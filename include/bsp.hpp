#ifndef ALH_BSP_H
#define ALH_BSP_H

#include <cassert>
#include <cstdint>
#include <cmath>
#include <vector>
#include <algorithm>

struct line_t;

struct vec2_t {
    float x, y;
    friend vec2_t operator+(vec2_t const& lhs, vec2_t const& rhs) {return {lhs.x+rhs.x, lhs.y+rhs.y}; }
    friend vec2_t operator-(vec2_t const& lhs, vec2_t const& rhs) {return {lhs.x-rhs.x, lhs.y-rhs.y}; }
    friend vec2_t operator*(vec2_t const& lhs, float const& rhs) {return {lhs.x*rhs, lhs.y*rhs}; }
    friend vec2_t operator/(vec2_t const& lhs, float const& rhs) {return {lhs.x/rhs, lhs.y/rhs}; }
    bool is_left_of(line_t l) const;
    bool is_exactly_on(line_t l) const;
};

float dist2(vec2_t const& p, vec2_t const& q) {
    vec2_t pq = q - p;
    return (pq.x*pq.x + pq.y*pq.y);
}

vec2_t get_normal(vec2_t p, vec2_t q) {
    vec2_t pq = q - p;
    float dist = sqrtf(pq.x*pq.x + pq.y*pq.y);
    assert(0.f != dist);
    return {-pq.y / dist, pq.x / dist};
}

struct line_t {
    vec2_t p, q;
    vec2_t normal;
    line_t &w_normal() { normal = get_normal(p, q); return *this; }
};

bool vec2_t::is_left_of(line_t l) const {
    vec2_t dl = l.q - l.p;
    return (x - l.p.x)*(-dl.y) + (y - l.p.y)*(dl.x) < 0;
}

bool vec2_t::is_exactly_on(line_t l) const {
    vec2_t dl = l.q - l.p;
    return (x - l.p.x)*(-dl.y) + (y - l.p.y)*(dl.x) == 0;
}

static constexpr size_t EMPTY_LEAF = (size_t)(-1);
static constexpr size_t SOLID_LEAF = (size_t)(-2);

struct node_t {
    line_t plane;
    size_t right, left;
};

bool line_intersect_gg3(line_t const& l1, line_t const& l2, float &alpha, float &beta) {
    vec2_t a, b, c;
    a = l1.q - l1.p;
    b = l2.p - l2.q;
    c = l1.p - l2.p;

    float num_a, denom_a;
    num_a = b.y*c.x - b.x*c.y;
    denom_a = a.y*b.x - a.x*b.y;

    if (denom_a == 0) {
        return false; // collinear segments
    } else {
        alpha = num_a / denom_a;
        beta = (a.x*c.y - a.y*c.x) / denom_a;
        return true;
    }
}

using bsp_t = typename std::vector<node_t>;

struct build_context_t {
    bsp_t nodes;
    std::vector<line_t> tmp;
};

// todo: maybe turn into method for line_t
bool split_line(line_t const& hyperplane, line_t const& subj, line_t &out1, line_t &out2) {
    // split subj with hyperplane
    float alpha, beta;
    bool b_result;
    
    if (b_result = line_intersect_gg3(hyperplane, subj, alpha, beta)) {
        vec2_t i_pq = subj.p + (subj.q - subj.p) * beta;

        float eps = 0.1;
        if (b_result &= (dist2(subj.p, i_pq) > eps*eps && dist2(subj.q, i_pq) > eps*eps)) {
            out1 = {subj.p, i_pq};
            out2 = {i_pq, subj.q};
        }
    }

    return b_result;
}

size_t build_helper(build_context_t &ctx, line_t hyperplane, size_t i_begin, size_t i_end) {

    // split line segment if p and q are on opposite sides of hyperplane
    for (size_t i=i_begin; i != i_end; i++) {
        if (ctx.tmp[i].p.is_left_of(hyperplane) != ctx.tmp[i].q.is_left_of(hyperplane)) {
            line_t out1, out2;
            if (split_line(hyperplane, ctx.tmp[i], out1, out2)) {
                ctx.tmp[i] = out1;
                ctx.tmp.push_back(out2);
            }
        }
    }

    // then sort according to midpoints
    auto right_side = std::partition(ctx.tmp.begin() + i_begin, ctx.tmp.end(), [=](line_t const& l) {
        vec2_t mid = (l.p + l.q) / 2.f;
        return mid.is_left_of(hyperplane);
    });
    size_t i_split = std::distance(ctx.tmp.begin(), right_side);

    // create node
    size_t i_self = ctx.nodes.size();
    ctx.nodes.emplace_back();
    ctx.nodes[i_self].plane = hyperplane.w_normal();

    // pop last line in tmp, (todo: select with heuristic and swap+pop)
    if (size_t i_last = ctx.tmp.size(); i_split < i_last) {
        line_t right_split = ctx.tmp.back(); ctx.tmp.pop_back();
        ctx.nodes[i_self].right = build_helper(ctx, right_split, i_split, i_last);
    } else {
        ctx.nodes[i_self].right = EMPTY_LEAF;
    }

    if (i_begin < i_split) {
        line_t left_split = ctx.tmp.back(); ctx.tmp.pop_back();
        ctx.nodes[i_self].left = build_helper(ctx, left_split, i_begin, i_split);
    } else {
        ctx.nodes[i_self].left = SOLID_LEAF;
    }

    return i_self;
}

// build bsp-tree from lines
bsp_t build(std::vector<line_t> lines) {
    build_context_t ctx;
    ctx.tmp = std::move(lines);

    // select root hyperplane
    line_t hyperplane = ctx.tmp.back(); ctx.tmp.pop_back();
    build_helper(ctx, hyperplane, 0, ctx.tmp.size());

    return ctx.nodes;
}

// split tree in-place, return indices to new tree's root (?)

//size_t split_tree(bsp_t &bsp, size_t nid, line_t const& hyperplane) {
//
//    // detect leaves
//    if (EMPTY_LEAF == nid || SOLID_LEAF == nid) return nid;
//
//    node_t &node = bsp[nid];
//
////    if ((bsp[nid].plane.p.is_left_of(hyperplane) || bsp[nid].plane.p.is_exactly_on(hyperplane))
////    && (bsp[nid].plane.q.is_left_of(hyperplane) || bsp[nid].plane.q.is_exactly_on(hyperplane))) {
////        node.left = split_tree(bsp, node.left, hyperplane);
////        node.right = split_tree(bsp, node.right, hyperplane);
////        return EMPTY_LEAF;
////    } else if ((!bsp[nid].plane.p.is_left_of(hyperplane) || bsp[nid].plane.p.is_exactly_on(hyperplane))
////    && (!bsp[nid].plane.q.is_left_of(hyperplane) || bsp[nid].plane.q.is_exactly_on(hyperplane))) {
////        node.left = split_tree(bsp, node.left, hyperplane);
////        node.right = split_tree(bsp, node.right, hyperplane);
////        return nid;
////    }
//    // traverse down the tree if hyperplane is on one side of current node
//    if (bsp[nid].plane.p.is_left_of(hyperplane) == bsp[nid].plane.q.is_left_of(hyperplane)
//    || bsp[nid].plane.p.is_exactly_on(hyperplane)
//    || bsp[nid].plane.q.is_exactly_on(hyperplane)) {
//        split_tree(bsp, node.left, hyperplane);
//        split_tree(bsp, node.right, hyperplane);
//        return nid;
//    } else {
//        // split the segment
//        line_t subj = bsp[nid].plane;
//        line_t out1, out2;
//        split_line(hyperplane, subj, out1, out2);
//
//        // assert that out1 and out2 each have a point exactly on hyperplane
//        assert(out1.p.is_exactly_on(hyperplane) || out1.q.is_exactly_on(hyperplane));
//        assert(out2.p.is_exactly_on(hyperplane) || out2.q.is_exactly_on(hyperplane));
//
//        //line_t ll, lr;
//        //ll = bsp[nid].plane.p.is_left_of(hyperplane) ? out1 : out2;
//        //lr = bsp[nid].plane.p.is_left_of(hyperplane) ? out2 : out1;
//
//        node.plane = out1;
//
//        bsp.push_back(node); // push copy
//        bsp.back().plane = out2;
//        
//        node.right = split_tree(bsp, bsp.size() - 1, hyperplane);
//        node.left = split_tree(bsp, node.left, hyperplane);
//        
//        return nid;
//    }
//}

//// merge two bsp trees
//// descend down the bsp_a, partition bsp_b along the way
//void merge(bsp_t const& bsp_a, size_t const& nid, bsp_t const& bsp_b) {
//    
//}

// write traversal functions
bool is_solid(bsp_t const& bsp, size_t const& nid, vec2_t point) {
    // recurse until leaf
    if (EMPTY_LEAF == nid) return false;
    if (SOLID_LEAF == nid) return true;

    node_t const& node = bsp[nid];
    if (point.is_left_of(node.plane)) {
        return is_solid(bsp, node.left, point);
    } else {
        return is_solid(bsp, node.right, point);
    }
}

bool sweep_impl(bsp_t const& bsp, size_t const& nid, vec2_t const& p1, vec2_t const& p2, float t1, float t2, line_t last_line, vec2_t &out, line_t &out_line) {
    if (EMPTY_LEAF == nid) return false;
    if (SOLID_LEAF == nid) {
        out_line = last_line;
        out = p1 + (p2 - p1) * (t1 - 1e-4);
        return true;
    }

    node_t const& node = bsp[nid];
    line_t hyperplane = node.plane;

    vec2_t p_s1 = p1 + (p2 - p1) * t1;
    vec2_t p_s2 = p1 + (p2 - p1) * t2;

    if (p_s1.is_left_of(hyperplane) == p_s2.is_left_of(hyperplane)) {
        size_t next = p_s1.is_left_of(hyperplane) ? node.left : node.right;
        return sweep_impl(bsp, next, p1, p2, t1, t2, last_line, out, out_line);
    } else { // split swept line (aka pass new t1 & t2)
        float t, s;
        line_intersect_gg3({p1, p2}, hyperplane, t, s);

        size_t first = p_s1.is_left_of(hyperplane) ? node.left : node.right;
        size_t second = p_s1.is_left_of(hyperplane) ? node.right : node.left;

        if (sweep_impl(bsp, first, p1, p2, t1, t, last_line, out, out_line)) return true;
        return sweep_impl(bsp, second, p1, p2, t, t2, hyperplane, out, out_line);
    }
}

bool sweep(bsp_t const& bsp, line_t const& line, vec2_t &intersection, line_t &intersected) {
    return sweep_impl(bsp, 0, line.p, line.q, 0.f, 1.f, line_t(), intersection, intersected);
}

vec2_t project_to(vec2_t const& p, line_t const& line) {
    float p_a_x, p_a_y, b_a_x, b_a_y;
    p_a_x = p.x - line.p.x;
    p_a_y = p.y - line.p.y;
    b_a_x = line.q.x - line.p.x;
    b_a_y = line.q.y - line.p.y;

    float num = p_a_x*b_a_x + p_a_y*b_a_y;
    float denom = b_a_x*b_a_x + b_a_y*b_a_y;
    float s = num/denom;
    
    return {line.p.x + b_a_x * s, line.p.y + b_a_y * s};
}

void dot_solve(bsp_t const& bsp, vec2_t const& p1, vec2_t &p2) {
    // repeat sweep and dot projection until p2 isn't solid
    vec2_t intersection;
    line_t line;
    while (sweep(bsp, {p1, p2}, intersection, line)) {
        p2 = project_to(p2, line) + line.normal * 0.1;
    }
}

// todo:
//  - implement merge function
//  - use merge to implement union operation

//  - have a look at new heuristic
//  - check if implementation can use TCO

#endif