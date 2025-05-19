#include "bsp.hpp"

namespace alh::bsp {

vec2_t edge_t::p(bsp_t const& bsp) const {
    vec2_t op = bsp.lines[lid].p;
    vec2_t oq = bsp.lines[lid].q;
    vec2_t p = op + (oq - op) * t1;
    return p;
}

vec2_t edge_t::q(bsp_t const& bsp) const {
    vec2_t op = bsp.lines[lid].p;
    vec2_t oq = bsp.lines[lid].q;
    vec2_t q = op + (oq - op) * t2;
    return q;
}

line_t edge_t::line(bsp_t const& bsp) const {
    return line_t{p(bsp), q(bsp)}.w_normal();
}

namespace {

struct build_context_t {
    bsp_t bsp;
    std::vector<edge_t> tmp;
};

bool split_line(bsp_t const& bsp, edge_t const& hyperplane, edge_t const& subj, edge_t &out1, edge_t &out2) {
    // split subj with hyperplane
    float alpha, beta;
    bool b_result;

    line_t const& l1 = bsp.lines[hyperplane.lid];
    line_t const& l2 = bsp.lines[subj.lid];

    if ((b_result = line_intersect_gg3(l1, l2, alpha, beta))) {
        if (b_result &= (subj.t1 < beta && beta < subj.t2)) {
            out1 = {subj.lid, subj.t1, beta};
            out2 = {subj.lid, beta, subj.t2};
        }
    }

    return b_result;
}

size_t build_impl(build_context_t &ctx, edge_t hyperplane, size_t i_begin, size_t i_end) {

    line_t lh = hyperplane.line(ctx.bsp);

    // split line segment if p and q are on opposite sides of hyperplane
    for (size_t i=i_begin; i != i_end && i != ctx.tmp.size(); i++) {
        if (ctx.tmp[i].p(ctx.bsp).is_left_of(lh) != ctx.tmp[i].q(ctx.bsp).is_left_of(lh)) {
            edge_t out1, out2;
            if (split_line(ctx.bsp, hyperplane, ctx.tmp[i], out1, out2)) {
                ctx.tmp[i] = out1;
                ctx.tmp.push_back(out2);
            }
        }
    }

    // then sort according to midpoints
    auto right_side = std::partition(ctx.tmp.begin() + i_begin, ctx.tmp.end(), [=](edge_t const& e) {
        vec2_t mid = (e.p(ctx.bsp) + e.q(ctx.bsp)) / 2.f;
        return mid.is_left_of(lh);
    });
    size_t i_split = std::distance(ctx.tmp.begin(), right_side);

    // create node
    size_t i_self = ctx.bsp.nodes.size();
    ctx.bsp.nodes.emplace_back();
    ctx.bsp.nodes[i_self].plane = hyperplane;

    // pop last line in tmp, (todo: select with heuristic and swap+pop)
    if (size_t i_last = ctx.tmp.size(); i_split < i_last) {
        edge_t right_split = ctx.tmp.back(); ctx.tmp.pop_back();
        ctx.bsp.nodes[i_self].right = build_impl(ctx, right_split, i_split, i_last);
    } else {
        ctx.bsp.nodes[i_self].right = EMPTY_LEAF;
    }

    if (i_begin < i_split) {
        edge_t left_split = ctx.tmp.back(); ctx.tmp.pop_back();
        ctx.bsp.nodes[i_self].left = build_impl(ctx, left_split, i_begin, i_split);
    } else {
        ctx.bsp.nodes[i_self].left = SOLID_LEAF;
    }

    return i_self;
}

bool sweep_impl(bsp_t const& bsp, size_t nid, line_t const& line, float t1, float t2, edge_t last_edge, vec2_t &out, line_t &out_line) {
    if (EMPTY_LEAF == nid) return false;
    if (SOLID_LEAF == nid) {
        out_line = line_t{last_edge.p(bsp), last_edge.q(bsp)}.w_normal();
        out = line.p + (line.q - line.p) * (t1 - 1e-4);
        return true;
    }

    node_t const& node = bsp.nodes[nid];
    line_t lh = node.plane.line(bsp);

    vec2_t p_t1 = line.p + (line.q - line.p) * t1;
    vec2_t p_t2 = line.p + (line.q - line.p) * t2;

    if (p_t1.is_left_of(lh) == p_t2.is_left_of(lh)) {
        size_t next = p_t1.is_left_of(lh) ? node.left : node.right;
        return sweep_impl(bsp, next, line, t1, t2, last_edge, out, out_line);
    } else { // split swept line (aka pass new t1 & t2)
        float t, s;
        line_intersect_gg3(line, lh, t, s);

        size_t first = p_t1.is_left_of(lh) ? node.left : node.right;
        size_t second = p_t1.is_left_of(lh) ? node.right : node.left;

        if (sweep_impl(bsp, first, line, t1, t, last_edge, out, out_line)) return true;
        return sweep_impl(bsp, second, line, t, t2, node.plane, out, out_line);
    }
}

// clip line against bsp tree recursively
bool clip_impl(clip_context_t &ctx, size_t nid, float t1, float t2) {
    
    if (SOLID_LEAF == nid) return ctx.on_solid(ctx, t1, t2, ctx.userdata);
    if (EMPTY_LEAF == nid) return ctx.on_empty(ctx, t1, t2, ctx.userdata);

    node_t const& node = (*ctx.bsp).nodes[nid];
    line_t lh = node.plane.line(*ctx.bsp);

    line_t l = (*ctx.bsp).lines[ctx.edge->lid];
    vec2_t p_t1 = l.p + (l.q - l.p) * t1;
    vec2_t p_t2 = l.p + (l.q - l.p) * t2;

    if (p_t1.is_left_of(lh) == p_t2.is_left_of(lh)) {
        size_t next = p_t1.is_left_of(lh) ? node.left : node.right;
        return clip_impl(ctx, next, t1, t2);
    } else { // split swept line (aka pass new t1 & t2)
        float t, s;
        line_intersect_gg3(l, lh, t, s);
        
        if (t1 < t && t < t2) {
            size_t first = p_t1.is_left_of(lh) ? node.left : node.right;
            size_t second = p_t1.is_left_of(lh) ? node.right : node.left;

            if (clip_impl(ctx, first, t1, t)) return true;
            return clip_impl(ctx, second, t, t2);
        } else {
            // use midpoint
            size_t next = ((p_t1 + p_t2) / 2.f).is_left_of(lh) ? node.left : node.right;
            return clip_impl(ctx, next, t1, t2);
        }
    }
}

bsp_t combine_buffers(bsp_t const& a, bsp_t const& b) {
    // place two bsp trees in a shared buffer (the trees are still disjoint)
    bsp_t c;
    c.lines.insert(c.lines.end(), a.lines.begin(), a.lines.end());
    c.lines.insert(c.lines.end(), b.lines.begin(), b.lines.end());
    c.nodes.insert(c.nodes.end(), a.nodes.begin(), a.nodes.end());
    for (node_t n : b.nodes) {
        n.plane.lid += a.lines.size();
        if (EMPTY_LEAF != n.left && SOLID_LEAF != n.left) n.left += a.nodes.size();
        if (EMPTY_LEAF != n.right && SOLID_LEAF != n.right) n.right += a.nodes.size();
        c.nodes.push_back(n);
    }
    return c;
}

bool cb_do_nothing(clip_context_t &, float, float, void *) { return false; }

bool cb_push_edge(clip_context_t &ctx, float t1, float t2, void *userdata) {
    std::vector<edge_t> &out = *(std::vector<edge_t> *)userdata;
    out.push_back({ctx.edge->lid, t1, t2});
    return false;
}

bool cb_push_flipped_edge(clip_context_t &ctx, float t1, float t2, void *userdata) {
    std::vector<edge_t> &out = *(std::vector<edge_t> *)userdata;
    out.push_back({ctx.edge->lid, t2, t1});
    return false;
}

} // namespace

bsp_t build(std::vector<line_t> lines, std::vector<edge_t> edges) {
    assert(lines.size() > 0);
    assert(edges.size() > 0);

    build_context_t ctx;
    ctx.tmp = std::move(edges);
    ctx.bsp.lines = std::move(lines);
    for (line_t &l : ctx.bsp.lines) l.w_normal();

    // select root hyperplane
    edge_t hyperplane = ctx.tmp.back(); ctx.tmp.pop_back();
    build_impl(ctx, hyperplane, 0, ctx.tmp.size());

    return ctx.bsp;
}

// build bsp-tree from lines only
bsp_t build(std::vector<line_t> lines) {
    assert(lines.size() > 0);

    // create edges
    std::vector<edge_t> edges;
    for (size_t i=0; i<lines.size(); i++)
        edges.push_back({i, 0.f, 1.f});

    return build(lines, edges);
}

bool is_solid(bsp_t const& bsp, size_t nid, vec2_t const& point) {
    // recurse until leaf
    if (EMPTY_LEAF == nid) return false;
    if (SOLID_LEAF == nid) return true;

    node_t const& node = bsp.nodes[nid];

    if (point.is_left_of(node.plane.line(bsp))) {
        return is_solid(bsp, node.left, point);
    } else {
        return is_solid(bsp, node.right, point);
    }
}

bool sweep(bsp_t const& bsp, line_t const& line, vec2_t &intersection, line_t &intersected) {
    return sweep_impl(bsp, 0, line, 0.f, 1.f, edge_t(), intersection, intersected);
}

//bool clip(clip_context_t &ctx, size_t root) {
//    return clip_impl(ctx, root, 0.f, 1.f);
//}

void dot_solve(bsp_t const& bsp, vec2_t const& p1, vec2_t &p2) {
    // repeat sweep and dot projection until p2 isn't solid
    vec2_t intersection;
    line_t line;
    while (sweep(bsp, {p1, p2}, intersection, line)) {
        p2 = project_to(p2, line) + line.normal * 0.1;
    }
}

bsp_t union_op(bsp_t const& a, bsp_t const& b) {
    
    assert(!(a.nodes.empty() && b.nodes.empty()));
    if (a.nodes.empty() ^ b.nodes.empty()) { // short circuit if empty operand
        return (a.nodes.empty()) ? b : a;
    }

    bsp_t c = combine_buffers(a, b);
    std::vector<edge_t> out;
    clip_context_t ctx;
    ctx.userdata = &out;
    ctx.on_solid = cb_do_nothing;
    ctx.on_empty = cb_push_edge;
    ctx.bsp = &c;
    size_t aroot = 0;
    size_t broot = a.nodes.size();

    for (auto it = c.nodes.begin() + a.nodes.size(); it != c.nodes.end(); ++it) {
        ctx.edge = &(*it).plane;
        clip_impl(ctx, aroot, (*it).plane.t1, (*it).plane.t2);
    }

    for (auto it = c.nodes.begin(); it != c.nodes.begin() + a.nodes.size(); ++it) {
        ctx.edge = &(*it).plane;
        clip_impl(ctx, broot, (*it).plane.t1, (*it).plane.t2);
    }
    
    return build(c.lines, out);
}

bsp_t intersect_op(bsp_t const& a, bsp_t const& b) {
    
    bsp_t c = combine_buffers(a, b);
    std::vector<edge_t> out;
    clip_context_t ctx;
    ctx.userdata = &out;
    ctx.on_solid = cb_push_edge;
    ctx.on_empty = cb_do_nothing;
    ctx.bsp = &c;
    size_t aroot = 0;
    size_t broot = a.nodes.size();

    for (auto it = c.nodes.begin() + a.nodes.size(); it != c.nodes.end(); ++it) {
        ctx.edge = &(*it).plane;
        clip_impl(ctx, aroot, (*it).plane.t1, (*it).plane.t2);
    }

    for (auto it = c.nodes.begin(); it != c.nodes.begin() + a.nodes.size(); ++it) {
        ctx.edge = &(*it).plane;
        clip_impl(ctx, broot, (*it).plane.t1, (*it).plane.t2);
    }
    
    return build(c.lines, out);
}

bsp_t difference_op(bsp_t const& a, bsp_t const& b) {
    
    bsp_t c = combine_buffers(a, b);
    std::vector<edge_t> out;
    clip_context_t ctx;
    ctx.userdata = &out;
    ctx.on_solid = cb_push_flipped_edge;
    ctx.on_empty = cb_do_nothing;
    ctx.bsp = &c;
    size_t aroot = 0;
    size_t broot = a.nodes.size();

    for (auto it = c.nodes.begin() + a.nodes.size(); it != c.nodes.end(); ++it) {
        ctx.edge = &(*it).plane;
        clip_impl(ctx, aroot, (*it).plane.t1, (*it).plane.t2);
    }

    ctx.on_solid = cb_do_nothing;
    ctx.on_empty = cb_push_edge;

    for (auto it = c.nodes.begin(); it != c.nodes.begin() + a.nodes.size(); ++it) {
        ctx.edge = &(*it).plane;
        clip_impl(ctx, broot, (*it).plane.t1, (*it).plane.t2);
    }
    
    return build(c.lines, out);
}

//bsp_t xor_op(bsp_t const& a, bsp_t const& b) {
//    bsp_t a_sub_b = difference_op(a, b);
//    bsp_t b_sub_a = difference_op(b, a);
//    return union_op(a_sub_b, b_sub_a);
//}

} // namespace alh::bsp