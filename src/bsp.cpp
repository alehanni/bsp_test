#include "bsp.hpp"

namespace alh::bsp {

namespace {

struct build_context_t {
    bsp_t nodes;
    std::vector<line_t> tmp;
};

// todo: maybe turn into method for line_t
bool split_line(line_t const& hyperplane, line_t const& subj, line_t &out1, line_t &out2) {
    // split subj with hyperplane
    float alpha, beta;
    bool b_result;
    
    if ((b_result = line_intersect_gg3(hyperplane, subj, alpha, beta))) {
        vec2_t i_pq = subj.p + (subj.q - subj.p) * beta;

        float eps = 0.1;
        if (b_result &= (dist2(subj.p, i_pq) > eps*eps && dist2(subj.q, i_pq) > eps*eps)) {
            out1 = {subj.p, i_pq, subj.userdata};
            out2 = {i_pq, subj.q, subj.userdata};
        }
    }

    return b_result;
}

size_t build_impl(build_context_t &ctx, line_t hyperplane, size_t i_begin, size_t i_end) {
    // split line segment if p and q are on opposite sides of hyperplane
    for (size_t i=i_begin; i != i_end && i != ctx.tmp.size(); i++) {
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
        ctx.nodes[i_self].right = build_impl(ctx, right_split, i_split, i_last);
    } else {
        ctx.nodes[i_self].right = EMPTY_LEAF;
    }

    if (i_begin < i_split) {
        line_t left_split = ctx.tmp.back(); ctx.tmp.pop_back();
        ctx.nodes[i_self].left = build_impl(ctx, left_split, i_begin, i_split);
    } else {
        ctx.nodes[i_self].left = SOLID_LEAF;
    }

    return i_self;
}

bool sweep_impl(bsp_t const& bsp, size_t nid, line_t const& line, float t1, float t2, line_t last_line, vec2_t &out, line_t &out_line) {
    if (EMPTY_LEAF == nid) return false;
    if (SOLID_LEAF == nid) {
        out_line = last_line;
        out = line.p + (line.q - line.p) * (t1 - 1e-4);
        return true;
    }

    node_t const& node = bsp[nid];
    line_t hyperplane = node.plane;

    vec2_t p_t1 = line.p + (line.q - line.p) * t1;
    vec2_t p_t2 = line.p + (line.q - line.p) * t2;

    if (p_t1.is_left_of(hyperplane) == p_t2.is_left_of(hyperplane)) {
        size_t next = p_t1.is_left_of(hyperplane) ? node.left : node.right;
        return sweep_impl(bsp, next, line, t1, t2, last_line, out, out_line);
    } else { // split swept line (aka pass new t1 & t2)
        float t, s;
        line_intersect_gg3(line, hyperplane, t, s);

        size_t first = p_t1.is_left_of(hyperplane) ? node.left : node.right;
        size_t second = p_t1.is_left_of(hyperplane) ? node.right : node.left;

        if (sweep_impl(bsp, first, line, t1, t, last_line, out, out_line)) return true;
        return sweep_impl(bsp, second, line, t, t2, hyperplane, out, out_line);
    }
}

// clip line against bsp tree recursively
bool clip_impl(clip_context_t &ctx, size_t nid, float t1, float t2) {
    
    if (SOLID_LEAF == nid) return ctx.on_solid(ctx, t1, t2, ctx.userdata);
    if (EMPTY_LEAF == nid) return ctx.on_empty(ctx, t1, t2, ctx.userdata);

    node_t const& node = (*ctx.bsp)[nid];
    line_t hyperplane = node.plane;

    vec2_t p_t1 = ctx.line->p + (ctx.line->q - ctx.line->p) * t1;
    vec2_t p_t2 = ctx.line->p + (ctx.line->q - ctx.line->p) * t2;

    if (p_t1.is_left_of(hyperplane) == p_t2.is_left_of(hyperplane)) {
        size_t next = p_t1.is_left_of(hyperplane) ? node.left : node.right;
        return clip_impl(ctx, next, t1, t2);
    } else { // split swept line (aka pass new t1 & t2)
        float t, s;
        line_intersect_gg3(*ctx.line, hyperplane, t, s);

        vec2_t p_t = ctx.line->p + (ctx.line->q - ctx.line->p) * t;
        
        float eps = 0.1;
        if (dist2(p_t, p_t1) > eps*eps && dist2(p_t, p_t2) > eps*eps) {
            size_t first = p_t1.is_left_of(hyperplane) ? node.left : node.right;
            size_t second = p_t1.is_left_of(hyperplane) ? node.right : node.left;

            if (clip_impl(ctx, first, t1, t)) return true;
            return clip_impl(ctx, second, t, t2);
        } else {
            // use midpoint
            size_t next = ((p_t1 + p_t2) / 2.f).is_left_of(hyperplane) ? node.left : node.right;
            return clip_impl(ctx, next, t1, t2);
        }
    }
}

bool cb_do_nothing(clip_context_t &, float, float, void *) { return false; }

bool cb_push_segment(clip_context_t &ctx, float t1, float t2, void *userdata) {
    vec2_t p_t1 = ctx.line->p + (ctx.line->q - ctx.line->p) * t1;
    vec2_t p_t2 = ctx.line->p + (ctx.line->q - ctx.line->p) * t2;
    std::vector<line_t> &out = *(std::vector<line_t> *)userdata;
    out.push_back({p_t1, p_t2, ctx.line->userdata});
    return false;
}

bool cb_push_flipped_segment(clip_context_t &ctx, float t1, float t2, void *userdata) {
    vec2_t p_t1 = ctx.line->p + (ctx.line->q - ctx.line->p) * t1;
    vec2_t p_t2 = ctx.line->p + (ctx.line->q - ctx.line->p) * t2;
    std::vector<line_t> &out = *(std::vector<line_t> *)userdata;
    out.push_back({p_t2, p_t1, ctx.line->userdata});
    return false;
}

} // namespace

// build bsp-tree from lines
bsp_t build(std::vector<line_t> lines) {
    assert(lines.size() > 0);

    build_context_t ctx;
    ctx.tmp = std::move(lines);

    // select root hyperplane
    line_t hyperplane = ctx.tmp.back(); ctx.tmp.pop_back();
    build_impl(ctx, hyperplane, 0, ctx.tmp.size());

    return ctx.nodes;
}

bool is_solid(bsp_t const& bsp, size_t nid, vec2_t const& point) {
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

bool sweep(bsp_t const& bsp, line_t const& line, vec2_t &intersection, line_t &intersected) {
    return sweep_impl(bsp, 0, line, 0.f, 1.f, line_t(), intersection, intersected);
}

bool clip(clip_context_t &ctx, size_t root) {
    return clip_impl(ctx, root, 0.f, 1.f);
}

void dot_solve(bsp_t const& bsp, vec2_t const& p1, vec2_t &p2) {
    // repeat sweep and dot projection until p2 isn't solid
    vec2_t intersection;
    line_t line;
    while (sweep(bsp, {p1, p2}, intersection, line)) {
        p2 = project_to(p2, line) + line.normal * 0.1;
    }
}

bsp_t union_op(bsp_t const& a, bsp_t const& b) {
    
    assert(!(a.empty() && b.empty()));
    if (a.empty() ^ b.empty()) { // short circuit if empty operand
        return (a.empty()) ? b : a;
    }

    std::vector<line_t> out;
    clip_context_t ctx;
    ctx.userdata = &out;
    ctx.on_solid = cb_do_nothing;
    ctx.on_empty = cb_push_segment;
    ctx.bsp = &a;

    for (node_t node : b) {
        ctx.line = &node.plane;
        clip_impl(ctx, 0, 0.f, 1.f);
    }

    ctx.bsp = &b;

    for (node_t node : a) {
        ctx.line = &node.plane;
        clip_impl(ctx, 0, 0.f, 1.f);
    }
    
    return build(out);
}

bsp_t intersect_op(bsp_t const& a, bsp_t const& b) {
    std::vector<line_t> out;

    clip_context_t ctx;
    ctx.userdata = &out;
    ctx.on_solid = cb_push_segment;
    ctx.on_empty = cb_do_nothing;
    ctx.bsp = &a;

    for (node_t node : b) {
        ctx.line = &node.plane;
        clip_impl(ctx, 0, 0.f, 1.f);
    }

    ctx.bsp = &b;

    for (node_t node : a) {
        ctx.line = &node.plane;
        clip_impl(ctx, 0, 0.f, 1.f);
    }
    
    return build(out);
}

bsp_t difference_op(bsp_t const& a, bsp_t const& b) {
    std::vector<line_t> out;

    clip_context_t ctx;
    ctx.userdata = &out;
    ctx.on_solid = cb_push_flipped_segment;
    ctx.on_empty = cb_do_nothing;
    ctx.bsp = &a;

    for (node_t node : b) {
        ctx.line = &node.plane;
        clip_impl(ctx, 0, 0.f, 1.f);
    }

    ctx.on_solid = cb_do_nothing;
    ctx.on_empty = cb_push_segment;
    ctx.bsp = &b;

    for (node_t node : a) {
        ctx.line = &node.plane;
        clip_impl(ctx, 0, 0.f, 1.f);
    }
    
    return build(out);
}

bsp_t xor_op(bsp_t const& a, bsp_t const& b) {
    bsp_t a_sub_b = difference_op(a, b);
    bsp_t b_sub_a = difference_op(b, a);
    return union_op(a_sub_b, b_sub_a);
}

} // namespace alh::bsp