#include "bsp.hpp"

namespace alh::bsp {

namespace {

struct build_context_t {
    uint32_t leaf_id_acc;
    bsp_t nodes;
    std::vector<paramline_t> tmp;
};

bool split_line(paramline_t const& hyperplane, paramline_t const& subj, paramline_t &out1, paramline_t &out2) {
    // split subj with hyperplane
    float alpha, beta;
    bool b_result;
    
    line_t const& l1 = hyperplane.line;
    line_t const& l2 = subj.line; // original line for finding new t1 & t2
    line_t const& l3 = subj.apply(); // scaled line for comparing absolute distance

    if ((b_result = line_intersect_gg3(l1, l2, alpha, beta))) {
        vec2_t i_pq = l2.p + (l2.q - l2.p) * beta;
        constexpr float eps2 = 0.1 * 0.1; // minimum allowed deviation in position
        if (b_result &= (subj.t1 < beta == beta < subj.t2 && dist2(l3.p, i_pq) > eps2 && dist2(l3.q, i_pq) > eps2)) {
            out1 = {subj.line, subj.t1, beta};
            out2 = {subj.line, beta, subj.t2};
        }
    }

    return b_result;
}

id_t build_impl(build_context_t &ctx, paramline_t hyperplane, id_t i_begin, id_t i_end) {
    
    line_t lh = hyperplane.apply();
    
    // split line segment if p and q are on opposite sides of hyperplane
    for (id_t i=i_begin; i != i_end && i != ctx.tmp.size(); i++) {
        line_t l = ctx.tmp[i].apply();
        if (l.p.is_left_of(lh) != l.q.is_left_of(lh)) {
            paramline_t out1, out2;
            if (split_line(hyperplane, ctx.tmp[i], out1, out2)) {
                ctx.tmp[i] = out1;
                ctx.tmp.push_back(out2);
            }
        }
    }

    // then sort according to midpoints
    auto right_side = std::partition(ctx.tmp.begin() + i_begin, ctx.tmp.end(), [=](paramline_t const& pl) {
        line_t l = pl.apply();
        vec2_t mid = (l.p + l.q) / 2.f;
        mid = mid - l.w_normal().normal; // yeah you could consider me a bit of a hacker
        return mid.is_left_of(lh);
    });
    id_t i_split = std::distance(ctx.tmp.begin(), right_side);

    // create node
    id_t i_self = ctx.nodes.size();
    ctx.nodes.push_back({hyperplane, 0, 0});

    // pop last line in tmp, (todo: select with heuristic and swap+pop)
    if (id_t i_last = ctx.tmp.size(); i_split < i_last) {
        paramline_t right_split = ctx.tmp.back(); ctx.tmp.pop_back();
        ctx.nodes[i_self].right = build_impl(ctx, right_split, i_split, i_last);
    } else {
        ctx.nodes[i_self].right = ((ctx.leaf_id_acc++) | IS_LEAF) & ~IS_SOLID;
    }

    if (i_begin < i_split) {
        paramline_t left_split = ctx.tmp.back(); ctx.tmp.pop_back();
        ctx.nodes[i_self].left = build_impl(ctx, left_split, i_begin, i_split);
    } else {
        ctx.nodes[i_self].left = ((ctx.leaf_id_acc++) | IS_LEAF) | IS_SOLID;
    }

    return i_self;
}

bool sweep_impl(bsp_t const& bsp, id_t nid, line_t const& line, float t1, float t2, paramline_t last_line, vec2_t &out, line_t &out_line) {
    if (empty_leaf(nid)) return false;
    if (solid_leaf(nid)) {
        out_line = last_line.apply().w_normal();
        out = line.p + (line.q - line.p) * (t1 - 1e-4);
        return true;
    }

    bsp_node_t const& node = bsp[nid];
    line_t lh = node.plane.apply();

    vec2_t p_t1 = line.p + (line.q - line.p) * t1;
    vec2_t p_t2 = line.p + (line.q - line.p) * t2;

    if (p_t1.is_left_of(lh) == p_t2.is_left_of(lh)) {
        id_t next = p_t1.is_left_of(lh) ? node.left : node.right;
        return sweep_impl(bsp, next, line, t1, t2, last_line, out, out_line);
    } else { // split swept line (aka pass new t1 & t2)
        float t, s;
        line_intersect_gg3(line, lh, t, s);

        id_t first = p_t1.is_left_of(lh) ? node.left : node.right;
        id_t second = p_t1.is_left_of(lh) ? node.right : node.left;

        if (sweep_impl(bsp, first, line, t1, t, last_line, out, out_line)) return true;
        return sweep_impl(bsp, second, line, t, t2, node.plane, out, out_line);
    }
}

// clip line against bsp tree recursively
bool clip_impl(clip_context_t &ctx, id_t nid, float t1, float t2) {
    
    if (solid_leaf(nid)) return ctx.on_solid(ctx, t1, t2, ctx.userdata);
    if (empty_leaf(nid)) return ctx.on_empty(ctx, t1, t2, ctx.userdata);

    bsp_node_t const& node = (*ctx.bsp)[nid];
    line_t lh = node.plane.apply();

    line_t l = (*ctx.paramline).line;
    vec2_t p_t1 = l.p + (l.q - l.p) * t1;
    vec2_t p_t2 = l.p + (l.q - l.p) * t2;

    if (p_t1.is_left_of(lh) == p_t2.is_left_of(lh)) {
        id_t next = p_t1.is_left_of(lh) ? node.left : node.right;
        return clip_impl(ctx, next, t1, t2);
    } else { // split swept line (aka pass new t1 & t2)
        float t, s;
        line_intersect_gg3(l, lh, t, s);

        vec2_t p_t = l.p + (l.q - l.p) * t;
        constexpr float eps2 = 0.1 * 0.1; // minimum allowed deviation in position

        if (t1 < t == t < t2 && dist2(p_t, p_t1) > eps2 && dist2(p_t, p_t2) > eps2) {
            assert(t != t1 && t != t2);
            id_t first = p_t1.is_left_of(lh) ? node.left : node.right;
            id_t second = p_t1.is_left_of(lh) ? node.right : node.left;

            if (clip_impl(ctx, first, t1, t)) return true;
            return clip_impl(ctx, second, t, t2);
        } else {
            // use midpoint
            id_t next = ((p_t1 + p_t2) / 2.f).is_left_of(lh) ? node.left : node.right;
            return clip_impl(ctx, next, t1, t2);
        }
    }
}

bool cb_do_nothing(clip_context_t &, float, float, void *) { return false; }

bool cb_push_segment(clip_context_t &ctx, float t1, float t2, void *userdata) {
    std::vector<paramline_t> &out = *(std::vector<paramline_t> *)userdata;
    out.push_back(paramline_t(ctx.paramline->line, t1, t2));
    return false;
}

bool cb_push_flipped_segment(clip_context_t &ctx, float t1, float t2, void *userdata) {
    std::vector<paramline_t> &out = *(std::vector<paramline_t> *)userdata;
    out.push_back(paramline_t(ctx.paramline->line, t2, t1));
    return false;
}

} // namespace

// build bsp-tree from lines
bsp_t build(std::vector<line_t> lines) {
    assert(lines.size() > 0);

    std::vector<paramline_t> paramlines;
    paramlines.insert(paramlines.end(), lines.begin(), lines.end());

    return build(paramlines);
}

bsp_t build(std::vector<paramline_t> paramlines) {
    assert(paramlines.size() > 0);

    build_context_t ctx;
    ctx.leaf_id_acc = 0;
    ctx.tmp = std::move(paramlines);

    // select root hyperplane
    paramline_t hyperplane = ctx.tmp.back(); ctx.tmp.pop_back();
    build_impl(ctx, hyperplane, 0, ctx.tmp.size());

    return ctx.nodes;
}

bool is_solid(bsp_t const& bsp, id_t nid, vec2_t const& point) {
    // recurse until leaf
    if (empty_leaf(nid)) return false;
    if (solid_leaf(nid)) return true;

    bsp_node_t const& node = bsp[nid];
    if (point.is_left_of(node.plane.apply())) {
        return is_solid(bsp, node.left, point);
    } else {
        return is_solid(bsp, node.right, point);
    }
}

id_t leaf_id(bsp_t const& bsp, id_t nid, vec2_t const& point) {
    if (is_leaf(nid)) return (nid & ~IS_LEAF) & ~IS_SOLID;

    bsp_node_t const& node = bsp[nid];
    if (point.is_left_of(node.plane.apply())) {
        return leaf_id(bsp, node.left, point);
    } else {
        return leaf_id(bsp, node.right, point);
    }
}

bool sweep(bsp_t const& bsp, line_t const& line, vec2_t &intersection, line_t &intersected) {
    return sweep_impl(bsp, 0, line, 0.f, 1.f, paramline_t(), intersection, intersected);
}

void dot_solve(bsp_t const& bsp, vec2_t const& p1, vec2_t &p2) {
    // repeat sweep and dot projection until p2 isn't solid
    vec2_t intersection;
    line_t line;
    while (sweep(bsp, {p1, p2}, intersection, line)) {
        if (p2.is_left_of(line))
            p2 = project_to(p2, line) + line.normal * 0.1;
        else
            p2 = project_to(p2, line) - line.normal * 0.1;
    }
}

bsp_t union_op(bsp_t const& a, bsp_t const& b) {
    
    assert(!(a.empty() && b.empty()));
    if (a.empty() ^ b.empty()) { // short circuit if empty operand
        return (a.empty()) ? b : a;
    }

    std::vector<paramline_t> out;
    clip_context_t ctx;
    ctx.userdata = &out;
    ctx.on_solid = cb_do_nothing;
    ctx.on_empty = cb_push_segment;
    ctx.bsp = &a;

    for (bsp_node_t n : b) {
        ctx.paramline = &n.plane;
        clip_impl(ctx, 0, n.plane.t1, n.plane.t2);
    }

    ctx.bsp = &b;

    for (bsp_node_t n : a) {
        ctx.paramline = &n.plane;
        clip_impl(ctx, 0, n.plane.t1, n.plane.t2);
    }
    
    return build(out);
}

bsp_t intersect_op(bsp_t const& a, bsp_t const& b) {

    std::vector<paramline_t> out;
    clip_context_t ctx;
    ctx.userdata = &out;
    ctx.on_solid = cb_push_segment;
    ctx.on_empty = cb_do_nothing;
    ctx.bsp = &a;

    for (bsp_node_t n : b) {
        ctx.paramline = &n.plane;
        clip_impl(ctx, 0, n.plane.t1, n.plane.t2);
    }

    ctx.bsp = &b;

    for (bsp_node_t n : a) {
        ctx.paramline = &n.plane;
        clip_impl(ctx, 0, n.plane.t1, n.plane.t2);
    }
    
    return build(out);
}

bsp_t difference_op(bsp_t const& a, bsp_t const& b) {

    std::vector<paramline_t> out;
    clip_context_t ctx;
    ctx.userdata = &out;
    ctx.on_solid = cb_push_flipped_segment;
    ctx.on_empty = cb_do_nothing;
    ctx.bsp = &a;

    for (bsp_node_t n : b) {
        ctx.paramline = &n.plane;
        clip_impl(ctx, 0, n.plane.t1, n.plane.t2);
    }

    ctx.on_solid = cb_do_nothing;
    ctx.on_empty = cb_push_segment;
    ctx.bsp = &b;

    for (bsp_node_t n : a) {
        ctx.paramline = &n.plane;
        clip_impl(ctx, 0, n.plane.t1, n.plane.t2);
    }
    
    return build(out);
}

bsp_t xor_op(bsp_t const& a, bsp_t const& b) {
    bsp_t a_sub_b = difference_op(a, b);
    bsp_t b_sub_a = difference_op(b, a);
    return union_op(a_sub_b, b_sub_a);
}

} // namespace alh::bsp