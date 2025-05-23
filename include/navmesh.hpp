#ifndef ALH_NAVMESH_H
#define ALH_NAVMESH_H

#include <cassert>
#include <limits>
#include <vector>

#include "alh.hpp"
#include "bsp.hpp"
#include <cstdio> // NOTE: REMOVE LATER
#include "dbg_shapes.hpp" // NOTE: REMOVE LATER

namespace alh::bsp::navmesh {

struct traj_t {
    uint64_t mask; // mask corresponds to depth
    uint64_t bits; // sequence of left/right traversal choices
};

void get_traj_impl(bsp_t const& bsp, size_t nid, size_t depth, vec2_t const& point, traj_t &tr) {
    assert(depth < 64);

    if (nid == EMPTY_LEAF || nid == SOLID_LEAF) return;
    assert(nid < bsp.size());

    tr.mask = ~(0xffffffffffffffff << (depth + 1));

    node_t const& node = bsp[nid];
    if (point.is_left_of(node.plane.apply())) {
        return get_traj_impl(bsp, node.left, depth + 1, point, tr);
    } else {
        tr.bits |= 1 << depth;
        return get_traj_impl(bsp, node.right, depth + 1, point, tr);
    }
}

traj_t get_traj(bsp_t const& bsp, vec2_t const& point) {
    traj_t tr = { 0 };
    get_traj_impl(bsp, 0, 0, point, tr);
    return tr;
}

void bsp_bb(bsp_t const& bsp, vec2_t &min, vec2_t &max) {
    static constexpr float inf = std::numeric_limits<float>::infinity();
    
    min = {inf, inf};
    max = {-inf, -inf};
    
    for (node_t node : bsp) {
        line_t l = node.plane.apply();
        min.x = std::min(std::min(min.x, l.p.x), l.q.x);
        min.y = std::min(std::min(min.y, l.p.y), l.q.y);
        max.x = std::max(std::max(max.x, l.p.x), l.q.x);
        max.y = std::max(std::max(max.y, l.p.y), l.q.y);
    }
}

bsp_t bsp_poly(bsp_t const& bsp, size_t root, traj_t tr) {

    bsp_t poly;
    node_t node = bsp[root];
    uint64_t mask = 1;

    size_t nid = root;

    while(tr.mask & mask) {

        // todo: replace scaling with bb intersection points
        paramline_t longplane = node.plane;
        longplane.t1 = -1000.f;
        longplane.t2 = 1000.f;
        assert(EMPTY_LEAF != nid);
        assert(SOLID_LEAF != nid);
        longplane.line.userdata = (void *)nid;

        if (tr.bits & mask) {
            if (1 == mask)
                poly = build({longplane.flip()});
            else
                poly = bsp::difference_op(poly, build({longplane}));
            
            if (EMPTY_LEAF == node.right || SOLID_LEAF == node.right) break;
            nid = node.right;
            node = bsp[node.right];
        } else {
            if (1 == mask)
                poly = build({longplane});
            else
                poly = bsp::difference_op(poly, build({longplane.flip()}));
            
            if (EMPTY_LEAF == node.left || SOLID_LEAF == node.left) break;
            nid = node.left;
            node = bsp[node.left];
        }

        mask = mask << 1;
    }

    return poly;
}

} // namespace alh::bsp::navmesh

#endif