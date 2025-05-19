#ifndef ALH_BSP_H
#define ALH_BSP_H

#include <cassert>
#include <cstdint>
#include <cmath>
#include <vector>
#include <algorithm>

#include "alh.hpp"

namespace alh::bsp {

    static constexpr size_t EMPTY_LEAF = (size_t)(-1);
    static constexpr size_t SOLID_LEAF = (size_t)(-2);

    struct bsp_t;

    struct edge_t {
        size_t lid;
        float t1, t2;
        vec2_t p(bsp_t const& bsp) const;
        vec2_t q(bsp_t const& bsp) const;
        line_t line(bsp_t const& bsp) const;
    };
    
    struct node_t {
        edge_t plane;
        size_t right, left;
    };

    struct bsp_t {
        std::vector<line_t> lines;
        std::vector<node_t> nodes;
    };

    struct clip_context_t;
    typedef bool (*leaf_callback)(clip_context_t &ctx, float t1, float t2, void *userdata);

    struct clip_context_t {
        const bsp_t *bsp;
        const edge_t *edge;
        leaf_callback on_empty;
        leaf_callback on_solid;
        void *userdata; // used for output, etc.
    };

    bsp_t build(std::vector<line_t> lines);
    
    bool is_solid(bsp_t const& bsp, size_t nid, vec2_t const& point);
    bool sweep(bsp_t const& bsp, line_t const& line, vec2_t &intersection, line_t &intersected);
    bool clip(clip_context_t &ctx, size_t root);
    void dot_solve(bsp_t const& bsp, vec2_t const& p1, vec2_t &p2);

    bsp_t union_op(bsp_t const& a, bsp_t const& b);
    bsp_t intersect_op(bsp_t const& a, bsp_t const& b);
    bsp_t difference_op(bsp_t const& a, bsp_t const& b);
    bsp_t xor_op(bsp_t const& a, bsp_t const& b);

} // namespace alh::bsp

// todo:
//  - write better split function
//  - have a look at new heuristic
//  - check if implementation can use TCO

#endif