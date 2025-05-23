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

    struct paramline_t {
        line_t line;
        float t1, t2;

        paramline_t() { }
        paramline_t(line_t l) : line(l), t1(0.f), t2(1.f) { }
        paramline_t(line_t l, float t1, float t2) : line(l), t1(t1), t2(t2) { assert(t1 != t2); }

        paramline_t &flip() {
            float tmp = t1;
            t1 = t2;
            t2 = tmp;
            return *this;
        }

        //explicit operator line_t() const {
        //    line_t l;
        //    l.p = line.p + (line.q - line.p) * t1;
        //    l.q = line.p + (line.q - line.p) * t2;
        //    l.userdata = line.userdata;
        //    return l;
        //}

        line_t apply() const {
            line_t l;
            l.p = line.p + (line.q - line.p) * t1;
            l.q = line.p + (line.q - line.p) * t2;
            assert(l.p.x != l.q.x || l.p.y != l.q.y);
            l.userdata = line.userdata;
            return l;
        }
    };

    struct node_t {
        paramline_t plane;
        size_t right, left;
    };

    using bsp_t = typename std::vector<node_t>;

    struct clip_context_t;
    typedef bool (*leaf_callback)(clip_context_t &ctx, float t1, float t2, void *userdata);

    struct clip_context_t {
        const bsp_t *bsp;
        const paramline_t *paramline;
        leaf_callback on_empty;
        leaf_callback on_solid;
        void *userdata; // used for output, etc.
    };

    bsp_t build(std::vector<line_t> lines);
    bsp_t build(std::vector<paramline_t> paramlines);
    
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