#ifndef ALH_BSP_HPP
#define ALH_BSP_HPP

#include <cassert>
#include <cstdint>
#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>

#include "alh.hpp"

namespace alh::bsp {

    using id_t = uint32_t;
    static constexpr id_t NULL_ID = (id_t)(-1);
    static constexpr id_t IS_LEAF = 0x80000000;
    static constexpr id_t IS_SOLID = 0x40000000;
    static inline bool is_leaf(id_t nid) { return (nid & IS_LEAF); };
    static inline bool solid_leaf(id_t nid) { return (is_leaf(nid) && nid & IS_SOLID); }
    static inline bool empty_leaf(id_t nid) { return (is_leaf(nid) && !(nid & IS_SOLID)); }

    struct paramline_t {
        line_t line;
        float t1, t2;

        paramline_t() { }
        paramline_t(line_t l) : line(l), t1(0.f), t2(1.f) { }
        paramline_t(line_t l, float t1, float t2) : line(l), t1(t1), t2(t2) { assert(fabs(t1 - t2) > 1e-4); }

        paramline_t &flip() {
            float tmp = t1;
            t1 = t2;
            t2 = tmp;
            return *this;
        }

        line_t apply() const {
            line_t l;
            l.p = line.p + (line.q - line.p) * t1;
            l.q = line.p + (line.q - line.p) * t2;
            assert(l.p.x != l.q.x || l.p.y != l.q.y);
            l.userdata = line.userdata;
            return l;
        }
    };

    struct bsp_node_t {
        paramline_t plane;
        id_t right, left;
    };

    using bsp_t = typename std::vector<bsp_node_t>;

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
    
    bool is_solid(bsp_t const& bsp, id_t nid, vec2_t const& point);
    id_t leaf_id(bsp_t const& bsp, id_t nid, vec2_t const& point);
    bool sweep(bsp_t const& bsp, line_t const& line, vec2_t &intersection, line_t &intersected);
    bool clip(clip_context_t &ctx, id_t root);
    void dot_solve(bsp_t const& bsp, vec2_t const& p1, vec2_t &p2);

    void serialize(bsp_t const& in);
    bsp_t deserialize(uint8_t *data, size_t len);

    bsp_t union_op(bsp_t const& a, bsp_t const& b);
    bsp_t intersect_op(bsp_t const& a, bsp_t const& b);
    bsp_t difference_op(bsp_t const& a, bsp_t const& b);
    bsp_t xor_op(bsp_t const& a, bsp_t const& b);

} // namespace alh::bsp

#endif