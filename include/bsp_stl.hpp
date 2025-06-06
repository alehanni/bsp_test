
#include <cstdint>
#include <cstring>

#include "alh.hpp"
#include "bsp.hpp"

#ifndef ALH_BSP_STL_HPP
#define ALH_BSP_STL_HPP

namespace alh::bsp {

    bsp_t from_stl(const uint8_t *data, size_t len) {
        // build bsp from the union of .stl triangles
        assert(sizeof(float) == 4);
        assert(len > 84);
        size_t off = 80; // skip header

        uint32_t n_tris = 0;
        std::memcpy(&n_tris, data + off, 4); // get number of triangles
        off += 4;
        assert(n_tris > 0);

        // turn triangles into lines
        std::vector<line_t> lines;
        lines.reserve(n_tris * 3);
        for(size_t n=0; n<n_tris && off<len; n++, off+=50) {
            vec2_t v1, v2, v3; // (skip z component)
            memcpy(&v1, data + off + 12, 8);
            memcpy(&v2, data + off + 24, 8);
            memcpy(&v3, data + off + 36, 8);

            line_t l1, l2, l3;
            if (cross(v1 - v3, v2 - v3) > 0.f) {
                l1 = {v1, v2};
                l2 = {v2, v3};
                l3 = {v3, v1};
            } else {
                l1 = {v1, v3};
                l2 = {v3, v2};
                l3 = {v2, v1};
            }

            lines.push_back(l1);
            lines.push_back(l2);
            lines.push_back(l3);
        }

        // remove diagonals
        for (size_t i=0; i<lines.size(); i++) {
            for (size_t j=i+1; j<lines.size(); j++) {
                line_t l1 = lines[i];
                line_t l2 = lines[j];
                if (l1.p == l2.q && l1.q == l2.p) {
                    lines[i] = lines[lines.size() - 1];
                    lines[j] = lines[lines.size() - 2];
                    lines.pop_back();
                    lines.pop_back();
                    i -= 1;
                    break;
                }
            }
        }

        bsp_t bsp = build(lines);
        return bsp;
    }

} // namespace alh::bsp

#endif