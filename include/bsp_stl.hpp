
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

        // iteratively build bsp
        bsp_t bsp;
        for(size_t n=0; n<n_tris && off<len; n++, off+=50) {
            vec2_t v1, v2, v3; // (skip z component)
            memcpy(&v1, data + off + 12, 8);
            memcpy(&v2, data + off + 24, 8);
            memcpy(&v3, data + off + 36, 8);

            // inflate triangle (to help union op)
            float eps = 1e-6;
            vec2_t avg = (v1 + v2 + v3) / 3.f;
            v1 = avg + (v1 - avg) * (1.f + eps);
            v2 = avg + (v2 - avg) * (1.f + eps);
            v3 = avg + (v3 - avg) * (1.f + eps);

            line_t l1, l2, l3;
            if (cross(v1 - v3, v2 - v3) < 0.f) {
                l1 = {v1, v2};
                l2 = {v2, v3};
                l3 = {v3, v1};
            } else {
                l1 = {v1, v3};
                l2 = {v3, v2};
                l3 = {v2, v1};
            }

            bsp = union_op(bsp, build((std::vector<line_t>){l1, l2, l3}));
        }

        return bsp;
    }

} // namespace alh::bsp

#endif