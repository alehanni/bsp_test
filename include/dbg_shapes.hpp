#ifndef ALH_DBG_SHAPES_HPP
#define ALH_DBG_SHAPES_HPP

#include <cassert>
#include <cstdint>
#include <cmath>
#include <vector>

#include "raylib.h"
#include "alh.hpp"

using namespace alh;

struct dbg_shapes {

    void queue_linestrip(std::vector<float> in_coords, uint32_t col) {
        assert(0 == in_coords.size() % 2);
        assert(in_coords.size() >= 4);
        // add linestrip to draw queue
        linestrip_t ls;
        ls.first = coords.size();
        coords.insert(coords.cend(), in_coords.begin(), in_coords.end());
        ls.last = coords.size() - 1;
        ls.color = col;
        linestrips.push_back(ls);
    }

    void queue_circle(float x, float y, float r, uint32_t col) {
        circles.push_back({x, y, r, col});
    }

    void draw() {
        for (linestrip_t ls : linestrips) {
            int count = (ls.last - ls.first + 1) / 2;
            DrawLineStrip((Vector2 *)&coords[ls.first], count, hex_to_color(ls.color));
        }

        for (circle_t c : circles) {
            DrawCircleLinesV({c.x, c.y}, c.r, hex_to_color(c.color));
        }

        circles.clear();
        linestrips.clear();
        coords.clear();
    }

private:

    Color hex_to_color(uint32_t in) {
        Color c;
        c.r = (in & 0xff0000) / (256 * 256);
        c.g = (in & 0x00ff00) / 256;
        c.b = (in & 0x0000ff);
        c.a = 0xff;
        return c;
    }

    struct linestrip_t {
        std::size_t first;
        std::size_t last;
        uint32_t color;
    };

    inline static std::vector<float> coords;
    inline static std::vector<linestrip_t> linestrips;
    
    struct circle_t {
        float x, y, r;
        uint32_t color;
    };

    inline static std::vector<circle_t> circles;
};

void dbg_rectangle(float x, float y, float w, float h, uint32_t col) {
    dbg_shapes().queue_linestrip({
        x, y,
        x, y + h,
        x + w, y + h,
        x + w, y,
        x, y
    }, col);
}

void dbg_rectangle(float x, float y, float w, float h, float px, float py, float a, uint32_t col) {
    // draw a rectangle rotated a radians around local pivot point (px, py)
    float nx = cosf(a);
    float ny = sinf(a);
    float s = px / w;
    float t = py / h;
    vec2_t tx0 = {x/2.f - s*w*nx, y/2.f - s*w*ny};
    vec2_t ty0 = {x/2.f + t*h*ny, y/2.f - t*h*nx};
    vec2_t wside = {w*nx, w*ny};
    vec2_t hside = {-h*ny, h*nx};
    dbg_shapes().queue_linestrip({
        tx0.x + ty0.x, tx0.y + ty0.y,
        tx0.x + ty0.x + hside.x, tx0.y + ty0.y + hside.y,
        tx0.x + ty0.x + wside.x + hside.x, tx0.y + ty0.y + wside.y + hside.y,
        tx0.x + ty0.x + wside.x, tx0.y + ty0.y + wside.y,
        tx0.x + ty0.x, tx0.y + ty0.y
    }, col);
}

void dbg_line(float x1, float y1, float x2, float y2, uint32_t col) {
    dbg_shapes().queue_linestrip({x1, y1, x2, y2}, col);
}

void dbg_circle(float x, float y, float r, uint32_t col) {
    dbg_shapes().queue_circle(x, y, r, col);
}

// todo: add re-ordering step based on depth property

#endif