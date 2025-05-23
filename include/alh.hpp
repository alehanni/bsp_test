#ifndef ALH_H
#define ALH_H

#include <cassert>
#include <cmath>

namespace alh {

struct vec2_t {
    float x, y;
    friend vec2_t operator+(vec2_t const& lhs, vec2_t const& rhs) {return {lhs.x+rhs.x, lhs.y+rhs.y}; }
    friend vec2_t operator-(vec2_t const& lhs, vec2_t const& rhs) {return {lhs.x-rhs.x, lhs.y-rhs.y}; }
    friend vec2_t operator*(vec2_t const& lhs, float const& rhs) {return {lhs.x*rhs, lhs.y*rhs}; }
    friend vec2_t operator/(vec2_t const& lhs, float const& rhs) {return {lhs.x/rhs, lhs.y/rhs}; }
    
    bool is_left_of(auto l) const { // todo: avoid using auto here
        vec2_t dl = l.q - l.p;
        return (x - l.p.x)*(-dl.y) + (y - l.p.y)*(dl.x) < 0;
    }

    bool is_exactly_on(auto l) const {
        vec2_t dl = l.q - l.p;
        return (x - l.p.x)*(-dl.y) + (y - l.p.y)*(dl.x) == 0;
    }
};

static float dist2(vec2_t const& p, vec2_t const& q) {
    vec2_t pq = q - p;
    return (pq.x*pq.x + pq.y*pq.y);
}

static vec2_t get_normal(vec2_t p, vec2_t q) {
    vec2_t pq = q - p;
    float dist = sqrtf(pq.x*pq.x + pq.y*pq.y);
    assert(0.f != dist);
    return {-pq.y / dist, pq.x / dist};
}

struct line_t {
    vec2_t p, q;
    void *userdata;
    
    vec2_t normal;
    line_t &w_normal() { normal = get_normal(p, q); return *this; }
};

static bool line_intersect_gg3(line_t const& l1, line_t const& l2, float &alpha, float &beta) {
    vec2_t a, b, c;
    a = l1.q - l1.p;
    b = l2.p - l2.q;
    c = l1.p - l2.p;

    float num_a, denom_a;
    num_a = b.y*c.x - b.x*c.y;
    denom_a = a.y*b.x - a.x*b.y;

    if (denom_a == 0) {
        return false; // collinear segments
    } else {
        alpha = num_a / denom_a;
        beta = (a.x*c.y - a.y*c.x) / denom_a;
        return true;
    }
}

static vec2_t project_to(vec2_t const& p, line_t const& line) {
    float p_a_x, p_a_y, b_a_x, b_a_y;
    p_a_x = p.x - line.p.x;
    p_a_y = p.y - line.p.y;
    b_a_x = line.q.x - line.p.x;
    b_a_y = line.q.y - line.p.y;

    float num = p_a_x*b_a_x + p_a_y*b_a_y;
    float denom = b_a_x*b_a_x + b_a_y*b_a_y;
    float s = num/denom;
    
    return {line.p.x + b_a_x * s, line.p.y + b_a_y * s};
}

}

#endif