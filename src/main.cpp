#include <cstdio>
#include <bitset>

#include "raylib.h"
#include "bsp.hpp"
#include "bsp_stl.hpp"
#include "navmesh.hpp"
#include "dbg_shapes.hpp"

#include "test.stl.h"

#if defined(PLATFORM_WEB)
    #include <emscripten/emscripten.h>
#endif

using namespace alh;

vec2_t player_pos = {240.f, 120.f};

static inline void draw_cross(vec2_t p, uint32_t col) {
    // draw a cute little cross
    p.x = floorf(p.x);
    p.y = floorf(p.y);
    vec2_t p1, p2, p3, p4;
    p1 = p - vec2_t{2.f, 2.f};
    p2 = p + vec2_t{2.f, 2.f};
    p3 = p - vec2_t{2.f, -2.f};
    p4 = p + vec2_t{2.f, -2.f};
    dbg_line(p1.x, p1.y, p2.x, p2.y, col);
    dbg_line(p3.x, p3.y, p4.x, p4.y, col);
}

void draw_tree(bsp::bsp_t &bsp, bsp::id_t nid, Color col) {
    // draw bsp tree by traversing nodes

    //Color colors[] = {LIGHTGRAY, GRAY, DARKGRAY, YELLOW, GOLD, ORANGE, PINK, RED, MAROON, GREEN, LIME, DARKGREEN, SKYBLUE, BLUE, DARKBLUE, PURPLE, VIOLET, DARKPURPLE, BEIGE, BROWN, DARKBROWN};
    //Color rcol = colors[GetRandomValue(0, 20)];

    bsp::bsp_node_t const& node = bsp[nid];
    line_t l = node.plane.apply();

    DrawLineV({l.p.x, l.p.y}, {l.q.x, l.q.y}, col);
    if (!bsp::empty_leaf(node.left) && !bsp::solid_leaf(node.left))
        draw_tree(bsp, node.left, col);
    if (!bsp::empty_leaf(node.right) && !bsp::solid_leaf(node.right))
        draw_tree(bsp, node.right, col);
}

void draw_cell(bsp::bsp_t const& bsp, bsp::id_t nid, vec2_t const& point) {

    if (bsp::empty_leaf(nid) || bsp::solid_leaf(nid)) return;

    bsp::bsp_node_t const& node = bsp[nid];
    line_t l = node.plane.apply();

    vec2_t p2 = l.p + (l.q - l.p) * -100.f;
    vec2_t q2 = l.p + (l.q - l.p) * 100.f;

    DrawLineV({p2.x, p2.y}, {q2.x, q2.y}, DARKGRAY);
    DrawLineV({l.p.x, l.p.y}, {l.q.x, l.q.y}, RED);
    if (point.is_left_of(l)) {
        draw_cell(bsp, node.left, point);
    } else {
        draw_cell(bsp, node.right, point);
    }
}

void draw_navmesh(bsp::navmesh::navmesh_t const& navmesh) {
    auto const& links = navmesh.links;
    for (auto const& n : navmesh.nodes) {
        draw_cross(n.position, 0x8080ff);
        //for (size_t i= n.links_start; i != n.links_end; i++) {
        //    auto const& n2 = navmesh.nodes[links[i].target];
        //    dbg_line(n.position.x, n.position.y, n2.position.x, n2.position.y, 0x404080);
        //}
    }
}

void draw_path(bsp::navmesh::navmesh_t const& navmesh, std::vector<size_t> const& path) {
    vec2_t p1, p2;
    alh::bsp::navmesh::nav_node_t n1, n2;
    
    for (size_t i=0; i<path.size()-1; i++) {
        n1 = navmesh.nodes[path[i]];
        n2 = navmesh.nodes[path[i+1]];
        p1 = n1.position;
        p2 = n2.position;
        dbg_line(p1.x, p1.y, p2.x, p2.y, 0xff0000);
        
        draw_cross(p1, 0xff0000);
        //dbg_line(p1.x, p1.y, n1.face.p.x, n1.face.p.y, 0x00a0ff);
        //dbg_line(p1.x, p1.y, n1.face.q.x, n1.face.q.y, 0xffa000);
    }

    draw_cross(p2, 0xff0000);
    //dbg_line(p2.x, p2.y, n2.face.p.x, n2.face.p.y, 0x00a0ff);
    //dbg_line(p2.x, p2.y, n2.face.q.x, n2.face.q.y, 0xffa000);
}

bsp::bsp_t g_bsp;
bsp::navmesh::navmesh_t g_navmesh;

void update_draw_frame() {

    // update
    vec2_t next_pos = player_pos;
    if (IsKeyDown(KEY_RIGHT) || IsKeyDown(KEY_D)) next_pos.x = player_pos.x + 1.f;
    if (IsKeyDown(KEY_LEFT) || IsKeyDown(KEY_A)) next_pos.x = player_pos.x - 1.f;
    if (IsKeyDown(KEY_DOWN) || IsKeyDown(KEY_S)) next_pos.y = player_pos.y + 1.f;
    if (IsKeyDown(KEY_UP) || IsKeyDown(KEY_W)) next_pos.y = player_pos.y - 1.f;

    bsp::dot_solve(g_bsp, player_pos, next_pos);
    player_pos = next_pos;

    Vector2 mpos = GetMousePosition();
    vec2_t target_pos = {mpos.x, mpos.y};
    
    {
        vec2_t result;
        line_t line;
        if (bsp::sweep(g_bsp, {player_pos, target_pos}, result, line))
            target_pos = result;
    }
    
    // draw
    BeginDrawing();
    
    ClearBackground(BLACK);

    // draw polygons for empty leaves
    auto ids = bsp::navmesh::empty_leaves(g_bsp);
    for (auto id : ids) {
        auto cell = bsp::navmesh::leaf_poly(g_bsp, id);
        draw_tree(cell, 0, DARKGRAY);

        // get center
        vec2_t acc = {0, 0};
        for (bsp::bsp_node_t n : cell)
            acc = acc + n.plane.apply().p;
        acc = acc / cell.size();

        draw_cross(acc, 0x505050);
    }

    SetRandomSeed(0xdeafbeef);
    draw_tree(g_bsp, 0, WHITE);

    bsp::id_t id_mpos = bsp::leaf_id(g_bsp, 0, {mpos.x, mpos.y});
    bsp::id_t id_player = bsp::leaf_id(g_bsp, 0, {player_pos.x, player_pos.y});
    
    if (!bsp::is_solid(g_bsp, 0, {mpos.x, mpos.y})) {
        auto points = bsp::navmesh::find_path(g_bsp,
                                              g_navmesh,
                                              {player_pos.x, player_pos.y},
                                              {mpos.x, mpos.y});
        for (size_t i=0; i<points.size()-1; i++) {
            vec2_t p1 = points[i];
            vec2_t p2 = points[i+1];
            dbg_line(p1.x, p1.y, p2.x, p2.y, 0xff00ff);
        }
    }

    line_t const& rootplane = g_bsp.front().plane.apply();
    DrawLineV({rootplane.p.x, rootplane.p.y}, {rootplane.q.x, rootplane.q.y}, GREEN);

    DrawLineV({player_pos.x, player_pos.y}, {target_pos.x, target_pos.y}, GREEN);
    DrawCircle(player_pos.x, player_pos.y, 3.f, bsp::is_solid(g_bsp, 0, player_pos) ? RED : BLUE);
    DrawCircle(mpos.x, mpos.y, 3.f, bsp::is_solid(g_bsp, 0, {mpos.x, mpos.y}) ? RED : BLUE);
    draw_navmesh(g_navmesh);

    dbg_shapes().draw();

    EndDrawing();
}

int main() {

    // build bsp from .stl
    g_bsp = bsp::from_stl(test_stl, test_stl_len);

//    std::vector<line_t> lines = {
//        {{60, 40}, {340, 40}},
//        {{340, 40}, {320, 220}},
//        {{320, 220}, {60, 240}},
//        {{60, 240}, {60, 40}},
//
//        {{200, 160}, {210, 140}},
//        {{210, 140}, {180, 150}},
//        {{180, 150}, {200, 160}},
//
//        {{220, 190}, {240, 190}},
//        {{240, 190}, {240, 160}},
//        {{240, 160}, {230, 160}},
//        {{230, 160}, {230, 180}},
//        {{230, 180}, {220, 180}},
//        {{220, 180}, {220, 190}},
//
//        {{140, 140}, {120, 120}},
//        {{120, 120}, {100, 140}},
//        {{100, 140}, {120, 160}},
//        {{120, 160}, {140, 140}}, // having two lines on the same hyperplane creates a zero-area leaf
//
//        {{200, 80}, {260, 90}},
//        {{260, 90}, {270, 80}},
//        {{270, 80}, {210, 70}},
//        {{210, 70}, {200, 80}},
//    };
//
//    g_bsp = bsp::build(lines);
    g_navmesh = bsp::navmesh::build(g_bsp);

    InitWindow(400, 300, "BSP test");
#if defined(PLATFORM_WEB)
    emscripten_set_main_loop(update_draw_frame, 0, 1);
#else

    SetTargetFPS(60);
    while (!WindowShouldClose()) {
        update_draw_frame();
    }
#endif

    CloseWindow();

    return 0;
}