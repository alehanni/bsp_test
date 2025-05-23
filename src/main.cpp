#include <cstdio>
#include <bitset>

#include "raylib.h"
#include "bsp.hpp"
#include "navmesh.hpp"
#include "dbg_shapes.hpp"

using namespace alh;

vec2_t player_pos = {240.f, 120.f};

void draw_tree(bsp::bsp_t &bsp, size_t nid, Color col) {
    // draw bsp tree by traversing nodes

    //Color colors[] = {LIGHTGRAY, GRAY, DARKGRAY, YELLOW, GOLD, ORANGE, PINK, RED, MAROON, GREEN, LIME, DARKGREEN, SKYBLUE, BLUE, DARKBLUE, PURPLE, VIOLET, DARKPURPLE, BEIGE, BROWN, DARKBROWN};
    //Color rcol = colors[GetRandomValue(0, 20)];

    bsp::node_t const& node = bsp[nid];
    line_t l = node.plane.apply();

    DrawLineV({l.p.x, l.p.y}, {l.q.x, l.q.y}, col);
    if (bsp::EMPTY_LEAF != node.left && bsp::SOLID_LEAF != node.left)
        draw_tree(bsp, node.left, col);
    if (bsp::EMPTY_LEAF != node.right && bsp::SOLID_LEAF != node.right)
        draw_tree(bsp, node.right, col);
}

void draw_cell(bsp::bsp_t const& bsp, size_t nid, vec2_t const& point) {

    if (bsp::EMPTY_LEAF == nid || bsp::SOLID_LEAF == nid) return;

    bsp::node_t const& node = bsp[nid];
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

int main() {

    std::vector<line_t> lines = {
        {{60, 40}, {340, 40}},
        {{340, 40}, {320, 220}},
        {{320, 220}, {60, 240}},
        {{60, 240}, {60, 40}},

        {{200, 160}, {210, 140}},
        {{210, 140}, {180, 150}},
        {{180, 150}, {200, 160}},

        {{220, 190}, {240, 190}},
        {{240, 190}, {240, 160}},
        {{240, 160}, {230, 160}},
        {{230, 160}, {230, 180}},
        {{230, 180}, {220, 180}},
        {{220, 180}, {220, 190}},

        {{140, 140}, {120, 120}},
        {{120, 120}, {100, 140}},
        {{100, 140}, {120, 160}},
        {{120, 160}, {140, 140}},

        {{200, 80}, {260, 90}},
        {{260, 90}, {270, 80}},
        {{270, 80}, {210, 70}},
        {{210, 70}, {200, 80}},
    };

//    for (size_t i=0; i<lines.size(); i++) {
//        lines[i].userdata = (void *)(i+1);
//    }

//    std::vector<line_t> lines = {
//        {{200, 150}, {180, 170}},
//        {{180, 170}, {220, 170}},
//        {{220, 170}, {200, 150}}
//    };

//    std::vector<line_t> lines = {
//        {{180, 170}, {220, 170}},
//        {{220, 170}, {220, 130}},
//        {{220, 130}, {180, 130}},
//        {{180, 130}, {180, 170}}
//    };

    auto bsp = bsp::build(lines);

//    std::vector<line_t> lines_a = {
//        {{200 - 10, 150 + 20}, {200 + 20 - 10, 150 + 20}},
//        {{200 - 10 + 20, 150 + 20}, {200 + 20 - 10, 150}},
//        {{200 - 10 + 20, 150}, {200 - 10, 150}},
//        {{200 - 10, 150}, {200 - 10, 150 + 20}}, //
//    };
//    bsp::bsp_t bsp_a = bsp::build(lines_a);
//
//    std::vector<line_t> lines_b = {
//        {{200, 150 + 20 + 5}, {200 + 20, 150 + 20 + 5}},
//        {{200 + 20, 150 + 20 + 5}, {200 + 20, 150 + 5}},
//        {{200 + 20, 150 + 5}, {200, 150 + 5}},
//        {{200, 150 + 5}, {200, 150 + 20 + 5}},
//    };
//    bsp::bsp_t bsp_b = bsp::build(lines_b);
//
//    auto bsp = bsp::union_op(bsp_a, bsp_b);

//    auto bsp_a = bsp::build({{{0, 150}, {400, 150}, (void *)0xdead}});
//    auto bsp_b = bsp::build({{{200, 0}, {200, 300}, (void *)0xbeef}});
//    auto bsp = bsp::difference_op(bsp_a, bsp_b);

    InitWindow(400, 300, "BSP test");
    SetTargetFPS(60);

    Color colors[] = {LIGHTGRAY, GRAY, DARKGRAY, YELLOW, GOLD, ORANGE, PINK, RED, MAROON, GREEN, LIME, DARKGREEN, SKYBLUE, BLUE, DARKBLUE, PURPLE, VIOLET, DARKPURPLE, BEIGE, BROWN, DARKBROWN};

    while(!WindowShouldClose()) {

        // update
        vec2_t next_pos = player_pos;
        if (IsKeyDown(KEY_RIGHT) || IsKeyDown(KEY_D)) next_pos.x = player_pos.x + 1.f;
        if (IsKeyDown(KEY_LEFT) || IsKeyDown(KEY_A)) next_pos.x = player_pos.x - 1.f;
        if (IsKeyDown(KEY_DOWN) || IsKeyDown(KEY_S)) next_pos.y = player_pos.y + 1.f;
        if (IsKeyDown(KEY_UP) || IsKeyDown(KEY_W)) next_pos.y = player_pos.y - 1.f;

        bsp::dot_solve(bsp, player_pos, next_pos);
        player_pos = next_pos;

        Vector2 mpos = GetMousePosition();
        vec2_t target_pos = {mpos.x, mpos.y};
        
        {
            vec2_t result;
            line_t line;
            if (bsp::sweep(bsp, {player_pos, target_pos}, result, line))
                target_pos = result;
        }

        // draw
        ClearBackground(BLACK);

        //SetRandomSeed(0x2);

        BeginDrawing();
        //for (auto it=bsp.begin(); it!=bsp.end(); ++it)
        //    DrawLineV({it->plane.p.x, it->plane.p.y}, {it->plane.q.x, it->plane.q.y}, colors[GetRandomValue(0, 20)]);
        
        SetRandomSeed(0xdeafbeef);
        draw_tree(bsp, 0, GRAY);
        //draw_cell(bsp, 0, {mpos.x, mpos.y});

        //traj_t traj = get_traj(bsp, {mpos.x, mpos.y});
        //std::bitset<64> s_bits(traj.bits);
        //std::bitset<64> s_mask(traj.mask);
        //printf("bits: %s, mask: %s\n", s_bits.to_string().c_str(), s_mask.to_string().c_str());

        //auto tr_mpos = get_traj(bsp, {mpos.x, mpos.y});
        //auto tr_player = get_traj(bsp, player_pos);
        //auto lcn = find_lcn(bsp, 0, tr_mpos, tr_player);
        //DrawLineV({lcn.plane.p.x, lcn.plane.p.y}, {lcn.plane.q.x, lcn.plane.q.y}, SKYBLUE);

//        auto tr_mpos = bsp::navmesh::get_traj(bsp, {mpos.x, mpos.y});
//        auto tr_player = bsp::navmesh::get_traj(bsp, player_pos);
//        bool b_neighbors = bsp::navmesh::neighbors(bsp, tr_mpos, tr_player);
//        printf("neighbor: %d\n", b_neighbors);
//
//        size_t depth=0, nid=0;
//        bsp::navmesh::find_lcn(bsp, 0, tr_mpos, tr_player, depth, nid);
//        printf("DEPTH: %d\n", depth);
        
//        {
//            auto l = bsp[nid].plane;
//            DrawLineV({l.p.x, l.p.y}, {l.q.x, l.q.y}, PINK);
//        }

        auto tr_mpos = bsp::navmesh::get_traj(bsp, {mpos.x, mpos.y});
        auto cell = bsp::navmesh::bsp_poly(bsp, 0, tr_mpos);
        draw_tree(cell, 0, PINK);

        for (bsp::node_t n : cell) {
            printf("%lu ", (size_t)n.plane.line.userdata);
        }
        printf("\n");

        //auto indices = bsp::navmesh::cell_indices(bsp, 0, tr_mpos);

        //std::vector<size_t> indices;
        //for (bsp::node_t node : cell) indices.push_back((size_t)node.plane.userdata);

        //for (size_t lu : indices) {
        //    printf("%lu ", lu);
        //}
        //printf("\n");

//        std::bitset<64> s_bits(tr_mpos.bits);
//        std::bitset<64> s_mask(tr_mpos.mask);
//        printf("bits: %s, mask: %s\n", s_bits.to_string().c_str(), s_mask.to_string().c_str());

        DrawLineV({player_pos.x, player_pos.y}, {target_pos.x, target_pos.y}, GREEN);
        DrawCircle(player_pos.x, player_pos.y, 3.f, bsp::is_solid(bsp, 0, player_pos) ? RED : BLUE);
        DrawCircle(mpos.x, mpos.y, 3.f, bsp::is_solid(bsp, 0, {mpos.x, mpos.y}) ? RED : BLUE);

        //DrawLineV({bsp[0].plane.line.p.x, bsp[0].plane.line.p.y}, {bsp[0].plane.line.q.x, bsp[0].plane.line.q.y}, PURPLE);

        //BeginBlendMode(BLEND_ADDITIVE);
        //dbg_shapes().draw();
        //EndBlendMode();

        EndDrawing();
    }

    CloseWindow();

    return 0;
}