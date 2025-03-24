
#include <cstdio>
#include "raylib.h"
#include "bsp.hpp"

vec2_t player_pos = {240.f, 150.f};

void draw_tree(bsp_t &bsp, size_t nid, Color col) {
    // draw bsp tree by traversing nodes
    node_t &node = bsp[nid];
    DrawLineV({node.plane.p.x, node.plane.p.y}, {node.plane.q.x, node.plane.q.y}, col);
    if (EMPTY_LEAF != node.left && SOLID_LEAF != node.left)
        draw_tree(bsp, node.left, col);
    if (EMPTY_LEAF != node.right && SOLID_LEAF != node.right)
        draw_tree(bsp, node.right, col);
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
    std::vector<node_t> bsp = build(lines);

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

        dot_solve(bsp, player_pos, next_pos);
        player_pos = next_pos;

        Vector2 mpos = GetMousePosition();
        vec2_t target_pos = {mpos.x, mpos.y};
        
        {
            vec2_t result;
            line_t line;
            if (sweep(bsp, {player_pos, target_pos}, result, line))
                target_pos = result;
        }

        // draw
        ClearBackground(BLACK);

        SetRandomSeed(0x2);

        BeginDrawing();
        for (auto it=bsp.begin(); it!=bsp.end(); ++it)
            DrawLineV({it->plane.p.x, it->plane.p.y}, {it->plane.q.x, it->plane.q.y}, colors[GetRandomValue(0, 20)]);
        
        //draw_tree(bsp, 0, BLUE);

        DrawLineV({player_pos.x, player_pos.y}, {target_pos.x, target_pos.y}, GREEN);
        DrawCircle(player_pos.x, player_pos.y, 3.f, is_solid(bsp, 0, player_pos) ? RED : BLUE);
        DrawCircle(mpos.x, mpos.y, 3.f, is_solid(bsp, 0, {mpos.x, mpos.y}) ? RED : BLUE);

        EndDrawing();
    }

    CloseWindow();

    return 0;
}