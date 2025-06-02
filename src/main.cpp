#include <cstdio>
#include <bitset>

#include "raylib.h"
#include "bsp.hpp"
#include "navmesh.hpp"
#include "dbg_shapes.hpp"

using namespace alh;

vec2_t player_pos = {240.f, 120.f};

void draw_cross(vec2_t p, uint32_t col) {
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
        //{{120, 160}, {140, 140}}, // having two lines on the same hyperplane creates a zero-area leaf

        {{200, 80}, {260, 90}},
        {{260, 90}, {270, 80}},
        {{270, 80}, {210, 70}},
        {{210, 70}, {200, 80}},
    };

    auto bsp = bsp::build(lines);
    auto navmesh = bsp::navmesh::build(bsp);

    //auto path = bsp::navmesh::dijkstra(navmesh, 0, 32);
    //auto points = bsp::navmesh::funnel(navmesh, path, {315.f, 64.f}, {76.f, 136.f});

    auto test = bsp::navmesh::dijkstra(navmesh, 10, 1);    

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

        //printf("%f, %f\n", mpos.x, mpos.y);

        //bsp::navmesh::dijkstra(navmesh, 1, 4);

        // draw
        ClearBackground(BLACK);

        //SetRandomSeed(0x2);

        BeginDrawing();
        //for (auto it=bsp.begin(); it!=bsp.end(); ++it)
        //    DrawLineV({it->plane.p.x, it->plane.p.y}, {it->plane.q.x, it->plane.q.y}, colors[GetRandomValue(0, 20)]);
        
        // draw polygons for empty leaves
//        auto trs = bsp::navmesh::empty_leaves(bsp);
//        for (auto tr : trs) {
//            auto cell = bsp::navmesh::bsp_poly(bsp, 0, tr);
//            draw_tree(cell, 0, DARKGRAY);
//
//            // get center
//            vec2_t acc = {0, 0};
//            for (bsp::bsp_node_t n : cell)
//                acc = acc + n.plane.apply().p;
//            acc = acc / cell.size();
//
//            draw_cross(acc, 0x505050);
//        }

        SetRandomSeed(0xdeafbeef);
        draw_tree(bsp, 0, WHITE);

        //auto traj = bsp::navmesh::traj_at_point(bsp, {mpos.x, mpos.y});
        //printf("{mask = %lu, bits = %lu}\n", traj.mask, traj.bits);

        //auto poly = bsp::navmesh::bsp_poly(bsp, 0, {127, 71});
        //draw_tree(poly, 0, PINK);

        //std::bitset<64> s_bits(traj.bits);
        //std::bitset<64> s_mask(traj.mask);
        //printf("bits: %s, mask: %s\n", s_bits.to_string().c_str(), s_mask.to_string().c_str());

        //auto tr_mpos = get_traj(bsp, {mpos.x, mpos.y});
        //auto tr_player = get_traj(bsp, player_pos);
        //auto lcn = find_lcn(bsp, 0, tr_mpos, tr_player);
        //DrawLineV({lcn.plane.p.x, lcn.plane.p.y}, {lcn.plane.q.x, lcn.plane.q.y}, SKYBLUE);

//        auto tr_mpos = bsp::navmesh::traj_at_point(bsp, {mpos.x, mpos.y});
//        auto cell = bsp::navmesh::bsp_poly(bsp, 0, tr_mpos);
//        draw_tree(cell, 0, PINK);

//        printf("old: ");
//        for (bsp::bsp_node_t n : cell) {
//            printf("%lu ", (size_t)n.plane.line.userdata);
//        }
//        printf("\n");

        bsp::id_t id_mpos = bsp::leaf_id(bsp, 0, {mpos.x, mpos.y});

        auto cell2 = bsp::navmesh::leaf_poly(bsp, id_mpos);
        draw_tree(cell2, 0, PINK);
        //printf("size: %lu\n", cell2.size());

        bsp::id_t id_player = bsp::leaf_id(bsp, 0, {player_pos.x, player_pos.y});

        
        if (!bsp::is_solid(bsp, 0, {mpos.x, mpos.y})) {
            auto path = bsp::navmesh::dijkstra(navmesh, id_player, id_mpos);
            
            printf("path: ");
            for (size_t i : path)
                printf("%lu ", i);
            printf("\n");

            draw_path(navmesh, path);
        }

//        printf("new: ");
//        for (bsp::bsp_node_t n : cell2) {
//            printf("%lu ", (size_t)n.plane.line.userdata);
//        }
//        printf("\n");
//
//        printf("id_mpos: %d\n", id_mpos);

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

        draw_navmesh(navmesh);
        //draw_path(navmesh, path);

        //for (size_t i=0; i<points.size()-1; i++) {
        //    vec2_t p1 = points[i];
        //    vec2_t p2 = points[i+1];
        //    dbg_line(p1.x, p1.y, p2.x, p2.y, 0xff00ff);
        //}

//        for (size_t i=navmesh.lookup[id_mpos].first; i<navmesh.lookup[id_mpos].second; i++) {
//            bsp::navmesh::nav_node_t const& n = navmesh.nodes[i];
//            vec2_t p = n.face.p + (n.face.q - n.face.p) * 0.5f;
//            draw_cross(p, 0xff0000);
//        }

        dbg_shapes().draw();

//        auto portals = bsp::navmesh::generate_portals(bsp);
//        for (auto portal : portals) {
//            line_t l = portal.paramline.apply();
//            DrawLineV({l.p.x, l.p.y}, {l.q.x, l.q.y}, ORANGE);
//        }

        EndDrawing();
    }

    CloseWindow();

    return 0;
}