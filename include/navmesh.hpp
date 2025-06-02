#ifndef ALH_NAVMESH_HPP
#define ALH_NAVMESH_HPP

#include <cassert>
#include <vector>
#include <numeric>
#include <algorithm>
#include <array>
#include <map>

#include "alh.hpp"
//#include "hashmap.hpp"
#include "bsp.hpp"

namespace alh::bsp::navmesh {

std::vector<id_t> parent_list(bsp_t const& bsp) {
    // identify index of parent for each index
    std::vector<id_t> parents(bsp.size(), NULL_ID);
    for (id_t i=0; i<bsp.size(); i++) {
        bsp_node_t const& node = bsp[i];
        if (!is_leaf(node.left)) parents[node.left] = i;
        if (!is_leaf(node.right)) parents[node.right] = i;
    }
    return parents;
}

std::vector<id_t> empty_leaves(bsp_t const& bsp) {
    // scan through nodes
    std::vector<id_t> ids;
    for (bsp_node_t node : bsp) {
        if (empty_leaf(node.left)) ids.push_back((node.left & ~IS_LEAF) & ~IS_SOLID);
        if (empty_leaf(node.right)) ids.push_back((node.right & ~IS_LEAF) & ~IS_SOLID);
    }
    return ids;
}

void bsp_bb(bsp_t const& bsp, vec2_t &min, vec2_t &max) {
    static constexpr float inf = std::numeric_limits<float>::infinity();
    
    min = {inf, inf};
    max = {-inf, -inf};
    
    for (bsp_node_t node : bsp) {
        line_t l = node.plane.apply();
        min.x = std::min(std::min(min.x, l.p.x), l.q.x);
        min.y = std::min(std::min(min.y, l.p.y), l.q.y);
        max.x = std::max(std::max(max.x, l.p.x), l.q.x);
        max.y = std::max(std::max(max.y, l.p.y), l.q.y);
    }
}

// get polygon for a specific leaf id
bsp_t leaf_poly(bsp_t const& bsp, id_t leaf_id) {

    // find parent node with leaf id
    bool is_left = false;
    id_t parent = 0;
    for (; parent<bsp.size(); parent++) {
        bsp_node_t const& n = bsp[parent];
        if ((is_left |= (n.left & ~IS_SOLID) == (leaf_id | IS_LEAF))) goto found;
        if ((n.right & ~IS_SOLID) == (leaf_id | IS_LEAF)) goto found;
    }
    assert(false); // leaf id not in tree!
found:

    // create first plane
    paramline_t plane = bsp[parent].plane;
    plane.t1 = -1000.f;
    plane.t2 = 1000.f;
    plane.line.userdata = (void *)(size_t) parent;
    if (!is_left) plane.flip();
    bsp_t poly = bsp::build({plane});

    // clip against parent planes
    id_t prev = parent;
    std::vector<id_t> parents = parent_list(bsp);
    while(NULL_ID != parents[parent]) {
        parent = parents[parent];

        paramline_t plane = bsp[parent].plane;
        plane.t1 = -1000.f;
        plane.t2 = 1000.f;
        plane.line.userdata = (void *)(size_t) parent;
        assert(bsp[parent].left == prev || bsp[parent].right == prev);
        if (bsp[parent].left == prev) plane.flip();
        poly = difference_op(poly, bsp::build({plane}));

        prev = parent;
    }

    return poly;
}

struct nav_node_t {
    vec2_t position;
    size_t links_start, links_end;
};

struct nav_link_t {
    size_t target;
    line_t portal;
    float weight;
};

struct navmesh_t {
    std::vector<nav_node_t> nodes;
    std::vector<nav_link_t> links;
};

// todo: build navmesh with nodes at leaf centerpoints instead
navmesh_t build(bsp_t const& bsp) {

    size_t n_leaves = 0;
    for (bsp_node_t node : bsp) {
        if (is_leaf(node.left)) n_leaves++;
        if (is_leaf(node.right)) n_leaves++;
    }

    bsp_t cells;
    std::vector<id_t> empty_ids = empty_leaves(bsp);
    std::vector<size_t> cell_start;
    cell_start.reserve(empty_ids.size() + 1);

    // generate polygons for empty leaves
    for (id_t id : empty_ids) {
        bsp_t poly = leaf_poly(bsp, id);
        cell_start.push_back(cells.size());
        cells.insert(cells.end(), poly.begin(), poly.end());
    }
    cell_start.push_back(cells.size());

    // get centerpoints
    std::vector<nav_node_t> nodes(n_leaves, {{0.f, 0.f}, 0, 0});
    for (size_t i=0; i<empty_ids.size(); i++) {
        vec2_t acc = {0.f, 0.f};
        for (size_t j=cell_start[i]; j<cell_start[i+1]; j++) {
            bsp_node_t const& n = cells[j];
            acc = acc + n.plane.apply().p;
        }
        acc = acc / (cell_start[i+1] - cell_start[i]);
        nodes[empty_ids[i]].position = std::move(acc);
    }

    // identify links and create portals
    std::vector<std::pair<id_t, paramline_t>> tagged_planes;
    tagged_planes.reserve(cells.size());
    for (id_t i=0; i<empty_ids.size(); i++) {
        for (size_t j=cell_start[i]; j<cell_start[i+1]; j++)
            tagged_planes.push_back(std::make_pair(empty_ids[i], cells[j].plane));
    }

    std::vector<std::pair<id_t, nav_link_t>> tagged_links;
    for (id_t i=0; i<tagged_planes.size(); i++) {
        for (id_t j=i+1; j<tagged_planes.size(); j++) {
            
            paramline_t a, b;
            a = tagged_planes[i].second;
            b = tagged_planes[j].second;

            if (a.line.userdata == b.line.userdata) {
                // same hyperplane => neighbors if overlap
                id_t a_lid, b_lid;
                a_lid = tagged_planes[i].first;
                b_lid = tagged_planes[j].first;

                float a_tmin, a_tmax, b_tmin, b_tmax;
                a_tmin = std::min(a.t1, a.t2);
                a_tmax = std::max(a.t1, a.t2);
                b_tmin = std::min(b.t1, b.t2);
                b_tmax = std::max(b.t1, b.t2);

                paramline_t portal;

                bool b_overlap = false;
                if ((b_overlap |= a_tmin <= b_tmin && b_tmin < a_tmax)) {
                    portal.t1 = b_tmin;
                } else if ((b_overlap |= b_tmin <= a_tmin && a_tmin < b_tmax)) {
                    portal.t1 = a_tmin;
                }

                if (b_overlap) {
                    portal.t2 = std::min(a_tmax, b_tmax);
                    tagged_links.push_back({a_lid, {b_lid, portal.apply(), 0.f}});
                    tagged_links.push_back({b_lid, {a_lid, portal.flip().apply(), 0.f}});
                }
            }
        }
    }

    // sort links using tags
    std::sort(tagged_links.begin(), tagged_links.end(), [](auto a, auto b){
        return (a.first < b.first);
    });

    // point node ranges to sorted links
    for (id_t id : empty_ids) {
        auto cond = [id](auto a){ return (a.first == id); };
        auto start = std::find_if(tagged_links.begin(), tagged_links.end(), cond);
        nodes[id].links_start = std::distance(tagged_links.begin(), start);
        auto end = std::find_if_not(start, tagged_links.end(), cond);
        nodes[id].links_end = std::distance(tagged_links.begin(), end);
    }

    // discard tags
    std::vector<nav_link_t> links;
    links.reserve(tagged_links.size());
    for (auto tagged_link : tagged_links) links.push_back(tagged_link.second);

    // set weights to euclidean distance
    for (nav_node_t const& n : nodes) {
        for (size_t i=n.links_start; i<n.links_end; i++) {
            nav_link_t &link = links[i];
            nav_node_t const& n2 = nodes[link.target];
            link.weight = dist(n.position, n2.position);
        }
    }

    // create navmesh
    navmesh_t navmesh;
    navmesh.nodes = std::move(nodes);
    navmesh.links = std::move(links);
    return navmesh;
}

using path_t = typename std::vector<size_t>;

path_t dijkstra(navmesh_t const& navmesh, size_t src, size_t dest) {

    // use sorting to create priority queue
    size_t len = navmesh.nodes.size();
    assert(src < len);
    assert(dest < len);

    std::vector<uint16_t> dists(len, 0xffff);
    std::vector<uint16_t> prev(len, 0xffff);
    std::vector<uint32_t> queue;
    queue.reserve(len);
    queue.push_back(0xffff0000 | src); // high bytes are bitflipped distance,
                                       // low bytes are node offset
    bool b_found = false;

    while (!queue.empty()) {
        std::sort(queue.begin(), queue.end());
        uint16_t ui = queue.back() & 0x0000ffff;
        
        assert(0xffff != ui);
        if ((b_found |= dest == ui))
            break;
        
        uint16_t ud = ~(queue.back() & 0xffff0000) >> 16;
        queue.pop_back();
        
        nav_node_t const& n = navmesh.nodes[ui];
        for (size_t i=n.links_start; i<n.links_end; i++) {
            nav_link_t const& link = navmesh.links[i];
            uint16_t vi = link.target;
            assert(ui != vi);

            uint16_t alt = ud + link.weight;

            if (alt < dists[vi]) {
                dists[vi] = alt;
                prev[vi] = ui;
                queue.push_back((~alt << 16) | vi);
            }
        }
    }

    // unroll prev
    std::vector<size_t> path;
    if (b_found) {
        size_t id = dest;
        for(;;) {
            path.push_back(id);
            if (id == src) break;
            id = prev[id];
        }
    }
    std::reverse(path.begin(), path.end());

    return path;
}

/*
std::vector<vec2_t> funnel(navmesh_t const& navmesh, path_t const& path, vec2_t start, vec2_t goal) {

    assert(path.size() > 0);

    vec2_t apex = start;
    vec2_t p = navmesh.nodes[path[0]].face.p; // right funnel point
    vec2_t q = navmesh.nodes[path[0]].face.q; // left funnel point
    std::vector<vec2_t> out;

    for(size_t i=1; i<path.size(); i++) {
        
        line_t portal = navmesh.nodes[path[i]].face;

        if (portal.p.is_left_of((line_t){apex, q})) {
            out.push_back(apex);
            apex = q;
            p = portal.p;
            q = portal.q;
            continue;
        }

        if (portal.p.is_left_of((line_t){apex, p}))
            p = portal.p;
        
        if (!portal.q.is_left_of((line_t){apex, p})) {
            out.push_back(apex);
            apex = p;
            p = portal.p;
            q = portal.q;
            continue;
        }

        if (!portal.q.is_left_of((line_t){apex, q}))
            q = portal.q;
    }
    
    out.push_back(apex);
    out.push_back(goal);

    return out;
}
*/

void find_path(navmesh_t const& navmesh, vec2_t start, vec2_t goal) {

    // find portals at start/goal trajectory
    // select closest/ do a pseudo-first-step of dijkstra
    // run regular dijkstra (with multiple goals)



}

} // namespace navmesh

#endif