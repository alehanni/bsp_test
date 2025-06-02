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

//struct traj_t {
//    uint64_t mask; // mask corresponds to depth
//    uint64_t bits; // sequence of left/right traversal choices
//    friend bool operator==(traj_t const& lhs, traj_t const& rhs) { return (lhs.mask==rhs.mask && lhs.bits==rhs.bits); }
//    friend bool operator!=(traj_t const& lhs, traj_t const& rhs) { return !(lhs == rhs); }
//};

//void traj_at_point_impl(bsp_t const& bsp, size_t nid, size_t depth, vec2_t const& point, traj_t &tr) {
//    assert(depth < 64);
//
//    if (empty_leaf(nid) || solid_leaf(nid)) return;
//    assert(nid < bsp.size());
//
//    tr.mask = ~(0xffffffffffffffff << (depth + 1));
//
//    bsp_node_t const& node = bsp[nid];
//    if (point.is_left_of(node.plane.apply())) {
//        return traj_at_point_impl(bsp, node.left, depth + 1, point, tr);
//    } else {
//        tr.bits |= 1 << depth;
//        return traj_at_point_impl(bsp, node.right, depth + 1, point, tr);
//    }
//}

//traj_t traj_at_point(bsp_t const& bsp, vec2_t const& point) {
//    traj_t tr = { 0 };
//    traj_at_point_impl(bsp, 0, 0, point, tr);
//    return tr;
//}

//void empty_leaves_impl(bsp_t const& bsp, size_t nid, size_t depth, traj_t tr, std::vector<traj_t> &leaves) {
//    assert(depth < 64);
//
//    if (solid_leaf(nid)) return;
//    if (empty_leaf(nid)) {
//        leaves.push_back(tr);
//        return;
//    }
//    assert(nid < bsp.size());
//
//    tr.mask = ~(0xffffffffffffffff << (depth + 1));
//
//    bsp_node_t const& node = bsp[nid];
//    empty_leaves_impl(bsp, node.left, depth + 1, tr, leaves);
//    tr.bits |= 1 << depth;
//    empty_leaves_impl(bsp, node.right, depth + 1, tr, leaves);
//}
//
//std::vector<traj_t> empty_leaves(bsp_t const& bsp) {
//    traj_t tr = { 0 };
//    std::vector<traj_t> leaves;
//    empty_leaves_impl(bsp, 0, 0, tr, leaves);
//    return leaves;
//}

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

//std::vector<id_t> leaf_polygons(bsp_t const& bsp) {
//    // generate a polygon for each leaf
//    std::vector<id_t> parents = parent_list(bsp);
//
//    unsigned n_leaves = 0;
//    for (bsp_node_t node : bsp) {
//        if (is_leaf(node.left)) n_leaves++;
//        if (is_leaf(node.right)) n_leaves++;
//    }
//
//    bsp_t polygons;
//    std::vector<id_t> poly_roots(n_leaves, 0);
//    
//    for (id_t i=0; i<bsp.size(); i++) {
//        bsp_node_t const& node = bsp[i];
//        if (is_leaf(node.left) || is_leaf(node.right)) {
//
//            // build stack of cutting planes
//            std::vector<paramline_t> stack;
//
//            id_t prev = i;
//            id_t parent = parents[i];
//            while(NULL_ID != parent) {
//                assert(bsp[parent].left == prev || bsp[parent].right == prev);
//
//                paramline_t plane = bsp[parent].plane;
//                plane.t1 = -1000.f;
//                plane.t2 = 1000.f;
//                plane.line.userdata = (void *)(size_t)parent;
//
//                if (bsp[parent].left == prev) {
//                    stack.push_back(plane.flip());
//                } else {
//                    stack.push_back(plane);
//                }
//
//                prev = parent;
//                parent = parents[parent];
//            }
//
//            // build cell polygons
//            paramline_t plane = node.plane;
//            plane.t1 = -1000.f;
//            plane.t2 = 1000.f;
//            plane.line.userdata = (void *)(size_t)i;
//
//            if (is_leaf(node.left)) {
//                bsp_t poly = bsp::build({plane});
//                for (paramline_t plane : stack) {
//                    poly = bsp::difference_op(poly, build({plane}));
//                }
//                
//                id_t id = (node.left & ~IS_LEAF) & ~IS_SOLID;
//                poly_roots[id] = polygons.size();
//                polygons.insert(polygons.end(), poly.begin(), poly.end());
//            }
//
//            if (is_leaf(node.right)) {
//                bsp_t poly = bsp::build({plane.flip()});
//                for (paramline_t plane : stack) {
//                    poly = bsp::difference_op(poly, build({plane}));
//                }
//
//                id_t id = (node.right & ~IS_LEAF) & ~IS_SOLID;
//                poly_roots[id] = polygons.size();
//                polygons.insert(polygons.end(), poly.begin(), poly.end());
//            }
//        }
//    }
//
//    return poly_roots;
//}

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

//bsp_t bsp_poly(bsp_t const& bsp, size_t root, traj_t tr) {
//
//    bsp_t poly;
//    bsp_node_t node = bsp[root];
//    uint64_t mask = 1;
//
//    size_t nid = root;
//
//    while(tr.mask & mask) {
//
//        // todo: replace scaling with bb intersection points
//        paramline_t longplane = node.plane;
//        longplane.t1 = -1000.f;
//        longplane.t2 = 1000.f;
//        assert(!empty_leaf(nid));
//        assert(!solid_leaf(nid));
//        longplane.line.userdata = (void *)nid;
//
//        if (tr.bits & mask) {
//            if (1 == mask)
//                poly = build({longplane.flip()});
//            else
//                poly = bsp::difference_op(poly, build({longplane}));
//            
//            if (empty_leaf(node.right) || solid_leaf(node.right)) break;
//            nid = node.right;
//            node = bsp[node.right];
//        } else {
//            if (1 == mask)
//                poly = build({longplane});
//            else
//                poly = bsp::difference_op(poly, build({longplane.flip()}));
//            
//            if (empty_leaf(node.left) || solid_leaf(node.left)) break;
//            nid = node.left;
//            node = bsp[node.left];
//        }
//
//        mask = mask << 1;
//    }
//
//    return poly;
//}

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

struct portal_t {
    paramline_t paramline;
    id_t from;
    id_t to;
};

std::vector<portal_t> generate_portals(const bsp_t &bsp) {

    std::vector<id_t> ids = empty_leaves(bsp);
    bsp_t cells;
    std::vector<size_t> cell_start;

    // generate leaf cell boundaries
    for (id_t id : ids) {
        bsp_t poly = leaf_poly(bsp, id);
        cell_start.push_back(cells.size());
        cells.insert(cells.end(), poly.begin(), poly.end());
    }
    cell_start.push_back(cells.size());

    // create portal pair for each neighbor pair
    std::vector<portal_t> portals;

    for (size_t i=0; i<ids.size(); i++) {
        for (size_t j=i+1; j<ids.size(); j++) {
            std::vector<bsp_node_t>::iterator a_begin, a_end, b_begin, b_end;
            a_begin = cells.begin() + cell_start[i];
            a_end = cells.begin() + cell_start[i+1];
            b_begin = cells.begin() + cell_start[j];
            b_end = cells.begin() + cell_start[j+1];

            for (auto a_it = a_begin; a_it != a_end; ++a_it) {
                line_t al = a_it->plane.line;
                for (auto b_it = b_begin; b_it != b_end; ++b_it) {
                    line_t bl = b_it->plane.line;
                    if (al.userdata == bl.userdata) {
                        assert(al.p == bl.p && al.q == bl.q);
                        paramline_t apl = a_it->plane;
                        paramline_t bpl = b_it->plane;

                        float a_tmin, a_tmax, b_tmin, b_tmax;
                        a_tmin = std::min(apl.t1, apl.t2);
                        a_tmax = std::max(apl.t1, apl.t2);
                        b_tmin = std::min(bpl.t1, bpl.t2);
                        b_tmax = std::max(bpl.t1, bpl.t2);

                        portal_t portal1, portal2; // great games
                        portal1 = { al, ids[i], ids[j] };
                        portal2 = { al, ids[j], ids[i] };
                        bool b_overlap = false;

                        if ((b_overlap = a_tmin <= b_tmin && b_tmin < a_tmax)) {
                            portal1.paramline.t1 = b_tmin;
                        } else if ((b_overlap = b_tmin <= a_tmin && a_tmin < b_tmax)) {
                            portal1.paramline.t1 = a_tmin;
                        }

                        if (b_overlap) {
                            portal1.paramline.t2 = std::min(a_tmax, b_tmax);
                            portal2.paramline.t1 = portal1.paramline.t2;
                            portal2.paramline.t2 = portal1.paramline.t1;

                            portals.push_back(portal1);
                            for (size_t i=0, j=1; j<portals.size(); i++, j++) {
                                size_t i2 = portals.size() - i - 1;
                                size_t j2 = portals.size() - j - 1;
                                
                                if (portals[j2].from > portals[i2].from) {
                                    portal_t tmp = portals[j2];
                                    portals[j2] = portals[i2];
                                    portals[i2] = tmp;
                                } else {
                                    break;
                                }
                            }
                            
                            portals.push_back(portal2);
                            for (size_t i=0, j=1; j<portals.size(); i++, j++) {
                                size_t i2 = portals.size() - i - 1;
                                size_t j2 = portals.size() - j - 1;
                                
                                if (portals[j2].from > portals[i2].from) {
                                    portal_t tmp = portals[j2];
                                    portals[j2] = portals[i2];
                                    portals[i2] = tmp;
                                } else {
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return portals;
}

// build graph from bsp portals, use a* to route paths

struct nav_node_t {
    line_t face;
    size_t links_start, links_end;
};

struct nav_link_t {
    nav_link_t(id_t tgt) : target(tgt) { }
    size_t target;
    float weight;
};

struct navmesh_t {
    std::vector<nav_node_t> nodes;
    std::vector<nav_link_t> links;
    std::vector<std::pair<size_t, size_t>> lookup;
};

navmesh_t build(bsp_t const& bsp) {
    std::vector<portal_t> portals = generate_portals(bsp);
    
    size_t n_leaves = 0;
    for (bsp_node_t node : bsp) {
        if (is_leaf(node.left)) n_leaves++;
        if (is_leaf(node.right)) n_leaves++;
    }

    // use a vector of vectors to build, concatenate later
    std::vector<std::vector<id_t>> buckets;
    buckets.resize(n_leaves, {});
    for (size_t i=0; i<portals.size(); i++) {
        auto const& p = portals[i];
        buckets[p.from].push_back(i);
    }

    navmesh_t graph;

    // build lookup table
    graph.lookup.resize(n_leaves);
    for (size_t i=0; i<portals.size(); i++) {
        const id_t id = portals[i].from;
        id_t start=0, end=portals.size()-1;
        
        for (; id != portals[start].from; start++);
        for (; id != portals[end].from; end--);
        for (; id == portals[i].from; i++);

        graph.lookup[id].first = start;
        graph.lookup[id].second = end+1;
    }

    // link all the portals/nodes in the same bucket
    for (size_t i=0; i<portals.size(); i++) {
        auto const& portal = portals[i];
        auto const& bucket = buckets[portal.to];

        nav_node_t node;
        node.face = portal.paramline.apply();
        node.links_start = graph.links.size();
        graph.links.insert(graph.links.end(), bucket.begin(), bucket.end());
        node.links_end = graph.links.size();
        graph.nodes.push_back(std::move(node));
    }

    // set weights to euclidean distance
    for (nav_node_t const& n : graph.nodes) {
        // todo: use line-to-line distance instead of midpoint distance??
        vec2_t p1 = n.face.p + (n.face.q - n.face.p) * 0.5f;
        for (size_t i=n.links_start; i<n.links_end; i++) {
            nav_link_t &link = graph.links[i];
            nav_node_t const& n2 = graph.nodes[link.target];
            vec2_t p2 = n2.face.p + (n2.face.q - n2.face.p) * 0.5f;
            link.weight = dist(p1, p2);
        }
    }

    return graph;
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
        do {
            path.push_back(id);
        } while(id != src && (id = prev[id]));
    }
    std::reverse(path.begin(), path.end());

    return path;
}

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

void shortest_path(navmesh_t const& navmesh, vec2_t start, vec2_t goal) {

    // find portals at start/goal trajectory
    // select closest/ do a pseudo-first-step of dijkstra
    // run regular dijkstra (with multiple goals)



}

} // namespace navmesh

#endif