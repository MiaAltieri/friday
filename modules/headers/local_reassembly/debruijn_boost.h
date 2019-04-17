//
// Created by Kishwar Shafin on 10/19/18.
//
/*
 * UPDATE: THE BOOST LIBRARY IS MUCH FASTER THAN MY PREVIOUS IMPLEMENTATION AND DeepVariant ALREADY USES IT.
 * SO REPLACING EVERYTHING TO THIS: https://github.com/google/deepvariant/blob/master/deepvariant/realigner/debruijn_graph.h
 * AS DEBRUIJN GRAPHS FAIL TO HANDLE LONG READS, WE INTEND TO REPLACE THIS WITH POA.
 * BUT FOR NOW I WILL FOLLOW DeepVariant.
 * DeepVariant LICENSE: third_party/deepvariant.LICENSE
 * */
#ifndef FRIDAY_DEBRUIJN_GRAPH_H
#define FRIDAY_DEBRUIJN_GRAPH_H

#include "../dataio/bam_handler.h"
#include <iostream>
#include <set>
#include <stack>
#include <queue>
#include <algorithm>
#include <string>
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/depth_first_search.hpp"
#include "boost/graph/graph_traits.hpp"
#include "boost/graph/graphviz.hpp"
#include "boost/graph/reverse_graph.hpp"
#include "absl/container/flat_hash_map.h"
#include "absl/strings/string_view.h"

using namespace std;

namespace DeBruijnGraph_options {
    static constexpr int MIN_EDGE_SUPPORT = 2;
    static constexpr int MIN_BASE_QUALITY = 15;
    static constexpr int MIN_MAP_QUALITY = 14;
    static constexpr int MAX_ALLOWED_PATHS = 256;
    static constexpr int MIN_K = 10;
    static constexpr int MAX_K = 100;
    static constexpr int STEP_K = 1;
};
struct VertexInfo {
    string kmer;
};

struct EdgeInfo {
    int weight;   // The # of multiedges this edge represents.
    bool is_ref;  // True iff this edge is reflected by the reference sequence.
};

class DeBruijnGraph {
public:
    using BoostGraph = boost::adjacency_list<
            boost::setS,            // Out edge list type.
            boost::listS,           // Vertex list type.
            boost::bidirectionalS,  // Directed graph.
            VertexInfo,             // Vertex label.
            EdgeInfo>;              // Edge label.

    using VertexIterator = boost::graph_traits<BoostGraph>::vertex_iterator;
    using EdgeIterator = boost::graph_traits<BoostGraph>::edge_iterator;
    using AdjacencyIterator = boost::graph_traits<BoostGraph>::adjacency_iterator;

    using Vertex = boost::graph_traits<BoostGraph>::vertex_descriptor;
    using Edge = boost::graph_traits<BoostGraph>::edge_descriptor;
    using Path = std::vector<Vertex>;

    using RawVertexIndexMap = std::map<Vertex, int>;
    using VertexIndexMap =
    boost::const_associative_property_map<RawVertexIndexMap>;

    // Convenience method for rebuilding a table usable as the vertex_index_t
    // property that algorithms require.
    void RebuildIndexMap();

    // Accessor for the vertex index table.
    VertexIndexMap IndexMap() const;

    // Ensure a vertex with label kmer is present--adding if necessary.
    Vertex EnsureVertex(absl::string_view kmer);

    // Look up the vertex with this kmer label.
    Vertex VertexForKmer(absl::string_view kmer) const;

    // Is this graph cyclic?
    bool HasCycle() const;

    // Private constructor.  Public interface via factory only allows access to
    // acyclic DeBruijn graphs.  Argument `k` is used to construct the graph;
    // filtering settings are taken from options.
    DeBruijnGraph(
            const string& ref,
            const vector<type_read>& reads,
            int k,
            long long read_start,
            long long read_end);

    // Add edge, implicitly adding the vertices if needed.  If such an edge is
    // already present, we merely increment its weight to reflect its "multiedge"
    // degree.
    Edge AddEdge(Vertex from_vertex, Vertex to_vertex, bool is_ref);

    // Adds kmers from bases starting at start and stopping at end. We add a kmer
    // at each i from start to end (inclusive), and edges between all sequential
    // kmers. Since the first kmer spans k bases starting at start, start + k must
    // be <= bases.size(). Since the last kmer we add starts at end and is k bases
    // long, end + k <= bases.size() as well. Note that this function tolerates
    // end < 0, which causes the code to return immediately.
    void AddKmersAndEdges(absl::string_view bases, int start_index, int end_index,
                          bool is_ref);

    // Add all the edges implied by the given reference string.
    void AddEdgesForReference(absl::string_view ref);

    // Add all the edges implied by the given read (and according to our edge
    // filtering criteria).
    void AddEdgesForRead(const type_read& read);

    // Returns candidate haplotype paths through the graph.  If more that
    // options.max_num_paths paths are found, this will return an empty vector.
    std::vector<Path> CandidatePaths() const;

    // Returns the string traced by a path through the graph.
    string HaplotypeForPath(const Path& path) const;

    // Removes low weight non-ref edges from the graph.
    void Prune();

public:
    // Gets all the candidate haplotypes defined by paths through the graph.  If
    // more than options.max_num_paths() haplotypes are identified, returns an
    // empty vector, to preempt excessive computation.
    vector<string> CandidateHaplotypes() const;

    // Gets a GraphViz representation of the graph.
    string GraphViz() const;

    // Gets the kmer size used in this graph.
    int KmerSize() const { return k_; }

    long long region_start;
    long long region_end;
    BoostGraph g_;
    int k_;
    Vertex source_;
    Vertex sink_;

    // N.B.: kmer strings are owned by VertexInfo objects;
    // map keys are merely pointers.
    absl::flat_hash_map<absl::string_view, Vertex> kmer_to_vertex_;
    RawVertexIndexMap vertex_index_map_;
};


class DeBruijnGraphHelper {
public:
    // We attempt to build acyclic graphs with increasing kmer size until we
    // achieve an acyclic graph---kmer size starts with options.min_k and goes
    // linearly up to options.max_k, stepping by options.step_k.  If we are able
    // to construct an acyclic DeBruijn graph in this manner, it is returned;
    vector<string> Build(
            const string& ref,
            const vector<type_read>& reads,
            long long region_start,
            long long region_end);

};


class CycleDetector : public boost::dfs_visitor<> {
public:
    explicit CycleDetector(bool* has_cycle) : has_cycle(has_cycle) {}

    template <class Edge, class Graph>
    void back_edge(Edge, const Graph&) {
        *has_cycle = true;
    }

private:
    bool* has_cycle;
};

using Vertex = DeBruijnGraph::Vertex;
using VertexIndexMap = DeBruijnGraph::VertexIndexMap;
using Edge = DeBruijnGraph::Edge;
using Path = DeBruijnGraph::Path;

using Read = type_read;

template <class BoostGraph>
class EdgeLabelWriter {
public:
    explicit EdgeLabelWriter(const BoostGraph& g) : g_(g) {}

    void operator()(std::ostream& out, const Edge e) const {
        EdgeInfo ei = g_[e];
        out << "[label=" << ei.weight << (ei.is_ref ? " color=red" : "") << "]";
    }

private:
    const BoostGraph& g_;
};

class ReachableVertexVisitor : public boost::dfs_visitor<> {
public:
    explicit ReachableVertexVisitor(std::set<Vertex>* reachable_vertices)
            : reachable_vertices(reachable_vertices) {}

    template <class Edge, class Graph>
    void tree_edge(Edge e, const Graph& g) {
        Vertex from = boost::source(e, g);
        if (reachable_vertices->find(from) != reachable_vertices->end()) {
            Vertex to = boost::target(e, g);
            reachable_vertices->insert(to);
        }
    }

private:
    std::set<Vertex>* reachable_vertices;
};

template <class BoostGraphT, class VertexIndexMapT>
std::set<Vertex> VerticesReachableFrom(
        Vertex v, const BoostGraphT& g, const VertexIndexMapT& vertex_index_map) {
    std::set<Vertex> reachable_vertices{v};
    ReachableVertexVisitor vis(&reachable_vertices);
    boost::depth_first_search(
            g, boost::visitor(vis).root_vertex(v).vertex_index_map(vertex_index_map));
    return reachable_vertices;
}

#endif //FRIDAY_DEBRUIJN_GRAPH_H
