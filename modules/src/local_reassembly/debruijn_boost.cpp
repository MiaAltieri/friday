//
// Created by Kishwar Shafin on 10/19/18.
//

#include "../../headers/local_reassembly/debruijn_boost.h"
using Vertex = DeBruijnGraph::Vertex;
using VertexIndexMap = DeBruijnGraph::VertexIndexMap;
using Edge = DeBruijnGraph::Edge;
using Path = DeBruijnGraph::Path;

using Read = type_read;

Vertex DeBruijnGraph::EnsureVertex(absl::string_view kmer) {
    Vertex v;
    auto vertex_find = kmer_to_vertex_.find(kmer);
    if (vertex_find != kmer_to_vertex_.end()) {
        v = (*vertex_find).second;
    } else {
        string kmer_copy(kmer);
        v = boost::add_vertex(VertexInfo{kmer_copy}, g_);
        // N.B.: must use the long-lived string in the map key as the referent of
        // the absl::string_view key.
        kmer_to_vertex_[absl::string_view(g_[v].kmer)] = v;
    }

    return v;
}

Vertex DeBruijnGraph::VertexForKmer(absl::string_view kmer) const {
    return kmer_to_vertex_.at(kmer);
}

void DeBruijnGraph::RebuildIndexMap() {
    std::map<Vertex, int> table;
    VertexIterator vi, vend;
    std::tie(vi, vend) = boost::vertices(g_);
    int index = 0;
    for (; vi != vend; ++vi) {
        table[*vi] = index;
        ++index;
    }
    vertex_index_map_ = table;
}

VertexIndexMap DeBruijnGraph::IndexMap() const {
    boost::const_associative_property_map<RawVertexIndexMap> vmap(
            vertex_index_map_);
    return vmap;
}

bool DeBruijnGraph::HasCycle() const {
    bool has_cycle = false;
    CycleDetector cycle_detector(&has_cycle);
    boost::depth_first_search(
            g_, boost::visitor(cycle_detector).vertex_index_map(IndexMap()));
    return has_cycle;
}

DeBruijnGraph::DeBruijnGraph(
        const string& ref,
        const vector<type_read>& reads,
        int k, long long _start, long long _end) {

    region_start = _start;
    region_end = _end;
    k_ = k;

    if( k <=0 || k >= ref.size()) {
        cerr<<"DeBruijnGraph error, invalid K-mer size: "<<k<<" Reference size: "<<ref.size()<<endl;
    }
    AddEdgesForReference(ref);
    source_ = VertexForKmer(ref.substr(0, k_));
    sink_ = VertexForKmer(ref.substr(ref.size() - k_, k_));

    for(auto &read: reads) {
        if(read.mapping_quality < DeBruijnGraph_options::MIN_MAP_QUALITY) continue;
        if(read.pos_end < region_start) continue;
        if(read.pos > region_end) continue;
        AddEdgesForRead(read);
    }
    RebuildIndexMap();
}


// Indicates that we couldn't find a minimum k that can be used.
constexpr int kBoundsNoWorkingK = -1;
struct KBounds {
    int min_k;  // Mininum k to consider (inclusive).
    int max_k;  // Maximum k to consider (inclusive).
};


KBounds KMinMaxFromReference(const absl::string_view ref) {
    KBounds bounds;
    bounds.min_k = kBoundsNoWorkingK;
    bounds.max_k  = std::min(DeBruijnGraph_options::MAX_K, static_cast<int>(ref.size()) - 1);

    for (int k = DeBruijnGraph_options::MIN_K; k <= bounds.max_k; k += DeBruijnGraph_options::STEP_K) {
        bool has_cycle = false;
        std::set<absl::string_view> kmers;

        for (int i = 0; i < ref.size() - k + 1; i++) {
            absl::string_view kmer = ref.substr(i, k);
            if (kmers.insert(kmer).second == false) {
                // No insertion took place because the kmer already exists. This implies
                // that there's a cycle in the graph.
                has_cycle = true;
                break;
            }
        }

        if (!has_cycle) {
            bounds.min_k = k;
            break;
        }
    }

    return bounds;
}

vector<string> DeBruijnGraphHelper::Build(
        const string& ref,
        const vector<type_read>& reads,
        long long region_start,
        long long region_end) {
    KBounds bounds = KMinMaxFromReference(ref);

    if (bounds.min_k == kBoundsNoWorkingK) return {};

    for (int k = bounds.min_k; k <= bounds.max_k; k += DeBruijnGraph_options::STEP_K) {

        std::unique_ptr<DeBruijnGraph> graph = std::unique_ptr<DeBruijnGraph>(
                new DeBruijnGraph(ref, reads, k, region_start, region_end));
        if (graph->HasCycle()) {
            continue;
        } else {
            graph->Prune();
            return graph->CandidateHaplotypes();
        }
    }
    return {};
}

Edge DeBruijnGraph::AddEdge(Vertex from_vertex, Vertex to_vertex, bool is_ref) {
    bool was_present;
    Edge edge;
    std::tie(edge, was_present) = boost::edge(from_vertex, to_vertex, g_);
    if (!was_present) {
        std::tie(edge, std::ignore) = boost::add_edge(from_vertex, to_vertex,
                                                      EdgeInfo{0, false}, g_);
    }
    EdgeInfo& ei = g_[edge];
    ei.weight++;
    ei.is_ref |= is_ref;
    return edge;
}

void DeBruijnGraph::AddKmersAndEdges(absl::string_view bases, int start_index, int end_index,
                                     bool is_ref) {
    if(start_index >= 0 && start_index + k_ <= bases.size() && end_index + k_ <= bases.size()) {
        // End can be less than 0, in which case we return without doing any work.
        if (end_index > 0) {
            Vertex vertex_prev = EnsureVertex(bases.substr(start_index, k_));
            for (int i = start_index + 1; i <= end_index; ++i) {
                Vertex vertex_cur = EnsureVertex(bases.substr(i, k_));
                AddEdge(vertex_prev, vertex_cur, is_ref);
                vertex_prev = vertex_cur;
            }
        }
    }
}

void DeBruijnGraph::AddEdgesForReference(absl::string_view ref) {
    AddKmersAndEdges(ref, 0, ref.size() - k_, true /* is_ref */);
}


void DeBruijnGraph::AddEdgesForRead(const type_read& read) {
    const string bases = read.sequence;
    const absl::string_view bases_view(bases);
    // Note that this SIGNED int type declaration is key to avoid
    // bases.size() - k_ underflowing.
    if(bases.size() <= k_) {
        return;
    }

    const int stop = bases.size() - k_;
    int i = 0;
    int current_index = 0;
    while (i < stop) {
        int next_bad_position = read.bad_indicies[current_index];
        AddKmersAndEdges(bases_view, i, next_bad_position - k_, false /* is_ref */);
        i = next_bad_position + 1;
        current_index += 1;
    }
}

std::vector<Path> DeBruijnGraph::CandidatePaths() const {
    vector<Path> terminated_paths;
    queue<Path> extendable_paths;

    if(boost::out_degree(source_, g_) <= 0) {
        cerr<<"WARNING: DeBruijnGraph path out degree is 0"<<endl;
        return terminated_paths;
    }
    extendable_paths.push({source_});

    // Inefficient.
    while (!extendable_paths.empty()) {
        // Some windows can have an extremely branchy graph.  Ideally windows would
        // be chosen to avoid this.  We give up if we encounter too many paths.
        int n_total_paths = terminated_paths.size() + extendable_paths.size();
        if (n_total_paths > DeBruijnGraph_options::MAX_ALLOWED_PATHS) {
            return {};
        }

        Path path = extendable_paths.front();
        extendable_paths.pop();
        Vertex last_v = path.back();
        // For each successor of last_v, add path::successor to the
        // appropriate queue.
        AdjacencyIterator vi, vend;
        std::tie(vi, vend) = boost::adjacent_vertices(last_v, g_);
        for (; vi != vend; ++vi) {
            Path extended_path(path);
            extended_path.push_back(*vi);
            if (*vi == sink_ || boost::out_degree(*vi, g_) == 0) {
                terminated_paths.push_back(extended_path);
            } else {
                extendable_paths.push(extended_path);
            }
        }
    }
    return terminated_paths;
}

string DeBruijnGraph::HaplotypeForPath(const Path& path) const {
    std::stringstream haplotype;
    for (Vertex v : path) {
        haplotype << g_[v].kmer[0];
    }
    if (!path.empty()) {
        haplotype << g_[path.back()].kmer.substr(1, k_ - 1);
    }
    return haplotype.str();
}

std::vector<string> DeBruijnGraph::CandidateHaplotypes() const {
    std::vector<string> haplotypes;
    for (const Path& path : CandidatePaths()) {
        haplotypes.push_back(HaplotypeForPath(path));
    }
    std::sort(haplotypes.begin(), haplotypes.end());
    return haplotypes;
}

string DeBruijnGraph::GraphViz() const {
    std::stringstream graphviz;
    auto vertex_label_writer = boost::make_label_writer(
            boost::get(&VertexInfo::kmer, g_));
    boost::write_graphviz(
            graphviz,
            g_,
            vertex_label_writer,
            EdgeLabelWriter<BoostGraph>(g_),
            boost::default_writer(),
            IndexMap());
    return graphviz.str();
}

void DeBruijnGraph::Prune() {
    // Remove low-weight edges not in the reference.
    boost::remove_edge_if(
            [this](const Edge& e) {
                return !g_[e].is_ref && g_[e].weight < DeBruijnGraph_options::MIN_EDGE_SUPPORT;
            },
            g_);

    // Remove vertices not reachable forward from src or backward from sink.
    VertexIterator vbegin, vend;
    std::tie(vbegin, vend) = boost::vertices(g_);
    std::unordered_set<Vertex> all_vertices(vbegin, vend);

    std::set<Vertex> fwd_reachable_vertices, rev_reachable_vertices;
    fwd_reachable_vertices = VerticesReachableFrom(
            source_, g_, IndexMap());
    rev_reachable_vertices = VerticesReachableFrom(
            sink_, boost::make_reverse_graph(g_), IndexMap());

    std::unordered_set<Vertex> reachable_vertices;
    std::set_intersection(
            fwd_reachable_vertices.begin(), fwd_reachable_vertices.end(),
            rev_reachable_vertices.begin(), rev_reachable_vertices.end(),
            std::inserter(reachable_vertices, reachable_vertices.end()));
    for (Vertex v : all_vertices) {
        if (reachable_vertices.find(v) == reachable_vertices.end()) {
            kmer_to_vertex_.erase(g_[v].kmer);
            boost::clear_vertex(v, g_);
            boost::remove_vertex(v, g_);
        }
    }
    RebuildIndexMap();
}
