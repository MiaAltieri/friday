from build import FRIDAY
from graphviz import Digraph


class DeBruijnGraphOptions(object):
    MIN_K = 10
    MAX_K = 100
    STEP_K = 1
    MIN_EDGE_SUPPORT = 2
    MAX_ALLOWED_PATHS = 256

    # base and map quality
    MIN_BASE_QUALITY = 20
    MIN_MAP_QUALITY = 20


class DeBruijnHaplotyper:
    def __init__(self, bam_file_path, fasta_file_path, contig, start, end):
        self.bam_file_path = bam_file_path
        self.fasta_file_path = fasta_file_path
        self.contig = contig
        self.region_start = start
        self.region_end = end

    @staticmethod
    def visualize(graph, output_filename, on_pruned_nodes=False):
        dot = Digraph(comment='dbruijn graph')
        for i in range(1, graph.current_hash_value):
            if on_pruned_nodes is True and i not in graph.good_nodes:
                continue
            kmer = graph.node_hash_int_to_str[i]
            dot.node(str(i), label=kmer)

        for i in range(1, graph.current_hash_value):
            if on_pruned_nodes is True and i not in graph.good_nodes:
                continue
            for j in range(len(graph.out_nodes[i])):
                node_a = i
                node_b = graph.out_nodes[i][j]
                if on_pruned_nodes is True and node_b not in graph.good_nodes:
                    continue
                weight = graph.edges[(node_a, node_b)]
                dot.edge(str(node_a), str(node_b), label=str(weight))

        print("GRAPH SAVED")
        # print(self.dot.source)
        dot.render('outputs/'+output_filename+'.sv')

    def find_haplotypes(self):
        # get the reference from the fasta file
        fasta_handler = FRIDAY.FASTA_handler(self.fasta_file_path)
        reference_sequence = fasta_handler.get_reference_sequence(self.contig, self.region_start, self.region_end)

        min_k, max_k = FRIDAY.DeBruijnGraph.find_min_k_from_ref(reference_sequence,
                                                                DeBruijnGraphOptions.MIN_K,
                                                                DeBruijnGraphOptions.MAX_K,
                                                                DeBruijnGraphOptions.STEP_K)
        # couldn't build ref without cycle
        if min_k == -1:
            return None

        # get the reads from the bam file
        bam_handler = FRIDAY.BAM_handler(self.bam_file_path)
        reads = bam_handler.get_reads(self.contig,
                                      self.region_start,
                                      self.region_end,
                                      DeBruijnGraphOptions.MIN_MAP_QUALITY,
                                      DeBruijnGraphOptions.MIN_BASE_QUALITY)

        # print(reference_sequence)
        for kmer_size in range(min_k, max_k+1, DeBruijnGraphOptions.STEP_K):
            dbg_graph = FRIDAY.DeBruijnGraph(self.region_start, self.region_end)
            haplotypes = dbg_graph.generate_haplotypes(reference_sequence, reads, kmer_size)
            # self.visualize(dbg_graph, 'cpp_dbg', True)
            if haplotypes:
                return reference_sequence, haplotypes

        return reference_sequence, []
