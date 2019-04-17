from build import FRIDAY
from graphviz import Digraph

from modules.python.Options import DeBruijnGraphOptions


class DeBruijnHaplotyper:
    def __init__(self, fasta_handler, contig, start, end):
        self.fasta_handler = fasta_handler
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
            if i not in graph.out_nodes:
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

    def find_haplotypes(self, reads):
        # get the reference from the fasta file
        reference_sequence = self.fasta_handler.get_reference_sequence(self.contig, self.region_start, self.region_end)

        haplotypes = FRIDAY.DebruijnGraphHelper().Build(reference_sequence, reads, self.region_start, self.region_end)

        if not haplotypes:
            return reference_sequence, []

        return reference_sequence, haplotypes
