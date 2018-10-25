//
// Created by Kishwar Shafin on 10/18/18.
//

#ifndef FRIDAY_PYBIND_API_H
#define FRIDAY_PYBIND_API_H

#include "dataio/fasta_handler.h"
#include "dataio/bam_handler.h"
#include "local_reassembly/active_region_finder.h"
#include "local_reassembly/debruijn_graph.h"
#include "local_reassembly/aligner.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

PYBIND11_MODULE(FRIDAY, m) {
        // Alignment CLASS
        py::class_<StripedSmithWaterman::Alignment>(m, "Alignment")
            .def(py::init<>())
            .def_readwrite("best_score", &StripedSmithWaterman::Alignment::sw_score)
            .def_readwrite("best_score2", &StripedSmithWaterman::Alignment::sw_score_next_best)
            .def_readwrite("reference_begin", &StripedSmithWaterman::Alignment::ref_begin)
            .def_readwrite("reference_end", &StripedSmithWaterman::Alignment::ref_end)
            .def_readwrite("query_begin", &StripedSmithWaterman::Alignment::query_begin)
            .def_readwrite("query_end", &StripedSmithWaterman::Alignment::query_end)
            .def_readwrite("ref_end_next_best", &StripedSmithWaterman::Alignment::ref_end_next_best)
            .def_readwrite("mismatches", &StripedSmithWaterman::Alignment::mismatches)
            .def_readwrite("cigar_string", &StripedSmithWaterman::Alignment::cigar_string)
            .def_readwrite("cigar", &StripedSmithWaterman::Alignment::cigar)
            .def("Clear", &StripedSmithWaterman::Alignment::Clear);

        // Filter Class
        py::class_<StripedSmithWaterman::Filter>(m, "Filter")
            .def_readwrite("report_begin_position", &StripedSmithWaterman::Filter::report_begin_position)
            .def_readwrite("report_cigar", &StripedSmithWaterman::Filter::report_cigar)
            .def_readwrite("score_filter", &StripedSmithWaterman::Filter::score_filter)
            .def_readwrite("distance_filter", &StripedSmithWaterman::Filter::distance_filter)
            .def(py::init<>())
            .def(py::init<const bool&, const bool&, const uint16_t&, const uint16_t&>());

        // Aligner Class
        py::class_<StripedSmithWaterman::Aligner>(m, "Aligner")
            .def(py::init<>())
            .def(py::init<const uint8_t&, const uint8_t&, const uint8_t&, const uint8_t&>())
            .def("SetReferenceSequence", &StripedSmithWaterman::Aligner::SetReferenceSequence)
            .def("Align_cpp", &StripedSmithWaterman::Aligner::Align_cpp);

        // SSW alignment class
        py::class_<ReadAligner>(m, "ReadAligner")
            .def(py::init<const int &, const int &, const string &>())
            .def("align_reads", &ReadAligner::align_reads);

        // Debruijn graph class
        py::class_<DeBruijnGraph>(m, "DeBruijnGraph")
            .def(py::init<const long long &, const long long &>())
            .def_readwrite("current_hash_value", &DeBruijnGraph::current_hash_value)
            .def_readwrite("node_hash_int_to_str", &DeBruijnGraph::node_hash_int_to_str)
            .def_readwrite("good_nodes", &DeBruijnGraph::good_nodes)
            .def_readwrite("out_nodes", &DeBruijnGraph::out_nodes)
            .def_readwrite("edges", &DeBruijnGraph::edges)

            .def("generate_haplotypes", &DeBruijnGraph::generate_haplotypes)
            .def("find_min_k_from_ref", &DeBruijnGraph::find_min_k_from_ref);

        // data structure for sequence name and their length
        py::class_<ActiveRegionFinder>(m, "ActiveRegionFinder")
            .def(py::init<const string &, const string &, long long &, long long&>())
            .def("find_active_region", &ActiveRegionFinder::find_active_region);

        // data structure for sequence name and their length
        py::class_<type_sequence>(m, "type_sequence")
            .def_readwrite("sequence_length", &type_sequence::sequence_length)
            .def_readwrite("sequence_name", &type_sequence::sequence_name);

        // data structure for CIGAR operation
        py::class_<CigarOp>(m, "CigarOp")
            .def_readwrite("cigar_op", &CigarOp::operation)
            .def_readwrite("cigar_len", &CigarOp::length);

        // data structure for read attributes aka read flags
        py::class_<type_read_flags>(m, "type_read_flags")
            .def_readwrite("is_paired", &type_read_flags::is_paired)
            .def_readwrite("is_proper_pair", &type_read_flags::is_proper_pair)
            .def_readwrite("is_unmapped", &type_read_flags::is_unmapped)
            .def_readwrite("is_mate_unmapped", &type_read_flags::is_mate_unmapped)
            .def_readwrite("is_reverse", &type_read_flags::is_reverse)
            .def_readwrite("is_mate_is_reverse", &type_read_flags::is_mate_is_reverse)
            .def_readwrite("is_read1", &type_read_flags::is_read1)
            .def_readwrite("is_read2", &type_read_flags::is_read2)
            .def_readwrite("is_secondary", &type_read_flags::is_secondary)
            .def_readwrite("is_qc_failed", &type_read_flags::is_qc_failed)
            .def_readwrite("is_duplicate", &type_read_flags::is_duplicate)
            .def_readwrite("is_supplementary", &type_read_flags::is_supplementary)
            .def(py::init());

        // data structure for read
        py::class_<type_read>(m, "type_read")
            .def_readwrite("pos", &type_read::pos)
            .def_readwrite("pos_end", &type_read::pos_end)
            .def_readwrite("query_name", &type_read::query_name)
            .def_readwrite("flags", &type_read::flags)
            .def_readwrite("sequence", &type_read::sequence)
            .def_readwrite("cigar_tuples", &type_read::cigar_tuples)
            .def_readwrite("mapping_quality", &type_read::mapping_quality)
            .def_readwrite("base_qualities", &type_read::base_qualities)
            .def_readwrite("bad_indicies", &type_read::bad_indicies);

        // bam handler API
        py::class_<BAM_handler>(m, "BAM_handler")
            .def(py::init<const string &>())
            .def("get_chromosome_sequence_names", &BAM_handler::get_chromosome_sequence_names)
            .def("get_chromosome_sequence_names_with_length", &BAM_handler::get_chromosome_sequence_names_with_length)
            .def("get_reads", &BAM_handler::get_reads);

        // FASTA handler API
        py::class_<FASTA_handler>(m, "FASTA_handler")
            .def(py::init<const string &>())
            .def("get_reference_sequence", &FASTA_handler::get_reference_sequence)
            .def("get_chromosome_sequence_length", &FASTA_handler::get_chromosome_sequence_length)
            .def("get_chromosome_names", &FASTA_handler::get_chromosome_names);
}
#endif //FRIDAY_PYBIND_API_H
