//
// Created by Kishwar Shafin on 10/18/18.
//

#ifndef FRIDAY_PYBIND_API_H
#define FRIDAY_PYBIND_API_H

#include "dataio/fasta_handler.h"
#include "dataio/bam_handler.h"
#include "dataio/vcf_handler.h"
#include "local_reassembly/active_region_finder.h"
#include "local_reassembly/aligner.h"
#include "image_generator/image_generator.h"
#include "candidate_finding/candidate_finder.h"
#include "local_reassembly/debruijn_boost.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
namespace py = pybind11;

PYBIND11_MODULE(FRIDAY, m) {
        py::class_<PileupImage>(m, "PileupImage")
            .def(py::init<>())
            .def(py::init<const string &, const int &, vector<vector<vector<uint8_t> > > & >())
            .def_readwrite("chromosome_name", &PileupImage::chromosome_name)
            .def_readwrite("start_pos", &PileupImage::start_pos)
            .def_readwrite("end_pos", &PileupImage::end_pos)
            .def_readwrite("image", &PileupImage::image)
            .def_readwrite("label", &PileupImage::label)
            .def_readwrite("name", &PileupImage::name)
            .def(py::pickle(
                    [](const PileupImage &p) { // __getstate__
                        /* Return a tuple that fully encodes the state of the object */
                        return py::make_tuple(p.name, p.label, p.image);
                    },
                    [](py::tuple t) { // __setstate__
                        if (t.size() != 3)
                            throw std::runtime_error("Invalid state!");

                        /* Create a new C++ instance */
                        PileupImage p(t[0].cast<string>(), t[1].cast<int>(), t[2].cast<vector<vector<vector<uint8_t> > > >());

                        /* Assign any additional state */
                        //dp.setExtra(t[1].cast<int>());

                        return p;
                    }
            ));

        py::class_<ImageGenerator>(m, "ImageGenerator")
            .def(py::init<const string &, const string &, long long &, long long&>())
            .def("create_image", &ImageGenerator::create_image);

//        py::class_<SummaryGenerator>(m, "SummaryGenerator")
//            .def(py::init<const string &, const string &, long long &, long long &>())
//            .def_readwrite("genomic_pos", &SummaryGenerator::genomic_pos)
//            .def_readwrite("labels", &SummaryGenerator::labels)
//            .def_readwrite("image", &SummaryGenerator::image)
//            .def("generate_train_summary", &SummaryGenerator::generate_train_summary)
//            .def("generate_summary", &SummaryGenerator::generate_summary);
        py::class_<PositionalCandidateRecord>(m, "PositionalCandidateRecord")
            .def(py::init<>())
            .def(py::init<const string &, long long &, long long&,
                          const string &, const vector<string> &,
                          const vector<int> &, const vector<double> &,
                          const int &, const vector<string> >())
            .def("set_genotype", &PositionalCandidateRecord::set_genotype)
            .def("add_image_name", &PositionalCandidateRecord::add_image_name)
            .def_readwrite("name", &PositionalCandidateRecord::name)
            .def_readwrite("chromosome_name", &PositionalCandidateRecord::chromosome_name)
            .def_readwrite("pos_start", &PositionalCandidateRecord::pos_start)
            .def_readwrite("pos_end", &PositionalCandidateRecord::pos_end)
            .def_readwrite("ref", &PositionalCandidateRecord::ref)
            .def_readwrite("alternate_alleles", &PositionalCandidateRecord::alternate_alleles)
            .def_readwrite("allele_depths", &PositionalCandidateRecord::allele_depths)
            .def_readwrite("allele_frequencies", &PositionalCandidateRecord::allele_frequencies)
            .def_readwrite("genotype", &PositionalCandidateRecord::genotype)
            .def_readwrite("read_support_alleles", &PositionalCandidateRecord::read_support_alleles)
            .def_readwrite("depth", &PositionalCandidateRecord::depth)
            .def_readwrite("image_names", &PositionalCandidateRecord::image_names)
            .def_readwrite("labeled", &PositionalCandidateRecord::labeled)
            .def(py::pickle(
                    [](const PositionalCandidateRecord &p) { // __getstate__
                        /* Return a tuple that fully encodes the state of the object */
                        return py::make_tuple(p.chromosome_name, p.pos_start, p.pos_end,
                                              p.ref, p.alternate_alleles,
                                              p.allele_depths, p.allele_frequencies,
                                              p.depth, p.image_names);
                    },
                    [](py::tuple t) { // __setstate__
                        if (t.size() != 9)
                            throw std::runtime_error("Invalid state!");

                        /* Create a new C++ instance */
                        PositionalCandidateRecord p(t[0].cast<string>(), t[1].cast<long long>(), t[2].cast<long long>(),
                                                    t[3].cast<string>(), t[4].cast< vector<string> >(),
                                                    t[5].cast<vector<int> >(), t[6].cast<vector<double> >(),
                                                    t[7].cast<int>(), t[8].cast<vector<string> >());

                        /* Assign any additional state */
                        //dp.setExtra(t[1].cast<int>());

                        return p;
                    }
            ));


        py::class_<CandidateAllele>(m, "CandidateAllele")
            .def(py::init<>())
            .def_readwrite("ref", &CandidateAllele::ref)
            .def_readwrite("alt", &CandidateAllele::alt)
            .def_readwrite("alt_type", &CandidateAllele::alt_type);

        py::class_<Candidate>(m, "Candidate")
            .def(py::init<>())
            .def("print", &Candidate::print)
            .def("set_genotype", &Candidate::set_genotype)
            .def_readwrite("pos", &Candidate::pos)
            .def_readwrite("pos_end", &Candidate::pos_end)
            .def_readwrite("genotype", &Candidate::genotype)
            .def_readwrite("supporting_read_ids", &Candidate::supporting_read_ids)
            .def_readwrite("allele", &Candidate::allele);


//        py::class_<PositionalCandidateRecord>(m, "PositionalCandidateRecord")
//            .def(py::init<>())
//            .def(py::init<const string &, long long &, long long&, const string &,
//                 const string &, const string &, const int &, const int &>())
//            .def("print", &PositionalCandidateRecord::print)
//            .def("set_genotype", &PositionalCandidateRecord::set_genotype)
//            .def("get_candidate_record", &PositionalCandidateRecord::get_candidate_record)
//            .def_readwrite("chromosome_name", &PositionalCandidateRecord::chromosome_name)
//            .def_readwrite("pos", &PositionalCandidateRecord::pos)
//            .def_readwrite("pos_end", &PositionalCandidateRecord::pos_end)
//            .def_readwrite("ref", &PositionalCandidateRecord::ref)
//            .def_readwrite("alt1", &PositionalCandidateRecord::alt1)
//            .def_readwrite("alt2", &PositionalCandidateRecord::alt2)
//            .def_readwrite("alt1_type", &PositionalCandidateRecord::alt1_type)
//            .def_readwrite("alt2_type", &PositionalCandidateRecord::alt2_type)
//            .def(py::pickle(
//                    [](const PositionalCandidateRecord &p) { // __getstate__
//                        /* Return a tuple that fully encodes the state of the object */
//                        return py::make_tuple(p.chromosome_name, p.pos, p.pos_end, p.ref, p.alt1, p.alt2, p.alt1_type, p.alt2_type);
//                    },
//                    [](py::tuple t) { // __setstate__
//                        if (t.size() != 8)
//                            throw std::runtime_error("Invalid state!");
//
//                        /* Create a new C++ instance */
//                        PositionalCandidateRecord p(t[0].cast<string>(), t[1].cast<long long>(), t[2].cast<long long>(), t[3].cast<string>(),
//                                                    t[4].cast<string>(), t[5].cast<string>(), t[6].cast<int>(), t[7].cast<int>());
//
//                        /* Assign any additional state */
//                        //dp.setExtra(t[1].cast<int>());
//
//                        return p;
//                    }
//            ));


        // Candidate finder
        py::class_<CandidateFinder>(m, "CandidateFinder")
            .def(py::init<const string &, const string &, long long &, long long&, long long&, long long&>())
            .def_readwrite("position_to_read_map", &CandidateFinder::position_to_read_map)
            .def("find_candidates", &CandidateFinder::find_candidates);

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
        py::class_<DeBruijnGraphHelper>(m, "DebruijnGraphHelper")
            .def(py::init<>())
            .def("Build", &DeBruijnGraphHelper::Build);

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
            .def("__lt__", &type_read::operator<, py::is_operator())
            .def_readwrite("pos", &type_read::pos)
            .def_readwrite("pos_end", &type_read::pos_end)
            .def_readwrite("query_name", &type_read::query_name)
            .def_readwrite("read_id", &type_read::read_id)
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
            .def("get_sample_names", &BAM_handler::get_sample_names)
            .def("get_chromosome_sequence_names_with_length", &BAM_handler::get_chromosome_sequence_names_with_length)
            .def("get_reads", &BAM_handler::get_reads);

        // FASTA handler API
        py::class_<FASTA_handler>(m, "FASTA_handler")
            .def(py::init<const string &>())
            .def("get_reference_sequence", &FASTA_handler::get_reference_sequence)
            .def("get_chromosome_sequence_length", &FASTA_handler::get_chromosome_sequence_length)
            .def("get_chromosome_names", &FASTA_handler::get_chromosome_names);

        // VCF handler API
        py::class_<VCF_handler>(m, "VCF_handler")
            .def(py::init<const string &>())
            .def("get_vcf_records", &VCF_handler::get_vcf_records)
            .def("get_positional_vcf_records", &VCF_handler::get_positional_vcf_records);

        // VCF handler API
        py::class_<type_vcf_record>(m, "type_vcf_record")
            .def_readwrite("chromosome_name", &type_vcf_record::chromosome_name)
            .def_readwrite("start_pos", &type_vcf_record::start_pos)
            .def_readwrite("end_pos", &type_vcf_record::end_pos)
            .def_readwrite("id", &type_vcf_record::id)
            .def_readwrite("qual", &type_vcf_record::qual)
            .def_readwrite("is_filter_pass", &type_vcf_record::is_filter_pass)
            .def_readwrite("sample_name", &type_vcf_record::sample_name)
            .def_readwrite("genotype", &type_vcf_record::genotype)
            .def_readwrite("filters", &type_vcf_record::filters)
            .def_readwrite("alleles", &type_vcf_record::alleles);


        py::class_<type_alt_allele>(m, "type_alt_allele")
            .def_readonly("ref", &type_alt_allele::ref)
            .def_readonly("alt_allele", &type_alt_allele::alt_allele)
            .def_readonly("alt_type", &type_alt_allele::alt_type)
            .def("get_ref", &type_alt_allele::get_ref)
            .def("get_alt_allele", &type_alt_allele::get_alt_allele)
            .def("get_alt_type", &type_alt_allele::get_alt_type);

        py::class_<type_positional_vcf_record>(m, "type_positional_vcf_record")
            .def_readonly("chromosome_name", &type_positional_vcf_record::chromosome_name)
            .def_readonly("start_pos", &type_positional_vcf_record::start_pos)
            .def_readonly("end_pos", &type_positional_vcf_record::end_pos)
            .def_readonly("id", &type_positional_vcf_record::id)
            .def_readonly("qual", &type_positional_vcf_record::qual)
            .def_readonly("is_filter_pass", &type_positional_vcf_record::is_filter_pass)
            .def_readonly("sample_name", &type_positional_vcf_record::sample_name)
            .def_readonly("genotype", &type_positional_vcf_record::genotype)
            .def_readonly("filters", &type_positional_vcf_record::filters)
            .def_readonly("alt_allele", &type_positional_vcf_record::alt_allele);

}
#endif //FRIDAY_PYBIND_API_H
