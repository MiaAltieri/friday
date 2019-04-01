//
// Created by Kishwar Shafin on 10/24/18.
//

#include "../../headers/candidate_finding/candidate_finder.h"


CandidateFinder::CandidateFinder(string reference_sequence,
                                 string chromosome_name,
                                 long long region_start,
                                 long long region_end,
                                 long long ref_start,
                                 long long ref_end) {
    this->reference_sequence = reference_sequence;
    this->region_start = region_start;
    this->region_end = region_end;
    this->chromosome_name = chromosome_name;
    this->ref_start = ref_start;
    this->ref_end = ref_end;
    AlleleMap.resize(region_end - region_start + 1);
}

void CandidateFinder::add_read_alleles(type_read &read, vector<int> &coverage) {
    int read_index = 0;
    long long ref_position = read.pos;
    int cigar_index = 0;
    int base_quality = 0;
    long long reference_index, region_index;

//    cout<<CIGAR_OPERATIONS::SOFT_CLIP<<" "<<CIGAR_OPERATIONS::IN<<endl;
//    cout<<read.pos<<" "<<read.pos_end<<endl;
//    for (auto &cigar: read.cigar_tuples) {
//
//        cout<<"("<<cigar.operation<<", "<<cigar.length<<"), ";
//    }
//    cout<<endl;
    for (int cigar_i=0; cigar_i<read.cigar_tuples.size(); cigar_i++) {
        CigarOp cigar = read.cigar_tuples[cigar_i];
        switch (cigar.operation) {
            case CIGAR_OPERATIONS::EQUAL:
            case CIGAR_OPERATIONS::DIFF:
            case CIGAR_OPERATIONS::MATCH:
                cigar_index = 0;
                if (ref_position < region_start) {
                    cigar_index = min(region_start - ref_position, (long long) cigar.length);
                    read_index += cigar_index;
                    ref_position += cigar_index;
                }
                for (int i = cigar_index; i < cigar.length; i++) {
                    reference_index = ref_position - ref_start;
                    region_index = ref_position - region_start;

                    if (ref_position >= region_start && ref_position <= region_end &&
                        reference_sequence[reference_index] != read.sequence[read_index] &&
                        read.base_qualities[read_index] >= CandidateFinder_options::min_base_quality) {
                        //look forward and make sure this is not an anchor base
                        if(i == cigar.length - 1 && cigar_i + 1 < read.cigar_tuples.size()) {
                            CigarOp next_cigar = read.cigar_tuples[cigar_i + 1];
                            if(next_cigar.operation == CIGAR_OPERATIONS::IN ||
                               next_cigar.operation == CIGAR_OPERATIONS::DEL) {
                                // this is an anchor base of a delete or an insert, don't process this.
                                coverage[region_index] += 1;
                            }
                        } else {
                            // process the SNP allele here
                            string ref(1, reference_sequence[reference_index]);
                            string alt(1, read.sequence[read_index]);
                            Candidate candidate_alt(ref_position, ref_position + 1, ref, alt, AlleleType::SNP_ALLELE);
                            if (AlleleFrequencyMap.find(candidate_alt) != AlleleFrequencyMap.end()) {
                                AlleleFrequencyMap[candidate_alt] += 1;
                            } else {
                                AlleleFrequencyMap[candidate_alt] = 1;
                            }

                            ReadSupportMap[candidate_alt].push_back(read.read_id);

                            if (AlleleMap[region_index].find(candidate_alt) == AlleleMap[region_index].end())
                                AlleleMap[region_index].insert(candidate_alt);
//                          cout<<read.pos<<" "<<ref_position<<" "<<region_start<<" "<<i<<" ";
//                          cout<<"SNP: "<<read.query_name<<" "<<ref_position<<" "<<ref<<" "<<alt<<" "<<AlleleFrequencyMap[candidate_alt]<<" "<<AlleleMap[region_index].size()<<endl;
                            coverage[region_index] += 1;
//                            allele_lengths[region_index] = max(allele_lengths[region_index], 1);
                        }
                    } else if (ref_position <= region_end &&
                               read.base_qualities[read_index] >= CandidateFinder_options::min_base_quality) {
                        coverage[region_index] += 1;
                    }
                    read_index += 1;
                    ref_position += 1;
                }
                break;
//            case CIGAR_OPERATIONS::SOFT_CLIP:
            case CIGAR_OPERATIONS::IN:
                base_quality = *std::min_element(read.base_qualities.begin() + read_index,
                                                 read.base_qualities.begin() + (read_index + cigar.length));
                reference_index = ref_position - ref_start - 1;
                region_index = ref_position - region_start - 1;


                if (ref_position - 1 >= region_start &&
                    ref_position - 1 <= region_end) {
                    // process insert allele here
                    string ref = reference_sequence.substr(reference_index, 1);
                    string alt;
                    if (read_index - 1 >= 0) alt = read.sequence.substr(read_index - 1, cigar.length + 1);
                    else alt = ref + read.sequence.substr(read_index, cigar.length);

                    Candidate candidate_alt(ref_position - 1, ref_position, ref, alt, AlleleType::INSERT_ALLELE);
                    if (AlleleFrequencyMap.find(candidate_alt) != AlleleFrequencyMap.end()) {
                        AlleleFrequencyMap[candidate_alt] += 1;
                    } else {
                        AlleleFrequencyMap[candidate_alt] = 1;
                    }

                    ReadSupportMap[candidate_alt].push_back(read.read_id);

                    if (AlleleMap[region_index].find(candidate_alt) == AlleleMap[region_index].end())
                        AlleleMap[region_index].insert(candidate_alt);
//                    allele_lengths[region_index] = max(allele_lengths[region_index], 1);

//                    cout<<"INSERT: "<<ref_position-1<<" "<<ref<<" "<<alt<<" "<<AlleleFrequencyMap[candidate_alt]<<" "<<base_quality<<endl;
                }
                read_index += cigar.length;
                break;

            case CIGAR_OPERATIONS::DEL:
//                base_quality = read.base_qualities[max(0, read_index)];
                reference_index = ref_position - ref_start - 1;
                region_index = ref_position - region_start - 1;

                if (ref_position - 1 >= region_start && ref_position - 1 <= region_end &&
                    ref_position + cigar.length < ref_end) {
                    // process delete allele here
                    string ref = reference_sequence.substr(ref_position - ref_start - 1, cigar.length + 1);
                    string alt;

                    if (read_index - 1 >= 0) alt = read.sequence.substr(read_index - 1, 1);
                    else alt = reference_sequence.substr(ref_position - ref_start - 1, 1);

                    Candidate candidate_alt(ref_position - 1, ref_position - 1 + cigar.length + 1, ref, alt,
                                            AlleleType::DELETE_ALLELE);

                    if (AlleleFrequencyMap.find(candidate_alt) != AlleleFrequencyMap.end()) {
                        AlleleFrequencyMap[candidate_alt] += 1;
                    } else {
                        AlleleFrequencyMap[candidate_alt] = 1;
                    }

                    ReadSupportMap[candidate_alt].push_back(read.read_id);

//                    cout<<"DEL: "<<ref_position<<" "<<ref<<" "<<alt<<" "<<AlleleFrequencyMap[candidate_alt]<<endl;

                    if (AlleleMap[region_index].find(candidate_alt) == AlleleMap[region_index].end())
                        AlleleMap[region_index].insert(candidate_alt);
//                    allele_lengths[region_index] = max(allele_lengths[region_index], cigar.length + 1);
                }

//                if (ref_position >= region_start && ref_position <= region_end) {
//                    for (long long pos = ref_position; pos <= min(region_end, ref_position + cigar.length); pos++)
//                        coverage[pos - region_start] += 1;
//                }

                ref_position += cigar.length;
                break;
            case CIGAR_OPERATIONS::SOFT_CLIP:
//                base_quality = *std::min_element(read.base_qualities.begin() + read_index,
//                                                 read.base_qualities.begin() + (read_index + cigar.length));
//                reference_index = ref_position - ref_start - 1;
//                region_index = ref_position - region_start - 1;
//
//                if(ref_position - 1 >= region_start &&
//                   ref_position - 1 <= region_end &&
//                   base_quality >= ActiveRegionFinder_options::min_base_quality) {
//                    // process insert allele here
//                    string ref = reference_sequence.substr(reference_index, 1);
//                    string alt;
//                    if(read_index - 1 >= 0) alt = read.sequence.substr(read_index - 1, cigar.length + 1);
//                    else alt = ref + read.sequence.substr(read_index, cigar.length);
//
//                    Candidate candidate_alt(ref_position - 1, ref_position, ref, alt, AlleleType::INSERT_ALLELE);
//                    if(AlleleFrequencyMap.find(candidate_alt) != AlleleFrequencyMap.end()) {
//                        AlleleFrequencyMap[candidate_alt] += 1;
//                    } else {
//                        AlleleFrequencyMap[candidate_alt] = 1;
//                    }
//
//                    if(AlleleMap[region_index].find(candidate_alt) ==  AlleleMap[region_index].end())
//                        AlleleMap[region_index].insert(candidate_alt);
//
//                    cout<<"SOFT CLIP: "<<ref_position<<" "<<ref<<" "<<alt<<" "<<AlleleFrequencyMap[candidate_alt]<<endl;
//                }
                read_index += cigar.length;
                break;
            case CIGAR_OPERATIONS::REF_SKIP:
            case CIGAR_OPERATIONS::PAD:
                ref_position += cigar.length;
                break;
            case CIGAR_OPERATIONS::HARD_CLIP:
                break;

        }
    }
}

vector<PositionalCandidateRecord> CandidateFinder::find_candidates(
        vector <type_read> reads) {

    map<long long, vector <Candidate> > all_positional_candidates;
    set<long long> filtered_candidate_positions;

    vector<int> coverage(region_end - region_start + 1, 0);
    vector<int> allele_ends(region_end - region_start + 1, 0);
    vector<PositionalCandidateRecord> all_records;
    for (auto &read:reads) {

//        if(read.read_id == "HWI-D00360:7:H88WKADXX:2:1103:3811:90435/0") {
//            cout<<read.read_id<<" "<<read.mapping_quality<<" "<<read.pos<<endl;
//            cout<<read.sequence<<endl;
//            cout<<read.flags.is_paired<<endl;
//            cout<<read.flags.is_proper_pair<<endl;
//            cout<<read.flags.is_unmapped<<endl;
//            cout<<read.flags.is_mate_unmapped<<endl;
//            cout<<read.flags.is_mate_is_reverse<<endl;
//            cout<<read.flags.is_secondary<<endl;
//            cout<<read.flags.is_qc_failed<<endl;
//            cout<<read.flags.is_duplicate<<endl;
//            cout<<read.flags.is_supplementary<<endl;
//            for(auto &cigar:read.cigar_tuples) {
//                cout<<cigar.operation<<" "<<cigar.length<<endl;
//            }
//            cout<<endl;
//        }

        add_read_alleles(read, coverage);
    }

    int ref_buffer = region_start - ref_start;
    // get all the positions that pass the threshold
    for (long long i = 0; i < coverage.size(); i++) {
        allele_ends[i] = 1;
        // first figure out the longest delete
        for (auto &candidate: AlleleMap[i]) {
            double freq_can = 0.0;
            if (coverage[i] > 0)
                freq_can = 100.0 * ((double) AlleleFrequencyMap[candidate] / (double) coverage[i]);

            if (freq_can >= CandidateFinder_options::freq_threshold &&
                AlleleFrequencyMap[candidate] >= CandidateFinder_options::min_count_threshold) {
                if(candidate.allele.alt_type == DELETE_TYPE) {
                    allele_ends[i] = max(allele_ends[i], (int) candidate.allele.ref.length());
                }
            }
        }

        PositionalCandidateRecord positional_record;
        positional_record.chromosome_name = this->chromosome_name;
        positional_record.pos_start = this->region_start + i;
        positional_record.pos_end = positional_record.pos_start + allele_ends[i];
        positional_record.depth = coverage[i];
        positional_record.name = chromosome_name + "_" +
                                 to_string(positional_record.pos_start) + "_" + to_string(positional_record.pos_end);
        if(positional_record.pos_end > this->ref_end) {
            // goes beyond the reference regions skip;
            continue;
        }
        positional_record.ref = reference_sequence.substr(ref_buffer + i, max(1, allele_ends[i]));

        bool candidate_found = false;
        for (auto &candidate: AlleleMap[i]) {
            int freq_can = 0;
            if (coverage[i] > 0)
                freq_can = 100.0 * ((double) AlleleFrequencyMap[candidate] / (double) coverage[i]);

            if (freq_can >= CandidateFinder_options::freq_threshold &&
                AlleleFrequencyMap[candidate] >= CandidateFinder_options::min_count_threshold) {
                candidate_found = true;
                filtered_candidate_positions.insert(i + this->region_start);

                string alt_allele = candidate.allele.alt;

                if(candidate.pos_end < positional_record.pos_end) {
                    int suffix_start = candidate.pos_end - ref_start;
                    int length = positional_record.pos_end - candidate.pos_end;
                    string ref_suffix = reference_sequence.substr(suffix_start, length);
//                    cout<<suffix_start<<" "<<length<<" "<<ref_suffix<<endl;
                    alt_allele = alt_allele + ref_suffix;
                }
                // allele, depth and frequency
                positional_record.alternate_alleles.push_back(alt_allele);
                positional_record.allele_depths.push_back(AlleleFrequencyMap[candidate]);
                positional_record.allele_frequencies.push_back(((double) AlleleFrequencyMap[candidate] /
                                                                max(1.0, (double) coverage[i])));
                positional_record.read_support_alleles.push_back(ReadSupportMap[candidate]);

//                positional_candidates.push_back(new_candidate);
//
//                cout << "CANDIDATE: " << i + this->region_start << " " << candidate.allele.ref << " " << candidate.pos_end << " "
//                     << alt_allele << " " << " " << AlleleFrequencyMap[candidate]
//                     << " " << coverage[i] << " "<< freq_can << endl;
            } else {
//                cout << "SKIPPED CANDIDATE: " << i + this->region_start << " " << candidate.allele.ref << " "
//                     << candidate.allele.alt << " " << candidate.allele.alt_type << " " << AlleleFrequencyMap[candidate]
//                     << " " << coverage[i] << " "<< freq_can << endl;
            }
        }

        if (!candidate_found) continue;
//        cout<< "{ " << "Chr: " << positional_record.chromosome_name<<
//            " , POS_START: "<<positional_record.pos_start<<
//            " , POS_END: "<<positional_record.pos_end<<
//            " , REF: "<<positional_record.ref<<" , DEPTH: "<<positional_record.depth<<endl;;
//
//        for(int ii=0; ii<positional_record.alternate_alleles.size(); ii++) {
//            cout<<"{ ALT: "<<positional_record.alternate_alleles[ii];
//            cout<<" , AD:  "<<positional_record.allele_depths[ii];
//            cout<<" , VAF: "<<positional_record.allele_frequencies[ii];
//            cout<<"\n SUPPORTED READS: "<<endl;
//            for(int jj=0;jj<positional_record.read_support_alleles[ii].size(); jj++) {
//                cout<<positional_record.read_support_alleles[ii][jj]<<endl;
//            }
//            cout<<" }"<<endl;
//        }
//        cout<<"}"<<endl;
        all_records.push_back(positional_record);

//        sort(positional_candidates.begin(), positional_candidates.end());
//
//        PositionalCandidateRecord pos_candidate;
//        pair<double, Candidate> candidates = positional_candidates.back();
//        positional_candidates.pop_back();
//        pos_candidate.set_positions(this->chromosome_name, candidates.second.pos, candidates.second.pos_end);
//        pos_candidate.set_reference(candidates.second.allele.ref);
//        pos_candidate.set_alt1(candidates.second.allele.alt, candidates.second.allele.alt_type);
//
//        if (!positional_candidates.empty()) {
//            candidates = positional_candidates.back();
//            positional_candidates.pop_back();
//            pos_candidate.set_alt2(candidates.second.allele.alt,
//                                   candidates.second.allele.alt_type,
//                                   candidates.second.pos_end);
//        }

//        all_positional_candidates[i + this->region_start] = positional_candidates;
    }
    return all_records;
}