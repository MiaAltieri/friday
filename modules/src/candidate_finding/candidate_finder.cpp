//
// Created by Kishwar Shafin on 10/24/18.
//

#include "../../headers/candidate_finding/candidate_finder.h"


CandidateFinder::CandidateFinder(string reference_sequence,
                                 string chromosome_name,
                                 long long region_start,
                                 long long region_end) {
    this->reference_sequence = reference_sequence;
    this->region_start = region_start;
    this->region_end = region_end;
    this->chromosome_name = chromosome_name;
}

void CandidateFinder::add_read_alleles(type_read &read, vector<int> &coverage) {
    int read_index = 0;
    long long ref_position = read.pos;
    int cigar_index = 0;
    int base_quality = 0;

    for(auto &cigar: read.cigar_tuples) {
        switch (cigar.operation) {
            case CIGAR_OPERATIONS::EQUAL:
            case CIGAR_OPERATIONS::DIFF:
            case CIGAR_OPERATIONS::MATCH:
                cigar_index = 0;
                if(ref_position < region_start) {
                    cigar_index = min(region_start - ref_position, (long long)cigar.length);
                    read_index += cigar_index;
                    ref_position += cigar_index;
                }
                for(int i=cigar_index; i < cigar.length ; i++) {
                    if(ref_position < region_end &&
                       reference_sequence[ref_position - region_start] != read.sequence[read_index] &&
                       read.base_qualities[ref_position - region_start] >= CandidateFinder_options::min_base_quality) {
                        // process the SNP allele here
                        string ref(1, reference_sequence[ref_position - region_start]);
                        string alt(1, read.sequence[read_index]);
                        Candidate candidate_alt(ref_position, ref_position, ref, alt, AlleleType::SNP_ALLELE);
                        if(AlleleFrequencyMap.find(candidate_alt) != AlleleFrequencyMap.end()) {
                            AlleleFrequencyMap[candidate_alt] += 1;
                        } else {
                            AlleleFrequencyMap[candidate_alt] = 1;
                        }
                        if(find(AlleleMap[ref_position].begin(), AlleleMap[ref_position].end(), candidate_alt) ==
                           AlleleMap[ref_position].end()) {
                            AlleleMap[ref_position].push_back(candidate_alt);
                        }
//                        cout<<"SNP: "<<ref_position<<" "<<ref<<" "<<alt<<" "<<AlleleFrequencyMap[candidate_alt]<<endl;
                        coverage[ref_position - region_start] += 1;
                    }
                    else if(ref_position < region_end &&
                            read.base_qualities[ref_position - region_start] >= CandidateFinder_options::min_base_quality) {
                        coverage[ref_position - region_start] += 1;
                    }
                    read_index += 1;
                    ref_position += 1;
                }
                break;

            case CIGAR_OPERATIONS::IN:
                base_quality = *std::min_element(read.base_qualities.begin() + read_index,
                                                 read.base_qualities.begin() + (read_index + cigar.length));
                if(ref_position - 1 >= region_start &&
                   ref_position - 1 < region_end &&
                   base_quality >= ActiveRegionFinder_options::min_base_quality) {
                    // process insert allele here
                    string ref = reference_sequence.substr(ref_position - region_start - 1, 1);
                    string alt = read.sequence.substr(read_index, cigar.length);
                    Candidate candidate_alt(ref_position - 1, ref_position - 1, ref, alt, AlleleType::INSERT_ALLELE);
                    if(AlleleFrequencyMap.find(candidate_alt) != AlleleFrequencyMap.end()) {
                        AlleleFrequencyMap[candidate_alt] += 1;
                    } else {
                        AlleleFrequencyMap[candidate_alt] = 1;
                    }
                    if(find(AlleleMap[ref_position - 1].begin(), AlleleMap[ref_position - 1].end(), candidate_alt) ==
                       AlleleMap[ref_position - 1].end()) {
                        AlleleMap[ref_position - 1].push_back(candidate_alt);
                    }
//                    cout<<"INSERT: "<<ref_position<<" "<<ref<<" "<<alt<<" "<<AlleleFrequencyMap[candidate_alt]<<endl;
                }
                read_index += cigar.length;
                break;

            case CIGAR_OPERATIONS::DEL:
                base_quality = read.base_qualities[max(0, read_index - 1)];
                if(ref_position >= region_start &&
                   ref_position + cigar.length < region_end &&
                   base_quality >= CandidateFinder_options::min_base_quality) {
                    // process delete allele here
                    string ref = reference_sequence.substr(ref_position - region_start - 1, 1);
                    string alt = ref + reference_sequence.substr(ref_position - region_start , cigar.length);
                    Candidate candidate_alt(ref_position - 1, ref_position - 1, ref, alt, AlleleType::DELETE_ALLELE);
                    if(AlleleFrequencyMap.find(candidate_alt) != AlleleFrequencyMap.end()) {
                        AlleleFrequencyMap[candidate_alt] += 1;
                    } else {
                        AlleleFrequencyMap[candidate_alt] = 1;
                    }
                    if(find(AlleleMap[ref_position - 1].begin(), AlleleMap[ref_position - 1].end(), candidate_alt) ==
                       AlleleMap[ref_position - 1].end()) {
                        AlleleMap[ref_position - 1].push_back(candidate_alt);
                    }
//                    cout<<"DELETE: "<<ref_position<<" "<<ref<<" "<<alt<<" "<<AlleleFrequencyMap[candidate_alt]<<endl;
                }

                for(long long pos = ref_position; pos <= min(region_end, ref_position + cigar.length); pos++)
                    coverage[pos - region_start] += 1;
                ref_position += cigar.length;
                break;
            case CIGAR_OPERATIONS::SOFT_CLIP:
                base_quality = *std::min_element(read.base_qualities.begin() + read_index,
                                                 read.base_qualities.begin() + (read_index + cigar.length));
                if(ref_position - 1 >= region_start &&
                   ref_position - 1 < region_end &&
                   base_quality >= CandidateFinder_options::min_base_quality) {
                    //process soft clip allele here
                    string ref = reference_sequence.substr(ref_position - region_start - 1, 1);
                    string alt = read.sequence.substr(read_index, cigar.length);
                    Candidate candidate_alt(ref_position, ref_position - 1, ref, alt, AlleleType::INSERT_ALLELE);
                    if(AlleleFrequencyMap.find(candidate_alt) != AlleleFrequencyMap.end()) {
                        AlleleFrequencyMap[candidate_alt] += 1;
                    } else {
                        AlleleFrequencyMap[candidate_alt] = 1;
                    }
                    if(find(AlleleMap[ref_position - 1].begin(), AlleleMap[ref_position - 1].end(), candidate_alt) ==
                       AlleleMap[ref_position - 1].end()) {
                        AlleleMap[ref_position - 1].push_back(candidate_alt);
                    }
                }
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

vector<PositionalCandidateRecord> CandidateFinder::find_candidates(vector<type_read> reads) {

    vector<PositionalCandidateRecord> all_positional_candidates;

    vector<int> coverage(region_end - region_start + 1, 0);
    for(auto &read:reads) {
        add_read_alleles(read, coverage);
    }

    vector<long long> positions;
    // get all the positions that pass the threshold
    for(int i=0; i < coverage.size(); i++) {
        vector< pair<double, Candidate> > positional_candidates;
        for(auto& candidate: AlleleMap[region_start + i]) {
            int freq = 0;
            if(coverage[i] > 0)
                freq = (int)ceil(100.0 *  ((double) AlleleFrequencyMap[candidate] / (double) coverage[i]));

            if(freq < CandidateFinder_options::freq_threshold)
                continue;

            positional_candidates.push_back(make_pair(freq, candidate));
        }
        if(positional_candidates.empty()) continue;

        sort(positional_candidates.begin(), positional_candidates.end());

        PositionalCandidateRecord pos_candidate;
        pair<double, Candidate> candidates = positional_candidates.back();
        positional_candidates.pop_back();
        pos_candidate.set_positions(this->chromosome_name, candidates.second.pos, candidates.second.pos_end);
        pos_candidate.set_reference(candidates.second.allele.ref);
        pos_candidate.set_alt1(candidates.second.allele.alt, candidates.second.allele.alt_type);

        if(!positional_candidates.empty()) {
            candidates = positional_candidates.back();
            positional_candidates.pop_back();
            pos_candidate.set_alt2(candidates.second.allele.alt, candidates.second.allele.alt_type);
        }

        all_positional_candidates.push_back(pos_candidate);
    }

    return all_positional_candidates;
}