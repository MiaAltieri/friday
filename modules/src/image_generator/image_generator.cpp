//
// Created by Kishwar Shafin on 11/1/18.
//

#include "../../headers/image_generator/image_generator.h"

ImageGenerator::ImageGenerator(string reference_sequence,
                               string chromosome_name,
                               long long ref_start,
                               long long ref_end,
                               map<long long, PositionalCandidateRecord> all_positional_candidates) {
    this->reference_sequence = reference_sequence;
    this->chromosome_name = chromosome_name;
    this->ref_start = ref_start;
    this->ref_end = ref_end;
    this->all_positional_candidates = all_positional_candidates;
    this->global_base_color = {{'A', 250}, {'C', 30}, {'G', 180}, {'T', 100}, {'.', 0}, {'*', 0}, {'N', 10}};
}

pair<vector< vector<uint8_t> >, pair<long long, long long> > ImageGenerator::read_to_image_row(type_read read) {
    long long read_start = -1;
    long long read_end = -1;
    vector< vector<uint8_t> > image_row;

    uint8_t base_color, base_qual_color, map_qual_color, strand_color, alt_color;
    map_qual_color = uint8_t(PileupPixels::MAX_COLOR_VALUE *
            (min(read.mapping_quality, PileupPixels::MAP_QUALITY_CAP) / PileupPixels::MAP_QUALITY_CAP));
    strand_color = read.flags.is_reverse ? 240 : 0;

    int read_index = 0;
    long long ref_position = read.pos;
    int cigar_index = 0;
    int base_quality = 0;
    long long reference_index;


    for(auto &cigar: read.cigar_tuples) {
        switch (cigar.operation) {
            case CIGAR_OPERATIONS::EQUAL:
            case CIGAR_OPERATIONS::DIFF:
            case CIGAR_OPERATIONS::MATCH:
                cigar_index = 0;
                if(ref_position < ref_start) {
                    cigar_index = min(ref_start - ref_position, (long long)cigar.length);
                    read_index += cigar_index;
                    ref_position += cigar_index;
                }
                for(int i=cigar_index; i < cigar.length ; i++) {
                    reference_index = ref_position - ref_start;
                    //read.base_qualities[read_index] base quality
                    if(ref_position >= ref_start && ref_position < ref_end &&
                       reference_sequence[reference_index] != read.sequence[read_index]) {
                        if(read_start == -1) {
                            read_start = ref_position;
                            read_end = ref_position;
                        }
                        read_start = min(read_start, ref_position);
                        read_end = max(read_end, ref_position);

                        base_color = global_base_color[read.sequence[read_index]];

                        base_qual_color = uint8_t(PileupPixels::MAX_COLOR_VALUE *
                                (min(read.base_qualities[read_index], PileupPixels::BASE_QUALITY_CAP)
                                 / PileupPixels::BASE_QUALITY_CAP));

                        image_row.push_back({ base_color, base_qual_color, map_qual_color, strand_color});

                        // process the SNP allele here
                        string ref(1, reference_sequence[reference_index]);
                        string alt(1, read.sequence[read_index]);
                        Candidate candidate_alt(ref_position, ref_position + 1, ref, alt, AlleleType::SNP_ALLELE);


//                        cout<<"SNP: "<<ref_position<<" "<<ref<<" "<<alt<<endl;
                    } else if(ref_position >= ref_start && ref_position < ref_end) {
                        if(read_start == -1) {
                            read_start = ref_position;
                            read_end = ref_position;
                        }
                        read_start = min(read_start, ref_position);
                        read_end = max(read_end, ref_position);

                        base_qual_color = uint8_t(PileupPixels::MAX_COLOR_VALUE *
                                                     (min(read.base_qualities[read_index], PileupPixels::BASE_QUALITY_CAP)
                                                      / PileupPixels::BASE_QUALITY_CAP));

                        base_color = global_base_color[read.sequence[read_index]];

                        image_row.push_back({ base_color, base_qual_color, map_qual_color, strand_color});
                    }
                    read_index += 1;
                    ref_position += 1;
                }
                break;
            case CIGAR_OPERATIONS::IN:
                base_quality = *std::min_element(read.base_qualities.begin() + read_index,
                                                 read.base_qualities.begin() + (read_index + cigar.length));
                reference_index = ref_position - ref_start - 1;

                if(ref_position - 1 >= ref_start &&
                   ref_position - 1 < ref_end) {
                    if(read_start == -1) {
                        read_start = ref_position - 1;
                        read_end = ref_position - 1;
                    }
                    read_start = min(read_start, ref_position - 1);
                    read_end = max(read_end, ref_position - 1);

                    // process insert allele here
                    string ref = reference_sequence.substr(reference_index, 1);
                    string alt;
                    if(read_index - 1 >= 0) alt = read.sequence.substr(read_index - 1, cigar.length + 1);
                    else alt = ref + read.sequence.substr(read_index, cigar.length);

                    base_qual_color = uint8_t(PileupPixels::MAX_COLOR_VALUE *
                                                 (min(base_quality, PileupPixels::BASE_QUALITY_CAP)
                                                  / PileupPixels::BASE_QUALITY_CAP));

                    base_color = global_base_color['*'];

                    if(!image_row.empty()) {
                        image_row.pop_back();
                    }

                    image_row.push_back({ base_color, base_qual_color, map_qual_color, strand_color});

//                    cout<<"INSERT: "<<ref_position<<" "<<ref<<" "<<alt<<" "<<endl;
                }
                read_index += cigar.length;
                break;
            case CIGAR_OPERATIONS::REF_SKIP:
            case CIGAR_OPERATIONS::PAD:
            case CIGAR_OPERATIONS::DEL:
//                base_quality = read.base_qualities[max(0, read_index)];
                reference_index = ref_position - ref_start - 1;

                if(ref_position - 1 >= ref_start && ref_position - 1 < ref_end &&
                   ref_position + cigar.length < ref_end) {
                    // process delete allele here
                    string ref = reference_sequence.substr(ref_position - ref_start - 1, 1);
                    string alt = ref + reference_sequence.substr(ref_position - ref_start , cigar.length);
                    Candidate candidate_alt(ref_position - 1, ref_position - 1 + cigar.length, ref, alt,
                                            AlleleType::DELETE_ALLELE);

                    base_color = global_base_color['.'];

                    if(!image_row.empty()) {
                        image_row.pop_back();
                    }

                    image_row.push_back({ base_color, 0, 0, strand_color});

//                    cout<<"DEL: "<<ref_position<<" "<<ref<<" "<<alt<<" "<<endl;
                }
                ref_position += cigar.length;
                break;
            case CIGAR_OPERATIONS::SOFT_CLIP:
                read_index += cigar.length;
                break;
            case CIGAR_OPERATIONS::HARD_CLIP:
                break;

        }
    }
    return make_pair(image_row, make_pair(read_start, read_end));
}