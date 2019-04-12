//
// Created by Kishwar Shafin on 11/1/18.
//

#include "../../headers/image_generator/image_generator.h"

ImageGenerator::ImageGenerator(string reference_sequence,
                               string chromosome_name,
                               long long ref_start,
                               long long ref_end) {
    this->reference_sequence = reference_sequence;
    this->chromosome_name = chromosome_name;
    this->ref_start = ref_start;
    this->ref_end = ref_end;
    this->global_base_color = {{'C', 50}, {'T', 100}, {'G', 150}, {'A', 200}, {'*', 250}, {'.', 10}, {'N', 10}};
}

string ImageGenerator::get_reference_sequence(long long st_pos, long long end_pos) {
    if(st_pos < ref_start) {
        cerr<<"REFERENCE FETCH ERROR ST_POS < REF_START"<<endl;
    }
    if(end_pos >= ref_end) {
        cerr<<"REFERENCE FETCH ERROR END_POS >= REF_END"<<endl;
    }
    long long length = end_pos - st_pos;
    string ref_seq = this->reference_sequence.substr(st_pos - this->ref_start, length);
    return ref_seq;
}

vector<vector<int> > ImageGenerator::read_to_image_row(type_read read,
                                                       long long &read_start,
                                                       long long &read_end,
                                                       bool supports_allele) {
    read_start = -1;
    read_end = -1;
    vector<vector<int> > image_row;

    int base_color, base_qual_color, map_qual_color, strand_color, support_color, mismatch_color;
    map_qual_color = (double) PileupPixels::MAX_COLOR_VALUE *
                     ((double) min(read.mapping_quality, PileupPixels::MAP_QUALITY_CAP) /
                      (double)PileupPixels::MAP_QUALITY_CAP);
    strand_color = read.flags.is_reverse ? PileupPixels::MAX_COLOR_VALUE : 70;
    support_color = supports_allele ? PileupPixels::MAX_COLOR_VALUE : 70;

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

                        base_qual_color = (double)PileupPixels::MAX_COLOR_VALUE *
                                          ((double)min(read.base_qualities[read_index], PileupPixels::BASE_QUALITY_CAP)
                                           / (double)PileupPixels::BASE_QUALITY_CAP);
                        mismatch_color = PileupPixels::MAX_COLOR_VALUE;
                        image_row.push_back({ base_color, base_qual_color, map_qual_color, strand_color,
                                              mismatch_color, support_color});
                    } else if(ref_position >= ref_start && ref_position < ref_end) {
                        if(read_start == -1) {
                            read_start = ref_position;
                            read_end = ref_position;
                        }
                        read_start = min(read_start, ref_position);
                        read_end = max(read_end, ref_position);

                        base_qual_color = (double)PileupPixels::MAX_COLOR_VALUE *
                                          ((double)min(read.base_qualities[read_index], PileupPixels::BASE_QUALITY_CAP)
                                           / (double)PileupPixels::BASE_QUALITY_CAP);

                        base_color = global_base_color[read.sequence[read_index]];
                        mismatch_color = 70;
                        image_row.push_back({ base_color, base_qual_color, map_qual_color, strand_color,
                                              mismatch_color, support_color});
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
                    base_qual_color = (double)PileupPixels::MAX_COLOR_VALUE *
                                      ((double)min(base_quality, PileupPixels::BASE_QUALITY_CAP)
                                       / (double)PileupPixels::BASE_QUALITY_CAP);

                    base_color = global_base_color['*'];

                    if(!image_row.empty()) {
                        image_row.pop_back();
                    }
                    mismatch_color = PileupPixels::MAX_COLOR_VALUE;

                    image_row.push_back({ base_color, base_qual_color, map_qual_color, strand_color,
                                          mismatch_color, support_color});

//                    cout<<"INSERT: "<<ref_position<<" "<<ref<<" "<<alt<<" "<<int(alt_color)<<endl;
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
                    base_color = global_base_color['.'];

                    if(!image_row.empty()) {
                        image_row.pop_back();
                    }
                    if (read_start == -1) {
                        read_start = ref_position - 1;
                        read_end = ref_position - 1;
                    }
                    mismatch_color = PileupPixels::MAX_COLOR_VALUE;

                    base_qual_color = (double)PileupPixels::MAX_COLOR_VALUE *
                                      ((double)min(20, PileupPixels::BASE_QUALITY_CAP)
                                       / (double)PileupPixels::BASE_QUALITY_CAP);

                    for(int i= -1; i<cigar.length; i++) {
                        read_start = min(read_start, ref_position + i);
                        read_end = max(read_end, ref_position + i);
                        if(i==-1) {
                            image_row.push_back({global_base_color['*'], base_qual_color, map_qual_color, strand_color,
                                                 mismatch_color, support_color});
                        }else {
                            image_row.push_back({global_base_color['.'], base_qual_color, map_qual_color, strand_color,
                                                 mismatch_color, support_color});
                        }
                    }

//                    cout<<"DEL: "<<ref_position<<" "<<ref<<" "<<alt<<" "<<int(alt_color)<<endl;
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
    return image_row;
}

vector<vector<int> > ImageGenerator::get_reference_row(string ref_seq, int left_pad, int right_pad) {
    vector<vector<int> > reference_row;
    reference_row.insert(reference_row.end(), left_pad, {0, 0, 0, 0, 0, 0});
    for(auto&base: ref_seq) {
        int base_color = global_base_color[base];
        reference_row.push_back({base_color, PileupPixels::MAX_COLOR_VALUE, PileupPixels::MAX_COLOR_VALUE,
                                 70, 70, PileupPixels::MAX_COLOR_VALUE});
    }
    reference_row.insert(reference_row.end(), right_pad, {0, 0, 0, 0, 0, 0});
    assert(reference_row.size() == 2 * PileupPixels::CONTEXT_SIZE);
    return reference_row;
}

long long ImageGenerator::overlap_length_between_ranges(pair<long long, long long> range_a,
                                                        pair<long long, long long> range_b) {
    return max((long long)0, (min(range_a.second, range_b.second) - max(range_a.first, range_b.first)));
}


void ImageGenerator::assign_read_to_image(PileupImage& pileup_image,
                                          vector<vector<int> >& image_row,
                                          long long read_start,
                                          long long read_end,
                                          int left_pad,
                                          int right_pad) {
    if(pileup_image.image.size() >= PileupPixels::IMAGE_HEIGHT) return;
    long long read_start_index = 0, read_end_index = image_row.size();
    int left_empties = 0, right_empties = 0;

    if(read_start < pileup_image.start_pos) {
        // read starts before the window
        read_start_index = pileup_image.start_pos - read_start;
    }
    else if(read_start > pileup_image.start_pos) {
        // read starts after the window, so we need to add empty pixels to the left
        left_empties = read_start - pileup_image.start_pos;
    }

    if(read_end >= pileup_image.end_pos) {
        // read goes beyond the window
        read_end_index = pileup_image.end_pos - read_start;
    } else if(read_end < pileup_image.end_pos) {
        // read ends well before the end position
        right_empties = pileup_image.end_pos - read_end - 1;
    }
    // core values from the read
    vector<vector<int> > core;
    core.insert(core.end(), image_row.begin() + read_start_index, image_row.begin() + read_end_index);

    vector<vector<int> > read_row;
    read_row.insert(read_row.end(), left_pad, {0, 0, 0, 0, 0, 0});
    read_row.insert(read_row.end(), left_empties, {0, 0, 0, 0, 0, 0});
    read_row.insert(read_row.end(), core.begin(), core.end());
    read_row.insert(read_row.end(), right_empties, {0, 0, 0, 0, 0, 0});
    read_row.insert(read_row.end(), right_pad, {0, 0, 0, 0, 0, 0});
    assert(read_row.size() == 2 * PileupPixels::CONTEXT_SIZE);

    pileup_image.image.push_back(read_row);
}


PileupImage ImageGenerator::create_image(PositionalCandidateRecord candidate,
                                                 vector< pair<type_read, bool> > reads,
                                                 int genotype) {
    // container for all the images we will generate
    PileupImage pileup_image;
    pileup_image.chromosome_name = candidate.chromosome_name;
    pileup_image.start_pos = max(ref_start, candidate.pos_start - PileupPixels::CONTEXT_SIZE);
    pileup_image.end_pos = min(ref_end, candidate.pos_start + PileupPixels::CONTEXT_SIZE);
    pileup_image.label = genotype;


    // calculate how much left padding required so the variant of interest can be in the middle column
    int left_band = candidate.pos_start - pileup_image.start_pos;
    int left_pad = PileupPixels::CONTEXT_SIZE - left_band;

    int right_band = pileup_image.end_pos - candidate.pos_start;
    int right_pad = PileupPixels::CONTEXT_SIZE - right_band;

    // initialize all the pileup images with reference sequence
    string ref_seq = get_reference_sequence(pileup_image.start_pos, pileup_image.end_pos);
    pileup_image.image.insert(pileup_image.image.end(), PileupPixels::REF_ROW_BAND,
                              get_reference_row(ref_seq, left_pad, right_pad));
//    // now iterate through each of the reads and add it to different windows if read overlaps

    for(int i=0; i<reads.size(); i++) {
        type_read read = reads[i].first;
        bool supported = reads[i].second;
//        cout<<read.read_id<<" "<<read.pos<<" "<<read.pos_end<<" "<<supported<<endl;
        long long read_start, read_end;
        // convert the read to a pileup row
        vector<vector<int> > image_row = read_to_image_row(read, read_start, read_end, supported);
        if(read_end < pileup_image.start_pos || read_end > pileup_image.end_pos) {
            continue;
        }
        assign_read_to_image(pileup_image, image_row, read_start, read_end, left_pad, right_pad);
//        cout<<start_ref_position<<" "<<end_ref_position<<endl;
    }
    vector<vector<int> > empty_row;
    empty_row.insert(empty_row.end(), PileupPixels::CONTEXT_SIZE * 2, {0, 0, 0, 0, 0, 0});

    int empties_needed = PileupPixels::IMAGE_HEIGHT - pileup_image.image.size();
    if(empties_needed > 0){
        pileup_image.image.insert(pileup_image.image.end(), empties_needed, empty_row);
    }
    assert(pileup_image.image.size() == PileupPixels::IMAGE_HEIGHT);
    return pileup_image;
}
