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
    this->global_base_color = {{'A', 250}, {'C', 30}, {'G', 180}, {'T', 100}, {'.', 10}, {'*', 20}, {'N', 5}};
}

void ImageGenerator::set_positional_vcf(map<long long, vector<type_positional_vcf_record> > pos_vcf) {
    this->pos_vcf = pos_vcf;
}

int ImageGenerator::get_which_allele(long long pos, string ref, string alt, int alt_type) {
    if(all_positional_candidates.find(pos) == all_positional_candidates.end()) return 0;

    PositionalCandidateRecord candidate = all_positional_candidates[pos];

    if(candidate.ref.compare(ref) == 0 && candidate.alt1.compare(alt) == 0 && candidate.alt1_type == alt_type) {
        return 1;
    }
    if(candidate.ref.compare(ref) == 0 && candidate.alt2.compare(alt) == 0 && candidate.alt2_type == alt_type) {
        return 2;
    }
    return 0;
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

void ImageGenerator::read_to_image_row(type_read read,
                                       long long &read_start,
                                       long long &read_end,
                                       vector<vector<int> >& image_row_alt1,
                                       vector<vector<int> >& image_row_alt2) {
    read_start = -1;
    read_end = -1;

    int base_color, base_qual_color, map_qual_color, strand_color, alt1_color, alt2_color;
    map_qual_color = (double) PileupPixels::MAX_COLOR_VALUE *
            ((double) min(read.mapping_quality, PileupPixels::MAP_QUALITY_CAP) /
             (double)PileupPixels::MAP_QUALITY_CAP);
    strand_color = read.flags.is_reverse ? 240 : 70;

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

                        // process the SNP allele here
                        string ref(1, reference_sequence[reference_index]);
                        string alt(1, read.sequence[read_index]);
                        int which_allele = get_which_allele(ref_position, ref, alt, AlleleType::SNP_ALLELE);

//                        cout<<"SNP: "<<ref_position<<" "<<ref<<" "<<alt<<" "<<which_allele<<" "<<endl;
                        alt1_color = 5;
                        alt2_color = 5;
                        if(which_allele == 1) alt1_color = 254;
                        if(which_allele == 2) alt2_color = 254;

                        image_row_alt1.push_back({base_color, base_qual_color, map_qual_color, strand_color, alt1_color});
                        image_row_alt2.push_back({base_color, base_qual_color, map_qual_color, strand_color, alt2_color});

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
                        alt1_color = 5;
                        alt2_color = 5;
                        image_row_alt1.push_back({base_color, base_qual_color, map_qual_color, strand_color, alt1_color});
                        image_row_alt2.push_back({base_color, base_qual_color, map_qual_color, strand_color, alt2_color});
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

                    base_qual_color = (double)PileupPixels::MAX_COLOR_VALUE *
                                      ((double)min(base_quality, PileupPixels::BASE_QUALITY_CAP)
                                       / (double)PileupPixels::BASE_QUALITY_CAP);

                    base_color = global_base_color['*'];

                    if(!image_row_alt1.empty()) image_row_alt1.pop_back();
                    if(!image_row_alt2.empty()) image_row_alt2.pop_back();

                    int which_allele = get_which_allele(ref_position - 1, ref, alt, AlleleType::INSERT_ALLELE);

                    alt1_color = 5;
                    alt2_color = 5;
                    if(which_allele == 1) alt1_color = 254;
                    if(which_allele == 2) alt2_color = 254;
                    image_row_alt1.push_back({base_color, base_qual_color, map_qual_color, strand_color, alt1_color});
                    image_row_alt2.push_back({base_color, base_qual_color, map_qual_color, strand_color, alt2_color});

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
                    string ref = reference_sequence.substr(ref_position - ref_start - 1, 1);
                    string alt = ref + reference_sequence.substr(ref_position - ref_start , cigar.length);
                    Candidate candidate_alt(ref_position - 1, ref_position - 1 + cigar.length, ref, alt,
                                            AlleleType::DELETE_ALLELE);
                    base_color = global_base_color['.'];

                    if(!image_row_alt1.empty()) image_row_alt1.pop_back();
                    if(!image_row_alt2.empty()) image_row_alt2.pop_back();

                    if (read_start == -1) {
                        read_start = ref_position - 1;
                        read_end = ref_position - 1;
                    }
                    int which_allele = get_which_allele(ref_position - 1, ref, alt, AlleleType::DELETE_ALLELE);
                    alt1_color = 5;
                    alt2_color = 5;
                    if(which_allele == 1) alt1_color = 254;
                    if(which_allele == 2) alt2_color = 254;
                    image_row_alt1.push_back({base_color, base_qual_color, map_qual_color, strand_color, alt1_color});
                    image_row_alt2.push_back({base_color, base_qual_color, map_qual_color, strand_color, alt2_color});


                    base_qual_color = (double)PileupPixels::MAX_COLOR_VALUE *
                                      ((double)min(20, PileupPixels::BASE_QUALITY_CAP)
                                       / (double)PileupPixels::BASE_QUALITY_CAP);
                    alt1_color = 5;
                    alt2_color = 5;
                    for(int i= 0; i<cigar.length; i++) {
                        read_start = min(read_start, ref_position + i);
                        read_end = max(read_end, ref_position + i);
                        image_row_alt1.push_back({base_color, base_qual_color, map_qual_color, strand_color, alt1_color});
                        image_row_alt2.push_back({base_color, base_qual_color, map_qual_color, strand_color, alt2_color});
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
}

vector<vector<int> > ImageGenerator::get_reference_row(string ref_seq) {
    vector<vector<int> > reference_row;
    for(auto&base: ref_seq) {
        int base_color = global_base_color[base];
        reference_row.push_back({base_color, PileupPixels::MAX_COLOR_VALUE, PileupPixels::MAX_COLOR_VALUE, 70, PileupPixels::MAX_COLOR_VALUE});
    }
    return reference_row;
}

long long ImageGenerator::overlap_length_between_ranges(pair<long long, long long> range_a,
                                                        pair<long long, long long> range_b) {
    return max((long long)0, (min(range_a.second, range_b.second) - max(range_a.first, range_b.first)));
}

void ImageGenerator::assign_read_to_window(PileupImage& pileup_image,
                                           vector<vector<int> > image_row_alt1,
                                           vector<vector<int> > image_row_alt2,
                                           long long read_start,
                                           long long read_end){
    if(pileup_image.image_alt1.size() >= PileupPixels::IMAGE_HEIGHT) return;
    long long read_start_index = 0, read_end_index = image_row_alt1.size();
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
    vector<vector<int> > core_alt1;
    vector<vector<int> > core_alt2;
    core_alt1.insert(core_alt1.end(), image_row_alt1.begin() + read_start_index, image_row_alt1.begin() + read_end_index);
    core_alt2.insert(core_alt2.end(), image_row_alt2.begin() + read_start_index, image_row_alt2.begin() + read_end_index);

    vector<vector<int> > window_row_alt1, window_row_alt2;
    window_row_alt1.insert(window_row_alt1.end(), left_empties, {0, 0, 0, 0, 0});
    window_row_alt1.insert(window_row_alt1.end(), core_alt1.begin(), core_alt1.end());
    window_row_alt1.insert(window_row_alt1.end(), right_empties, {0, 0, 0, 0, 0});

    window_row_alt2.insert(window_row_alt2.end(), left_empties, {0, 0, 0, 0, 0});
    window_row_alt2.insert(window_row_alt2.end(), core_alt2.begin(), core_alt2.end());
    window_row_alt2.insert(window_row_alt2.end(), right_empties, {0, 0, 0, 0, 0});

    pileup_image.image_alt1.push_back(window_row_alt1);
    pileup_image.image_alt2.push_back(window_row_alt2);
}


int ImageGenerator::get_image_label(int gt1, int gt2) {
    if(gt1 == Genotype::HOM && gt2 == Genotype::HOM) return 0;
    if(gt1 == Genotype::HET && gt2 == Genotype::HOM) return 1;
    if(gt1 == Genotype::HOM_ALT && gt2 == Genotype::HOM) return 2;
    if(gt1 == Genotype::HOM && gt2 == Genotype::HET) return 3;
    if(gt2 == Genotype::HOM_ALT && gt1 == Genotype::HOM) return 4;
    if(gt1 == Genotype::HET && gt2 == Genotype::HET) return 5;

    return 0;
}

void ImageGenerator::set_window_labels(pair<long long, long long> window,
                                       vector<int>& label_alt1,
                                       vector<int>& label_alt2) {
    for(long long pos = window.first; pos < window.second; pos++) {
        if(all_positional_candidates.find(pos) == all_positional_candidates.end()) {
            label_alt1.push_back(0);
            label_alt2.push_back(0);
            continue;
        }

        PositionalCandidateRecord candidate = all_positional_candidates[pos];

        int gt1 = Genotype::HOM, gt2 = Genotype::HOM;
        for(auto &vcf_rec:pos_vcf[pos]) {
            for(auto &alt: vcf_rec.alt_allele) {
                if(alt.ref.compare(candidate.ref) == 0
                   && alt.alt_allele.compare(candidate.alt1) == 0
                   && alt.alt_type == candidate.alt1_type) {
                    if(vcf_rec.genotype[0] == 1 && vcf_rec.genotype[1] == 1) {
                        gt1 = Genotype::HOM_ALT;
                    }else if(vcf_rec.genotype[0] == 1 || vcf_rec.genotype[1] == 1) {
                        gt1 = Genotype::HET;
                    }
                }

                if(alt.ref.compare(candidate.ref) == 0
                   && alt.alt_allele.compare(candidate.alt2) == 0
                   && alt.alt_type == candidate.alt2_type) {
                    if(vcf_rec.genotype[0] == 2 && vcf_rec.genotype[1] == 2) {
                        gt2 = Genotype::HOM_ALT;
                    }else if(vcf_rec.genotype[0] == 2 || vcf_rec.genotype[1] == 2) {
                        gt2 = Genotype::HET;
                    }
                }
            }
        }
        label_alt1.push_back(gt1);
        label_alt2.push_back(gt2);
    }
}


vector<PileupImage> ImageGenerator::create_window_pileups(vector<pair<long long, long long> > windows,
                                                          vector<type_read> reads, bool train_mode) {
    // container for all the images we will generate
    vector<PileupImage> pileup_images(windows.size());
    int inferred_window_size = windows[0].second - windows[0].first;

    // initialize all the pileup images with reference sequence
    for(int i=0; i<windows.size(); i++) {
        pileup_images[i].set_values(this->chromosome_name, windows[i].first, windows[i].second);
        string ref_seq = get_reference_sequence(windows[i].first, windows[i].second);
        pileup_images[i].image_alt1.insert(pileup_images[i].image_alt1.end(),
                                           PileupPixels::REF_ROW_BAND,
                                           get_reference_row(ref_seq));
        pileup_images[i].image_alt2.insert(pileup_images[i].image_alt2.end(),
                                           PileupPixels::REF_ROW_BAND,
                                           get_reference_row(ref_seq));
        if(train_mode) {
            set_window_labels(windows[i], pileup_images[i].label_alt1, pileup_images[i].label_alt2);
        }
    }

    // now iterate through each of the reads and add it to different windows if read overlaps
    for(int i=0; i<reads.size(); i++) {
        reads[i].set_read_id(i);
        vector<int> windows_indices;
        for(int j=0; j < windows.size(); j++) {
            if(overlap_length_between_ranges(make_pair(reads[i].pos, reads[i].pos_end), windows[j])){
                windows_indices.push_back(j);
            }
        }
        // if the read doesn't overlap with any of the window, then don't convert it
        if(windows_indices.empty()) continue;
        long long read_start, read_end;
        // convert the read to a pileup row
        vector<vector<int> > image_row_alt1, image_row_alt2;
        read_to_image_row(reads[i], read_start, read_end, image_row_alt1, image_row_alt2);
        for(auto&index: windows_indices) {
            if(read_end < windows[index].first || read_start > windows[index].second) continue;
            // assign the read to each of the window it overlaps with
            assign_read_to_window(pileup_images[index], image_row_alt1, image_row_alt2, read_start, read_end);
        }
    }
    vector<vector<int> > empty_row;
    empty_row.insert(empty_row.end(), inferred_window_size, {0, 0, 0, 0, 0});
    for(auto&pileup: pileup_images) {
        int empties_needed = PileupPixels::IMAGE_HEIGHT - pileup.image_alt1.size();
        if(empties_needed > 0){
            pileup.image_alt1.insert(pileup.image_alt1.end(), empties_needed, empty_row);
            pileup.image_alt2.insert(pileup.image_alt2.end(), empties_needed, empty_row);
        }
    }
    return pileup_images;
}