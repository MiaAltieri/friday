//
// Created by Kishwar Shafin on 11/1/18.
//

#ifndef FRIDAY_IMAGE_GENERATOR_H
#define FRIDAY_IMAGE_GENERATOR_H

#include <algorithm>
using namespace std;
#include "../dataio/bam_handler.h"
#include "../candidate_finding/candidate_finder.h"

namespace PileupPixels {
    static constexpr int MAX_COLOR_VALUE = 254;
    static constexpr int BASE_QUALITY_CAP = 40;
    static constexpr int MAP_QUALITY_CAP = 60;
    static constexpr int REF_ROW_BAND = 5;
    static constexpr int IMAGE_HEIGHT = 100;
};

struct PileupImage {
    string chromosome_name;
    long long start_pos;
    long long end_pos;
    vector<vector<vector<uint8_t> > >image;

    void set_values(string chromosome_name, long long start_pos, long long end_pos) {
        this->chromosome_name = chromosome_name;
        this->start_pos = start_pos;
        this->end_pos = end_pos;
    }
};

class ImageGenerator {
    long long ref_start;
    long long ref_end;
    string chromosome_name;
    string reference_sequence;
    map<long long, PositionalCandidateRecord> all_positional_candidates;
    map<char, uint8_t> global_base_color;
public:
    ImageGenerator(string reference_sequence,
                   string chromosome_name,
                   long long ref_start,
                   long long ref_end,
                   map<long long, PositionalCandidateRecord> all_positional_candidates);

    string get_reference_sequence(long long st_pos, long long end_pos);
    vector<vector<uint8_t> > read_to_image_row(type_read read, long long &read_start, long long &read_end);
    vector<vector<uint8_t> > get_reference_row(string ref_seq);
    vector<uint8_t> get_window_labels(pair<long long, long long> window,
                                      map<long long, vector<type_positional_vcf_record> > pos_vcf);
    vector<PileupImage> create_window_pileups(vector<pair<long long, long long> > windows,
                                              vector<type_read> reads);
    int get_which_allele(long long pos, string ref, string alt, int alt_type);
    long long overlap_length_between_ranges(pair<long long, long long> range_a,
                                            pair<long long, long long> range_b);
    void assign_read_to_window(PileupImage& pileup_image,
                               vector<vector<uint8_t> >& image_row,
                               long long read_start,
                               long long read_end);
};


#endif //FRIDAY_IMAGE_GENERATOR_H
