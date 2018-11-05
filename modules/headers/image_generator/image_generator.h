//
// Created by Kishwar Shafin on 11/1/18.
//

#ifndef FRIDAY_IMAGE_GENERATOR_H
#define FRIDAY_IMAGE_GENERATOR_H

#include <algorithm>
using namespace std;
#include "../dataio/bam_handler.h"

namespace PileupPixels {
    static constexpr int MAX_COLOR_VALUE = 254;
    static constexpr int BASE_QUALITY_CAP = 40;
    static constexpr int MAP_QUALITY_CAP = 60;
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
    pair<vector< vector<uint8_t> >, pair<long long, long long> > read_to_image_row(type_read read);
};


#endif //FRIDAY_IMAGE_GENERATOR_H
