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
    static constexpr int IMAGE_HEIGHT = 250;
    static constexpr int CONTEXT_SIZE = 125;
};


namespace Genotypes {
    static constexpr int HOM = 0;
    static constexpr int HET = 1;
    static constexpr int HOM_ALT = 2;
};

struct PileupImage {
    string chromosome_name;
    long long start_pos;
    long long end_pos;
    vector<vector<vector<uint8_t> > >image;
    int label;
    string name;
    PileupImage(string name, int label, vector<vector<vector<uint8_t> > > image) {
        this->name = name;
        this->label = label;
        this->image = image;
    }
    PileupImage() {
    }
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
    map<char, uint8_t> global_base_color;
    map< string,  vector<vector<uint8_t> > > read_row_cache;
    map<string, long long> read_start_cache;
    map<string, long long> read_end_cache;
public:
    ImageGenerator(string reference_sequence,
                   string chromosome_name,
                   long long ref_start,
                   long long ref_end);

    string get_reference_sequence(long long st_pos, long long end_pos);
    vector<vector<uint8_t> > read_to_image_row(type_read &read,
                                           long long &read_start,
                                           long long &read_end);
    vector<vector<uint8_t> > get_reference_row(string ref_seq, int left_pad, int right_pad);
    PileupImage create_image(PositionalCandidateRecord candidate,
                             vector<type_read> &reads,
                             vector< pair<int, bool> > read_index_and_support,
                             int genotype);
    void assign_read_to_image(PileupImage& pileup_image,
                              vector<vector<uint8_t> >& image_row,
                              long long read_start,
                              long long read_end,
                              int left_pad,
                              int right_pad);
    long long overlap_length_between_ranges(pair<long long, long long> range_a,
                                            pair<long long, long long> range_b);
};


#endif //FRIDAY_IMAGE_GENERATOR_H
