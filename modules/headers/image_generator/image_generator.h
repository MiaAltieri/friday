//
// Created by Kishwar Shafin on 11/1/18.
//

#ifndef FRIDAY_IMAGE_GENERATOR_H
#define FRIDAY_IMAGE_GENERATOR_H

#include <algorithm>
using namespace std;
#include "../dataio/bam_handler.h"

class ImageGenerator {
public:
    void generate_image(pair<long long, long long> windows, vector<type_read> reads, string ref);
};


#endif //FRIDAY_IMAGE_GENERATOR_H
