//
// Created by Kishwar Shafin on 11/1/18.
//

#include "../../headers/image_generator/image_generator.h"

void ImageGenerator::generate_image(pair<long long, long long> window, vector<type_read> reads, string ref) {
        cout<<window.first<<" "<<window.second<<endl;
        cout<<ref<<endl;
        long long window_length = window.second - window.first;
        for(auto &read: reads){
            if(read.pos > window.second or read.pos_end < window.first) continue;


            long long read_start_index = max((long long) 0, window.first - read.pos);
            long long read_length = window_length;

            if(read.pos_end < window.second) read_length = read.pos_end - max(read.pos, window.first);
            if(read.pos > window.first)
                for(int i=window.first; i<read.pos; i++)cout<<" ";
//            cout<<window.first<<" "<<window.second<<endl;
//            cout<<read.pos<<" "<<read.pos_end<<endl;
//            cout<<read_start_index<<" "<<read_length<<endl;
            string read_seq = read.sequence.substr(read_start_index, read_length);
            cout<<read_seq<<endl;
        }
}