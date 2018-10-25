//
// Created by Kishwar Shafin on 10/24/18.
//

#ifndef FRIDAY_CANDIDATE_FINDER_H
#define FRIDAY_CANDIDATE_FINDER_H

#include <cmath>
#include "../dataio/bam_handler.h"
using namespace std;

namespace CandidateFinder_options {
    static constexpr int min_mapping_quality = 14;
    static constexpr int min_base_quality = 15;
    static constexpr int freq_threshold = 12;
};

namespace AlleleType {
    static constexpr int SNP_ALLELE = 1;
    static constexpr int INSERT_ALLELE = 2;
    static constexpr int DELETE_ALLELE = 3;
};


struct PositionalCandidateRecord{
    string chromosome_name;
    long long pos;
    long long pos_end;
    string ref;
    string alt1;
    string alt2;
    int alt1_type;
    int alt2_type;

    PositionalCandidateRecord() {
        this->alt1_type = 0;
        this->alt2_type = 0;
        this->alt1 = '.';
        this->alt2 = '.';
    }

    void set_positions(string chromosme_name, long long pos, long long pos_end) {
        this->chromosome_name = chromosme_name;
        this->pos = pos;
        this->pos_end = pos_end;
    }

    void set_reference(string ref) {
        this->ref = ref;
    }

    void set_alt1(string alt1, int alt1_type) {
        this->alt1 = alt1;
        this->alt1_type = alt1_type;
    }

    void set_alt2(string alt2, int alt2_type) {
        this->alt2 = alt2;
        this->alt2_type = alt2_type;
    }
};

struct CandidateAllele{
    string ref;
    string alt;
    int alt_type;

    CandidateAllele(string ref, string alt, int alt_type) {
        this->ref = ref;
        this->alt = alt;
        this->alt_type = alt_type;
    }

    CandidateAllele() {}
};

struct Candidate{
    long long pos;
    long long pos_end;
    CandidateAllele allele;

    Candidate(long long pos_start, long long pos_end, string ref, string alt, int alt_type) {
        this->pos = pos_start;
        this->pos_end = pos_end;
        this->allele.alt = alt;
        this->allele.ref = ref;
        this->allele.alt_type = alt_type;
    }
    bool operator< (const Candidate& that ) const {
        if(this->pos != that.pos) return this->pos < that.pos;
        if(this->pos_end != that.pos_end) return this->pos_end < that.pos_end;
        if(this->allele.alt != that.allele.alt) return this->allele.alt < that.allele.alt;
        if(this->allele.ref != that.allele.ref) return this->allele.ref < that.allele.ref;
        if(this->allele.alt_type != that.allele.alt_type) return this->allele.alt_type < that.allele.alt_type;
        return this->pos < that.pos;
    }

    bool operator==(const Candidate& that ) const {
        if(this->pos == that.pos &&
           this->pos_end == that.pos_end &&
           this->allele.ref == that.allele.ref &&
           this->allele.alt == that.allele.alt &&
           this->allele.alt_type == that.allele.alt_type)
            return true;
        return false;
    }

};

class CandidateFinder {
    long long region_start;
    long long region_end;
    string chromosome_name;
    string reference_sequence;
    map<Candidate, int> AlleleFrequencyMap;
    map< long long, vector<Candidate> > AlleleMap;
public:
    CandidateFinder(string reference_sequence,
                    string chromosome_name,
                    long long region_start,
                    long long region_end);
    void add_read_alleles(type_read &read, vector<int> &coverage);
    vector<PositionalCandidateRecord>  find_candidates(vector<type_read> reads);
};



#endif //FRIDAY_CANDIDATE_FINDER_H
