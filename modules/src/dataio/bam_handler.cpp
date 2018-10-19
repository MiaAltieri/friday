//
// Created by Kishwar Shafin on 6/14/18.
//
#include "../../headers/dataio/bam_handler.h"

BAM_handler::BAM_handler(string path) {
    this->hts_file = sam_open(path.c_str(), "r");
    // Used to check if valid bam file
    if(this->hts_file == NULL){
        cerr<<"INVALID BAM FILE. PLEASE CHECK IF PATH IS CORRECT: "<<path<<endl;
        exit (EXIT_FAILURE);
    }

    // Setting up the indexing
    this->idx = sam_index_load(this->hts_file, path.c_str());
    if(this->idx == NULL){
        cerr<<"INVALID BAM INDEX FILE. PLEASE CHECK IF FILE IS INDEXED: "<<path<<endl;
        exit (EXIT_FAILURE);
    }
}

set<string> BAM_handler::get_sample_names() {
    bam_hdr_t* header;
    header = sam_hdr_read(this->hts_file);
    if(header == 0){
        cerr<<"INVALID BAM FILE. PLEASE CHECK IF PATH IS CORRECT"<<endl;
        exit (EXIT_FAILURE);
    }

    int l_text = header->l_text;
    char *text = header->text;
    set<string> samples;
    stringstream full_header_tokenizer(text);
    string line;
    while(getline(full_header_tokenizer, line, '\n')){
        istringstream header_line_tokenizer(line);
        string token;
        getline(header_line_tokenizer, token, '\t');
        if(token == "@RG") {
            while(getline(header_line_tokenizer, token, '\t')){
                istringstream tag_tokenizer(token);
                string tag;
                getline(tag_tokenizer, tag, ':');
                if(tag == "SM"){
                    string sample_name;
                    getline(tag_tokenizer, sample_name, ':');
                    samples.insert(sample_name);
                }
            }
        }
    }
    bam_hdr_destroy(header);
    return samples;
}

type_read_flags BAM_handler::get_read_flags(int flag) {
    type_read_flags flags;
    if ( flag&BAM_FPAIRED ) flags.is_paired = 1;
    if ( flag&BAM_FPROPER_PAIR ) flags.is_proper_pair = 1;
    if ( flag&BAM_FUNMAP ) flags.is_unmapped = 1;
    if ( flag&BAM_FMUNMAP ) flags.is_mate_unmapped = 1;
    if ( flag&BAM_FREVERSE ) flags.is_reverse = 1;
    if ( flag&BAM_FMREVERSE ) flags.is_mate_is_reverse = 1;
    if ( flag&BAM_FREAD1 ) flags.is_read1 = 1;
    if ( flag&BAM_FREAD2 ) flags.is_read2 = 1;
    if ( flag&BAM_FSECONDARY ) flags.is_secondary = 1;
    if ( flag&BAM_FQCFAIL ) flags.is_qc_failed = 1;
    if ( flag&BAM_FDUP ) flags.is_duplicate = 1;
    if ( flag&BAM_FSUPPLEMENTARY ) flags.is_supplementary = 1;
    return flags;
}

vector<type_sequence> BAM_handler::get_chromosome_sequence_names_with_length() {
    // Setting up the header
    bam_hdr_t* header;
    header = sam_hdr_read(this->hts_file);
    if(header == 0){
        cerr<<"INVALID BAM FILE. PLEASE CHECK IF PATH IS CORRECT"<<endl;
        exit (EXIT_FAILURE);
    }


    // Get all the sequence names. These are the chromosome names from the bed header file.
    vector<type_sequence> sequence_names;
    int total_targets = header->n_targets;
    for (int i=0; i < total_targets; i++){
        type_sequence sq_info;
        sq_info.sequence_name = header->target_name[i];
        sq_info.sequence_length = header->target_len[i];
        sequence_names.push_back(sq_info);
    }

    bam_hdr_destroy(header);

    return sequence_names;
}

vector<string> BAM_handler::get_chromosome_sequence_names() {
    // Setting up the header
    bam_hdr_t* header;
    header = sam_hdr_read(this->hts_file);
    if(header == 0){
        cerr<<"INVALID BAM FILE. PLEASE CHECK IF PATH IS CORRECT"<<endl;
        exit (EXIT_FAILURE);
    }
    // Get all the sequence names. These are the chromosome names from the bed header file.
    vector<string> sequence_names;
    int total_targets = header->n_targets;
    for (int i=0; i < total_targets; i++){
        string seq_name = header->target_name[i];
        sequence_names.push_back(seq_name);
    }

    return sequence_names;
}

vector<type_read> BAM_handler::get_reads(string chromosome, long long start, long long stop) {
    vector <type_read> reads;

    // Setting up the header
    bam_hdr_t* header;
    header = sam_hdr_read(this->hts_file);
    if(header == 0){
        cerr<<"INVALID BAM FILE. PLEASE CHECK IF PATH IS CORRECT"<<endl;
        exit (EXIT_FAILURE);
    }

    // get the id of the chromosome
    const int tid = bam_name2id(header, chromosome.c_str());

    // get the iterator
    hts_itr_t *iter  = sam_itr_queryi(this->idx, tid, start, stop);

    // initialize an alignment
    bam1_t* alignment = bam_init1();

    while(sam_itr_next(this->hts_file, iter, alignment) >= 0) {
        //get read flags
        type_read_flags read_flags = get_read_flags(alignment->core.flag);
        if(read_flags.is_supplementary || read_flags.is_qc_failed || read_flags.is_duplicate || read_flags.is_secondary
           || read_flags.is_unmapped){
            continue;
        }

        // mapping quality
        if(alignment->core.qual <= 0){
            continue;
        }

        // get query name
        string query_name = bam_get_qname(alignment);

        // get the position where it gets aligned
        int32_t pos = alignment->core.pos;
        uint32_t len = alignment->core.l_qseq;

        // get the base qualities and sequence bases
        uint8_t *seqi = bam_get_seq(alignment);
        uint8_t *qual = bam_get_qual(alignment);

        vector<int> base_qualities;
        string read_seq;

        for (int i = 0; i < len; i++) {
            int base_quality = (int) qual[i];
            base_qualities.push_back(base_quality);

            read_seq += seq_nt16_str[bam_seqi(seqi, i)];
        }

        // get the cigar operations of the alignment
        uint32_t *cigar = bam_get_cigar(alignment);
        string str_cigar;
        vector <type_cigar> cigar_tuples;

        for(int k = 0; k < alignment->core.n_cigar; k++) {
            int cigar_op = bam_cigar_op(cigar[k]);
            int cigar_len = bam_cigar_oplen(cigar[k]);

            type_cigar cigar_instance;

            cigar_instance.cigar_op = cigar_op;
            cigar_instance.cigar_len = cigar_len;

            cigar_tuples.push_back(cigar_instance);
        }

        // mapping quality
        int map_quality = alignment->core.qual;

        // set all fetched attributes
        type_read read;
        read.query_name = query_name;
        read.pos = pos;
        read.sequence = read_seq;
        read.flags = read_flags;
        read.mapping_quality = map_quality;
        read.base_qualities = base_qualities;
        read.cigar_tuples = cigar_tuples;
        reads.push_back(read);
    }
    bam_destroy1(alignment);
    hts_itr_destroy(iter);
    bam_hdr_destroy(header);


    return reads;
}

BAM_handler::~BAM_handler() {
    // destroy everything
    hts_idx_destroy( this->idx );
    sam_close( this->hts_file );
}