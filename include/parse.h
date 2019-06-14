#ifndef PARSE_H
#define PARSE_H

#include <fstream>
#include "util.h"

#define DEF_POP_SIZE		    20
#define DEF_BREED_POP_SIZE	6
#define DEF_NUM_GENS		    10
#define DEF_COUPLING_BITS	  16
#define DEF_MAX_COUPLING	  100
#define DEF_COUP_VAR		    0.25
#define DEF_COUP_MEAN		    0.0
#define DEF_MUTATE_PROB		  0.1
#define DEF_HYPER_THRESH    0.8
#define DEF_NUM_CROSSOVERS  1

//flags
#define WAIT_CON		1
#define VERBOSE			2

#define BUF_SIZE		50

namespace Genetics {

class ArgStore {
  private:
    std::mt19937 generator;
//REMOVAL CANDIDATE
    std::binomial_distribution<unsigned char>* long_bin;
    std::binomial_distribution<unsigned char>* short_bin;
    std::bernoulli_distribution* bern;
    unsigned int short_bin_n = 0;
//END REMOVAL CANDIDATE

    ArgStore(ArgStore const&){};
    ArgStore& operator=(ArgStore const& arg) {return *this;};

    unsigned flags;
    size_t pop_size;
    size_t breed_pop_size;
    size_t num_gens;
    int num_crossovers;
    double init_coup_var;
    double init_coup_mean;
    double mutate_prob;
    double hypermutation_threshold;
    std::string out_fname;

    bool activate = true;
    int seed = 0;

  public: 
    ArgStore();

    void initialize_from_args(size_t argc, char** argv);
    void initialize_from_file(const char* fname);
    void initialize();
    void read_file(std::string fname);
    std::mt19937& get_generator() { return generator; }
//REMOVAL CANDIDATE
    unsigned int sample_binomial(unsigned int n);
    bool sample_bernoulli();
    void print_data();
//END REMOVAL CANDIDATE

    size_t get_pop_size()	{ return pop_size; }
    size_t get_survivors()	{ return breed_pop_size; }
    size_t get_num_gens() 	{ return num_gens; }
    int get_num_crossovers() 	{ return num_crossovers; }
    double get_init_coup_var()	{ return init_coup_var; }
    double get_init_coup_mean() { return init_coup_mean; }
    double get_mutate_prob()	{ return mutate_prob; }
    double get_hypermutation_threshold()	{ return hypermutation_threshold; }
    bool wait_for_con()		{ return flags & WAIT_CON; }
    bool verbose()		{ return flags & VERBOSE; }
    std::string get_out_fname() { return out_fname; }
};
    
}

#endif //PARSE_H
