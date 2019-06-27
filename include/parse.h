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
#define DEF_CROSSOVER_PROB	0.8
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
    std::bernoulli_distribution* bern_mut;
    std::bernoulli_distribution* bern_cross;
    unsigned int short_bin_n = 0;
//END REMOVAL CANDIDATE

    unsigned flags;
    size_t pop_size;
    size_t breed_pop_size;
    size_t num_gens;
    int num_crossovers;
    double init_coup_var;
    double init_coup_mean;
    double crossover_prob;
    double mutate_prob;
    double hypermutation_threshold;
    std::string out_fname;

    bool activate = true;
    int seed = 0;

  public: 
    ArgStore();
    ArgStore(const ArgStore& o);
    ArgStore(ArgStore&& o);
    ~ArgStore();

    void initialize_from_args(size_t argc, char** argv);
    void initialize_from_file(const char* fname);
    void initialize();
    void read_file(std::string fname);
    std::mt19937& get_generator() { return generator; }
//REMOVAL CANDIDATE
    unsigned int sample_binomial(unsigned int n);
    bool random_mutation();
    bool random_crossover();
    void print_data();
//END REMOVAL CANDIDATE

    size_t get_pop_size()			{ return pop_size; }
    void set_pop_size(size_t n)			{ pop_size = n; }
    size_t get_survivors()			{ return breed_pop_size; }
    void set_survivors(size_t n)		{ breed_pop_size = n; }
    size_t get_num_gens() 			{ return num_gens; }
    void set_num_gens(size_t n) 		{ num_gens = n; }
    int get_num_crossovers() 			{ return num_crossovers; }
    void set_num_crossovers(int n) 		{ num_crossovers = n; }
    double get_init_coup_var()			{ return init_coup_var; }
    void set_init_coup_var(double x)		{ init_coup_var = x; }
    double get_init_coup_mean() 		{ return init_coup_mean; }
    void set_init_coup_mean(double x) 		{ init_coup_mean = x; }
    double get_crossover_prob()			{ return crossover_prob; }
    void set_crossover_prob(double x)		{ crossover_prob = x; }
    double get_mutate_prob()			{ return mutate_prob; }
    void set_mutate_prob(double x)		{ mutate_prob = x; }
    double get_hypermutation_threshold()	{ return hypermutation_threshold; }
    void set_hypermutation_threshold(double x)	{ hypermutation_threshold = x; }
    bool wait_for_con()				{ return flags & WAIT_CON; }
    void set_wait_for_con(bool b = true)	{ flags |= WAIT_CON * ((_uint)b); }
    bool verbose()				{ return flags & VERBOSE; }
    void set_verbose(bool b = true)		{ flags |= VERBOSE * ((_uint)b); }
    std::string get_out_fname() 		{ return out_fname; }
};
    
}

#endif //PARSE_H
