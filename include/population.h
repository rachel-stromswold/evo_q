#ifndef POPULATION_H
#define POPULATION_H

#include "organism.h"
#include <vector>
#include <fstream>
#include <memory>

#define NUM_GENES	10
#define NUM_CARRY	1

#define MAX_NUM_GENS	100
#define OUT_BUF_SIZE    10

#define DEF_SORT_PARAM	-3

#define FLAG_NONE_SET	0
#define FLAG_STATS_SET	1
#define FLAG_DIST_SET	2
#define FLAG_BEST_FOUND	4
#define FLAG_FRONTS	8
#define VALID_BEST	16

#if __cplusplus >= 201402L
#define DEPRECATED(msg) [[ deprecated(msg) ]]
#elif defined(__GNUC__)
#define DEPRECATED(msg) __attribute__ ((deprecated(msg)))
#elif defined(_MSC_VER)
#define DEPRECATED(msg) __declspec(deprecated(msg))
#else
#define DEPRECATED(msg) 
#endif

namespace Genetics {

struct FitnessStats {
  double mean;
  double var;
  double max;
  double min;
};

class ConvergenceCriteria {
  public:
    virtual bool evaluate_convergence(_uint N_OBJS, FitnessStats* stats) = 0;
};

class Population {
  private:
    _uint N_BITS;
    _uint N_PARAMS;
    _uint N_OBJS;
    _uint generation = 0;
    _uchar calculated_flags = 0;
  protected:
    size_t sort_org_calls = 0;
    size_t carryover_num;//How many of the best individuals carry over to the next generation 

    //OWNED POINTERS
    FitnessStats* pop_stats = NULL;
    //EXTERNALLY MANAGED POINTERS
    std::shared_ptr<PhenotypeMap> map;
    
    ArgStore args;
    //all offspring from the previous generation
    size_t offspring_num;
    std::vector<std::shared_ptr<Organism>> offspring;
    std::vector<std::shared_ptr<Organism>> old_gen;
    size_t min_penalty_ind, max_penalty_ind;
    //which offspring will survive to enter the next breeding round
    size_t survivors_num;
    std::vector<std::shared_ptr<Organism>> survivors;
    //guarantee that the best organism appears in the next generation
//    size_t best_organism_ind;
    Organism best_organism;
    //labels for generating data output
//    Vector<String> var_labels;
//    Vector<String> obj_labels;
    char** var_labels;
    char** obj_labels;

    //cull in place is slightly faster but less accurate than the standard cull method
    bool cull_in_place();
    //cull first sorts the organisms and selects them based on the ratio of their relative fitness to the total relative fitness
    bool cull(); 
    void breed_shuffle();
    void breed();
    void tournament_selection();
    void find_best_organism();
    void calculate_distances();
    void hypermutate();
    void set_best_organism(_uint i);

  public:
    Population(_uint pn_bits, _uint pn_objs, PhenotypeMap* p_map, String conf_fname = "");
    Population(_uint pn_bits, _uint pn_objs, Organism* tmplt, PhenotypeMap* p_map, String conf_fname = "");
    Population(_uint pn_bits, _uint pn_objs, std::shared_ptr<PhenotypeMap> p_map, String conf_fname = "");
    Population(_uint pn_bits, _uint pn_objs, Organism* tmplt, std::shared_ptr<PhenotypeMap> p_map, String conf_fname = "");
    ~Population();
    Population(Population& o);
    Population& operator=(Population& o);
    Population(Population&& o);

    void set_convergence_type(ConvergenceCriteria* conv);
    void resize_population(_uint new_size);
    void set_n_survivors(_uint new_size);
    void evaluate(Problem* prob);
    void evaluate_async(Problem* prob);
    bool iterate(ConvergenceCriteria* conv = NULL);

    std::shared_ptr<Organism> get_best_organism(size_t i = 0);
    std::shared_ptr<Organism> get_organism(size_t i);
    std::shared_ptr<Organism> get_child(size_t i);

    void run(Problem* prob);
    void swap_orgs(int i, int j);
    int partition(_uint ind, std::vector<std::shared_ptr<Organism>>* work_arr, int s, int e);
    void sort_orgs(unsigned int ind, std::vector<std::shared_ptr<Organism>>* arr, int s = DEF_SORT_PARAM, int e = DEF_SORT_PARAM);
    Vector<String> get_header();
    Vector<String> get_pop_data();

    FitnessStats get_pop_stats(_uint i = 0) { return pop_stats[i]; }
    DEPRECATED("get_min_fitness is deprecated, use get_pop_stats instead") double get_min_fitness(_uint i = 0) {
      return pop_stats[i].min;
    }
    DEPRECATED("get_max_fitness is deprecated, use get_pop_stats instead") double get_max_fitness(_uint i = 0) {
      return pop_stats[i].max;
    }

    void set_var_label(_uint ind, String val); 
    void set_obj_label(_uint ind, String val);

    size_t get_offspring_num() { return offspring_num; }
    size_t get_survivors_num() { return survivors_num; }
    _uint get_n_bits() { return N_BITS; }
    _uint get_n_objs() { return N_OBJS; }
    ArgStore& get_args() { return args; }
};

class Population_NSGAII : public Population {
  private:
    _uchar calculated_flags = FLAG_NONE_SET;
    size_t survivors_num;
    //the ngsa alternative to elitism
    std::vector<std::vector<std::shared_ptr<Organism>>> pareto_fronts;
    void hypermutate();
    void update_fitness_ranges(Organism* org, size_t i);
    _uint generation = 0;
    void make_fronts(std::vector<std::shared_ptr<Organism>>* cmb_arr);

  public:
    Population_NSGAII(_uint pn_bits, _uint pn_objs, PhenotypeMap* p_map, String conf_fname = "");
    Population_NSGAII(_uint pn_bits, _uint pn_objs, Organism* tmplt, PhenotypeMap* p_map, String conf_fname = "");
    Population_NSGAII(_uint pn_bits, _uint pn_objs, std::shared_ptr<PhenotypeMap> p_map, String conf_fname = "");
    Population_NSGAII(_uint pn_bits, _uint pn_objs, Organism* tmplt, std::shared_ptr<PhenotypeMap> p_map, String conf_fname = "");
    ~Population_NSGAII();

    void evaluate(Problem* prob);
    std::shared_ptr<PhenotypeMap> get_map() { return this->map; }

    //cull in place is slightly faster but less accurate than the standard cull method
    bool cull_in_place();
    //cull first sorts the organisms and selects them based on the ratio of their relative fitness to the total relative fitness
    void cull();
    void breed();
    void calculate_pop_stats();
    bool iterate(ConvergenceCriteria* conv = NULL);

    _uint get_n_pareto_fronts() {
      return pareto_fronts.size();
    }
    std::vector<std::shared_ptr<Organism>> get_pareto_front(_uint i);

    Vector<String> get_header();
    Vector<String> get_pop_data();

    std::shared_ptr<Organism> get_organism(size_t i) { return (i < offspring_num) ? this->old_gen[i] : this->offspring[i - offspring_num]; }

    size_t get_offspring_num() { return this->offspring_num; }
    size_t get_survivors_num() { return this->survivors_num; }
};

}

#endif //POPULATION_H
