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

namespace Genetics {

struct FitnessStats {
  double mean;
  double var;
  double max;
  double min;
};

class ConvergenceCriteria {
  public:
    virtual bool evaluate_convergence(FitnessStats* stats) = 0;
};

class Population {
  private:
    _uint N_BITS;
    _uint N_PARAMS;
    _uint N_OBJS;
    _uint generation = 0;
  protected:
    size_t sort_org_calls = 0;
    size_t carryover_num;//How many of the best individuals carry over to the next generation 

    //OWNED POINTERS
    FitnessStats* pop_stats = NULL;
    //EXTERNALLY MANAGED POINTERS
    PhenotypeMap* map = NULL;
    
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
    size_t best_organism_ind;
    std::shared_ptr<Organism> best_organism;
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
    void find_best_organism();

  public:
    Population(_uint pn_bits, _uint pn_objs, PhenotypeMap* p_map, String conf_fname = "");
    Population(_uint pn_bits, _uint pn_objs, Organism* tmplt, PhenotypeMap* p_map, String conf_fname = "");
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

    double get_min_fitness(_uint i = 0) { return pop_stats[i].min; }
    double get_max_fitness(_uint i = 0) { return pop_stats[i].max; }

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
    ~Population_NSGAII();

    void evaluate(Problem* prob);
    PhenotypeMap* get_map() { return this->map; }

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
    std::vector<std::shared_ptr<Organism>> get_pareto_front(_uint i) {
      if (pareto_fronts.size() == 0) { error(1, "The population has not yet been evaluated."); }
      return pareto_fronts[i];
    }

    Vector<String> get_header();
    Vector<String> get_pop_data();

    std::shared_ptr<Organism> get_organism(size_t i) { return (i < offspring_num) ? this->old_gen[i] : this->offspring[i - offspring_num]; }

    size_t get_offspring_num() { return this->offspring_num; }
    size_t get_survivors_num() { return this->survivors_num; }
};

}

#endif //POPULATION_H