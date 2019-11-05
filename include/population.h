#ifndef POPULATION_H
#define POPULATION_H

#include "compare.h"
#include <vector>
#include <fstream>
#include <memory>

//check for the presence of libomp
#ifdef USE_LIBOMP
#include <omp.h>
#endif

#define NUM_GENES	10
#define NUM_CARRY	1

#define MAX_NUM_GENS	100
#define OUT_BUF_SIZE    50

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
    virtual ~ConvergenceCriteria() = default;
};

template <class FitType, class SelectType, class=void>
class Population {
private:
  _uint N_BITS;
  _uint N_PARAMS;
  _uint N_OBJS;
  _uint generation = 0;
  _uchar calculated_flags = 0;
  void evaluate_best(Problem<FitType>* prob, double forget_weight=0.0) {
    best_organism->evaluate_fitness(prob);
    for (_uint i = 0; i < N_OBJS; ++i) {
      pop_stats[i].max = best_organism->get_fitness(i);
    }
  }

protected:
  static_assert( std::is_base_of<Selector<FitType>, SelectType>::value, "SelectType must be derived from Selector<FitType>" );
  SelectType sel;
  size_t sort_org_calls = 0;
  size_t carryover_num;//How many of the best individuals carry over to the next generation 

  //OWNED POINTERS
  FitnessStats* pop_stats = NULL;
  //EXTERNALLY MANAGED POINTERS
  std::shared_ptr<PhenotypeMap> map;
  
  ArgStore args;
  //all offspring from the previous generation
  size_t offspring_num;
  std::vector<std::shared_ptr<Organism<FitType>>> offspring;
  std::vector<std::shared_ptr<Organism<FitType>>> old_gen;
  size_t min_penalty_ind, max_penalty_ind;
  //which offspring will survive to enter the next breeding round
  size_t survivors_num;
  std::vector<std::shared_ptr<Organism<FitType>>> survivors;
  //guarantee that the best organism appears in the next generation
  _uint best_organism_ind = 0;
  std::shared_ptr<Organism<FitType>> best_organism;
  std::shared_ptr<Organism<FitType>> alltime_best_organism;
  //labels for generating data output
  char** var_labels;
  char** obj_labels;
  std::vector<bool> is_obj_cost;
  int print_penalties = 0;

  //cull in place is slightly faster but less accurate than the standard cull method
  bool cull_in_place() {
    size_t j = 0;

    double difference = pop_stats[0].max - pop_stats[0].min;
    //avoid divide by 0
    if (difference == 0) {
      error(CODE_WARN, "All organisms have the same fitness, exiting");
      return true;
    }
    std::uniform_real_distribution<double> dist(0, difference);

    if (survivors.size() < survivors_num) {
      survivors.resize(survivors_num);
    }
    for (size_t i = 0; i < this->offspring_num && j < survivors_num; ++i) {
      /* if M_f is the maximum fitness and m_f is the minimum the minimum, while x is the
       * fitness of a given organism, then the probability of survival is x/(M_f-m_f) or 1
       * if M_f = m_f */
      if (dist(args.get_generator()) < old_gen[i]->get_fitness(0) - pop_stats[0].min) {
        survivors[j] = old_gen[i];
        j++;
      }
    }

    std::uniform_int_distribution<int> ind_dist(0, this->offspring_num - 1);
    while (j < survivors_num) {
      survivors[j] = old_gen[ind_dist( args.get_generator() )];
      j++;
    }
    return false;
  }

  /**
   * \brief An implementation of simple roulette selection. This function first sorts the organisms and selects them based on the ratio of their relative fitness to the total relative fitness.
   *
   * \returns True if all organisms have the same fitness (results have converged).
   */
  bool cull() {
    sort_orgs(0, &old_gen);
    //if all the organisms have the same fitness then reinitialize the population
    if (old_gen[0]->get_fitness(0) - old_gen[this->offspring_num-1]->get_fitness(0) < 0.001 ) {
      return true;
    }
    double min_fit = old_gen[this->offspring_num-1]->get_fitness(0);
    double total_fit = 0;
    for (size_t i = 0; i < this->offspring_num; ++i) {
      total_fit += old_gen[i]->get_fitness(0) - min_fit;
    }
    std::uniform_real_distribution<double> dist(0, total_fit);

    //maintain a list of organisms that have already been added so no organism appears twice
    size_t* banned = (size_t*)malloc(sizeof(size_t)*survivors_num);
    for (size_t i = 0; i < survivors_num; ++i) { banned[i] = -1; }

    for (size_t i = 0; i < survivors_num; ++i) {
      double val = dist( args.get_generator() );
      size_t j = 0;
      while (val > 0.0) {
        val -= old_gen[j]->get_fitness(0) - min_fit;
        j++;
      }
      j -= 1;
      while (contains<size_t>(banned, survivors_num, j)) {
        j++;
        if (j >= survivors_num) {
          j = 0;
        }
      }
      survivors[i] = old_gen[j];
      banned[i] = j;
    }
    if (args.verbose()) {
      std::cerr << "Added orgs to survs:";
      for (size_t i = 0; i < survivors_num; ++i) {
        std::cerr << " " << banned[i];
      }
      std::cerr << std::endl;
    }
    free(banned);
  //  best_organism_ind = 0;
    return false;
  }
  void breed_shuffle() {
    std::uniform_int_distribution<size_t> dist_surv0(0, survivors_num - 1);
    std::vector<Organism<FitType>*> children;
    Organism<FitType>** shuffled_inds = (Organism<FitType>**)malloc(sizeof(Organism<FitType>*)*offspring_num);
    for (size_t i = 0; i < survivors_num; ++i) {
      shuffled_inds[i] = survivors[i].get();
    }
    for (size_t i = 0; i < survivors_num; ++i) {
      size_t ind_o = dist_surv0( args.get_generator() );
      Organism<FitType>* tmp = shuffled_inds[i];
      shuffled_inds[i] = shuffled_inds[ind_o];
      shuffled_inds[ind_o] = tmp;
    }
    size_t last_org_ind = offspring_num;
    for (size_t i = survivors_num; i < offspring_num; ++i) {
      size_t org_ind= dist_surv0( args.get_generator() );
      //ensure that we don't see the same organism breeding with itself
      if (org_ind == last_org_ind) {
        org_ind = (org_ind + 1) % offspring_num;
      }
      shuffled_inds[i] = survivors[org_ind].get();
      last_org_ind = org_ind;
    }
    if (this->offspring_num % 2 == 1) {
      //elitist algorithm, make the first individual in the next generation the previous most fit
      //offspring[0] = std::make_shared<Organism<FitType>>(best_organism);
      offspring[0] = best_organism;

      for (size_t i = 1; 2*i < this->offspring_num; i++) {
        children = shuffled_inds[2*i - 1]->breed(&args, shuffled_inds[2*i]);
        offspring[2*i] = std::shared_ptr<Organism<FitType>>(children[0]);
        offspring[2*i - 1] = std::shared_ptr<Organism<FitType>>(children[1]);
      }
    } else {
      for (size_t i = 0; 2*i + 1 < this->offspring_num; i++) {
        children = shuffled_inds[2*i]->breed(&args, shuffled_inds[2*i + 1]);
        offspring[2*i] = std::shared_ptr<Organism<FitType>>(children[0]);
        offspring[2*i + 1] = std::shared_ptr<Organism<FitType>>(children[1]);
      }
    }

    offspring.swap(old_gen);
  }
  /**
   * \brief Apply soft penalties to organisms, this gaurantees that a penalized organism will always have a lower fitness than an unpenalized organism
   *
   * \note Penalized organisms all have a nonzero penalty weight. Let this penalty weight be given by $P$. Let $m$ and $M$ be the minimum and maximum fitness among all unpenalized organisms respectively. Given an organism $O$ with an unpenalized fitness $F(O)$, after penalties are applied,
   * \f[
   * F_new(O)=m - F(O) - max(|M|, |m|)*P
   * \f]
   * By default $F(O)$ is set to zero unless the user specified function evaluate_fitness"("O")" calls O.update"()". 
   * Users may or may not want to assign fitnesses to penalized organisms depending on whether such a fitness is well defined.
   */
  void apply_penalties(Problem<FitType>* prob) {
    //calculate penalties based on the range of fitnesses
    //double penalty_fact = min(abs(pop_stats[0].max - pop_stats[0].min), 1);
    double penalty_fact = abs(pop_stats[0].max - pop_stats[0].min);
    if (penalty_fact == 0) { penalty_fact = 1; }
    if ( pop_stats[0].max == std::numeric_limits<double>::infinity() ) {
      penalty_fact = abs(pop_stats[0].min);
    } else if ( pop_stats[0].min == -std::numeric_limits<double>::infinity() ) {
      penalty_fact = abs(pop_stats[0].max);
    }
    bool penalties_applied = false;
#ifdef USE_LIBOMP
#pragma omp parallel for
    for (size_t i = 0; i < this->offspring_num; ++i) {
      if (old_gen[i]->penalized()) {
        double new_fit = pop_stats[0].min - max(pop_stats[0].max - old_gen[i]->get_fitness(0), 0);
        if (old_gen[i]->get_fitness > pop_stats[0].max) {
          new_fit = pop_stats[0].min
        }
        new_fit -= penalty_fact*old_gen[i]->get_penalty();
        old_gen[i]->set_fitness(new_fit);
        penalties_applied = true;
      }
    }
#else
    for (size_t i = 0; i < this->offspring_num; ++i) {
      if (old_gen[i]->penalized()) {
        double new_fit = pop_stats[0].min - max(pop_stats[0].max - old_gen[i]->get_fitness(0), 0);
        new_fit -= penalty_fact*old_gen[i]->get_penalty();
        old_gen[i]->set_fitness(new_fit);
        penalties_applied = true;
      }
    }
#endif
    if (!penalties_applied) {
      calculated_flags |= FLAG_STATS_SET | FLAG_BEST_FOUND;
    } else {
      calculated_flags = FLAG_NONE_SET;
    }

    if (best_organism->get_fitness(0) > alltime_best_organism->get_fitness(0)) { 
      for (_uint j = 0; j < args.noise_compensate(); ++j) {
        evaluate_best(prob, args.forget_weight);
      }
      //alltime_best_organism = best_organism->copy();
      alltime_best_organism = best_organism;
      //alltime_best_organism->set_fitness(0, best_organism->get_fitness(0));
    }
    //check to see if there has been a decrease in fitness
    if ( args.noise_compensate() &&
         alltime_best_organism->valid() && 
         alltime_best_organism != best_organism &&
         alltime_best_organism->get_fitness(0) > best_organism->get_fitness(0) ) {
      //run more evaluations to see if the new organism is actually better
      for (_uint i = 0; i < args.noise_compensate(); ++i) {
        alltime_best_organism->evaluate_fitness(prob);
        evaluate_best(prob, args.forget_weight);
      }
      if ( alltime_best_organism->get_fitness(0) > best_organism->get_fitness(0) ) {
        //best_organism->swap(alltime_best_organism);
        best_organism = alltime_best_organism;
        pop_stats[0].max = best_organism->get_fitness(0);
      }
    }
  }
  void find_best_organism() {
    for (int j = 0; j < N_OBJS; ++j) {
      pop_stats[j].max = best_organism->get_fitness(j);
      pop_stats[j].min = best_organism->get_fitness(j);
      pop_stats[j].mean = best_organism->get_fitness(j) / offspring_num;
      for (size_t i = 0; i < offspring_num; ++i) {
        double fitness_i = old_gen[i]->get_fitness(j);
        if (fitness_i > pop_stats[j].max) {
          //TODO: make this usefully track multiple objectives
          if (j == 0) {
            set_best_organism(i);
          }
        }
        if (fitness_i < pop_stats[j].min) {
    pop_stats[j].min = fitness_i;
        }
        pop_stats[j].mean += fitness_i / offspring_num;
      }
      //calculate the variance
      pop_stats[j].var = 0;
      for (size_t i = 0; i < offspring_num; ++i) {
        double fitness_i = old_gen[i]->get_fitness(j);
        pop_stats[j].var += (fitness_i - pop_stats[j].mean)*(fitness_i - pop_stats[j].mean);
      }
    }
    
    calculated_flags |= FLAG_STATS_SET | FLAG_BEST_FOUND;
  }
  void breed() {
    find_best_organism();
    std::uniform_int_distribution<size_t> dist_surv0(0, survivors_num - 1);
    std::uniform_int_distribution<size_t> dist_surv1(0, survivors_num - 2);
    std::vector<Organism<FitType>*> children;
    size_t* shuffled_inds = (size_t*)malloc(sizeof(size_t)*survivors_num);
    if (this->offspring_num % 2 == 1) {
      //offspring[0] = std::make_shared<Organism<FitType>>(best_organism);
      offspring[0] = best_organism;

      for (size_t i = 1; 2*i < this->offspring_num; i++) {
        size_t par1_i = dist_surv0( args.get_generator() );
        //use the survivors_num-1 distribution to guarantee different parents
        size_t par2_i = dist_surv1( args.get_generator() );
        if (par2_i >= par1_i) {
    par2_i++;
        }
        children = survivors[par1_i].get()->breed(&args, survivors[par2_i].get());
        offspring[2*i] = std::shared_ptr<Organism<FitType>>(children[0]);
        offspring[2*i - 1] = std::shared_ptr<Organism<FitType>>(children[1]);
      }
    } else {
      for (size_t i = 0; 2*i + 1 < this->offspring_num; i++) {
        size_t par1_i = dist_surv0( args.get_generator() );
        //use the survivors_num-1 distribution to guarantee different parents
        size_t par2_i = dist_surv1( args.get_generator() );
        if (par2_i >= par1_i) {
    par2_i++;
        }
        children = survivors[par1_i].get()->breed(&args, survivors[par2_i].get());
        offspring[2*i] = std::shared_ptr<Organism<FitType>>(children[0]);
        offspring[2*i + 1] = std::shared_ptr<Organism<FitType>>(children[1]);
      }
    }

    free(shuffled_inds);
    offspring.swap(old_gen);
  }
  void tournament_selection() {
    _uint arena_size = args.get_survivors();
    SampleDraw sampler(offspring_num, arena_size, args.tournament_replacement());
    std::uniform_int_distribution<_uint> selector(0, offspring_num - 1);
    std::vector<Organism<FitType>*> children;

    for (size_t i = 0; 2*i + 1 < offspring_num; ++i) {
      std::shared_ptr< Organism<FitType> > first_parent, second_parent;
      std::vector<_uint> t1 = sampler( args.get_generator() );
      std::vector<_uint> t2 = sampler( args.get_generator() );
      first_parent = old_gen[t1[0]];
      second_parent = old_gen[t2[0]];
      for (size_t j = 1; j < t1.size(); ++j) {
        //check whether the fitness is an improvement and use variance as a tiebreaker
        if ( sel.compare(old_gen[t1[j]], first_parent) ) {
          first_parent = old_gen[t1[0]];
        }
        if ( sel.compare(old_gen[t2[j]], second_parent) ) {
          second_parent = old_gen[t2[0]];
        }
      }
      //ensure that we don't use the same parent twice with reasonable probability
      if (first_parent == second_parent) {
        second_parent = old_gen[selector( args.get_generator() )];
      }

      children = first_parent->breed(&args, second_parent.get());
      offspring[2*i] = std::shared_ptr<Organism<FitType>>(children[0]);
      offspring[2*i + 1] = std::shared_ptr<Organism<FitType>>(children[1]);
    }
    if (offspring_num % 2 == 1) {
      //offspring[offspring_num - 1] = std::make_shared<Organism<FitType>>(best_organism);
      offspring[offspring_num - 1] = best_organism;
    } else {
      //offspring[0] = std::make_shared<Organism<FitType>>(best_organism);
      offspring[0] = best_organism;
    }
    old_gen.swap(offspring);
    calculated_flags &= !FLAG_FRONTS;
  }
  void calculate_distances() {
    for (_uint i = 0; i < offspring_num; ++i) {
      old_gen[i]->get_fitness_info().distance = 0;
    }
    for (_uint i = 0; i < N_OBJS; ++i) {
      sort_orgs(i, &old_gen);
      for (_uint j = 1; j < offspring_num - 1; ++j) {
        double tmp_dist = old_gen[j + 1]->get_fitness(i) - old_gen[j - 1]->get_fitness(i);
        old_gen[j]->get_fitness_info().distance += tmp_dist*tmp_dist;
      }
      old_gen[0]->get_fitness_info().distance = std::numeric_limits<double>::infinity();
      old_gen[offspring_num - 1]->get_fitness_info().distance = std::numeric_limits<double>::infinity();
    }
    calculated_flags |= FLAG_DIST_SET;
  }
  void hypermutate() {
    if ( !(calculated_flags & FLAG_DIST_SET) ) {
      calculate_distances();
    }
    if ( !(calculated_flags & FLAG_BEST_FOUND) ) {
      find_best_organism();
    }
    //sort by distance
    sort_orgs(get_n_objs(), &old_gen);
    //select the half of the most crowded individuals in the first front
    for (_uint i = offspring_num - 1; i > args.get_replacement_fraction()*offspring_num; --i) {
  //    if (old_gen[i] != best_organism_current) {
        old_gen[i]->randomize(&args);
  //    }
    }
  }
  void set_best_organism(_uint i, bool force=false) {
    pop_stats[0].max = old_gen[i]->get_fitness(0);
    /*Organism<FitType> tmp_org = old_gen[i]->copy();
    if (force || tmp_org > best_organism) {
      best_organism = tmp_org;
      for (_uint j = 0; j < N_OBJS; ++j) {
        best_organism->set_fitness( j, old_gen[i]->get_fitness(j) );
      }
      best_organism_ind = i;
    }*/
    if ( force || !best_organism || sel.compare(old_gen[i], best_organism) ) {
      best_organism = old_gen[i];
      best_organism_ind = i;
    }
    calculated_flags |= VALID_BEST;
  }

public:
//    Population(_uint pn_bits, _uint pn_objs, PhenotypeMap* p_map, ArgStore p_args);
//    Population(_uint pn_bits, _uint pn_objs, Organism<FitType>* tmplt, PhenotypeMap* p_map, ArgStore p_args);
  Population(_uint pn_bits, _uint pn_objs, std::shared_ptr<PhenotypeMap> p_map, bool latin=true) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs),
  is_obj_cost(N_OBJS, false)
  {
    map = p_map;
    createOrganisms(NULL, latin);
  }
  Population(_uint pn_bits, _uint pn_objs, Organism<FitType>* tmplt, std::shared_ptr<PhenotypeMap> p_map) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs),
  is_obj_cost(N_OBJS, false)
  {
    map = p_map;
    createOrganisms(tmplt, false);
  }
  Population(_uint pn_bits, _uint pn_objs, std::shared_ptr<PhenotypeMap> p_map, ArgStore p_args, bool latin=true) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs),
  args(p_args),
  is_obj_cost(N_OBJS, false)
  {
    map = p_map;
    createOrganisms(NULL, latin);
  }
  Population(_uint pn_bits, _uint pn_objs, Organism<FitType>* tmplt, std::shared_ptr<PhenotypeMap> p_map, ArgStore p_args) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs),
  args(p_args),
  is_obj_cost(N_OBJS, false)
  {
    map = p_map;
    createOrganisms(tmplt, false);
  }
  void createOrganisms(Organism<FitType>* tmplt, bool latin) {
    pop_stats = (FitnessStats*)malloc(sizeof(FitnessStats)*N_OBJS);
    this->survivors_num = args.get_survivors();
    this->offspring_num = args.get_pop_size();
    if (this->offspring_num % 2 == 0) {
      this->offspring_num++;
    }
    //we need to keep the old and new generation in separate arrays to avoid overwriting data
    this->offspring.insert(this->offspring.end(), this->offspring_num, std::shared_ptr<Organism<FitType>>(NULL));
    if (tmplt) {
      this->old_gen.push_back(std::make_shared<Organism<FitType>>(*tmplt));
      //initally fill up the offspring randomly
      for (size_t i = 1; i < this->offspring_num; ++i) {
        this->old_gen.push_back( std::make_shared<Organism<FitType>>(N_BITS, N_OBJS, map) );
        this->old_gen.back()->randomize(&args, tmplt);
      }
    } else {
      if (latin) {
        Shuffle samp(offspring_num);
        std::uniform_real_distribution<double> in_cube_dist(0, 1);
        Vector<_uint> row;
        
        this->old_gen.reserve(this->offspring_num);
        for (size_t i = 0; i < this->offspring_num; ++i) {
          this->old_gen.push_back( std::make_shared<Organism<FitType>>(N_BITS, N_OBJS, map) );
        }
        for (size_t i = 0; i < map->get_num_params(); ++i) {
          row = samp( args.get_generator() ); 
          double row_width = ( map->get_range_max(i) - map->get_range_min(i) )/offspring_num;
          for (size_t j = 0; j < offspring_num; ++j) {
            double x = in_cube_dist( args.get_generator() );
            x = (x + row[j])*row_width + map->get_range_min(i);
            this->old_gen[j]->set_real(i, x);
          }
        }
      } else {
        this->old_gen.reserve(this->offspring_num);
        //initally fill up the offspring randomly
        for (size_t i = 0; i < this->offspring_num; ++i) {
          this->old_gen.push_back( std::make_shared<Organism<FitType>>(N_BITS, N_OBJS, map) );
          this->old_gen[i]->randomize(&args);
        }
      }
    }
    for (_uint i = 0; i < N_OBJS; ++i) {
      this->pop_stats[i].max = -std::numeric_limits<double>::infinity();
      this->pop_stats[i].min = std::numeric_limits<double>::infinity();
    }
    char buf[OUT_BUF_SIZE];
  //  var_labels.resize(map->get_num_params());
    N_PARAMS = map->get_num_params();
    var_labels = (char**)malloc(sizeof(char*)*N_PARAMS);
    for (_uint j = 0; j < map->get_num_params(); ++j) {
      var_labels[j] = (char*)malloc(sizeof(char)*OUT_BUF_SIZE);
      snprintf(var_labels[j], OUT_BUF_SIZE, "x_%d", j);
    }
    obj_labels = (char**)malloc(sizeof(char*)*N_OBJS);
    for (_uint j = 0; j < N_OBJS; ++j) {
      obj_labels[j] = (char*)malloc(sizeof(char)*OUT_BUF_SIZE);
      snprintf(obj_labels[j], OUT_BUF_SIZE - 1, "f_%d(x)", j);
    }

    calculated_flags = FLAG_NONE_SET;
  }
  ~Population() {
    for (size_t i = 0; i < this->offspring_num; ++i) {
      old_gen[i].reset();
    }
    if (pop_stats) {
      free(pop_stats);
    }
    if (var_labels) {
      for (_uint i = 0; i < N_PARAMS; ++i) {
        free(var_labels[i]);
      }
      free(var_labels);
    }
    if (obj_labels) {
      for (_uint i = 0; i < N_OBJS; ++i) {
        free(obj_labels[i]);
      }
      free(obj_labels);
    }
  }
  Population(Population& o) :
  args(o.args),
  best_organism(o.best_organism)
  {
    N_BITS = o.get_n_bits();
    N_PARAMS = o.N_PARAMS;
    N_OBJS = o.N_OBJS;
    map = o.map;
    old_gen = o.old_gen;
    offspring = o.offspring;
    pop_stats = (FitnessStats*)malloc(sizeof(FitnessStats)*N_OBJS);
    var_labels = (char**)malloc(sizeof(char*)*N_PARAMS);
    obj_labels = (char**)malloc(sizeof(char*)*N_OBJS);

    for (_uint i = 0; i < N_OBJS; ++i) {
      pop_stats[i].min = o.pop_stats[i].min;
      pop_stats[i].max = o.pop_stats[i].max;
      pop_stats[i].var = o.pop_stats[i].var;
      size_t len = strlen(o.obj_labels[i]) + 1;
      obj_labels[i] = (char*)malloc( sizeof(char)*len);
      snprintf(obj_labels[i], len, "%s", o.obj_labels[i]);
    }
    for (_uint i = 0; i < N_PARAMS; ++i) {
      size_t len = strlen(o.var_labels[i]) + 1;
      var_labels[i] = (char*)malloc(sizeof(char)*len);
      snprintf(obj_labels[i], len, "%s", o.obj_labels[i]);
    }
    calculated_flags = FLAG_NONE_SET;
  }
  Population& operator=(Population& o) {
    int tmp_N_BITS = N_BITS;
    int tmp_N_OBJS = N_OBJS;
    std::shared_ptr<PhenotypeMap> tmp_map = map;
    FitnessStats* tmp_pop_stats = pop_stats;
    char** tmp_var_labels = var_labels;
    char** tmp_obj_labels = obj_labels;
    N_BITS = o.N_BITS;
    N_OBJS = o.N_OBJS;
    map = o.map;
    pop_stats = tmp_pop_stats;
    var_labels = o.var_labels;
    obj_labels = o.obj_labels;
    o.N_BITS = tmp_N_BITS;
    o.N_OBJS = tmp_N_OBJS;
    o.map = tmp_map;
    o.pop_stats = tmp_pop_stats;
    o.var_labels = tmp_var_labels;
    o.obj_labels = tmp_obj_labels;

    old_gen = o.old_gen;
    offspring = o.offspring;
    best_organism = o.best_organism;

    calculated_flags = FLAG_NONE_SET;
    
    return *this;
  }
  Population(Population&& o) :
  map(o.map),
  args(std::move(o.args)),
  best_organism(std::move(o.best_organism))
  {
    N_BITS = o.get_n_bits();
    for (size_t i = 0; i < this->offspring_num; ++i) {
      old_gen[i].reset();
    }
    old_gen = std::move(o.old_gen);
    offspring = std::move(o.offspring);
    pop_stats = o.pop_stats;
    var_labels = o.var_labels;
    obj_labels = o.obj_labels;
    o.pop_stats = NULL;
    o.var_labels = NULL;
    o.obj_labels = NULL;

    calculated_flags = FLAG_NONE_SET;
  }

  void set_convergence_type(ConvergenceCriteria* conv);
  void set_penalty_printing(bool val = true) { print_penalties = (val)? 1 : 0; }
  void resize_population(_uint new_size) {
    size_t old_size = old_gen.size();
    if (new_size > old_size) {
      offspring.insert( offspring.end(), new_size - old_size, std::shared_ptr<Organism<FitType>>(NULL) );
      old_gen.reserve(new_size);
      for (size_t i = old_size; i < old_size; ++i) {
        old_gen.push_back( std::make_shared<Organism<FitType>>(N_BITS, N_OBJS, map) );
        old_gen.back()->randomize(&args);
      }
    } else {
      offspring.resize(new_size);
      old_gen.resize(new_size);
    }
    args.set_pop_size(new_size);

    calculated_flags = FLAG_NONE_SET;
  }
  void set_n_survivors(_uint new_size) {
    size_t old_size = survivors_num;
    if (new_size > old_size) {
      survivors.insert( offspring.end(), new_size - old_size, std::shared_ptr<Organism<FitType>>(NULL) );
    } else {
      survivors.resize(new_size);
    }
    args.set_survivors(new_size);
  }
#ifdef USE_LIBOMP
  void evaluate_async(Problem* prob) {
    if (N_OBJS == 1) { 
      for (_uint i = 0; i < offspring_num; ++i) {
        old_gen[i]->apply_penalty(0);
      }
      //calculate averages for organisms that appear twice in the population
      if ( args.average_multiples() ) {
#pragma omp parallel for
        for (_uint i = 0; i < offspring_num; ++i) {
          Vector<_uint> identical_set;
          double avg_fit = 0.0;
          bool apply_averages = true;
          for (_uint j = 0; j < offspring_num; ++j) {
            if ( i == j || *(old_gen[j]) == *(old_gen[i]) ) {
              identical_set.push_back(j);
              //ensure that we only calculate the identical set once
              if (i < j) {
                apply_averages = false;
              }
            }
          }

#pragma omp parallel for
          for (_uint j = 0; j < identical_set.size(); ++j) {
            for (_uint k = 0; k < args.noise_compensate() + 1; ++k) {
              old_gen[i]->evaluate_fitness(prob, args.forget_weight);
            }
          }
          
          //don't recalculate if we don't have to
          if (apply_averages) {
            for (auto it = identical_set.begin(); it != identical_set.end(); ++it) {
              old_gen[*it]->update(0, old_gen[i]->get_fitness(0));
            }
            if (old_gen[i]->get_fitness(0) > best_organism->get_fitness(0) && !old_gen[i]->penalized()) {
                set_best_organism(i);
            }
          }
        }
      } else {
        Vector<_uint*> skip_set;
        for (_uint i = 0; i < this->offspring_num; ++i) {
          if (args.skip_multiples()) {
            for (size_t j = 0; j < i; ++j) {
              if (*(old_gen[j]) == *(old_gen[i])) {
                _uint* tmp = (_uint*)malloc(sizeof(_uint)*2);
                tmp[0] = j; tmp[1] = i;
                skip_set.push_back(tmp);
              }
            }
          } else if (args.perturb_multiples()) {
            for (size_t j = 0; j < i; ++j) {
              if (*(old_gen[j]) == *(old_gen[i])) {
                old_gen[j]->mutate(&args);
              }
            }
          }
        }

#pragma omp parallel for
        for (_uint i = 0; i < offspring_num; ++i) {
          _uint j = 0;
          for (; j < skip_set.size(); ++j) {
            if ( skip_set[j][0] == i ) { break; }
          }
          if (j >= skip_set.size()) {
            if ( args.verbose() ) {
              std::cout << "Now evaluating organism " << i << std::endl;
            }
            old_gen[i]->evaluate_fitness(prob);
            for (_uint k = 0; k < args.noise_compensate(); ++k) {
              old_gen[i]->evaluate_fitness(prob, args.forget_weight);
            }
          }
        }

        _uint first_valid_i = 0;
        if ( best_organism->valid() ) {
          if ( args.noise_compensate() ) {
            evaluate_best(prob, args.forget_weight);
          }
          pop_stats[0].max = best_organism->get_fitness(0);
          pop_stats[0].min = best_organism->get_fitness(0);
        } else {
          //iterate until we find an organism that isn't penalized and set it to be the best
          do {
            if (first_valid_i == offspring_num) {
              error(CODE_MISC, "All organisms in population had applied penalty.");
            }
            ++first_valid_i;
          } while( old_gen[first_valid_i]->penalized() );
          set_best_organism(first_valid_i - 1);
          //alltime_best_organism = best_organism->copy();
    alltime_best_organism = best_organism;
          pop_stats[0].max = old_gen[first_valid_i]->get_fitness(0);
          pop_stats[0].min = old_gen[first_valid_i]->get_fitness(0);
        }

        //this can't be parallelized easily
        for (_uint i = 0; i < offspring_num; ++i) {
          _uint j = 0;
          for (; j < skip_set.size(); ++j) {
            if ( skip_set[j][0] == i ) { break; }
          }
          //if we are in the skip set, then set fitness accordingly
          if (j < skip_set.size()) {
            uint prev_ind = skip_set[j][1];
            old_gen[i]->update(0, old_gen[prev_ind]->get_fitness(0));
            old_gen[i]->apply_penalty(old_gen[prev_ind]->get_penalty());
          } else if (old_gen[i]->get_fitness(0) > best_organism->get_fitness(0) && !old_gen[i]->penalized()) {
            //check the organism again to make sure this isn't a fluke
            for (_uint j = 0; j < args.noise_compensate(); ++j) {
              old_gen[i]->evaluate_fitness(prob, args.forget_weight);
              evaluate_best(prob, args.forget_weight);
            }
            set_best_organism(i);
            if (old_gen[i]->get_fitness(0) < pop_stats[0].min) {
              pop_stats[0].min = old_gen[i]->get_fitness(0);
            }
          }
        }

        //free allocated memory
        for (_uint j = 0; j < skip_set.size(); ++j) {
          free(skip_set[j]);
        }
      }
    } else {
      //TODO: figure out what the default behavior should be
    }
    double penalty_fact;
    if (pop_stats[0].max > 0) {
      penalty_fact = pop_stats[0].min;
    } else {
      penalty_fact = -pop_stats[0].min;
    }
    bool penalties_applied = false;
#pragma omp parallel for
    for (size_t i = 0; i < this->offspring_num; ++i) {
      if (old_gen[i]->penalized()) {
        old_gen[i]->update(pop_stats[0].min - penalty_fact*old_gen[i]->get_penalty());
        penalties_applied = true;
      }
    }
    if (!penalties_applied) {
      calculated_flags |= FLAG_STATS_SET | FLAG_BEST_FOUND;
    } else {
      calculated_flags = FLAG_NONE_SET;
    }
    apply_penalties(prob);
  }
#endif
  //function for population where organism FitType has a member average_fitness
  template <typename T = FitType> inline
  typename enable_if_c< has_average_fitness<T, void(T&)>::value, void >::type
  evaluate(Problem<T>* prob) {
    std::cout << "Called evaluate with " << typeid(T).name() << " which does implement average_fitness.\n";
    size_t start_i = 0;
    if ( this->best_organism && this->best_organism->valid() ) {
      if ( this->args.noise_compensate() ) {
        this->evaluate_best(prob, this->args.forget_weight);
      }
      this->pop_stats[0].max = this->best_organism->get_fitness(0);
    } else {
      //iterate until we find an organism that isn't penalized and set it to be the best
      do {
        if (start_i == this->offspring_num) {
          error(CODE_MISC, "All organisms in population had an applied penalty.");
        }
        this->old_gen[start_i]->apply_penalty(0);
        if ( this->args.noise_compensate() ) {
          this->old_gen[start_i]->reset_fitness();
          this->old_gen[start_i]->evaluate_fitness(prob);
        } else {
          this->old_gen[start_i]->reset_fitness();
          this->old_gen[start_i]->evaluate_fitness(prob);
        }
        ++start_i;
      } while( this->old_gen[start_i - 1]->penalized() );
      this->set_best_organism(start_i - 1, this->args.noise_compensate());
      //alltime_best_organism = best_organism->copy();
      this->alltime_best_organism = this->best_organism;
    }
    this->pop_stats[0].min = this->best_organism->get_fitness(0);
   
    for (_uint i = start_i; i < this->offspring_num; ++i) {
      this->old_gen[i]->apply_penalty(0);
      if ( args.perturb_multiples() ) {
        for (_uint j = 0; j < i; ++j) {
          if ( *(this->old_gen[i]) == *(this->old_gen[j]) ) {
            this->old_gen[i]->mutate( &(this->args) );
          }
        }
      }
      this->old_gen[i]->evaluate_fitness(prob);
      //average fitnesses with previous organisms if appropriate
      if ( this->args.average_multiples() ) {
        for (_uint j = 0; j < i; ++j) {
          if ( *(this->old_gen[i]) == *(this->old_gen[i]) ) {
            this->old_gen[i]->get_fitness_info().average_fitness( this->old_gen[j]->get_fitness_info() );
          }
        }
      }
      // update the max and min fitnesses if we need to
      if ( this->old_gen[i]->get_fitness(0) > this->best_organism->get_fitness(0)
      && !(this->old_gen[i]->penalized()) ) {
        //check the organism again to make sure this isn't a fluke
        for (_uint j = 0; j < this->args.noise_compensate(); ++j) {
          this->old_gen[i]->evaluate_fitness(prob);
          //evaluate_best(prob, args.forget_weight);
        }
        if (!(this->old_gen[i]->penalized())
         &&  (this->old_gen[i]->get_fitness(0) > this->best_organism->get_fitness(0)
           || this->old_gen[i]->get_fitness(0) > pop_stats[0].max)) {
          this->set_best_organism(i);
        }
      }
      if (this->old_gen[i]->get_fitness(0) < this->pop_stats[0].min) {
        this->pop_stats[0].min = this->old_gen[i]->get_fitness(0);
      }
    }
    this->apply_penalties(prob);
    //TODO: make this more flexible and allow for custom comparators
    for (_uint i = 0; i < this->offspring_num; ++i) {
      if (!(this->old_gen[i]->penalized())
       &&  (this->old_gen[i]->get_fitness(0) > this->best_organism->get_fitness(0)
         || this->old_gen[i]->get_fitness(0) > pop_stats[0].max)) {
        this->set_best_organism(i, true);
      }
    }
    this->pop_stats[0].max = this->best_organism->get_fitness();
  }
  //function for population where organism FitType does not have a member average_fitness
  template <typename T = FitType> inline
  typename enable_if_c< !has_average_fitness<T, void(T&)>::value, void >::type
  evaluate(Problem<T>* prob) {
    //std::cout << "Called evaluate with " << typeid(T).name() << " which does not implement average_fitness.\n";
    size_t start_i = 0;
    if ( !(this->best_organism) || !(this->best_organism->valid()) ) {
      //iterate until we find an organism that isn't penalized and set it to be the best
      do {
        if (start_i == this->offspring_num) {
          error(CODE_MISC, "All organisms in population had an applied penalty.");
        }
        this->old_gen[start_i]->apply_penalty(0);
        this->old_gen[start_i]->evaluate_fitness(prob);
        ++start_i;
      } while( this->old_gen[start_i - 1]->penalized() );
      this->set_best_organism(start_i - 1, true);
      //alltime_best_organism = best_organism->copy();
      this->alltime_best_organism = this->best_organism;
    }
    this->pop_stats[0].min = this->best_organism->get_fitness(0);
   
    for (_uint i = start_i; i < this->offspring_num; ++i) {
      this->old_gen[i]->apply_penalty(0);
    }
    //calculate averages for organisms that appear twice in the population
    for (_uint i = start_i; i < this->offspring_num; ++i) {
      if ( this->args.verbose() ) {
        std::cout << "Now evaluating organism " << i << std::endl;
      }

      bool found_identical = false;
      //look for duplicates of the current organism
      for (size_t j = 0; j < i; ++j) {
        //handle them
        if ( *(this->old_gen[j]) == *(this->old_gen[i]) ) {
          this->old_gen[i]->update( 0, this->old_gen[j]->get_fitness() );
          found_identical = true;
          break;
        }
      }
      if (!found_identical) {
        this->old_gen[i]->evaluate_fitness(prob);
      }
      // update the max and min fitnesses if we need to
      if (this->old_gen[i]->get_fitness(0) > this->best_organism->get_fitness(0) && !(this->old_gen[i]->penalized())) {
        this->set_best_organism(i);
      }

      if (this->old_gen[i]->get_fitness(0) < this->pop_stats[0].min) {
        this->pop_stats[0].min = this->old_gen[i]->get_fitness(0);
      }
    }
    this->pop_stats[0].max = this->best_organism->get_fitness();

    this->apply_penalties(prob);
  }
  //void evaluate(Problem<FitType>* prob) { evaluate_imp(prob, NULL); }
  bool iterate(ConvergenceCriteria* conv = NULL) {
    find_best_organism();
    //check for hypermutation
    for (_uint i = 0; i < N_OBJS; ++i) {
      double range_ratio = (pop_stats[i].max - pop_stats[i].min)/ pop_stats[i].max;
      if (range_ratio < 0) {
        range_ratio *= -1;
      }
      if (1.0 - range_ratio > args.get_hypermutation_threshold()) {
        hypermutate();
        break;
      }
    }
    if (args.use_tournament()) {
      tournament_selection();
    } else {
      if (cull_in_place()) {
        return true;
      }
      breed();
    } 
    generation++;
    if (conv) {
    //TODO: implement a binary private member variable that tracks whether sorting needs to be performed for the sake of efficiency
      return conv->evaluate_convergence(1, pop_stats);
    } else {
      return (generation < args.get_num_gens());
    }
    calculated_flags = FLAG_NONE_SET;
  }
  void run(Problem<FitType>* prob) {
    evaluate(prob);
    if ( this->args.wait_for_con() ) {
      size_t i = 1;
      unsigned streak = 0;
      double prev_ftns = this->get_best_organism()->get_fitness(0);
      while (i < MAX_NUM_GENS) {
        if (this->iterate()) {
          break;
        }
#ifdef USE_LIBOMP
        if (args.async()) {
          evaluate_async(prob);
        } else {
          evaluate(prob);
        }
#else
        evaluate(prob);
#endif

        if (this->get_best_organism()->get_fitness(0) > prev_ftns) {
          prev_ftns = this->get_best_organism()->get_fitness(0);
          streak = 0;
        }

        streak++;
        i++;
        if (streak > this->args.get_num_gens()) {
          break;
        }
      }

      if (i >= MAX_NUM_GENS) {
        std::cout << "failed to converge after " << i << " generations." << std::endl;
      } else {
        std::cout << "converged to result after " << i << " generations." << std::endl;
      }
    } else {
      if (this->args.verbose()) {
        std::cout << "Now evaluating generation 0..." << std::endl;
      }
      for (size_t i = 1; i < this->args.get_num_gens() + 1; ++i) {
        if (this->args.verbose()) {
          std::cout << "Now evaluating generation " << i << "..." << std::endl;
        }
        //produce the next generation in the population
        if (this->iterate()) {
          break;
        }
        this->evaluate(prob);
      } 
    }
  }
  std::shared_ptr< Organism<FitType> > get_best_organism(size_t i = 0) {
    if ( (calculated_flags & FLAG_BEST_FOUND) == 0 ) {
      find_best_organism();
    }
    if (i == 0) {
      std::shared_ptr<Organism<FitType>> tmp_org = std::make_shared<Organism<FitType>>( best_organism->copy() );
      /*for (_uint i = 0; i < N_OBJS; ++i) {
        tmp_org->update(i, pop_stats[i].max);
      }*/
      return tmp_org;
    } else {
      sort_orgs(0, &old_gen);
      return old_gen[i];
    }
  }
  std::shared_ptr< Organism<FitType> > get_organism(size_t i) {
    if (i >= old_gen.size())
      error(CODE_ARG_RANGE, "Attempt to access invalid index %d when the maximum allowed is %d.", i, old_gen.size());
    return old_gen[i];
  }
  std::shared_ptr< Organism<FitType> > get_child(size_t i) {
    if (i >= offspring.size())
      error(CODE_ARG_RANGE, "Attempt to access invalid index %d when the maximum allowed is %d.", i, offspring.size());
    return offspring[i];
  }
  void swap_orgs(int i, int j) {
    std::shared_ptr<Organism<FitType>> tmp = old_gen[j];
    old_gen[j] = old_gen[i];
    old_gen[i] = tmp;
  }
  int partition(_uint fit_ind, std::vector<std::shared_ptr< Organism<FitType> >>* work_arr, int s, int e) {
    double p = (*work_arr)[e]->get_fitness(fit_ind);
    int i = s;
    for (int j = s; j < e; ++j) {
      if ((*work_arr)[j]->get_fitness(fit_ind) > p) {
        (*work_arr)[i].swap( (*work_arr)[j] );
        i++;
      }
    }
    if (i != e) {
      std::shared_ptr<Organism<FitType>> tmp = (*work_arr)[i];
      (*work_arr)[i] = (*work_arr)[e];
      (*work_arr)[e] = tmp;
    }
    return i;
  }
  void sort_orgs(unsigned int fit_ind, std::vector<std::shared_ptr< Organism<FitType> >>* work_arr, int s = DEF_SORT_PARAM, int e = DEF_SORT_PARAM) {
    sort_org_calls += 1;
    if (sort_org_calls > work_arr->size()*5) {
      error(CODE_WARN, "Suspiciously many calls (%u) to sort org have been made", sort_org_calls);
    }
    if (s == DEF_SORT_PARAM || e == DEF_SORT_PARAM) {
      sort_org_calls = 0;
      sort_orgs(fit_ind, work_arr, 0, work_arr->size() - 1);
    } else if (s < e) {
      int p = partition(fit_ind, work_arr, s, e);
      sort_orgs(fit_ind, work_arr, s, p-1);
      sort_orgs(fit_ind, work_arr, p+1, e);
    }
  }
  Vector<String> get_best_header() {
    String def;
    Vector<String> ret;
    ret.reserve( N_PARAMS + print_penalties + N_OBJS );
    char buf[OUT_BUF_SIZE];

    //print information about the best organism
    for (_uint j = 0; j < N_PARAMS; ++j) {
      snprintf(buf, OUT_BUF_SIZE - 1, "best %s", var_labels[j]);
      ret.push_back( String(buf) );
    }
    if (print_penalties != 0) {
      ret.push_back( String("best penalty") );
    }
    for (_uint j = 0; j < N_OBJS; ++j) {
      snprintf(buf, OUT_BUF_SIZE - 1, "best %s", obj_labels[j]);
      ret.push_back( String(buf) );
    }
    
    return ret;
  }
  Vector<String> get_header() {
    String def;
    Vector<String> ret;
    ret.reserve( old_gen.size()*(N_PARAMS + print_penalties + N_OBJS) );

    //print information about every other organism in the population
    for (_uint i = 0; i < old_gen.size(); ++i) {
      for (_uint j = 0; j < N_PARAMS; ++j) {
        ret.push_back( String(var_labels[j]) );
      }

      if (print_penalties != 0) {
        ret.push_back( String("Penalized") );
      }

      for (_uint j = 0; j < N_OBJS; ++j) {
        ret.push_back( String(obj_labels[j]) );
      }
    }
    
    return ret;
  }
  Vector<String> get_best_data() {
    sort_orgs(0, &old_gen);
    char buf[OUT_BUF_SIZE];
    String def;
    Vector<String> ret(N_OBJS + print_penalties + map->get_num_params(), def);

    size_t param_o = 0;
    size_t fitness_o = map->get_num_params() + print_penalties;
    //print information about the best organism
    for (_uint j = 0; j < map->get_num_params(); ++j) {
      ret[j] = best_organism->get_chromosome_string(j);
    }
    if (print_penalties) {
      snprintf(buf, OUT_BUF_SIZE - 1, "%f", best_organism->get_penalty());
      ret[map->get_num_params()] = buf;
    }
    for (_uint j = 0; j < N_OBJS; ++j) {
      if ( args.noise_compensate() ) {
        double std_dev = sqrt( best_organism->get_fitness_info().get_uncertainty() );
        if (is_obj_cost[j]) {
          snprintf(buf, OUT_BUF_SIZE - 1, "%f\u00B1%f", best_organism->get_cost(j), std_dev);
        } else {
          snprintf(buf, OUT_BUF_SIZE - 1, "%f\u00B1%f", best_organism->get_fitness(j), std_dev);
        }
      } else {
        if (is_obj_cost[j]) {
          snprintf(buf, OUT_BUF_SIZE - 1, "%f", best_organism->get_cost(j));
        } else {
          snprintf(buf, OUT_BUF_SIZE - 1, "%f", best_organism->get_fitness(j));
        }
      }
      ret[fitness_o + j] = buf;
    }

    return ret;
  }
  Vector<String> get_pop_data() {
    sort_orgs(0, &old_gen);
    _uint span = N_OBJS + print_penalties + map->get_num_params();
    String def;
    Vector<String> ret(span*offspring_num, def);
    char buf[OUT_BUF_SIZE];

    //print information about the rest of the population
    for (_uint i = 0; i < old_gen.size(); ++i) {
      size_t param_o = i*span;
      size_t fitness_o = param_o + map->get_num_params() + print_penalties;

      //print out the parameters
      for (_uint j = 0; j < map->get_num_params(); ++j) {
        ret[param_o + j] = old_gen[i]->get_chromosome_string(j);
      }
      //print out the penalty applied to the organism
      if (print_penalties) {
        snprintf(buf, OUT_BUF_SIZE - 1, "%f", old_gen[i]->get_penalty());
        ret[param_o + map->get_num_params()] = buf;
      }
      //print out the fitness value(s)
      for (_uint j = 0; j < N_OBJS; ++j) {
        if ( args.noise_compensate() ) {
          double std_dev = sqrt( old_gen[i]->get_fitness_info().get_uncertainty() );
          if (is_obj_cost[j]) {
            snprintf(buf, OUT_BUF_SIZE - 1, "%f\u00B1%f", old_gen[i]->get_cost(j), std_dev);
          } else {
            snprintf(buf, OUT_BUF_SIZE - 1, "%f\u00B1%f", old_gen[i]->get_fitness(j), std_dev);
          }
        } else {
          if (is_obj_cost[j]) {
            snprintf(buf, OUT_BUF_SIZE - 1, "%f", old_gen[i]->get_cost(j));
          } else {
            snprintf(buf, OUT_BUF_SIZE - 1, "%f", old_gen[i]->get_fitness(j));
          }
        }
        
        ret[fitness_o + j] = buf;
      }
    }

    return ret;
  }

  FitnessStats get_pop_stats(_uint i = 0) { return pop_stats[i]; }
  DEPRECATED("get_min_fitness is deprecated, use get_pop_stats instead") double get_min_fitness(_uint i = 0) {
    return pop_stats[i].min;
  }
  DEPRECATED("get_max_fitness is deprecated, use get_pop_stats instead") double get_max_fitness(_uint i = 0) {
    return pop_stats[i].max;
  }

  /**
   * \brief Sets the objective referenced by index ind to use a fitness (corresponding to a maximization problem).
   *
   * \param ind	The index of the parameter to set
   * \seealso set_cost
   */
  void update(_uint ind) { is_obj_cost[ind] = false; }
  /**
   * \brief Sets the objective referenced by index ind to use a cost (corresponding to a minimization problem).
   *
   * \param ind	The index of the parameter to set
   * \seealso update
   */
  void set_cost(_uint ind) { is_obj_cost[ind] = true; }
  void set_var_label(_uint ind, String val) {
    if (ind >= N_PARAMS) {
      error(CODE_WARN, "Invalid parameter index %u provided for set_var_label. The parameter index must be less than %u.", ind, N_PARAMS);
    } else {
      size_t vs = val.size() + 1;
      if (vs > OUT_BUF_SIZE) {
        free(var_labels[ind]);
        var_labels[ind] = (char*)malloc(sizeof(char)*vs);
      }
      for (int j = 0; j < vs; ++j) {
        var_labels[ind][j] = val[j];
      }
      var_labels[ind][vs - 1] = 0;
    }
  }
  void set_obj_label(_uint ind, String val) {
    if (ind >= N_OBJS) {
      error(CODE_WARN, "Invalid objective index %u provided for set_obj_label. The objective index must be less than %u.", ind, N_OBJS);
    } else {
      size_t vs = val.size() + 1;
      if (vs > OUT_BUF_SIZE) {
        free(obj_labels[ind]);
        obj_labels[ind] = (char*)malloc(sizeof(char)*vs);
      }
      for (int j = 0; j < vs; ++j) {
        obj_labels[ind][j] = val[j];
      }
      obj_labels[ind][vs - 1] = 0;
    }
  }

  size_t get_offspring_num() { return offspring_num; }
  size_t get_survivors_num() { return survivors_num; }
  _uint get_n_bits() { return N_BITS; }
  _uint get_n_objs() { return N_OBJS; }
  ArgStore& get_args() { return args; }
};

template <class FitType, class SelectType>
class Population_NSGAII : public Population<FitType, SelectType> {
private:
  _uchar calculated_flags = FLAG_NONE_SET;
  size_t survivors_num;
  _uint generation = 0;
  //the ngsa alternative to elitism
  std::vector<std::vector<std::shared_ptr< Organism<FitType> >>> pareto_fronts;
  void hypermutate() {
    this->sort_orgs(this->get_n_objs(), &(pareto_fronts[0]));
    //select the half of the most crowded individuals in the first front
    for (_uint i = pareto_fronts[0].size() - 1; i > pareto_fronts[0].size()/2; --i) {
      pareto_fronts[0][i]->randomize(&(this->args));
    }
  }
  void update_fitness_ranges(Organism<FitType>* org, size_t i) {
    if (org->penalized()) {
      this->max_penalty_ind = i + 1;
      if (this->min_penalty_ind == this->offspring_num) {
        this->min_penalty_ind = i;
      }
    } else {
      for (_uint j = 0; j < this->get_n_objs(); ++j) {
        if (org->get_fitness(j) > this->pop_stats[j].max) {
          this->pop_stats[j].max = org->get_fitness(j);
        }
        if (org->get_fitness(j) < this->pop_stats[j].min) {
          this->pop_stats[j].min = org->get_fitness(j);
        }
      }
    }
  }
  void make_fronts(std::vector<std::shared_ptr< Organism<FitType> >>* cmb_arr) {
    pareto_fronts.clear();
    std::vector<std::shared_ptr<Organism<FitType>>> empty;
    pareto_fronts.push_back(empty);
  //  std::vector<std::vector<std::shared_ptr<Organism<FitType>>>> dominating(cmb_arr.size(), empty);

    for (size_t i = 0; i < cmb_arr->size(); ++i) {
      (*cmb_arr)[i]->get_fitness_info().n_dominations = 0;
      for (size_t j = 0; j < cmb_arr->size(); ++j) {
        if (i != j) {
          //if the ith solution dominates the jth add the jth entry to the list of dominated solutions, otherwise increment the number of dominating solutions
          if ( this->sel.compare( (*cmb_arr)[i], (*cmb_arr)[j] ) ) {
      //TODO: figure out if we need to use the dominating array for anything
  //          dominating[i].push_back( (*cmb_arr)[j] );
          } else if ( this->sel.compare( (*cmb_arr)[j], (*cmb_arr)[i] ) ) {
            (*cmb_arr)[i]->get_fitness_info().n_dominations++;
          }
        }
      }
      if ( (*cmb_arr)[i]->get_fitness_info().n_dominations == 0) {
        (*cmb_arr)[i]->get_fitness_info().rank = 0;
        pareto_fronts[0].push_back( (*cmb_arr)[i] );
      }
    }

    size_t i = 0;
    while (i < pareto_fronts.size() && pareto_fronts[i].size() != 0) {
      pareto_fronts.push_back(empty);
      for (size_t j = 0; j < pareto_fronts[i].size(); ++j) {
        for (size_t k = 0; k < cmb_arr->size(); ++k) {
          (*cmb_arr)[k]->get_fitness_info().n_dominations--;
          if ((*cmb_arr)[k]->get_fitness_info().n_dominations == 0) {
            (*cmb_arr)[k]->get_fitness_info().rank = i + 1;
            pareto_fronts[i + 1].push_back( (*cmb_arr)[k] );
          }
        }
      }
      i++;
    }
    calculated_flags |= FLAG_FRONTS;
  }

public:
  //Population_NSGAII(_uint pn_bits, _uint pn_objs, PhenotypeMap* p_map, String conf_fname = "");
  //Population_NSGAII(_uint pn_bits, _uint pn_objs, Organism<FitType>* tmplt, PhenotypeMap* p_map, String conf_fname = "");
  Population_NSGAII(_uint pn_bits, _uint pn_objs, std::shared_ptr<PhenotypeMap> p_map) :
  Population<FitType, SelectType>(pn_bits, pn_objs, p_map)
  {
    if (this->offspring_num % 2 == 1) {
      this->offspring_num--;
    }
    this->old_gen.resize(this->offspring_num);
    this->offspring.resize(this->offspring_num);
  }
  Population_NSGAII(_uint pn_bits, _uint pn_objs, Organism<FitType>* tmplt, std::shared_ptr<PhenotypeMap> p_map) :
  Population<FitType, SelectType>(pn_bits, pn_objs, tmplt, p_map)
  {
    if (this->offspring_num % 2 == 1) {
      this->offspring_num--;
    }
    this->old_gen.resize(this->offspring_num);
    this->offspring.resize(this->offspring_num);
  }
  Population_NSGAII(_uint pn_bits, _uint pn_objs, std::shared_ptr<PhenotypeMap> p_map, ArgStore p_args) :
  Population<FitType, SelectType>(pn_bits, pn_objs, p_map, p_args)
  {
    if (this->offspring_num % 2 == 1) {
      this->offspring_num--;
    }
    this->old_gen.resize(this->offspring_num);
    this->offspring.resize(this->offspring_num);
  }
  Population_NSGAII(_uint pn_bits, _uint pn_objs, Organism<FitType>* tmplt, std::shared_ptr<PhenotypeMap> p_map, ArgStore p_args) :
  Population<FitType, SelectType>(pn_bits, pn_objs, tmplt, p_map, p_args)
  {
    if (this->offspring_num % 2 == 1) {
      this->offspring_num--;
    }
    this->old_gen.resize(this->offspring_num);
    this->offspring.resize(this->offspring_num);
  }
  ~Population_NSGAII() {
    for (size_t i = 0; i < this->offspring_num; ++i) {
      if (this->old_gen[i] != NULL) {
        this->old_gen[i].reset();
      }
      if (this->offspring[i] != NULL) {
        this->offspring[i].reset();
      }
    }
  //  free(pop_stats);
  }

  void evaluate(Problem<FitType>* prob) {
    std::vector<std::shared_ptr<Organism<FitType>>> cmb_arr = this->old_gen;
    //initialize max and min fitness values
    for (_uint j = 0; j < this->get_n_objs(); ++j) {
      this->pop_stats[j].max = -std::numeric_limits<double>::infinity();
      this->pop_stats[j].min = std::numeric_limits<double>::infinity();
    }
    this->min_penalty_ind = this->offspring_num, this->max_penalty_ind = 0;
    for (size_t i = 0; i < this->offspring_num; ++i) {
      //fitnesses for the old generation are leftover and already have applied penalties, ensure we don't apply these penalties again
      this->old_gen[i]->apply_penalty(0);
      if (this->offspring[i] != NULL) {
        this->offspring[i]->apply_penalty(0);
        this->offspring[i]->evaluate_fitness(prob);
        this->offspring[i]->get_fitness_info().distance = 0;
        //to maintain elitism we look at both the parent and offspring generations
        cmb_arr.push_back(this->offspring[i]);
        update_fitness_ranges(this->offspring[i].get(), this->offspring_num + i);
      }
      this->old_gen[i]->evaluate_fitness(prob);
      //TODO: make sure that this actually works with references
      this->old_gen[i]->get_fitness_info().distance = 0;//initialize for later crowding calculations
      update_fitness_ranges(this->old_gen[i].get(), i);
    }
    //apply penalties to each organism in the combined array
    for (size_t i = this->min_penalty_ind; i < this->max_penalty_ind; ++i) {
      if (cmb_arr[i]->penalized()) {
        for (size_t j = 0; j < this->get_n_objs(); ++j) {
          double n_fit = this->pop_stats[j].min;
          if (this->pop_stats[j].max > 0) {
            n_fit -= this->pop_stats[j].max*(cmb_arr[i]->get_penalty());
          } else {
            n_fit += this->pop_stats[j].min*(cmb_arr[i]->get_penalty());
          }
          cmb_arr[i]->update(n_fit);
        }
      }
    }
    make_fronts(&cmb_arr);
  }
  std::shared_ptr<PhenotypeMap> get_map() { return this->map; }

  //cull in place is slightly faster but less accurate than the standard cull method
  bool cull_in_place() {
    size_t j = 0;

    double difference = this->pop_stats[0].max - this->pop_stats[0].min;
    //avoid divide by 0
    if (difference == 0) {
      error(CODE_WARN, "All organisms have the same fitness, exiting");
      return true;
    }
    std::uniform_real_distribution<double> dist(0, difference);

    if (this->survivors.size() < survivors_num) {
      this->survivors.resize(survivors_num);
    }
    for (size_t i = 0; i < this->offspring_num && j < survivors_num; ++i) {
      /* if M_f is the maximum fitness and m_f is the minimum the minimum, while x is the
       * fitness of a given organism, then the probability of survival is x/(M_f-m_f) or 1
       * if M_f = m_f */
      if (dist(this->args.get_generator()) < this->old_gen[i]->get_fitness(0) - this->pop_stats[0].min) {
        this->survivors[j] = this->old_gen[i];
        j++;
      }
    }

    std::uniform_int_distribution<int> ind_dist(0, this->offspring_num - 1);
    while (j < survivors_num) {
      this->survivors[j] = this->old_gen[ind_dist( this->args.get_generator() )];
      j++;
    }
    return false;
  }
  /**
   * \brief An implementation of simple roulette selection. Cull first sorts the organisms and selects them based on the ratio of their relative fitness to the total relative fitness.
   *
   * \returns True if all organisms have the same fitness (results have converged).
   */
  bool cull() {
    sort_orgs(0, &(this->old_gen));
    //if all the organisms have the same fitness then reinitialize the population
    if (this->old_gen[0]->get_fitness(0) - this->old_gen[this->offspring_num-1]->get_fitness(0) < 0.001 ) {
      return true;
    }
    double min_fit = this->old_gen[this->offspring_num-1]->get_fitness(0);
    double total_fit = 0;
    for (size_t i = 0; i < this->offspring_num; ++i) {
      total_fit += this->old_gen[i]->get_fitness(0) - min_fit;
    }
    std::uniform_real_distribution<double> dist(0, total_fit);

    //maintain a list of organisms that have already been added so no organism appears twice
    size_t* banned = (size_t*)malloc(sizeof(size_t)*survivors_num);
    for (size_t i = 0; i < survivors_num; ++i) { banned[i] = -1; }

    for (size_t i = 0; i < survivors_num; ++i) {
      double val = dist( this->args.get_generator() );
      size_t j = 0;
      while (val > 0.0) {
        val -= this->old_gen[j]->get_fitness(0) - min_fit;
        j++;
      }
      j -= 1;
      while (contains<size_t>(banned, survivors_num, j)) {
        j++;
        if (j >= survivors_num) {
          j = 0;
        }
      }
      this->survivors[i] = this->old_gen[j];
      banned[i] = j;
    }
    if (this->args.verbose()) {
      std::cerr << "Added orgs to survs:";
      for (size_t i = 0; i < survivors_num; ++i) {
        std::cerr << " " << banned[i];
      }
      std::cerr << std::endl;
    }
    free(banned);
  //  best_organism_ind = 0;
    return false;
  }
  void breed() {
    _uint arena_size = this->args.get_survivors();
    SampleDraw sampler(this->offspring_num, arena_size);
    std::uniform_int_distribution<_uint> selector(0, this->offspring_num - 1);
    std::vector<Organism<FitType>*> children;

    for (size_t i = 0; 2*i + 1 < this->offspring_num; ++i) {
      Organism<FitType> *first_parent, *second_parent;
      std::vector<_uint> t1 = sampler( this->args.get_generator() );
      std::vector<_uint> t2 = sampler( this->args.get_generator() );
      first_parent = this->old_gen[t1[0]].get();
      second_parent = this->old_gen[t2[0]].get();
      for (size_t j = 1; j < t1.size(); ++j) {
        if (this->old_gen[t1[j]]->get_fitness_info().rank < first_parent->get_fitness_info().rank) {
          first_parent = this->old_gen[t1[0]].get();
        }
        if (this->old_gen[t2[j]]->get_fitness_info().rank < second_parent->get_fitness_info().rank) {
          second_parent = this->old_gen[t2[0]].get();
        }
      }
      //ensure that we don't use the same parent twice with reasonable probably
      if (first_parent == second_parent) {
        second_parent = this->old_gen[selector( this->args.get_generator() )].get();
      }

      children = first_parent->breed(&(this->args), second_parent);
      this->offspring[2*i] = std::shared_ptr<Organism<FitType>>(children[0]);
      this->offspring[2*i + 1] = std::shared_ptr<Organism<FitType>>(children[1]);
    }
    if (this->offspring_num % 2 == 1) {
      this->offspring_num--;
      this->offspring.pop_back();
    }
    this->old_gen.swap(this->offspring);
    calculated_flags &= !FLAG_FRONTS;
  }
  void calculate_pop_stats() {
    for (_uint i = 0; i < this->get_n_objs(); ++i) {
      this->pop_stats[i].max = -std::numeric_limits<double>::infinity();
      this->pop_stats[i].min = std::numeric_limits<double>::infinity();
      this->pop_stats[i].mean = 0;
      this->pop_stats[i].var = 0;
      size_t n_orgs = 0;
      //calculate the maximum, minimum and mean fitness values
      for (_uint j = 0; j < pareto_fronts.size(); ++j) {
        for (_uint k = 0; k < pareto_fronts[j].size(); ++k) {
    double fit_ijk = pareto_fronts[j][k]->get_fitness(i);
    if (fit_ijk > this->pop_stats[i].max) {
      this->pop_stats[i].max = fit_ijk;
    }
    if (fit_ijk < this->pop_stats[i].min) {
      this->pop_stats[i].min = fit_ijk;
    }
    this->pop_stats[i].mean += fit_ijk;
    n_orgs++;
        }
      }
      if (n_orgs > 0) {
        this->pop_stats[i].mean /= n_orgs;
      }
      if (n_orgs > 1) {
        //calculate the variance
        for (_uint j = 0; j < pareto_fronts.size(); ++j) {
          for (_uint k = 0; k < pareto_fronts[j].size(); ++k) {
            double fit_ijk = pareto_fronts[j][k]->get_fitness(i);
            this->pop_stats[i].var += (fit_ijk - this->pop_stats[i].mean)*(fit_ijk - this->pop_stats[i].mean)
                                    / (n_orgs - 1);
          }
        }
      }
    }
  }
  bool iterate(ConvergenceCriteria* conv = NULL) {
    if (pareto_fronts.size() == 0 && !(calculated_flags & FLAG_FRONTS)) {
      std::vector<std::shared_ptr<Organism<FitType>>> cmb_arr = this->old_gen;
      cmb_arr.reserve(this->old_gen.size() + this->offspring.size());
      for (_uint i = 0; i < this->offspring.size(); ++i) {
        if (this->offspring[i]) {
          cmb_arr.push_back(this->offspring[i]);
        }
      }
      make_fronts(&cmb_arr);
    }
    //if we exceed the threshold, start hypermutation
    if ((double)(pareto_fronts[0].size())/(this->offspring_num) > this->args.get_hypermutation_threshold()) {
      hypermutate();
    }
    size_t i = 0;
    std::vector<std::shared_ptr<Organism<FitType>>> tmp(this->offspring_num, NULL);
    _uint k = 0;
    while (i < pareto_fronts.size() && (k + pareto_fronts[i].size()) <= this->offspring_num) {
      for (size_t j = 0; j < pareto_fronts[i].size(); ++j) {
        tmp[k] = pareto_fronts[i][j];
        ++k;
      }
      ++i;
    }
    //fill in the last elements from the remaining pareto front ranked according to crowding
    if (k < this->offspring_num) {
      for (size_t ii = 0; ii < pareto_fronts[i].size(); ++ii) {
        pareto_fronts[i][ii]->get_fitness_info().distance = 0;
      }
      //sort by each objective function for crowding evaluation
      for (size_t j = 0; j < this->get_n_objs(); ++j) {
        this->sort_orgs(j, &(pareto_fronts[i]));

        //the highest and lowest values should be considered to have no crowding
        pareto_fronts[i].front()->get_fitness_info().distance = std::numeric_limits<double>::infinity();
        pareto_fronts[i].back()->get_fitness_info().distance = std::numeric_limits<double>::infinity();
        double range = pareto_fronts[i].front()->get_fitness(j) - pareto_fronts[i].back()->get_fitness(j);
        if (range == 0) {
          for (size_t ii = 0; ii < pareto_fronts[i].size(); ++ii) {
            pareto_fronts[i][ii]->get_fitness_info().distance = 0;
          }
        } else {
          for (size_t ii = 1; ii < pareto_fronts[i].size() - 1; ++ii) {
            if (pareto_fronts[i][ii]->get_fitness_info().distance != std::numeric_limits<double>::infinity()) {
              double sj_h = pareto_fronts[i][ii-1]->get_fitness(j);
              double sj_l = pareto_fronts[i][ii+1]->get_fitness(j);
              double d_norm = (sj_h - sj_l) / range;
              pareto_fronts[i][ii]->get_fitness_info().distance += d_norm;
            }
          }
        }
      }
      this->sort_orgs(this->get_n_objs(), &(pareto_fronts[i]));
      size_t j = 0;
      //select the least crowded individuals
      size_t p_i_size = pareto_fronts[i].size();
      while (j < p_i_size && k < this->offspring_num) {
        tmp[k] = pareto_fronts[i][j];
        ++j;
        ++k;
      }
      //we don't need the rest of the elements so they should be deleted
      while (j < pareto_fronts[i].size()) {
        if (pareto_fronts[i][j] != NULL) {
          pareto_fronts[i][j].reset();
          pareto_fronts[i][j] = NULL;
        }
        ++j;
      }
      ++i;
    }
    while (i < pareto_fronts.size()) {
      for (size_t j = 0; j < pareto_fronts[i].size(); ++j) {
        if (pareto_fronts[i][j] != NULL) {
          pareto_fronts[i][j].reset();
          pareto_fronts[i][j] = NULL;
        }
      }
      ++i;
    }
    this->old_gen = tmp;
    //ensure that we aren't keeping null pointers around for safety
    pareto_fronts.clear();

    breed();
    if (conv) {
      calculate_pop_stats();
      return conv->evaluate_convergence(this->get_n_objs(), this->pop_stats);
    } else {
      generation++;
      return (generation > this->args.get_num_gens());
    }
  }

  _uint get_n_pareto_fronts() {
    return pareto_fronts.size();
  }
  std::vector<std::shared_ptr< Organism<FitType> >> get_pareto_front(_uint i) {
    if (pareto_fronts.size() == 0) {
error(CODE_MISC, "The population has not yet been evaluated.");
    }
    if ( i >= pareto_fronts.size() ) {
      error(CODE_ARG_RANGE, "Attempt to access invalid index %d out of %d fronts.",
      i, pareto_fronts.size());
    }
    return pareto_fronts[i];
  }

  Vector<String> get_best_header() {
    //TODO: create an actual implementation
    Vector<String> ret;
    return ret;
  }
  Vector<String> get_header() {
    Vector<String> ret;
    char buf[OUT_BUF_SIZE];
    for (_uint i = 0; i < this->old_gen.size(); ++i) {
      snprintf(buf, OUT_BUF_SIZE, "rank_%d", i);
      ret.push_back( String(buf) );

      for (_uint j = 0; j < this->map->get_num_params(); ++j) {
        ret.push_back( String(this->var_labels[j]) );
      }

      for (_uint j = 0; j < this->get_n_objs(); ++j) {
        ret.push_back( String(this->obj_labels[j]) );
      }
    }
    return ret;
  }
  Vector<String> get_best_data() {
    //TODO: create an actual implementation
    Vector<String> ret;
    return ret;
  }
  Vector<String> get_pop_data() {
    Vector<String> ret;
    bool make_front = true;
    for (_uint i = 0; i < this->offspring_num; ++i) {
      if (this->offspring[i].get() == NULL) {
        make_front = false;
        break;
      }
    }
    if (make_front) {
      make_fronts( &(this->offspring) );
      String tmp_str;
      char buf[OUT_BUF_SIZE];
      for (_uint n = 0; n < pareto_fronts.size(); ++n) {
        for (_uint i = 0; i < pareto_fronts[n].size(); ++i) {
          snprintf(buf, OUT_BUF_SIZE - 1, "%d", n);
          tmp_str = buf;
          ret.push_back(tmp_str);

          for (_uint j = 0; j < this->map->get_num_params(); ++j) {
            ret.push_back( pareto_fronts[n][i]->get_chromosome_string(j) );
          }

          for (_uint j = 0; j < this->get_n_objs(); ++j) {
            snprintf(buf, OUT_BUF_SIZE - 1, "%f", pareto_fronts[n][i]->get_fitness(j));
            ret.push_back(String(buf));
          }
        }
      }
    }
    return ret;
  }

  std::shared_ptr< Organism<FitType> > get_organism(size_t i) {
    return (i < this->offspring_num) ? this->old_gen[i] : this->offspring[i - this->offspring_num];
  }

  size_t get_offspring_num() { return this->offspring_num; }
  size_t get_survivors_num() { return this->survivors_num; }
};

}

#endif //POPULATION_H
