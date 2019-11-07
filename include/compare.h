#ifndef COMPARE_H
#define COMPARE_H

#define DEF_SORT_PARAM	-3

#define FLAG_FRONTS	8

#include "organism.h"

namespace Genetics {

struct FitnessStats {
  double mean;
  double var;
  double max;
  double min;
};

/*struct OrganismPair {
  std::shared_ptr<Organism> first, second;
};*/

template <typename FitType>
class Comparator {
public:
  static int compare(std::shared_ptr< Organism<FitType> > a, std::shared_ptr< Organism<FitType> > b) {
    FitType a_info = a->get_fitness_info();
    FitType b_info = a->get_fitness_info();
    if ( a->get_n_objs() != b->get_n_objs() ) {
      error(CODE_MISC, "Comparison of fitness values with a different number of objectives.");
    }
    int ret = 0;
    for (_uint i = 0; i < a_info.get_n_objs(); ++i) {
      if (a->get_fitness() < b->get_fitness()) {
        --ret;
      } else if (a->get_fitness() > b->get_fitness()) {
        ++ret;
      }
    }
    return true;
  }
};

template <typename FitType>
class NSGAII_Comparator : public Comparator<FitType> {
public:
  static_assert( std::is_base_of<MultiFitness, FitType>::value, "FitType must be derived from MultiFitness for NSGAII comparator" );
  static int compare(std::shared_ptr< Organism<FitType> > a, std::shared_ptr< Organism<FitType> > b) {
    for (_uint i = 0; i < a->get_fitness_info().get_n_objs(); ++i) {
      if (a->get_fitness(i) <= b->get_fitness(i)) {
        return 0;
      }
    }
    return 1;
  }
};

typedef std::pair<_uint, _uint> ParentIndSet;

template <typename FitType, typename Comparator=Comparator<FitType>>
class Selector {
protected:
  int partition(_uint fit_ind, std::vector<std::shared_ptr< Organism<FitType> >>& work_arr, int s, int e) {
    //double p = (*work_arr)[e]->get_fitness(fit_ind);
    int i = s;
    for (int j = s; j < e; ++j) {
      if ( Comparator::compare( work_arr[j], work_arr[e] ) > 0) {
        std::shared_ptr<Organism<FitType>> tmp = work_arr[i];
        work_arr[i] = work_arr[j];
        work_arr[j] = tmp;
        ++i;
      }
    }
    if (i != e) {
      std::shared_ptr<Organism<FitType>> tmp = work_arr[i];
      work_arr[i] = work_arr[e];
      work_arr[e] = tmp;
    }
    return i;
  }
public:
  void sort_orgs(unsigned int fit_ind, std::vector<std::shared_ptr< Organism<FitType> >>& work_arr, int s = DEF_SORT_PARAM, int e = DEF_SORT_PARAM) {
    if (s == DEF_SORT_PARAM || e == DEF_SORT_PARAM) {
      sort_orgs(fit_ind, work_arr, 0, work_arr.size() - 1);
    } else if (s < e) {
      int p = partition(fit_ind, work_arr, s, e);
      sort_orgs(fit_ind, work_arr, s, p-1);
      sort_orgs(fit_ind, work_arr, p+1, e);
    }
  }
  Selector() {
    static_assert( std::is_base_of<Fitness, FitType>::value, "FitType must be derived from Fitness" );
  }
  //virtual OrganismPair select(Population& pop);
  virtual ~Selector() = default;
  //returns >=1 if a > b, 0 if a = b and -1 if a <= b
  //errors if the two fitness values are not comprable (they do not have the same number of objectives)
  
  virtual Vector<ParentIndSet> select(ArgStore& args, Vector<std::shared_ptr<Organism<FitType>>>& old_gen, Vector<std::shared_ptr<Organism<FitType>>>& offspring) = 0;
};

template <class FitType>
class TournamentSelector : public Selector<FitType> {
public:
  typedef Comparator<FitType> Comparator;

  Vector<ParentIndSet> select(ArgStore& args, Vector<std::shared_ptr<Organism<FitType>>>& old_gen, Vector<std::shared_ptr<Organism<FitType>>>& offspring) {
    _uint arena_size = args.read_custom_double("arena_size", 2);
    if (arena_size < 2) { arena_size = 2; }
    bool tournament_replacement = (args.get_custom_parameter("tournament_replacement") != "");
    _uint offspring_num = old_gen.size();
    SampleDraw sampler(offspring_num, arena_size, tournament_replacement);
    std::uniform_int_distribution<_uint> selector(0, offspring_num - 1);
    std::vector<Organism<FitType>*> children;
    Vector<ParentIndSet> ret( divideup(offspring_num, 2) );

    for (size_t i = 0; 2*i + 1 < offspring_num; ++i) {
      std::shared_ptr< Organism<FitType> > first_parent, second_parent;
      std::vector<_uint> t1 = sampler( args.get_generator() );
      std::vector<_uint> t2 = sampler( args.get_generator() );
      ret[i].first  = t1[0];
      ret[i].second = t2[0];
      for (size_t j = 1; j < t1.size(); ++j) {
        //check whether the fitness is an improvement and use variance as a tiebreaker
        if (Comparator::compare(old_gen[t1[j]], old_gen[ret[i].first]) > 0) {
          ret[i].first  = t1[j];
        }
        if (Comparator::compare(old_gen[t2[j]], old_gen[ret[i].second]) > 0) {
          ret[i].second = t2[j];
        }
      }
      //ensure that we don't use the same parent twice with reasonable probability
      if (ret[i].first == ret[i].second) {
        ret[i].second = selector( args.get_generator() );
      }
    }
    return ret;
  }
};

template <class FitType>
class SurvivalSelector : public Selector<FitType> {
public:
  typedef Comparator<FitType> Comparator;

  Vector<ParentIndSet> select(ArgStore& args, Vector<std::shared_ptr<Organism<FitType>>>& old_gen, Vector<std::shared_ptr<Organism<FitType>>>& offspring) {
    sort_orgs(0, old_gen);
    /*TODO: determine whether we actually need to figure out a way to keep this in place
     * //if all the organisms have the same fitness then reinitialize the population
    if (old_gen[0]->get_fitness(0) - old_gen[offspring_num-1]->get_fitness(0) < 0.001 ) {
      return true;
    }*/
    size_t offspring_num = old_gen.size();
    double min_fit = old_gen[offspring_num-1]->get_fitness(0);
    double total_fit = 0;
    for (size_t i = 0; i < offspring_num; ++i) {
      total_fit += old_gen[i]->get_fitness(0) - min_fit;
    }
    std::uniform_real_distribution<double> dist(0, total_fit);

    _uint survivors_num = read_double(args.get_custom_parameter("num_survivors"), offspring_num/2);
    std::vector< std::shared_ptr<Organism<FitType>> > survivors(survivors_num);
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

    std::uniform_int_distribution<size_t> dist_surv0(0, survivors_num - 1);
    std::uniform_int_distribution<size_t> dist_surv1(0, survivors_num - 2);
    Vector<ParentIndSet> ret( divideup(offspring_num, 2) );
    for (_uint i = 0; 2*i + 1< offspring_num; ++i) {
      ret[i].first  = dist_surv0( args.get_generator() );
      ret[i].second = dist_surv1( args.get_generator() );
      if (ret[i].second >= ret[i].first) {
        ++ret[i].second;
      }
    }
    return ret;
  }
};

template <typename FitType>
class NSGAII_TournamentSelector : public Selector<FitType, NSGAII_Comparator<FitType>> {
public:
  typedef std::shared_ptr< Organism<FitType> > OrgPtr;
  typedef NSGAII_Comparator<FitType> Comparator;
private:  
  _uint n_objs = 1;
  Vector<Vector< std::shared_ptr<Organism<FitType>> >> pareto_fronts;
  Vector<FitnessStats> pop_stats;
  _uchar calculated_flags = 0;
  size_t min_penalty_ind, max_penalty_ind;
  void make_fronts(std::vector<OrgPtr>& cmb_arr) {
    pareto_fronts.clear();
    std::vector<OrgPtr> empty;
    pareto_fronts.push_back(empty);

    for (size_t i = 0; i < cmb_arr.size(); ++i) {
      cmb_arr[i]->get_fitness_info().n_dominations = 0;
      for (_uint j = 0; j < n_objs; ++j) {
        if (cmb_arr[i]->get_fitness(j) > pop_stats[j].max) {
          pop_stats[j].max = cmb_arr[i]->get_fitness(j);
        }
        if (cmb_arr[i]->get_fitness(j) < pop_stats[j].min) {
          pop_stats[j].min = cmb_arr[i]->get_fitness(j);
        }
      }
      for (size_t j = 0; j < cmb_arr.size(); ++j) {
        if (i != j) {
          //if the ith solution dominates the jth add the jth entry to the list of dominated solutions, otherwise increment the number of dominating solutions
          if ( Comparator::compare(cmb_arr[j], cmb_arr[i]) > 0 ) {
            cmb_arr[i]->get_fitness_info().n_dominations++;
          }
        }
      }
      if ( cmb_arr[i]->get_fitness_info().n_dominations == 0) {
        cmb_arr[i]->get_fitness_info().rank = 0;
        pareto_fronts[0].push_back(cmb_arr[i]);
      }
    }

    size_t i = 0;
    while (i < pareto_fronts.size() && pareto_fronts[i].size() != 0) {
      pareto_fronts.push_back(empty);
      for (size_t j = 0; j < pareto_fronts[i].size(); ++j) {
        for (size_t k = 0; k < cmb_arr.size(); ++k) {
          cmb_arr[k]->get_fitness_info().n_dominations--;
          if (cmb_arr[k]->get_fitness_info().n_dominations == 0) {
            cmb_arr[k]->get_fitness_info().rank = i + 1;
            pareto_fronts[i + 1].push_back(cmb_arr[k]);
          }
        }
      }
      i++;
    }
    calculated_flags |= FLAG_FRONTS;
  }
 
  Vector<ParentIndSet> gen_breed_pairs(ArgStore& args, Vector<OrgPtr>& old_gen, Vector<OrgPtr>& offspring) {
    _uint arena_size = args.read_custom_double("arena_size", 2);
    size_t offspring_num = old_gen.size();
    SampleDraw sampler(offspring_num, arena_size);
    std::uniform_int_distribution<_uint> selector(0, offspring_num - 1);
    Vector<ParentIndSet> ret( divideup(offspring_num, 2) );

    for (size_t i = 0; 2*i + 1 < offspring_num; ++i) {
      std::vector<_uint> t1 = sampler( args.get_generator() );
      std::vector<_uint> t2 = sampler( args.get_generator() );
      ret[i].first = t1[0];
      ret[i].second = t2[0];
      for (size_t j = 1; j < t1.size(); ++j) {
        if (old_gen[t1[j]]->get_fitness_info().rank < old_gen[ret[i].first]->get_fitness_info().rank) {
          ret[i].first = t1[j];
        }
        if (old_gen[t2[j]]->get_fitness_info().rank < old_gen[ret[i].second]->get_fitness_info().rank) {
          ret[i].second = t2[j];
        }
      }
      //ensure that we don't use the same parent twice with reasonable probably
      if (ret[i].first == ret[i].second) {
        ret[i].second = selector( args.get_generator() );
      }
    }
    return ret;
  }
public:
  NSGAII_TournamentSelector() : Selector<FitType, NSGAII_Comparator<FitType>>() {
    static_assert( std::is_base_of<MultiFitness, FitType>::value, "FitType must be derived from MultiFitness for the NSGAII selector" );
  }
  
  Vector<ParentIndSet> select(ArgStore& args, Vector<OrgPtr>& old_gen, Vector<OrgPtr>& offspring) {
    if (n_objs != old_gen[0]->get_fitness_info().get_n_objs()) {
      n_objs = old_gen[0]->get_fitness_info().get_n_objs();
      pop_stats.resize(n_objs);
    }
    if (pareto_fronts.size() == 0 && !(calculated_flags & FLAG_FRONTS)) {
      std::vector<OrgPtr> cmb_arr = old_gen;
      cmb_arr.reserve(old_gen.size() + offspring.size());
      for (_uint i = 0; i < offspring.size(); ++i) {
        if (offspring[i]) {
          cmb_arr.push_back(offspring[i]);
        }
      }
      make_fronts(cmb_arr);
    }
    size_t offspring_num = old_gen.size();
    size_t i = 0;
    std::vector<OrgPtr> tmp(offspring_num, NULL);
    _uint k = 0;
    while (i < pareto_fronts.size() && (k + pareto_fronts[i].size()) <= offspring_num) {
      for (size_t j = 0; j < pareto_fronts[i].size(); ++j) {
        tmp[k] = pareto_fronts[i][j];
        ++k;
      }
      ++i;
    }
    //fill in the last elements from the remaining pareto front ranked according to crowding
    if (k < offspring_num) {
      for (size_t ii = 0; ii < pareto_fronts[i].size(); ++ii) {
        std::cout << "i=" << i << ", ii=" << ii << ", k=" << k << std::endl;
        pareto_fronts[i][ii]->get_fitness_info().distance = 0;
      }
      //sort by each objective function for crowding evaluation
      for (size_t j = 0; j < n_objs; ++j) {
        this->sort_orgs(j, pareto_fronts[i]);

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
      this->sort_orgs(n_objs, pareto_fronts[i]);
      size_t j = 0;
      //select the least crowded individuals
      size_t p_i_size = pareto_fronts[i].size();
      while (j < p_i_size && k < offspring_num) {
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
    old_gen = tmp;
    //ensure that we aren't keeping null pointers around for safety
    pareto_fronts.clear();

    return gen_breed_pairs(args, old_gen, offspring);
  }
};

}

#endif
