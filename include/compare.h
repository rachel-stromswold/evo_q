#ifndef COMPARE_H
#define COMPARE_H

#include "organism.h"

namespace Genetics {

/*struct OrganismPair {
  std::shared_ptr<Organism> first, second;
};*/

template <class FitType>
class Selector {
public:
  Selector() { static_assert( std::is_base_of<Fitness, FitType>::value, "FitType must be derived from Fitness" ); }
  //virtual OrganismPair select(Population& pop);
  virtual ~Selector() = default;
  virtual bool compare(std::shared_ptr< Organism<FitType> > a, std::shared_ptr< Organism<FitType> > b) = 0;
  virtual void select(ArgStore* args, Vector<std::shared_ptr<Organism<FitType>>>& old_gen, Vector<std::shared_ptr<Organism<FitType>>>& offspring) = 0;
};

template <class FitType>
class TournamentSelector : public Selector<FitType> {
public:
  bool compare(std::shared_ptr< Organism<FitType> > a, std::shared_ptr< Organism<FitType> > b) {
    FitType a_info = a->get_fitness_info();
    FitType b_info = a->get_fitness_info();
    if ( a_info.get_n_objs() != b_info.get_n_objs() ) {
      return false;
    }
    for (_uint i = 0; i < a_info.get_n_objs(); ++i) {
      if (a->get_fitness() < b->get_fitness()) {
        return false;
      }
    }
    return true;
  }

  void select(ArgStore* args, Vector<std::shared_ptr<Organism<FitType>>>& old_gen, Vector<std::shared_ptr<Organism<FitType>>>& offspring) {
    _uint arena_size = read_double(args->get_custom_parameter("arena_size"), 2);
    if (arena_size < 2) { arena_size = 2; }
    bool tournament_replacement = (args->get_custom_parameter("tournament_replacement") != "");
    _uint offspring_num = old_gen.size();
    if (offspring.size() != old_gen.size()) {
      error(CODE_MISC, "size of offspring and old generation differ in call to select");
    }
    SampleDraw sampler(offspring_num, arena_size, tournament_replacement);
    std::uniform_int_distribution<_uint> selector(0, offspring_num - 1);
    std::vector<Organism<FitType>*> children;

    for (size_t i = 0; 2*i + 1 < offspring_num; ++i) {
      std::shared_ptr< Organism<FitType> > first_parent, second_parent;
      std::vector<_uint> t1 = sampler( args->get_generator() );
      std::vector<_uint> t2 = sampler( args->get_generator() );
      first_parent = old_gen[t1[0]];
      second_parent = old_gen[t2[0]];
      for (size_t j = 1; j < t1.size(); ++j) {
        //check whether the fitness is an improvement and use variance as a tiebreaker
        if ( compare(old_gen[t1[j]], first_parent) ) {
          first_parent = old_gen[t1[0]];
        }
        if ( compare(old_gen[t2[j]], second_parent) ) {
          second_parent = old_gen[t2[0]];
        }
      }
      //ensure that we don't use the same parent twice with reasonable probability
      if (first_parent == second_parent) {
        second_parent = old_gen[selector( args->get_generator() )];
      }

      children = first_parent->breed(args, second_parent.get());
      offspring[2*i] = std::shared_ptr<Organism<FitType>>(children[0]);
      offspring[2*i + 1] = std::shared_ptr<Organism<FitType>>(children[1]);
    }
  }
};

}

#endif
