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
};

}

#endif
