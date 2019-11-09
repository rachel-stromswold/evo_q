#include "compare.h"

namespace Genetics {

template <class FitType>
bool TournamentSelector<FitType>::compare(Organism<FitType>& a, Organism<FitType>& b) {
  FitType fit_a = a.get_fitness_stats();
  FitType fit_b = a.get_fitness_stats();
  for (_uint i = 0; i < fit_a.get_n_objs(); ++i) {
    if ( fit_a.get_fitness(i) < fit_b.get_fitness(i) ) {
      return false;
    }
  }
  return true;
}

}
