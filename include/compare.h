#ifndef COMPARE_H
#define COMPARE_H

namespace Genetics {

struct OrganismPair {
  std::shared_ptr<Organism> first, second;
};

class Selector {
public:
  virtual OrganismPair select(Population& pop);
  virtual ~Selector() = default;
};

class TournamentSelector : public Selector {

};

}

#endif
