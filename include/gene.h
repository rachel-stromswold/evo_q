#ifndef GENE_H
#define GENE_H

#include <stdlib.h>
#include <iostream>
#include <random>
#include <climits>
#include <type_traits>

#include "parse.h"
#include "phenotype.h"

#define MUTATE_PROB	0.1

namespace Genetics {

class Chromosome {
private:
  _uint N_BITS;
  size_t N_BYTES;
  size_t N;

protected:
  static const _uint bin_size = sizeof(unsigned long)*8;
//  unsigned long genes[(N_BYTES+sizeof(unsigned long)-1)/sizeof(unsigned long)];
  Vector<unsigned long> genes;
  size_t getBitStream (size_t n, size_t k, size_t x);

public:
  Chromosome(_uint pn_bits);
  Chromosome(_uint pn_bits, Chromosome* o);
  Chromosome(Chromosome& other);
  Chromosome(Chromosome&& other);
  void exchange(Chromosome* other, size_t k);
  void exchange_uniform(ArgStore* args, Chromosome* other);

  unsigned int get_N() { return N; }
  unsigned int get_n_bits() { return N_BITS; }
  Chromosome& operator=(Chromosome& other);

  void reset();

  unsigned char operator[](unsigned int i);
  //randomly mutate each bit in the gene
  void mutate(ArgStore* args);
  void slow_mutate(ArgStore* args);
  //set the gene to a new completely random value
  void randomize(ArgStore* args);
  //sets the gene to encode the value specified by min, max
  void set_to_num(PhenotypeMap* al, _uint ind, double value);
  void set_to_int(PhenotypeMap* al, _uint ind, int value);
  void set_to_ulong(PhenotypeMap* al, _uint ind, unsigned long value);
  //returns the corresponding integer for the gene
  unsigned long gene_to_ulong(PhenotypeMap* al, _uint ind);
  int gene_to_int(PhenotypeMap* al, _uint ind);
  //returns a double value corresponding to the gene, it will have a value between max and min
  double gene_to_num(PhenotypeMap* al, _uint ind);
  String get_string(PhenotypeMap* al, _uint ind);
};

}

#endif //GENE_H
