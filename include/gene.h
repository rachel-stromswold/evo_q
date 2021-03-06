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
#define REAL_ENABLED	1
#define REAL_ACTIVE	2
#define DISC_DISABLED	4

namespace Genetics {

class Chromosome {
private:
  _uint N_BITS;
  _uint N_BYTES;
  size_t N;
  double* real_vals = NULL;
  size_t real_vals_size = 0;

protected:
  static const _uint bin_size = sizeof(unsigned long)*8;
//  unsigned long genes[(N_BYTES+sizeof(unsigned long)-1)/sizeof(unsigned long)];
  unsigned long* genes = NULL;
  size_t getBitStream (size_t n, size_t k, size_t x);
  _uchar use_real = 0;
  size_t get_real_vals_size() { return real_vals_size; }

public:
  Chromosome(_uint pn_bits);
  Chromosome(_uint pn_bits, _uchar real_mode);
  Chromosome(_uint pn_bits, Chromosome* o);
  Chromosome(Chromosome& other);
  Chromosome(Chromosome&& other);
  ~Chromosome();
  void exchange(Chromosome* other, size_t k);
  void exchange_uniform(ArgStore& args, Chromosome* other);

  unsigned int get_N() { return N; }
  unsigned int get_n_bits() { return N_BITS; }
  void swap(Chromosome& other);
  Chromosome& operator=(Chromosome& other);
  Chromosome& operator=(Chromosome&& other);
  bool operator==(Chromosome& other);

  void reset();

  unsigned char operator[](unsigned int i);
  //randomly mutate each bit in the gene
  bool real_space_mutate(ArgStore& args);
  void mutate(ArgStore& args);
  void slow_mutate(ArgStore& args);
  //set the gene to a new completely random value
  void randomize(PhenotypeMap* al, ArgStore& args);
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
  Vector<double> get_real_vector(PhenotypeMap* al);
};

}

#endif //GENE_H
