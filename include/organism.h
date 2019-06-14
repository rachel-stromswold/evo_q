#ifndef ORGANISM_H
#define ORGANISM_H

#include "gene.h"
//#include "util.h"
#include <fstream>
#include <sstream>
#include <string.h>

#include <unistd.h>
#include <sys/wait.h>

#define COUPLING_BITS	16
#define M_COUP		10//The maximum value that a coupling parameter may take
#define BUF_SIZE	50
#define N_CROSSOVERS  1

#define COUP_VAR	20.0
#define COUP_MEAN	0.0

namespace Genetics {

class Organism;

struct Result {
  size_t index;
  Vector<double> fit_vals;
  String misc_str;
};

class Problem {
public:
  _uint N_BITS, N_PARAMS, N_OBJS;
  PhenotypeMap map;
  Vector<Result> result_list;
//  Problem() {std::cout << "Initializing problem...\n"; }
  Problem(unsigned n_bits, unsigned n_params, int n_objs) : N_BITS(n_bits), N_PARAMS(n_params), N_OBJS(n_objs), map(n_bits) {}

  virtual void evaluate_fitness(Organism* org) {}
  virtual void evaluate_fitness_async(size_t index, Chromosome genes) {}
};

class Organism {
private:
  _uint N_BITS;
  _uint N_OBJS;

  char output_stream[BUF_SIZE];
  Vector<double> fitness;
  size_t output_len;
  PhenotypeMap* al;

protected:
  Chromosome* genes;
  size_t n_nodes;

public:
  //this is not used internally, but can be set when evaluating the fitness
  String misc_data;
  double coupling_range;
  double coupling_prec;

  int n_dominations;
  int rank;
  double distance;

  Organism(int N_BITS, int N_OBJS, PhenotypeMap* p_al);
  Organism(int N_BITS, int N_OBJS, Chromosome p_genes, PhenotypeMap* p_al);
  Organism(const Organism &obj);
  Organism(Organism&& obj);
  ~Organism();

  Organism& operator=(Organism& obj);

  std::vector<Organism*> breed(ArgStore* args, Organism* par1);
  void reset();
  void randomize(ArgStore* args);
  void randomize(ArgStore* args, Organism* orgtmp);

  void evaluate_fitness(Problem* prob);

  double get_fitness(_uint i = 0);
  void set_fitness(double val);
  void set_fitness(_uint i, double val);

  void set_int(_uint i, int value);
  void set_uint(_uint i, int value);
  void set_real(_uint i, double value);
  double read_real(_uint i);
  int read_int(_uint i);
  bool dominates(Organism* other);
  String get_chromosome_string(_uint i) { return genes->get_string(al, i); }
  char* get_output_stream() { return output_stream; }
  size_t get_output_len() {return output_len; }
  int get_rank() { return rank; }
  _uint get_n_bits() { return N_BITS; }
  _uint get_n_params() { return al->get_num_params(); }
  _uint get_n_objs() { return N_BITS; }
};

}

#endif //ORGANISM_H
