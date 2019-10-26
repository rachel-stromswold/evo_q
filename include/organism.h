#ifndef ORGANISM_H
#define ORGANISM_H

#include "gene.h"
//#include "util.h"
#include <fstream>
#include <sstream>
#include <string.h>
#include <memory>

#include <unistd.h>
#include <sys/wait.h>

#define COUPLING_BITS	16
#define M_COUP		10//The maximum value that a coupling parameter may take
#define BUF_SIZE	50
#define N_CROSSOVERS  1

#define COUP_VAR	20.0
#define COUP_MEAN	0.0

namespace Genetics {

template <class FitType>
class Organism;

template <class FitType>
class Problem {
public:
  _uint N_BITS, N_PARAMS, N_OBJS;
  std::shared_ptr<PhenotypeMap> map;
  //Vector<Result> result_list;
//  Problem() {std::cout << "Initializing problem...\n"; }
  Problem(unsigned n_bits, unsigned n_params, int n_objs) : N_BITS(n_bits), N_PARAMS(n_params), N_OBJS(n_objs), map(std::make_shared<PhenotypeMap>(n_bits)) {}

  virtual void evaluate_fitness(Organism<FitType>* org) {}
  virtual void evaluate_fitness_async(size_t index, Chromosome genes) {}
};

/**
 * \brief	The abstract class FitnessStats is used by implementations of the Selector abstract class to select which organism is more fit
 */
class FitnessStats {
protected:
  _uint N_OBJS = 1;
public:
  FitnessStats(_uint pn_objs = 1) { N_OBJS = pn_objs; }
  virtual double get_fitness(_uint i = 0);
  double get_cost() { return -get_fitness(); }
  virtual ~FitnessStats() = default;
  virtual void update(double val, _uint i = 0);
  void reset() {}
  _uint get_n_objs() { return N_OBJS; }
};

/**
 * \brief	An implementation of FitnessStats designed for single-objective optimization in the noise-free case
 */
class SingleFitness : public FitnessStats {
protected:
  double fitness;
public:
  SingleFitness();
  double get_fitness(_uint i = 0);
  virtual void update(double val, _uint i = 0);
};

/**
 * \brief	An implementation of FitnessStats designed for multi-objective optimization in the noise-free case
 */
class MultiFitness : public FitnessStats {
protected:
  Vector<double> fitness;
  _uint n_dominations = 0;
  _uint rank = 0;
  double distance = 0;

public:
  MultiFitness(_uint pn_objs);
  double get_fitness(_uint i);
  void update(double val, _uint i);
};

/**
 * \brief	An implementation of FitnessStats designed for single-objective optimization in the noisy case
 */
class NoisyFitness : public FitnessStats {
protected:
  double fitness = 0, variance = 0;

  _uint n_evaluations = 0;

public:
  NoisyFitness();
  double get_fitness(_uint i = 0);
  virtual void update(double val, _uint i = 0);
};

/**
 * \brief	An implementation of FitnessStats designed for single-objective optimization in the noisy case. This implementation uses a parameter forget_weight that biases results to more heavily weigh recent observations
 */
class NoisyFitnessForgetful : public FitnessStats {
protected:
  double fitness = 0, variance = 0;
  double forget_weight;
  bool evaluated = false;
public:
  NoisyFitnessForgetful(double p_forget_weight);
  double get_fitness(_uint i = 0);
  void update(double val, _uint i = 0);
};

/**
 * \brief	An implementation of FitnessStats designed for single-objective optimization in the noisy case
 */
class NoisyMultiFitness : public FitnessStats {
protected:
  Vector<double> fitness;
  Vector<double> variances;
  _uint N_OBJS;
  _uint n_evaluations;

public:
  NoisyMultiFitness(_uint pn_objs);
  double get_fitness(_uint i);
  void update(double val, _uint i);
};

template <class FitType>
class Organism {
  static_assert( std::is_base_of<FitnessStats, FitType>::value, "FitType must be derived from FitnessStats" );
private:
  _uint N_BITS;
  _uint N_OBJS;

  char output_stream[BUF_SIZE];
  double penalty = 0.0;
  size_t output_len;
  std::shared_ptr<PhenotypeMap> al;

protected:
  Chromosome genes;
  size_t n_nodes;
  FitType fit;
  //Vector<double> fitness;
  //Vector<double> fit_vars;
  //int n_evaluations = 0;

public:
  //this is not used internally, but can be set when evaluating the fitness
  String misc_data;
  double coupling_range;
  double coupling_prec; 

  Organism();
  Organism(int N_BITS, int N_OBJS, PhenotypeMap* p_al);
  Organism(int N_BITS, int N_OBJS, Chromosome p_genes, PhenotypeMap* p_al);
  Organism(int N_BITS, int N_OBJS, std::shared_ptr<PhenotypeMap> p_al);
  Organism(int N_BITS, int N_OBJS, Chromosome p_genes, std::shared_ptr<PhenotypeMap> p_al);
  /**
   * \brief	Update the FitnessStats object to match the template temp. This should be used for parameter setting before any evaluation calls have been made.
   */
  void set_fitness_stats(FitType temp);
  //Organism(const Organism &obj);
  //Organism(Organism&& obj);
  //~Organism();
  Organism copy();

  //Organism& operator=(Organism& obj);
  bool operator==(Organism& obj);
  bool operator!=(Organism& obj);
  bool operator>(Organism& obj);
  bool operator<(Organism& obj);

  void swap(Organism& obj);
  bool valid() { return (al != NULL && N_OBJS > 0 && N_BITS > 0); }

  std::vector<Organism*> breed(ArgStore* args, Organism* par1);
  void mutate(ArgStore* args);
  void reset();
  void randomize(ArgStore* args);
  void randomize(ArgStore* args, Organism* orgtmp);

  //void evaluate_fitness_noisy(Problem<FitType>* prob, double forget_weight=0);
  void evaluate_fitness(Problem<FitType>* prob);

  FitType get_fitness_stats();
  double get_fitness(_uint i = 0);
  double get_cost(_uint i = 0);
  void set_fitness(double val);
  void set_cost(double val);
  void apply_penalty(double val) { penalty = val; }

  double get_penalty() { return penalty; }
  bool penalized() { return penalty != 0; }
  void set_fitness(_uint i, double val);
  void set_cost(_uint i, double val);
  //_uint get_n_evaluations() { return n_evaluations; }
  //average fitness with another organism if they both have the same genotype
  void average_fitness(Organism* other);
  void copy_fitness_data(Organism* other);

  void set_int(_uint i, int value);
  void set_uint(_uint i, int value);
  void set_real(_uint i, double value);
  double read_real(_uint i);
  int read_int(_uint i);
  _uint read_uint(_uint i);
  bool dominates(Organism* other);
  String get_chromosome_string(_uint i);
  char* get_output_stream() { return output_stream; }
  size_t get_output_len() {return output_len; }
  //int get_rank() { return rank; }
  _uint get_n_bits() { return N_BITS; }
  _uint get_n_params() { return al->get_num_params(); }
  _uint get_n_objs() { return N_OBJS; }
};

}

#endif //ORGANISM_H
