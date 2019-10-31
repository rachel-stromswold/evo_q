#ifndef ORGANISM_H
#define ORGANISM_H

#include "gene.h"
//#include "util.h"
#include <fstream>
#include <sstream>
#include <string.h>
#include <memory>
#include <type_traits>

#include <unistd.h>
#include <sys/wait.h>

#define COUPLING_BITS	16
#define M_COUP		10//The maximum value that a coupling parameter may take
#define BUF_SIZE	50
#define N_CROSSOVERS  1

#define COUP_VAR	20.0
#define COUP_MEAN	0.0

namespace Genetics {

/**
 * \brief	The abstract class
 */
/*class Fitness {
protected:
  _uint N_OBJS;
public:
  void reset() {}
  _uint get_n_objs() { return N_OBJS; }
  virtual double get_cost() { return -get_fitness(); }
  double get_uncertainty(_uint i = 0) { return 0.0; }

  virtual ~Fitness() = default;
  virtual double get_fitness(_uint i = 0) = 0;
  virtual void update(double val, _uint i = 0) = 0;
};*/


/**
 * \brief	An implementation of Fitness designed for single-objective optimization in the noise-free case. Other single objective problems can use Fitness classes that inherit from SingleFitness for generation of comparisons.
 */
class Fitness {
protected:
  double fitness;
public:
  double distance = 0;

  Fitness() { fitness = 0; }

  void reset() {}
  void update(double val, _uint i = 0) { fitness = val; }
  _uint get_n_objs() { return 1; }
  double get_fitness(_uint i = 0) { return fitness; }
  double get_cost(_uint i = 0) { return -get_fitness(i); }
  double get_uncertainty(_uint i = 0) { return 0.0; }
};
typedef Fitness SingleFitness;

/**
 * \brief	The class MultiFitness is used by implementations of the Selector abstract class to select which organism is more fit when multiple objectives are to be considered
 */
class MultiFitness : public Fitness {
protected:
  _uint N_OBJS;
  Vector<double> fitness;
 
public:
  _uint n_dominations = 0;
  _uint rank = 0;

  MultiFitness(_uint pn_objs = 1) : fitness(pn_objs, 0.0) { N_OBJS = pn_objs; }

  void reset() {}
  void update(double val, _uint i);

  _uint get_n_objs() { return N_OBJS; }
  double get_fitness(_uint i);
  double get_cost(_uint i) { return -get_fitness(i); }
  double get_uncertainty(_uint i) { return 0.0; }
  _uint get_rank() { return rank; }
  double get_distance() { return distance; }
  _uint get_n_dominations() { return n_dominations; }
};

/**
 * \brief	An implementation of Fitness designed for single-objective optimization in the noisy case
 */
class NoisyFitness : public SingleFitness {
protected:
  double fitness = 0, variance = 0;
  _uint n_evaluations = 0;

public:
  _uint get_n_evaluations() { return n_evaluations; }
  double get_uncertainty(_uint i = 0) { return sqrt(variance); }
  void update(double val, _uint i = 0);
  void average_fitness(NoisyFitness* other);
};

/**
 * \brief	An implementation of Fitness designed for single-objective optimization in the noisy case. This implementation uses a parameter forget_weight that biases results to more heavily weigh recent observations
 */
class NoisyFitnessForgetful : public NoisyFitness {
protected:
  double fitness = 0, variance = 0;
  double forget_weight;
  bool evaluated = false;
public:
  NoisyFitnessForgetful(double p_forget_weight = 1.0);
  void update(double val, _uint i = 0);
};

/**
 * \brief	An implementation of Fitness designed for single-objective optimization in the noisy case
 */
class NoisyMultiFitness : public Fitness {
protected:
  Vector<double> fitness;
  Vector<double> variances;
  _uint N_OBJS;
  _uint n_evaluations;

public:
  NoisyMultiFitness(_uint pn_objs = 1);
  double get_fitness(_uint i);
  void update(double val, _uint i);
  double get_uncertainty(_uint i) { return sqrt(variances[i]); }
};

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

  virtual void evaluate_fitness(Organism<FitType>* org) = 0;
#ifdef USE_LIBOMP
  void evaluate_fitness_async(Organism<FitType>* org, _uint i = 0) { evaluate_fitess(org); }
#endif
};

template <class FitType>
class Organism {
  static_assert( std::is_base_of<Fitness,FitType>::value, "FitType must be derived from SingleFitness or MultiFitness" );
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

  Organism() : genes(0) { N_BITS = 0;N_OBJS = 0; }
  Organism(int pn_bits, int pn_objs, PhenotypeMap* p_al) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs),
  genes(N_BITS),
  al(p_al)
  {
    memset(output_stream, 0, BUF_SIZE);
    reset_fitness();
  }
  Organism(int pn_bits, int pn_objs, Chromosome p_genes, PhenotypeMap* p_al) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs),
  genes(p_genes),
  al(p_al)
  {
    reset_fitness();
  }
  Organism(int pn_bits, int pn_objs, std::shared_ptr<PhenotypeMap> p_al) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs),
  genes(N_BITS),
  al(p_al)
  {
    memset(output_stream, 0, BUF_SIZE);
    reset_fitness();
  }
  Organism(int pn_bits, int pn_objs, Chromosome p_genes, std::shared_ptr<PhenotypeMap> p_al) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs),
  genes(p_genes),
  al(p_al)
  {
    reset_fitness();
  }
  /**
   * \brief	Update the Fitness object to match the template temp. This should be used for parameter setting before any evaluation calls have been made.
   */
  void set_fitness_stats(FitType temp) { fit = temp; }
  //Organism(const Organism &obj);
  //Organism(Organism&& obj);
  //~Organism();
  Organism copy() {
    Organism<FitType> ret(N_BITS, N_OBJS, genes, al);
    ret.fit = fit;
    
    return ret;
  }

  bool operator==(Organism& obj) {
    for (_uint i = 0; i < al->get_num_params(); ++i) {
      Type t = al->get_type(i);
      if (t == t_real) {
        if (read_real(i) != obj.read_real(i)) {
          return false;
        }
      } else if (t == t_int) {
        if (read_int(i) != obj.read_int(i)) {
          return false;
        }
      } else {
        if (read_uint(i) != obj.read_uint(i)) {
          return false;
        }
      }
    }
    return true;
  }
  bool operator!=(Organism& obj) {
    for (_uint i = 0; i < al->get_num_params(); ++i) {
      Type t = al->get_type(i);
      if (t == t_real) {
        if (read_real(i) != obj.read_real(i)) {
          return true;
        }
      } else if (t == t_int) {
        if (read_int(i) != obj.read_int(i)) {
          return true;
        }
      } else {
        if (read_uint(i) != obj.read_uint(i)) {
          return true;
        }
      }
    }
    return false;
  }
  //bool operator>(Organism& obj) { return !obj.dominates(this); }
  //bool operator<(Organism& obj) { return obj.dominates(this); }

  void swap(Organism& obj) {
    FitType tmp = fit;
    fit = obj.fit;
    obj.fit = tmp;
  }
  bool valid() { return (al != NULL && N_OBJS > 0 && N_BITS > 0); }

  std::vector<Organism*> breed(ArgStore* args, Organism* o) {
    if (get_n_bits() != o->get_n_bits()) {
      error(CODE_MISC, "Cannot breed organsims with a differing number of bits, %d and %d.", get_n_bits(), o->get_n_bits());
    }
    std::vector<Organism*> children(2);

    memset(output_stream, 0, BUF_SIZE);
    Chromosome gene0(genes);
    Chromosome gene1(o->genes);

    if (args->random_crossover()) {
      if (args->get_num_crossovers() <= 0) {
        gene0.exchange_uniform(args, &gene1);
      } else {
        std::uniform_int_distribution<size_t> rint( 0, gene0.get_n_bits() - 1 );
        for (int n = 0; n < args->get_num_crossovers(); ++n) {
          size_t exch_bit = rint( args->get_generator() );
          gene0.exchange(&gene1, exch_bit);
        }
      }
      children[0] = new Organism<FitType>(N_BITS, N_OBJS, gene0, al);
      children[1] = new Organism<FitType>(N_BITS, N_OBJS, gene1, al);
    } else {
      children[0] = new Organism<FitType>(*this);
      children[1] = new Organism<FitType>(*o);
    }
#ifdef MUT_SLOW
    children[0]->genes.slow_mutate(args);
    children[1]->genes.slow_mutate(args);
#else
    children[0]->genes.mutate(args);
    children[1]->genes.mutate(args);
#endif
    return children;
  }
  void mutate(ArgStore* args) {
#ifdef MUT_SLOW
    genes.slow_mutate(args);
#else
    genes.mutate(args);
#endif
  }
  void reset_fitness() {
    fit.reset();
    misc_data = "";
  }
  void randomize(ArgStore* args) { genes.randomize(al.get(), args); }
  void randomize(ArgStore* args, Organism* orgtmp) {
    //make most of the genes similar, with one gene more wildly varied
    std::uniform_int_distribution<size_t> chrom(0, al->get_num_params() - 1);
    size_t high_ind = chrom( args->get_generator() );

    double var = args->get_init_coup_var();
    double lvar = var/al->get_num_params();
    genes.reset();
    for (size_t i = 0; i < al->get_num_params(); i++) {
      Type t = al->get_type(i);
      if (i != high_ind) {
        double mean;
        if (t == t_real) {
          mean = orgtmp->read_real(i);
          double range = al->get_range_max(i) - al->get_range_min(i);
          std::normal_distribution<double> norm(mean, lvar*range);
          double x = norm( args->get_generator() );
          genes.set_to_num(al.get(), i, x);
        } else {
          _uint max_possible = 1 << al->get_block_length(i);
          mean = (double)(max_possible - genes.gene_to_int(al.get(), i))/2;
          //scale lvar to lvar*max_possible/2 and set n and p to produce the according mean and variance
          double p = 1 - lvar*max_possible/(2*mean);
          int n = (int)mean/p;
          std::binomial_distribution<int> dist(n, p);
          int x = dist( args->get_generator() )*2 + genes.gene_to_int(al.get(), i);
          genes.set_to_int(al.get(), i, x);
        }
      }
    }
    if (al->get_type(high_ind) == t_real) {
      //set the genome representation for the high variance index
      double mean = orgtmp->read_real(high_ind);
      double range = al->get_range_max(high_ind) - al->get_range_min(high_ind);
      std::normal_distribution<double> norm(mean, var*range);
      double x = norm( args->get_generator() );
      genes.set_to_num(al.get(), high_ind, x);
    } else {
      int tmpx = orgtmp->read_int(high_ind);
      _uint max_possible = 1 << al->get_block_length(high_ind);
      std::normal_distribution<double> dist((double)tmpx, var*max_possible);
      int x = (int)dist( args->get_generator() );
      genes.set_to_int(al.get(), high_ind, x);
      std::cout << " max_possible = " << max_possible << " x = " << x << " orgtmp_x = " << tmpx << "\n";
    }
  }

  //void evaluate_fitness_noisy(Problem<FitType>* prob, double forget_weight=0);
  void evaluate_fitness(Problem<FitType>* prob) { prob->evaluate_fitness(this); }

  FitType& get_fitness_info() { return fit; }

  double get_fitness(_uint i = 0) { return fit.get_fitness(i); }
  double get_cost(_uint i = 0) { return fit.get_cost(i); }
  void set_fitness(double val) { fit.update(val, 0); }
  void set_cost(double val) { fit.update(-val, 0); }
  void apply_penalty(double val) { penalty = val; }

  double get_penalty() { return penalty; }
  bool penalized() { return penalty != 0; }
  void set_fitness(_uint i, double val) {
    if (i >= fit.get_n_objs()) {
      error(CODE_ARG_RANGE, "Attempt to modify invalid fitness index %d, size is %d.", i, fit.get_n_objs());
    } else {
      if (std::is_base_of<MultiFitness, FitType>::value) {
        fit.update(val, i);
      }
    }
  }
  void set_cost(_uint i, double val) { set_fitness(i, -val); }
  //_uint get_n_evaluations() { return n_evaluations; }
  //average fitness with another organism if they both have the same genotype
  /*void average_fitness(Organism* other);
  void copy_fitness_data(Organism* other);*/

  void set_int(_uint i, int value) { genes.set_to_int(al.get(), i, value); }
  void set_uint(_uint i, _uint value) { genes.set_to_ulong(al.get(), i, value); }
  void set_real(_uint i, double value) { genes.set_to_num(al.get(), i, value); }
  double read_real(_uint i) { return genes.gene_to_num(al.get(), i); }
  int read_int(_uint i) { return genes.gene_to_int(al.get(), i); }
  _uint read_uint(_uint i) { return genes.gene_to_ulong(al.get(), i); }
  //bool dominates(Organism* other);
  String get_chromosome_string(_uint i) {
    if (N_BITS == 0 || !al) {
      error(1, "Attempt to access string for uninitialized organism.");
    }
    if ( i >= al->get_num_params() ) {
      error(1, "Attempt to access invalid parameter with index %d.", i);
    }
    return genes.get_string(al.get(), i);
  }
  char* get_output_stream() { return output_stream; }
  size_t get_output_len() {return output_len; }
  //int get_rank() { return rank; }
  _uint get_n_bits() { return N_BITS; }
  _uint get_n_params() { return al->get_num_params(); }
  _uint get_n_objs() { return N_OBJS; }
};

}

#endif //ORGANISM_H
