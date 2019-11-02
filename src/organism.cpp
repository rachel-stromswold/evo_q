#include "organism.h"

namespace Genetics {

// ==================== COMPARATORS ====================

void MultiFitness::update(double val, _uint i) {
  if (i < N_OBJS) {
    fitness[i] = val;
  }
}
double MultiFitness::get_fitness(_uint i) {
  if (i < N_OBJS) {
    return fitness[i];
  } else if (i == N_OBJS) {
    //TODO: this should probably be replaced with -distance?
    return distance;
  } else if (i == N_OBJS + 1) {
    return -(double)rank;
  } else {
    return -(double)n_dominations;
  }
}

void NoisyFitness::update(double val, _uint i) {
  if (n_evaluations == 0) {
    variance = 0.0;
    fitness = val;
    n_evaluations = 1;
  } else {
    double old_mu = fitness;
    double new_mu = (n_evaluations*fitness + val) / (n_evaluations + 1);
    variance = old_mu*(old_mu - 2*new_mu) + pow(new_mu, 2) + (n_evaluations - 1)*variance/n_evaluations + pow(val - new_mu, 2)/n_evaluations;
    fitness = new_mu;
    ++n_evaluations;
  }
}

void NoisyFitness::average_fitness(NoisyFitness& other) { 
  _uint my_n = n_evaluations;
  _uint their_n = other.n_evaluations; 

  double my_mu_1 = fitness;
  double their_mu_1 = other.get_fitness();
  //these two functions combine the sample means and sample variances from two different measurements
  fitness = (my_n*fitness + their_n*other.get_fitness()) / (my_n + their_n);
  variance = my_n*( my_mu_1*(my_mu_1 - 2*fitness) + pow(fitness, 2) )/(my_n + their_n - 1) + 
             their_n*( their_mu_1*(their_mu_1 - 2*fitness) + pow(fitness, 2) )/(my_n + their_n - 1) +
          ( (my_n - 1)*variance + (their_n - 1)*other.variance )/(my_n + their_n - 1);

  other.variance = variance;
  other.fitness = fitness;

  n_evaluations = my_n + their_n;
  other.n_evaluations = n_evaluations;
}

NoisyFitnessForgetful::NoisyFitnessForgetful(double p_forget_weight) { forget_weight = p_forget_weight; }
void NoisyFitnessForgetful::update(double val, _uint i) { 
  if (!evaluated) {
    fitness = val;
    variance = 0;
    evaluated = true;
  } else {
    //TODO: figure out how to drop previous evaluations (I'm 99.9% sure that you need to keep track of the history)
    double mu = (fitness + forget_weight*val) / (forget_weight + 1);
    //TODO: double check whether this is actually correct
    variance = ( pow(fitness - mu, 2) + pow(forget_weight*val - mu, 2) ) / forget_weight;
    fitness = mu;
  }
}
void NoisyFitnessForgetful::average_fitness(NoisyFitnessForgetful& other) {
  variance = pow((other.fitness - fitness) / 2, 2);
  fitness = (other.fitness + fitness) / 2;
  other.variance = variance;
  other.fitness = fitness;
}

NoisyMultiFitness::NoisyMultiFitness(_uint pn_objs) { N_OBJS = pn_objs; }
double NoisyMultiFitness::get_fitness(_uint i) { return fitness[i]; }
void NoisyMultiFitness::update(double val, _uint i) { 
  if (n_evaluations == 0) {
    variances[i] = 0.0;
    n_evaluations = 1;
  } else {
    double old_mu = val;
    double new_mu = (n_evaluations*fitness[i] + val) / (n_evaluations + 1);
    variances[i] = old_mu*(old_mu - 2*new_mu) + pow(new_mu, 2) + (n_evaluations - 1)*variances[i]/n_evaluations + pow(val - new_mu, 2)/n_evaluations;
    fitness[i] = new_mu;
    //fit_vars[i] = (fit_vars[i]*(n_evaluations - 1) + old_mu*(old_mu - 2*mu) + mu*mu)/n_evaluations;
    ++n_evaluations;
  }
}

}

// ==================== ORGANISM ====================
/*
template <class FitType>
Organism<FitType>::Organism() : genes(0) {
  N_BITS = 0;
  N_OBJS = 0;
  //for (auto it = fitness.begin(); it != fitness.end(); ++it) {
    *it = -std::numeric_limits<double>::infinity();
  }
//  genes = NULL;
}

template <class FitType>
Organism<FitType>::Organism(int pn_bits, int pn_objs, PhenotypeMap* p_al) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs),
  genes(N_BITS),
  al(p_al)
{
//  genes = new Chromosome(N_BITS);
  memset(output_stream, 0, BUF_SIZE);
  reset();
}

template <class FitType>
Organism<FitType>::Organism(int pn_bits, int pn_objs, Chromosome p_genes, PhenotypeMap* p_al) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs),
  genes(p_genes),
  al(p_al)
{
//  genes = new Chromosome(N_BITS);
//  *genes = p_genes;
  reset();
}

template <class FitType>
Organism<FitType>::Organism(int pn_bits, int pn_objs, std::shared_ptr<PhenotypeMap> p_al) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs),
  genes(N_BITS),
  al(p_al)
{
//  genes = new Chromosome(N_BITS);
  memset(output_stream, 0, BUF_SIZE);
  reset();
}

template <class FitType>
Organism<FitType>::Organism(int pn_bits, int pn_objs, Chromosome p_genes, std::shared_ptr<PhenotypeMap> p_al) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs),
  genes(p_genes),
  al(p_al)
{
//  genes = new Chromosome(N_BITS);
//  *genes = p_genes;
  reset();
}


template <class FitType>
void Organism<FitType>::set_fitness_stats(FitType temp) { fit = temp; }

template <class FitType>
void Organism<FitType>::swap(Organism& obj) {
  FitType tmp = fit;
  fit = obj.fit;
  obj.fit = tmp;
}

template <class FitType>
Organism<FitType> Organism<FitType>::copy() {
  Organism<FitType> ret(N_BITS, N_OBJS, genes, al);
  ret.fit = fit;
  
  return ret;
}

template <class FitType>
bool Organism<FitType>::operator==(Organism& obj) {
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

template <class FitType>
bool Organism<FitType>::operator!=(Organism& obj) {
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

template <class FitType>
bool Organism<FitType>::operator>(Organism& obj) {
  return !obj.dominates(this);
}

template <class FitType>
bool Organism<FitType>::operator<(Organism& obj) {
  return obj.dominates(this);
}

template <class FitType>
void Organism<FitType>::reset_fitness() {
  fit.reset();
  misc_data = "";
}

template <class FitType>
std::vector<Organism<FitType>*> Organism<FitType>::breed(ArgStore* args, Organism* o) {
  if (get_n_bits() != o->get_n_bits()) {
    error(CODE_MISC, "Cannot breed organsims with a differing number of bits, %d and %d.", get_n_bits(), o->get_n_bits());
  }
//    N_OBJS = N_OBJS*(N_OBJS + 1)/2;
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

template <class FitType>
void Organism<FitType>::mutate(ArgStore* args) {
#ifdef MUT_SLOW
  genes.slow_mutate(args);
#else
  genes.mutate(args);
#endif
}

template <class FitType>
void Organism<FitType>::randomize(ArgStore* args) {
  genes.randomize(al.get(), args);
}

template <class FitType>
void Organism<FitType>::randomize(ArgStore* args, Organism* orgtmp) {
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

template <class FitType>
void Organism<FitType>::evaluate_fitness(Problem<FitType>* prob) {
  prob->evaluate_fitness(this);
}

template <class FitType>
double Organism<FitType>::get_fitness(unsigned int i) {
  return fit.get_fitness(i);
}

template <class FitType>
double Organism<FitType>::get_cost(unsigned int i) {
  return fit.get_cost(i);
}

template <class FitType>
void Organism<FitType>::set_fitness(double val) {
  fit.set_fitness(val, 0);
}

template <class FitType>
void Organism<FitType>::set_cost(double val) {
  fit.set_fitness(-val, 0);
}

template <class FitType>
void Organism<FitType>::set_fitness(_uint i, double val) {
  if (i >= fit.get_n_objs()) {
    error(CODE_ARG_RANGE, "Attempt to modify invalid fitness index %d, size is %d.", i, fit.get_n_objs());
  } else {
    fit.set_fitness(val, i);
  }
}

template <class FitType>
void Organism<FitType>::set_int(_uint i, int value) {
  genes.set_to_int(al.get(), i, value);
}

template <class FitType>
void Organism<FitType>::set_real(_uint i, double value) {
  genes.set_to_num(al.get(), i, value);
}

template <class FitType>
double Organism<FitType>::read_real(_uint i) {
  return genes.gene_to_num(al.get(), i);
}

template <class FitType>
int Organism<FitType>::read_int(_uint i) {
  return genes.gene_to_int(al.get(), i);
}

template <class FitType>
_uint Organism<FitType>::read_uint(_uint i) {
  return genes.gene_to_ulong(al.get(), i);
}

template <class FitType>
String Organism<FitType>::get_chromosome_string(_uint i) {
  if (N_BITS == 0 || !al) {
    error(1, "Attempt to access string for uninitialized organism.");
  }
  if ( i >= al->get_num_params() ) {
    error(1, "Attempt to access invalid parameter with index %d.", i);
  }
  return genes.get_string(al.get(), i);
}

}*/
