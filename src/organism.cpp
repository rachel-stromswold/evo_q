#include "organism.h"

namespace Genetics {

Organism::Organism() : genes(0), fitness(1) {
  N_BITS = 0;
  N_OBJS = 0;
  for (auto it = fitness.begin(); it != fitness.end(); ++it) {
    *it = -std::numeric_limits<double>::infinity();
  }
//  genes = NULL;
}

Organism::Organism(int pn_bits, int pn_objs, PhenotypeMap* p_al) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs),
  fitness(N_OBJS, (double)0.0),
  genes(N_BITS),
  al(p_al)
{
//  genes = new Chromosome(N_BITS);
  memset(output_stream, 0, BUF_SIZE);
  reset();
}

Organism::Organism(int pn_bits, int pn_objs, Chromosome p_genes, PhenotypeMap* p_al) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs),
  fitness(N_OBJS, (double)0.0),
  genes(p_genes),
  al(p_al)
{
//  genes = new Chromosome(N_BITS);
//  *genes = p_genes;
  reset();
}

Organism::Organism(int pn_bits, int pn_objs, std::shared_ptr<PhenotypeMap> p_al) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs),
  fitness(N_OBJS, (double)0.0),
  genes(N_BITS),
  al(p_al)
{
//  genes = new Chromosome(N_BITS);
  memset(output_stream, 0, BUF_SIZE);
  reset();
}

Organism::Organism(int pn_bits, int pn_objs, Chromosome p_genes, std::shared_ptr<PhenotypeMap> p_al) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs),
  fitness(N_OBJS, (double)0.0),
  genes(p_genes),
  al(p_al)
{
//  genes = new Chromosome(N_BITS);
//  *genes = p_genes;
  reset();
}

/*Organism::~Organism() {
  if (genes) {
    delete genes;
  }
}*/

/*Organism::Organism(const Organism &obj) :
  fitness(obj.fitness),
  misc_data(obj.misc_data)
{
  N_BITS = obj.N_BITS;
  N_OBJS = obj.N_OBJS;
  genes = new Chromosome(*(obj.genes));
  al = obj.al;
}

Organism Organism::copy() {
  Organism ret(N_BITS, N_OBJS, al);
  ret.genes = new Chromosome(*genes);
  return ret;
}

Organism::Organism(Organism&& obj) :
  fitness(std::move(obj.fitness)),
  misc_data(std::move(obj.misc_data))
{
  genes = obj.genes;
  al = obj.al;
  obj.genes = NULL;
}

Organism& Organism::operator=(Organism& obj) {
  Chromosome* tmp_genes = genes;
  genes = obj.genes;
  fitness = obj.fitness;
  misc_data = obj.misc_data;
  obj.genes = tmp_genes;
  al = obj.al;
  return *this;
}*/

void Organism::swap(Organism& obj) {
  _uint tmp_n_evaluations = n_evaluations;
  n_evaluations = obj.n_evaluations;
  obj.n_evaluations = n_evaluations;
  fitness.swap(obj.fitness);
}

/*Organism& Organism::operator=(Organism& obj) {
  swap(obj);
  return *this;
}*/

Organism Organism::copy() {
  Organism ret(N_BITS, N_OBJS, genes, al);
  for (_uint i = 0; i < N_OBJS; ++i) {
    ret.set_fitness(i, get_fitness(i));
  }
  ret.n_evaluations = n_evaluations;
  return ret;
}

bool Organism::operator==(Organism& obj) {
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

bool Organism::operator!=(Organism& obj) {
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

bool Organism::operator>(Organism& obj) {
  return !obj.dominates(this);
}

bool Organism::operator<(Organism& obj) {
  return obj.dominates(this);
}

void Organism::reset() {
  for (_uint i = 0; i < fitness.size(); ++i) {
    fitness[i] = 0;
  }
  misc_data = "";
}

std::vector<Organism*> Organism::breed(ArgStore* args, Organism* o) {
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
    children[0] = new Organism(N_BITS, N_OBJS, gene0, al);
    children[1] = new Organism(N_BITS, N_OBJS, gene1, al);
  } else {
    children[0] = new Organism(*this);
    children[1] = new Organism(*o);
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

void Organism::mutate(ArgStore* args) {
#ifdef MUT_SLOW
  genes.slow_mutate(args);
#else
  genes.mutate(args);
#endif
}

void Organism::randomize(ArgStore* args) {
  genes.randomize(al.get(), args);
}

void Organism::randomize(ArgStore* args, Organism* orgtmp) {
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

void Organism::evaluate_fitness_noisy(Problem* prob) {
  double* prev_fitness = (double*)malloc(sizeof(double)*N_OBJS);
  for (_uint i = 0; i < N_OBJS; ++i) {
    prev_fitness[i] = fitness[i];
  }
  prob->evaluate_fitness(this);
  for (_uint i = 0; i < N_OBJS; ++i) {
    //the commented line does not forget old evaluations, the new one does
    //fitness[i] = (n_evaluations*prev_fitness[i] + fitness[i]) / (n_evaluations + 1);
    fitness[i] = (prev_fitness[i] + fitness[i]) / 2;
  }
  ++n_evaluations;
}

void Organism::evaluate_fitness(Problem* prob) {
  prob->evaluate_fitness(this);
  ++n_evaluations;
}

double Organism::get_fitness(unsigned int i) {
  if (i == N_OBJS)
    return distance;
  else if (i == N_OBJS + 1)
    return -(double)rank;
  else if (i == N_OBJS + 2)
    return -(double)n_dominations;
  else
    return fitness[i];
}

double Organism::get_cost(unsigned int i) {
  if (i == N_OBJS)
    return -distance;
  else if (i == N_OBJS + 1)
    return (double)rank;
  else if (i == N_OBJS + 2)
    return (double)n_dominations;
  else
    return -fitness[i];
}

void Organism::set_fitness(double val) {
  fitness[0] = val;
}

void Organism::set_cost(double val) {
  fitness[0] = -val;
}

void Organism::set_fitness(_uint i, double val) {
  if (i >= fitness.size()) {
    if (i < N_OBJS) {
      fitness.resize(N_OBJS);
    } else {
      error(CODE_ARG_RANGE, "Attempt to modify invalid fitness index %d, size is %d.", i, fitness.size());
    }
  }
  fitness[i] = val;
}

void Organism::set_cost(_uint i, double val) {
  if (i >= fitness.size()) {
    if (i < N_OBJS) {
      fitness.resize(N_OBJS);
    } else {
      error(CODE_ARG_RANGE, "Attempt to modify invalid fitness index %d, size is %d.", i, fitness.size());
    }
  }
  fitness[i] = -val;
}

void Organism::set_int(_uint i, int value) {
  genes.set_to_int(al.get(), i, value);
}

void Organism::set_real(_uint i, double value) {
  genes.set_to_num(al.get(), i, value);
}

double Organism::read_real(_uint i) {
  return genes.gene_to_num(al.get(), i);
}

int Organism::read_int(_uint i) {
  return genes.gene_to_int(al.get(), i);
}

_uint Organism::read_uint(_uint i) {
  return genes.gene_to_ulong(al.get(), i);
}

bool Organism::dominates(Organism* other) {
  bool ret = false;
  for (_uint i = 0; i < fitness.size(); ++i) {
    if ( fitness[i] > other->get_fitness(i) ) {
      ret = true;
    }
    if ( other->get_fitness(i) > fitness[i]) {
      return false;
    }
  }
  return ret;
}

}
