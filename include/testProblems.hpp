#ifndef TEST_PROBLEM_H
#define TEST_PROBLEM_H

#include "../include/convergence.h"
#include "../include/testProblems.hpp"
#include <chrono>
#include <thread>
#include <fstream>

#define DEFAULT_LCG_SEED 0xA1A3A5A7A9ABAFA5

#define NUM_BITS	16
#define NUM_BITS_LARGE  64
#define NUM_OBJS  	2
#define NUM_CHROMS	3

#define NOISY_DOMAIN  1
#define NOISY_VAR     0.25
#define SLEEP_TIME	1000

#define INV_ERR_N 5

typedef Genetics::Organism<Genetics::NoisyFitness> OrganismNoisy;
typedef Genetics::Organism<Genetics::MultiFitness> OrganismMulti;
typedef Genetics::TournamentSelector<Genetics::NoisyFitness> TournamentNoisy;
typedef Genetics::Population<Genetics::NoisyFitness, TournamentNoisy> PopulationNoisy;
typedef Genetics::Organism<Genetics::SingleFitness> OrganismSingle;
typedef Genetics::TournamentSelector<Genetics::SingleFitness> TournamentSingle;
typedef Genetics::Population<Genetics::SingleFitness, TournamentSingle> PopulationSingle;

//for random number generation
class LCG {
private:
  //a-1 is divisible by 4 (and obviously 2 (the only prime factor of 2^64)
  static const unsigned long a = 3*( ((unsigned long)1 << 24) + 5 );
  static const unsigned long c = 170859375;//=15^7 which is relatively prime to m = 2^64
  static const unsigned long x0 = DEFAULT_LCG_SEED;
  static const unsigned long high_mask = ULONG_MAX << 32;
  double c_k[INV_ERR_N];

  // the modulo 64 is implicit
  unsigned long state; void update_state() { state = (state*a + c);}
public:
  LCG(unsigned long seed=x0) : state(seed) {
    //calculate the c_k terms in the taylor series for the inverse error function (from wikipedia)
    c_k[0] = 1;
    for (int k = 1; k < INV_ERR_N; ++k) {
      c_k[k] = 0;
      for (int m = 0; m < k; ++m) {
        c_k[k] += (c_k[m]*c_k[k-1-m])/( (m+1)*(2*m+1) );
      }
    }
    for (int k = 0; k < INV_ERR_N; ++k) {
      //sqrt(pi)
      c_k[k] *= pow(0.88622693, 2*k + 1)/(2*k + 1);
    }
  }
  unsigned long random_ulong() {
    unsigned long ret = state >> 32;
    update_state();
    ret = ret | (state & high_mask);
    return ret;
  }
  double random_int() {
    unsigned long val = random_ulong();
    if ( val & (~(ULONG_MAX >> 1)) ) {
      return -1*((int)val);
    } else {
      return (int)val;
    }
  }
  double random_real() { return (double)random_ulong()/ULONG_MAX; }
  double gaussian_0_1() {
    double uniform = 2*random_real() - 1;
    double ret = 0;
    for (int k = 0; k < INV_ERR_N; ++k) {
      ret += c_k[k]*pow(uniform, 2*k + 1);
    }
    return ret;
  }
};

class TestProblemSingle : public Genetics::Problem<Genetics::SingleFitness> {
public:
  TestProblemSingle();
  void evaluate_fitness(OrganismSingle* org);
};

class TestProblemMulti : public Genetics::Problem<Genetics::MultiFitness> {
public:
  TestProblemMulti();
  void evaluate_fitness(OrganismMulti* org);
};

class TestProblemSlow : public Genetics::Problem<Genetics::SingleFitness> {
public:
  TestProblemSlow();
  void evaluate_fitness(OrganismSingle* org);
};

class TestProblemNoisy : public Genetics::Problem<Genetics::NoisyFitness> {
private:
  LCG gen;
  double domain, penalty_domain, variance;

public:
  TestProblemNoisy(double p_domain = NOISY_DOMAIN, double p_penalty_domain = NOISY_DOMAIN/8, double p_variance = NOISY_VAR);
  double evaluate_fitness_noiseless(OrganismNoisy* org);
  void evaluate_fitness(OrganismNoisy* org);
};

#endif
