#include "organism.h"

namespace Genetics {

// ==================== COMPARATORS ====================

void MultiFitness::update(double val, _uint i) {
  if (i < N_OBJS) {
    fitness[i] = val;
  }
}
void MultiFitness::set_fitness(double val, _uint i) {
  update(val, i);
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
//double NoisyMultiFitness::get_fitness(_uint i) { return fitness[i]; }
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

double NoisyMultiFitness::get_fitness(_uint i) {
  return fitness[i];
}

}
