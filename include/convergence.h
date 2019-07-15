#ifndef CONVERGENCE_H
#define CONVERGENCE_H

#include "population.h"

#ifndef MIN_CUTOFF
#define MIN_CUTOFF	0.01
#endif

namespace Genetics {

class Conv_VarianceCutoff : public ConvergenceCriteria {
  private:
    double cutoff;

    Vector<double> max_variance;//used for normalization
  public:
    Conv_VarianceCutoff(double p_cutoff);
    bool evaluate_convergence(_uint N_OBJS, FitnessStats* stats);
};

class Conv_RangeCutoff : public ConvergenceCriteria {
  private:
    double cutoff;

    Vector<double> max_range;//used for normalization
  public:
    Conv_RangeCutoff(double p_cutoff);
    bool evaluate_convergence(_uint N_OBJS, FitnessStats* stats);
};

class Conv_Plateau : public ConvergenceCriteria {
  private:
    double fitness_threshold;
    _uint generation_cutoff;
    _uint gens_without_improvement;

    Vector<double> max_fitness;//used for normalization
    Vector<double> prev_fitness;
  public:
    Conv_Plateau(double p_fitness_threshold, _uint p_generation_cutoff);
    bool evaluate_convergence(_uint N_OBJS, FitnessStats* stats);
};

class uConv_VarianceCutoff : public ConvergenceCriteria {
  private:
    Vector<double> cutoffs;
  public:
    uConv_VarianceCutoff(Vector<double> p_cutoffs);
    bool evaluate_convergence(_uint N_OBJS, FitnessStats* stats);
};

class uConv_RangeCutoff : public ConvergenceCriteria {
  private:
    Vector<double> cutoffs;
  public:
    uConv_RangeCutoff(Vector<double> p_cutoffs);
    bool evaluate_convergence(_uint N_OBJS, FitnessStats* stats);
};

class uConv_Plateau : public ConvergenceCriteria {
  private:
    _uint generation_cutoff;
    _uint gens_without_improvement;

    Vector<double> fitness_thresholds;
    Vector<double> prev_fitness;
  public:
    uConv_Plateau(Vector<double> p_fitness_thresholds, _uint p_generation_cutoff);
    bool evaluate_convergence(_uint N_OBJS, FitnessStats* stats);
};

}

#endif//CONVERGENCE_H
