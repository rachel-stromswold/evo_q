#ifndef CONVERGENCE_H
#define CONVERGENCE_H

namespace Genetics {

class Conv_VarianceCutoff : public ConvergenceCriteria {
  private:
    _uint N_OBJS;
    double cutoff;

    Vector<double> max_variance;//used for normalization
  public:
    Conv_VarianceCutoff(_uint N_OBJS, double p_cutoff);
    bool evaluate_convergence(FitnessStats* stats);
};

class Conv_RangeCutoff : public ConvergenceCriteria {
  private:
    _uint N_OBJS;
    double cutoff;

    Vector<double> max_range;//used for normalization
  public:
    Conv_RangeCutoff(_uint N_OBJS, double p_cutoff);
    bool evaluate_convergence(FitnessStats* stats);
};

class Conv_Plateau : public ConvergenceCriteria {
  private:
    _uint N_OBJS;
    double fitness_threshold;
    _uint generation_cutoff;
    _uint gens_without_improvement;

    Vector<double> max_fitness;//used for normalization
    Vector<double> prev_fitness;
  public:
    Conv_Plateau(_uint N_OBJS, double p_fitness_threshold, _uint p_generation_cutoff);
    bool evaluate_convergence(FitnessStats* stats);
};

class uConv_VarianceCutoff : public ConvergenceCriteria {
  private:
    _uint N_OBJS;
    Vector<double> cutoffs;
  public:
    uConv_VarianceCutoff(_uint N_OBJS, Vector<double> p_cutoffs);
    bool evaluate_convergence(FitnessStats* stats);
};

class uConv_RangeCutoff : public ConvergenceCriteria {
  private:
    _uint N_OBJS;
    Vector<double> cutoffs;
  public:
    uConv_RangeCutoff(_uint N_OBJS, Vector<double> p_cutoffs);
    bool evaluate_convergence(FitnessStats* stats);
};

class uConv_Plateau : public ConvergenceCriteria {
  private:
    _uint N_OBJS;
    _uint generation_cutoff;
    _uint gens_without_improvement;

    Vector<double> fitness_thresholds;
    Vector<double> prev_fitness;
  public:
    uConv_Plateau(_uint N_OBJS, Vector<double> p_fitness_thresholds, _uint p_generation_cutoff);
    bool evaluate_convergence(FitnessStats* stats);
};

}

#endif//CONVERGENCE_H
