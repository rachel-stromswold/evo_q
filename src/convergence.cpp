#include "convergence.h"

namespace Genetics {

Conv_VarianceCutoff::Conv_VarianceCutoff(_uint N_OBJS, double p_cutoff) :
  max_variance(N_OBJS, -std::numeric_limits<double>::infinity())
{
  this->N_OBJS = N_OBJS;

  cutoff = p_cutoff;
  if (cutoff < MIN_CUTOFF) {
    cutoff = MIN_CUTOFF;
  }
}

bool Conv_VarianceCutoff::evaluate_convergence(FitnessStats* stats) {
  bool ret = true;
  for (_uint i = 0; i < N_OBJS; ++i) {
    if (stats[i].var > max_variance[i]) {
      max_variance[i] = stats[i].var;
    }
    if (max_variance[i] != 0) {
      stats[i].var /= max_variance[i];
    }
    if (stats[i].var > cutoff) {
      ret = false;
    }
  }
  return ret;
}

Conv_RangeCutoff::Conv_RangeCutoff(_uint N_OBJS, double p_cutoff) :
  max_range(N_OBJS, -std::numeric_limits<double>::infinity())
{
  this->N_OBJS = N_OBJS;

  cutoff = p_cutoff;
  if (cutoff < MIN_CUTOFF) {
    cutoff = MIN_CUTOFF;
  }
}

bool Conv_RangeCutoff::evaluate_convergence(FitnessStats* stats) {
  bool ret = true;
  for (_uint i = 0; i < N_OBJS; ++i) {
    double range = stats[i].max - stats[i].min;
    if (range > max_range[i]) {
      max_range[i] = range;
    }
    if (max_range[i] != 0) {
      range /= max_range[i];
    }
    if (range > cutoff) {
      ret = false;
    }
  }
  return ret;
}

Conv_Plateau::Conv_Plateau(_uint N_OBJS, double p_fitness_threshold, _uint p_generation_cutoff) :
  max_fitness(N_OBJS, -std::numeric_limits<double>::infinity()),
  prev_fitness(N_OBJS, -std::numeric_limits<double>::infinity())
{
  this->N_OBJS = N_OBJS;

  fitness_threshold = p_fitness_threshold;
  generation_cutoff = p_generation_cutoff;
  gens_without_improvement = 0;
  if (fitness_threshold < MIN_CUTOFF) {
    fitness_threshold = MIN_CUTOFF;
  }
  for (_uint i = 0; i < N_OBJS; ++i) {
    prev_fitness[i] = -std::numeric_limits<double>::infinity();
    max_fitness[i] = -std::numeric_limits<double>::infinity();
  }
}

bool Conv_Plateau::evaluate_convergence(FitnessStats* stats) {
  for (_uint i = 0; i < N_OBJS; ++i) {
    double max_rel = stats[i].max;
    if (stats[i].max > max_fitness[i]) {
      max_fitness[i] = stats[i].max;
    }
    if (prev_fitness[i] != 0) {
      max_rel /= prev_fitness[i];
    } else {
      max_rel += 1.0;
    }
    if (max_rel - 1.0 < fitness_threshold &&
	1.0 - max_rel < fitness_threshold) {
      gens_without_improvement++;
    } else {
      gens_without_improvement = 0;
      prev_fitness[i] = stats[i].max;
    }
  }
  if (gens_without_improvement > generation_cutoff) {
    return true;
  }
  return false;
}

uConv_VarianceCutoff::uConv_VarianceCutoff(_uint N_OBJS, Vector<double> p_cutoffs) :
  cutoffs(p_cutoffs)
{
  this->N_OBJS = N_OBJS;

  for (_uint i = 0; i < N_OBJS; ++i) {
    if (cutoffs[i] <= 0) {
      cutoffs[i] = std::numeric_limits<double>::infinity();
    }
  }
}

bool uConv_VarianceCutoff::evaluate_convergence(FitnessStats* stats) {
  for (_uint i = 0; i < N_OBJS; ++i) {
    if (stats[i].var > cutoffs[i]) {
      return false;
    }
  }
  return true;
}

uConv_RangeCutoff::uConv_RangeCutoff(_uint N_OBJS, Vector<double> p_cutoffs) :
  cutoffs(p_cutoffs)
{
  this->N_OBJS = N_OBJS;

  for (_uint i = 0; i < N_OBJS; ++i) {
    if (cutoffs[i] <= 0) {
      cutoffs[i] = std::numeric_limits<double>::infinity();
    }
  }
}

bool uConv_RangeCutoff::evaluate_convergence(FitnessStats* stats) {
  for (_uint i = 0; i < N_OBJS; ++i) {
    double range = stats[i].max - stats[i].min;
    if (range > cutoffs[i]) {
      return false;
    }
  }
  return true;
}

uConv_Plateau::uConv_Plateau(_uint N_OBJS, Vector<double> p_fitness_thresholds, _uint p_generation_cutoff) :
  fitness_thresholds(p_fitness_thresholds),
  prev_fitness(N_OBJS, -std::numeric_limits<double>::infinity())
{
  this->N_OBJS = N_OBJS;

  generation_cutoff = p_generation_cutoff;
  gens_without_improvement = 0;
  for (_uint i = 0; i < N_OBJS; ++i) {
    if (fitness_thresholds[i] <= 0) {
      fitness_thresholds[i] = std::numeric_limits<double>::infinity();
    }
  }
}

bool uConv_Plateau::evaluate_convergence(FitnessStats* stats) {
  for (_uint i = 0; i < N_OBJS; ++i) {
    if (stats[i].max - prev_fitness[i] < fitness_thresholds[i] &&
	prev_fitness[i] - stats[i].max < fitness_thresholds[i]) {
      gens_without_improvement++;
    } else {
      gens_without_improvement = 0;
      prev_fitness[i] = stats[i].max;
    }
  }
  return (gens_without_improvement > generation_cutoff);
}

}
