#include "convergence.h"

namespace Genetics {

Conv_VarianceCutoff::Conv_VarianceCutoff(double p_cutoff) {
  cutoff = p_cutoff;
  if (cutoff < MIN_CUTOFF) {
    cutoff = MIN_CUTOFF;
  }
}

bool Conv_VarianceCutoff::evaluate_convergence(Vector<FitnessStats> stats) {
  if (stats.size() > max_variance.size()) {
    max_variance.insert(max_variance.end(), stats.size() - max_variance.size(), -std::numeric_limits<double>::infinity());
  }
  bool ret = true;
  for (_uint i = 0; i < stats.size(); ++i) {
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

Conv_RangeCutoff::Conv_RangeCutoff(double p_cutoff) {
  cutoff = p_cutoff;
  if (cutoff < MIN_CUTOFF) {
    cutoff = MIN_CUTOFF;
  }
}

bool Conv_RangeCutoff::evaluate_convergence(Vector<FitnessStats> stats) {
  if (stats.size() > max_range.size()) {
    max_range.insert(max_range.end(), stats.size() - max_range.size(), -std::numeric_limits<double>::infinity());
  }
  bool ret = true;
  for (_uint i = 0; i < stats.size(); ++i) {
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

Conv_Plateau::Conv_Plateau(double p_fitness_threshold, _uint p_generation_cutoff) {
  fitness_threshold = p_fitness_threshold;
  generation_cutoff = p_generation_cutoff;
  gens_without_improvement = 0;
  if (fitness_threshold < MIN_CUTOFF) {
    fitness_threshold = MIN_CUTOFF;
  }
}

bool Conv_Plateau::evaluate_convergence(Vector<FitnessStats> stats) {
  if ( stats.size() > max_fitness.size() || stats.size() > prev_fitness.size() ) {
    max_fitness.insert(max_fitness.end(), stats.size() - max_fitness.size(), -std::numeric_limits<double>::infinity());
    prev_fitness.insert(prev_fitness.end(), stats.size() - prev_fitness.size(), -std::numeric_limits<double>::infinity());
  }
  bool improved = false;
  for (_uint i = 0; i < stats.size(); ++i) {
    double max_rel = stats[i].max;
    if (stats[i].max > max_fitness[i]) {
      max_fitness[i] = stats[i].max;
    }
    if (prev_fitness[i] != 0) {
      //ensure that max_rel is greater than or equal to 1
      if (stats[i].max < 0) {
	max_rel = prev_fitness[i] / stats[i].max;
      } else {
	max_rel /= prev_fitness[i];
      }
    } else {
      max_rel = 1.0;
    }
    if (max_rel - 1.0 >= fitness_threshold) {
      improved = true;
      gens_without_improvement = 0;
      prev_fitness[i] = stats[i].max;
    }
  }
  if (!improved) {
    gens_without_improvement++;
  }
  if (gens_without_improvement > generation_cutoff) {
    return true;
  }
  return false;
}

uConv_VarianceCutoff::uConv_VarianceCutoff(Vector<double> p_cutoffs) :
cutoffs(p_cutoffs)
{
  for (_uint i = 0; i < cutoffs.size(); ++i) {
    if (cutoffs[i] <= 0) {
      cutoffs[i] = std::numeric_limits<double>::infinity();
    }
  }
}

bool uConv_VarianceCutoff::evaluate_convergence(Vector<FitnessStats> stats) {
  if (stats.size() > cutoffs.size()) {
    error(CODE_ARG_RANGE, "population range convergence out of bounds, %d objects were supplied, but bounds were only specified for %d variables.", stats.size(), cutoffs.size());
  }
  for (_uint i = 0; i < stats.size(); ++i) {
    if (stats[i].var > cutoffs[i]) {
      return false;
    }
  }
  return true;
}

uConv_RangeCutoff::uConv_RangeCutoff(Vector<double> p_cutoffs) :
cutoffs(p_cutoffs)
{
  for (_uint i = 0; i < cutoffs.size(); ++i) {
    if (cutoffs[i] <= 0) {
      cutoffs[i] = std::numeric_limits<double>::infinity();
    }
  }
}

bool uConv_RangeCutoff::evaluate_convergence(Vector<FitnessStats> stats) {
  if (stats.size() > cutoffs.size()) {
    error(CODE_ARG_RANGE, "population range convergence out of bounds, %d objects were supplied, but bounds were only specified for %d variables.", stats.size(), cutoffs.size());
  }
  for (_uint i = 0; i < stats.size(); ++i) {
    double range = stats[i].max - stats[i].min;
    if (range > cutoffs[i]) {
      return false;
    }
  }
  return true;
}

uConv_Plateau::uConv_Plateau(Vector<double> p_fitness_thresholds, _uint p_generation_cutoff) :
  fitness_thresholds(p_fitness_thresholds)
{
  generation_cutoff = p_generation_cutoff;
  gens_without_improvement = 0;
  prev_fitness.insert(prev_fitness.end(), p_fitness_thresholds.size(), -std::numeric_limits<double>::infinity() );

  for (_uint i = 0; i < fitness_thresholds.size(); ++i) {
    if (fitness_thresholds[i] <= 0) {
      fitness_thresholds[i] = std::numeric_limits<double>::infinity();
    }
  }
}

bool uConv_Plateau::evaluate_convergence(Vector<FitnessStats> stats) {
  if ( stats.size() > fitness_thresholds.size() || stats.size() > prev_fitness.size() ) {
    fitness_thresholds.insert(fitness_thresholds.end(), stats.size() - fitness_thresholds.size(), -std::numeric_limits<double>::infinity());
    prev_fitness.insert(prev_fitness.end(), stats.size() - prev_fitness.size(), -std::numeric_limits<double>::infinity());
  }

  bool improved = false;
  for (_uint i = 0; i < stats.size(); ++i) {
    if (stats[i].max - prev_fitness[i] > fitness_thresholds[i]) {
      improved = true;
      break;
    }
  }
  if (!improved) {
    gens_without_improvement++;
  } else {
    gens_without_improvement = 0;
    for (_uint i = 0; i < stats.size(); ++i) {
      prev_fitness[i] = stats[i].max;
    }
  }
  return (gens_without_improvement > generation_cutoff);
}

}
