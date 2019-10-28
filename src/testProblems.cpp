#include "testProblems.hpp"

// =============================== SINGLE ===============================

TestProblemSingle::TestProblemSingle() : Genetics::Problem(NUM_BITS, NUM_CHROMS, 1) {
  map->initialize(1, Genetics::t_real);
}
void TestProblemSingle::evaluate_fitness(Genetics::Organism* org) {
  double x = org->read_real(0);
  org->set_fitness(0, -(x*x));
}

// =============================== MULTI ===============================

TestProblemMulti::TestProblemMulti() : Genetics::Problem(NUM_BITS, NUM_CHROMS, NUM_OBJS) {
  map->initialize(NUM_CHROMS, Genetics::t_real);
}
void TestProblemMulti::evaluate_fitness(Genetics::Organism* org) {
  double x = org->read_real(0);
  org->set_fitness(0, -x*x);
  org->set_fitness(1, -(x - 2)*(x - 2));
  org->apply_penalty(0.0);
}

// =============================== SLOW ===============================

TestProblemSlow::TestProblemSlow() : Genetics::Problem(NUM_BITS, NUM_CHROMS, 1) {
  map->initialize(1, Genetics::t_real);
}
void TestProblemSlow::evaluate_fitness(Genetics::Organism* org) {
  std::this_thread::sleep_for( std::chrono::milliseconds(SLEEP_TIME) );
  double x = org->read_real(0);
  org->set_fitness(0, -(x*x));
}

// =============================== NOISY ===============================

TestProblemNoisy::TestProblemNoisy(double p_domain, double p_penalty_domain, double p_variance) :
Genetics::Problem(NUM_BITS, NUM_CHROMS, 1),
domain(p_domain),
penalty_domain(p_penalty_domain),
variance(p_variance) {
  Genetics::Vector<Genetics::VarContainer> tmp_vars;
  for (int i = 0; i < NUM_CHROMS; ++i) {
    tmp_vars.emplace_back(0, -domain, domain, Genetics::t_real);
  }
  map->initialize(tmp_vars);
  penalty_domain = abs(penalty_domain);
}

double TestProblemNoisy::evaluate_fitness_noiseless(Genetics::Organism* org) {
  double ret = 0.0;
  double penalty = 0;
  for (unsigned i = 0; i < NUM_CHROMS; ++i) {
    double x_i = org->read_real(i);
    ret += x_i*x_i;
    if (x_i < penalty_domain && x_i > -penalty_domain) {
      penalty += 1;
    }
  }
  org->apply_penalty(penalty);
  return ret;
}

void TestProblemNoisy::evaluate_fitness(Genetics::Organism* org) {
  double val = evaluate_fitness_noiseless(org);
  val += variance*gen.gaussian_0_1();
  org->set_cost(0, val);
}
