#ifndef EVO_Q_HPP
#define EVO_Q_HPP

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include "population.h"
#include "convergence.h"

#define POP_SINGLE	1
#define POP_NSGAII	2

typedef unsigned int _uint;
namespace py = pybind11;
namespace ge = Genetics;

struct ProblemData {
  ge::Chromosome gene;
  ge::Vector<double> fit_vals;
  bool evaluation_complete;
  bool evaluation_ready;
  ProblemData(int pn_bits, int pn_objs) : gene(pn_bits), fit_vals(pn_objs, 0.0) {}
};

class PythonException : public std::exception {
public:
    explicit PythonException(std::string m) : message(m) {}
    virtual const char * what() const noexcept override {return message.c_str();}
private:
    std::string message = "";
};

class ConvergenceWrapper : public ge::ConvergenceCriteria {
  //using ConvergenceCriteria::ConvergenceCriteria;
  bool evaluate_convergence(_uint N_OBJS, ge::FitnessStats* stats) {
    PYBIND11_OVERLOAD_PURE(
	bool,			/** return type **/
	ge::ConvergenceCriteria,/** parent class**/
	evaluate_convergence,	/** function to override**/
	N_OBJS,			/** argument(s) **/
	stats
    );
  }
};

class PythonFitness;

class PythonFitnessSettings {
protected:
  bool function_set = false;
  py::function update_func;

  bool noise_compensate;
  double forget_weight;

public:
  void set_n_objs() { n_objs = pn_objs; }
  void set_update_function(py::function f) { update_func = f;function_set = true; }
  friend class PythonFitness
};

class PythonFitness : public ge::NoisyMultiFitness {
private:
  PythonFitnessSettings& fit_sets;

public:
  void update(double val, _uint i = 0);
};

class OrganismWrapper {
private:
  std::shared_ptr<ge::Organism> org;
  ge::Problem* prob;

public:
  OrganismWrapper(std::shared_ptr<ge::Organism> p_org, ge::Problem* p_prob);
  OrganismWrapper(ge::Organism* p_org, ge::Problem* p_prob);

  bool check_validity();
  int get_rank();
  double get_fitness();
  py::list get_fitness_list();
  double read_real(_uint ind);
  int read_int(_uint ind);
  int read_uint(_uint ind);
  py::list get_phenotype();
  void set_fitness(double val);
  void set_fitness(py::list fit_vals);
};

class PopulationWrapper {
private:
  ge::Population* pop;
  ge::Problem* prob;
  size_t parent_ind;
  size_t child_ind;
//  unsigned char pop_type;

public:
  PopulationWrapper(_uint N_BITS, _uint N_OBJS, ge::Problem* p_prob, std::shared_ptr<ge::PhenotypeMap> map, ge::String conf_file, bool latin);
  PopulationWrapper(_uint N_BITS, _uint N_OBJS, ge::Problem* p_prob, ge::Organism* tmplt, std::shared_ptr<ge::PhenotypeMap> map, ge::String conf_file);
  void evaluate();
  void iterate();
  py::dict run(ge::ConvergenceCriteria* convp, bool store_intermediate);
  OrganismWrapper* get_parent(int ind);
  py::list get_best(_uint i = 0);
  py::list get_parents();
  py::list get_children();
  double get_max_fitness(_uint i = 0) { return pop->get_pop_stats(i).max; }
  double get_min_fitness(_uint i = 0) { return pop->get_pop_stats(i).min; }
  double get_fitness_range(_uint i = 0) { return get_max_fitness() - get_min_fitness(); }
  double get_fitness_var(_uint i = 0)  { return pop->get_pop_stats(i).var; }

  void set_cost(_uint i = 0) { pop->set_cost(i); }
  size_t get_pop_size()			{ return pop->get_args().get_pop_size(); }
  void set_pop_size(size_t n)		{ pop->resize_population(n); }
  size_t get_survivors()		{ return pop->get_args().get_survivors(); }
  void set_survivors(size_t n)		{ pop->set_n_survivors(n); }
  size_t get_num_gens() 		{ return pop->get_args().get_num_gens(); }
  void set_num_gens(size_t n) 		{ pop->get_args().set_num_gens(n); }
  int get_num_crossovers() 		{ return pop->get_args().get_num_crossovers(); }
  void set_num_crossovers(int n) 	{ pop->get_args().set_num_crossovers(n); }
  double get_init_coup_var()		{ return pop->get_args().get_init_coup_var(); }
  void set_init_coup_var(double x)	{ pop->get_args().set_init_coup_var(x); }
  double get_init_coup_mean() 		{ return pop->get_args().get_init_coup_mean(); }
  void set_init_coup_mean(double x) 	{ pop->get_args().set_init_coup_mean(x); }
  double get_mutate_prob()		{ return pop->get_args().get_mutate_prob(); }
  void set_mutate_prob(double x)	{ pop->get_args().set_mutate_prob(x); }
  double get_crossover_prob()		{ return pop->get_args().get_crossover_prob(); }
  void set_crossover_prob(double x)	{ pop->get_args().set_crossover_prob(x); }
  double get_hypermutation_threshold()	{ return pop->get_args().get_hypermutation_threshold(); }
  void set_hypermutation_threshold(double x)	{ pop->get_args().set_hypermutation_threshold(x); }
  void set_noise_compensation(_uint val) { pop->get_args().set_noise_compensation(val); }
  void set_selection(std::string str);
//  std::shared_ptr<OrganismWrapper>

//  const PopulationWrapper* get_begin() { parent_ind = 0; return this; }
//  std::shared_ptr<OrganismWrapper> get_next() { parent_ind += 1; }
};

class PythonProblem : public ge::Problem {
protected:
  bool template_set;
  ge::Organism tmplt_org;
  bool evaluation_set = false;
  py::function evaluation_func;

public:
  PythonProblem(int n_bits, int n_params, int n_objs) : 
    ge::Problem(n_bits, n_params, n_objs),
    tmplt_org(n_bits, n_objs, map)
  {
    template_set = false;
  }

  void evaluate_fitness(ge::Organism* org);
  void set_phenotype_parameters(py::list param_list);

  void set_parameter_range(unsigned param_ind, double min, double max) { map->set_range(param_ind, min, max); }
  void set_fitness_function(py::function f) { evaluation_func = f; evaluation_set = true; }

  void set_template_parameter(unsigned param_ind, py::object o);

  PopulationWrapper* initialize_population(ge::String conf_file, bool latin);
};

#endif //EVO_Q_HPP
