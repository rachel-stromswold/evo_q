#ifndef EVO_Q_HPP
#define EVO_Q_HPP

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include "population.h"

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

class PythonProblem;

class OrganismWrapper {
private:
  std::shared_ptr<ge::Organism> org;
  ge::Problem* prob;

public:
  OrganismWrapper(std::shared_ptr<ge::Organism> p_org, ge::Problem* p_prob);
  OrganismWrapper(ge::Organism* p_org, ge::Problem* p_prob);

  int get_rank() { return org->get_rank(); }
  double get_fitness() { return org->get_fitness(); }
  double read_real(_uint ind) { return org->read_real(ind); }
  int read_int(_uint ind) { return org->read_int(ind); }
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
  unsigned char pop_type;

public:
  PopulationWrapper(_uint N_BITS, _uint N_OBJS, ge::Problem* p_prob, ge::PhenotypeMap* map, ge::String conf_file);
  PopulationWrapper(_uint N_BITS, _uint N_OBJS, ge::Problem* p_prob, ge::Organism* tmplt, ge::PhenotypeMap* map, ge::String conf_file);
  void evaluate() { pop->evaluate(prob); }
  OrganismWrapper* get_parent(int ind);
  py::list get_best(_uint i = 0);
  py::list get_parents();
  py::list get_children();

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
  double get_hypermutation_threshold()	{ return pop->get_args().get_hypermutation_threshold(); }
//  std::shared_ptr<OrganismWrapper>

//  const PopulationWrapper* get_begin() { parent_ind = 0; return this; }
//  std::shared_ptr<OrganismWrapper> get_next() { parent_ind += 1; }
  void iterate() { pop->iterate(); }
};

class PythonProblem : public ge::Problem {
protected:
  bool template_set;
  ge::Organism tmplt_org;
  bool evaluation_set = false;
  py::function evaluation_func;

public:
  PythonProblem(std::string module_name, int n_bits, int n_params, int n_objs) : 
    ge::Problem(n_bits, n_params, n_objs),
    tmplt_org(n_bits, n_objs, &map)
  {
    template_set = false;
  }

  void evaluate_fitness(ge::Organism* org);
  void set_phenotype_parameters(py::list param_list);

  void set_parameter_range(unsigned param_ind, double min, double max) { map.set_range(param_ind, min, max); }
  void set_fitness_function(py::function f) { evaluation_func = f; evaluation_set = true; }

  void set_template_parameter(unsigned param_ind, py::object o);

  PopulationWrapper* initialize_population(ge::String conf_file);
};

#endif //EVO_Q_HPP