#include "python.hpp"

int str_starts(const char* test_str, const char* match_str) { 
  for (size_t i = 0; match_str[i] != 0; ++i) {
    if (match_str[i] == 0 || test_str[i] < match_str[i]) {
      return -1;
    }
    if (match_str[i] > test_str[i]) {
      return 1;
    }
  }
  return 0;
}

// ============================= ORGANISM_WRAPPER =============================

OrganismWrapper::OrganismWrapper(std::shared_ptr<ge::Organism> p_org, ge::Problem* p_prob) :
org(p_org) {
  prob = p_prob;
}

OrganismWrapper::OrganismWrapper(ge::Organism* p_org, ge::Problem* p_prob) {
  org = std::make_shared<ge::Organism>(*p_org);
  prob = p_prob;
}

bool OrganismWrapper::check_validity() {
  if(!org) {
    throw std::runtime_error("Attempt to access uninitialized organism");
    return false;
  }
  return true;
}

int OrganismWrapper::get_rank() {
  check_validity();
  return org->get_rank();
}
double OrganismWrapper::get_fitness() {
  check_validity();
  return org->get_fitness();
}
py::list OrganismWrapper::get_fitness_list() {
  py::list ret;
  for (_uint i = 0; i < org->get_n_objs(); ++i) {
    ret.append(org->get_fitness(i));
  }
  return ret;
}
double OrganismWrapper::read_real(_uint ind) {
  check_validity();
  return org->read_real(ind);
}
int OrganismWrapper::read_int(_uint ind) {
  check_validity();
  return org->read_int(ind);
}
int OrganismWrapper::read_uint(_uint ind) {
  check_validity();
  return org->read_uint(ind);
}

void OrganismWrapper::set_fitness(double val) {
  check_validity();
  org->set_fitness(val);
}

void OrganismWrapper::set_fitness(py::list fit_vals) {
  check_validity();
  for (size_t i = 0; i < fit_vals.size() && i < org->get_n_objs(); ++i) {
    org->set_fitness(i, fit_vals[i].cast<double>());
  }
}

py::list OrganismWrapper::get_phenotype() {
  check_validity();
  py::list ret;
  for (size_t i = 0; i < org->get_n_params(); ++i) {
    Genetics::String s = org->get_chromosome_string(i);
    ret.append( std::string(s.c_str()) );
  }
  return ret;
}

// ============================= POPULATION_WRAPPER =============================

PopulationWrapper::PopulationWrapper(_uint N_BITS, _uint N_OBJS, Genetics::Problem* p_prob, std::shared_ptr<Genetics::PhenotypeMap> map, ge::String conf_file) {
  if (N_OBJS <= 1) {
    pop = new ge::Population(N_BITS, N_OBJS, map, conf_file);
    pop_type = POP_SINGLE;
  } else {
    pop = new ge::Population_NSGAII(N_BITS, N_OBJS, map, conf_file);
    pop_type = POP_NSGAII;
  }
  prob = p_prob;
}

PopulationWrapper::PopulationWrapper(_uint N_BITS, _uint N_OBJS, ge::Problem* p_prob, ge::Organism* tmplt, std::shared_ptr<ge::PhenotypeMap> map, ge::String conf_file)
{
  if (N_OBJS <= 1) {
    pop = new ge::Population(N_BITS, N_OBJS, map, conf_file);
    pop_type = POP_SINGLE;
  } else {
    pop = new ge::Population_NSGAII(N_BITS, N_OBJS, map, conf_file);
    pop_type = POP_NSGAII;
  }
  prob = p_prob;
}

OrganismWrapper* PopulationWrapper::get_parent(int ind) {
  return new OrganismWrapper(pop->get_organism(ind), prob);
}

py::list PopulationWrapper::get_best(_uint i) {
  py::list ret;
  if (pop_type == POP_SINGLE) {
    std::shared_ptr<Genetics::Organism> best = pop->get_best_organism();
    if (!best) {
      throw std::runtime_error("Couldn't find best organism, perhaps no evaluations have been performed?");
    }
    ret.append( OrganismWrapper(best, prob) );
  } else {
    std::vector<std::shared_ptr<ge::Organism>> par = ((ge::Population_NSGAII*)pop)->get_pareto_front(i);
    for (_uint i = 0; i < par.size(); ++i) {
      ret.append( OrganismWrapper(par[i], prob) );
    }
  }

  return ret;
}

py::list PopulationWrapper::get_parents() {
  py::list ret;
  for (size_t i = 0; i < pop->get_offspring_num(); ++i) {
    ret.append( OrganismWrapper(pop->get_organism(i), prob) );
  }
//  for (int i = 0; i < pop->get_children_size(); ++i) {
//    std::cout << "c++: " << get_chromosome_string(0) << "\n";
//  }
  return ret;
}

py::list PopulationWrapper::get_children() {
  py::list ret;
  for (size_t i = 0; i < pop->get_offspring_num(); ++i) {
    ret.append( OrganismWrapper(pop->get_child(i), prob) );
  }
  return ret;
}

void PopulationWrapper::evaluate() {
  if (pop->get_n_objs() > 1) {
    ((Genetics::Population_NSGAII*)pop)->evaluate(prob);
  } else {
    pop->evaluate(prob);
  }
}

void PopulationWrapper::iterate() {
  if (pop->get_n_objs() > 1) {
    ((Genetics::Population_NSGAII*)pop)->iterate();
  } else {
    pop->iterate();
  }
}

// ============================= PYTHON_PROBLEM =============================

void PythonProblem::evaluate_fitness(Genetics::Organism* org) {
  if (!evaluation_set) {
    throw std::runtime_error("fitness function must be set before evaluations can be performed");
  }
  py::list fit_vals_list = evaluation_func( OrganismWrapper(org, (Genetics::Problem*)this) );
  if (fit_vals_list.size() < org->get_n_objs()) {
    std::stringstream iss;
    iss << "list returned from evaluation function is shorter (length " << fit_vals_list.size() << ") than the number of objectives (" << org->get_n_objs() << ")";
    throw PythonException( iss.str() );
  }
  for (size_t i = 0; i < org->get_n_objs(); ++i) {
    org->set_fitness(i, fit_vals_list[i].cast<double>());
  }
}

void PythonProblem::set_phenotype_parameters(py::list param_list) {
  if (param_list.size() != N_PARAMS) {
    Genetics::error(1, "Expected %d parameters in the list, %d were provided.", N_PARAMS, param_list.size());
  }
  Genetics::Vector<Genetics::VarContainer> vc_list;
  for (size_t i = 0; i < param_list.size(); ++i) {
    std::string type = param_list[i].cast<std::string>();
    const char* type_str = type.data();
#ifdef PRINT_DEBUG_DATA
    std::cout << "orig = " << type << ", c_str = " << type_str << "\n";
#endif
    if (str_starts(type_str, "real") == 0) {
      const char* range_dat = strchr(type_str, (int)'(');
      vc_list.push_back( Genetics::VarContainer(0, 0, 0, Genetics::t_real) );
      //parse the range data and set values accordingly
      if (range_dat != NULL) {
	//copy over the portion of the string with (<min_value>, <max_value>)
	char* min_str = strdup(range_dat + 1);
	char* max_str = (char*)strchr(type_str, (int)',');//cast to silence gcc complaint
	//Don't you just love hackish parsing and manipulation of c-strings?
	if (max_str) {
	  *max_str = 0;
	  ++max_str;

	  char* end_max = (char*)strchr(max_str, (int)')');//cast to silence gcc complaint
	  if (end_max == NULL) {
	    Genetics::error(0, "Expected terminating ')' in %s.", range_dat);
	  } else {
	    *end_max = 0;
	  }
#ifdef PRINT_DEBUG_DATA
	  std::cout << "min = " << min_str << ", max = " << max_str;
#endif
	  vc_list.back().range_lo = atof(min_str);
	  vc_list.back().range_hi = atof(max_str);
	} else {
	  Genetics::error(1, "Expected both an upper and lower bound of the form (<min>, <max>) in %s", range_dat);
	}
	free(min_str);
      }
    } else if (strcmp(type_str, "uint") == 0) {
      vc_list.push_back( Genetics::VarContainer(0, 0, 0, Genetics::t_uint) );
    } else if (strcmp(type_str, "int") == 0) {
      vc_list.push_back( Genetics::VarContainer(0, 0, 0, Genetics::t_int) );
    } else if (str_starts(type_str, "bitstream") == 0) {
      const char* size_str = strchr(type_str, (int)'_');
      vc_list.push_back( Genetics::VarContainer(0, 0, 0, Genetics::t_bitstream) );
      if (size_str[0]) {
	++size_str;
	unsigned int len = atoi(size_str);
	vc_list.back().loc = len;
      }
    } else {
      Genetics::error(0, "Unrecognized type specifier %s.", type_str);
    }
  }
  map->initialize(vc_list);
#ifdef PRINT_DEBUG_DATA
  std::cout << "map_size = " << map->get_num_params() << "\n";
#endif
}

void PythonProblem::set_template_parameter(unsigned param_ind, py::object o) {
  template_set = true;
  Genetics::Type t = map->get_type(param_ind);
  if (t == Genetics::t_int || t == Genetics::t_bitstream || t == Genetics::t_uint) {
    tmplt_org.set_int(param_ind, o.cast<int>());
  } else if (t == Genetics::t_real) {
    tmplt_org.set_real(param_ind, o.cast<double>());
  }
}

PopulationWrapper* PythonProblem::initialize_population(ge::String conf_file) {
  PopulationWrapper* pop;
  if (template_set) {
    pop = new PopulationWrapper(N_BITS, N_OBJS, (Genetics::Problem*)this, &tmplt_org, map, conf_file);
  } else {
    pop = new PopulationWrapper(N_BITS, N_OBJS, (Genetics::Problem*)this, map, conf_file);
  }
#ifdef PRINT_DEBUG_DATA
  std::cout << "map_size = " << map->get_num_params() << "\n";
#endif
  return pop;
}

PYBIND11_MODULE(evo_p, m) {
  m.doc() = "pybind11 example plugin";
  py::class_<OrganismWrapper>(m, "Organism")
      .def("get_rank", &OrganismWrapper::get_rank)
      .def("get_fitness",&OrganismWrapper::get_fitness)
      .def("get_fitness_list",&OrganismWrapper::get_fitness_list)
      .def("read_real", &OrganismWrapper::read_real)
      .def("read_int", &OrganismWrapper::read_int)
      .def("read_uint", &OrganismWrapper::read_uint)
      .def("set_fitness", (void (OrganismWrapper::*)(double)) &OrganismWrapper::set_fitness)
      .def("set_fitness", (void (OrganismWrapper::*)(py::list)) &OrganismWrapper::set_fitness)
      .def("get_phenotype", &OrganismWrapper::get_phenotype);
  py::class_<PopulationWrapper>(m, "Population")
      .def("evaluate", &PopulationWrapper::evaluate)
      .def("iterate", &PopulationWrapper::iterate, py::return_value_policy::take_ownership)
      .def("get_best", &PopulationWrapper::get_best, py::return_value_policy::reference_internal, py::arg("i") = 0)
      .def("get_max_fitness", &PopulationWrapper::get_max_fitness, py::arg("i") = 0)
      .def("get_min_fitness", &PopulationWrapper::get_min_fitness, py::arg("i") = 0)
      .def("get_fitness_range", &PopulationWrapper::get_fitness_range, py::arg("i") = 0)
      .def("get_fitness_var", &PopulationWrapper::get_fitness_var, py::arg("i") = 0)
      .def("get_parent", &PopulationWrapper::get_parent)
      .def("get_parents", &PopulationWrapper::get_parents)
      .def("get_children", &PopulationWrapper::get_children)
      .def("get_pop_size", &PopulationWrapper::get_pop_size)
      .def("set_pop_size", &PopulationWrapper::set_pop_size)
      .def("get_n_survivors", &PopulationWrapper::get_survivors)
      .def("set_n_survivors", &PopulationWrapper::set_survivors)
      .def("get_num_gens", &PopulationWrapper::get_num_gens)
      .def("set_num_gens", &PopulationWrapper::set_num_gens)
      .def("get_num_crossovers", &PopulationWrapper::get_num_crossovers)
      .def("set_num_crossovers", &PopulationWrapper::set_num_crossovers)
      .def("get_param_var", &PopulationWrapper::get_init_coup_mean)
      .def("set_param_var", &PopulationWrapper::set_init_coup_mean)
      .def("get_mutate_prob", &PopulationWrapper::get_mutate_prob)
      .def("set_mutate_prob", &PopulationWrapper::set_mutate_prob)
      .def("get_crossover_prob", &PopulationWrapper::get_crossover_prob)
      .def("set_crossover_prob", &PopulationWrapper::set_crossover_prob)
      .def("get_hypermutation_threshold", &PopulationWrapper::get_hypermutation_threshold)
      .def("set_hypermutation_threshold", &PopulationWrapper::set_hypermutation_threshold);
  py::class_<PythonProblem>(m, "Problem")
      .def(py::init<int, int, int>())
      .def("initialize_population", &PythonProblem::initialize_population, py::return_value_policy::take_ownership, py::arg("conf_file") = "")
      .def("set_template_parameter", &PythonProblem::set_template_parameter)
      .def("set_parameter_range", &PythonProblem::set_parameter_range)
      .def("set_phenotype_parameters", &PythonProblem::set_phenotype_parameters)
      .def("set_fitness_function", &PythonProblem::set_fitness_function);
}
