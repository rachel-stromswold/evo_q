#include "population.h"

namespace Genetics {
 
Population::Population(_uint pn_bits, _uint pn_objs, PhenotypeMap* p_map, String conf_fname) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs)
{
  if (conf_fname != "") {
    args.initialize_from_file(conf_fname.c_str());
  }
  map = p_map;
  pop_stats = (FitnessStats*)malloc(sizeof(FitnessStats)*N_OBJS);
  this->survivors_num = args.get_survivors();
  this->offspring_num = args.get_pop_size();
  if (this->offspring_num % 2 == 0) {
    this->offspring_num++;
  }
  //we need to keep the old and new generation in separate arrays to avoid overwriting data
  this->offspring.insert(this->offspring.end(), this->offspring_num, std::shared_ptr<Organism>(NULL)); 

  //initally fill up the offspring randomly
  for (size_t i = 0; i < this->offspring_num; ++i) {
    this->old_gen.push_back(std::make_shared<Organism>(N_BITS, N_OBJS, map));
    this->old_gen[i]->randomize(&args);
  }

  for (_uint i = 0; i < N_OBJS; ++i) {
    this->pop_stats[i].max = -std::numeric_limits<double>::infinity();
    this->pop_stats[i].min = std::numeric_limits<double>::infinity();
  }

  N_PARAMS = map->get_num_params();
  var_labels = (char**)malloc(sizeof(char*)*N_PARAMS);
  for (_uint j = 0; j < map->get_num_params(); ++j) {
    var_labels[j] = (char*)malloc(sizeof(char)*OUT_BUF_SIZE);
    snprintf(var_labels[j], OUT_BUF_SIZE, "x_%d", j);
  }
  obj_labels = (char**)malloc(sizeof(char*)*N_OBJS);
  for (_uint j = 0; j < N_OBJS; ++j) {
    obj_labels[j] = (char*)malloc(sizeof(char)*OUT_BUF_SIZE);
    snprintf(obj_labels[j], OUT_BUF_SIZE - 1, "f_%d(x)", j);
  }
}

Population::Population(_uint pn_bits, _uint pn_objs, Organism* tmplt, PhenotypeMap* p_map, String conf_fname) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs)
{
  if (conf_fname != "") {
    args.initialize_from_file(conf_fname.c_str());
  }
  map = p_map;
  pop_stats = (FitnessStats*)malloc(sizeof(FitnessStats)*N_OBJS);
  this->survivors_num = args.get_survivors();
  this->offspring_num = args.get_pop_size();
  if (this->offspring_num % 2 == 0) {
    this->offspring_num++;
  }
  //we need to keep the old and new generation in separate arrays to avoid overwriting data
  this->offspring.insert(this->offspring.end(), this->offspring_num, NULL);

  this->old_gen.push_back(std::make_shared<Organism>(*tmplt));
  //initally fill up the offspring randomly
  for (size_t i = 1; i < this->offspring_num; ++i) {
    this->old_gen.push_back( std::make_shared<Organism>(N_BITS, N_OBJS, map) );
    this->old_gen.back()->randomize(&args, tmplt);
  }

  for (_uint i = 0; i < N_OBJS; ++i) {
    this->pop_stats[i].max = -std::numeric_limits<double>::infinity();
    this->pop_stats[i].min = std::numeric_limits<double>::infinity();
  }

  char buf[OUT_BUF_SIZE];
//  var_labels.resize(map->get_num_params());
  N_PARAMS = map->get_num_params();
  var_labels = (char**)malloc(sizeof(char*)*N_PARAMS);
  for (_uint j = 0; j < map->get_num_params(); ++j) {
    var_labels[j] = (char*)malloc(sizeof(char)*OUT_BUF_SIZE);
    snprintf(var_labels[j], OUT_BUF_SIZE, "x_%d", j);
  }
  obj_labels = (char**)malloc(sizeof(char*)*N_OBJS);
  for (_uint j = 0; j < N_OBJS; ++j) {
    obj_labels[j] = (char*)malloc(sizeof(char)*OUT_BUF_SIZE);
    snprintf(obj_labels[j], OUT_BUF_SIZE - 1, "f_%d(x)", j);
  }
}

Population::~Population() {
  for (size_t i = 0; i < this->offspring_num; ++i) {
    old_gen[i].reset();
  }
  if (pop_stats) {
    free(pop_stats);
  }
  if (var_labels) {
    for (_uint i = 0; i < N_PARAMS; ++i) {
      free(var_labels[i]);
    }
    free(var_labels);
  }
  if (obj_labels) {
    for (_uint i = 0; i < N_OBJS; ++i) {
      free(obj_labels[i]);
    }
    free(obj_labels);
  }
}

Population::Population(Population& o) : args(o.args) {
  N_BITS = o.get_n_bits();
  N_PARAMS = o.N_PARAMS;
  N_OBJS = o.N_OBJS;
  map = o.map;
  old_gen = o.old_gen;
  offspring = o.offspring;
  pop_stats = (FitnessStats*)malloc(sizeof(FitnessStats)*N_OBJS);
  var_labels = (char**)malloc(sizeof(char*)*N_PARAMS);
  obj_labels = (char**)malloc(sizeof(char*)*N_OBJS);

  for (_uint i = 0; i < N_OBJS; ++i) {
    pop_stats[i].min = o.pop_stats[i].min;
    pop_stats[i].max = o.pop_stats[i].max;
    pop_stats[i].var = o.pop_stats[i].var;
    size_t len = strlen(o.obj_labels[i]) + 1;
    obj_labels[i] = (char*)malloc( sizeof(char)*len);
    snprintf(obj_labels[i], len, "%s", o.obj_labels[i]);
  }
  for (_uint i = 0; i < N_PARAMS; ++i) {
    size_t len = strlen(o.var_labels[i]) + 1;
    var_labels[i] = (char*)malloc(sizeof(char)*len);
    snprintf(obj_labels[i], len, "%s", o.obj_labels[i]);
  }
}

Population& Population::operator=(Population& o) {
  int tmp_N_BITS = N_BITS;
  int tmp_N_OBJS = N_OBJS;
  PhenotypeMap* tmp_map = map;
  FitnessStats* tmp_pop_stats = pop_stats;
  char** tmp_var_labels = var_labels;
  char** tmp_obj_labels = obj_labels;
  N_BITS = o.N_BITS;
  N_OBJS = o.N_OBJS;
  map = o.map;
  pop_stats = tmp_pop_stats;
  var_labels = o.var_labels;
  obj_labels = o.obj_labels;
  o.N_BITS = tmp_N_BITS;
  o.N_OBJS = tmp_N_OBJS;
  o.map = tmp_map;
  o.pop_stats = tmp_pop_stats;
  o.var_labels = tmp_var_labels;
  o.obj_labels = tmp_obj_labels;

  old_gen = o.old_gen;
  offspring = o.offspring;
  
  return *this;
}

Population::Population(Population&& o) : map(o.map), args(std::move(o.args)) {
  N_BITS = o.get_n_bits();
  for (size_t i = 0; i < this->offspring_num; ++i) {
    old_gen[i].reset();
  }
  old_gen = std::move(o.old_gen);
  offspring = std::move(o.offspring);
  pop_stats = o.pop_stats;
  var_labels = o.var_labels;
  obj_labels = o.obj_labels;
  o.pop_stats = NULL;
  o.var_labels = NULL;
  o.obj_labels = NULL;
}

void Population::resize_population(_uint new_size) {
  size_t old_size = old_gen.size();
  if (new_size > old_size) {
    offspring.insert( offspring.end(), new_size - old_size, std::shared_ptr<Organism>(NULL) );
    old_gen.reserve(new_size);
    for (size_t i = old_size; i < old_size; ++i) {
      old_gen.push_back( std::make_shared<Organism>(N_BITS, N_OBJS, map) );
      old_gen.back()->randomize(&args);
    }
  } else {
    offspring.resize(new_size);
    old_gen.resize(new_size);
  }
  args.set_pop_size(new_size);
}

void Population::set_n_survivors(_uint new_size) {
  size_t old_size = survivors_num;
  if (new_size > old_size) {
    survivors.insert( offspring.end(), new_size - old_size, std::shared_ptr<Organism>(NULL) );
  } else {
    survivors.resize(new_size);
  }
  args.set_survivors(new_size);
}

void Population::evaluate(Problem* prob) {
  if (N_OBJS == 1) {
    // we need to do the first loop once to initialize the minimum
    old_gen[0]->evaluate_fitness(prob);
    best_organism_ind = 0;

    pop_stats[0].max = old_gen[0]->get_fitness(0);
    pop_stats[0].min = pop_stats[0].max;
   
    for (size_t i = 1; i < this->offspring_num; ++i) {
      old_gen[i]->apply_penalty(0);
      if ( args.verbose() ) {
        std::cout << "Now evaluating organism " << i << std::endl;
      }
      if (old_gen.size() == 0) {
        std::cout << "index = " << i;
        error(1, "Tried to evaluate fitness for null organism");
      }
      old_gen[i]->evaluate_fitness(prob);

      // update the max and min fitnesses if we need to
      if (old_gen[i]->get_fitness(0) > pop_stats[0].max && !old_gen[i]->penalized()) {
        pop_stats[0].max = old_gen[i]->get_fitness(0);
        best_organism = old_gen[i];
        best_organism_ind = i;
      }

      if (old_gen[i]->get_fitness(0) < pop_stats[0].min) {
        pop_stats[0].min = old_gen[i]->get_fitness(0);
      }
    }
  } else {
    //TODO: figure out what the default behavior should be
  }
  double penalty_fact;
  if (pop_stats[0].max > 0) {
    penalty_fact = pop_stats[0].min;
  } else {
    penalty_fact = -pop_stats[0].min;
  }
  for (size_t i = 0; i < this->offspring_num; ++i) {
    if (old_gen[i]->penalized()) {
      old_gen[i]->set_fitness(pop_stats[0].min - penalty_fact*old_gen[i]->get_penalty());
    }
  }
}

void Population::find_best_organism() {
  for (int j = 0; j < N_OBJS; ++j) {
    pop_stats[j].max = old_gen[0]->get_fitness(j);
    pop_stats[j].min = old_gen[0]->get_fitness(j);
    pop_stats[j].mean = old_gen[0]->get_fitness(j) / offspring_num;
    best_organism_ind = 0;
    for (size_t i = 1; i < offspring_num; ++i) {
      double fitness_i = old_gen[i]->get_fitness(j);
      if (fitness_i > pop_stats[j].max) {
	pop_stats[j].max = fitness_i;
	best_organism_ind = i;
      }
      if (fitness_i < pop_stats[j].min) {
	pop_stats[j].min = fitness_i;
      }
      pop_stats[j].mean += fitness_i / offspring_num;
    }
    //calculate the variance
    pop_stats[j].var = 0;
    for (size_t i = 0; i < offspring_num; ++i) {
      double fitness_i = old_gen[i]->get_fitness(j);
      pop_stats[j].var += (fitness_i - pop_stats[j].mean)*(fitness_i - pop_stats[j].mean);
    }
  }
}

void Population::evaluate_async(Problem* prob) {
}

void Population::swap_orgs(int i, int j) {
  std::shared_ptr<Organism> tmp = old_gen[j];
  old_gen[j] = old_gen[i];
  old_gen[i] = tmp;
}

int Population::partition(unsigned int fit_ind,
					  std::vector<std::shared_ptr<Organism>>* work_arr,
					  int s, int e) {
  double p = (*work_arr)[e]->get_fitness(fit_ind);
  int i = s;
  for (int j = s; j < e; ++j) {
    if ((*work_arr)[j]->get_fitness(fit_ind) > p) {
      std::shared_ptr<Organism> tmp = (*work_arr)[i];
      (*work_arr)[i] = (*work_arr)[j];
      (*work_arr)[j] = tmp;
      i++;
    }
  }
  if (i != e) {
    std::shared_ptr<Organism> tmp = (*work_arr)[i];
    (*work_arr)[i] = (*work_arr)[e];
    (*work_arr)[e] = tmp;
  }
  return i;
}

void Population::sort_orgs(unsigned int fit_ind,
					   std::vector<std::shared_ptr<Organism>>* work_arr,
					   int s, int e) {
  sort_org_calls += 1;
  if (sort_org_calls > work_arr->size()*5) {
    error(0, "Suspiciously many calls (%u) to sort org have been made", sort_org_calls);
  }
  if (s == DEF_SORT_PARAM || e == DEF_SORT_PARAM) {
    sort_org_calls = 0;
    sort_orgs(fit_ind, work_arr, 0, work_arr->size() - 1);
  } else if (s < e) {
    int p = partition(fit_ind, work_arr, s, e);
    sort_orgs(fit_ind, work_arr, s, p-1);
    sort_orgs(fit_ind, work_arr, p+1, e);
  }
}

bool Population::cull_in_place() {
  size_t j = 0;

  double difference = pop_stats[0].max - pop_stats[0].min;
  //avoid divide by 0
  if (difference == 0) {
    return true;
  }
  std::uniform_real_distribution<double> dist(0, difference);

  if (survivors.size() < survivors_num) {
    survivors.resize(survivors_num);
  }
  for (size_t i = 0; i < this->offspring_num && j < survivors_num; ++i) {
    /* if M_f is the maximum fitness and m_f is the minimum the minimum, while x is the
     * fitness of a given organism, then the probability of survival is x/(M_f-m_f) or 1
     * if M_f = m_f */
    if (dist(args.get_generator()) < old_gen[i]->get_fitness(0) - pop_stats[0].min) {
      survivors[j] = old_gen[i];
      j++;
    }
  }

  std::uniform_int_distribution<int> ind_dist(0, this->offspring_num - 1);
  while (j < survivors_num) {
    survivors[j] = old_gen[ind_dist( args.get_generator() )];
    j++;
  }
  return false;
}

bool Population::cull() {
  sort_orgs(0, &old_gen);
  //if all the organisms have the same fitness then reinitialize the population
  if (old_gen[0]->get_fitness(0) - old_gen[this->offspring_num-1]->get_fitness(0) < 0.001 ) {
    return true;
  }
  double min_fit = old_gen[this->offspring_num-1]->get_fitness(0);
  double total_fit = 0;
  for (size_t i = 0; i < this->offspring_num; ++i) {
    total_fit += old_gen[i]->get_fitness(0) - min_fit;
  }
  std::uniform_real_distribution<double> dist(0, total_fit);

  //maintain a list of organisms that have already been added so no organism appears twice
  size_t* banned = (size_t*)malloc(sizeof(size_t)*survivors_num);
  for (size_t i = 0; i < survivors_num; ++i) { banned[i] = -1; }

  for (size_t i = 0; i < survivors_num; ++i) {
    double val = dist( args.get_generator() );
    size_t j = 0;
    while (val > 0.0) {
      val -= old_gen[j]->get_fitness(0) - min_fit;
      j++;
    }
    j -= 1;
    while (contains<size_t>(banned, survivors_num, j)) {
      j++;
      if (j >= survivors_num) {
        j = 0;
      }
    }
    survivors[i] = old_gen[j];
    banned[i] = j;
  }
  if (args.verbose()) {
    std::cout << "Added orgs to survs:";
    for (size_t i = 0; i < survivors_num; ++i) {
      std::cout << " " << banned[i];
    }
    std::cout << std::endl;
  }
  free(banned);
  best_organism_ind = 0;
  return false;
}

void Population::breed_shuffle() {
  std::uniform_int_distribution<size_t> dist_surv0(0, survivors_num - 1);
  std::vector<Organism*> children;
  Organism** shuffled_inds = (Organism**)malloc(sizeof(Organism*)*offspring_num);
  for (size_t i = 0; i < survivors_num; ++i) {
    shuffled_inds[i] = survivors[i].get();
  }
  for (size_t i = 0; i < survivors_num; ++i) {
    size_t ind_o = dist_surv0( args.get_generator() );
    Organism* tmp = shuffled_inds[i];
    shuffled_inds[i] = shuffled_inds[ind_o];
    shuffled_inds[ind_o] = tmp;
  }
  size_t last_org_ind = offspring_num;
  for (size_t i = survivors_num; i < offspring_num; ++i) {
    size_t org_ind= dist_surv0( args.get_generator() );
    //ensure that we don't see the same organism breeding with itself
    if (org_ind == last_org_ind) {
      org_ind = (org_ind + 1) % offspring_num;
    }
    shuffled_inds[i] = survivors[org_ind].get();
    last_org_ind = org_ind;
  }
  if (this->offspring_num % 2 == 1) {
    //elitist algorithm, make the first individual in the next generation the previous most fit
    offspring[0] = old_gen[best_organism_ind];

    for (size_t i = 1; 2*i < this->offspring_num; i++) {
      children = shuffled_inds[2*i - 1]->breed(&args, shuffled_inds[2*i]);
      offspring[2*i] = std::shared_ptr<Organism>(children[0]);
      offspring[2*i - 1] = std::shared_ptr<Organism>(children[1]);
    }
  } else {
    for (size_t i = 0; 2*i + 1 < this->offspring_num; i++) {
      children = shuffled_inds[2*i]->breed(&args, shuffled_inds[2*i + 1]);
      offspring[2*i] = std::shared_ptr<Organism>(children[0]);
      offspring[2*i + 1] = std::shared_ptr<Organism>(children[1]);
    }
  }

  offspring.swap(old_gen);
}

void Population::breed() {
  find_best_organism();
  std::uniform_int_distribution<size_t> dist_surv0(0, survivors_num - 1);
  std::uniform_int_distribution<size_t> dist_surv1(0, survivors_num - 2);
  std::vector<Organism*> children;
  size_t* shuffled_inds = (size_t*)malloc(sizeof(size_t)*survivors_num);
  if (this->offspring_num % 2 == 1) {
    offspring[0] = old_gen[best_organism_ind];

    for (size_t i = 1; 2*i < this->offspring_num; i++) {
      size_t par1_i = dist_surv0( args.get_generator() );
      //use the survivors_num-1 distribution to guarantee different parents
      size_t par2_i = dist_surv1( args.get_generator() );
      if (par2_i >= par1_i) {
	par2_i++;
      }
      children = survivors[par1_i].get()->breed(&args, survivors[par2_i].get());
      offspring[2*i] = std::shared_ptr<Organism>(children[0]);
      offspring[2*i - 1] = std::shared_ptr<Organism>(children[1]);
    }
  } else {
    for (size_t i = 0; 2*i + 1 < this->offspring_num; i++) {
      size_t par1_i = dist_surv0( args.get_generator() );
      //use the survivors_num-1 distribution to guarantee different parents
      size_t par2_i = dist_surv1( args.get_generator() );
      if (par2_i >= par1_i) {
	par2_i++;
      }
      children = survivors[par1_i].get()->breed(&args, survivors[par2_i].get());
      offspring[2*i] = std::shared_ptr<Organism>(children[0]);
      offspring[2*i + 1] = std::shared_ptr<Organism>(children[1]);
    }
  }

  free(shuffled_inds);
  offspring.swap(old_gen);
}

bool Population::iterate(ConvergenceCriteria* conv) {
  if (cull_in_place()) {
    return true;
  }
  breed(); 
  generation++;
  if (conv) {
    return conv->evaluate_convergence(pop_stats);
  } else {
    return (generation < args.get_num_gens());
  }
}

std::shared_ptr<Organism> Population::get_best_organism(size_t i) {
  if (best_organism_ind >= old_gen.size()) {
    find_best_organism();
  }
  if (i == 0) {
    return old_gen[best_organism_ind];
  } else {
    //TODO: implement a binary private member variable that tracks whether sorting needs to be performed for the sake of efficiency
    sort_orgs(0, &old_gen);
    return old_gen[i];
  }
}

std::shared_ptr<Organism> Population::get_organism(size_t i) {
  if (i > old_gen.size())
    error(1, "Attempt to access invalid index %d when the maximum allowed is %d.", i, old_gen.size());
  return old_gen[i];
}

std::shared_ptr<Organism> Population::get_child(size_t i) {
  if (i > offspring.size())
    error(1, "Attempt to access invalid index %d when the maximum allowed is %d.", i, offspring.size());
  return offspring[i];
}

void Population::run(Problem* prob) {
  evaluate(prob);
  if ( args.wait_for_con() ) {
    size_t i = 1;
    unsigned streak = 0;
    double prev_ftns = get_best_organism()->get_fitness(0);
    while (i < MAX_NUM_GENS) {
      if (iterate()) {
        break;
      }
      evaluate(prob);

      if (get_best_organism()->get_fitness(0) > prev_ftns) {
        prev_ftns = get_best_organism()->get_fitness(0);
        streak = 0;
      }

      streak++;
      i++;
      if (streak > args.get_num_gens()) {
        break;
      }
    }

    if (i >= MAX_NUM_GENS) {
      std::cout << "failed to converge after " << i << " generations." << std::endl;
    } else {
      std::cout << "converged to result after " << i << " generations." << std::endl;
    }
  } else {
    if (args.verbose()) {
      std::cout << "Now evaluating generation 0..." << std::endl;
    }
    for (size_t i = 1; i < args.get_num_gens() + 1; ++i) {
      if (args.verbose()) {
        std::cout << "Now evaluating generation " << i << "..." << std::endl;
      }
      //produce the next generation in the population
      if (iterate()) {
        break;
      }
      evaluate(prob);
    } 
  }
}

void Population::set_var_label(_uint ind, String val) {
  if (ind >= N_PARAMS) {
    error(0, "Invalid parameter index %u provided for set_var_label. The parameter index must be less than %u.", ind, N_PARAMS);
  } else {
    size_t vs = val.size() + 1;
    if (vs > OUT_BUF_SIZE) {
      free(var_labels[ind]);
      var_labels[ind] = (char*)malloc(sizeof(char)*vs);
    }
    for (int j = 0; j < vs; ++j) {
      var_labels[ind][j] = val[j];
    }
    var_labels[ind][vs - 1] = 0;
  }
}

void Population::set_obj_label(_uint ind, String val) {
  if (ind >= N_OBJS) {
    error(0, "Invalid objective index %u provided for set_obj_label. The objective index must be less than %u.", ind, N_OBJS);
  } else {
    size_t vs = val.size() + 1;
    if (vs > OUT_BUF_SIZE) {
      free(obj_labels[ind]);
      obj_labels[ind] = (char*)malloc(sizeof(char)*vs);
    }
    for (int j = 0; j < vs; ++j) {
      obj_labels[ind][j] = val[j];
    }
    obj_labels[ind][vs - 1] = 0;
  }
}

Vector<String> Population::get_header() {
  String def;
  Vector<String> ret;
  ret.reserve(old_gen.size()*(N_PARAMS + N_OBJS));

  for (_uint i = 0; i < old_gen.size(); ++i) {
    for (_uint j = 0; j < N_PARAMS; ++j) {
      ret.push_back(String(var_labels[j]));
    }

    for (_uint j = 0; j < N_OBJS; ++j) {
      ret.push_back(String(obj_labels[j]));
    }
  }
  return ret;
}

Vector<String> Population::get_pop_data() {
  sort_orgs(0, &old_gen);
  _uint span = N_OBJS + map->get_num_params();
  String def;
  Vector<String> ret(offspring_num*span, def);
  char buf[OUT_BUF_SIZE];

  for (_uint i = 0; i < old_gen.size(); ++i) {
    for (_uint j = 0; j < Population::map->get_num_params(); ++j) {
      ret[i*span + j] = old_gen[i]->get_chromosome_string(j);
    }

    for (_uint j = 0; j < N_OBJS; ++j) {
      snprintf(buf, OUT_BUF_SIZE - 1, "%f", old_gen[i]->get_fitness(j));
      ret[i*span + map->get_num_params() + j] = buf;
    }
  }

  return ret;
}

// =============================== POPULATION_NGSAII ===============================

Population_NSGAII::~Population_NSGAII() {
  for (size_t i = 0; i < Population::offspring_num; ++i) {
    if (Population::old_gen[i] != NULL) {
      Population::old_gen[i].reset();
    }
    if (Population::offspring[i] != NULL) {
      Population::offspring[i].reset();
    }
  }
}

void Population_NSGAII::hypermutate() {
  sort_orgs(get_n_objs(), &(pareto_fronts[0]));
  //select the half of the most crowded individuals in the first front
  for (_uint i = pareto_fronts[0].size() - 1; i > pareto_fronts[0].size()/2; --i) {
    pareto_fronts[0][i]->randomize(&args);
  }
}

void Population_NSGAII::update_fitness_ranges(Organism* org, size_t i) {
  if (org->penalized()) {
    max_penalty_ind = i + 1;
    if (min_penalty_ind == Population::offspring_num) {
      min_penalty_ind = i;
    }
  } else {
    for (_uint j = 0; j < this->get_n_objs(); ++j) {
      if (org->get_fitness(j) > Population::pop_stats[j].max) {
        Population::pop_stats[j].max = org->get_fitness(j);
      }
      if (org->get_fitness(j) < Population::pop_stats[j].min) {
        Population::pop_stats[j].min = org->get_fitness(j);
      }
    }
  }
}

void Population_NSGAII::evaluate(Problem* prob) {
  std::vector<std::shared_ptr<Organism>> empty;
  pareto_fronts.push_back(empty);

  std::vector<std::shared_ptr<Organism>> cmb_arr = Population::old_gen;
  //initialize max and min fitness values
  for (_uint j = 0; j < this->get_n_objs(); ++j) {
    pop_stats[j].max = -std::numeric_limits<double>::infinity();
    pop_stats[j].min = std::numeric_limits<double>::infinity();
  }
  min_penalty_ind = Population::offspring_num, max_penalty_ind = 0;
  for (size_t i = 0; i < Population::offspring_num; ++i) {
    //fitnesses for the old generation are leftover and already have applied penalties, ensure we don't apply these penalties again
    Population::old_gen[i]->apply_penalty(0);
    if (Population::offspring[i] != NULL) {
      Population::offspring[i]->apply_penalty(0);
      Population::offspring[i]->evaluate_fitness(prob);
      Population::offspring[i]->distance = 0;
      //to maintain elitism we look at both the parent and offspring generations
      cmb_arr.push_back(Population::offspring[i]);
      update_fitness_ranges(Population::offspring[i].get(), Population::offspring_num + i);
    } else {
      Population::old_gen[i]->evaluate_fitness(prob);
      Population::old_gen[i]->distance = 0;//initialize for later crowding calculations
      update_fitness_ranges(Population::old_gen[i].get(), i);
    }
  }
  //apply penalties to each organism in the combined array
  for (size_t i = min_penalty_ind; i < max_penalty_ind; ++i) {
    if (cmb_arr[i]->penalized()) {
      for (size_t j = 0; j < this->get_n_objs(); ++j) {
        double n_fit = pop_stats[j].min;
        if (pop_stats[j].max > 0) {
          n_fit -= pop_stats[j].max*(cmb_arr[i]->get_penalty());
        } else {
          n_fit += pop_stats[j].min*(cmb_arr[i]->get_penalty());
        }
        cmb_arr[i]->set_fitness(n_fit);
      }
    }
  }
  std::vector<std::vector<std::shared_ptr<Organism>>> dominating(cmb_arr.size(), empty);

  for (size_t i = 0; i < cmb_arr.size(); ++i) {
    cmb_arr[i]->n_dominations = 0;
    for (size_t j = 0; j < cmb_arr.size(); ++j) {
      if (i != j) {
        //if the ith solution dominates the jth add the jth entry to the list of dominated solutions, otherwise increment the number of dominating solutions
        if ( cmb_arr[i].get()->dominates(cmb_arr[j].get()) ) {
          dominating[i].push_back(cmb_arr[j]);
        } else if ( cmb_arr[j].get()->dominates(cmb_arr[i].get()) ) {
          cmb_arr[i]->n_dominations++;
        }
      }
    }
    if (cmb_arr[i]->n_dominations == 0) {
      cmb_arr[i]->rank = 0;
      pareto_fronts[0].push_back(cmb_arr[i]);
    }
  }

  size_t i = 0;
  while (i < pareto_fronts.size() && pareto_fronts[i].size() != 0) {
    pareto_fronts.push_back(empty);
    for (size_t j = 0; j < pareto_fronts[i].size(); ++j) {
      for (size_t k = 0; k < cmb_arr.size(); ++k) {
        cmb_arr[k]->n_dominations--;
        if (cmb_arr[k]->n_dominations == 0) {
          cmb_arr[k]->rank = i + 1;
          pareto_fronts[i + 1].push_back(cmb_arr[k]);
        }
      }
    }
    i++;
  }
}

bool Population_NSGAII::iterate() {
  //if we exceed the threshold, start hypermutation
  if ((double)(pareto_fronts[0].size())/offspring_num > args.get_hypermutation_threshold()) {
    hypermutate();
  }
  size_t i = 0;
  std::vector<std::shared_ptr<Organism>> tmp(Population::offspring_num, NULL);
  _uint k = 0;
  while (i < pareto_fronts.size() && (k + pareto_fronts[i].size()) <= Population::offspring_num) {
    for (size_t j = 0; j < pareto_fronts[i].size(); ++j) {
      tmp[k] = pareto_fronts[i][j];
      ++k;
    }
    ++i;
  }
  //fill in the last elements from the remaining pareto front ranked according to crowding
  if (k < Population::offspring_num) {
    for (size_t ii = 0; ii < pareto_fronts[i].size(); ++ii) {
      pareto_fronts[i][ii]->distance = 0;
    }
    //sort by each objective function for crowding evaluation
    for (size_t j = 0; j < get_n_objs(); ++j) {
      Population::sort_orgs(j, &(pareto_fronts[i]));

      //the highest and lowest values should be considered to have no crowding
      pareto_fronts[i].front()->distance = std::numeric_limits<double>::infinity();
      pareto_fronts[i].back()->distance = std::numeric_limits<double>::infinity();
      double range = pareto_fronts[i].front()->get_fitness(j) - pareto_fronts[i].back()->get_fitness(j);
      if (range == 0) {
        for (size_t ii = 0; ii < pareto_fronts[i].size(); ++ii) {
          pareto_fronts[i][ii]->distance = 0;
        }
      } else {
        for (size_t ii = 1; ii < pareto_fronts[i].size() - 1; ++ii) {
          if (pareto_fronts[i][ii]->distance != std::numeric_limits<double>::infinity()) {
            double sj_h = pareto_fronts[i][ii-1]->get_fitness(j);
            double sj_l = pareto_fronts[i][ii+1]->get_fitness(j);
            double d_norm = (sj_h - sj_l) / range;
            pareto_fronts[i][ii]->distance += d_norm;
          }
        }
      }
    }
    Population::sort_orgs(get_n_objs(), &(pareto_fronts[i]));
    size_t j = 0;
    //select the least crowded individuals
    size_t p_i_size = pareto_fronts[i].size();
    while (j < p_i_size && k < Population::offspring_num) {
      tmp[k] = pareto_fronts[i][j];
      ++j;
      ++k;
    }
    //we don't need the rest of the elements so they should be deleted
    while (j < pareto_fronts[i].size()) {
      if (pareto_fronts[i][j] != NULL) {
        pareto_fronts[i][j].reset();
        pareto_fronts[i][j] = NULL;
      }
      ++j;
    }
    ++i;
  }
  while (i < pareto_fronts.size()) {
    for (size_t j = 0; j < pareto_fronts[i].size(); ++j) {
      if (pareto_fronts[i][j] != NULL) {
        pareto_fronts[i][j].reset();
        pareto_fronts[i][j] = NULL;
      }
    }
    ++i;
  }
  Population::old_gen = tmp;
  //ensure that we aren't keeping null pointers around for safety
  pareto_fronts.clear();

  breed();
  return false;
}

void Population_NSGAII::breed() {
  _uint arena_size = args.get_survivors();
  SampleDraw sampler(Population::offspring_num, arena_size);
  std::uniform_int_distribution<_uint> selector(0, Population::offspring_num - 1);
  std::vector<Organism*> children;

  for (size_t i = 0; i < Population::offspring_num / 2; ++i) {
    Organism *first_parent, *second_parent;
    std::vector<_uint> t1 = sampler( args.get_generator() );
    std::vector<_uint> t2 = sampler( args.get_generator() );
    first_parent = Population::old_gen[t1[0]].get();
    second_parent = Population::old_gen[t2[0]].get();
    for (size_t j = 1; j < t1.size(); ++j) {
      if (t1[j] >= Population::old_gen.size()) {
        error(1, "Invalid index created from random sampling, %d", t1[j]);
      }
      if (Population::old_gen[t1[j]]->get_rank() < first_parent->get_rank()) {
        first_parent = Population::old_gen[t1[0]].get();
      }
      if (Population::old_gen[t2[j]]->get_rank() < second_parent->get_rank()) {
        second_parent = Population::old_gen[t2[0]].get();
      }
    }
    //ensure that we don't use the same parent twice with reasonable probably
    if (first_parent == second_parent) {
      second_parent = Population::old_gen[selector( args.get_generator() )].get();
    }

    children = first_parent->breed(&args, second_parent);
    Population::offspring[2*i] = std::shared_ptr<Organism>(children[0]);
    if (2*i + 1 < Population::offspring_num) {
      Population::offspring[2*i + 1] = std::shared_ptr<Organism>(children[1]);
    }
  }
  Population::old_gen.swap(Population::offspring);
}

Vector<String> Population_NSGAII::get_header() {
  Vector<String> ret;
  char buf[OUT_BUF_SIZE];
  for (_uint i = 0; i < 2*Population::old_gen.size(); ++i) {
    snprintf(buf, OUT_BUF_SIZE, "rank_%d", i);
    ret.push_back( String(buf) );

    for (_uint j = 0; j < Population::map->get_num_params(); ++j) {
      ret.push_back( String(Population::var_labels[j]) );
    }

    for (_uint j = 0; j < get_n_objs(); ++j) {
      ret.push_back( String(Population::obj_labels[j]) );
    }
  }
  return ret;
}

Vector<String> Population_NSGAII::get_pop_data() {
  String tmp_str;
  Vector<String> ret;
  char buf[OUT_BUF_SIZE];
  for (_uint n = 0; n < pareto_fronts.size(); ++n) {
    for (_uint i = 0; i < pareto_fronts[n].size(); ++i) {
      snprintf(buf, OUT_BUF_SIZE - 1, "%d", n);
      tmp_str = buf;
      ret.push_back(tmp_str);

      for (_uint j = 0; j < Population::map->get_num_params(); ++j) {
        ret.push_back( pareto_fronts[n][i]->get_chromosome_string(j) );
      }

      for (_uint j = 0; j < get_n_objs(); ++j) {
        snprintf(buf, OUT_BUF_SIZE - 1, "%f", pareto_fronts[n][i]->get_fitness(j));
        ret.push_back(String(buf));
      }
    }
  }
  return ret;
}

}
