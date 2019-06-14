#include "population.h"

namespace Genetics {
  
Population::Population(_uint pn_bits, _uint pn_objs, PhenotypeMap* p_map, String conf_fname) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs)
{
  map = p_map;
  min_fitness = (double*)malloc(sizeof(double)*N_OBJS);
  max_fitness = (double*)malloc(sizeof(double)*N_OBJS);
  this->survivors_num = args.get_survivors();
  this->offspring_num = args.get_pop_size();
  //we need to keep the old and new generation in separate arrays to avoid overwriting data
  this->offspring.insert(this->offspring.end(), this->offspring_num, std::shared_ptr<Organism>(NULL)); 

  //initally fill up the offspring randomly
  for (size_t i = 0; i < this->offspring_num; ++i) {
    this->old_gen.push_back(std::make_shared<Organism>(N_BITS, N_OBJS, map));
    this->old_gen[i]->randomize(&args);
  }

  for (_uint i = 0; i < N_OBJS; ++i) {
    this->max_fitness[i] = 0;
  }

  var_labels.resize( map->get_num_params() );
  obj_labels.resize( N_OBJS );
#ifdef USE_CUSTOM_CONTAINERS
  for (_uint j = 0; j < map->get_num_params(); ++j) {
    var_labels[j] << "x_" << j;
  }
  for (_uint j = 0; j < N_OBJS; ++j) {
    obj_labels[j] << "f_" << j << "(x)";
  }
#else
  char* buf = (char*)malloc(sizeof(char)*(3 + map->get_num_params() + N_OBJS));
  for (_uint j = 0; j < map->get_num_params(); ++j) {
    sprintf(buf, "x_%d", j);
    var_labels[j] = buf;
  }
  for (_uint j = 0; j < N_OBJS; ++j) {
    snprintf(buf, 2 + map->get_num_params() + N_OBJS, "f_%d(x)", j);
    obj_labels[j] = buf;
  }
  free(buf);
#endif
  if (conf_fname != "") {
    args.initialize_from_file(conf_fname.c_str());
  }
}

Population::Population(_uint pn_bits, _uint pn_objs, Organism* tmplt, PhenotypeMap* p_map, String conf_fname) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs)
{
  map = p_map;
  min_fitness = (double*)malloc(sizeof(double)*N_OBJS);
  max_fitness = (double*)malloc(sizeof(double)*N_OBJS);
  this->survivors_num = args.get_survivors();
  this->offspring_num = args.get_pop_size();
  //we need to keep the old and new generation in separate arrays to avoid overwriting data
  this->offspring.insert(this->offspring.end(), this->offspring_num, NULL);

  this->old_gen.push_back(std::make_shared<Organism>(*tmplt));
  //initally fill up the offspring randomly
  for (size_t i = 1; i < this->offspring_num; ++i) {
    this->old_gen.push_back( std::make_shared<Organism>(N_BITS, N_OBJS, map) );
    this->old_gen.back()->randomize(&args, tmplt);
  }

  for (_uint i = 0; i < N_OBJS; ++i) {
    this->max_fitness[i] = 0;
  }

  char buf[OUT_BUF_SIZE];
  var_labels.resize(map->get_num_params());
  for (_uint j = 0; j < map->get_num_params(); ++j) {
    snprintf(buf, OUT_BUF_SIZE - 1, "x_%d", j);
    var_labels[j] = String(buf);
  }
  obj_labels.resize( N_OBJS );
  for (_uint j = 0; j < N_OBJS; ++j) {
    snprintf(buf, OUT_BUF_SIZE - 1, "f_%d(x)", j);
    obj_labels[j] = String(buf);
  }
  if (conf_fname != "") {
    args.initialize_from_file(conf_fname.c_str());
  }
}

Population::~Population() {
  for (size_t i = 0; i < this->offspring_num; ++i) {
    old_gen[i].reset();
  }
  if (min_fitness) {
    free(min_fitness);
  }
  if (max_fitness) {
    free(max_fitness);
  }
}

Population::Population(Population& o) {
  N_BITS = o.get_n_bits();
  N_OBJS = o.N_OBJS;
  map = o.map;
  old_gen = o.old_gen;
  offspring = o.offspring;
  min_fitness = (double*)malloc(sizeof(double)*N_OBJS);
  max_fitness = (double*)malloc(sizeof(double)*N_OBJS);
  for (_uint i = 0; i < N_OBJS; ++i) {
    min_fitness[i] = o.min_fitness[i];
    max_fitness[i] = o.max_fitness[i];
  }
}

Population& Population::operator=(Population& o) {
  int tmp_N_BITS = N_BITS;
  int tmp_N_OBJS = N_OBJS;
  PhenotypeMap* tmp_map = map;
  double* tmp_min_fit = min_fitness;
  double* tmp_max_fit = max_fitness;
  N_BITS = o.N_BITS;
  N_OBJS = o.N_OBJS;
  map = o.map;
  min_fitness = o.min_fitness;
  max_fitness = o.max_fitness;
  o.N_BITS = tmp_N_BITS;
  o.N_OBJS = tmp_N_OBJS;
  o.map = tmp_map;
  o.min_fitness = tmp_min_fit;
  o.max_fitness = tmp_max_fit;

  old_gen = o.old_gen;
  offspring = o.offspring;
  
  return *this;
}

Population::Population(Population&& o) : map(o.map) {
  N_BITS = o.get_n_bits();
  for (size_t i = 0; i < this->offspring_num; ++i) {
    old_gen[i].reset();
  }
  old_gen = std::move(o.old_gen);
  offspring = std::move(o.offspring);
  if (N_OBJS != o.get_n_objs()) {
    N_OBJS = o.get_n_objs();

    for (_uint i = 0; i < N_OBJS; ++i) {
      min_fitness[i] = o.min_fitness[i];
      max_fitness[i] = o.max_fitness[i];
    }
  }
  min_fitness = o.min_fitness;
  max_fitness = o.max_fitness;
}

void Population::evaluate(Problem* prob) {
  if (N_OBJS == 1) {
    // we need to do the first loop once to initialize the minimum
    old_gen[0]->evaluate_fitness(prob);
    best_organism_ind = 0;

    max_fitness[0] = old_gen[0]->get_fitness(0);
    min_fitness[0] = max_fitness[0];
   
    for (size_t i = 1; i < this->offspring_num; ++i) {
      if ( args.verbose() ) {
        std::cout << "Now evaluating organism " << i << std::endl;
      }
      if (old_gen.size() == 0) {
        std::cout << "index = " << i;
        error(1, "Tried to evaluate fitness for null organism");
      }
      old_gen[i]->evaluate_fitness(prob);

      // update the max and min fitnesses if we need to
      if (old_gen[i]->get_fitness(0) > max_fitness[0]) {
        max_fitness[0] = old_gen[i]->get_fitness(0);
        best_organism = old_gen[i];
        best_organism_ind = i;
      }

      if (old_gen[i]->get_fitness(0) < min_fitness[0]) {
        min_fitness[0] = old_gen[i]->get_fitness(0);
      }
    }
  } else {
    //TODO: figure out what the default behavior should be
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

  double difference = max_fitness - min_fitness;
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
    if (dist(args.get_generator()) < old_gen[i]->get_fitness(0)-min_fitness[0]) {
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

void Population::breed() {
  std::uniform_int_distribution<size_t> dist_surv0(0, survivors_num - 1);
  std::uniform_int_distribution<size_t> dist_surv1(0, survivors_num - 2);
  std::vector<Organism*> children;
  if (this->offspring_num % 2 == 1) {
    //elitist algorithm, make the first individual in the next generation the previous most fit
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
  //get rid of the old generation
  size_t i = 0;
  for (; i < best_organism_ind; ++i) {
    old_gen[i].reset();
  }
  i++;//skip over the best organism
  for (; i < this->offspring_num; ++i) {
    old_gen[i].reset();
  }

  offspring.swap(old_gen);
}

bool Population::iterate() {
  if (cull_in_place()) {
    return true;
  }
  breed();
  return false;
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

Vector<String> Population::get_header() {
  String def;
  Vector<String> ret;

  for (_uint i = 0; i < old_gen.size(); ++i) {
    for (_uint j = 0; j < map->get_num_params(); ++j) {
      ret.push_back(var_labels[j]);
    }

    for (_uint j = 0; j < N_OBJS; ++j) {
      ret.push_back(obj_labels[j]);
    }
  }
  return ret;
}

Vector<String> Population::get_pop_data() {
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

void Population_NSGAII::evaluate(Problem* prob) {
  std::vector<std::shared_ptr<Organism>> empty;
  pareto_fronts.push_back(empty);

  std::vector<std::shared_ptr<Organism>> cmb_arr = Population::old_gen;
  for (size_t i = 0; i < Population::offspring_num; ++i) {
    Population::old_gen[i]->evaluate_fitness(prob);
    Population::old_gen[i]->distance = 0;//initialize for later crowding calculations
    if (Population::offspring[i] != NULL) {
      Population::offspring[i]->evaluate_fitness(prob);
      Population::offspring[i]->distance = 0;
      //to maintain elitism we look at both the parent and offspring generations
      cmb_arr.push_back(Population::offspring[i]);
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
}

Vector<String> Population_NSGAII::get_header() {
  Vector<String> ret;
  char buf[OUT_BUF_SIZE];
  for (_uint i = 0; i < 2*Population::old_gen.size(); ++i) {
    snprintf(buf, OUT_BUF_SIZE - 1, "rank_%d", i);
    ret.push_back(String(buf));

    for (_uint j = 0; j < Population::map->get_num_params(); ++j) {
      ret.push_back(Population::var_labels[j]);
    }

    for (_uint j = 0; j < get_n_objs(); ++j) {
      ret.push_back(Population::obj_labels[j]);
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
