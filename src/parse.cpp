#include "parse.h"

namespace Genetics {

ArgStore::ArgStore() : generator() {
  pop_size = DEF_POP_SIZE;
  breed_pop_size = DEF_BREED_POP_SIZE;
  num_gens = DEF_NUM_GENS;
  num_crossovers = DEF_NUM_CROSSOVERS;
  init_coup_var = DEF_COUP_VAR;
  init_coup_mean = DEF_COUP_MEAN;
  crossover_prob = DEF_CROSSOVER_PROB;
  mutate_prob = DEF_MUTATE_PROB;
  hypermutation_threshold = DEF_HYPER_THRESH;
  replacement_fraction = DEF_REPLACE_FRAC;
  long_bin = NULL;
  short_bin = NULL;
  bern_mut = NULL;
  bern_cross = NULL;
}

ArgStore::ArgStore(const ArgStore& o) {
  pop_size = o.pop_size;
  breed_pop_size = o.breed_pop_size;
  num_gens = o.num_gens;
  num_crossovers = o.num_crossovers;
  init_coup_var = o.init_coup_var;
  init_coup_mean = o.init_coup_mean;
  crossover_prob = o.crossover_prob;
  mutate_prob = o.mutate_prob;
  hypermutation_threshold = o.hypermutation_threshold;
  replacement_fraction = DEF_REPLACE_FRAC;
  if (o.long_bin) {
    long_bin = new std::binomial_distribution<unsigned char>(*o.long_bin);
  } else {
    long_bin = NULL;
  }
  if (o.short_bin) {
    short_bin = new std::binomial_distribution<unsigned char>(*o.short_bin);
    short_bin_n = o.short_bin_n;
  } else {
    short_bin = NULL;
    short_bin_n = 0;
  }
  if (o.bern_mut) {
    bern_mut = new std::bernoulli_distribution(*o.bern_mut);
  } else {
    bern_mut = NULL;
  }
  if (o.bern_cross) {
    bern_cross = new std::bernoulli_distribution(*o.bern_cross);
  } else {
    bern_cross = NULL;
  }
}

ArgStore::ArgStore(ArgStore&& o) {
  pop_size = o.pop_size;
  breed_pop_size = o.breed_pop_size;
  num_gens = o.num_gens;
  num_crossovers = o.num_crossovers;
  init_coup_var = o.init_coup_var;
  init_coup_mean = o.init_coup_mean;
  crossover_prob = o.crossover_prob;
  mutate_prob = o.mutate_prob;
  hypermutation_threshold = o.hypermutation_threshold;
  replacement_fraction = o.replacement_fraction;
  //managed pointers dumb copying
  long_bin = o.long_bin;
  short_bin = o.short_bin;
  bern_mut = o.bern_mut;
  bern_cross = o.bern_cross;
  o.long_bin = NULL;
  o.short_bin = NULL;
  o.bern_mut = NULL;
  o.bern_cross = NULL;
}

ArgStore::~ArgStore() {
  if (long_bin) {
    delete long_bin;
  }
  if (short_bin) {
    delete short_bin;
  }
  if (bern_mut) {
    delete bern_mut;
  }
  if (bern_cross) {
    delete bern_cross;
  }
}

void ArgStore::initialize_from_file(const char* fname) {
  pop_size = DEF_POP_SIZE;
  breed_pop_size = DEF_BREED_POP_SIZE;
  num_gens = DEF_NUM_GENS;
  num_crossovers = DEF_NUM_CROSSOVERS;
  init_coup_var = DEF_COUP_VAR;
  init_coup_mean = DEF_COUP_MEAN;
  crossover_prob = DEF_CROSSOVER_PROB;
  mutate_prob = DEF_MUTATE_PROB;
  out_fname = "output.csv";
  flags = 0;
  seed = 0;

  FILE* fp = fopen(fname, "r");
  if (!fp) {
    error(CODE_ARG_INVALID, "Could not open configuration file %s for reading.", fname);
  }
/*  char str[BUF_SIZE];
  double val;
  char val_str[BUF_SIZE];
  long pos = ftell(fp);
  int retval = fscanf(fp, "%s %lf\n", str, &val);
  int retval = read(fp, BUF_SIZE);*/
  char *str = NULL, *token = NULL, *val = NULL, *saveptr;
  size_t strlen;
  int retval = getline(&str, &strlen, fp);
  while ( retval > 0) {
    token = strtok_r(str, ":=", &saveptr);
    val = strtok_r(NULL, ":=", &saveptr);
    if (strcmp(token, "population_size") == 0) {
      pop_size = atoi(val);
    } else if (strcmp(token, "tournament_size") == 0 || strcmp(str, "breed_pop_size:") == 0){
      breed_pop_size = atoi(val);
    } else if (strcmp(token, "num_generations") == 0) {
      num_gens = atoi(val);
    } else if (strcmp(token, "num_crossovers") == 0) {
      num_crossovers = atoi(val);
    } else if (strcmp(token, "parameter_variance") == 0) {
      init_coup_var = atof(val);
    } else if (strcmp(token, "parameter_mean") == 0) {
      init_coup_mean = atof(val);
    } else if (strcmp(token, "crossover_probability") == 0) {
      crossover_prob = atof(val);
    } else if (strcmp(token, "mutation_probability") == 0) {
      mutate_prob = atof(val);
    } else if (strcmp(token, "hypermutation_threshold") == 0) {
      hypermutation_threshold = atof(val);
    } else if (strcmp(token, "replacement_fraction") == 0) {
      replacement_fraction = atof(val);
    } else if (strcmp(token, "output_file") == 0) {
      out_fname = val;
    } else if (strcmp(token, "seed") == 0) {
      seed = atoi(val);
    } else if (strcmp(token, "verbose") == 0) {
      flags = flags | VERBOSE;
    } else if (strcmp(token, "wait") == 0) {
      flags = flags | WAIT_CON;
    }
    if (str) {
      free(str);
      str = NULL;
    } else {
      break;
    }
    retval = getline(&str, &strlen, fp);
  }
  if (str) {
    free(str);
  }
  if (seed == 0) {
    std::random_device rd;
    seed = rd();
  }
  generator.seed(seed);
  if (long_bin) { delete long_bin; }
  if (bern_mut) { delete bern_mut; }
  if (bern_cross) { delete bern_cross; }
  long_bin = new std::binomial_distribution<unsigned char>(sizeof(unsigned long), mutate_prob);
  bern_mut = new std::bernoulli_distribution(mutate_prob);
  bern_cross = new std::bernoulli_distribution(crossover_prob);
  fclose(fp);
}

void ArgStore::initialize() {
  pop_size = DEF_POP_SIZE;
  breed_pop_size = DEF_BREED_POP_SIZE;
  num_gens = DEF_NUM_GENS;
  num_crossovers = DEF_NUM_CROSSOVERS;
  init_coup_var = DEF_COUP_VAR;
  init_coup_mean = DEF_COUP_MEAN;
  crossover_prob = DEF_CROSSOVER_PROB;
  mutate_prob = DEF_MUTATE_PROB;
  flags = 0;
  out_fname = "output.csv";
  long_bin = new std::binomial_distribution<unsigned char>(sizeof(unsigned long)*8, mutate_prob);
  bern_mut = new std::bernoulli_distribution(mutate_prob);
  bern_cross = new std::bernoulli_distribution(crossover_prob);
}

void ArgStore::initialize_from_args(size_t argc, char** argv) {
  pop_size = DEF_POP_SIZE;
  breed_pop_size = DEF_BREED_POP_SIZE;
  num_gens = DEF_NUM_GENS;
  num_crossovers = DEF_NUM_CROSSOVERS;
  init_coup_var = DEF_COUP_VAR;
  init_coup_mean = DEF_COUP_MEAN;
  crossover_prob = DEF_CROSSOVER_PROB;
  mutate_prob = DEF_MUTATE_PROB;
  flags = 0;
  out_fname = "output.csv";

  for (size_t i = 1; i < argc; ++i) {
    if (argv[i] != NULL && argv[i][0] == '-') {
      if (i == argc - 1) {
        if (argv[i][1] != 'w' && argv[i][1] != 'v') {
          std::cout << "ERROR: specify a parameter to use with " << argv[i] << std::endl;
          exit(0);
        }
      }

      bool test_again = true;
      char* namestr = strtok(argv[i], "=");
      char* valstr = strtok(NULL, "=");
      if (namestr && !valstr && argv[i][2] == '-') {
        std::cout << "ERROR: correct syntax is ./opt " << namestr << "=<value>" << std::endl;
        exit(0);
      } else if (namestr && valstr) {
        if (strcmp(namestr, "--pop_size") == 0) {
          pop_size = atoi( valstr );
          test_again = false;
        } else if (strcmp(namestr, "--survivors") == 0) {
          breed_pop_size = atoi( valstr );
          test_again = false;
        } else if (strcmp(namestr, "--generations") == 0) {
          num_gens = atoi( valstr );
          test_again = false;
        } else if (strcmp(namestr, "--variance") == 0 || strcmp(namestr, "--var") == 0) {
          init_coup_var = atoi( valstr );
          test_again = false;
        } else if (strcmp(namestr, "--crossover-prob") == 0) {
          crossover_prob = atof( valstr );
          test_again = false;
        } else if (strcmp(namestr, "--mutate-prob") == 0) {
          mutate_prob = atof( valstr );
          test_again = false;
        } else if (strcmp(namestr, "--hyper-mutate-threshold") == 0) {
          hypermutation_threshold = atof( valstr );
          test_again = false;
        } else if (strcmp(namestr, "--replacement-fraction") == 0) {
          replacement_fraction = atof( valstr );
          test_again = false;
        }
      }

      if (test_again) {
        if (i == argc - 1) {
          if (argv[i][1] != 'w' && argv[i][1] != 'v') {
            std::cout << "ERROR: specify a parameter to use with " << argv[i] << std::endl;
            exit(0);
          }
        }
        switch (argv[i][1]) {
          case 'p': pop_size = atoi( argv[i+1] );
              i++; break;
          case 'u': breed_pop_size = atoi( argv[i+1] );
              i++; break;
          case 'g': num_gens = atoi( argv[i+1] );
              i++; break;
          case 'c': num_crossovers = atoi( argv[i+1] );
              i++; break;
          case 'a': init_coup_var = atof( argv[i+1] );
              i++; break;
          case 'm': init_coup_mean = atof( argv[i+1] );
              i++; break;
          case 'r': crossover_prob = atof( argv[i+1] );
              i++; break;
          case 't': mutate_prob = atof( argv[i+1] );
              i++; break;
          case 'o': out_fname = argv[i+1];
              i++; break;
          case 'w': flags = flags | WAIT_CON;
              break;
          case 'v': flags = flags | VERBOSE;
              break;
          case 's': seed = atoi(argv[i+1]); generator.seed(seed);
              activate = false; i++; break;
        }
      }
    }
  }

  //entropy seed the random generator if no seed has been provided
  if (activate) {
    std::random_device rd;
    seed = rd();
    generator.seed(seed);
  }

  //print out warnings
  if (crossover_prob > 1.0 || crossover_prob < 0.0) {
    std::cout << "invalid probability of mutation specified by the -t or --crossover-prob parameter"<< std::endl;
    error(CODE_MATH_ERROR, "Enter a valid probability in the range [0, 1]");
  }
  if (mutate_prob > 1.0 || mutate_prob < 0.0) {
    std::cout << "invalid probability of mutation specified by the -t or --mutate-prob parameter"<< std::endl;
    error(CODE_MATH_ERROR, "Enter a valid probability in the range [0, 1]");
  }
  if (pop_size < 2) {
    error(CODE_ARG_INVALID, "The population size must be greater than 1.");
  }
  if (breed_pop_size < 2) {
    error(CODE_ARG_INVALID, "The number of survivors in each generation must be greater than 1.");
  }
  if (breed_pop_size > pop_size) {
    error(CODE_ARG_INVALID, "The number of surviving individuals cannot be larger than the size of the population.");
  }
  if (long_bin) { delete long_bin; }
  if (bern_mut) { delete bern_mut; }
  if (bern_cross) { delete bern_cross; }
  long_bin = new std::binomial_distribution<unsigned char>(sizeof(unsigned long)*8, mutate_prob);
  bern_mut = new std::bernoulli_distribution(mutate_prob);
  bern_cross = new std::bernoulli_distribution(crossover_prob);
}

void ArgStore::print_data() {
  //print out data
  std::cout << "Now optimizing using the following parameters:" << std::endl
	    << "\ttotal pop. size: " << get_pop_size() << std::endl
	    << "\tnumber survivors/gen: " << get_survivors() << std::endl
	    << "\tnumber gens: " << get_num_gens() << std::endl
	    << "\tvar: " << get_init_coup_var() << std::endl
	    << "\tmean: " << get_init_coup_mean() << std::endl
	    << "\tnumber of crossovers: " << get_num_crossovers() << std::endl
	    << "\tcrossover probability: " << get_crossover_prob() << std::endl
	    << "\tmutation probability: " << get_mutate_prob() << std::endl
	    << "\thyper-mutation threshold: " << get_hypermutation_threshold() << std::endl
	    << "\toutput file: " << get_out_fname() << std::endl;
	    if (activate) { std::cout << "\tusing entropy seed " << seed << std::endl; }
	    else { std::cout << "\tusing user provided seed " << seed << std::endl; }
}

unsigned int ArgStore::sample_binomial(unsigned int n) {
  if (n == sizeof(unsigned long)*CHAR_BIT) {
    if (!long_bin) {
      long_bin = new std::binomial_distribution<unsigned char>(sizeof(unsigned long)*CHAR_BIT);
    }
    return (*long_bin)(generator);
  } else if (n != short_bin_n) {
    if (short_bin != NULL) {
      delete short_bin;
    }
    short_bin_n = n;
    short_bin = new std::binomial_distribution<unsigned char>(n, mutate_prob);
  }
  return (*short_bin)(generator);
}

bool ArgStore::random_crossover() {
  if (!bern_cross) {
    bern_cross = new std::bernoulli_distribution(crossover_prob);
  }
  return (*bern_cross)(generator);
}

bool ArgStore::random_mutation() {
  if (!bern_mut) {
    bern_mut = new std::bernoulli_distribution(mutate_prob);
  }
  return (*bern_mut)(generator);
}

}
