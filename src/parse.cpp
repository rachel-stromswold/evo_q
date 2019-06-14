#include "parse.h"

namespace Genetics {

ArgStore::ArgStore() : generator() {
  pop_size = DEF_POP_SIZE;
  breed_pop_size = DEF_BREED_POP_SIZE;
  num_gens = DEF_NUM_GENS;
  num_crossovers = DEF_NUM_CROSSOVERS;
  init_coup_var = DEF_COUP_VAR;
  init_coup_mean = DEF_COUP_MEAN;
  mutate_prob = DEF_MUTATE_PROB;
  hypermutation_threshold = DEF_HYPER_THRESH;
  long_bin = NULL;
  short_bin = NULL;
  bern = NULL;
}

void ArgStore::initialize_from_file(const char* fname) {
  pop_size = DEF_POP_SIZE;
  breed_pop_size = DEF_BREED_POP_SIZE;
  num_gens = DEF_NUM_GENS;
  num_crossovers = DEF_NUM_CROSSOVERS;
  init_coup_var = DEF_COUP_VAR;
  init_coup_mean = DEF_COUP_MEAN;
  mutate_prob = DEF_MUTATE_PROB;
  out_fname = "output.csv";
  flags = 0;
  seed = 0;

  FILE* fp = fopen(fname, "r");
  if (!fp) {
    error(1, "Could not open configuration file %s for reading.", fname);
  }
  char str[BUF_SIZE];
  double val;
  char val_str[BUF_SIZE];
  long pos = ftell(fp);
  int retval = fscanf(fp, "%s %lf\n", str, &val);
  while (retval != EOF && retval != 0) {
    if (strcmp(str, "population_size:") == 0) {
      pop_size = (_uint)val;
    } else if (strcmp(str, "tournament_size:") == 0 || strcmp(str, "breed_pop_size:") == 0){
      breed_pop_size = (_uint)val;
    } else if (strcmp(str, "num_generations:") == 0) {
      num_gens = (_uint)val;
    } else if (strcmp(str, "num_crossovers:") == 0) {
      num_crossovers = (_uint)val;
    } else if (strcmp(str, "parameter_variance:") == 0) {
      init_coup_var = val;
    } else if (strcmp(str, "parameter_mean:") == 0) {
      init_coup_mean = val;
    } else if (strcmp(str, "mutatation_probability:") == 0) {
      mutate_prob = val;
    } else if (strcmp(str, "hypermutation_threshold:") == 0) {
      hypermutation_threshold = val;
    } else if (strcmp(str, "output_file:") == 0) {
      fseek(fp, pos, SEEK_SET);
      retval = fscanf(fp, "%s %s\n", str, val_str);
      out_fname = val_str;
    } else if (strcmp(str, "seed:") == 0) {
      seed = (int)val;
    } else if (strcmp(str, "verbose") == 0 || (strcmp(str, "verbose:") == 0 && val != 0.0)) {
      flags = flags | VERBOSE;
    } else if (strcmp(str, "wait") == 0 || (strcmp(str, "wait:") == 0 && val != 0.0)) {
      flags = flags | WAIT_CON;
    }
    pos = ftell(fp);
    retval = fscanf(fp, "%s %lf\n", str, &val);
  }
  if (seed == 0) {
    std::random_device rd;
    seed = rd();
  }
  generator.seed(seed);
  print_data();
  if (long_bin) { delete long_bin; }
  if (bern) { delete bern; }
  long_bin = new std::binomial_distribution<unsigned char>(sizeof(unsigned long), mutate_prob);
  bern = new std::bernoulli_distribution(mutate_prob);
}

void ArgStore::initialize() {
  pop_size = DEF_POP_SIZE;
  breed_pop_size = DEF_BREED_POP_SIZE;
  num_gens = DEF_NUM_GENS;
  num_crossovers = DEF_NUM_CROSSOVERS;
  init_coup_var = DEF_COUP_VAR;
  init_coup_mean = DEF_COUP_MEAN;
  mutate_prob = DEF_MUTATE_PROB;
  flags = 0;
  out_fname = "output.csv";
  long_bin = new std::binomial_distribution<unsigned char>(sizeof(unsigned long)*8, mutate_prob);
  bern = new std::bernoulli_distribution(mutate_prob);
}

void ArgStore::initialize_from_args(size_t argc, char** argv) {
  pop_size = DEF_POP_SIZE;
  breed_pop_size = DEF_BREED_POP_SIZE;
  num_gens = DEF_NUM_GENS;
  num_crossovers = DEF_NUM_CROSSOVERS;
  init_coup_var = DEF_COUP_VAR;
  init_coup_mean = DEF_COUP_MEAN;
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
        } else if (strcmp(namestr, "--mutate-prob") == 0) {
          mutate_prob = atof( valstr );
          test_again = false;
        } else if (strcmp(namestr, "--hyper-mutate-threshold") == 0) {
          hypermutation_threshold = atof( valstr );
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
  if (mutate_prob > 1.0 || mutate_prob < 0.0) {
    std::cout << "invalid probability of mutation specified by the -t or --mutate-prob parameter"<< std::endl;
    error(1, "Enter a valid probability in the range [0, 1]");
  }
  if (pop_size < 2) {
    error(1, "The population size must be greater than 1.");
  }
  if (breed_pop_size < 2) {
    error(1, "The number of survivors in each generation must be greater than 1.");
  }
  if (breed_pop_size > pop_size) {
    error(1, "The number of surviving individuals cannot be larger than the size of the population.");
  }
  print_data();
  if (long_bin) { delete long_bin; }
  if (bern) { delete bern; }
  long_bin = new std::binomial_distribution<unsigned char>(sizeof(unsigned long)*8, mutate_prob);
  bern = new std::bernoulli_distribution(mutate_prob);
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
	    << "\tmutation probability: " << get_mutate_prob() << std::endl
	    << "\thyper-mutation threshold: " << get_hypermutation_threshold() << std::endl
	    << "\toutput file: " << get_out_fname() << std::endl;
	    if (activate) { std::cout << "\tusing entropy seed " << seed << std::endl; }
	    else { std::cout << "\tusing user provided seed " << seed << std::endl; }
}

unsigned int ArgStore::sample_binomial(unsigned int n) {
  if (n == sizeof(unsigned long)*8) {
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

bool ArgStore::sample_bernoulli() {
  return (*bern)(generator);
}

}
