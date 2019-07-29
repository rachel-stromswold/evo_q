//#include "../include/evo_q.hpp"
#include "../include/convergence.h"
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#define CATCH_CONFIG_MAIN
#include "../include/catch.hpp"
#define N_TRIALS 100
#include <chrono>
#include <thread>

#define NUM_BITS	16
#define NUM_OBJS  	2
#define NUM_CHROMS	1

#define M_RANGE		2.0
#define TEST_CONV_GEN	10
#define SLEEP_TIME	1000

#define APPROX(a,b) ((a-b < 0.01 && a-b >= 0) || (b-a < 0.01 && b-a >= 0))

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

//for random number generation
class LCG {
private:
  //a-1 is divisible by 4 (and obviously 2 (the only prime factor of 2^64)
  static const unsigned long a = 3*( ((unsigned long)1 << 24) + 5 );
  static const unsigned long c = 170859375;//=15^7 which is relatively prime to m = 2^64
  static const unsigned long x0 = 0xA1A3A5A7A9ABAFA5;
  static const unsigned long high_mask = ULONG_MAX << 32;

  unsigned long state; void update_state() { state = (state*a + c)/* the modulo 64 is implicit */;}
public:
  LCG() : state(x0) {}
  unsigned long random_ulong() {
    unsigned long ret = state >> 32;
    update_state();
    ret = ret | (state & high_mask);
    return ret;
  }
  double random_int() {
    unsigned long val = random_ulong();
    if ( val & (~(ULONG_MAX >> 1)) ) {
      return -1*((int)val);
    } else {
      return (int)val;
    }
  }
  double random_real() { return (double)random_ulong()/ULONG_MAX; }
};

class ChromosomeTestFixture {
protected:
  Genetics::ArgStore args;
  Genetics::Chromosome chrom_96_bits_A, chrom_96_bits_B;//this should be large enough to require the use of two words
  Genetics::Chromosome chrom_32_bits_A, chrom_32_bits_B;
  Genetics::PhenotypeMap map_96;
  Genetics::PhenotypeMap map_32;
public:
  ChromosomeTestFixture() : chrom_96_bits_A(96), chrom_96_bits_B(96), chrom_32_bits_A(32), chrom_32_bits_B(32), map_96(96), map_32(32) {}
};

class TestProblemPenalties : public Genetics::Problem {
private:
  bool apply_penalty;
  double penalty_amt;
public:
  TestProblemPenalties(double p_penalty_amt) : Genetics::Problem(NUM_BITS, NUM_CHROMS, 1) {
    map->initialize(1, Genetics::t_real);
    apply_penalty = true;
    penalty_amt = p_penalty_amt;
  }
  void evaluate_fitness(Genetics::Organism* org) {
    org->set_fitness(0, 1);
    if (apply_penalty) {
      org->apply_penalty(penalty_amt);
      apply_penalty = false;
    }
  }
};

class TestProblemSingle : public Genetics::Problem {
public:
  TestProblemSingle() : Genetics::Problem(NUM_BITS, NUM_CHROMS, 1) {
    map->initialize(1, Genetics::t_real);
  }
  void evaluate_fitness(Genetics::Organism* org) {
    double x = org->read_real(0);
    org->set_fitness(0, -(x*x));
  }
};

class TestProblemMulti : public Genetics::Problem {
public:
  TestProblemMulti() : Genetics::Problem(NUM_BITS, NUM_CHROMS, NUM_OBJS) {
    map->initialize(NUM_CHROMS, Genetics::t_real);
  }
  void evaluate_fitness(Genetics::Organism* org) {
    double x = org->read_real(0);
    org->set_fitness(0, -x*x);
    org->set_fitness(1, -(x - 2)*(x - 2));
    org->apply_penalty(0.0);
  }
};

class TestProblemSlow : public Genetics::Problem {
public:
  TestProblemSlow() : Genetics::Problem(NUM_BITS, NUM_CHROMS, 1) {
    map->initialize(1, Genetics::t_real);
  }
  void evaluate_fitness(Genetics::Organism* org) {
    std::this_thread::sleep_for( std::chrono::milliseconds(SLEEP_TIME) );
    double x = org->read_real(0);
    org->set_fitness(0, -(x*x));
  }
};

// =============================== CHROMOSOME TEST CASES ===============================

TEST_CASE_METHOD( ChromosomeTestFixture, "Ensure that chromosomes are correctly encoded and decoded for integers", "[chromosomes]" ) {
  REQUIRE(chrom_96_bits_A.get_n_bits() == 96);
  REQUIRE(chrom_96_bits_B.get_n_bits() == 96);
  REQUIRE(chrom_32_bits_A.get_n_bits() == 32);
  REQUIRE(chrom_32_bits_B.get_n_bits() == 32);
  map_96.initialize(2, Genetics::t_uint);
  map_32.initialize(2, Genetics::t_uint);
  
  unsigned long mask_96 = ( (unsigned long)1 << (96/2) ) - 1;
  unsigned long mask_32 = ( (unsigned long)1 << (32/2) ) - 1;

  LCG generator;
  unsigned long val = 1;
  for (int i = 0; i < N_TRIALS; ++i) {
    val = generator.random_ulong() & mask_32;
    chrom_32_bits_A.set_to_ulong(&map_32, 0, val);
    INFO("Setting 32 bit genome: i = " << i)
    REQUIRE(chrom_32_bits_A.gene_to_ulong(&map_32, 0) == val);
  }
  for (int i = 0; i < N_TRIALS; ++i) {
    val = generator.random_ulong() & mask_32;
    chrom_32_bits_A.set_to_ulong(&map_32, 1, val);
    INFO("Setting 32 bit genome: i = " << i)
    REQUIRE(chrom_32_bits_A.gene_to_ulong(&map_32, 1) == val);
  }
  for (int i = 0; i < N_TRIALS; ++i) {
    val = generator.random_ulong() & mask_96;
    chrom_96_bits_A.set_to_ulong(&map_96, 0, val);
    INFO("Setting 96 bit genome: i = " << i)
    REQUIRE(chrom_96_bits_A.gene_to_ulong(&map_96, 0) == val);
  }
  for (int i = 0; i < N_TRIALS; ++i) {
    val = generator.random_ulong() & mask_96;
    chrom_96_bits_A.set_to_ulong(&map_96, 1, val);
    INFO("Setting 96 bit genome: i = " << i)
    REQUIRE(chrom_96_bits_A.gene_to_ulong(&map_96, 1) == val);
  }
}

TEST_CASE_METHOD( ChromosomeTestFixture, "Ensure that chromosomes are correctly encoded and decoded for reals", "[chromosomes]" ) {
  REQUIRE(chrom_96_bits_A.get_n_bits() == 96);
  REQUIRE(chrom_96_bits_B.get_n_bits() == 96);
  REQUIRE(chrom_32_bits_A.get_n_bits() == 32);
  REQUIRE(chrom_32_bits_B.get_n_bits() == 32);
  //initialize the 96 bit map
  map_96.initialize(2, Genetics::t_real);
  map_96.set_range(0, 0, 100);
  map_96.set_range(1, 0, 1);
  //initialize the 32 bit map
  map_32.initialize(2, Genetics::t_real);
  map_32.set_range(0, 0, 100);
  map_32.set_range(1, 0, 1);

  LCG generator;
  double val = 1.0;
  for (int i = 0; i < N_TRIALS; ++i) {
    val = generator.random_real()*100.0;
    chrom_32_bits_A.set_to_num(&map_32, 0, val);
    INFO("Setting 32 bit genome: i = " << i << " val = " << val << " rep = " << chrom_32_bits_A.gene_to_num(&map_32, 0))
    REQUIRE(APPROX(chrom_32_bits_A.gene_to_num(&map_32, 0), val));
  }
  for (int i = 0; i < N_TRIALS; ++i) {
    val = generator.random_real();
    chrom_32_bits_A.set_to_num(&map_32, 1, val);
    INFO("Setting 32 bit genome: i = " << i << " val = " << val << " rep = " << chrom_32_bits_A.gene_to_num(&map_32, 1))
    REQUIRE(APPROX(chrom_32_bits_A.gene_to_num(&map_32, 1), val));
  }
  for (int i = 0; i < N_TRIALS; ++i) {
    val = generator.random_real()*100.0;
    chrom_96_bits_A.set_to_num(&map_96, 0, val);
    INFO("Setting 32 bit genome: i = " << i << " val = " << val << " rep = " << chrom_32_bits_A.gene_to_num(&map_32, 0))
    REQUIRE(APPROX(chrom_96_bits_A.gene_to_num(&map_96, 0), val));
  }
  for (int i = 0; i < N_TRIALS; ++i) {
    val = generator.random_real();
    chrom_96_bits_A.set_to_num(&map_96, 1, val);
    INFO("Setting 32 bit genome: i = " << i << " val = " << val << " rep = " << chrom_32_bits_A.gene_to_num(&map_32, 1))
    REQUIRE(APPROX(chrom_96_bits_A.gene_to_num(&map_96, 1), val));
  }
}

TEST_CASE_METHOD( ChromosomeTestFixture, "Ensure that chromosome exchange works properly", "[chromosomes]" ) {
  REQUIRE(chrom_96_bits_A.get_n_bits() == 96);
  REQUIRE(chrom_96_bits_B.get_n_bits() == 96);
  REQUIRE(chrom_32_bits_A.get_n_bits() == 32);
  REQUIRE(chrom_32_bits_B.get_n_bits() == 32);
  map_96.initialize(1, Genetics::t_uint);
  map_32.initialize(1, Genetics::t_uint); 

  for (int i = 0; i < 64; ++i) {
    chrom_96_bits_A.set_to_ulong(&map_96, 0, ULONG_MAX);
    REQUIRE( chrom_96_bits_A.gene_to_ulong(&map_96, 0) == ULONG_MAX );
    chrom_96_bits_B.set_to_ulong(&map_96, 0, 0);
    //perform the exchange
    chrom_96_bits_A.exchange(&chrom_96_bits_B, i);
    unsigned long tmp_long = chrom_96_bits_A.gene_to_ulong(&map_96, 0);
    INFO("i = " << i)
    REQUIRE(tmp_long == (ULONG_MAX << i));
  }
  unsigned long mask = ((unsigned long)1 << 32) - 1;
  for (int i = 0; i < 32; ++i) {
    chrom_32_bits_A.set_to_ulong(&map_32, 0, ULONG_MAX);
    REQUIRE( chrom_32_bits_A.gene_to_ulong(&map_32, 0) == (ULONG_MAX & mask) );
    chrom_32_bits_B.set_to_ulong(&map_32, 0, 0);
    //perform the exchange
    chrom_32_bits_A.exchange(&chrom_32_bits_B, i);

    unsigned long tmp_long = chrom_32_bits_A.gene_to_ulong(&map_32, 0);
    INFO("i = " << i)
    REQUIRE( tmp_long == ((ULONG_MAX << i) & mask) );
  }
}

TEST_CASE( "Ensure that phenotype mappings are set up properly", "[PhenotypeMapping]") {
  SECTION ( "Phenotype maps correctly handle 1 type initialization" ) {
    Genetics::PhenotypeMap map_int(32);
    map_int.initialize(4, Genetics::t_int);
    REQUIRE(map_int.get_block_length(0) == 8);
    REQUIRE(map_int.get_block_length(1) == 8);
    REQUIRE(map_int.get_block_length(2) == 8);
    REQUIRE(map_int.get_block_length(3) == 8);
    REQUIRE(map_int.is_int(0));
    REQUIRE(map_int.is_int(1));
    REQUIRE(map_int.is_int(2));
    REQUIRE(map_int.is_int(3));

    Genetics::PhenotypeMap map_real(32);
    map_real.initialize(4, Genetics::t_real);
    REQUIRE(map_real.is_real(0));
    REQUIRE(map_real.is_real(1));
    REQUIRE(map_real.is_real(2));
    REQUIRE(map_real.is_real(3));
    REQUIRE(map_real.get_block_length(0) == 8);
    REQUIRE(map_real.get_block_length(1) == 8);
    REQUIRE(map_real.get_block_length(2) == 8);
    REQUIRE(map_real.get_block_length(3) == 8);
  }

  SECTION ( "Phenotype maps correctly handle more complex initialization" ) {
    Genetics::PhenotypeMap map_mixed(32);
    Genetics::Vector<Genetics::VarContainer> vc;
    vc.push_back(Genetics::VarContainer(0, 0.0, 1.0, Genetics::t_real));
    vc.push_back(Genetics::VarContainer(0, 0.0, 0.0, Genetics::t_int));
    vc.push_back(Genetics::VarContainer(16, 0.0, 0.0, Genetics::t_bitstream));
    map_mixed.initialize(vc);
    REQUIRE(map_mixed.is_real(0));
    REQUIRE(map_mixed.is_int(1));
    REQUIRE(map_mixed.is_bitstream(2));
    REQUIRE(map_mixed.get_block_length(0) == 10);
    REQUIRE(map_mixed.get_block_length(1) == 5);
    REQUIRE(map_mixed.get_block_length(2) == 16);
  }

  SECTION ( "Ensure reading from a list works" ) {
    Genetics::Vector<Genetics::String> param_list;
    param_list.push_back("real_(0, 1.0)");
    param_list.push_back("int");
    param_list.push_back("bitstream_16");
    Genetics::PhenotypeMap map(32);

    Genetics::Vector<Genetics::VarContainer> vc_list;
    for (size_t i = 0; i < param_list.size(); ++i) {
      Genetics::String type(param_list[i]);
      const char* type_str = type.c_str();
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
    map.initialize(vc_list);
#ifdef PRINT_DEBUG_DATA
    std::cout << "map_size = " << map.get_num_params() << "\n";
#endif
    REQUIRE(map.is_real(0));
    REQUIRE(map.is_int(1));
    REQUIRE(map.is_bitstream(2));
    REQUIRE(map.get_block_length(0) == 10);
    REQUIRE(map.get_block_length(1) == 5);
    REQUIRE(map.get_block_length(2) == 16);
  }
}

// =============================== POPULATION TEST CASES ===============================

TEST_CASE ( "Organisms are correctly created and decoded", "[organisms]" ) {
  std::shared_ptr<Genetics::PhenotypeMap> map_16 = std::make_shared<Genetics::PhenotypeMap>(16);
  std::shared_ptr<Genetics::PhenotypeMap> map_32 = std::make_shared<Genetics::PhenotypeMap>(32);
  Genetics::Organism template_16_bits_1_params(16, 1, map_16);
  Genetics::Organism template_32_bits_2_params(32, 2, map_32);
  LCG generator;

  SECTION ( "The organism correctly sets and reads unsigned integers" ) {
    map_16->initialize(1, Genetics::t_uint);
    map_32->initialize(2, Genetics::t_uint);
    for (int i = 0; i < N_TRIALS; ++i) {
      int val = i - N_TRIALS/2;
      template_16_bits_1_params.set_int(0, val);
      REQUIRE(template_16_bits_1_params.read_int(0) == val);
    }
    for (int i = 0; i < N_TRIALS; ++i) {
      int val1 = i - N_TRIALS;
      template_32_bits_2_params.set_int(0, val1);
      int val2 = i;
      template_32_bits_2_params.set_int(1, val2);
      REQUIRE(template_32_bits_2_params.read_int(0) == val1);
      REQUIRE(template_32_bits_2_params.read_int(1) == val2);
    }
  }
  SECTION ( "The organism correctly sets and read reals" ) {
    map_16->initialize(1, Genetics::t_real);
    map_32->initialize(2, Genetics::t_real);
    for (int i = 0; i < N_TRIALS; ++i) {
      double val = generator.random_real();
      template_16_bits_1_params.set_real(0, val);
      REQUIRE( APPROX(template_16_bits_1_params.read_real(0), val) );
    }
    for (int i = 0; i < N_TRIALS; ++i) {
      double val1 = generator.random_real();
      template_32_bits_2_params.set_real(0, val1);
      double val2 = generator.random_real();
      template_32_bits_2_params.set_real(1, val2);
      INFO("res = " << template_32_bits_2_params.read_real(0) << " val = " << val1)
      REQUIRE( APPROX(template_32_bits_2_params.read_real(0), val1) );
      INFO("res = " << template_32_bits_2_params.read_real(1) << " val = " << val2)
      REQUIRE( APPROX(template_32_bits_2_params.read_real(1), val2) );
    }
  }
}

TEST_CASE ("ArgStore successfully parses a file") {
  TestProblemSingle prob;
  Genetics::Conv_Plateau plat_cut(0.05, TEST_CONV_GEN / 2);
  Genetics::ArgStore args;
  args.initialize_from_file("ga.conf");
  Genetics::Population pop( NUM_BITS, 1, prob.map, args);

  //ga.conf must be as follows for this to work:
  /*
   * population_size: 50
   * tournament_size: 2
   * num_generations: 50
   * num_crossovers: 2
   * parameter_variance: 0.3
   * parameter_mean: 0.0
   * mutation_probability: 0.016
   * crossover_probability: 0.9
   * hypermutation_threshold: 1.0
   * output_file: output.csv
   */
  REQUIRE( pop.get_args().get_pop_size() == 50 );
  REQUIRE( pop.get_args().get_survivors() == 2 );
  REQUIRE( pop.get_args().get_num_gens() == 50 );
  REQUIRE( pop.get_args().get_num_crossovers() == 2 );
  REQUIRE( pop.get_args().get_init_coup_var() == 0.3 );
  REQUIRE( pop.get_args().get_crossover_prob() == 0.9 );
  REQUIRE( pop.get_args().get_mutate_prob() == 0.016 );
  REQUIRE( pop.get_args().get_hypermutation_threshold() == 1.0 );
  REQUIRE( pop.get_args().get_out_fname() == "output.csv" );
}

TEST_CASE ("Populations are correctly created and data is successfully read", "[populations]") {
//  Genetics::ArgStore inst;
  Genetics::ArgStore args;
  args.initialize_from_file("ga.conf");
//  inst.initialize_from_file(conf_file.c_str());
  std::shared_ptr<Genetics::PhenotypeMap> map = std::make_shared<Genetics::PhenotypeMap>(NUM_BITS);
  map->initialize(NUM_CHROMS, Genetics::t_real);
  map->set_range(0, -10, 10);
  Genetics::Organism tmplt(NUM_BITS, NUM_OBJS, map);
  tmplt.set_real(0, 4.6);
  Genetics::Population pop_none(NUM_BITS, NUM_OBJS, &tmplt, map, args);
  Genetics::Population pop_tmplt(NUM_BITS, NUM_OBJS, map, args);

  SECTION ( "The header is read correctly" ) {
    Genetics::Vector<Genetics::String> pop_dat = pop_none.get_header();
    for (unsigned i = 0; i < pop_dat.size() / 3; ++i) {
      REQUIRE( pop_dat[i*3] == "x_0");
      REQUIRE( pop_dat[i*3 + 1] == "f_0(x)" );
      REQUIRE( pop_dat[i*3 + 2] == "f_1(x)" );
    }
#ifdef PRINT_GARBAGE
    for (Genetics::Vector<Genetics::String>::iterator it = pop_dat.begin(); it != pop_dat.end(); ++it) {
      printf("%s ", it->c_str());
    }
    printf("\n");
    Genetics::Vector<Genetics::String> pop_dat2 = pop_tmplt.get_header();
    for (Genetics::Vector<Genetics::String>::iterator it = pop_dat2.begin(); it != pop_dat2.end(); ++it) {
      printf("%s ", it->c_str());
    }
    printf("\n");
#endif
  }
  SECTION ( "The population data is read correctly" ) {
    Genetics::Vector<Genetics::String> pop_dat = pop_none.get_pop_data();
    Genetics::Vector<Genetics::String> pop_dat2 = pop_tmplt.get_pop_data();
    REQUIRE( pop_dat.size() % 3 == 0 );
    REQUIRE( pop_dat.size() == pop_dat2.size() );

#ifdef PRINT_GARBAGE
    for (Genetics::Vector<Genetics::String>::iterator it = pop_dat.begin(); it != pop_dat.end(); ++it) {
      printf("%s ", it->c_str());
    }
    printf("\n");
    for (Genetics::Vector<Genetics::String>::iterator it = pop_dat2.begin(); it != pop_dat2.end(); ++it) {
      printf("%s ", it->c_str());
    }
    printf("\n");
#endif
  }
}

TEST_CASE ("Penalties are applied properly", "[populations]") {
  double penalty_str = 4.0;
  TestProblemPenalties prob(penalty_str);
  Genetics::ArgStore args;
  args.initialize_from_file("ga.conf");
  Genetics::Population pop( NUM_BITS, 1, prob.map, args);
  pop.evaluate(&prob);

  double fit_0 = pop.get_organism(0)->get_fitness(0);
  for (unsigned i = 1; i < pop.get_offspring_num(); ++i) {
    double fit_i = pop.get_organism(i)->get_fitness(0);
    INFO( "i = " << i )
    REQUIRE( fit_i == 1.0 );
    REQUIRE( fit_0 < fit_i );
  }
}

TEST_CASE ("Simple evolution of a single objective converges to roughly appropriate result", "[populations]") {
  TestProblemSingle prob;
  Genetics::ArgStore args;
  args.initialize_from_file("ga.conf");
  Genetics::Population pop( NUM_BITS, 1, prob.map, args);
  for (size_t i = 0; i < 100; ++i) {
    pop.evaluate(&prob);
    pop.iterate();
  }
  pop.evaluate(&prob);
  REQUIRE( APPROX(pop.get_best_organism()->get_fitness(0), 0) );
}

TEST_CASE ("Simple evolution of a multi objective converges to roughly appropriate result", "[populations]") {
  TestProblemMulti prob;
  Genetics::ArgStore args;
  args.initialize_from_file("ga.conf");
  Genetics::Population_NSGAII pop( NUM_BITS, NUM_OBJS, prob.map, args);
  bool converged = false;
  while (!converged) {
    pop.evaluate(&prob);
    converged = pop.iterate();
  }
  int modulo = NUM_OBJS + prob.map->get_num_params() + 1;
  pop.evaluate(&prob);
  Genetics::Vector<Genetics::String> pop_dat = pop.get_pop_data();
  REQUIRE( pop_dat.size() == pop.get_header().size() );
#ifdef PRINT_GARBAGE
  size_t num_params = prob.map->get_num_params();
  for (size_t i = 0; modulo*(i + 1) - 1 < pop_dat.size(); ++i) {
    std::cout << "organism " << i
	      << ", rank = " << pop_dat[modulo*i];
    for (size_t j = 0; j < num_params; ++j) {
      std::cout << ", x_" << j << " = " << pop_dat[modulo*i + j + 1];
    }
    for (size_t j = 0; j < prob.N_OBJS; ++j) {
      std::cout << ", fitness_" << j << " = " << pop_dat[modulo*i + j + num_params + 1];
    }
    std::cout << std::endl;
  }
#endif
}

TEST_CASE ("Convergence functions work", "[populations]") {
  Genetics::Conv_VarianceCutoff var_cut(0.011);
  Genetics::Conv_RangeCutoff range_cut(0.11);
  Genetics::Conv_Plateau plat_cut(0.1, TEST_CONV_GEN / 2);
  Genetics::FitnessStats s[2];
  //note that 1/n - 1/(n+1) = 1/(n^2 + n) < 1/10 => n^2 + n > 10 => n >= 3
  for (unsigned i = 1; i < TEST_CONV_GEN + 1; ++i) {
    s[0].var = 1.0 / (i*i);
    s[0].mean = 1.0 - 1.5 / i;
    s[0].max = 1.0 - 1.0 / i;
    s[0].min = 1.0 - 2.0 / i;
    s[1].var = 1.0 / (i*i*i*i);
    s[1].mean = 1.0 - 1.5 / (i*i);
    s[1].max = 1.0 - 1.0 / (i*i);
    s[1].min = 1.0 - 2.0 / (i*i);
    if (i < TEST_CONV_GEN) {
      INFO("i = " << i << " var_1 = " << s[0].var << " var_2 = " << s[1].var)
      REQUIRE( !var_cut.evaluate_convergence(2, s) );
      INFO("i = " << i << " range_1 = " << s[0].max - s[0].min << " range_2 = " << s[1].max - s[1].min)
      REQUIRE( !range_cut.evaluate_convergence(2, s) );
    } else {
      INFO("i = " << i << " var_1 = " << s[0].var << " var_2 = " << s[1].var)
      REQUIRE( var_cut.evaluate_convergence(2, s) );
      INFO("i = " << i << " range_1 = " << s[0].max - s[0].min << " range_2 = " << s[1].max - s[1].min)
      REQUIRE( range_cut.evaluate_convergence(2, s) );
    }
  }
}

TEST_CASE ("Combine convergence checking and evolution", "[populations]") {
  TestProblemSingle prob;
  Genetics::Conv_Plateau plat_cut(0.05, TEST_CONV_GEN / 2);
  Genetics::ArgStore args;
  args.initialize_from_file("ga.conf");
  Genetics::Population pop( NUM_BITS, 1, prob.map, args);

  bool converged = false;
  int generation = 0;
  while (!converged) {
    pop.evaluate(&prob);
    converged = pop.iterate(&plat_cut);
    generation++;
  }
  std::cout << "Converged after " << generation << " generations\n";
  pop.evaluate(&prob);
  double mean_x, mean_fit, var_x, var_fit;
  mean_x = 0.0;mean_fit = 0.0;var_x = 0.0;var_fit = 0.0;
  Genetics::Vector<Genetics::String> pop_dat = pop.get_pop_data();
  for (size_t i = 0; 2*i + 1 < pop_dat.size(); ++i) {
    double f = atof( pop_dat[2*i + 1].c_str() );
    double x = atof( pop_dat[2*i].c_str() );
    mean_fit += f / pop_dat.size();
    mean_x += x / pop_dat.size();
  }
  for (size_t i = 0; 2*i + 1 < pop_dat.size(); ++i) {
    double f = atof( pop_dat[2*i + 1].c_str() );
    double x = atof( pop_dat[2*i].c_str() );
    var_fit += (f - mean_fit)*(f - mean_fit)/ (pop_dat.size() + 1);
    var_x += (x - mean_x)*(x - mean_x) / (pop_dat.size() + 1);
#ifdef PRINT_GARBAGE
    std::cout << "organism " << i
	      << " x = " << pop_dat[2*i]
	      << " fitness = " << pop_dat[2*i + 1] << std::endl;
#endif
  }
  std::cout << "x1 = " << mean_x << "\u00B1" << var_x
	    << ", fitness = " << mean_fit << "\u00B1" << var_fit << std::endl;
}

TEST_CASE ("timing info with and without parallelization", "[populations]") {
  TestProblemSingle prob;
  Genetics::Conv_Plateau plat_cut(0.05, TEST_CONV_GEN / 2);
  Genetics::ArgStore args_async;
  Genetics::ArgStore args_serial;
  args_async.set_async(true);
  args_serial.set_async(false);

  double serial_mean = 0.0;double async_mean = 0.0;
  double serial_var = 0.0;double async_var = 0.0;
  double serial_times[N_TRIALS];
  double async_times[N_TRIALS];
  for (unsigned i = 0; i < N_TRIALS; ++i) {
    //evaluate serial
    auto t_start = std::chrono::high_resolution_clock::now();
      Genetics::Population pop_serial( NUM_BITS, 1, prob.map, args_serial);
      pop_serial.run(&prob);
    auto t_end = std::chrono::high_resolution_clock::now();
    serial_times[i] = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    serial_mean += serial_times[i] / N_TRIALS;

    //evaluate async
    t_start = std::chrono::high_resolution_clock::now();
      Genetics::Population pop_async( NUM_BITS, 1, prob.map, args_async);
      pop_async.run(&prob);
    t_end = std::chrono::high_resolution_clock::now();
    async_times[i] = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    async_mean += async_times[i] / N_TRIALS;
  }

  for (unsigned i = 0; i < N_TRIALS; ++i) {
    serial_var += (serial_times[i] - serial_mean)*(serial_times[i] - serial_mean) / (N_TRIALS - 1);
    async_var += (async_times[i] - async_mean)*(async_times[i] - async_mean) / (N_TRIALS - 1);
  }
  std::cout << "Benchmark\t\tnumber of trials\telapsed time\n";
  std::cout << std::fixed << std::setprecision(2) << "Serial evaluation\t"
	    << N_TRIALS << "\t\t\t" << serial_mean << "\u00B1" << sqrt(serial_var)
	    << "\n";
  std::cout << std::fixed << std::setprecision(2) << "Parallel evaluation\t"
	    << N_TRIALS << "\t\t\t" << async_mean << "\u00B1" << sqrt(async_var)
	    << "\n";
}
