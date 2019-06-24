#include "../include/evo_q.hpp"
#define CATCH_CONFIG_MAIN
#include "../include/catch.hpp"
#define N_TRIALS 100

#define NUM_BITS	16
#define NUM_OBJS  	2
#define NUM_CHROMS	1

#define M_RANGE		2.0

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

class TestProblemSingle : public Genetics::Problem {
public:
  TestProblemSingle() : Genetics::Problem(NUM_BITS, NUM_CHROMS, 1) {
    map.initialize(1, Genetics::t_real);
  }
  void evaluate_fitness(Genetics::Organism* org) {
    double x = org->read_real(0);
    org->set_fitness(0, x*x);
  }
};

class TestProblemMulti : public Genetics::Problem {
public:
  TestProblemMulti() : Genetics::Problem(NUM_BITS, NUM_CHROMS, NUM_OBJS) {
    map.initialize(NUM_CHROMS, Genetics::t_real);
  }
  void evaluate_fitness(Genetics::Organism* org) {
    double x = org->read_real(0);
    org->set_fitness(0, x*x);
    org->set_fitness(1, (x - 2)*(x - 2));
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
  Genetics::PhenotypeMap map_16(16);
  Genetics::PhenotypeMap map_32(32);
  Genetics::Organism template_16_bits_1_params(16, 1, &map_16);
  Genetics::Organism template_32_bits_2_params(32, 2, &map_32);
  LCG generator;

  SECTION ( "The organism correctly sets and reads unsigned integers" ) {
    map_16.initialize(1, Genetics::t_uint);
    map_32.initialize(2, Genetics::t_uint);
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
    map_16.initialize(1, Genetics::t_real);
    map_32.initialize(2, Genetics::t_real);
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

TEST_CASE ("Populations are correctly created and data is successfully read", "[populations]") {
//  Genetics::ArgStore inst;
  Genetics::String conf_file("ga.conf");
//  inst.initialize_from_file(conf_file.c_str());
  Genetics::PhenotypeMap map(NUM_BITS);
  map.initialize(NUM_CHROMS, Genetics::t_real);
  map.set_range(0, -10, 10);
  Genetics::Organism tmplt(NUM_BITS, NUM_OBJS, &map);
  tmplt.set_real(0, 4.6);
  Genetics::Population pop_none(NUM_BITS, NUM_OBJS, &tmplt, &map, conf_file);
  Genetics::Population pop_tmplt(NUM_BITS, NUM_OBJS, &map, conf_file);

  SECTION ( "The header is read correctly" ) {
    Genetics::Vector<Genetics::String> pop_dat = pop_none.get_header();
    for (Genetics::Vector<Genetics::String>::iterator it = pop_dat.begin(); it != pop_dat.end(); ++it) {
      printf("%s ", it->c_str());
    }
    printf("\n");
    Genetics::Vector<Genetics::String> pop_dat2 = pop_tmplt.get_header();
    for (Genetics::Vector<Genetics::String>::iterator it = pop_dat2.begin(); it != pop_dat2.end(); ++it) {
      printf("%s ", it->c_str());
    }
    printf("\n");
  }
  SECTION ( "The population data is read correctly" ) {
    Genetics::Vector<Genetics::String> pop_dat = pop_none.get_pop_data();
    for (Genetics::Vector<Genetics::String>::iterator it = pop_dat.begin(); it != pop_dat.end(); ++it) {
      printf("%s ", it->c_str());
    }
    printf("\n");
    Genetics::Vector<Genetics::String> pop_dat2 = pop_tmplt.get_pop_data();
    for (Genetics::Vector<Genetics::String>::iterator it = pop_dat2.begin(); it != pop_dat2.end(); ++it) {
      printf("%s ", it->c_str());
    }
    printf("\n");
  }
}

TEST_CASE ("Simple evolution of a single objective converges to roughly appropriate result", "[populations]") {
  TestProblemSingle prob;
  Genetics::String conf_file("ga.conf");
  Genetics::Population pop( NUM_BITS, NUM_OBJS, &(prob.map), conf_file);
  for (size_t i = 0; i < 100; ++i) {
    pop.evaluate(&prob);
    pop.iterate();
  }
  pop.evaluate(&prob);
  REQUIRE( APPROX(pop.get_best_organism()->get_fitness(0), 0) );
}

TEST_CASE ("Simple evolution of a multi objective converges to roughly appropriate result", "[populations]") {
  TestProblemMulti prob;
  Genetics::String conf_file("ga.conf");
  Genetics::Population pop( NUM_BITS, NUM_OBJS, &(prob.map), conf_file);
  for (size_t i = 0; i < 100; ++i) {
    pop.evaluate(&prob);
    pop.iterate();
  }
  pop.evaluate(&prob);
  Genetics::Vector<Genetics::String> pop_dat = pop.get_pop_data();
  for (size_t i = 0; i < pop_dat.size(); ++i) {
    std::cout << "organism " << i << " = " << pop_dat[i] << std::endl;
  }
}
