//#include "../include/evo_q/evo_q.hpp"

#define CATCH_CONFIG_ENABLE_BENCHMARKING
#define CATCH_CONFIG_MAIN

#include "testProblems.hpp"
#include "../include/catch.hpp"

#define N_TRIALS 100

#define M_RANGE		2.0
#define TEST_CONV_GEN	10

#define EPSILON 0.01
#define APPROX(a,b) ((a-b < EPSILON && a-b >= 0) || (b-a < EPSILON && b-a >= 0))

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

#define INV_ERR_N 5
template <class PopType>
class PopulationPrinter {
private:
  std::ofstream ofs;
  PopType* pop;
  int gen = -1;
public:
  PopulationPrinter(PopType* p_pop, std::string fname) :
  ofs(fname, std::ofstream::out)
  {
    pop = p_pop;
  }

  void print_line(bool include_best = true) {
    if (gen < 0) {
      ofs << "gen";
      if (include_best) {
        //ofs << ",n_evals";
        Genetics::Vector<Genetics::String> best_data = pop->get_best_header();
        for (auto it = best_data.begin(); it != best_data.end(); ++it) {
          ofs << "," << *it;
        }
      }
      Genetics::Vector<Genetics::String> data = pop->get_header();
      for (auto it = data.begin(); it != data.end(); ++it) {
        ofs << "," << *it;
      }
      ofs << "\n";
      gen = 0;
    } else {
      ofs << gen;
      if (include_best) {
        Genetics::Vector<Genetics::String> best_data = pop->get_best_data();
        //ofs << "," << pop->get_best_organism()->get_fitness_info().get_n_evaluations();
        for (auto it = best_data.begin(); it != best_data.end(); ++it) {
          ofs << "," << *it;
        }
      }
      Genetics::Vector<Genetics::String> data = pop->get_pop_data();
      for (auto it = data.begin(); it != data.end(); ++it) {
        ofs << "," << *it;
      }
      ofs << "\n";
      ++gen;
    }
  }
};

class ChromosomeTestFixture {
protected:
  Genetics::ArgStore args;
  Genetics::Chromosome chrom_96_bits_A, chrom_96_bits_B;//this should be large enough to require the use of two words
  Genetics::Chromosome chrom_64_bits_A, chrom_64_bits_B;
  Genetics::Chromosome chrom_32_bits_A, chrom_32_bits_B;
  Genetics::PhenotypeMap map_96;
  Genetics::PhenotypeMap map_64;
  Genetics::PhenotypeMap map_32;
public:
  ChromosomeTestFixture() : chrom_96_bits_A(96), chrom_96_bits_B(96), chrom_64_bits_A(64), chrom_64_bits_B(64), chrom_32_bits_A(32), chrom_32_bits_B(32), map_96(96), map_64(64), map_32(32) {}
};

class TestProblemPenalties : public Genetics::Problem<Genetics::SingleFitness> {
private:
  bool apply_penalty;
  bool use_costs = false;
  double penalty_amt;
public:
  TestProblemPenalties(double p_penalty_amt, bool p_use_costs) :
  Genetics::Problem<Genetics::SingleFitness>(NUM_BITS, NUM_CHROMS, 1) {
    map->initialize(1, Genetics::t_real);
    apply_penalty = true;
    penalty_amt = p_penalty_amt;
    use_costs = p_use_costs;
  }
  void evaluate_fitness(OrganismSingle* org) {
    if (use_costs) {
      org->set_cost(1);
    } else {
      org->update(1);
    }
    
    if (apply_penalty) {
      org->apply_penalty(penalty_amt);
      apply_penalty = false;
    }
  }
};

class TestProblemAvgs : public Genetics::Problem<Genetics::NoisyFitness> {
public:
  LCG gen;
  Genetics::Vector<double> fit_vals;

public:
  TestProblemAvgs(unsigned long seed = DEFAULT_LCG_SEED) :
  Genetics::Problem<Genetics::NoisyFitness>(NUM_BITS, NUM_CHROMS, 1), gen(seed) {}

  void evaluate_fitness(OrganismNoisy* org) {
    double val = gen.random_real();
    fit_vals.push_back(val);
    org->update(val);
  }

  double get_val(size_t i) { return fit_vals[i]; }
  void get_mean_var(double* mean, double* var) {
    *mean = 0.0;
    if ( fit_vals.size() > 0 ) {
      for (auto it = fit_vals.begin(); it != fit_vals.end(); ++it) {
        *mean += (*it);
      }
      *mean /= fit_vals.size();
      *var = 0.0;
      if ( fit_vals.size() > 1 ) {
        for (auto it = fit_vals.begin(); it != fit_vals.end(); ++it) {
          *var += (*it - *mean)*(*it - *mean);
        }
        *var /= (fit_vals.size() - 1);
      }
    }
  }
};

// =============================== CHROMOSOME TEST CASES ===============================

TEST_CASE_METHOD( ChromosomeTestFixture, "Ensure that chromosomes are correctly encoded and decoded for integers", "[chromosomes]" ) {
  REQUIRE(chrom_96_bits_A.get_n_bits() == 96);
  REQUIRE(chrom_96_bits_B.get_n_bits() == 96);
  REQUIRE(chrom_64_bits_A.get_n_bits() == 64);
  REQUIRE(chrom_64_bits_B.get_n_bits() == 64);
  REQUIRE(chrom_32_bits_A.get_n_bits() == 32);
  REQUIRE(chrom_32_bits_B.get_n_bits() == 32);
  map_96.initialize(2, Genetics::t_uint);
  map_64.initialize(2, Genetics::t_uint);
  map_32.initialize(2, Genetics::t_uint);
  
  unsigned long mask_96 = ( (unsigned long)1 << (96/2) ) - 1;
  unsigned long mask_64 = ( (unsigned long)1 << (64/2) ) - 1;
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
    val = generator.random_ulong() & mask_64;
    chrom_64_bits_A.set_to_ulong(&map_64, 0, val);
    INFO("Setting 64 bit genome: i = " << i)
    REQUIRE(chrom_64_bits_A.gene_to_ulong(&map_64, 0) == val);
  }
  for (int i = 0; i < N_TRIALS; ++i) {
    val = generator.random_ulong() & mask_64;
    chrom_64_bits_A.set_to_ulong(&map_64, 1, val);
    INFO("Setting 64 bit genome: i = " << i)
    REQUIRE(chrom_64_bits_A.gene_to_ulong(&map_64, 1) == val);
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
  REQUIRE(chrom_64_bits_A.get_n_bits() == 64);
  REQUIRE(chrom_64_bits_B.get_n_bits() == 64);
  REQUIRE(chrom_32_bits_A.get_n_bits() == 32);
  REQUIRE(chrom_32_bits_B.get_n_bits() == 32);
  //initialize the 96 bit map
  map_96.initialize(2, Genetics::t_real);
  map_96.set_range(0, 0, 100);
  map_96.set_range(1, 0, 1);
  //initialize the 64 bit map
  map_64.initialize(2, Genetics::t_real);
  map_64.set_range(0, 0, 100);
  map_64.set_range(1, 0, 1);
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
    chrom_64_bits_A.set_to_num(&map_64, 0, val);
    INFO("Setting 64 bit genome: i = " << i << " val = " << val << " rep = " << chrom_64_bits_A.gene_to_num(&map_64, 0))
    REQUIRE( APPROX(chrom_64_bits_A.gene_to_num(&map_64, 0), val) );
  }
  for (int i = 0; i < N_TRIALS; ++i) {
    val = generator.random_real();
    chrom_64_bits_A.set_to_num(&map_64, 1, val);
    INFO("Setting 64 bit genome: i = " << i << " val = " << val << " rep = " << chrom_64_bits_A.gene_to_num(&map_64, 1))
    REQUIRE(APPROX(chrom_64_bits_A.gene_to_num(&map_64, 1), val));
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
  REQUIRE(chrom_64_bits_A.get_n_bits() == 64);
  REQUIRE(chrom_64_bits_B.get_n_bits() == 64);
  REQUIRE(chrom_32_bits_A.get_n_bits() == 32);
  REQUIRE(chrom_32_bits_B.get_n_bits() == 32);
  map_64.initialize(1, Genetics::t_uint);
  map_96.initialize(2, Genetics::t_uint);
  map_32.initialize(1, Genetics::t_uint); 
  int single_size_96 = 96 / 2;// = 48

  for (int i = 0; i < 64; ++i) {
    chrom_64_bits_A.set_to_ulong(&map_64, 0, ULONG_MAX);
    REQUIRE( chrom_64_bits_A.gene_to_ulong(&map_64, 0) == ULONG_MAX );
    chrom_64_bits_B.set_to_ulong(&map_64, 0, 0);
    //perform the exchange
    chrom_64_bits_A.exchange(&chrom_64_bits_B, i);
    unsigned long tmp_long = chrom_64_bits_A.gene_to_ulong(&map_64, 0);
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
  mask = ((unsigned long)1 << single_size_96) - 1;
  for (int i = 0; i < single_size_96; ++i) {
    chrom_96_bits_A.set_to_ulong(&map_96, 1, ULONG_MAX);
    REQUIRE( chrom_96_bits_A.gene_to_ulong(&map_96, 1) == (ULONG_MAX & mask) );
    chrom_96_bits_B.set_to_ulong(&map_96, 1, 0);
    //perform the exchange
    chrom_96_bits_A.exchange(&chrom_96_bits_B, i + single_size_96);

    unsigned long tmp_long = chrom_96_bits_A.gene_to_ulong(&map_96, 1);
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
    Genetics::PhenotypeMap map_mixed(31);
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
    Genetics::PhenotypeMap map(31);

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

// =============================== ORGANISM TEST CASES ===============================

TEST_CASE ("Accumulated averages and standard deviations work") {
  TestProblemAvgs avg;
  TestProblemAvgs avg2(0xA1B3C5D7E9FBABA5);
  OrganismNoisy org1(1, 1, std::make_shared<Genetics::PhenotypeMap>(1));
  OrganismNoisy org2(1, 1, std::make_shared<Genetics::PhenotypeMap>(1));
  double mean, var;
  double org_mean, org_stdev;
  for (int i = 0; i < N_TRIALS; ++i) {
    org1.evaluate_fitness(&avg);
    org2.evaluate_fitness(&avg2);
    avg.get_mean_var(&mean, &var);
    org_mean = org1.get_fitness();
    org_stdev = org1.get_fitness_info().get_uncertainty();

    REQUIRE( org1.get_fitness_info().get_n_evaluations() == i+1 );
    INFO( "i=" << i << " accumulated fitness: " << org_mean << ", actual: " << mean << ", accumulated variance: " << org_stdev << ", actual: " << var)
    REQUIRE( APPROX(org_mean, mean) );
    INFO( "i=" << i << " accumulated fitness: " << org_mean << ", actual: " << mean << ", accumulated variance: " << org_stdev << ", actual: " << var)
    REQUIRE( APPROX(org_stdev*org_stdev, var) );
  }

  //use the more memory intensive but sure-fire way of calculating mean and variance
  mean = 0;
  var = 0;
  for (int i = 0; i < N_TRIALS; ++i) {
    mean += avg.get_val(i);
    mean += avg2.get_val(i);
  }
  mean /= 2*N_TRIALS;
  for (int i = 0; i < N_TRIALS; ++i) {
    var += (avg.get_val(i) - mean)*(avg.get_val(i) - mean);
    var += (avg2.get_val(i) - mean)*(avg2.get_val(i) - mean);
  }
  var /= (2*N_TRIALS - 1);
  
  Genetics::NoisyFitness org1_fit_val = org1.get_fitness_info();
  Genetics::NoisyFitness org2_fit_val = org2.get_fitness_info();
  org1_fit_val.average_fitness(org2_fit_val);
  org_mean = org1_fit_val.get_fitness();
  org_stdev = org1_fit_val.get_uncertainty();
  std::cout << "Using accumulated fitnesses with " << 2*N_TRIALS << " trials produced a value of " << org_mean << " and a variance " << org_stdev << ". The accepted values were " << mean << " and " << var << " respectively. Abs_mu=" << (org_mean - mean) << ", %_mu=" << 100*(org_mean - mean)/mean << ", Abs_var=" << (org_stdev - var) << ", %_var=" << 100*(org_stdev - var)/var << std::endl;
  INFO( "accumulated fitness: " << org_mean << ", actual: " << mean << ", accumulated variance: " << org_stdev << ", actual: " << var)
  REQUIRE( APPROX(org_mean, mean) );
  INFO( "accumulated fitness: " << org_mean << ", actual: " << mean << ", accumulated variance: " << org_stdev << ", actual: " << var)
  REQUIRE( APPROX(org_stdev, sqrt(var)) );
}

TEST_CASE ( "Organisms are correctly created and decoded", "[organisms]" ) {
  std::shared_ptr<Genetics::PhenotypeMap> map_16 = std::make_shared<Genetics::PhenotypeMap>(16);
  std::shared_ptr<Genetics::PhenotypeMap> map_32 = std::make_shared<Genetics::PhenotypeMap>(32);
  OrganismSingle template_16_bits_1_params(16, 1, map_16);
  OrganismSingle template_32_bits_2_params(32, 2, map_32);
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

// =============================== POPULATION TEST CASES ===============================

#define POP_TEST_MIN	-10
#define POP_TEST_MAX	10
TEST_CASE ("Populations are correctly created and data is successfully read", "[populations]") {
//  Genetics::ArgStore inst;
  Genetics::ArgStore args;
  args.initialize_from_file("ga.conf");
//  inst.initialize_from_file(conf_file.c_str());
  std::shared_ptr<Genetics::PhenotypeMap> map = std::make_shared<Genetics::PhenotypeMap>(NUM_BITS);
  map->initialize(NUM_CHROMS, Genetics::t_real);
  map->set_range(0, POP_TEST_MIN, POP_TEST_MAX);
  OrganismSingle tmplt(NUM_BITS, NUM_OBJS, map);
  tmplt.set_real(0, 4.6);
  PopulationSingle pop_none(NUM_BITS, NUM_OBJS, &tmplt, map, args);
  PopulationSingle pop_tmplt(NUM_BITS, NUM_OBJS, map, args);

  SECTION ( "The header is read correctly" ) {
    Genetics::Vector<Genetics::String> pop_dat = pop_none.get_header();
    unsigned span = NUM_CHROMS + NUM_OBJS;
    for (unsigned i = 0; i < pop_dat.size() / span; ++i) {
      for (unsigned j = 0; j < NUM_CHROMS; ++j) {
        std::stringstream oss;
        oss << "x_" << j;
        REQUIRE( pop_dat[i*span + j] == oss.str() );
      }
      for (unsigned j = 0; j < NUM_OBJS; ++j) {
        std::stringstream oss;
        oss << "f_" << j << "(x)";
        REQUIRE( pop_dat[i*span + NUM_CHROMS + j] == oss.str() );
      }
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
    REQUIRE( pop_dat.size() % (NUM_CHROMS + NUM_OBJS) == 0 );
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

  SECTION ( "Latin-Hypercube initialization works" ) {
    unsigned n_orgs = (POP_TEST_MAX - POP_TEST_MIN) / 2;
    args.set_pop_size(n_orgs);
    for (unsigned dims = 1; dims < 5; ++dims) {
      std::shared_ptr<Genetics::PhenotypeMap> map_multi = std::make_shared<Genetics::PhenotypeMap>(NUM_BITS_LARGE);
      map_multi->initialize(dims, Genetics::t_real);
      for (unsigned i = 0; i < dims; ++i) {
        map_multi->set_range(i, POP_TEST_MIN, POP_TEST_MAX);
      }
      PopulationSingle pop(NUM_BITS, NUM_OBJS, map_multi, args, true);
      for (unsigned k = 0; k < dims; ++k) {
        for (unsigned i = 0; i < pop.get_offspring_num(); ++i) {
          unsigned row_i = (unsigned)floor( pop.get_offspring_num()*(pop.get_organism(i)->read_real(k) - POP_TEST_MIN) / (POP_TEST_MAX - POP_TEST_MIN) );
          for (unsigned j = i+1; j < pop.get_offspring_num(); ++j) {
            unsigned row_j = (unsigned)floor( pop.get_offspring_num()*(pop.get_organism(j)->read_real(k) - POP_TEST_MIN) / (POP_TEST_MAX - POP_TEST_MIN) );
            INFO ( "i=" << i << " ( " << pop.get_organism(i)->read_real(k) << " ), j=" << j << " ( " << pop.get_organism(j)->read_real(k) << " ) , collision in x_" << k << " " << POP_TEST_MIN << "<x<" << POP_TEST_MAX << ", n_orgs=" << pop.get_offspring_num() )
            REQUIRE( row_i != row_j );
          }
        }
      }
    }
  }
}

TEST_CASE ("Penalties are applied properly", "[populations]") {
  double penalty_str = 4.0;
  TestProblemPenalties prob_fit(penalty_str, false);
  TestProblemPenalties prob_cost(penalty_str, true);
  Genetics::ArgStore args;
  args.initialize_from_file("ga.conf");
  PopulationSingle pop_fit( NUM_BITS, 1, prob_fit.map, args);
  PopulationSingle pop_cost( NUM_BITS, 1, prob_cost.map, args);
  pop_fit.evaluate(&prob_fit);
  pop_cost.evaluate(&prob_cost);

  double fit_0 = pop_fit.get_organism(0)->get_fitness(0);
  double cost_0 = pop_cost.get_organism(0)->get_cost(0);
  for (unsigned i = 1; i < pop_fit.get_offspring_num() && i < pop_cost.get_offspring_num(); ++i) {
    double fit_i = pop_fit.get_organism(i)->get_fitness(0);
    double cost_i = pop_cost.get_organism(i)->get_cost(0);
    INFO( "i = " << i )
    REQUIRE( fit_i == 1.0 );
    REQUIRE( fit_0 < fit_i );
    REQUIRE( cost_i == 1.0 );
    REQUIRE( cost_0 > fit_i );
  }
}

// =============================== ARG_STORE TEST CASES ===============================

TEST_CASE ("ArgStore successfully parses a file") {
  TestProblemSingle prob;
  Genetics::Conv_Plateau plat_cut(0.05, TEST_CONV_GEN / 2);
  Genetics::ArgStore args;
  args.initialize_from_file("ga.conf");
  PopulationSingle pop( NUM_BITS, 1, prob.map, args);

  //ga.conf must be as follows for this to work:
  // population_size: 50
  // arena_size: 2
  // num_generations: 50
  // num_crossovers: 2
  // parameter_variance: 0.3
  // parameter_mean: 0.0
  // mutation_probability: 0.016
  // crossover_probability: 0.9
  // hypermutation_threshold: 1.0
  // output_file: output.csv
  REQUIRE( pop.get_args().get_pop_size() == 50 );
  REQUIRE( pop.get_args().read_custom_double("arena_size", 2) );
  REQUIRE( pop.get_args().get_num_gens() == 50 );
  REQUIRE( pop.get_args().get_num_crossovers() == 2 );
  REQUIRE( pop.get_args().get_init_param_var() == 0.3 );
  REQUIRE( pop.get_args().get_crossover_prob() == 0.9 );
  REQUIRE( pop.get_args().get_mutate_prob() == 0.016 );
  REQUIRE( pop.get_args().get_hypermutation_threshold() == 1.0 );
  REQUIRE( pop.get_args().get_out_fname() == "output.csv" );
}

#define N_VAR_TRIALS  5
TEST_CASE ("Noisy population evaluations don't break") {
  PopulationPrinter<PopulationNoisy>* out;
  char out_name[50];
  for (int n = 0; n < N_VAR_TRIALS; ++n) {
    double penalty_domain = NOISY_DOMAIN*(1 - pow(0.9, n));
    std::cout << "Now running noisy test with domain: " << NOISY_DOMAIN << " and penalty domain: " << penalty_domain << std::endl;
    TestProblemNoisy prob(NOISY_DOMAIN, penalty_domain, NOISY_VAR);
    Genetics::ArgStore args;
    args.initialize_from_file("ga_noisy.conf");
    PopulationNoisy pop( NUM_BITS, 1, prob.map, args);
    pop.set_cost(0);
    snprintf(out_name, 50, "noisy_output_%d.csv", n);
    out = new PopulationPrinter<PopulationNoisy>(&pop, out_name);
    out->print_line();

    //evaluate for 10 generations
    for (int gen = 0; gen < 20; ++gen) {
      pop.evaluate(&prob);
      std::shared_ptr<OrganismNoisy> best_org = pop.get_best_organism();
      double best_fitness = best_org->get_fitness(0);
      Genetics::FitnessStats pop_stats = pop.get_pop_stats();
      INFO( "gen=" << gen )
      REQUIRE( pop_stats.max == best_fitness );

      bool best_found = false;
      bool best_in_pop = false;
      for (int i = 0; i < pop.get_offspring_num(); ++i) {
        std::shared_ptr<OrganismNoisy> org_i = pop.get_organism(i);
        //check for information about the fitness relative to the best
        INFO( "i=" << i << ", best fitness=" << best_fitness << ", current observed=" << org_i->get_fitness(0) )
        REQUIRE( (org_i->get_fitness(0) < best_fitness || APPROX(org_i->get_fitness(0), best_fitness)) );
        if ( org_i->get_fitness(0) == best_fitness ) {
          best_found = true;
        }
        //check whether this organism has the same genotype as the best organism
        if ( *best_org == *org_i ) {
          std::cout << "\ti = " << i << ", value = " << org_i->read_real(0) << std::endl;
          best_in_pop = true;
        }
        //if the organism is penalized ensure that it is less fit than every other organism
        if (org_i->penalized()) {
          for (int j = 0; j < pop.get_offspring_num(); ++j) {
            if ( !(pop.get_organism(j)->penalized()) ) {
              INFO("j=" << j)
              REQUIRE( org_i->get_fitness() < pop.get_organism(j)->get_fitness() );
            }
          }
        }
      }
      INFO( "gen=" << gen << ", best_org.value[i]=" << best_org->read_real(0) << ", fitness=" << best_org->get_fitness() )
      REQUIRE( best_in_pop );

      REQUIRE( !(pop.get_best_organism()->penalized()) );
      if (gen == 0) {
        REQUIRE( best_found );
      }
      if (best_found) {
        std::cout << "\timprovement in generation " << gen << ". New fitness = " << best_fitness << std::endl;
      }
      out->print_line();
    }
    std::vector<std::pair<std::shared_ptr<OrganismNoisy>, unsigned>> species = pop.get_species_list();
    std::cout << "total number of organisms: " << pop.get_offspring_num() << ", number of species: " << species.size();
    for (unsigned i = 0; i < species.size(); ++i) {
      std::cout << "\t species_" << i << ": " << 100*(double)(species[i].second)/pop.get_offspring_num() << "%\n";
    }
    delete out;
  }
}

// =============================== CONVERGENCE TEST CASES ===============================

TEST_CASE ("Simple evolution of a single objective converges to roughly appropriate result", "[populations]") {
  TestProblemSingle prob;
  Genetics::ArgStore args;
  args.initialize_from_file("ga.conf");
  PopulationSingle pop( NUM_BITS, 1, prob.map, args);
  PopulationPrinter<PopulationSingle> out(&pop, "simple_evo.csv");
  out.print_line();
  for (size_t i = 0; i < 100; ++i) {
    pop.evaluate(&prob);
    pop.iterate();
    out.print_line();
  }
  pop.evaluate(&prob);
  REQUIRE( APPROX(pop.get_best_organism()->get_fitness(0), 0) );
}

typedef Genetics::Population<Genetics::MultiFitness, Genetics::NSGAII_TournamentSelector<Genetics::MultiFitness>> Population_NSGAII;
TEST_CASE ("Simple evolution of a multi objective converges to roughly appropriate result", "[populations]") {
  TestProblemMulti prob;
  Genetics::ArgStore args;
  args.initialize_from_file("ga.conf");

  Population_NSGAII pop( NUM_BITS, NUM_OBJS, prob.map, args);
  pop.set_cost(0);
  pop.set_cost(1);
  PopulationPrinter<Population_NSGAII> out(&pop, "multi_objective.csv");
  out.print_line(false);
  bool converged = false;
  int gen = 0;
  while (!converged) {
    pop.evaluate(&prob);
    for (int i = 0; i < pop.get_offspring_num(); ++i) {
      //the fitness function we are using makes this impossible
      INFO("organism " << i << " in generation " << gen << " has invalid fitness (x-2)^2=x^2=0")
      REQUIRE( (!APPROX(pop.get_organism(i)->get_fitness(0), 0) || !APPROX(pop.get_organism(i)->get_fitness(1), 0)) );
    }
    out.print_line(false);
    converged = pop.iterate();
    ++gen;
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
  Genetics::Vector<Genetics::FitnessStats> s(2);
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
      REQUIRE( !var_cut.evaluate_convergence(s) );
      INFO("i = " << i << " range_1 = " << s[0].max - s[0].min << " range_2 = " << s[1].max - s[1].min)
      REQUIRE( !range_cut.evaluate_convergence(s) );
    } else {
      INFO("i = " << i << " var_1 = " << s[0].var << " var_2 = " << s[1].var)
      REQUIRE( var_cut.evaluate_convergence(s) );
      INFO("i = " << i << " range_1 = " << s[0].max - s[0].min << " range_2 = " << s[1].max - s[1].min)
      REQUIRE( range_cut.evaluate_convergence(s) );
    }
  }
}

TEST_CASE ("Combine convergence checking and evolution", "[populations]") {
  TestProblemSingle prob;
  Genetics::Conv_Plateau plat_cut(0.05, TEST_CONV_GEN / 2);
  Genetics::ArgStore args;
  args.initialize_from_file("ga.conf");
  PopulationSingle pop( NUM_BITS, 1, prob.map, args);
  PopulationPrinter<PopulationSingle> out(&pop, "convergence_test.csv");
  out.print_line();

  bool converged = false;
  int generation = 0;
  while (!converged) {
    pop.evaluate(&prob);
    converged = pop.iterate(&plat_cut);
    generation++;
    out.print_line();
  }
  std::cout << "Converged after " << generation << " generations\n";
  pop.evaluate(&prob);
  out.print_line();
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
  //unicode for \pm
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
      PopulationSingle pop_serial( NUM_BITS, 1, prob.map, args_serial);
      pop_serial.run(&prob);
    auto t_end = std::chrono::high_resolution_clock::now();
    serial_times[i] = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    serial_mean += serial_times[i] / N_TRIALS;

    //evaluate async
    t_start = std::chrono::high_resolution_clock::now();
      PopulationSingle pop_async( NUM_BITS, 1, prob.map, args_async);
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
