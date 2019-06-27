#ifndef PHENOTYPE_H
#define PHENOTYPE_H

#include <climits>
#include "util.h"

#define R_FACT    2//by default, reals are allocated twice as much storage space as integers
#define MAX_SIZE  64//doubles and ints only have a maximum of 64 bits anyways, so allocating more memory is pointless

namespace Genetics {

/**
 * \brief	This struct is designed as a wrapper for parameters used by a fitness function.
 * 		Contained are the type, which may be t_int, t_uint, t_real, or t_bitstream for
 * 		integers, unsigned integers, real values or bitstreams respectively. The
 * 		internal variables range_lo, range_hi, and factor are only used for real valued
 * 		parameters where they store the maximum and minimum values that a real value may
 * 		take and a conversion factor respectively. See the methods get_range_min,
 * 		get_range_high and get_factor of the class PhenotypeMap for more details.
 */
struct VarContainer {
  _uint loc;
  double range_hi;
  double range_lo;
  double factor;
  Type type;
  VarContainer(_uint p_loc, double p_lo, double p_hi, Type p_type) {loc = p_loc;range_lo = p_lo;range_hi = p_hi;type = p_type;}
  VarContainer() { loc = 0;range_hi = 0; range_lo = 0; factor = 0; type = t_real; }
};

/**
 * \brief	A class for managing streams of raw bits and converting them into useful
 * 		parameters that may be evaluated by a fitness function.
 *
 * \note	This class does not actually store the data used by organisms. That task is left
 * 		to the Chromosome class.
 */
class PhenotypeMap {
private:
  int N_BITS;

  _uint n_dims = 0;
  std::vector<VarContainer> vars;
//  std::vector<_uint> var_locs;
//  std::vector<std::pair<double, double>> ranges;
//  std::vector<Type> var_types;
  void allocate_locations(_uint n_ints, _uint n_reals, std::vector<_uint> bstream_lens);
public:
  static const _uint sto_size = sizeof(unsigned long)*8;

  PhenotypeMap(int pn_bits);
  PhenotypeMap(const PhenotypeMap &obj);

  void parse_string(std::string str);
  void initialize(int n, Type t);
  void initialize(Vector<VarContainer> vc);
  void initialize(VarContainer* vc, size_t n_elements);

  void set_range(_uint ind, double min, double max);
  double get_range_min(_uint ind);
  double get_range_max(_uint ind);
  _uint get_block_location(_uint ind);
  _uint get_block_length(_uint ind);
  void get_block(_uint ind, _uint* loc, _uint* len);

  void get_masks(_uint ind, _uint* loc, _uint* p_off, unsigned long* lmask, unsigned long* hmask);
  unsigned long get_low_mask(_uint ind);
  unsigned long get_high_mask(_uint ind);
  double get_factor(_uint ind);
  bool is_real(_uint ind);
  bool is_int(_uint ind);
  bool is_uint(_uint ind);
  bool is_bitstream(_uint ind);
  Type get_type(_uint ind) { return vars[ind].type; }
  _uint get_num_params() { if (vars.size() == 0) { return 0; } return vars.size() - 1; }
  int get_n_bits() { return N_BITS; }
};

}

#endif //PHENOTYPE_H
