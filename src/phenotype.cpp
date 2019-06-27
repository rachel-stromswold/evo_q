#include "phenotype.h"

namespace Genetics {
  
/**
 * \brief constructor
 *
 * \param n_bits	The number of bits available for use by the phenotype map
 */
PhenotypeMap::PhenotypeMap(int pn_bits) : N_BITS(pn_bits) {
}

/**
 * \brief 	copy constructor
 *
 * \param obj The PhenotypeMap reference to be copied
 */
PhenotypeMap::PhenotypeMap(const PhenotypeMap &obj) {
  N_BITS = obj.N_BITS;
  n_dims = obj.n_dims;
  vars = obj.vars;
}

/**
 * \brief Initializes the Phenotype map to store parameters based on the input string
 *
 * \details 	This function parses a string of the form
 * 		<n type 1>[type 1]<n type 2>[type 2]... where <n type m> is an integer and
 * 		[type 1] is a character code corresponding to a Type enum. Available character
 * 		codes are 'd' or 'i' for integer parameters, 'r' for real valued parameters or
 * 		'b' for generic bitstrings. For example, the input string "1d3r" will initialize 
 * 		the PhenotypeMap object to store 4 parameters, one of which is an integer while
 * 		the other three parameters may take on real values. When specifying real values,
 * 		the user should call the set_range() method to specify the range of values for
 * 		the given parameter.
 *
 * \param str	A string in the format specified above.
 */
void PhenotypeMap::parse_string(std::string str) {
  char t;
  std::vector<_uint> bstream_lens;
  _uint n_ints = 0, n_reals = 0;
  for (std::string::iterator it = str.begin(); it != str.end(); ++it) {
    t = *it;
    int n = 1;
    //if t is a digit, then convert to a number and store the result in n
    if (t >= '0' && t <= '9') {
      std::string n_tmp = read_number(&it);
      n = std::stoi(n_tmp);
    }
    //read the type to be inserted
    if (t == 'd' || t == 'i') {
      vars.insert(vars.end(), n, VarContainer(0, 0.0, 1.0, t_int));
      n_ints++;
    } else if (t == 'r') {
      double hi, lo;
      //when inserting a double we need to specify a range
      std::string tmp = read_number(&it);
      lo = std::stod(tmp);
      tmp = read_number(&it);
      hi = std::stod(tmp);
      vars.insert(vars.end(), n, VarContainer(0, lo, hi, t_real));
      n_reals++;
    } else if (t == 'b') {
      vars.insert(vars.end(), n, VarContainer(0, 0.0, 1.0, t_bitstream));
      std::string tmp = read_number(&it);
      bstream_lens.push_back( (_uint)std::stoul(tmp) );
    }
    n_dims += n;
  }
  allocate_locations(n_ints, n_reals, bstream_lens);
}

/**
 * \brief Initialize the phenotype map using a vector of VarContainer objects
 *
 * \details	This method initializes the phenotype map to contain the types specified in
 * 		the argument vc. If vc cointains any real valued parameters, then the range will
 * 		be set appropriately without the need for a separate call to the set_range()
 * 		method.
 * \param vc	A vector of VarContainer objects used to initialize the map
 */
void PhenotypeMap::initialize(Vector<VarContainer> vc) {
  std::vector<_uint> bstream_lens(0, 0);
  _uint n_real = 0, n_int = 0;
  for (_uint i = 0; i < vc.size(); ++i) {
    if (vc[i].type == t_bitstream) {
      vars.push_back(VarContainer(0, 0, 0, vc[i].type));
      bstream_lens.push_back(vc[i].loc);
    } else if (vc[i].type == t_int || vc[i].type == t_uint) {
      vars.push_back(VarContainer(0, 0.0, 0.0, vc[i].type));
      ++n_int;
    } else if (vc[i].type == t_real) {
      vars.push_back(VarContainer(0, vc[i].range_lo, vc[i].range_hi, vc[i].type));
      ++n_real;
    }
  }
  allocate_locations(n_int, n_real, bstream_lens);
}

/**
 * \brief Initialize the phenotype map using an array of VarContainer objects
 *
 * \details	This method initializes the phenotype map to contain the types specified in
 * 		the argument vc. If vc cointains any real valued parameters, then the range will
 * 		be set appropriately without the need for a separate call to the set_range()
 * 		method.
 * \param vc	A pointer to the first element in the array of VarContainer objects
 * \param n_elemnts	The number of elements contained in the array
 */
void PhenotypeMap::initialize(VarContainer* vc, size_t n_elements) {
  std::vector<_uint> bstream_lens(0, 0);
  _uint n_real = 0, n_int = 0;
  for (_uint i = 0; i < n_elements; ++i) {
    if (vc[i].type == t_bitstream) {
      vars.push_back(VarContainer(0, 0, 0, vc[i].type));
      bstream_lens.push_back(vc[i].loc);
    } else if (vc[i].type == t_int || vc[i].type == t_uint) {
      vars.push_back(VarContainer(0, 0.0, 0.0, vc[i].type));
      ++n_int;
    } else if (vc[i].type == t_real) {
      vars.push_back(VarContainer(0, vc[i].range_lo, vc[i].range_hi, vc[i].type));
      ++n_real;
    }
  }
  allocate_locations(n_int, n_real, bstream_lens);
}

/**
 * \brief Initialize the phenotype map to use n objects of type t
 *
 * \param n	The number of parameters to be used by the PhenotypeMap
 * \param t	The type to be used, this may be t_int, t_uint, t_real or t_bitstream
 */
void PhenotypeMap::initialize(int n, Type t) {
  vars.insert(vars.end(), n, VarContainer(0, 0.0, 1.0, t));
  std::vector<_uint> bstream_lens(0, 0);
  if (t == t_bitstream) {
    bstream_lens.insert(bstream_lens.end(), n, N_BITS/n);
    allocate_locations(0, 0, bstream_lens);
  } else if (t == t_real) {
    allocate_locations(0, n, bstream_lens);
  } else if (t == t_int || t == t_uint) {
    allocate_locations(n, 0, bstream_lens);
  }
}

/**
 * \brief A private helper function called by the various initializers to allocate parameters
 */
void PhenotypeMap::allocate_locations(_uint n_ints, _uint n_reals, std::vector<_uint> bstream_lens) {
  _uint n_avail = N_BITS;
  for (std::vector<_uint>::iterator it = bstream_lens.begin(); it != bstream_lens.end(); ++it) {
    n_avail -= *it;
  }
  if (n_avail <= 0 && (n_ints > 0 || n_reals > 0)) {
    error(1, "Insufficient number of bits in genome representation to store all data.");
  }
  _uint rsize = 2*n_avail/(n_ints+2*n_reals);
  _uint isize = n_avail/(n_ints+2*n_reals);
  if (rsize > MAX_SIZE) {
    rsize = MAX_SIZE;
  }
  if (isize > MAX_SIZE) {
    isize = MAX_SIZE;
  }
  //since it is possible for the ints and reals to not fill up the allocated space we need to add filler to pad out the genome
  if (isize*n_ints + rsize*n_reals < n_avail) {
    error(1, "The genome does not make full use of available data, consider reducing the number of allocated bits.");
  }
  unsigned long real_mask = ( ((unsigned long)0x01 << rsize) - 1 );
  if (rsize == sto_size) {
    real_mask = ULONG_MAX;
  }
  _uint i = 0;
  _uint bstream_ind = 0;
  for (std::vector<VarContainer>::iterator it = vars.begin(); it != vars.end(); ++it) {
    it->loc = i;
    if (it->type == t_int || it->type == t_uint) {
      i += isize;
    } else if (it->type == t_real) {
      i += rsize;
      it->factor = (it->range_hi - it->range_lo)/( (double)real_mask );
    } else if (it->type == t_bitstream) {
      i += bstream_lens[bstream_ind];
      bstream_ind++;
    }
  }
  n_dims = vars.size();
  VarContainer filler(i, 0.0, 1.0, t_terminator);
  vars.push_back(filler);
}

/**
 * \brief 	If the parameter at index ind is of type t_real, this function sets this
 * 		parameter to take a value between min and max.
 *
 * \param ind	The index of the parameter to set
 * \param min	The minimum value the parameter may take
 * \param max	The maximum value the parameter may take
 *
 * \note	Bitstreams are converted into real values by the formula
 * 		gray_decode(n)/(max - min) + min
 * 		where gray_decode performs the gray decoding of the bitstream n into an integer.
 * 		This introduces an inherent and inevitable discretization error with a behavior
 * 		that differs from standard floating point representations. This problem may be
 * 		alleviated by using smaller ranges or allocating more bits to the given
 * 		parameter. The former option is more desirable if possible, as longer bitstreams
 * 		may slow convergence
 *
 * \warning	An error will be thrown if the parameter specified by ind is not real or the
 * 		minimum value is greater than or equal to the maximum
 */
void PhenotypeMap::set_range(_uint ind, double min, double max) {
  if (!is_real(ind)) {
    error(1, "Index %d does not have a real value.");
  }
  if (min >= max) {
    error(1, "The maximum value in the range must be strictly greater than the minimum. ind=%d, range=[%f,%f]", ind, min, max);
  }
  vars[ind].range_lo = min;
  vars[ind].range_hi = max;
  _uint rsize = get_block_length(ind);
  unsigned long real_mask = ( ((unsigned long)0x01 << rsize) - 1 );
  if (rsize == sto_size) {
    real_mask = ULONG_MAX;
  }
  vars[ind].factor = (max - min)/( (double)real_mask );
}

/**
 * \param ind	The index of the parameter
 *
 * \returns	The upper range that the parameter specified by ind may take
 *
 * \warning	An error is thrown if ind is larger than the number of available parameters or
 * 		if the parameter specified by ind is not real
 */
double PhenotypeMap::get_range_min(_uint ind) {
  if (ind > n_dims) {
    error(1, "Attempted to access invalid range minimum index=%d", ind);
  }
  if (!is_real(ind)) {
    error(1, "Index %d does not have a real value.");
  }
  return vars[ind].range_lo;
}

/**
 * \param ind	The index of the parameter
 *
 * \returns	The upper range that the parameter specified by ind may take
 *
 * \warning	An error is thrown if ind is larger than the number of available parameters or
 * 		if the parameter specified by ind is not real
 */
double PhenotypeMap::get_range_max(_uint ind) {
  if (ind > n_dims) {
    error(1, "Attempted to access invalid range maximum index=%d", ind);
  }
  if (!is_real(ind)) {
    error(1, "Index %d does not have a real value.");
  }
  return vars[ind].range_hi;
}

/**
 * \param[in] ind	the index of the parameter to fetch data for
 * \param[out] loc	the low bit index of the parameters representation
 * \param[out] len	the length allocated to the given parameter in bits
 *
 * \warning	This method is designed as a helper for the Chromosome class. In most
 * 		circumstances, the user will probably not want to call this method directly.
 *
 * \returns	The block of bits allocated to a given parameter
 */
void PhenotypeMap::get_block(_uint ind, _uint* loc, _uint* len) {
  if (ind > n_dims) {
    error(1, "Attempted to access invalid block index=%d", ind);
  }
  if (loc != NULL) {
    *loc = vars[ind].loc;
  }
  if (len != NULL) {
    *len = vars[ind+1].loc - vars[ind].loc;
  }
}

/**
 * \param[in] ind	the index of the parameter to fetch data for
 *
 * \warning	This method is designed as a helper for the Chromosome class. In most
 * 		circumstances, the user will probably not want to call this method directly.
 *
 * \returns	The starting bit of the block allocated to a given parameter
 */
_uint PhenotypeMap::get_block_location(_uint ind) {
  if (ind > n_dims) {
    error(1, "Attempted to access invalid block index=%d", ind);
  }

  return vars[ind].loc;
}

/**
 * \param[in] ind	the index of the parameter to fetch data for
 *
 * \warning	This method is designed as a helper for the Chromosome class. In most
 * 		circumstances, the user will probably not want to call this method directly.
 *
 * \returns	The length of the block allocated to a given parameter
 */
unsigned int PhenotypeMap::get_block_length(_uint ind) {
  if (ind > n_dims) {
    error(1, "Attempted to access invalid block index=%d", ind);
  }

  return vars[ind+1].loc - vars[ind].loc;
}

/**
 * \details 	This method is designed as a utility to make writing to and reading from 
 * 		bitstreams easier. Each parameter is stored in a string of bits that may span
 * 		across multiple words. As such, when reading or writing values to memory, we need
 * 		a mask for the bits stored in the high word and the low word.
 *
 * \param[in] ind	the index of the parameter to fetch data for
 * \param[out] loc	a pointer to store the location of the given parameter block in
 * \param[out] p_off	a pointer to store the offset (the location mod word_size) of the given
 * 			parameter block in
 * \param[out] lmask	a pointer to store the bitmask that will be mapped onto the low word
 * \param[out] hmask	a pointer to store the bitmask that will be mapped onto the high word
 *
 * \warning	This method is designed as a helper for the Chromosome class. In most
 * 		circumstances, the user will probably not want to call this method directly.
 */
void PhenotypeMap::get_masks(_uint ind, _uint* loc, _uint* p_off, unsigned long* lmask, unsigned long* hmask) {
  _uint len;
  get_block(ind, loc, &len);

  _uint off = *loc % sto_size;
  *lmask = (off == 0)? ULONG_MAX : ((unsigned long)1 << (sto_size - off)) - 1;
  if (len < sto_size) {
    *lmask &= ((unsigned long)1 << len) - 1; 
  }

  if (len < sto_size - off) {
    *hmask = 0;
  } else {
    *hmask = ( (unsigned long)1 << (len - (sto_size - off)) ) - 1;
  }

  *p_off = off;
}

/**
 * \details 	This function is similar to get_masks, but only returns the low mask. The value
 * 		x can be written to memory by setting
 * 			low = (low & ~get_low_mask(ind)) | (x & get_low_mask(ind))
 * 			high = (high & ~get_high_mask(ind)) | (x & get_high_mask(ind))
 * 		where low is the low word in memory and high is the high word in memory
 *
 * \param ind	The index of the parameter we are interested in.
 *
 * \returns	The low mask associated with the given parameter.
 */
unsigned long PhenotypeMap::get_low_mask(_uint ind) {
  _uint loc, len;
  get_block(ind, &loc, &len);

  _uint off = loc % sto_size;
  unsigned long lmask = ((unsigned long)1 << (sto_size - off)) - 1;
  if (len < sto_size) {
    lmask &= ((unsigned long)1 << len) - 1;
  }
  return lmask;
}

/**
 * \details 	This function is similar to get_masks, but only returns the high mask. The value
 * 		x can be written to memory by setting
 * 			low = (low & ~get_low_mask(ind)) | (x & get_low_mask(ind))
 * 			high = (high & ~get_high_mask(ind)) | (x & get_high_mask(ind))
 * 		where low is the low word in memory and high is the high word in memory
 *
 * \param ind	The index of the parameter we are interested in.
 *
 * \returns	The high mask associated with the given parameter.
 */
unsigned long PhenotypeMap::get_high_mask(_uint ind) {
  _uint loc, len;
  get_block(ind, &loc, &len);

  _uint off = loc % sto_size;
  if (len < sto_size - off) {
    return 0;
  }
  unsigned long mask = ((unsigned long)1 << off) - 1;
  if (len < sto_size) {
    mask &= ( (unsigned long)1 << (len - (sto_size - off)) ) - 1;
  }
  return mask;
}

/**
 * \details	Strings of bits are converted into real values using the formula
 * 			x = gray_decode(bitstream)*(largest_possible_int)/(max - min) + min
 * 		where gray_decode converts a gray coded bitstream into an integer, min and max
 * 		are the largest real values that x may assume respectively, and 
 * 		largest_possible_int is the maximum integer value that may be stowed in n bits
 * 		(2^n - 1). This function returns the factor largest_possible_int/(max - min) for
 * 		the sake of easing conversions.
 *
 * \returns	The factor (largest_possible_int)/(max - min) as described above.
 */
double PhenotypeMap::get_factor(_uint ind) {
  if (!is_real(ind)) {
    error(1, "Attempt to access real conversion factor for index that is not real valued. index = %d", ind);
    return 1.0;
  }
  return vars[ind].factor;
}

/** 
 * \returns	True if the parameter at index ind is real valued (t_real) and false otherwise.
 */
bool PhenotypeMap::is_real(_uint ind) {
  if (ind > n_dims) {
    error(1, "Invalid index %d", ind);
    return false;
  }
  return vars[ind].type == t_real;
}

/** 
 * \returns	True if the parameter at index ind is integer valued (t_int) and false otherwise.
 */
bool PhenotypeMap::is_int(_uint ind) {
  if (ind > n_dims) {
    error(1, "Invalid index %d", ind);
    return false;
  }
  return vars[ind].type == t_int;
}

/** 
 * \returns	True if the parameter at index ind stores an unsigned integer (t_uint) and false otherwise.
 */
bool PhenotypeMap::is_uint(_uint ind) {
  if (ind > n_dims) {
    error(1, "Invalid index %d", ind);
    return false;
  }
  return vars[ind].type == t_uint;
}

/** 
 * \returns	True if the parameter at index ind stores a bitstream (t_bitstream) and false otherwise.
 */
bool PhenotypeMap::is_bitstream(_uint ind) {
  if (ind > n_dims) {
    error(1, "Invalid index %d", ind);
    return false;
  }
  return vars[ind].type == t_bitstream;
}

}
