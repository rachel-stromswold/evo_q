#include "phenotype.h"

namespace Genetics {
  
PhenotypeMap::PhenotypeMap(int pn_bits) : N_BITS(pn_bits) {
}

PhenotypeMap::PhenotypeMap(const PhenotypeMap &obj) {
  N_BITS = obj.N_BITS;
  n_dims = obj.n_dims;
  vars = obj.vars;
}

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
    error(0, "The genome does not make full use of available data, consider reducing the number of allocated bits.");
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

void PhenotypeMap::set_range(_uint ind, double min, double max) {
  if (!is_real(ind)) {
    error(0, "Index %d does not have a real value.");
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

double PhenotypeMap::get_range_min(_uint ind) {
  if (ind > n_dims) {
    error(1, "Attempted to access invalid range minimum index=%d", ind);
  }
  if (!is_real(ind)) {
    error(0, "Index %d does not have a real value.");
  }
  return vars[ind].range_lo;
}

double PhenotypeMap::get_range_max(_uint ind) {
  if (ind > n_dims) {
    error(1, "Attempted to access invalid range maximum index=%d", ind);
  }
  if (!is_real(ind)) {
    error(0, "Index %d does not have a real value.");
  }
  return vars[ind].range_hi;
}

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

_uint PhenotypeMap::get_block_location(_uint ind) {
  if (ind > n_dims) {
    error(1, "Attempted to access invalid block index=%d", ind);
  }

  return vars[ind].loc;
}

unsigned int PhenotypeMap::get_block_length(_uint ind) {
  if (ind > n_dims) {
    error(1, "Attempted to access invalid block index=%d", ind);
  }

  return vars[ind+1].loc - vars[ind].loc;
}

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

double PhenotypeMap::get_factor(_uint ind) {
  if (!is_real(ind)) {
    error(0, "Attempt to access real conversion factor for index that is not real valued. index = %d", ind);
    return 1.0;
  }
  return vars[ind].factor;
}

bool PhenotypeMap::is_real(_uint ind) {
  if (ind > n_dims) {
    error(0, "Invalid index %d", ind);
    return false;
  }
  return vars[ind].type == t_real;
}

bool PhenotypeMap::is_int(_uint ind) {
  if (ind > n_dims) {
    error(0, "Invalid index %d", ind);
    return false;
  }
  return vars[ind].type == t_int;
}

bool PhenotypeMap::is_uint(_uint ind) {
  if (ind > n_dims) {
    error(0, "Invalid index %d", ind);
    return false;
  }
  return vars[ind].type == t_uint;
}

bool PhenotypeMap::is_bitstream(_uint ind) {
  if (ind > n_dims) {
    error(0, "Invalid index %d", ind);
    return false;
  }
  return vars[ind].type == t_bitstream;
}

}
