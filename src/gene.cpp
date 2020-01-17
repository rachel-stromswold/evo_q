#include "gene.h"

namespace Genetics {

Chromosome::Chromosome(_uint pn_bits) :
  N_BITS(pn_bits),
  N_BYTES( (N_BITS+7)/8 ),
  N( (N_BYTES+sizeof(unsigned long)-1)/sizeof(unsigned long) )
{
  genes = (unsigned long*)malloc(sizeof(unsigned long)*N);
  for (unsigned i = 0; i < N; ++i) {
    genes[i] = 0;
  }
  use_real = 0;
}

Chromosome::Chromosome(_uint pn_bits, _uchar real_mode) :
  N_BITS(pn_bits),
  N_BYTES( (N_BITS+7)/8 ),
  N( (N_BYTES+sizeof(unsigned long)-1)/sizeof(unsigned long) )
{
  genes = (unsigned long*)malloc(sizeof(unsigned long)*N);
  for (unsigned i = 0; i < N; ++i) {
    genes[i] = 0;
  }
  use_real = real_mode;
}

Chromosome::Chromosome(_uint pn_bits, Chromosome* other) :
  N_BITS(pn_bits),
  N_BYTES( (N_BITS+7)/8 ),
  N( (N_BYTES+sizeof(unsigned long)-1)/sizeof(unsigned long) )
{
  genes = (unsigned long*)malloc(sizeof(unsigned long)*N);
  for (unsigned i = 0; i < N; ++i) {
    genes[i] = other->genes[i];
  }
  use_real = other->use_real;

  real_vals_size = other->get_real_vals_size();
  if (real_vals_size > 0) {
    real_vals = (double*)malloc(sizeof(double)*real_vals_size);
  } else {
    real_vals = NULL;
  }
}

Chromosome::Chromosome(Chromosome& other) {
  N_BITS = other.N_BITS;
  N_BYTES = other.N_BYTES;
  N = other.N;
  use_real = other.use_real;
  genes = (unsigned long*)malloc(sizeof(unsigned long)*N);
  for (unsigned i = 0; i < N; ++i) {
    genes[i] = other.genes[i];
  }
  real_vals_size = other.real_vals_size;
  if (real_vals_size > 0) {
    real_vals = (double*)malloc(sizeof(double)*real_vals_size);
  } else {
    real_vals = NULL;
  }
}

Chromosome::Chromosome(Chromosome&& other) :
  real_vals(std::move(other.real_vals))
{
  N_BITS = other.N_BITS;
  N_BYTES = other.N_BYTES;
  N = other.N;
  use_real = other.use_real;
  genes = other.genes;
  other.genes = NULL;
  real_vals = other.real_vals;
  real_vals_size = other.real_vals_size;
  other.real_vals = NULL;
  other.real_vals_size = 0;
}

Chromosome::~Chromosome() {
  if (genes) {
    free(genes);
  }
  if (real_vals) {
    free(real_vals);
  }
}

void Chromosome::swap(Chromosome& other) {
  _uint tmp_n_bits = N_BITS;
  _uint tmp_n_bytes = N_BYTES;
  _uint tmp_n = N;
  N_BITS = other.N_BITS;
  N_BYTES = other.N_BYTES;
  N = other.N;
  other.N_BITS = tmp_n_bits;
  other.N_BYTES = tmp_n_bytes;
  other.N = tmp_n;
  unsigned long* tmp = genes;
  genes = other.genes;
  other.genes = tmp;
}

Chromosome& Chromosome::operator=(Chromosome& other) {
  swap(other);
  return *this;
}

Chromosome& Chromosome::operator=(Chromosome&& other) {
  other.swap(*this);
  return *this;
}

bool Chromosome::operator==(Chromosome& other) {
  if (N_BITS == other.N_BITS) {
    for (_uint i = 0; i < N_BITS / bin_size; ++i) {
      if (genes[i] != other.genes[i]) {
	return false;
      }
    }
    if (N_BITS % bin_size != 0) {
      unsigned long mask = (1 << (N_BITS % bin_size)) - 1;
      return (mask & genes[N_BITS / bin_size]) == (mask & other.genes[N_BITS / bin_size]);
    }
    return true;
  }
  return false;
}

void Chromosome::exchange(Chromosome* other, size_t k) {
  if (N_BITS != other->get_n_bits()) {
    error(CODE_ARG_INVALID, "Cannot exchange chromosomes with a differing number of bits, %d and %d.", get_n_bits(), other->get_n_bits());
  }
  size_t j = 0;
  size_t k_ind = k/bin_size;

  if (k < N_BITS) {
    //this should be evaluated at compile time
    if (N == 1) {
      //evil bitwise hacking to swap genes
      unsigned long myg = genes[0];
      unsigned long og = other->genes[0];
      genes[0] = og & (((unsigned long)1 << k) - 1);
      genes[0] = genes[0] | ( myg & ~(((unsigned long)1 << k) - 1) );
      other->genes[0] = myg & (((unsigned long)1 << k) - 1);
      other->genes[0] = other->genes[0] | ( og & ~(((unsigned long)1 << k) - 1) );
    } else {
      for (; j < k_ind; j++) {
        unsigned long tmp = genes[j];
        genes[j] = other->genes[j];
        other->genes[j] = tmp;
      }
      //evil bitwise hacking to swap genes
      unsigned long myg = genes[k_ind];
      _uint l = k % bin_size;
      size_t one = (size_t)1;
      size_t o_mask = ((one << l) - 1);
      genes[k_ind] &= ~o_mask;
      genes[k_ind] |= other->genes[k_ind] & o_mask;
      other->genes[k_ind] &= o_mask;
      other->genes[k_ind] |= (myg & o_mask);
    }
  } else {
    error(CODE_MISC, "tried to copy too many bytes");
  }
}

void Chromosome::exchange_uniform(ArgStore& args, Chromosome* other) {
  if (N_BITS != other->get_n_bits()) {
    error(CODE_MISC, "Cannot exchange chromosomes with a differing number of bits, %d and %d.", get_n_bits(), other->get_n_bits());
  }
  std::binomial_distribution<> n_flipped( bin_size, 0.5);
  for (size_t i = 0; i < N; i++) {
    unsigned char num_ones = n_flipped( args.get_generator() );

    //pick a random bitmask with num_ones bits flipped
    std::uniform_int_distribution<unsigned long> dist( 0, nChoosek(bin_size, num_ones) - 1 );
    unsigned long x = dist( args.get_generator() );
    unsigned long flip_mask = getBitStream( bin_size, num_ones, x );

    unsigned long myg = genes[i];
    unsigned long og = other->genes[i];
    genes[i] = (og & flip_mask) | (myg & ~flip_mask);
    other->genes[i] = (myg & flip_mask) | (og & ~flip_mask) ;
  }
}

void Chromosome::reset() {
  for (size_t i = 0 ; i < N; ++i) {
    genes[i] = 0;
  }
}

unsigned char Chromosome::operator[](unsigned int i) {
  if (i < N_BITS) {
    return ( genes[i/(bin_size)] >> (i%(bin_size)) ) & 1;
  }
  return 0;
}

//this function is a bijection mapping from the set of integers between 0 and n choose k to the set of n-bit integers with binary representation containing k 1 bits. Generating one random number and using this function requires fewer calls to the prng and is thus more efficient than calling the Bernoulli distribution for each bit
size_t Chromosome::getBitStream (size_t n, size_t k, size_t x) {
  if (k == 0) { return 0; }
  if (k == n) { return (0x01 << n)-1; }
	
  if (x < nChoosek(n-1, k-1) ) {
    return 0x01 | (getBitStream(n-1, k-1, x) << 1);
  } else {
    return getBitStream( n-1, k, x-nChoosek(n-1, k-1) ) << 1;
  }
}

bool Chromosome::real_space_mutate(ArgStore& args) {
  if (use_real & REAL_ACTIVE) {
    std::normal_distribution<double> norm(0, args.get_init_param_var());
    for (size_t i = 0; i < real_vals_size; ++i) {
      real_vals[i] += norm( args.get_generator() );
      //ensure that the value is between 0 and 1
      if (real_vals[i] > 1.0) {
	real_vals[i] = 1.0;
      }
      if (real_vals[i] < 0.0) {
	real_vals[i] = 0.0;
      }
    }
    return false;
  }
  return true;
}

void Chromosome::mutate(ArgStore& args) {
  bool perform_bit_mutation = true;
  if (use_real & REAL_ENABLED) {
    perform_bit_mutation = real_space_mutate(args);
  }
  if (perform_bit_mutation) {
    for (size_t i = 0; i < N; i++) {
      unsigned char num_ones = args.sample_binomial( bin_size );

      //pick a random bitmask with num_ones bits flipped
      std::uniform_int_distribution<unsigned long> dist( 0, nChoosek(bin_size, num_ones) - 1);
      unsigned long x = dist(args.get_generator());

      genes[i] = genes[i] ^ getBitStream( bin_size, num_ones, x);
    }
  }
}

void Chromosome::slow_mutate(ArgStore& args) {
  bool perform_bit_mutation = true;
  if (use_real & REAL_ENABLED) {
    perform_bit_mutation = real_space_mutate(args);
  }
  if (perform_bit_mutation) {
    for (unsigned i = 0; i < N - 1; ++i) {
      for (unsigned j = 0; j < bin_size; ++j) {
	if ( args.random_mutation() ) {
	  genes[i] = genes[i] ^ (0x01 << j);
	}
      }
    }
    //we have to do something special for the last long because it isn't filled
    for (unsigned j = 0; j < N_BITS%(bin_size); ++j) {
      if ( args.random_mutation() ) {
	genes[N-1] = genes[N-1] ^ (0x01 << j);
      }
    }
  }
}

void Chromosome::randomize(PhenotypeMap* al, ArgStore& args) {
  std::uniform_int_distribution<unsigned long> dist(0, ULONG_MAX);
  for (size_t i = 0; i < N; i++) {
    genes[i] = dist(args.get_generator());
  }
  /*std::uniform_real_distribution<double> dist(0, 1.0);

  for (size_t i = 0; i < al->get_num_params(); i++) {
    Type t = al->get_type(i);
    double val = dist(args.get_generator());
    if (t == t_bitstream || t == t_uint || t == t_int) {
      _uint l = al->get_block_length(i);
      std::uniform_int_distribution<unsigned int> int_dist(0, 1 << l);
      set_to_ulong( al, i, int_dist( args.get_generator() ) );
    } else if (t == t_real) {
      double max = al->get_range_max(i);
      double min = al->get_range_min(i);
      double val = dist( args.get_generator() )*(max - min) + min;
      set_to_num(al, i, val);
    }
  }*/
}

void Chromosome::set_to_ulong(PhenotypeMap* al, _uint ind, unsigned long value) {
  _uint loc = 0, len = 0, off = 0;
  al->get_block(ind, &loc, &len);

  //get the high and low masks
  unsigned long low_mask = 0, high_mask = 0;
  al->get_masks(ind, &loc, &off, &low_mask, &high_mask);

  _uint i2 = loc / (bin_size);
  if (i2 > N) {
    error(CODE_ARG_INVALID, "Something has gone wrong, invalid location returned.");
  }
  genes[i2] &= ~(low_mask << off);
  genes[i2] |= (value & low_mask) << off;
  if (high_mask) {
    genes[i2 + 1] &= ~(high_mask);
    genes[i2 + 1] |= (value >> (bin_size - off)) & high_mask;
  }
}

void Chromosome::set_to_num(PhenotypeMap* al, _uint ind, double value) {
  double min = al->get_range_min(ind);
  double max = al->get_range_max(ind);
  if (value <= min) {
    value = min;
    set_to_ulong(al, ind, 0);
  } else if (value >= max) {
    //the value 2^(n-1) gets mapped to (2^n) - 1 using the Gray code
    set_to_ulong(al, ind, (unsigned long)1 << (al->get_block_length(ind) - 1));
  } else {
    //the ratio of the difference to the largest possible value
    _uint loc = 0, len = 0;
    al->get_block(ind, &loc, &len);

    double ratio = (value - min)/al->get_factor(ind);
    unsigned long val = encodeGray( (unsigned long)(ratio) );

    set_to_ulong(al, ind, val);
  }
}

void Chromosome::set_to_int (PhenotypeMap* al, _uint ind, int value) {
  _uint loc = 0, len = 0;
  al->get_block(ind, &loc, &len);
  unsigned long val;
  if (value >= 0) {
    val = value;
    val = val << 1;
  } else {
    val = -value;
    val = (val << 1) | 1;
  }
  set_to_ulong(al, ind, val);
}

unsigned long Chromosome::gene_to_ulong(PhenotypeMap* al, _uint ind) {
  _uint loc = 0, off = 0;
  unsigned long low_mask = 0, high_mask = 0;
  al->get_masks(ind, &loc, &off, &low_mask, &high_mask);

  _uint i2 = loc / (bin_size);
  unsigned long ret = (genes[i2] >> off) & low_mask;
  if (high_mask) {
    ret |= (genes[i2 + 1] & high_mask) << (bin_size - off);
  }
  return ret;
}

int Chromosome::gene_to_int(PhenotypeMap* al, _uint ind) {
  unsigned long val = gene_to_ulong(al, ind);
  int ret = val >> 1;
  if (val & 1) {
    ret *= -1;
  }
  return ret;
}

double Chromosome::gene_to_num(PhenotypeMap* al, _uint ind) {
  if ( (use_real & REAL_ACTIVE) && al->is_real(ind) ) {
    double val = (al->get_range_max(ind) - al->get_range_min(ind))*real_vals[ind] + al->get_range_min(ind);
  } else {
    double raw = (double)decodeGray( gene_to_ulong(al, ind) );
    if (!al->is_real(ind)) {
      return raw;
    }
    double min = al->get_range_min(ind);

    return al->get_factor(ind)*raw + min;
  }
  return 0;
}

String Chromosome::get_string(PhenotypeMap* al, _uint ind) {
#ifdef USE_CUSTOM_CONTAINERS
  String iss;
#else
  std::stringstream iss;
#endif
  if (al->is_real(ind)) {
    iss << gene_to_num(al, ind);
  } else if (al->is_int(ind)) {
    iss << gene_to_int(al, ind);
  } else if (al->is_uint(ind)) {
    iss << gene_to_ulong(al, ind);
  } else if (al->is_bitstream(ind)) {
    _uint off, loc, len;
    unsigned long lmask, hmask;
    al->get_masks(ind, &off, &loc, &lmask, &hmask);
    al->get_block(ind, NULL, &len);
    _uint i = loc / bin_size;
    unsigned long val = (genes[i] >> off) & lmask;
    if (i + 1 < N) {
      val |= (genes[i + 1] & hmask) << (bin_size - off);
    }

    _uint j = 0;
    for (; j < (len/bin_size) - 1; ++j) {
      for (_uint k = 0; k < bin_size; ++k) {
        iss << ((genes[i+j] >> k) & 1);
      }
    }
    _uint mlen = len % bin_size;
    lmask = ((unsigned long)1 << mlen) - 1;
    if (mlen == bin_size) {
      lmask = ULONG_MAX;
    }
    val = genes[i+j] & lmask;
    for (_uint k = 0; k < mlen; ++k) {
      iss << ((val >> k) & 1);
    }
  }
#ifdef USE_CUSTOM_CONTAINERS
  return iss;
#else
  return iss.str();
#endif
}

}
