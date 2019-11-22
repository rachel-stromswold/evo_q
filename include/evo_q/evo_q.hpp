#ifndef GENETICS_H
#define GENETICS_H

#include <stdio.h>
#include <cstring>
#include <iostream>
#include <random>
#include <math.h>
#include <stdarg.h>
#include <climits>
#include <memory>

//#ifndef USE_CUSTOM_CONTAINERS
#include <vector>
#include <unordered_map>
#include <string>
#include <sstream>
//#endif

#define PI 3.1415926535897932

#define DEFAULT_STRING_SIZE 8
#define DEFAULT_PRINT_SIZE  24

#define CODE_MISC		5
#define CODE_ARG_INVALID	4
#define CODE_ARG_RANGE  	3
#define CODE_MATH_ERROR  	2
#define CODE_WARN		1
#define CODE_NONE		0

//UTIL_H
#define DEF_POP_SIZE            20
#define DEF_BREED_POP_SIZE      6
#define DEF_NUM_GENS            10
#define DEF_COUPLING_BITS       16
#define DEF_MAX_COUPLING        100
#define DEF_COUP_VAR            0.25
#define DEF_COUP_MEAN           0.0
#define DEF_MUTATE_PROB         0.1
#define DEF_HYPER_THRESH        0.8
#define DEF_NUM_CROSSOVERS      1

//flags
#define WAIT_CON		1
#define VERBOSE			2
#define NOISE_COMPENSATE	4
#define MULTIPLES_NONE		7
#define MULTIPLES_SKIP		8
#define MULTIPLES_PERTURB	16
#define MULTIPLES_AVG		24
#define ASYNC_EVAL		32

//selection types
#define SELECT_ROULETTE		0
#define SELECT_TOURNAMENT	1
#define SELECT_USE_REPLACE	3
#define SELECT_ROULETTE_POOL	4

//POPULATION_H
#define NUM_GENES	10
#define NUM_CARRY	1

#define MAX_NUM_GENS	100
#define OUT_BUF_SIZE    50

#define FLAG_NONE_SET	0
#define FLAG_STATS_SET	1
#define FLAG_DIST_SET	2
#define FLAG_BEST_FOUND	4
#define VALID_BEST	16
#define BUF_SIZE      50
#define DEF_SORT_PARAM    -3
#define FLAG_NONE_SET	0
#define FLAG_STATS_SET	1
#define FLAG_DIST_SET	2
#define FLAG_BEST_FOUND	4
#define FLAG_FRONTS	8
#if __cplusplus >= 201402L
#define DEPRECATED(msg) [[ deprecated(msg) ]]
#elif defined(__GNUC__)
#define DEPRECATED(msg) __attribute__ ((deprecated(msg)))
#elif defined(_MSC_VER)
#define DEPRECATED(msg) __declspec(deprecated(msg))
#else
#define DEPRECATED(msg) 
#endif

namespace Genetics {
//UTIL_H
    template <bool condition, typename T>
    struct enable_if_c {};

    template <typename T>
    struct enable_if_c<true, T> { typedef T type; };

    //error message if we don't have average_fitness with the correct type
    template<typename, typename T>
    struct has_average_fitness {
        static_assert(
            std::integral_constant<T, false>::value,
            "Second template parameter needs to be of function type.");
    };

    // specialization that does the checking
    template<typename T, typename Ret, typename... Args>
    struct has_average_fitness<T, Ret(Args...)> {
    private:
        template<typename U>
        static constexpr auto check(U*) ->
        typename std::is_same<
                decltype( std::declval<U>().average_fitness( std::declval<Args>()... ) ),
                Ret    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            >::type;   // attempt to call it and see if the return type is correct

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<T>(0)) type;

    public:
        static constexpr bool value = type::value;
    };
    
    typedef unsigned int _uint;
    typedef unsigned char _uchar;
    
    typedef double (*FIT_FUNC_PTR)(size_t, char*);
    
    typedef enum Type{t_int, t_uint, t_real, t_bitstream, t_terminator} Type;
    
    void error (int fatal, std::string msg, ...);
    
    unsigned long encodeGray (unsigned long number);
    unsigned long decodeGray (unsigned long code);
    
    size_t getlen_and_clean (char* str);
    
    unsigned int nChoosek( unsigned int n, unsigned int k );
    
    _uint factorial(_uint n);
    inline double max(double a, double b) {
      if (a > b) {
        return a;
      } else {
        return b;
      }
    }
    
    std::string read_number(std::string::iterator* it);
    
    //returns whether the array
    template <class T>
    bool contains(T* arr, size_t len, T val) {
        for (size_t i = 0; i < len; ++i) {
            if (arr[i] == val) {
                return true;
            }
        }
        return false;
    }
    
#ifdef USE_CUSTOM_CONTAINERS
    template<class T>
    class Vector {
    protected:
        T* buf;
        size_t buf_capacity;
        size_t buf_size;
        size_t checksum = 0;
        
        void erase_all() {
            for (size_t i = 0; i < buf_size; ++i) {
                buf[i].~T();
            }
        }
        
        void grow(size_t new_capacity = 0) {
            if (new_capacity <= buf_capacity) {
                if (buf_capacity == 0) { buf_capacity = 1; }
                new_capacity = buf_capacity*2;
            }
            T* new_buf = (T*)malloc(sizeof(T)*new_capacity);
            for (size_t i = 0; i < buf_size; ++i) {
                new( &(new_buf[i]) ) T(std::move(buf[i]));
            }
            free(buf);
            buf = new_buf;
            buf_capacity = new_capacity;
        }
        
    public:
        void resize(_uint new_size) {
            T* new_buf = (T*)malloc(sizeof(T)*new_size);
            for (_uint i = 0; i < buf_size && i < new_size; ++i) {
                new_buf[i] = buf[i];
            }
            for (_uint i = buf_size; i < new_size; ++i) {
                new( &(new_buf[i]) ) T();
            }
            for (_uint i = new_size; i < buf_size; ++i) {
                buf[i].~T();
            }
            free(buf);
            buf = new_buf;
            
            buf_capacity = new_size;
            buf_size = new_size;
        }
        void resize(_uint new_size, T& dat) {
            T* new_buf = (T*)malloc(sizeof(T)*new_size);
            for (_uint i = 0; i < buf_size && i < new_size; ++i) {
                new( &(new_buf[i]) ) T( std::move(buf[i]) );
            }
            for (_uint i = buf_size; i < new_size; ++i) {
                new( &(new_buf[i]) ) T(dat);
            }
            for (_uint i = new_size; i < buf_size; ++i) {
                buf[i].~T();
            }
            free(buf);
            buf = new_buf;
            
            buf_capacity = new_size;
            buf_size = new_size;
        }
        void resize_erase(_uint new_size) {
            for (_uint i = 0; i < buf_size; ++i) {
                buf[i].~T();
            }
            free(buf);
            buf = (T*)malloc(sizeof(T)*new_size);
            for (_uint i = 0; i < new_size; ++i) {
                new( &(buf[i]) ) T();
            }
            buf_size = new_size;
            buf_capacity = new_size;
        }
        void shrink_to_fit() {
            buf_capacity = buf_size;
            T* new_buf = (T*)malloc(sizeof(T)*buf_size);
            for (_uint i = 0; i < buf_size; ++i) {
                new( &(new_buf[i]) ) T( std::move(buf[i]) );
            }
            free(buf);
            buf = new_buf;
        }
        void set(T* dat, _uint n = 0) {
            if (buf_capacity < n) {
                resize_erase(n);
            }
            for (size_t i = 0; i < buf_size; ++i) {
                buf[i] = dat[i];
            }
            for (size_t i = buf_size; i < n; ++i) {
                new( &(buf[i]) ) T(dat[i]);
            }
            buf_size = n;
        }
        typedef T* iterator;
        Vector() {
            buf = (T*)malloc(sizeof(T)*DEFAULT_STRING_SIZE);
            buf_capacity = DEFAULT_STRING_SIZE;
            buf_size = 0;
        }
        Vector(_uint n, T* dat) {
            buf = (T*)malloc(sizeof(T)*n);
            buf_size = 0;
            buf_capacity = n;
            set(dat, n);
        }
        Vector(_uint n, const T& dat) {
            buf = (T*)malloc(sizeof(T)*n);
            buf_size = n;
            buf_capacity = n;
            for (_uint i = 0; i < n; ++i) {
                new( &(buf[i]) ) T(dat);
            }
        }
        Vector(const Vector<T>& obj) {
            if (this != &obj) {
                buf = (T*)malloc(sizeof(T)*(obj.buf_capacity));
                buf_size = obj.buf_size;
                buf_capacity = obj.buf_capacity;
                for (_uint i = 0; i < buf_size; ++i) {
                    new( &(buf[i]) ) T(obj.buf[i]);
                }
            }
        }
        Vector(Vector<T>&& obj) {
            if (this != &obj) {
                buf = obj.buf;
                buf_capacity = obj.buf_capacity;
                buf_size = obj.buf_size;
                obj.buf_capacity = 0;
                obj.buf_size = 0;
                obj.buf = NULL;
            }
        }
        Vector& operator=(Vector<T>&& obj) {
            if (this != &obj) {
                size_t tmp_capacity = buf_capacity;
                size_t tmp_size = buf_size;
                T* tmp_buf = buf;
                buf_capacity = obj.buf_capacity;
                buf_size = obj.buf_size;
                buf = obj.buf;
                obj.buf_capacity = tmp_capacity;
                obj.buf_size = tmp_size;
                obj.buf = tmp_buf;
            }
            return *this;
        }
        void swap(Vector<T>& obj) {
            Vector<T> tmp = std::move(obj);
            obj = std::move(*this);
            *this = std::move(tmp);
        }
        Vector<T>& operator=(Vector<T>& obj) {
            swap(obj);
            return *this;
        }
        T& operator[](size_t index) {
            if (index >= buf_size) {
                error(1, "Attempt to access invalid index %d.", index);
            }
            return buf[index];
        }
        ~Vector() {
            erase_all();
            free(buf);
            buf=NULL;
        }
        size_t capacity() {
            return buf_capacity;
        }
        size_t size() {
            return buf_size;
        }
        iterator begin() {
            return &(buf[0]);
        }
        iterator end() {
            return buf + buf_size;
        }
        void insert(iterator it, _uint n, T& tmplt) {
            if (buf_size + n >= buf_capacity) {
                grow(buf_size + 2*n);
            }
            size_t start_ind = it - begin();
            if (start_ind != buf_size) {
                for (size_t i = buf_size - 1; i >= start_ind; --i) {
                    new( &(buf[i+n]) ) T(std::move(buf[i]));
                }
            }
            for (size_t i = 0; i < n; ++i) {
                new( &(buf[i + start_ind]) ) T(tmplt);
            }
            buf_size += n;
        }
        void push_back(const T& tmplt) {
            if (buf_size >= buf_capacity) {
                grow(2*buf_capacity);
            }
            new( &(buf[buf_size]) ) T(tmplt);
            buf_size++;
        }
        T pop_back() {
            buf_size--;
            T ret(buf[buf_size]);
            buf[buf_size].~T();
            return ret;
        }
        
        T& back() {
            return buf[buf_size - 1];
        }
        
        T pop(int index) {
            T tmp = buf[index];
            for (int i = index; i < buf_size - 1; --i) {
                buf[i] = buf[i + 1];
            }
            buf_size--;
            buf[buf_size].~T();
            return tmp;
        }
        
        void erase_back() {
            buf_size--;
            buf[buf_size].~T();
        }
        T* contents() {
            return buf;
        }
        bool operator==(Vector<T>& obj) {
            if (buf_size != obj.size()) {
                return false;
            }
            for (size_t i = 0; i < buf_size; ++i) {
                if (obj[i] != buf[i]) {
                    return false;
                }
            }
            return true;
        }
    };
    
    class String : public Vector<char> {
    public:
        String() {
            Vector<char>::buf_size = 1;
            Vector<char>::buf[0] = 0;
        }
        void set(const char* dat) {
            Vector<char>::buf_size = 0;
            char t;
            do {
                t = dat[Vector<char>::buf_size];
                if (Vector<char>::buf_size == Vector<char>::buf_capacity) {
                    grow(Vector<char>::buf_capacity*2);
                }
                this->buf[Vector<char>::buf_size] = t;
                Vector<char>::buf_size++;
            } while (t != 0);
        }
        void append(const char* str) {
            //we need to remove the trailing 0
            erase_back();
            char t = str[0];
            _uint i = 0;
            while (t != 0) {
                t = str[i];
                push_back(t);
                i++;
            }
            buf[buf_capacity - 1] = 0;
        }
        String(const char* str) {
            buf_size = 1;
            append(str);
        }
        String& operator<<(int n) {
            char tmp[DEFAULT_PRINT_SIZE];
            snprintf(tmp, DEFAULT_PRINT_SIZE - 1, "%d", n);
            tmp[DEFAULT_PRINT_SIZE - 1] = 0;
            append(tmp);
            return *this;
        }
        String& operator<<(long n) {
            char tmp[DEFAULT_PRINT_SIZE];
            snprintf(tmp, DEFAULT_PRINT_SIZE - 1, "%ld", n);
            tmp[DEFAULT_PRINT_SIZE - 1] = 0;
            append(tmp);
            return *this;
        }
        String& operator<<(unsigned long n) {
            char tmp[DEFAULT_PRINT_SIZE];
            snprintf(tmp, DEFAULT_PRINT_SIZE - 1, "%lu", n);
            tmp[DEFAULT_PRINT_SIZE - 1] = 0;
            append(tmp);
            return *this;
        }
        String& operator<<(const char* str) {
            char tmp[DEFAULT_PRINT_SIZE];
            snprintf(tmp, DEFAULT_PRINT_SIZE - 1, "%s", str);
            tmp[DEFAULT_PRINT_SIZE - 1] = 0;
            append(tmp);
            return *this;
        }
        String& operator<<(_uint n) {
            char tmp[DEFAULT_PRINT_SIZE];
            snprintf(tmp, DEFAULT_PRINT_SIZE - 1, "%u", n);
            tmp[DEFAULT_PRINT_SIZE - 1] = 0;
            append(tmp);
            return *this;
        }
        String& operator<<(double x) {
            char tmp[DEFAULT_PRINT_SIZE];
            snprintf(tmp, DEFAULT_PRINT_SIZE - 1, "%f", x);
            tmp[DEFAULT_PRINT_SIZE - 1] = 0;
            append(tmp);
            return *this;
        }
        String& operator<<(char* str) {
            append(str);
            return *this;
        }
        String& operator<<(char c) {
            buf[Vector<char>::buf_size - 1] = c;
            push_back(0);
            return *this;
        }
        String& operator=(const char* str) {
            set(str);
            return *this;
        }
        char* c_str() {
            return Vector<char>::buf;
        }
    };
#else
    template<typename T>
    using Vector = std::vector<T>;
    typedef std::string String;
#endif

    template<typename T>
    bool contains(Vector<T>& vec, T target) {
      for (auto it = vec.begin(); it != vec.end(); ++it) {
        if (*it == target) {
          return true;
        }
      }
      return false;
    }

#ifndef USE_EXCEPTIONS
    bool has_error();
    String get_error();
#endif

int divideup(int numerator, int denominator);
    
//draw k elements from the integer range from 0 to n useful for sampling from arrays
class SampleDraw {
private:
  _uint n_;
  _uint k_;
  std::uniform_int_distribution<_uint> dist;
  bool replace;

public:
  SampleDraw(_uint p_n, _uint p_k, bool replace=false);

  void reset() { dist.reset(); }
  _uint n() { return n_; }
  _uint k() { return k_; }

  template <class Generator>
  Vector<_uint> operator()(Generator& g) {
    Vector<_uint> ret(k_, 0);
    for (_uint i = 0; i < k_; ++i) {
      if (replace) {
        ret[i] = dist( g) % n_;
      } else {
        //TODO: this method biases results since UINT_MAX is not necessarily divisible by n-i
        _uint x = dist( g ) % (n_ - i);
        bool append = true;
        for (_uint j = 0; j < i; ++j) {
          if (x >= ret[j]) {
            ++x;
          } else {
            //move everything up to maintain a proper sorting
            for (_uint l = i; l > j; --l) {
              ret[l] = ret[l - 1];
            }
            ret[j] = x;
            append = false;
            break;
          }
        }
        if (append) {
          ret[i] = x;
        }
      }
    }
    return ret;
  }
};

//rearrange n elements into a random order
class Shuffle {
private:
  _uint n_;
  std::uniform_int_distribution<_uint> dist;
  bool replace;

public:
  Shuffle(_uint p_n);

  void reset() { dist.reset(); }
  _uint n() { return n_; }

  template <class Generator>
  Vector<_uint> operator()(Generator& g) {
    Vector<_uint> ret(n_, 0);
    for (_uint i = 0; i < n_; ++i) {
      _uint x = dist( g ) % (n_ - i);
      bool contained = true;
      while (contained) {
        contained = false;
        for (_uint j = 0; j < i; ++j) {
          if (ret[j] == x) { contained = true;break; }
        }
        if (!contained) {
          ret[i] = x;
        } else {
          ++x;
        }
      }
    }
    return ret;
  }
};

//PHENOTYPE_H
    
    struct VarContainer {
        _uint loc;
        double range_hi;
        double range_lo;
        double factor;
        Type type;
        VarContainer(_uint p_loc, double p_lo, double p_hi, Type p_type) {loc = p_loc;range_lo = p_lo;range_hi = p_hi;type = p_type;}
        VarContainer() { loc = 0;range_hi = 0; range_lo = 0; factor = 0; type = t_real; }
    };
    
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
    
//PARSE_H
    
class ArgStore {
  private:
    std::mt19937 generator;
//REMOVAL CANDIDATE
    std::binomial_distribution<unsigned char>* long_bin;
    std::binomial_distribution<unsigned char>* short_bin;
    std::bernoulli_distribution* bern_mut;
    std::bernoulli_distribution* bern_cross;
    unsigned int short_bin_n = 0;
//END REMOVAL CANDIDATE

    unsigned flags;
    size_t pop_size;
    size_t breed_pop_size;
    size_t num_gens;
    int num_crossovers;
    double init_coup_var;
//REMOVAL_CANDIDATE
    double init_coup_mean;
//END REMOVAL CANDIDATE
    double crossover_prob;
    double mutate_prob;
    double hypermutation_threshold;
    double replacement_fraction;
    String out_fname;

    bool activate = true;
    bool async_evaluation = true;
    _uchar selection_type = SELECT_ROULETTE;
    int seed = 0;
    _uint noise_compensation_runs = 0;
    std::unordered_map<String, String> custom_parameters;

  public: 
    double forget_weight = 0;

    ArgStore();
    ArgStore(const ArgStore& o);
    ArgStore(ArgStore&& o);
    ~ArgStore();

    void initialize_from_args(size_t argc, char** argv);
    void initialize_from_file(const char* fname);
    void initialize();
    void read_file(std::string fname);
    std::mt19937& get_generator() { return generator; }
//REMOVAL CANDIDATE
    unsigned int sample_binomial(unsigned int n);
    bool random_mutation();
    bool random_crossover();
    void print_data();
//END REMOVAL CANDIDATE

    bool skip_multiples() { return flags & MULTIPLES_SKIP; }
    bool average_multiples() { return flags & MULTIPLES_AVG; }
    bool perturb_multiples() { return flags & MULTIPLES_PERTURB; }
    _uint noise_compensate() { return noise_compensation_runs; }
    void set_noise_compensation(_uint val);

    void set_selection_type(_uchar val);
    //bool use_roulette() { return selection_type == SELECT_ROULETTE; }
    //bool use_tournament() { return selection_type & SELECT_TOURNAMENT; }
    //bool use_roulette_pool() { return selection_type == SELECT_ROULETTE_POOL; }
    //bool tournament_replacement() { return selection_type & SELECT_USE_REPLACE; }
    bool async() { return async_evaluation; }
    void set_async(bool val) { async_evaluation = val; }

    size_t get_pop_size()			{ return pop_size; }
    void set_pop_size(size_t n)			{ pop_size = n; }
    //size_t get_survivors()			{ return breed_pop_size; }
    //void set_survivors(size_t n)		{ breed_pop_size = n; }
    size_t get_num_gens() 			{ return num_gens; }
    void set_num_gens(size_t n) 		{ num_gens = n; }
    int get_num_crossovers() 			{ return num_crossovers; }
    void set_num_crossovers(int n) 		{ num_crossovers = n; }
    double get_init_coup_var()			{ return init_coup_var; }
    void set_init_coup_var(double x)		{ init_coup_var = x; }
    double get_init_coup_mean() 		{ return init_coup_mean; }
    void set_init_coup_mean(double x) 		{ init_coup_mean = x; }
    double get_crossover_prob()			{ return crossover_prob; }
    void set_crossover_prob(double x)		{ crossover_prob = x; }
    double get_mutate_prob()			{ return mutate_prob; }
    void set_mutate_prob(double x)		{ mutate_prob = x; }
    double get_hypermutation_threshold()	{ return hypermutation_threshold; }
    void set_hypermutation_threshold(double x)	{ hypermutation_threshold = x; }
    double get_replacement_fraction()		{ return replacement_fraction; }
    void set_replacement_fraction(double x)	{ replacement_fraction = x; }
    bool wait_for_con()				{ return flags & WAIT_CON; }
    void set_wait_for_con(bool b = true)	{ flags |= WAIT_CON * ((_uint)b); }
    bool verbose()				{ return flags & VERBOSE; }
    void set_verbose(bool b = true)		{ flags |= VERBOSE * ((_uint)b); }
    String get_custom_parameter(String val);
    int read_custom_int(String val, int default_val);
    double read_custom_double(String val, double default_val);
    String get_out_fname() 			{ return out_fname; }
};
    
//CHROMOSOME_H
    
class Chromosome {
private:
  _uint N_BITS;
  _uint N_BYTES;
  size_t N;
  double* real_vals = NULL;
  size_t real_vals_size = 0;

protected:
  static const _uint bin_size = sizeof(unsigned long)*8;
//  unsigned long genes[(N_BYTES+sizeof(unsigned long)-1)/sizeof(unsigned long)];
  unsigned long* genes = NULL;
  size_t getBitStream (size_t n, size_t k, size_t x);
  _uchar use_real = 0;
  size_t get_real_vals_size() { return real_vals_size; }

public:
  Chromosome(_uint pn_bits);
  Chromosome(_uint pn_bits, _uchar real_mode);
  Chromosome(_uint pn_bits, Chromosome* o);
  Chromosome(Chromosome& other);
  Chromosome(Chromosome&& other);
  ~Chromosome();
  void exchange(Chromosome* other, size_t k);
  void exchange_uniform(ArgStore& args, Chromosome* other);

  unsigned int get_N() { return N; }
  unsigned int get_n_bits() { return N_BITS; }
  void swap(Chromosome& other);
  Chromosome& operator=(Chromosome& other);
  Chromosome& operator=(Chromosome&& other);
  bool operator==(Chromosome& other);

  void reset();

  unsigned char operator[](unsigned int i);
  //randomly mutate each bit in the gene
  bool real_space_mutate(ArgStore& args);
  void mutate(ArgStore& args);
  void slow_mutate(ArgStore& args);
  //set the gene to a new completely random value
  void randomize(PhenotypeMap* al, ArgStore& args);
  //sets the gene to encode the value specified by min, max
  void set_to_num(PhenotypeMap* al, _uint ind, double value);
  void set_to_int(PhenotypeMap* al, _uint ind, int value);
  void set_to_ulong(PhenotypeMap* al, _uint ind, unsigned long value);
  //returns the corresponding integer for the gene
  unsigned long gene_to_ulong(PhenotypeMap* al, _uint ind);
  int gene_to_int(PhenotypeMap* al, _uint ind);
  //returns a double value corresponding to the gene, it will have a value between max and min
  double gene_to_num(PhenotypeMap* al, _uint ind);
  String get_string(PhenotypeMap* al, _uint ind);
};
    
//ORGANISM_H
    
/**
 * \brief	An implementation of Fitness designed for single-objective optimization in the noise-free case. Other single objective problems can use Fitness classes that inherit from SingleFitness for generation of comparisons.
 */
class Fitness {
protected:
  double fitness;
public:
  double distance = 0;
  Fitness() { fitness = 0; }

  void reset() {}
  //For this simple implementation of a fitness tracker, set_fitness() and update() behave identically. However, other trackers, such as NoisyFitness keep track of an average fitness which may behave differently.
  void set_fitness(double val, _uint i = 0) { fitness = val; }
  void update(double val, _uint i = 0) { fitness = val; }
  _uint get_n_objs() { return 1; }
  double get_fitness(_uint i = 0) { return fitness; }
  double get_cost(_uint i = 0) { return -get_fitness(i); }
  double get_uncertainty(_uint i = 0) { return 0.0; }
};
typedef Fitness SingleFitness;

/**
 * \brief	The class MultiFitness is used by implementations of the Selector abstract class to select which organism is more fit when multiple objectives are to be considered
 */
class MultiFitness : public Fitness {
protected:
  _uint N_OBJS;
  Vector<double> fitness;
 
public:
  _uint n_dominations = 0;
  _uint rank = 0;

  MultiFitness(_uint pn_objs = 1) : fitness(pn_objs, 0.0) { N_OBJS = pn_objs; }

  void reset() {}
  void set_fitness(double val, _uint i = 0);
  void update(double val, _uint i);
  _uint get_n_objs() { return N_OBJS; }
  double get_fitness(_uint i);
  double get_cost(_uint i) { return -get_fitness(i); }
  double get_uncertainty(_uint i) { return 0.0; }
  _uint get_rank() { return rank; }
  double get_distance() { return distance; }
  _uint get_n_dominations() { return n_dominations; }
};

/**
 * \brief	An implementation of Fitness designed for single-objective optimization in the noisy case
 */
class NoisyFitness : public SingleFitness {
protected:
  double variance = 0;
  _uint n_evaluations = 0;

public:
  _uint get_n_evaluations() { return n_evaluations; }
  double get_uncertainty(_uint i = 0) { return sqrt(variance); }
  void update(double val, _uint i = 0);
  void average_fitness(NoisyFitness& other);
};

/**
 * \brief	An implementation of Fitness designed for single-objective optimization in the noisy case. This implementation uses a parameter forget_weight that biases results to more heavily weigh recent observations
 */
class NoisyFitnessForgetful : public NoisyFitness {
protected:
  double fitness = 0, variance = 0;
  double forget_weight;
  bool evaluated = false;
public:
  NoisyFitnessForgetful(double p_forget_weight = 1.0);
  void update(double val, _uint i = 0);
  void average_fitness(NoisyFitnessForgetful& other);
};

/**
 * \brief	An implementation of Fitness designed for single-objective optimization in the noisy case
 */
class NoisyMultiFitness : public MultiFitness {
protected:
  Vector<double> variances;
  _uint n_evaluations;

public:
  NoisyMultiFitness(_uint pn_objs = 1);
  double get_fitness(_uint i);
  void update(double val, _uint i);
  double get_uncertainty(_uint i) { return sqrt(variances[i]); }
};

template <class FitType>
class Organism;

template <class FitType>
class Problem {
public:
  _uint N_BITS, N_PARAMS, N_OBJS;
  std::shared_ptr<PhenotypeMap> map;
  //Vector<Result> result_list;
//  Problem() {std::cout << "Initializing problem...\n"; }
  Problem(unsigned n_bits, unsigned n_params, int n_objs) : N_BITS(n_bits), N_PARAMS(n_params), N_OBJS(n_objs), map(std::make_shared<PhenotypeMap>(n_bits)) {}

  virtual void evaluate_fitness(Organism<FitType>* org) = 0;
#ifdef USE_LIBOMP
  void evaluate_fitness_async(Organism<FitType>* org, _uint i = 0) { evaluate_fitess(org); }
#endif
};

template <class FitType>
class Organism {
  static_assert( std::is_base_of<Fitness,FitType>::value, "FitType must be derived from SingleFitness or MultiFitness" );
private:
  typedef Organism<FitType> MyType;
  _uint N_BITS;
  _uint N_OBJS;

  char output_stream[BUF_SIZE];
  double penalty = 0.0;
  size_t output_len;
  std::shared_ptr<PhenotypeMap> al;

protected:
  Chromosome genes;
  size_t n_nodes;
  FitType fit;
  template <typename T = FitType>
  auto set_fit_n_objs(_uint n_objs) -> typename enable_if_c<std::is_base_of<MultiFitness,T>::value, void>::type {
    FitType tempFit(n_objs);
    fit = tempFit;
  }
  template <typename T = FitType>
  auto set_fit_n_objs(_uint n_objs) -> typename enable_if_c<!std::is_base_of<MultiFitness,T>::value, void>::type {}
  //Vector<double> fitness;
  //Vector<double> fit_vars;
  //int n_evaluations = 0;

public:
  //this is not used internally, but can be set when evaluating the fitness
  String misc_data;
  double coupling_range;
  double coupling_prec; 

  Organism() : genes(0) { N_BITS = 0;N_OBJS = 0; }
  Organism(int pn_bits, int pn_objs, PhenotypeMap* p_al) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs),
  genes(N_BITS),
  al(p_al)
  {
    set_fit_n_objs(pn_objs);
    memset(output_stream, 0, BUF_SIZE);
    reset_fitness();
  }
  Organism(int pn_bits, int pn_objs, Chromosome p_genes, PhenotypeMap* p_al) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs),
  genes(p_genes),
  al(p_al)
  {
    set_fit_n_objs(pn_objs);
    reset_fitness();
  }
  Organism(int pn_bits, int pn_objs, std::shared_ptr<PhenotypeMap> p_al) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs),
  genes(N_BITS),
  al(p_al)
  {
    set_fit_n_objs(pn_objs);
    memset(output_stream, 0, BUF_SIZE);
    reset_fitness();
  }
  Organism(int pn_bits, int pn_objs, Chromosome p_genes, std::shared_ptr<PhenotypeMap> p_al) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs),
  genes(p_genes),
  al(p_al)
  {
    set_fit_n_objs(pn_objs);
    reset_fitness();
  }
  /**
   * \brief	Update the Fitness object to match the template temp. This should be used for parameter setting before any evaluation calls have been made.
   */
  void set_fitness_stats(FitType temp) { fit = temp; }
  //Organism(const Organism &obj);
  //Organism(Organism&& obj);
  //~Organism();
  Organism copy() {
    MyType ret(N_BITS, N_OBJS, genes, al);
    ret.fit = fit;
    
    return ret;
  }

  bool operator==(Organism& obj) {
    for (_uint i = 0; i < al->get_num_params(); ++i) {
      Type t = al->get_type(i);
      if (t == t_real) {
        if (read_real(i) != obj.read_real(i)) {
          return false;
        }
      } else if (t == t_int) {
        if (read_int(i) != obj.read_int(i)) {
          return false;
        }
      } else {
        if (read_uint(i) != obj.read_uint(i)) {
          return false;
        }
      }
    }
    return true;
  }
  bool operator!=(Organism& obj) {
    for (_uint i = 0; i < al->get_num_params(); ++i) {
      Type t = al->get_type(i);
      if (t == t_real) {
        if (read_real(i) != obj.read_real(i)) {
          return true;
        }
      } else if (t == t_int) {
        if (read_int(i) != obj.read_int(i)) {
          return true;
        }
      } else {
        if (read_uint(i) != obj.read_uint(i)) {
          return true;
        }
      }
    }
    return false;
  }
  //bool operator>(Organism& obj) { return !obj.dominates(this); }
  //bool operator<(Organism& obj) { return obj.dominates(this); }

  void swap(Organism& obj) {
    FitType tmp = fit;
    fit = obj.fit;
    obj.fit = tmp;
  }
  bool valid() { return (al != NULL && N_OBJS > 0 && N_BITS > 0); }

  typedef std::shared_ptr< MyType > OrgPtr;
  std::pair<OrgPtr, OrgPtr> breed(ArgStore& args, OrgPtr o) {
    if (get_n_bits() != o->get_n_bits()) {
      error(CODE_MISC, "Cannot breed organsims with a differing number of bits, %d and %d.", get_n_bits(), o->get_n_bits());
    }
    std::pair<OrgPtr, OrgPtr> children;

    memset(output_stream, 0, BUF_SIZE);
    Chromosome gene0(genes);
    Chromosome gene1(o->genes);

    if (args.random_crossover()) {
      if (args.get_num_crossovers() <= 0) {
        gene0.exchange_uniform(args, &gene1);
      } else {
        std::uniform_int_distribution<size_t> rint( 0, gene0.get_n_bits() - 1 );
        for (int n = 0; n < args.get_num_crossovers(); ++n) {
          size_t exch_bit = rint( args.get_generator() );
          gene0.exchange(&gene1, exch_bit);
        }
      }
      children.first  = std::make_shared<MyType>(N_BITS, N_OBJS, gene0, al);
      children.second = std::make_shared<MyType>(N_BITS, N_OBJS, gene1, al);
    } else {
      children.first  = std::make_shared<MyType>(*this);
      children.second = std::make_shared<MyType>(*o);
    }
    children.first->mutate(args);
    children.second->mutate(args);

    return children;
  }
  void mutate(ArgStore& args) {
#ifdef MUT_SLOW
    genes.slow_mutate(args);
#else
    genes.mutate(args);
#endif
  }
  void reset_fitness() {
    fit.reset();
    misc_data = "";
  }
  void randomize(ArgStore& args) { genes.randomize(al.get(), args); }
  void randomize(ArgStore& args, Organism* orgtmp) {
    //make most of the genes similar, with one gene more wildly varied
    std::uniform_int_distribution<size_t> chrom(0, al->get_num_params() - 1);
    size_t high_ind = chrom( args.get_generator() );

    double var = args.get_init_coup_var();
    double lvar = var/al->get_num_params();
    genes.reset();
    for (size_t i = 0; i < al->get_num_params(); i++) {
      Type t = al->get_type(i);
      if (i != high_ind) {
        double mean;
        if (t == t_real) {
          mean = orgtmp->read_real(i);
          double range = al->get_range_max(i) - al->get_range_min(i);
          std::normal_distribution<double> norm(mean, lvar*range);
          double x = norm( args.get_generator() );
          genes.set_to_num(al.get(), i, x);
        } else {
          _uint max_possible = 1 << al->get_block_length(i);
          mean = (double)(max_possible - genes.gene_to_int(al.get(), i))/2;
          //scale lvar to lvar*max_possible/2 and set n and p to produce the according mean and variance
          double p = 1 - lvar*max_possible/(2*mean);
          int n = (int)mean/p;
          std::binomial_distribution<int> dist(n, p);
          int x = dist( args.get_generator() )*2 + genes.gene_to_int(al.get(), i);
          genes.set_to_int(al.get(), i, x);
        }
      }
    }
    if (al->get_type(high_ind) == t_real) {
      //set the genome representation for the high variance index
      double mean = orgtmp->read_real(high_ind);
      double range = al->get_range_max(high_ind) - al->get_range_min(high_ind);
      std::normal_distribution<double> norm(mean, var*range);
      double x = norm( args.get_generator() );
      genes.set_to_num(al.get(), high_ind, x);
    } else {
      int tmpx = orgtmp->read_int(high_ind);
      _uint max_possible = 1 << al->get_block_length(high_ind);
      std::normal_distribution<double> dist((double)tmpx, var*max_possible);
      int x = (int)dist( args.get_generator() );
      genes.set_to_int(al.get(), high_ind, x);
      std::cout << " max_possible = " << max_possible << " x = " << x << " orgtmp_x = " << tmpx << "\n";
    }
  }

  //void evaluate_fitness_noisy(Problem<FitType>* prob, double forget_weight=0);
  void evaluate_fitness(Problem<FitType>* prob) { penalty = 0;prob->evaluate_fitness(this); }

  FitType& get_fitness_info() { return fit; }

  double get_fitness(_uint i = 0) { return fit.get_fitness(i); }
  double get_cost(_uint i = 0) { return fit.get_cost(i); }
  void apply_penalty(double val) { penalty = val; }

  double get_penalty() { return penalty; }
  bool penalized() { return penalty != 0; }
  //set fitness functions with arguments
  void set_fitness(_uint i, double val) {
    if (i >= fit.get_n_objs()) {
      error(CODE_ARG_RANGE, "Attempt to modify invalid fitness index %d, size is %d.", i, fit.get_n_objs());
    } else {
      fit.set_fitness(val, i);
    }
  }
  void update(_uint i, double val) {
    if (i >= fit.get_n_objs()) {
      error(CODE_ARG_RANGE, "Attempt to modify invalid fitness index %d, size is %d.", i, fit.get_n_objs());
    } else {
      fit.update(val, i);
    }
  }
  void set_cost(_uint i, double val) { set_fitness(i, -val); }
  void update_cost(_uint i, double val) { update(i, -val); }
  //set fitness functions without arguments
  void set_fitness(double val) { fit.set_fitness(val, 0); }
  void set_cost(double val) { fit.set_fitness(-val, 0); }
  void update(double val) { fit.update(val, 0); }
  void update_cost(double val) { fit.update(-val, 0); }

  void set_int(_uint i, int value) { genes.set_to_int(al.get(), i, value); }
  void set_uint(_uint i, _uint value) { genes.set_to_ulong(al.get(), i, value); }
  void set_real(_uint i, double value) { genes.set_to_num(al.get(), i, value); }
  double read_real(_uint i) { return genes.gene_to_num(al.get(), i); }
  int read_int(_uint i) { return genes.gene_to_int(al.get(), i); }
  _uint read_uint(_uint i) { return genes.gene_to_ulong(al.get(), i); }
  String get_chromosome_string(_uint i) {
    if (N_BITS == 0 || !al) {
      error(1, "Attempt to access string for uninitialized organism.");
    }
    if ( i >= al->get_num_params() ) {
      error(1, "Attempt to access invalid parameter with index %d.", i);
    }
    return genes.get_string(al.get(), i);
  }
  char* get_output_stream() { return output_stream; }
  size_t get_output_len() {return output_len; }
  _uint get_n_bits() { return N_BITS; }
  _uint get_n_params() { return al->get_num_params(); }
  _uint get_n_objs() { return N_OBJS; }
};

//COMPARATOR_H

struct FitnessStats {
  double mean;
  double var;
  double max;
  double min;
};

template <typename FitType>
class Comparator {
public:
  static int compare(std::shared_ptr< Organism<FitType> > a, std::shared_ptr< Organism<FitType> > b) {
    FitType a_info = a->get_fitness_info();
    FitType b_info = a->get_fitness_info();
    if ( a->get_n_objs() != b->get_n_objs() ) {
      error(CODE_MISC, "Comparison of fitness values with a different number of objectives.");
    }
    int ret = 0;
    for (_uint i = 0; i < a_info.get_n_objs(); ++i) {
      if (a->get_fitness(i) < b->get_fitness(i)) {
        --ret;
      } else if (a->get_fitness(i) > b->get_fitness(i)) {
        ++ret;
      }
    }
    return ret;
  }
};

template <typename FitType>
class NSGAII_Comparator : public Comparator<FitType> {
public:
  static_assert( std::is_base_of<MultiFitness, FitType>::value, "FitType must be derived from MultiFitness for NSGAII comparator" );
  static int compare(std::shared_ptr< Organism<FitType> > a, std::shared_ptr< Organism<FitType> > b) {
    for (_uint i = 0; i < a->get_fitness_info().get_n_objs(); ++i) {
      if (a->get_fitness(i) <= b->get_fitness(i)) {
        return 0;
      }
    }
    return 1;
  }
};

typedef std::pair<_uint, _uint> ParentIndSet;

template <typename FitType, typename Comp=Comparator<FitType>>
class Selector {
protected:
  int partition(_uint fit_ind, std::vector<std::shared_ptr< Organism<FitType> >>& work_arr, int s, int e) {
    //double p = (*work_arr)[e]->get_fitness(fit_ind);
    int i = s;
    for (int j = s; j < e; ++j) {
      if ( Comp::compare( work_arr[j], work_arr[e] ) > 0) {
        std::shared_ptr<Organism<FitType>> tmp = work_arr[i];
        work_arr[i] = work_arr[j];
        work_arr[j] = tmp;
        ++i;
      }
    }
    if (i != e) {
      std::shared_ptr<Organism<FitType>> tmp = work_arr[i];
      work_arr[i] = work_arr[e];
      work_arr[e] = tmp;
    }
    return i;
  }
public:
  static const bool use_offspring=false;
  void sort_orgs(unsigned int fit_ind, std::vector<std::shared_ptr< Organism<FitType> >>& work_arr, int s = DEF_SORT_PARAM, int e = DEF_SORT_PARAM) {
    if (s == DEF_SORT_PARAM || e == DEF_SORT_PARAM) {
      sort_orgs(fit_ind, work_arr, 0, work_arr.size() - 1);
    } else if (s < e) {
      int p = partition(fit_ind, work_arr, s, e);
      sort_orgs(fit_ind, work_arr, s, p-1);
      sort_orgs(fit_ind, work_arr, p+1, e);
    }
  }
  Selector() {
    static_assert( std::is_base_of<Fitness, FitType>::value, "FitType must be derived from Fitness" );
  }
  //virtual OrganismPair select(Population& pop);
  virtual ~Selector() = default;
  //returns >=1 if a > b, 0 if a = b and -1 if a <= b
  //errors if the two fitness values are not comprable (they do not have the same number of objectives)
  
  virtual Vector<ParentIndSet> select(ArgStore& args, Vector<std::shared_ptr<Organism<FitType>>>& old_gen, Vector<std::shared_ptr<Organism<FitType>>>& offspring) = 0;

  //return the index of the best organism found in the population
  static _uint find_best_organism(Vector<std::shared_ptr<Organism<FitType>>>& orgs, Vector<FitnessStats>& pop_stats) {
    _uint ret = 0;
    if (orgs.size() > 0) {
      for (int j = 0; j < pop_stats.size(); ++j) {
        pop_stats[j].max = orgs[0]->get_fitness(j);
        pop_stats[j].min = orgs[0]->get_fitness(j);
        pop_stats[j].mean = orgs[0]->get_fitness(j) / orgs.size();
        for (size_t i = 0; i < orgs.size(); ++i) {
          double fitness_i = orgs[i]->get_fitness(j);
          if (fitness_i > pop_stats[j].max) {
            ret = i;
            pop_stats[j].max = fitness_i;
          }
          if (fitness_i < pop_stats[j].min) {
            pop_stats[j].min = fitness_i;
          }
          pop_stats[j].mean += fitness_i / orgs.size();
        }
        //calculate the variance
        pop_stats[j].var = 0;
        for (size_t i = 0; i < orgs.size(); ++i) {
          double fitness_i = orgs[i]->get_fitness(j);
          pop_stats[j].var += (fitness_i - pop_stats[j].mean)*(fitness_i - pop_stats[j].mean);
        }
        if (orgs.size() > 1) {
          pop_stats[j].var /= (orgs.size() - 1);
        }
      }
    }
    
    return ret;
  }
};

template <class FitType, typename MyComp=Comparator<FitType>>
class TournamentSelector : public Selector<FitType, MyComp> {
public:
  typedef MyComp Comp;
  Vector<ParentIndSet> select(ArgStore& args, Vector<std::shared_ptr<Organism<FitType>>>& old_gen, Vector<std::shared_ptr<Organism<FitType>>>& offspring) {
    _uint arena_size = args.read_custom_double("arena_size", 2);
    if (arena_size < 2) { arena_size = 2; }
    bool tournament_replacement = (args.get_custom_parameter("tournament_replacement") != "");
    _uint offspring_num = old_gen.size();
    SampleDraw sampler(offspring_num, arena_size, tournament_replacement);
    std::uniform_int_distribution<_uint> selector(0, offspring_num - 1);
    std::vector<Organism<FitType>*> children;
    Vector<ParentIndSet> ret( divideup(offspring_num, 2) );

    for (size_t i = 0; 2*i + 1 < offspring_num; ++i) {
      std::shared_ptr< Organism<FitType> > first_parent, second_parent;
      std::vector<_uint> t1 = sampler( args.get_generator() );
      std::vector<_uint> t2 = sampler( args.get_generator() );
      ret[i].first  = t1[0];
      ret[i].second = t2[0];
      for (size_t j = 1; j < t1.size(); ++j) {
        //check whether the fitness is an improvement and use variance as a tiebreaker
        if (Comp::compare(old_gen[t1[j]], old_gen[ret[i].first]) > 0) {
          ret[i].first  = t1[j];
        }
        if (Comp::compare(old_gen[t2[j]], old_gen[ret[i].second]) > 0) {
          ret[i].second = t2[j];
        }
      }
      //ensure that we don't use the same parent twice with reasonable probability
      if (ret[i].first == ret[i].second) {
        ret[i].second = selector( args.get_generator() );
      }
    }
    return ret;
  }
};

template <class FitType, typename MyComp=Comparator<FitType>>
class SurvivalSelector : public Selector<FitType, MyComp> {
public:
  typedef MyComp Comp;
  Vector<ParentIndSet> select(ArgStore& args, Vector<std::shared_ptr<Organism<FitType>>>& old_gen, Vector<std::shared_ptr<Organism<FitType>>>& offspring) {
    sort_orgs(0, old_gen);
    /*TODO: determine whether we actually need to figure out a way to keep this in place
     * //if all the organisms have the same fitness then reinitialize the population
    if (old_gen[0]->get_fitness(0) - old_gen[offspring_num-1]->get_fitness(0) < 0.001 ) {
      return true;
    }*/
    size_t offspring_num = old_gen.size();
    double min_fit = old_gen[offspring_num-1]->get_fitness(0);
    double total_fit = 0;
    for (size_t i = 0; i < offspring_num; ++i) {
      total_fit += old_gen[i]->get_fitness(0) - min_fit;
    }
    std::uniform_real_distribution<double> dist(0, total_fit);

    _uint survivors_num = args.read_custom_double("num_survivors", offspring_num/2);
    std::vector< std::shared_ptr<Organism<FitType>> > survivors(survivors_num);
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
      std::cerr << "Added orgs to survs:";
      for (size_t i = 0; i < survivors_num; ++i) {
        std::cerr << " " << banned[i];
      }
      std::cerr << std::endl;
    }
    free(banned);

    std::uniform_int_distribution<size_t> dist_surv0(0, survivors_num - 1);
    std::uniform_int_distribution<size_t> dist_surv1(0, survivors_num - 2);
    Vector<ParentIndSet> ret( divideup(offspring_num, 2) );
    for (_uint i = 0; 2*i + 1< offspring_num; ++i) {
      ret[i].first  = dist_surv0( args.get_generator() );
      ret[i].second = dist_surv1( args.get_generator() );
      if (ret[i].second >= ret[i].first) {
        ++ret[i].second;
      }
    }
    return ret;
  }
};

template <class FitType, typename MyComp=Comparator<FitType>>
class DominanceTournamentSelector : public TournamentSelector<FitType, MyComp> {
public:
  typedef MyComp Comp;

  //return the index of the best organism found in the population
  static _uint find_best_organism(Vector<std::shared_ptr<Organism<FitType>>>& orgs, Vector<FitnessStats>& pop_stats) {
    _uint ret = 0;
    if (orgs.size() > 0) {
      //set fitness to be given by the number of dominations
      Vector<int> n_dominations(orgs.size(), 0);
      for (_uint i = 0; i < orgs.size(); ++i) {
        for (_uint j = i + 1; j < orgs.size(); ++j) {
          int comp_val = MyComp::compare(orgs[i], orgs[j]);
          n_dominations[i] += comp_val;
          n_dominations[j] -= comp_val;
        }
        orgs[i]->set_fitness(0, n_dominations[i]);
      }

      pop_stats[0].max = n_dominations[0];
      pop_stats[0].min = n_dominations[0];
      pop_stats[0].mean = n_dominations[0] / orgs.size();
      for (size_t i = 0; i < orgs.size(); ++i) {
        double fitness_i = n_dominations[i];
        if (fitness_i > pop_stats[0].max) {
          ret = i;
        }
        if (fitness_i < pop_stats[0].min) {
          pop_stats[0].min = fitness_i;
        }
        pop_stats[0].mean += fitness_i / orgs.size();
      }
      //calculate the variance
      pop_stats[0].var = 0;
      for (size_t i = 0; i < orgs.size(); ++i) {
        double fitness_i = n_dominations[i];
        pop_stats[0].var += (fitness_i - pop_stats[0].mean)*(fitness_i - pop_stats[0].mean);
      }
      if (orgs.size() > 1) {
        pop_stats[0].var /= (orgs.size() - 1);
      }
    }
    
    return ret;
  }
};

template <typename FitType>
class NSGAII_TournamentSelector : public Selector<FitType, NSGAII_Comparator<FitType>> {
public:
  typedef std::shared_ptr< Organism<FitType> > OrgPtr;
  typedef NSGAII_Comparator<FitType> Comp;
  static const bool use_offspring=true;
private:  
  _uint n_objs = 1;
  Vector<Vector< std::shared_ptr<Organism<FitType>> >> pareto_fronts;
  Vector<FitnessStats> pop_stats;
  _uchar calculated_flags = 0;
  size_t min_penalty_ind, max_penalty_ind;
  void make_fronts(std::vector<OrgPtr>& cmb_arr) {
    pareto_fronts.clear();
    std::vector<OrgPtr> empty;
    pareto_fronts.push_back(empty);

    for (size_t i = 0; i < cmb_arr.size(); ++i) {
      cmb_arr[i]->get_fitness_info().n_dominations = 0;
      for (_uint j = 0; j < n_objs; ++j) {
        if (cmb_arr[i]->get_fitness(j) > pop_stats[j].max) {
          pop_stats[j].max = cmb_arr[i]->get_fitness(j);
        }
        if (cmb_arr[i]->get_fitness(j) < pop_stats[j].min) {
          pop_stats[j].min = cmb_arr[i]->get_fitness(j);
        }
      }
      for (size_t j = 0; j < cmb_arr.size(); ++j) {
        if (i != j) {
          //if the ith solution dominates the jth add the jth entry to the list of dominated solutions, otherwise increment the number of dominating solutions
          if ( Comp::compare(cmb_arr[j], cmb_arr[i]) > 0 ) {
            cmb_arr[i]->get_fitness_info().n_dominations++;
          }
        }
      }
      if ( cmb_arr[i]->get_fitness_info().n_dominations == 0) {
        cmb_arr[i]->get_fitness_info().rank = 0;
        pareto_fronts[0].push_back(cmb_arr[i]);
      }
    }

    size_t i = 0;
    while (i < pareto_fronts.size() && pareto_fronts[i].size() != 0) {
      pareto_fronts.push_back(empty);
      for (size_t j = 0; j < pareto_fronts[i].size(); ++j) {
        for (size_t k = 0; k < cmb_arr.size(); ++k) {
          cmb_arr[k]->get_fitness_info().n_dominations--;
          if (cmb_arr[k]->get_fitness_info().n_dominations == 0) {
            cmb_arr[k]->get_fitness_info().rank = i + 1;
            pareto_fronts[i + 1].push_back(cmb_arr[k]);
          }
        }
      }
      i++;
    }
    calculated_flags |= FLAG_FRONTS;
  }
 
  Vector<ParentIndSet> gen_breed_pairs(ArgStore& args, Vector<OrgPtr>& old_gen, Vector<OrgPtr>& offspring) {
    _uint arena_size = args.read_custom_double("arena_size", 2);
    size_t offspring_num = old_gen.size();
    SampleDraw sampler(offspring_num, arena_size);
    std::uniform_int_distribution<_uint> selector(0, offspring_num - 1);
    Vector<ParentIndSet> ret( divideup(offspring_num, 2) );

    for (size_t i = 0; 2*i + 1 < offspring_num; ++i) {
      std::vector<_uint> t1 = sampler( args.get_generator() );
      std::vector<_uint> t2 = sampler( args.get_generator() );
      ret[i].first = t1[0];
      ret[i].second = t2[0];
      for (size_t j = 1; j < t1.size(); ++j) {
        if (old_gen[t1[j]]->get_fitness_info().rank < old_gen[ret[i].first]->get_fitness_info().rank) {
          ret[i].first = t1[j];
        }
        if (old_gen[t2[j]]->get_fitness_info().rank < old_gen[ret[i].second]->get_fitness_info().rank) {
          ret[i].second = t2[j];
        }
      }
      //ensure that we don't use the same parent twice with reasonable probably
      if (ret[i].first == ret[i].second) {
        ret[i].second = selector( args.get_generator() );
      }
    }
    return ret;
  }
public:
  NSGAII_TournamentSelector() : Selector<FitType, NSGAII_Comparator<FitType>>() {
    static_assert( std::is_base_of<MultiFitness, FitType>::value, "FitType must be derived from MultiFitness for the NSGAII selector" );
  }
  
  Vector<ParentIndSet> select(ArgStore& args, Vector<OrgPtr>& old_gen, Vector<OrgPtr>& offspring) {
    if (n_objs != old_gen[0]->get_fitness_info().get_n_objs()) {
      n_objs = old_gen[0]->get_fitness_info().get_n_objs();
      pop_stats.resize(n_objs);
    }
    std::vector<OrgPtr> cmb_arr = old_gen;
    cmb_arr.reserve(old_gen.size() + offspring.size());
    for (_uint i = 0; i < offspring.size(); ++i) {
      if (offspring[i]) {
        cmb_arr.push_back(offspring[i]);
      }
    }
    make_fronts(cmb_arr);
    size_t offspring_num = old_gen.size();
    size_t i = 0;
    std::vector<OrgPtr> tmp(offspring_num, NULL);
    _uint k = 0;
    while (i < pareto_fronts.size() && (k + pareto_fronts[i].size()) <= offspring_num) {
      for (size_t j = 0; j < pareto_fronts[i].size(); ++j) {
        tmp[k] = pareto_fronts[i][j];
        ++k;
      }
      ++i;
    }
    //fill in the last elements from the remaining pareto front ranked according to crowding
    if (k < offspring_num) {
      for (size_t ii = 0; ii < pareto_fronts[i].size(); ++ii) {
        pareto_fronts[i][ii]->get_fitness_info().distance = 0;
      }
      //sort by each objective function for crowding evaluation
      for (size_t j = 0; j < n_objs; ++j) {
        this->sort_orgs(j, pareto_fronts[i]);

        //the highest and lowest values should be considered to have no crowding
        pareto_fronts[i].front()->get_fitness_info().distance = std::numeric_limits<double>::infinity();
        pareto_fronts[i].back()->get_fitness_info().distance = std::numeric_limits<double>::infinity();
        double range = pareto_fronts[i].front()->get_fitness(j) - pareto_fronts[i].back()->get_fitness(j);
        if (range == 0) {
          for (size_t ii = 0; ii < pareto_fronts[i].size(); ++ii) {
            pareto_fronts[i][ii]->get_fitness_info().distance = 0;
          }
        } else {
          for (size_t ii = 1; ii < pareto_fronts[i].size() - 1; ++ii) {
            if (pareto_fronts[i][ii]->get_fitness_info().distance != std::numeric_limits<double>::infinity()) {
              double sj_h = pareto_fronts[i][ii-1]->get_fitness(j);
              double sj_l = pareto_fronts[i][ii+1]->get_fitness(j);
              double d_norm = (sj_h - sj_l) / range;
              pareto_fronts[i][ii]->get_fitness_info().distance += d_norm;
            }
          }
        }
      }
      this->sort_orgs(n_objs, pareto_fronts[i]);
      size_t j = 0;
      //select the least crowded individuals
      size_t p_i_size = pareto_fronts[i].size();
      while (j < p_i_size && k < offspring_num) {
        tmp[k] = pareto_fronts[i][j];
        ++j;
        ++k;
      }
      ++i;
    }
    old_gen.swap(tmp);
    //ensure that we aren't keeping null pointers around for safety
    pareto_fronts.clear();

    return gen_breed_pairs(args, old_gen, offspring);
  }
};

//POPULATION_H    
class ConvergenceCriteria {
  public:
    virtual bool evaluate_convergence(Vector<FitnessStats> stats) = 0;
    virtual ~ConvergenceCriteria() = default;
};

template <class FitType, class SelectType, class=void>
class Population {
public:
  typedef typename SelectType::Comp Comp;
private:
  _uint N_BITS;
  _uint N_PARAMS;
  _uint N_OBJS;
  _uint generation = 0;
  _uchar calculated_flags = 0;
  void evaluate_best(Problem<FitType>* prob, double forget_weight=0.0) {
    best_organism->evaluate_fitness(prob);
    for (_uint i = 0; i < N_OBJS; ++i) {
      pop_stats[i].max = best_organism->get_fitness(i);
    }
  }

protected:
  static_assert( std::is_base_of<Selector<FitType, typename SelectType::Comp>, SelectType>::value, "SelectType must be derived from Selector<FitType, Comp>" );
  SelectType sel;
  typedef std::shared_ptr< Organism<FitType> > OrgPtr;
  size_t carryover_num;//How many of the best individuals carry over to the next generation 

  //EXTERNALLY MANAGED POINTERS
  std::shared_ptr<PhenotypeMap> map;
  
  ArgStore args;
  //all offspring from the previous generation
  size_t offspring_num;
  Vector<FitnessStats> pop_stats;
  std::vector<std::shared_ptr<Organism<FitType>>> offspring;
  std::vector<std::shared_ptr<Organism<FitType>>> old_gen;
  size_t min_penalty_ind, max_penalty_ind;
  //which offspring will survive to enter the next breeding round
  size_t survivors_num;
  std::vector<std::shared_ptr<Organism<FitType>>> survivors;
  //guarantee that the best organism appears in the next generation
  _uint best_organism_ind = 0;
  std::shared_ptr<Organism<FitType>> best_organism;
  std::shared_ptr<Organism<FitType>> alltime_best_organism;
  //labels for generating data output
  char** var_labels;
  char** obj_labels;
  std::vector<bool> is_obj_cost;
  int print_penalties = 0;

  //cull in place is slightly faster but less accurate than the standard cull method
  /*bool cull_in_place() {
    size_t j = 0;
    double difference = pop_stats[0].max - pop_stats[0].min;
    //avoid divide by 0
    if (difference == 0) {
      error(CODE_WARN, "All organisms have the same fitness, exiting");
      return true;
    }
    std::uniform_real_distribution<double> dist(0, difference);
    if (survivors.size() < survivors_num) {
      survivors.resize(survivors_num);
    }
    for (size_t i = 0; i < this->offspring_num && j < survivors_num; ++i) {
      // if M_f is the maximum fitness and m_f is the minimum the minimum, while x is the
      // fitness of a given organism, then the probability of survival is x/(M_f-m_f) or 1
      // if M_f = m_f
      if (dist(args.get_generator()) < old_gen[i]->get_fitness(0) - pop_stats[0].min) {
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
  }*/

  /**
   * \brief An implementation of simple roulette selection. This function first sorts the organisms and selects them based on the ratio of their relative fitness to the total relative fitness.
   *
   * \returns True if all organisms have the same fitness (results have converged).
   */
  /*bool cull() {
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
      std::cerr << "Added orgs to survs:";
      for (size_t i = 0; i < survivors_num; ++i) {
        std::cerr << " " << banned[i];
      }
      std::cerr << std::endl;
    }
    free(banned);
  //  best_organism_ind = 0;
    return false;
  }*/
  /*void breed_shuffle() {
    std::uniform_int_distribution<size_t> dist_surv0(0, survivors_num - 1);
    std::vector<Organism<FitType>*> children;
    Organism<FitType>** shuffled_inds = (Organism<FitType>**)malloc(sizeof(Organism<FitType>*)*offspring_num);
    for (size_t i = 0; i < survivors_num; ++i) {
      shuffled_inds[i] = survivors[i].get();
    }
    for (size_t i = 0; i < survivors_num; ++i) {
      size_t ind_o = dist_surv0( args.get_generator() );
      Organism<FitType>* tmp = shuffled_inds[i];
      shuffled_inds[i] = shuffled_inds[ind_o];
      shuffled_inds[ind_o] = tmp;
    }
    size_t last_org_ind = offspring_num;
    for (size_t i = survivors_num; i < offspring_num; ++i) {
      size_t org_ind= dist_surv0( args.get_generator() );
      //ensure that we don't see the same organism breeding with itself
      if (org_ind == last_org_ind) {
        org_ind = (org_ind + 1) % offspring_num;
      }
      shuffled_inds[i] = survivors[org_ind].get();
      last_org_ind = org_ind;
    }
    if (this->offspring_num % 2 == 1) {
      //elitist algorithm, make the first individual in the next generation the previous most fit
      //offspring[0] = std::make_shared<Organism<FitType>>(best_organism);
      offspring[0] = best_organism;
      for (size_t i = 1; 2*i < this->offspring_num; i++) {
        children = shuffled_inds[2*i - 1]->breed(&args, shuffled_inds[2*i]);
        offspring[2*i] = std::shared_ptr<Organism<FitType>>(children[0]);
        offspring[2*i - 1] = std::shared_ptr<Organism<FitType>>(children[1]);
      }
    } else {
      for (size_t i = 0; 2*i + 1 < this->offspring_num; i++) {
        children = shuffled_inds[2*i]->breed(&args, shuffled_inds[2*i + 1]);
        offspring[2*i] = std::shared_ptr<Organism<FitType>>(children[0]);
        offspring[2*i + 1] = std::shared_ptr<Organism<FitType>>(children[1]);
      }
    }
    offspring.swap(old_gen);
  }*/
  /**
   * \brief Apply soft penalties to organisms, this gaurantees that a penalized organism will always have a lower fitness than an unpenalized organism
   *
   * \note Penalized organisms all have a nonzero penalty weight. Let this penalty weight be given by $P$. Let $m$ and $M$ be the minimum and maximum fitness among all unpenalized organisms respectively. Given an organism $O$ with an unpenalized fitness $F(O)$, after penalties are applied,
   * \f[
   * F_new(O)=m - F(O) - max(|M|, |m|)*P
   * \f]
   * By default $F(O)$ is set to zero unless the user specified function evaluate_fitness"("O")" calls O.update"()". 
   * Users may or may not want to assign fitnesses to penalized organisms depending on whether such a fitness is well defined.
   */
  void apply_penalties(Problem<FitType>* prob) {
    //calculate penalties based on the range of fitnesses
    //double penalty_fact = min(abs(pop_stats[0].max - pop_stats[0].min), 1);
    double penalty_fact = abs(pop_stats[0].max - pop_stats[0].min);
    if (penalty_fact == 0) { penalty_fact = 1; }
    if ( pop_stats[0].max == std::numeric_limits<double>::infinity() ) {
      penalty_fact = abs(pop_stats[0].min);
    } else if ( pop_stats[0].min == -std::numeric_limits<double>::infinity() ) {
      penalty_fact = abs(pop_stats[0].max);
    }
    bool penalties_applied = false;
#ifdef USE_LIBOMP
#pragma omp parallel for
    for (size_t i = 0; i < this->offspring_num; ++i) {
      if (old_gen[i]->penalized()) {
        double new_fit = pop_stats[0].min - max(pop_stats[0].max - old_gen[i]->get_fitness(0), 0);
        if (old_gen[i]->get_fitness(0) > pop_stats[0].max) {
          new_fit = pop_stats[0].min;
        }
        new_fit -= penalty_fact*old_gen[i]->get_penalty();
        old_gen[i]->set_fitness(new_fit);
        penalties_applied = true;
      }
    }
#else
    for (size_t i = 0; i < this->offspring_num; ++i) {
      if (old_gen[i]->penalized()) {
        double new_fit = pop_stats[0].min - max(pop_stats[0].max - old_gen[i]->get_fitness(0), 0);
        new_fit -= penalty_fact*old_gen[i]->get_penalty();
        old_gen[i]->set_fitness(new_fit);
        penalties_applied = true;
      }
    }
#endif
    if (!penalties_applied) {
      calculated_flags |= FLAG_STATS_SET | FLAG_BEST_FOUND;
    } else {
      calculated_flags = FLAG_NONE_SET;
    }

    if (best_organism->get_fitness(0) > alltime_best_organism->get_fitness(0)) { 
      for (_uint j = 0; j < args.noise_compensate(); ++j) {
        evaluate_best(prob, args.forget_weight);
      }
      //alltime_best_organism = best_organism->copy();
      alltime_best_organism = best_organism;
      //alltime_best_organism->set_fitness(0, best_organism->get_fitness(0));
    }
    //check to see if there has been a decrease in fitness
    if ( args.noise_compensate() &&
         alltime_best_organism->valid() && 
         alltime_best_organism != best_organism &&
         alltime_best_organism->get_fitness(0) > best_organism->get_fitness(0) ) {
      //run more evaluations to see if the new organism is actually better
      for (_uint i = 0; i < args.noise_compensate(); ++i) {
        alltime_best_organism->evaluate_fitness(prob);
        evaluate_best(prob, args.forget_weight);
      }
      if ( alltime_best_organism->get_fitness(0) > best_organism->get_fitness(0) ) {
        //best_organism->swap(alltime_best_organism);
        best_organism = alltime_best_organism;
        pop_stats[0].max = best_organism->get_fitness(0);
      }
    }
  }
  void find_best_organism() {
    /*for (int j = 0; j < N_OBJS; ++j) {
      pop_stats[j].max = best_organism->get_fitness(j);
      pop_stats[j].min = best_organism->get_fitness(j);
      pop_stats[j].mean = best_organism->get_fitness(j) / offspring_num;
      for (size_t i = 0; i < offspring_num; ++i) {
        double fitness_i = old_gen[i]->get_fitness(j);
        if (fitness_i > pop_stats[j].max) {
          //TODO: make this usefully track multiple objectives
          if (j == 0) {
            set_best_organism(i);
          }
        }
        if (fitness_i < pop_stats[j].min) {
          pop_stats[j].min = fitness_i;
        }
        pop_stats[j].mean += fitness_i / offspring_num;
      }
      //calculate the variance
      pop_stats[j].var = 0;
      for (size_t i = 0; i < offspring_num; ++i) {
        double fitness_i = old_gen[i]->get_fitness(j);
        pop_stats[j].var += (fitness_i - pop_stats[j].mean)*(fitness_i - pop_stats[j].mean);
      }
    }*/

    set_best_organism( SelectType::find_best_organism(old_gen, pop_stats) );
    
    calculated_flags |= FLAG_STATS_SET | FLAG_BEST_FOUND;
  }
  void breed(Vector<ParentIndSet>&& parents) {
    if (2*parents.size() + 1 < offspring_num) {
      error(CODE_MISC, "Too few parents supplied in the parents array.");
    }
    for (_uint i = 0; 2*i + 1< offspring_num; ++i) {
      _uint par1_ind = parents[i].first;
      _uint par2_ind = parents[i].second;
      std::pair<OrgPtr, OrgPtr> children = old_gen[par1_ind]->breed(args, old_gen[par2_ind]);
      offspring[2*i] = children.first;
      offspring[2*i + 1] = children.second;
    }
    offspring[offspring_num - 1] = best_organism;//elitism
    old_gen.swap(offspring);
  }
  /*void breed() {
    find_best_organism();
    std::uniform_int_distribution<size_t> dist_surv0(0, survivors_num - 1);
    std::uniform_int_distribution<size_t> dist_surv1(0, survivors_num - 2);
    std::vector<Organism<FitType>*> children;
    size_t* shuffled_inds = (size_t*)malloc(sizeof(size_t)*survivors_num);
    if (this->offspring_num % 2 == 1) {
      //offspring[0] = std::make_shared<Organism<FitType>>(best_organism);
      offspring[0] = best_organism;
      for (size_t i = 1; 2*i < this->offspring_num; i++) {
        size_t par1_i = dist_surv0( args.get_generator() );
        //use the survivors_num-1 distribution to guarantee different parents
        size_t par2_i = dist_surv1( args.get_generator() );
        if (par2_i >= par1_i) {
    par2_i++;
        }
        children = survivors[par1_i].get()->breed(&args, survivors[par2_i].get());
        offspring[2*i] = std::shared_ptr<Organism<FitType>>(children[0]);
        offspring[2*i - 1] = std::shared_ptr<Organism<FitType>>(children[1]);
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
        offspring[2*i] = std::shared_ptr<Organism<FitType>>(children[0]);
        offspring[2*i + 1] = std::shared_ptr<Organism<FitType>>(children[1]);
      }
    }
    free(shuffled_inds);
    offspring.swap(old_gen);
  }*/
  
  void calculate_distances() {
    for (_uint i = 0; i < offspring_num; ++i) {
      old_gen[i]->get_fitness_info().distance = 0;
    }
    for (_uint i = 0; i < N_OBJS; ++i) {
      sel.sort_orgs(i, old_gen);
      for (_uint j = 1; j < offspring_num - 1; ++j) {
        double tmp_dist = old_gen[j + 1]->get_fitness(i) - old_gen[j - 1]->get_fitness(i);
        old_gen[j]->get_fitness_info().distance += tmp_dist*tmp_dist;
      }
      old_gen[0]->get_fitness_info().distance = std::numeric_limits<double>::infinity();
      old_gen[offspring_num - 1]->get_fitness_info().distance = std::numeric_limits<double>::infinity();
    }
    calculated_flags |= FLAG_DIST_SET;
  }
  void hypermutate() {
    if ( !(calculated_flags & FLAG_DIST_SET) ) {
      calculate_distances();
    }
    if ( !(calculated_flags & FLAG_BEST_FOUND) ) {
      find_best_organism();
    }
    //sort by distance
    sel.sort_orgs(get_n_objs(), old_gen);
    //select the half of the most crowded individuals in the first front
    for (_uint i = offspring_num - 1; i > args.get_replacement_fraction()*offspring_num; --i) {
  //    if (old_gen[i] != best_organism_current) {
        old_gen[i]->randomize(args);
  //    }
    }
  }
  void set_best_organism(_uint i, bool force=false, _uint j=0) {
    if (j < N_OBJS) {
      pop_stats[j].max = old_gen[i]->get_fitness(j);
      /*Organism<FitType> tmp_org = old_gen[i]->copy();
      if (force || tmp_org > best_organism) {
        best_organism = tmp_org;
        for (_uint j = 0; j < N_OBJS; ++j) {
          best_organism->set_fitness( j, old_gen[i]->get_fitness(j) );
        }
        best_organism_ind = i;
      }*/
      if ( force || !best_organism || Comp::compare(old_gen[i], best_organism) > 0 ) {
        best_organism = old_gen[i];
        best_organism_ind = i;
      }
      calculated_flags |= VALID_BEST;
    }
  }
  size_t find_first_unpenalized(Problem<FitType>* prob) {
    size_t start_i = 0;
    if ( best_organism && this->best_organism->valid() ) {
      if ( args.noise_compensate() ) {
        evaluate_best(prob, args.forget_weight);
      }
      pop_stats[0].max = best_organism->get_fitness(0);
    } else {
      do {
        if (start_i == offspring_num) {
          error(CODE_MISC, "All organisms in population had an applied penalty.");
        }
        old_gen[start_i]->apply_penalty(0);
        old_gen[start_i]->reset_fitness();
        old_gen[start_i]->evaluate_fitness(prob);
        ++start_i;
      } while( this->old_gen[start_i - 1]->penalized() );
      set_best_organism(start_i - 1, args.noise_compensate());
      //alltime_best_organism = best_organism->copy();
      alltime_best_organism = best_organism;
      --start_i;
    }
    pop_stats[0].min = best_organism->get_fitness(0);
    return start_i;
  }
  void handle_multiples() {
    if ( args.perturb_multiples() ) {
      for (_uint i = 0; i < offspring_num; ++i) {
        for (_uint j = 0; j < i; ++j) {
          if ( *(old_gen[i]) == *(old_gen[j]) ) {
            old_gen[i]->mutate(args);
          }
        }
      }
    }
  }

public:
//    Population(_uint pn_bits, _uint pn_objs, PhenotypeMap* p_map, ArgStore p_args);
//    Population(_uint pn_bits, _uint pn_objs, Organism<FitType>* tmplt, PhenotypeMap* p_map, ArgStore p_args);
  Population(_uint pn_bits, _uint pn_objs, std::shared_ptr<PhenotypeMap> p_map, bool latin=true) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs),
  is_obj_cost(N_OBJS, false)
  {
    map = p_map;
    createOrganisms(NULL, latin);
  }
  Population(_uint pn_bits, _uint pn_objs, Organism<FitType>* tmplt, std::shared_ptr<PhenotypeMap> p_map) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs),
  is_obj_cost(N_OBJS, false)
  {
    map = p_map;
    createOrganisms(tmplt, false);
  }
  Population(_uint pn_bits, _uint pn_objs, std::shared_ptr<PhenotypeMap> p_map, ArgStore p_args, bool latin=true) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs),
  args(p_args),
  is_obj_cost(N_OBJS, false)
  {
    map = p_map;
    createOrganisms(NULL, latin);
  }
  Population(_uint pn_bits, _uint pn_objs, Organism<FitType>* tmplt, std::shared_ptr<PhenotypeMap> p_map, ArgStore p_args) :
  N_BITS(pn_bits),
  N_OBJS(pn_objs),
  args(p_args),
  is_obj_cost(N_OBJS, false)
  {
    map = p_map;
    createOrganisms(tmplt, false);
  }
  void createOrganisms(Organism<FitType>* tmplt, bool latin) {
    pop_stats.resize(N_OBJS);
    //this->survivors_num = args.get_survivors();
    this->offspring_num = args.get_pop_size();
    /*if (this->offspring_num % 2 == 0) {
      this->offspring_num++;
    }*/
    //we need to keep the old and new generation in separate arrays to avoid overwriting data
    this->offspring.insert(this->offspring.end(), this->offspring_num, std::shared_ptr<Organism<FitType>>(NULL));
    if (tmplt) {
      this->old_gen.push_back(std::make_shared<Organism<FitType>>(*tmplt));
      //initally fill up the offspring randomly
      for (size_t i = 1; i < this->offspring_num; ++i) {
        this->old_gen.push_back( std::make_shared<Organism<FitType>>(N_BITS, N_OBJS, map) );
        this->old_gen.back()->randomize(args, tmplt);
      }
    } else {
      if (latin) {
        Shuffle samp(offspring_num);
        std::uniform_real_distribution<double> in_cube_dist(0, 1);
        Vector<_uint> row;
        
        this->old_gen.reserve(this->offspring_num);
        for (size_t i = 0; i < this->offspring_num; ++i) {
          this->old_gen.push_back( std::make_shared<Organism<FitType>>(N_BITS, N_OBJS, map) );
        }
        for (size_t i = 0; i < map->get_num_params(); ++i) {
          row = samp( args.get_generator() ); 
          double row_width = ( map->get_range_max(i) - map->get_range_min(i) )/offspring_num;
          for (size_t j = 0; j < offspring_num; ++j) {
            double x = in_cube_dist( args.get_generator() );
            x = (x + row[j])*row_width + map->get_range_min(i);
            this->old_gen[j]->set_real(i, x);
          }
        }
      } else {
        this->old_gen.reserve(this->offspring_num);
        //initally fill up the offspring randomly
        for (size_t i = 0; i < this->offspring_num; ++i) {
          this->old_gen.push_back( std::make_shared<Organism<FitType>>(N_BITS, N_OBJS, map) );
          this->old_gen[i]->randomize(args);
        }
      }
    }
    for (_uint i = 0; i < N_OBJS; ++i) {
      this->pop_stats[i].max = -std::numeric_limits<double>::infinity();
      this->pop_stats[i].min = std::numeric_limits<double>::infinity();
    }
    char buf[OUT_BUF_SIZE];
  //  var_labels.resize(map->get_num_params());
    N_PARAMS = map->get_num_params();
    var_labels = (char**)malloc(sizeof(char*)*N_PARAMS);
    for (_uint j = 0; j < map->get_num_params(); ++j) {
      var_labels[j] = (char*)malloc(sizeof(char)*OUT_BUF_SIZE);
      snprintf(var_labels[j], OUT_BUF_SIZE, "x_%d", j);
    }
    obj_labels = (char**)malloc(sizeof(char*)*N_OBJS);
    for (_uint j = 0; j < N_OBJS; ++j) {
      obj_labels[j] = (char*)malloc(sizeof(char)*OUT_BUF_SIZE);
      snprintf(obj_labels[j], OUT_BUF_SIZE - 1, "f_%d(x)", j);
    }

    calculated_flags = FLAG_NONE_SET;
  }
  ~Population() {
    for (size_t i = 0; i < this->offspring_num; ++i) {
      old_gen[i].reset();
    }
    if (var_labels) {
      for (_uint i = 0; i < N_PARAMS; ++i) {
        free(var_labels[i]);
      }
      free(var_labels);
    }
    if (obj_labels) {
      for (_uint i = 0; i < N_OBJS; ++i) {
        free(obj_labels[i]);
      }
      free(obj_labels);
    }
  }
  Population(Population& o) :
  args(o.args),
  best_organism(o.best_organism)
  {
    N_BITS = o.get_n_bits();
    N_PARAMS = o.N_PARAMS;
    N_OBJS = o.N_OBJS;
    map = o.map;
    old_gen = o.old_gen;
    offspring = o.offspring;
    pop_stats.resize(N_OBJS);
    var_labels = (char**)malloc(sizeof(char*)*N_PARAMS);
    obj_labels = (char**)malloc(sizeof(char*)*N_OBJS);

    for (_uint i = 0; i < N_OBJS; ++i) {
      pop_stats[i].min = o.pop_stats[i].min;
      pop_stats[i].max = o.pop_stats[i].max;
      pop_stats[i].var = o.pop_stats[i].var;
      size_t len = strlen(o.obj_labels[i]) + 1;
      obj_labels[i] = (char*)malloc( sizeof(char)*len);
      snprintf(obj_labels[i], len, "%s", o.obj_labels[i]);
    }
    for (_uint i = 0; i < N_PARAMS; ++i) {
      size_t len = strlen(o.var_labels[i]) + 1;
      var_labels[i] = (char*)malloc(sizeof(char)*len);
      snprintf(obj_labels[i], len, "%s", o.obj_labels[i]);
    }
    calculated_flags = FLAG_NONE_SET;
  }
  Population& operator=(Population& o) {
    N_BITS = o.N_BITS;
    N_OBJS = o.N_OBJS;
    map = o.map;
    pop_stats = o.pop_stats;
    //free existing pointers and reallocate array
    for (_uint i = 0; i < N_PARAMS; ++i) { if(var_labels[i]) { free(var_labels[i]); } }
    for (_uint i = 0; i < N_OBJS; ++i) { if(obj_labels[i]) { free(obj_labels[i]); } }
    var_labels = (char**)realloc(var_labels, sizeof(char*)*N_PARAMS);
    obj_labels = (char**)realloc(obj_labels, sizeof(char*)*N_OBJS);
    //deep copy
    for (_uint i = 0; i < N_PARAMS; ++i) {
      _uint len = strlen(o.var_labels[i]) + 1;
      var_labels[i] = (char*)malloc(sizeof(char)*len);
      strncpy(var_labels[i], o.var_labels[i], len);
    }
    for (_uint i = 0; i < N_OBJS; ++i) {
      _uint len = strlen(o.obj_labels[i]) + 1;
      obj_labels[i] = (char*)malloc(sizeof(char)*len);
      strncpy(obj_labels[i], o.obj_labels[i], len);
    }
    old_gen = o.old_gen;
    offspring = o.offspring;
    best_organism = o.best_organism;

    calculated_flags = FLAG_NONE_SET;
    
    return *this;
  }
  Population(Population&& o) :
  map(o.map),
  args(std::move(o.args)),
  best_organism(std::move(o.best_organism))
  {
    N_BITS = o.get_n_bits();
    for (size_t i = 0; i < this->offspring_num; ++i) {
      old_gen[i].reset();
    }
    old_gen = std::move(o.old_gen);
    offspring = std::move(o.offspring);
    pop_stats = o.pop_stats;
    var_labels = o.var_labels;
    obj_labels = o.obj_labels;
    o.var_labels = NULL;
    o.obj_labels = NULL;

    calculated_flags = FLAG_NONE_SET;
  }

  void set_convergence_type(ConvergenceCriteria* conv);
  void set_penalty_printing(bool val = true) { print_penalties = (val)? 1 : 0; }
  void resize_population(_uint new_size) {
    size_t old_size = old_gen.size();
    if (new_size > old_size) {
      offspring.insert( offspring.end(), new_size - old_size, std::shared_ptr<Organism<FitType>>(NULL) );
      old_gen.reserve(new_size);
      for (size_t i = old_size; i < old_size; ++i) {
        old_gen.push_back( std::make_shared<Organism<FitType>>(N_BITS, N_OBJS, map) );
        old_gen.back()->randomize(args);
      }
    } else {
      offspring.resize(new_size);
      old_gen.resize(new_size);
    }
    args.set_pop_size(new_size);

    calculated_flags = FLAG_NONE_SET;
  }
  /*void set_n_survivors(_uint new_size) {
    size_t old_size = survivors_num;
    if (new_size > old_size) {
      survivors.insert( offspring.end(), new_size - old_size, std::shared_ptr<Organism<FitType>>(NULL) );
    } else {
      survivors.resize(new_size);
    }
    args.set_survivors(new_size);
  }*/
#ifdef USE_LIBOMP
  template <typename T=FitType>
  void evaluate_async(Problem<T>* prob) {
    if (N_OBJS == 1) { 
      for (_uint i = 0; i < offspring_num; ++i) {
        old_gen[i]->apply_penalty(0);
      }
      //calculate averages for organisms that appear twice in the population
      if ( args.average_multiples() ) {
#pragma omp parallel for
        for (_uint i = 0; i < offspring_num; ++i) {
          Vector<_uint> identical_set;
          double avg_fit = 0.0;
          bool apply_averages = true;
          for (_uint j = 0; j < offspring_num; ++j) {
            if ( i == j || *(old_gen[j]) == *(old_gen[i]) ) {
              identical_set.push_back(j);
              //ensure that we only calculate the identical set once
              if (i < j) {
                apply_averages = false;
              }
            }
          }

#pragma omp parallel for
          for (_uint j = 0; j < identical_set.size(); ++j) {
            for (_uint k = 0; k < args.noise_compensate() + 1; ++k) {
              old_gen[i]->evaluate_fitness(prob);
            }
          }
          
          //don't recalculate if we don't have to
          if (apply_averages) {
            for (auto it = identical_set.begin(); it != identical_set.end(); ++it) {
              old_gen[*it]->update(0, old_gen[i]->get_fitness(0));
            }
            if (old_gen[i]->get_fitness(0) > best_organism->get_fitness(0) && !old_gen[i]->penalized()) {
                set_best_organism(i);
            }
          }
        }
      } else {
        Vector<_uint*> skip_set;
        for (_uint i = 0; i < this->offspring_num; ++i) {
          if (args.skip_multiples()) {
            for (size_t j = 0; j < i; ++j) {
              if (*(old_gen[j]) == *(old_gen[i])) {
                _uint* tmp = (_uint*)malloc(sizeof(_uint)*2);
                tmp[0] = j; tmp[1] = i;
                skip_set.push_back(tmp);
              }
            }
          } else if (args.perturb_multiples()) {
            for (size_t j = 0; j < i; ++j) {
              if (*(old_gen[j]) == *(old_gen[i])) {
                old_gen[j]->mutate(args);
              }
            }
          }
        }

#pragma omp parallel for
        for (_uint i = 0; i < offspring_num; ++i) {
          _uint j = 0;
          for (; j < skip_set.size(); ++j) {
            if ( skip_set[j][0] == i ) { break; }
          }
          if (j >= skip_set.size()) {
            if ( args.verbose() ) {
              std::cout << "Now evaluating organism " << i << std::endl;
            }
            old_gen[i]->evaluate_fitness(prob);
            for (_uint k = 0; k < args.noise_compensate(); ++k) {
              old_gen[i]->evaluate_fitness(prob);
            }
          }
        }

        _uint first_valid_i = 0;
        if ( best_organism->valid() ) {
          if ( args.noise_compensate() ) {
            evaluate_best(prob, args.forget_weight);
          }
          pop_stats[0].max = best_organism->get_fitness(0);
          pop_stats[0].min = best_organism->get_fitness(0);
        } else {
          //iterate until we find an organism that isn't penalized and set it to be the best
          do {
            if (first_valid_i == offspring_num) {
              error(CODE_MISC, "All organisms in population had applied penalty.");
            }
            ++first_valid_i;
          } while( old_gen[first_valid_i]->penalized() );
          set_best_organism(first_valid_i - 1);
          //alltime_best_organism = best_organism->copy();
    alltime_best_organism = best_organism;
          pop_stats[0].max = old_gen[first_valid_i]->get_fitness(0);
          pop_stats[0].min = old_gen[first_valid_i]->get_fitness(0);
        }

        //this can't be parallelized easily
        for (_uint i = 0; i < offspring_num; ++i) {
          _uint j = 0;
          for (; j < skip_set.size(); ++j) {
            if ( skip_set[j][0] == i ) { break; }
          }
          //if we are in the skip set, then set fitness accordingly
          if (j < skip_set.size()) {
            uint prev_ind = skip_set[j][1];
            old_gen[i]->update(0, old_gen[prev_ind]->get_fitness(0));
            old_gen[i]->apply_penalty(old_gen[prev_ind]->get_penalty());
          } else if (old_gen[i]->get_fitness(0) > best_organism->get_fitness(0) && !old_gen[i]->penalized()) {
            //check the organism again to make sure this isn't a fluke
            for (_uint j = 0; j < args.noise_compensate(); ++j) {
              old_gen[i]->evaluate_fitness(prob);
              evaluate_best(prob, args.forget_weight);
            }
            set_best_organism(i);
            if (old_gen[i]->get_fitness(0) < pop_stats[0].min) {
              pop_stats[0].min = old_gen[i]->get_fitness(0);
            }
          }
        }

        //free allocated memory
        for (_uint j = 0; j < skip_set.size(); ++j) {
          free(skip_set[j]);
        }
      }
    } else {
      //TODO: figure out what the default behavior should be
    }
    double penalty_fact;
    if (pop_stats[0].max > 0) {
      penalty_fact = pop_stats[0].min;
    } else {
      penalty_fact = -pop_stats[0].min;
    }
    bool penalties_applied = false;
#pragma omp parallel for
    for (size_t i = 0; i < this->offspring_num; ++i) {
      if (old_gen[i]->penalized()) {
        old_gen[i]->update(pop_stats[0].min - penalty_fact*old_gen[i]->get_penalty());
        penalties_applied = true;
      }
    }
    if (!penalties_applied) {
      calculated_flags |= FLAG_STATS_SET | FLAG_BEST_FOUND;
    } else {
      calculated_flags = FLAG_NONE_SET;
    }
    apply_penalties(prob);
  }
#endif
  //function for population where organism FitType has a member average_fitness
  template <typename T = FitType> inline
  typename enable_if_c< has_average_fitness<T, void(T&)>::value, void >::type
  evaluate(Problem<T>* prob) {
    size_t start_i = find_first_unpenalized(prob);
    handle_multiples();

    for (_uint i = start_i; i < this->offspring_num; ++i) {
      this->old_gen[i]->evaluate_fitness(prob);
      if (SelectType::use_offspring) {
        if (i < offspring.size() && offspring[i]) { offspring[i]->evaluate_fitness(prob); }
      }
      //average fitnesses with previous organisms if appropriate
      if ( this->args.average_multiples() ) {
        for (_uint j = 0; j < i; ++j) {
          if ( *(this->old_gen[i]) == *(this->old_gen[i]) ) {
            this->old_gen[i]->get_fitness_info().average_fitness( this->old_gen[j]->get_fitness_info() );
          }
        }
      }
      // update the max and min fitnesses if we need to
      for (_uint j = 0; j < N_OBJS; ++j) {
        if ( this->old_gen[i]->get_fitness(j) > this->best_organism->get_fitness(j)
        && !(this->old_gen[i]->penalized()) ) {
          //check the organism again to make sure this isn't a fluke
          for (_uint j = 0; j < this->args.noise_compensate(); ++j) {
            this->old_gen[i]->evaluate_fitness(prob);
          }
          if (!(this->old_gen[i]->penalized())
           &&  (this->old_gen[i]->get_fitness(j) > this->best_organism->get_fitness(j)
             || this->old_gen[i]->get_fitness(j) > pop_stats[j].max)) {
            this->set_best_organism(i);
          }
        }
        if (this->old_gen[i]->get_fitness(j) < this->pop_stats[j].min) {
          this->pop_stats[j].min = this->old_gen[i]->get_fitness(j);
        }
      }
    }
    this->apply_penalties(prob);
    for (_uint i = 0; i < this->offspring_num; ++i) {
      if (!(this->old_gen[i]->penalized()) &&  this->old_gen[i]->get_fitness(0) > this->best_organism->get_fitness(0)) {
        this->set_best_organism(i, true);
      }
    }
    this->pop_stats[0].max = this->best_organism->get_fitness();
  }
  //function for population where organism FitType does not have a member average_fitness
  template <typename T = FitType> inline
  typename enable_if_c< !has_average_fitness<T, void(T&)>::value, void >::type
  evaluate(Problem<T>* prob) {
    size_t start_i = find_first_unpenalized(prob);
    handle_multiples();
   
    //calculate averages for organisms that appear twice in the population
    for (_uint i = start_i; i < this->offspring_num; ++i) {
      bool found_identical = false;
      //look for duplicates of the current organism
      for (size_t j = 0; j < i; ++j) {
        //handle them
        if (!args.perturb_multiples() && *(this->old_gen[j]) == *(this->old_gen[i]) ) {
          for (_uint k = 0; k < N_OBJS; ++k) { this->old_gen[i]->update( k, this->old_gen[j]->get_fitness(k) ); }
          found_identical = true;
          break;
        }
      }
      if (!found_identical) {
        old_gen[i]->evaluate_fitness(prob);
      }
      if (SelectType::use_offspring) {
        if (i < offspring.size() && offspring[i]) { offspring[i]->evaluate_fitness(prob); }
      }
      // update the max and min fitnesses if we need to
      for (_uint j = 0; j < N_OBJS; ++j) {
        if ( !(this->old_gen[i]->penalized()) && Comp::compare(old_gen[i], best_organism) > 0 ) {
          this->set_best_organism(i, false, j);
        }
        if (this->old_gen[i]->get_fitness(j) < this->pop_stats[j].min) {
          this->pop_stats[j].min = this->old_gen[i]->get_fitness(j);
        }
      }
    }
    this->pop_stats[0].max = this->best_organism->get_fitness();
    this->apply_penalties(prob);
  }
  //void evaluate(Problem<FitType>* prob) { evaluate_imp(prob, NULL); }
  bool iterate(ConvergenceCriteria* conv = NULL) {
    find_best_organism();
    //check for hypermutation
    for (_uint i = 0; i < N_OBJS; ++i) {
      double range_ratio = (pop_stats[i].max - pop_stats[i].min)/ pop_stats[i].max;
      if (range_ratio < 0) {
        range_ratio *= -1;
      }
      if (1.0 - range_ratio > args.get_hypermutation_threshold()) {
        hypermutate();
        break;
      }
    }
    breed( sel.select(args, old_gen, offspring) );
    calculated_flags &= !FLAG_FRONTS;
    generation++;
    if (conv) {
      return conv->evaluate_convergence(pop_stats);
    } else {
      return (generation > args.get_num_gens());
    }
    calculated_flags = FLAG_NONE_SET;
  }
  void run(Problem<FitType>* prob) {
    evaluate(prob);
    if ( this->args.wait_for_con() ) {
      size_t i = 1;
      unsigned streak = 0;
      double prev_ftns = this->get_best_organism()->get_fitness(0);
      while (i < MAX_NUM_GENS) {
        if (this->iterate()) {
          break;
        }
#ifdef USE_LIBOMP
        if (args.async()) {
          evaluate_async(prob);
        } else {
          evaluate(prob);
        }
#else
        evaluate(prob);
#endif

        if (this->get_best_organism()->get_fitness(0) > prev_ftns) {
          prev_ftns = this->get_best_organism()->get_fitness(0);
          streak = 0;
        }

        streak++;
        i++;
        if (streak > this->args.get_num_gens()) {
          break;
        }
      }

      if (i >= MAX_NUM_GENS) {
        std::cout << "failed to converge after " << i << " generations." << std::endl;
      } else {
        std::cout << "converged to result after " << i << " generations." << std::endl;
      }
    } else {
      if (this->args.verbose()) {
        std::cout << "Now evaluating generation 0..." << std::endl;
      }
      for (size_t i = 1; i < this->args.get_num_gens() + 1; ++i) {
        if (this->args.verbose()) {
          std::cout << "Now evaluating generation " << i << "..." << std::endl;
        }
        //produce the next generation in the population
        if (this->iterate()) {
          break;
        }
        this->evaluate(prob);
      } 
    }
  }
  std::shared_ptr< Organism<FitType> > get_best_organism(size_t i = 0) {
    if ( pop_stats[i].max != best_organism->get_fitness(i) || (calculated_flags & FLAG_BEST_FOUND) == 0 ) {
      find_best_organism();
    }
    return best_organism;
    /*if (i == 0) {
      std::shared_ptr<Organism<FitType>> tmp_org = std::make_shared<Organism<FitType>>( best_organism->copy() );
      //for (_uint i = 0; i < N_OBJS; ++i) {
      //tmp_org->update(i, pop_stats[i].max);
      //}
      return tmp_org;
    } else {
      sel.sort_orgs(0, old_gen);
      return old_gen[i];
    }*/
  }

  std::shared_ptr< Organism<FitType> > get_organism(size_t i) {
    if (i >= old_gen.size())
      error(CODE_ARG_RANGE, "Attempt to access invalid index %d when the maximum allowed is %d.", i, old_gen.size());
    return old_gen[i];
  }
  std::shared_ptr< Organism<FitType> > get_child(size_t i) {
    if (i >= offspring.size())
      error(CODE_ARG_RANGE, "Attempt to access invalid index %d when the maximum allowed is %d.", i, offspring.size());
    return offspring[i];
  }
  
  Vector<String> get_best_header() {
    String def;
    Vector<String> ret;
    ret.reserve( N_PARAMS + print_penalties + N_OBJS );
    char buf[OUT_BUF_SIZE];

    //print information about the best organism
    for (_uint j = 0; j < N_PARAMS; ++j) {
      snprintf(buf, OUT_BUF_SIZE - 1, "best %s", var_labels[j]);
      ret.push_back( String(buf) );
    }
    if (print_penalties != 0) {
      ret.push_back( String("best penalty") );
    }
    for (_uint j = 0; j < N_OBJS; ++j) {
      snprintf(buf, OUT_BUF_SIZE - 1, "best %s", obj_labels[j]);
      ret.push_back( String(buf) );
    }
    
    return ret;
  }
  Vector<String> get_header() {
    String def;
    Vector<String> ret;
    ret.reserve( old_gen.size()*(N_PARAMS + print_penalties + N_OBJS) );

    //print information about every other organism in the population
    for (_uint i = 0; i < old_gen.size(); ++i) {
      for (_uint j = 0; j < N_PARAMS; ++j) {
        ret.push_back( String(var_labels[j]) );
      }

      if (print_penalties != 0) {
        ret.push_back( String("Penalized") );
      }

      for (_uint j = 0; j < N_OBJS; ++j) {
        ret.push_back( String(obj_labels[j]) );
      }
    }
    
    return ret;
  }
  Vector<String> get_best_data() {
    sel.sort_orgs(0, old_gen);
    char buf[OUT_BUF_SIZE];
    String def;
    Vector<String> ret(N_OBJS + print_penalties + map->get_num_params(), def);

    size_t param_o = 0;
    size_t fitness_o = map->get_num_params() + print_penalties;
    //print information about the best organism
    for (_uint j = 0; j < map->get_num_params(); ++j) {
      ret[j] = best_organism->get_chromosome_string(j);
    }
    if (print_penalties) {
      snprintf(buf, OUT_BUF_SIZE - 1, "%f", best_organism->get_penalty());
      ret[map->get_num_params()] = buf;
    }
    for (_uint j = 0; j < N_OBJS; ++j) {
      if ( args.noise_compensate() ) {
        double std_dev = sqrt( best_organism->get_fitness_info().get_uncertainty(j) );
        if (is_obj_cost[j]) {
          snprintf(buf, OUT_BUF_SIZE - 1, "%f\u00B1%f", best_organism->get_cost(j), std_dev);
        } else {
          snprintf(buf, OUT_BUF_SIZE - 1, "%f\u00B1%f", best_organism->get_fitness(j), std_dev);
        }
      } else {
        if (is_obj_cost[j]) {
          snprintf(buf, OUT_BUF_SIZE - 1, "%f", best_organism->get_cost(j));
        } else {
          snprintf(buf, OUT_BUF_SIZE - 1, "%f", best_organism->get_fitness(j));
        }
      }
      ret[fitness_o + j] = buf;
    }

    return ret;
  }
  Vector<String> get_pop_data() {
    sel.sort_orgs(0, old_gen);
    _uint span = N_OBJS + print_penalties + map->get_num_params();
    String def;
    Vector<String> ret(span*offspring_num, def);
    char buf[OUT_BUF_SIZE];

    //print information about the rest of the population
    for (_uint i = 0; i < old_gen.size(); ++i) {
      size_t param_o = i*span;
      size_t fitness_o = param_o + map->get_num_params() + print_penalties;

      //print out the parameters
      for (_uint j = 0; j < map->get_num_params(); ++j) {
        ret[param_o + j] = old_gen[i]->get_chromosome_string(j);
      }
      //print out the penalty applied to the organism
      if (print_penalties) {
        snprintf(buf, OUT_BUF_SIZE - 1, "%f", old_gen[i]->get_penalty());
        ret[param_o + map->get_num_params()] = buf;
      }
      //print out the fitness value(s)
      for (_uint j = 0; j < N_OBJS; ++j) {
        if ( args.noise_compensate() ) {
          double std_dev = sqrt( old_gen[i]->get_fitness_info().get_uncertainty(j) );
          if (is_obj_cost[j]) {
            snprintf(buf, OUT_BUF_SIZE - 1, "%f\u00B1%f", old_gen[i]->get_cost(j), std_dev);
          } else {
            snprintf(buf, OUT_BUF_SIZE - 1, "%f\u00B1%f", old_gen[i]->get_fitness(j), std_dev);
          }
        } else {
          if (is_obj_cost[j]) {
            snprintf(buf, OUT_BUF_SIZE - 1, "%f", old_gen[i]->get_cost(j));
          } else {
            snprintf(buf, OUT_BUF_SIZE - 1, "%f", old_gen[i]->get_fitness(j));
          }
        }
        
        ret[fitness_o + j] = buf;
      }
    }

    return ret;
  }

  Vector<std::pair<std::shared_ptr<Organism<FitType>>, _uint>> get_species_list(double tolerance=0.1, _uint dimension_threshold=1) {
    Vector<std::pair<std::shared_ptr<Organism<FitType>>, _uint>> ret;
    if (best_organism) { ret.emplace_back(best_organism, 0); }
    for (_uint i = 0; i < old_gen.size(); ++i) {
      bool found = false;
      for (_uint j = 0; j < ret.size(); ++j) {
        _uint num_differences = 0;
        for (_uint k = 0; k < map->get_num_params(); ++k) {
          if ( abs( old_gen[i]->read_real(k) - ret[j].first->read_real(k) ) > tolerance*(map->get_range_max(k) - map->get_range_min(k)) ) {
            ++num_differences;
            if (num_differences > dimension_threshold) { break; }
          }
        }
        if (num_differences < dimension_threshold) {
          ++ret[j].second;
          found = true;
          break;
        }
      }
      if (!found) {
        ret.emplace_back(old_gen[i], 1);
      }
    }
    return ret;
  }

  FitnessStats get_pop_stats(_uint i = 0) { return pop_stats[i]; }
  DEPRECATED("get_min_fitness is deprecated, use get_pop_stats instead") double get_min_fitness(_uint i = 0) {
    return pop_stats[i].min;
  }
  DEPRECATED("get_max_fitness is deprecated, use get_pop_stats instead") double get_max_fitness(_uint i = 0) {
    return pop_stats[i].max;
  }

  /**
   * \brief Sets the objective referenced by index ind to use a fitness (corresponding to a maximization problem).
   *
   * \param ind	The index of the parameter to set
   * \seealso set_cost
   */
  void update(_uint ind) { is_obj_cost[ind] = false; }
  /**
   * \brief Sets the objective referenced by index ind to use a cost (corresponding to a minimization problem).
   *
   * \param ind	The index of the parameter to set
   * \seealso update, set_fitness
   */
  void set_cost(_uint ind) { is_obj_cost[ind] = true; }
  /**
   * \brief Sets the objective referenced by index ind to use a cost (corresponding to a minimization problem).
   *
   * \param ind	The index of the parameter to set
   * \seealso update
   */
  void set_fitness(_uint ind) { is_obj_cost[ind] = false; }
  void set_var_label(_uint ind, String val) {
    if (ind >= N_PARAMS) {
      error(CODE_WARN, "Invalid parameter index %u provided for set_var_label. The parameter index must be less than %u.", ind, N_PARAMS);
    } else {
      size_t vs = val.size() + 1;
      if (vs > OUT_BUF_SIZE) {
        free(var_labels[ind]);
        var_labels[ind] = (char*)malloc(sizeof(char)*vs);
      }
      for (int j = 0; j < vs; ++j) {
        var_labels[ind][j] = val[j];
      }
      var_labels[ind][vs - 1] = 0;
    }
  }
  void set_obj_label(_uint ind, String val) {
    if (ind >= N_OBJS) {
      error(CODE_WARN, "Invalid objective index %u provided for set_obj_label. The objective index must be less than %u.", ind, N_OBJS);
    } else {
      size_t vs = val.size() + 1;
      if (vs > OUT_BUF_SIZE) {
        free(obj_labels[ind]);
        obj_labels[ind] = (char*)malloc(sizeof(char)*vs);
      }
      for (int j = 0; j < vs; ++j) {
        obj_labels[ind][j] = val[j];
      }
      obj_labels[ind][vs - 1] = 0;
    }
  }

  size_t get_offspring_num() { return offspring_num; }
  //size_t get_survivors_num() { return survivors_num; }
  _uint get_n_bits() { return N_BITS; }
  _uint get_n_objs() { return N_OBJS; }
  ArgStore& get_args() { return args; }
};
    
}

#endif //GENETICS_H
