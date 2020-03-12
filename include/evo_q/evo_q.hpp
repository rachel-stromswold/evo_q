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
    double init_param_var;
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
    double get_init_param_var()			{ return init_param_var; }
    void set_init_param_var(double x)		{ init_param_var = x; }
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
  Vector<double> get_real_vector(PhenotypeMap* al);
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

template <typename FitType>
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

template <typename FitType>
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

    double var = args.get_init_param_var();
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
  //_uint get_n_evaluations() { return n_evaluations; }
  //average fitness with another organism if they both have the same genotype
  /*void average_fitness(Organism* other);
  void copy_fitness_data(Organism* other);*/

  void set_int(_uint i, int value) { genes.set_to_int(al.get(), i, value); }
  void set_uint(_uint i, _uint value) { genes.set_to_ulong(al.get(), i, value); }
  void set_real(_uint i, double value) { genes.set_to_num(al.get(), i, value); }
  double read_real(_uint i) { return genes.gene_to_num(al.get(), i); }
  int read_int(_uint i) { return genes.gene_to_int(al.get(), i); }
  _uint read_uint(_uint i) { return genes.gene_to_ulong(al.get(), i); }
  Vector<double> read_real_vector() { return genes.get_real_vector(al.get()); }

  //bool dominates(Organism* other);
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
  //int get_rank() { return rank; }
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
  //compare should return >0 in the case where fitness(a) > fitness(b), <0 in the case where fitness(a) < fitness(b) or 0 in cases where comparison is equal or ill-defined
  static int compare(FitType& a, FitType& b) {
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
  static int compare(FitType& a, FitType& b) {
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
  int partition(_uint fit_ind, std::vector<std::shared_ptr< Organism<FitType> >>& work_arr, int s, int e);
public:
  static const bool use_offspring=false;
  void sort_orgs(unsigned int fit_ind, std::vector<std::shared_ptr< Organism<FitType> >>& work_arr, int s = DEF_SORT_PARAM, int e = DEF_SORT_PARAM);
  Selector();
  //virtual OrganismPair select(Population& pop);
  virtual ~Selector() = default;
  //returns >=1 if a > b, 0 if a = b and -1 if a <= b
  //errors if the two fitness values are not comprable (they do not have the same number of objectives)
  
  virtual Vector<ParentIndSet> select(ArgStore& args, Vector<std::shared_ptr<Organism<FitType>>>& old_gen, Vector<std::shared_ptr<Organism<FitType>>>& offspring) = 0;

  //return the index of the best organism found in the population
  static _uint find_best_organism(Vector<std::shared_ptr<Organism<FitType>>>& orgs, Vector<FitnessStats>& pop_stats);
};

template <class FitType, typename MyComp=Comparator<FitType>>
class TournamentSelector : public Selector<FitType, MyComp> {
public:
  typedef MyComp Comp;
  Vector<ParentIndSet> select(ArgStore& args, Vector<std::shared_ptr<Organism<FitType>>>& old_gen, Vector<std::shared_ptr<Organism<FitType>>>& offspring);
};

template <class FitType, typename MyComp=Comparator<FitType>>
class SurvivalSelector : public Selector<FitType, MyComp> {
public:
  typedef MyComp Comp;
  Vector<ParentIndSet> select(ArgStore& args, Vector<std::shared_ptr<Organism<FitType>>>& old_gen, Vector<std::shared_ptr<Organism<FitType>>>& offspring);
};

template <class FitType, typename MyComp=Comparator<FitType>>
class DominanceTournamentSelector : public TournamentSelector<FitType, MyComp> {
public:
  typedef MyComp Comp;

  //return the index of the best organism found in the population
  static _uint find_best_organism(Vector<std::shared_ptr<Organism<FitType>>>& orgs, Vector<FitnessStats>& pop_stats);
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
  void make_fronts(std::vector<OrgPtr>& cmb_arr);
 
  Vector<ParentIndSet> gen_breed_pairs(ArgStore& args, Vector<OrgPtr>& old_gen, Vector<OrgPtr>& offspring);
public:
  NSGAII_TournamentSelector() : Selector<FitType, NSGAII_Comparator<FitType>>() {
    static_assert( std::is_base_of<MultiFitness, FitType>::value, "FitType must be derived from MultiFitness for the NSGAII selector" );
  }
  
  Vector<ParentIndSet> select(ArgStore& args, Vector<OrgPtr>& old_gen, Vector<OrgPtr>& offspring);
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
  void evaluate_best(Problem<FitType>* prob, double forget_weight=0.0);

protected:
  static_assert( std::is_base_of<Selector<FitType, typename SelectType::Comp>, SelectType>::value, "SelectType must be derived from Selector<FitType, Comp>" );
  SelectType sel;
  typedef std::shared_ptr< Organism<FitType> > OrgPtr;
  size_t carryover_num;//How many of the best individuals carry over to the next generation 
  double penalty_fact;//Used to guarantee that penalized

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
  //bool cull_in_place();

  /**
   * \brief An implementation of simple roulette selection. This function first sorts the organisms and selects them based on the ratio of their relative fitness to the total relative fitness.
   *
   * \returns True if all organisms have the same fitness (results have converged).
   */
  //bool cull();
  //void breed_shuffle();
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
  void apply_penalties(Problem<FitType>* prob);
  void find_best_organism();
  void breed(Vector<ParentIndSet>&& parents);
  
  void calculate_distances();
  void hypermutate();
  void set_best_organism(_uint i, bool force=false, _uint j=0);
  size_t find_first_unpenalized(Problem<FitType>* prob);
  void handle_multiples();

public:
//    Population(_uint pn_bits, _uint pn_objs, PhenotypeMap* p_map, ArgStore p_args);
//    Population(_uint pn_bits, _uint pn_objs, Organism<FitType>* tmplt, PhenotypeMap* p_map, ArgStore p_args);
  Population(_uint pn_bits, _uint pn_objs, std::shared_ptr<PhenotypeMap> p_map, bool latin=true);
  Population(_uint pn_bits, _uint pn_objs, Organism<FitType>* tmplt, std::shared_ptr<PhenotypeMap> p_map);
  Population(_uint pn_bits, _uint pn_objs, std::shared_ptr<PhenotypeMap> p_map, ArgStore p_args, bool latin=true);
  Population(_uint pn_bits, _uint pn_objs, Organism<FitType>* tmplt, std::shared_ptr<PhenotypeMap> p_map, ArgStore p_args);
  void createOrganisms(Organism<FitType>* tmplt, bool latin);
  ~Population();
  Population(Population& o);
  Population& operator=(Population& o);
  Population(Population&& o);

  void set_convergence_type(ConvergenceCriteria* conv);
  void set_penalty_printing(bool val = true) { print_penalties = (val)? 1 : 0; }
  void resize_population(_uint new_size);
  void reset_population_fitness();
#ifdef USE_LIBOMP
  template <typename T=FitType>
  void evaluate_async(Problem<T>* prob);
#endif
  //function for population where organism FitType has a member average_fitness
  template <typename T = FitType> inline
  typename enable_if_c< has_average_fitness<T, void(T&)>::value, void >::type
  evaluate(Problem<T>* prob);
  //function for population where organism FitType does not have a member average_fitness
  template <typename T = FitType> inline
  typename enable_if_c< !has_average_fitness<T, void(T&)>::value, void >::type
  evaluate(Problem<T>* prob);
  //void evaluate(Problem<FitType>* prob) { evaluate_imp(prob, NULL); }
  bool iterate(ConvergenceCriteria* conv = NULL);
  void run(Problem<FitType>* prob);
  std::shared_ptr< Organism<FitType> > get_best_organism(size_t i = 0);

  std::shared_ptr< Organism<FitType> > get_organism(size_t i);
  std::shared_ptr< Organism<FitType> > get_child(size_t i);
  
  Vector<String> get_best_header();
  Vector<String> get_header();
  Vector<String> get_best_data();
  Vector<String> get_pop_data();

  Vector<std::pair<std::shared_ptr<Organism<FitType>>, _uint>> get_species_list(double tolerance=0.1, _uint dimension_threshold=1);

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
