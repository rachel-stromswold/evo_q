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
#include <string>
#include <sstream>
//#endif

#define PI 3.1415926535897932

#define DEFAULT_STRING_SIZE 8
#define DEFAULT_PRINT_SIZE  24

#define CODE_ERROR  2
#define CODE_FATAL  1
#define CODE_WARN   0

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
#define WAIT_CON      1
#define VERBOSE       2

//POPULATION_H
#define BUF_SIZE      50
#define DEF_SORT_PARAM    -3

namespace Genetics {
//UTIL_H
    
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
    
    //draw k elements from the integer range from 0 to n useful for sampling from arrays
    class SampleDraw {
    private:
        _uint n_;
        _uint k_;
        std::uniform_int_distribution<_uint> dist;
        
    public:
        SampleDraw(_uint p_n, _uint p_k) : dist(0, nChoosek(p_n-1, p_k-1)*factorial(p_k-1)) {
            n_ = p_n;
            k_ = p_k;
            if (k_ > n_) {
                error(1, "Cannot initialize sample draw with more samples than the population size.");
            }
        }
        
        void reset() { dist.reset(); }
        _uint n() { return n_; }
        _uint k() { return k_; }
        
        template <class Generator>
        std::vector<_uint> operator()(Generator& g) {
            std::vector<_uint> ret(k_, 0);
            for (_uint i = 0; i < k_; ++i) {
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
        std::bernoulli_distribution* bern;
        unsigned int short_bin_n = 0;
        //END REMOVAL CANDIDATE
        
        ArgStore(ArgStore const&){};
        ArgStore& operator=(ArgStore const& arg) {return *this;};
        
        unsigned flags;
        size_t pop_size;
        size_t breed_pop_size;
        size_t num_gens;
        int num_crossovers;
        double init_coup_var;
        double init_coup_mean;
        double mutate_prob;
        double hypermutation_threshold;
        std::string out_fname;
        
        bool activate = true;
        int seed = 0;
        
    public:
        ArgStore();
        
        void initialize_from_args(size_t argc, char** argv);
        void initialize_from_file(const char* fname);
        void initialize();
        void read_file(std::string fname);
        std::mt19937& get_generator() { return generator; }
        //REMOVAL CANDIDATE
        unsigned int sample_binomial(unsigned int n);
        bool sample_bernoulli();
        void print_data();
        //END REMOVAL CANDIDATE
        
        size_t get_pop_size()    { return pop_size; }
        size_t get_survivors()    { return breed_pop_size; }
        size_t get_num_gens()     { return num_gens; }
        int get_num_crossovers()     { return num_crossovers; }
        double get_init_coup_var()    { return init_coup_var; }
        double get_init_coup_mean() { return init_coup_mean; }
        double get_mutate_prob()    { return mutate_prob; }
        double get_hypermutation_threshold()    { return hypermutation_threshold; }
        bool wait_for_con()        { return flags & WAIT_CON; }
        bool verbose()        { return flags & VERBOSE; }
        std::string get_out_fname() { return out_fname; }
    };
    
//CHROMOSOME_H
    
    class Chromosome {
    private:
        _uint N_BITS;
        size_t N_BYTES;
        size_t N;
        
    protected:
        static const _uint bin_size = sizeof(unsigned long)*8;
        //  unsigned long genes[(N_BYTES+sizeof(unsigned long)-1)/sizeof(unsigned long)];
        Vector<unsigned long> genes;
        size_t getBitStream (size_t n, size_t k, size_t x);
        
    public:
        Chromosome(_uint pn_bits);
        Chromosome(_uint pn_bits, Chromosome* o);
        Chromosome(Chromosome& other);
        Chromosome(Chromosome&& other);
        void exchange(Chromosome* other, size_t k);
        void exchange_uniform(ArgStore* args, Chromosome* other);
        
        unsigned int get_N() { return N; }
        unsigned int get_n_bits() { return N_BITS; }
        Chromosome& operator=(Chromosome& other);
        
        void reset();
        
        unsigned char operator[](unsigned int i);
        //randomly mutate each bit in the gene
        void mutate(ArgStore* args);
        void slow_mutate(ArgStore* args);
        //set the gene to a new completely random value
        void randomize(ArgStore* args);
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
    
    class Organism;
    
    struct Result {
        size_t index;
        Vector<double> fit_vals;
        String misc_str;
    };
    
    class Problem {
    public:
        _uint N_BITS, N_PARAMS, N_OBJS;
        PhenotypeMap map;
        Vector<Result> result_list;
        //  Problem() {std::cout << "Initializing problem...\n"; }
        Problem(unsigned n_bits, unsigned n_params, int n_objs) : N_BITS(n_bits), N_PARAMS(n_params), N_OBJS(n_objs), map(n_bits) {}
        
        virtual void evaluate_fitness(Organism* org) {}
        virtual void evaluate_fitness_async(size_t index, Chromosome genes) {}
    };
    
    class Organism {
    private:
        _uint N_BITS;
        _uint N_OBJS;
        
        char output_stream[BUF_SIZE];
        Vector<double> fitness;
        size_t output_len;
        PhenotypeMap* al;
        
    protected:
        Chromosome* genes;
        size_t n_nodes;
        
    public:
        //this is not used internally, but can be set when evaluating the fitness
        String misc_data;
        double coupling_range;
        double coupling_prec;
        
        int n_dominations;
        int rank;
        double distance;
        
        Organism(int N_BITS, int N_OBJS, PhenotypeMap* p_al);
        Organism(int N_BITS, int N_OBJS, Chromosome p_genes, PhenotypeMap* p_al);
        Organism(const Organism &obj);
        Organism(Organism&& obj);
        ~Organism();
        
        Organism& operator=(Organism& obj);
        
        std::vector<Organism*> breed(ArgStore* args, Organism* par1);
        void reset();
        void randomize(ArgStore* args);
        void randomize(ArgStore* args, Organism* orgtmp);
        
        void evaluate_fitness(Problem* prob);
        
        double get_fitness(_uint i = 0);
        void set_fitness(double val);
        void set_fitness(_uint i, double val);
        
        void set_int(_uint i, int value);
        void set_uint(_uint i, int value);
        void set_real(_uint i, double value);
        double read_real(_uint i);
        int read_int(_uint i);
        bool dominates(Organism* other);
        String get_chromosome_string(_uint i) { return genes->get_string(al, i); }
        char* get_output_stream() { return output_stream; }
        size_t get_output_len() {return output_len; }
        int get_rank() { return rank; }
        _uint get_n_bits() { return N_BITS; }
        _uint get_n_params() { return al->get_num_params(); }
        _uint get_n_objs() { return N_BITS; }
    };
    
//POPULATION_H
    
    class Population {
    private:
        _uint N_BITS;
        _uint N_OBJS;
    protected:
        size_t sort_org_calls = 0;
        size_t carryover_num;//How many of the best individuals carry over to the next generation
        //OWNED POINTERS
        double* max_fitness = NULL;
        double* min_fitness = NULL;
        //EXTERNALLY MANAGED POINTERS
        PhenotypeMap* map = NULL;
        
        ArgStore args;
        //all offspring from the previous generation
        size_t offspring_num;
        std::vector<std::shared_ptr<Organism>> offspring;
        std::vector<std::shared_ptr<Organism>> old_gen;
        //which offspring will survive to enter the next breeding round
        size_t survivors_num;
        std::vector<std::shared_ptr<Organism>> survivors;
        //guarantee that the best organism appears in the next generation
        size_t best_organism_ind;
        std::shared_ptr<Organism> best_organism;
        //labels for generating data output
        Vector<String> var_labels;
        Vector<String> obj_labels;
        
        //cull in place is slightly faster but less accurate than the standard cull method
        bool cull_in_place();
        //cull first sorts the organisms and selects them based on the ratio of their relative fitness to the total relative fitness
        bool cull();
        void breed();
        
    public:
        Population(_uint pn_bits, _uint pn_objs, PhenotypeMap* p_map, String conf_fname = "");
        Population(_uint pn_bits, _uint pn_objs, Organism* tmplt, PhenotypeMap* p_map, String conf_fname = "");
        ~Population();
        Population(Population& o);
        Population& operator=(Population& o);
        Population(Population&& o);
        
        void evaluate(Problem* prob);
        void evaluate_async(Problem* prob);
        bool iterate();
        
        std::shared_ptr<Organism> get_best_organism() { return old_gen[best_organism_ind]; }
        std::shared_ptr<Organism> get_organism(size_t i) {
            if (i > old_gen.size())
                error(1, "Attempt to access invalid index %d when the maximum allowed is %d.", i, old_gen.size());
            return old_gen[i];
        }
        std::shared_ptr<Organism> get_child(size_t i) {
            if (i > offspring.size())
                error(1, "Attempt to access invalid index %d when the maximum allowed is %d.", i, offspring.size());
            return offspring[i];
        }
        
        void run(Problem* prob);
        void swap_orgs(int i, int j);
        int partition(_uint ind, std::vector<std::shared_ptr<Organism>>* work_arr, int s, int e);
        void sort_orgs(unsigned int ind, std::vector<std::shared_ptr<Organism>>* arr, int s = DEF_SORT_PARAM, int e = DEF_SORT_PARAM);
        Vector<String> get_header();
        Vector<String> get_pop_data();
        
        double get_min_fitness(_uint i = 0) { return min_fitness[i]; }
        double get_max_fitness(_uint i = 0) { return max_fitness[i]; }
        
        void set_var_label(_uint ind, String val) { var_labels[ind] = val; }
        void set_obj_label(_uint ind, String val) { obj_labels[ind] = val; }
        
        size_t get_offspring_num() { return offspring_num; }
        size_t get_survivors_num() { return survivors_num; }
        inline _uint get_n_bits() { return N_BITS; }
        inline _uint get_n_objs() { return N_OBJS; }
        ArgStore& get_args() { return args; }
    };
    
    class Population_NSGAII : public Population {
    private:
        size_t survivors_num;
        //the ngsa alternative to elitism
        std::vector<std::vector<std::shared_ptr<Organism>>> pareto_fronts;
        void hypermutate();
        
    public:
        Population_NSGAII(_uint pn_bits, _uint pn_objs, PhenotypeMap* p_map, String conf_fname) :
        Population(pn_bits, pn_objs, p_map, conf_fname) {}
        Population_NSGAII(_uint pn_bits, _uint pn_objs, Organism* tmplt, PhenotypeMap* p_map, String conf_fname) :
        Population(pn_bits, pn_objs, tmplt, p_map, conf_fname) {}
        ~Population_NSGAII();
        
        void evaluate(Problem* prob);
        PhenotypeMap* get_map() { return this->map; }
        
        //cull in place is slightly faster but less accurate than the standard cull method
        bool cull_in_place();
        //cull first sorts the organisms and selects them based on the ratio of their relative fitness to the total relative fitness
        void cull();
        void breed();
        bool iterate();
        
        _uint get_n_pareto_fronts() {
            return pareto_fronts.size();
        }
        std::vector<std::shared_ptr<Organism>> get_pareto_front(_uint i) {
            if (pareto_fronts.size() == 0) { error(1, "The population has not yet been evaluated."); }
            return pareto_fronts[i];
        }
        
        Vector<String> get_header();
        Vector<String> get_pop_data();
        
        std::shared_ptr<Organism> get_organism(size_t i) { return (i < offspring_num) ? this->old_gen[i] : this->offspring[i - offspring_num]; }
        
        size_t get_offspring_num() { return this->offspring_num; }
        size_t get_survivors_num() { return this->survivors_num; }
    };
    
}

#endif //GENETICS_H
