#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <cstring>
#include <iostream>
#include <random>
#include <math.h>
#include <stdarg.h>
#include <climits>
#include <type_traits>

#ifdef USE_EXCEPTIONS
#include <stdexcept>
#undef USE_CUSTOM_CONTAINERS
#endif

#ifndef USE_CUSTOM_CONTAINERS
#include <vector>
#include <string>
#include <sstream>
#endif

#define PI 3.1415926535897932

#define DEFAULT_STRING_SIZE 8
#define DEFAULT_PRINT_SIZE  24

#define CODE_MISC		5
#define CODE_ARG_INVALID	4
#define CODE_ARG_RANGE  	3
#define CODE_MATH_ERROR  	2
#define CODE_WARN		1
#define CODE_NONE		0

namespace Genetics {

template <bool condition, typename T>
struct enable_if_c {};

template <typename T>
struct enable_if_c<true, T> { typedef T type; };

/*template <typename Ts>
using void_t = void;



template <typename T>
struct has_average_fitness {
private:
  struct true_val{ char x[1]; };
  struct false_val{ char x[2]; };

  template <typename U>
  static true_val test( const U& a, decltype(&U::average_fitness) p = 0 ) {}
  //template <typename U>
  //static false_val test( ... ) {}
public:
  static const bool value = ( sizeof( test<T>(NULL) ) == sizeof( true_val ) );
};*/

// Primary template with a static assertion
// for a meaningful error message
// if it ever gets instantiated.
// We could leave it undefined if we didn't care.

template<typename, typename T>
struct has_average_fitness {
    static_assert(
        std::integral_constant<T, false>::value,
        "Second template parameter needs to be of function type.");
};

// specialization that does the checking

template<typename C, typename Ret, typename... Args>
struct has_average_fitness<C, Ret(Args...)> {
private:
    template<typename T>
    static constexpr auto check(T*)
    -> typename
        std::is_same<
            decltype( std::declval<T>().average_fitness( std::declval<Args>()... ) ),
            Ret    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        >::type;  // attempt to call it and see if the return type is correct

    template<typename>
    static constexpr std::false_type check(...);

    typedef decltype(check<C>(0)) type;

public:
    static constexpr bool value = type::value;
};

//end template checkers

typedef unsigned int _uint;
typedef unsigned char _uchar;

typedef double (*FIT_FUNC_PTR)(size_t, char*);

typedef enum Type{t_int, t_uint, t_real, t_bitstream, t_terminator} Type;

unsigned long encodeGray (unsigned long number);
unsigned long decodeGray (unsigned long code);

size_t getlen_and_clean (char* str);

unsigned int nChoosek( unsigned int n, unsigned int k );

_uint factorial(_uint n);

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
  size_t size() {
    return buf_size - 1;//minus the null terminating character
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
void error (int fatal, String msg, ...);
String read_number(String::iterator* it);

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

char* clean_c_str(char* src);

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

}

#endif
