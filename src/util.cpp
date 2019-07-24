#include "util.h"

namespace Genetics {

#ifndef USE_EXCEPTIONS
String global_err;
_uchar global_code = CODE_NONE;
#endif
  
void error (int code_class, String msg, ...) {
#ifndef USE_EXCEPTIONS
  global_code = code_class;
#endif
  va_list vl;
  va_start(vl, msg);
  char* message = strdup(msg.c_str());
  char* last_str = message;
#ifdef USE_CUSTOM_CONTAINERS
  String iss;
#else
  std::stringstream iss;
#endif
  if (code_class == CODE_WARN) {
    iss << "Warning: ";
  } else if (code_class == CODE_MATH_ERROR) {
    iss << "Math Error: ";
  } else if (code_class == CODE_ARG_RANGE) {
    iss << "Out of Range Error: ";
  } else if (code_class == CODE_ARG_INVALID) {
    iss << "Invalid Argument Error: ";
  } else {
    iss << "Error: ";
  }
  for (_uint i = 0; message[i] != 0; ++i) {
    if (message[i] == '%') {
      message[i] = 0;
      iss << last_str;
      //fprintf(stderr, "%s", last_str);
      i++;
      if (message[i] == 's') {
        char* tmp = va_arg(vl, char*);
	iss << tmp;
        //fprintf(stderr, "%s", tmp);
      } else if (message[i] == 'd' || message[i] == 'i') {
        int tmp = va_arg(vl, int);
	iss << tmp;
        //fprintf(stderr, "%d", tmp);
      } else if (message[i] == 'o') {
        unsigned int tmp = va_arg(vl, unsigned int);
	iss << "0o";
        //fprintf(stderr, "0o");
        for (int k = sizeof(unsigned int)*8/3; k >= 0; --k) {
	  iss << ( (tmp >> (k*3)) & 7 );
          //fprintf(stderr, "%d", (tmp >> (k*3)) & 7);
        }
      } else if (message[i] == 'u') {
        unsigned int tmp = va_arg(vl, unsigned int);
	iss << tmp;
        //fprintf(stderr, "%d", tmp);
      } else if (message[i] == 'x' || message[i] == 'X') {
        unsigned int tmp = va_arg(vl, unsigned int);
	iss << "0x";
        //fprintf(stderr, "0x");
        for (int k = sizeof(unsigned int)*2 - 1; k >= 0; --k) {
	  iss << ( (tmp >> (k*4)) & 15 );
          //fprintf(stderr, "%d", (tmp >> (k*4)) & 15);
        }
      } else if (message[i] == 'f' || message[i] == 'F'
              || message[i] == 'g' || message[i] == 'G'
              || message[i] == 'e' || message[i] == 'E') {
        double tmp = va_arg(vl, double);
	iss << tmp;
        //fprintf(stderr, "%g", tmp);
      } else if (message[i] == 'a' || message[i] == 'A') {
        double tmp = va_arg(vl, double);
	iss << tmp;
        //fprintf(stderr, "%a", tmp);
      } else if (message[i] == '%') {
	iss << '%';
        //fprintf(stderr, "%%");
      }
      i++;
      last_str = &(message[i]);
    }
  }
  iss << last_str << "\n";
  //fprintf(stderr, "%s", last_str);
  va_end(vl);
  //fprintf(stderr, "\n");
  free(message);

#ifdef USE_EXCEPTIONS
  if (code_class == CODE_WARN) {
    std::cerr << iss.str();
  } else if (code_class == CODE_MATH_ERROR) {
    throw std::domain_error(iss.str());
  } else if (code_class == CODE_ARG_RANGE) {
    throw std::out_of_range(iss.str());
  } else if (code_class == CODE_ARG_INVALID) {
    throw std::invalid_argument(iss.str());
  } else {
    throw std::runtime_error(iss.str());
  }
#else
#ifdef USE_CUSTOM_CONTAINERS
  global_err = iss;
#endif
#endif
}

#ifndef USE_EXCEPTIONS
bool has_error() {
  return global_code != CODE_NONE;
}

String get_error() {
  if (global_code) {
    String tmp = global_err;
    global_err = "";
    global_code = CODE_NONE;
    return tmp;
  }

  String ret("");
  return ret;
}
#endif

char* clean_c_str(char* str) {
  size_t l = strlen(str);
  while (l > 0 && (str[l - 1] == ' ' || str[l - 1] == '\t' || str[l - 1] == '\n')) {
    str[l - 1] = 0;
    --l;
  }
  size_t i = 0;
  while (l > 0 && str[i] != 0 && (str[i] == ' ' || str[i] == '\t')) {
    ++i;
    --l;
  }
  return str + i;
}

unsigned long decodeGray (unsigned long code) {
  unsigned long mask = code >> 1;
  unsigned long ret = code;

  while (mask != 0) {
    ret = ret ^ mask;
    mask = mask >> 1;
  }
  return ret;
}

unsigned long encodeGray (unsigned long number) {
  return number ^ (number >> 1);
}

size_t getlen_and_clean (char* str) {
  size_t i = 0;
  while (str[i] != 0) {
    if (str[i] == '\n') {
      str[i] = 0;
      return i;
    }
    i++;
  }
  return i;
}

unsigned int nChoosek( unsigned int n, unsigned int k ) {
  if (k > n) return 0;
  if (k * 2 > n) k = n-k;
  if (k == 0) return 1;

  unsigned int result = n;
  for( unsigned int i = 2; i <= k; ++i ) {
    result *= (n-i+1);
    result /= i;
  }
  return result;
}

_uint factorial(_uint n) {
  _uint ret = 1;
  for (_uint i = 1; i <= n; ++i) {
    ret *= i;
  }
  return ret;
}

String read_number(String::iterator* it) {
  char t = *(*it);
  if (t < '0' || t > '9') {
    ++(*it);
  }
  std::string ret = "";
  while ( t >= '0' && t <= '9') {
    ret += t;
    ++(*it);
    t = *(*it);
  }
  return ret;
}

SampleDraw::SampleDraw(_uint p_n, _uint p_k, bool replace) : dist(0, nChoosek(p_n-1, p_k-1)*factorial(p_k-1)) {
  this->replace = replace;
  n_ = p_n;
  k_ = p_k;
  if (k_ > n_) {
    error(1, "Cannot initialize sample draw with more samples than the population size.");
  }
}

}
