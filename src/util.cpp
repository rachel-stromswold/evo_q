#include "util.h"

namespace Genetics {
  
void error (int fatal, std::string msg, ...) {
  va_list vl;
  va_start(vl, msg);
  char* message = strdup(msg.c_str());
  char* last_str = message;
  if (fatal == CODE_FATAL) {
    fprintf(stderr, "Fatal Error: ");
  } else if (fatal == CODE_WARN) {
    fprintf(stderr, "Warning: ");
  } else if (fatal == CODE_ERROR) {
    fprintf(stderr, "Error: ");
    fatal = CODE_WARN;
  }
  for (_uint i = 0; message[i] != 0; ++i) {
    if (message[i] == '%') {
      message[i] = 0;
      fprintf(stderr, "%s", last_str);
      i++;
      if (message[i] == 's') {
        char* tmp = va_arg(vl, char*);
        fprintf(stderr, "%s", tmp);
      } else if (message[i] == 'd' || message[i] == 'i') {
        int tmp = va_arg(vl, int);
        fprintf(stderr, "%d", tmp);
      } else if (message[i] == 'o') {
        unsigned int tmp = va_arg(vl, unsigned int);
        fprintf(stderr, "0o");
        for (int k = sizeof(unsigned int)*8/3; k >= 0; --k) {
          fprintf(stderr, "%d", (tmp >> (k*3)) & 7);
        }
      } else if (message[i] == 'u') {
        unsigned int tmp = va_arg(vl, unsigned int);
        fprintf(stderr, "%d", tmp);
      } else if (message[i] == 'x' || message[i] == 'X') {
        unsigned int tmp = va_arg(vl, unsigned int);
        fprintf(stderr, "0x");
        for (int k = sizeof(unsigned int)*2 - 1; k >= 0; --k) {
          fprintf(stderr, "%d", (tmp >> (k*4)) & 15);
        }
      } else if (message[i] == 'f' || message[i] == 'F'
              || message[i] == 'g' || message[i] == 'G'
              || message[i] == 'e' || message[i] == 'E') {
        double tmp = va_arg(vl, double);
        fprintf(stderr, "%g", tmp);
      } else if (message[i] == 'a' || message[i] == 'A') {
        double tmp = va_arg(vl, double);
        fprintf(stderr, "%a", tmp);
      } else if (message[i] == '%') {
        fprintf(stderr, "%%");
      }
      i++;
      last_str = &(message[i]);
    }
    if (message[i] == '\\') {
      message[i] = 0;
      fprintf(stderr, "%s", last_str);
      char tmp[5];
      tmp[0] = '\\';
      if (message[i+1] >= '0' && message[i+1] <= '9') {
        int k = 1;
        while (k < 4 && message[i+k] >= '0' && message[i+k] <= '9') {
          tmp[k] = message[i+k];
          k++;
        }
        tmp[k] = 0;
        i += k;
      } else {
        tmp[1] = message[i+1];
        tmp[2] = 0;
        i += 2;
      }
      fprintf(stderr, "%s", tmp);
      last_str = &(message[i]);
    }
  }
  fprintf(stderr, "%s", last_str);
  va_end(vl);
  fprintf(stderr, "\n");
  free(message);
  if (fatal) {
    exit(1);
  }
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

std::string read_number(std::string::iterator* it) {
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

}
