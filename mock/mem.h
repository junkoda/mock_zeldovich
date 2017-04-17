#ifndef MEM_H
#define MEM_H 1

#include <string>

class Mem {
 public:
  Mem(const char name[]);
  Mem(const char name[], const size_t size);
  ~Mem();

  void reserve(size_t size, char const * const msg);
  void alloc();
  void* use_from_zero(size_t size);
  void* use_remaining(size_t size);
    
  std::string* name;
  void* buf;
  size_t size_alloc;
  size_t size_using;
};

#endif
