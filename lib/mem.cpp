//
// Memory (RAM) management module
//

#include <string>
#include <stdlib.h>
#include <assert.h>
#include "config.h"
#include "msg.h"
#include "error.h"
#include "util.h"
#include "mem.h"
#include "fft.h"

//
// step 1: Mem* mem= new Mem("name");
// step 2: mem->mem_reserve(mem, size1, "usage1");
//         mem->mem_reserve(mem, size2, "usage2");
// step 3: mem->mem_alloc();
//

using namespace std;

Mem::Mem(const char name_[]) :
  buf(0), size_alloc(0), size_using(0)
{
  name= new string(name_);

  msg_printf(msg_verbose, "Memory %s initilised.\n", name->c_str());
}

Mem::Mem(const char name_[], const size_t size) :
  buf(0), size_alloc(0), size_using(0)
{
  // Raises: MemoryError
  name= new string(name_);
  
  reserve(size, 0);
  alloc();
}

Mem::~Mem()
{
  free(buf);
  delete name;
}

void Mem::reserve(size_t size, char const * const msg)
{
  if(size % ALGN != 0)
    size += ALGN - (size % ALGN);
  assert(size % ALGN == 0);
  
  if(size > size_using)
    size_using= size;   // this is the amount going to be allocated

  if(msg)
    msg_printf(msg_info, "%s requested %lu MB for %s\n",
	       msg, mbytes(size), name->c_str());
}

void Mem::alloc()
{
  free(buf);

  buf= fft_malloc(size_using);

  if(buf == 0) {
    msg_printf(msg_fatal, "Error: Unable to allocate %lu MB for %s\n",
	       mbytes(size_using), name->c_str());
    throw MemoryError();
  }
  else
    msg_printf(msg_info, "%lu MB allocated for mem %s\n",
	       mbytes(size_using), name->c_str());

  size_alloc= size_using;
  size_using= 0;
}


void* Mem::use_from_zero(size_t size)
{
  size= size_align(size); 

  if(size > size_alloc) {
    msg_printf(msg_error,
	      "Error: Unable to use %lu MB in Mem %s (only %lu MB allocated)\n",
	       mbytes(size), name->c_str(), mbytes(size_alloc));
    throw MemoryError();
  }

  size_using= size;
  return buf;
}

void* Mem::use_remaining(size_t size)
{
  size= size_align(size);

  if(size + size_using > size_alloc) {
    msg_printf(msg_error,
	       "Error: Unable to use %lu MB in Mem %s; %lu MB allocated, "
	       "%lu remaining.\n",
	       mbytes(size), name->c_str(), mbytes(size_alloc),
	       mbytes(size_alloc - size_using));
    throw MemoryError();
  }
  
  assert(size_using % sizeof(float) == 0);
  float* const p= (float*) buf;
  const size_t n= size_using / sizeof(float);

  size_using += size;
  msg_printf(msg_verbose, "Using %lu of %lu in memory %s\n",
	     size_using, size_alloc, name->c_str());

  return p+n;
}
