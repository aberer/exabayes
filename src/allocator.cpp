// #include <new>			// 

// // #include "config.h"
// // #include "common.h"


// using namespace std; 


// // we do not actually need any of this as long as we do not use memory saving features 

// #if 0 


// #include <iostream>

// #if HAVE_PLL != 0
// extern "C"
// {
//   void *rax_malloc(size_t size);
//   void rax_free(void *p);
// }
// #endif

// #if 0 

// // TODO only compile when used with pll 


// // #define REPORT_MEMORY

// #if  HAVE_PLL != 0  &&  defined(_GLIBCXX_THROW) // BAD BAD hack 



// static void* allocate(size_t s)
// {
// #ifdef REPORT_MEMORY
//   cout  << "allocated " << s << " bytes for you "  << endl; 
// #endif
//   void *p = (void*)exa_malloc(s); 
//   return p; 
// }

// static void doFree(void *p)
// {
// #ifdef REPORT_MEMORY
//   cout << "freeing " << p<< endl; 
// #endif
//   exa_free(p) ; 
// }



// void* operator new(size_t s) _GLIBCXX_THROW  (bad_alloc) 
// {
//   void *p = allocate(s); 
//   if(p == NULL)
//     throw new  bad_alloc; 
//   return p ;     
// }

// void* operator new[](size_t s) _GLIBCXX_THROW  (bad_alloc) 
// {
//   void *p = allocate(s); 
//   if(p == NULL)
//     throw new  bad_alloc; 
//   return p ;     
// }


// void* operator new(size_t s, const nothrow_t&) _GLIBCXX_USE_NOEXCEPT
// {
//     void *p = allocate(s); 
//   return p ;     
// }

// void* operator new[](size_t s, const nothrow_t&) _GLIBCXX_USE_NOEXCEPT
// {
//   void *p = allocate(s); 
//   return p ;     
// }


// void operator delete(void* p, const nothrow_t&)     _GLIBCXX_USE_NOEXCEPT
// {
//   doFree(p); 
// }



// void operator delete[](void* p, const nothrow_t&)     _GLIBCXX_USE_NOEXCEPT
// {
//   doFree(p); 
// }


// void operator delete(void* p)  _GLIBCXX_USE_NOEXCEPT
// {
//   doFree(p); 
// }


// void operator delete[](void* p)  _GLIBCXX_USE_NOEXCEPT
// {
//   doFree(p); 
// }


// #endif

// #endif
// #endif
