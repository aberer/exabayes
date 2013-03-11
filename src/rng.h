#ifndef _RNG_H
#define _RNG_H

/* TODO
   
   we could mess around here a lot with 64-bit or the even cooler
   alternative to threefry (aesni). For the time being that's hardly
   worth it, threefry is already much better than the default RNG.
 */

#include <Random123/threefry.h>
#include <Random123/u01.h>

#define exa_rand(c,k) threefry2x32(c,k)

typedef threefry2x32_key_t randKey_t; 
typedef threefry2x32_ctr_t randCtr_t; 

#endif 
