/*===============================================================================/
/  Random number generator code for Monte Carlo alloy sim.                       /
/  Uses the Mersenne Twister algorithm.                                          /
/                                                                                /
/  C. D. Woodgate, Bristol                                                 2025  /
/                                                                                /
/ This file includes code derived from the original Mersenne Twister             /
/ code by Makoto Matsumoto and Takuji Nishimura.                                 /
/                                                                                /
/ The original code can be found at the following link:                          /
/ http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html         /
/                                                                                /
/ The code is subject to their original copyright notice copied below:           /
/                                                                                /
/ COPYRIGHT NOTICE FOR MERSENNE TWISTER CODE                                     /
/ Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,              /
/ All rights reserved.                                                           /
/                                                                                /
/ Redistribution and use in source and binary forms, with or without             /
/ modification, are permitted provided that the following conditions             /
/ are met:                                                                       /
/                                                                                /
/ 1. Redistributions of source code must retain the above copyright              /
/ notice, this list of conditions and the following disclaimer.                  /
/                                                                                /
/ 2. Redistributions in binary form must reproduce the above copyright           /
/ notice, this list of conditions and the following disclaimer in the            /
/ documentation and/or other materials provided with the distribution.           /
/                                                                                /
/ 3. The names of its contributors may not be used to endorse or promote         /
/ products derived from this software without specific prior written             /
/ permission.                                                                    /
/                                                                                /
/ THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS            /
/ "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT              /
/ LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR          /
/ A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR /
/ CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,          /
/ EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,            /
/ PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR             /
/ PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF         /
/ LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING           /
/ NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS             /
/ SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                   /
/                                                                                /
/===============================================================================*/

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

#include <time.h>
#include <stdio.h>

unsigned long mt[N];   /* the array for the state vector  */
int mti=N+1;           /* mti==N+1 means mt[N] is not initialized */


/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{

  mti=N+1;
  mt[0]= s & 0xffffffffUL;
  for (mti=1; mti<N; mti++) {
    mt[mti] =
      (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
    /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
    /* In the previous versions, MSBs of the seed affect   */
    /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
    mt[mti] &= 0xffffffffUL;
    /* for >32 bit machines */
  }
}


/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
  unsigned long y;
  /* static unsigned long mag01[2]={0x0UL, MATRIX_A}; */
  unsigned long mag01[2]={0x0UL, MATRIX_A};
  /* mag01[x] = x * MATRIX_A  for x=0,1 */
  
  if (mti >= N) { /* generate N words at one time */
    int kk;
    
    if (mti == N+1)   /* if init_genrand() has not been called, */
      init_genrand(5489UL); /* a default initial seed is used */
    
    for (kk=0;kk<N-M;kk++) {
      y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
      mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    for (;kk<N-1;kk++) {
      y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
      mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
    mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];
    
    mti = 0;
  }
  
  y = mt[mti++];
  
  /* Tempering */
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);
  
  return y;
}

/* generates a random number on [0,1)-real-interval */
double genrand(void)
{
  return genrand_int32()*(1.0/4294967296.0);
  /* divided by 2^32 */
}

unsigned long f90_init_genrand(int seedtime, int my_rank){
  unsigned long seed;
  /*printf(" Recieved seedtime = %d\n", seedtime);*/
  if(seedtime) {
    seed = time(NULL);
    seed += 11*my_rank;
    if(seed%2 ==0) { seed += 1;}
    /*printf(" # Using time-based random seed %ld\n",seed);*/
  }
  else {
    seed = 110179+my_rank;
  }
  init_genrand(seed);
  return seed;
}
