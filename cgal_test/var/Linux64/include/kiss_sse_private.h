/*---------------------------------------------------------------------------

Copyright (c) Dave Ritchie, 2008

  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
  * Neither the author nor the names of any contributors may be used to
    endorse or promote products derived from this software without
    specific prior written permission.
  
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.

-----------------------------------------------------------------------------*/

#ifndef KISS_SSE_PRIVATE_H
#define KISS_SSE_PRIVATE_H

#ifdef __cplusplus
extern "C" {
#endif

/*---------------------------------------------------------------------------*/

typedef float v4sf __attribute__ ((vector_size(16)));  /* mode(V4SF) */

typedef union _v4f {
   v4sf  vec         __attribute__ ((aligned(16)));
   float val[4]      __attribute__ ((aligned(16)));
} v4f; 

typedef struct _cv4f {
   union {
      v4sf r         __attribute__ ((aligned(16)));
      float rval[4]  __attribute__ ((aligned(16)));

   };
   union {
      v4sf i         __attribute__ ((aligned(16)));
      float ival[4]  __attribute__ ((aligned(16)));

   };
} cv4f;

typedef cv4f kiss_fft_cpx;


/*---------------------------------------------------------------------------*/

static __inline__ cv4f c_mul_cv4f(cv4f p, cv4f q) {  

   cv4f result;

   v4sf ac = __builtin_ia32_mulps(p.r, q.r);
   v4sf bd = __builtin_ia32_mulps(p.i, q.i);
   v4sf bc = __builtin_ia32_mulps(p.i, q.r);
   v4sf ad = __builtin_ia32_mulps(p.r, q.i);

   result.r = __builtin_ia32_subps(ac, bd);
   result.i = __builtin_ia32_addps(bc, ad);

   return result;
}

/*---------------------------------------------------------------------------*/

static __inline__ cv4f c_mulpi_cv4f(cv4f p) {  

   int   i;
   cv4f result;

   for (i=0; i<4; i++) result.rval[i] = p.ival[i];

// result.r = - p.i;
   result.i =  p.r;

   return result;
}

/*---------------------------------------------------------------------------*/

static __inline__ cv4f c_mulmi_cv4f(cv4f p) {  

   int   i;
   cv4f result;

   result.r =   p.i;
   result.i = - p.r;
   for (i=0; i<4; i++) result.ival[i] = -p.rval[i];

   return result;
}

/*---------------------------------------------------------------------------*/

static __inline__ cv4f c_div2_cv4f(cv4f p) {  

   cv4f result;

   v4sf half = {0.5, 0.5, 0.5, 0.5};

   result.r = __builtin_ia32_mulps(p.r, half);
   result.i = __builtin_ia32_mulps(p.i, half);

   return result;
}

/*---------------------------------------------------------------------------*/

static __inline__ cv4f c_sub_cv4f(cv4f p, cv4f q) {  

   cv4f result;

   result.r = __builtin_ia32_subps(p.r, q.r);
   result.i = __builtin_ia32_subps(p.i, q.i);

   return result;
}

/*---------------------------------------------------------------------------*/

static __inline__ cv4f c_add_cv4f(cv4f p, cv4f q) {  

   cv4f result;

   result.r = __builtin_ia32_addps(p.r, q.r);
   result.i = __builtin_ia32_addps(p.i, q.i);

   return result;
}

/*---------------------------------------------------------------------------*/

static __inline__ cv4f c_conj_cv4f(cv4f p) {  

   int   i;
   cv4f result;

   result.r =   p.r;
// result.i = - p.i;

   for (i=0; i<4; i++) result.ival[i] = -p.ival[i];


   return result;
}

/*---------------------------------------------------------------------------*/

static __inline__ cv4f c_set1_cv4f(float re, float im) {  

   int  i;
   cv4f result;

   for (i=0; i<4; i++) {
      result.rval[i] = re;
      result.ival[i] = im;
   }
/*
   result.r = (v4sf) __builtin_ia32_setps1(re);
   result.i = (v4sf) __builtin_ia32_setps1(im);
*/

   return result;
}

/*---------------------------------------------------------------------------*/
#ifdef __cplusplus
}
#endif

#endif /* KISS_SSE_PRIVATE_H */
/*---------------------------------------------------------------------------*/
