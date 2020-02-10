/* One-sided MPI implementation of Libcaf

Copyright (c) 2012-2016, Sourcery, Inc.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the Sourcery, Inc., nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL SOURCERY, INC., BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  */

#ifndef LIBCAF_GFORTRAN_DESCRIPTOR_H
#define LIBCAF_GFORTRAN_DESCRIPTOR_H

#include "libcaf-version-def.h"

/* GNU Fortran's array descriptor.  Keep in sync with libgfortran.h.  To be
   replaced by TS29113's ISO_Fortran_binding.h with CFI_cdesc_t.  */

enum
{ BT_UNKNOWN = 0, BT_INTEGER, BT_LOGICAL, BT_REAL, BT_COMPLEX,
  BT_DERIVED, BT_CHARACTER, BT_CLASS, BT_PROCEDURE, BT_HOLLERITH, BT_VOID,
  BT_ASSUMED
};

typedef struct descriptor_dimension
{
  ptrdiff_t _stride;
  ptrdiff_t lower_bound;
  ptrdiff_t _ubound;
}
descriptor_dimension;

#ifdef GCC_GE_8
  typedef struct dtype_type
  {
    size_t elem_len;
    int version;
    signed char rank;
    signed char type;
    signed short attribute;
  }
  dtype_type;
#endif

typedef struct gfc_descriptor_t {
  void *base_addr;
  size_t offset;
#ifdef GCC_GE_8
  dtype_type dtype;
  ptrdiff_t span;
#else
  ptrdiff_t dtype;
#endif
  descriptor_dimension dim[];
} gfc_descriptor_t;

#ifdef GCC_GE_8

#define GFC_MAX_DIMENSIONS 15
#define GFC_DTYPE_RANK_MASK 0x0F
#define GFC_DTYPE_TYPE_SHIFT 4
#define GFC_DTYPE_TYPE_MASK 0x70
#define GFC_DTYPE_SIZE_SHIFT 7

#define GFC_DESCRIPTOR_RANK(desc) (desc)->dtype.rank
#define GFC_DESCRIPTOR_TYPE(desc) (desc)->dtype.type
#define GFC_DESCRIPTOR_SIZE(desc) (desc)->dtype.elem_len
#define GFC_DTYPE_TYPE_SIZE(desc) (( ((desc)->dtype.type << GFC_DTYPE_TYPE_SHIFT) \
    | ((desc)->dtype.elem_len << GFC_DTYPE_SIZE_SHIFT) ) & GFC_DTYPE_TYPE_SIZE_MASK)

#else

#define GFC_MAX_DIMENSIONS 7
#define GFC_DTYPE_RANK_MASK 0x07
#define GFC_DTYPE_TYPE_SHIFT 3
#define GFC_DTYPE_TYPE_MASK 0x38
#define GFC_DTYPE_SIZE_SHIFT 6

#define GFC_DESCRIPTOR_RANK(desc) ((desc)->dtype & GFC_DTYPE_RANK_MASK)
#define GFC_DESCRIPTOR_TYPE(desc) (((desc)->dtype & GFC_DTYPE_TYPE_MASK) \
                                   >> GFC_DTYPE_TYPE_SHIFT)
#define GFC_DESCRIPTOR_SIZE(desc) ((desc)->dtype >> GFC_DTYPE_SIZE_SHIFT)
#define GFC_DTYPE_TYPE_SIZE(desc) ((desc)->dtype & GFC_DTYPE_TYPE_SIZE_MASK)

#endif

#define GFC_DTYPE_SIZE_MASK \
  ( ~((ptrdiff_t)(1 << GFC_DTYPE_SIZE_SHIFT) - 1)) // least significant bits to 0
#define GFC_DTYPE_TYPE_SIZE_MASK (GFC_DTYPE_SIZE_MASK | GFC_DTYPE_TYPE_MASK)

#define GFC_DTYPE_INTEGER_1 ((BT_INTEGER << GFC_DTYPE_TYPE_SHIFT) \
   | (sizeof(int8_t) << GFC_DTYPE_SIZE_SHIFT))
#define GFC_DTYPE_INTEGER_2 ((BT_INTEGER << GFC_DTYPE_TYPE_SHIFT) \
   | (sizeof(int16_t) << GFC_DTYPE_SIZE_SHIFT))
#define GFC_DTYPE_INTEGER_4 ((BT_INTEGER << GFC_DTYPE_TYPE_SHIFT) \
   | (sizeof(int32_t) << GFC_DTYPE_SIZE_SHIFT))
#define GFC_DTYPE_INTEGER_8 ((BT_INTEGER << GFC_DTYPE_TYPE_SHIFT) \
   | (sizeof(int64_t) << GFC_DTYPE_SIZE_SHIFT))
#if HAVE_INT128_T
#define GFC_DTYPE_INTEGER_16 ((BT_INTEGER << GFC_DTYPE_TYPE_SHIFT) \
   | (sizeof(__int128_t) << GFC_DTYPE_SIZE_SHIFT))
#endif

#define GFC_DTYPE_LOGICAL_4 ((BT_LOGICAL << GFC_DTYPE_TYPE_SHIFT) \
   | (sizeof(int) << GFC_DTYPE_SIZE_SHIFT))

#if 0
#define GFC_DTYPE_LOGICAL_1 ((BT_LOGICAL << GFC_DTYPE_TYPE_SHIFT) \
   | (sizeof(GFC_LOGICAL_1) << GFC_DTYPE_SIZE_SHIFT))
#define GFC_DTYPE_LOGICAL_2 ((BT_LOGICAL << GFC_DTYPE_TYPE_SHIFT) \
   | (sizeof(GFC_LOGICAL_2) << GFC_DTYPE_SIZE_SHIFT))
#define GFC_DTYPE_LOGICAL_8 ((BT_LOGICAL << GFC_DTYPE_TYPE_SHIFT) \
   | (sizeof(double) << GFC_DTYPE_SIZE_SHIFT))
#define GFC_DTYPE_LOGICAL_16 ((BT_LOGICAL << GFC_DTYPE_TYPE_SHIFT)\
   | (sizeof(GFC_LOGICAL_16) << GFC_DTYPE_SIZE_SHIFT))
#endif

#define GFC_DTYPE_REAL_4 ((BT_REAL << GFC_DTYPE_TYPE_SHIFT) \
   | (sizeof(float) << GFC_DTYPE_SIZE_SHIFT))
#define GFC_DTYPE_REAL_8 ((BT_REAL << GFC_DTYPE_TYPE_SHIFT) \
   | (sizeof(double) << GFC_DTYPE_SIZE_SHIFT))
#if 0
#ifdef HAVE_GFC_REAL_10
#define GFC_DTYPE_REAL_10  ((BT_REAL << GFC_DTYPE_TYPE_SHIFT) \
   | (sizeof(GFC_REAL_10) << GFC_DTYPE_SIZE_SHIFT))
#endif
#ifdef HAVE_GFC_REAL_16
#define GFC_DTYPE_REAL_16 ((BT_REAL << GFC_DTYPE_TYPE_SHIFT) \
   | (sizeof(GFC_REAL_16) << GFC_DTYPE_SIZE_SHIFT))
#endif
#endif

#define GFC_DTYPE_COMPLEX_4 ((BT_COMPLEX << GFC_DTYPE_TYPE_SHIFT) \
   | (sizeof(_Complex float) << GFC_DTYPE_SIZE_SHIFT))
#define GFC_DTYPE_COMPLEX_8 ((BT_COMPLEX << GFC_DTYPE_TYPE_SHIFT) \
   | (sizeof(_Complex double) << GFC_DTYPE_SIZE_SHIFT))
#if 0
#ifdef HAVE_GFC_COMPLEX_10
#define GFC_DTYPE_COMPLEX_10 ((BT_COMPLEX << GFC_DTYPE_TYPE_SHIFT) \
   | (sizeof(GFC_COMPLEX_10) << GFC_DTYPE_SIZE_SHIFT))
#endif
#ifdef HAVE_GFC_COMPLEX_16
#define GFC_DTYPE_COMPLEX_16 ((BT_COMPLEX << GFC_DTYPE_TYPE_SHIFT) \
   | (sizeof(GFC_COMPLEX_16) << GFC_DTYPE_SIZE_SHIFT))
#endif
#endif

/* FIXME: Hardwiring these values to what the mpi_caf.c macro GFC_DTYPE_TYPE_SIZE(desc)
    receives in the dtype component its gf_descriptor_t argument for character(kind=c_char)
    and logical(kind=c_bool) data:
*/

#ifdef GCC_GE_8

#define GFC_DTYPE_CHARACTER ((BT_CHARACTER << GFC_DTYPE_TYPE_SHIFT) \
   | (sizeof(char) << GFC_DTYPE_SIZE_SHIFT))

#else
#define GFC_DTYPE_CHARACTER 48
#endif


#endif  /* LIBCAF_GFORTRAN_DESCRIPTOR_H.  */
