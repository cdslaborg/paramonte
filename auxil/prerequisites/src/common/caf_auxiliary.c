/* Auxiliary functions for all of GNU Fortran libcaf implementations.

Copyright (c) 2012-2014, Sourcery, Inc.
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

#include "libcaf.h"


/* Check whether the array section is contiguous.  There are two possibilities
   either the stride matches always the extent of that dimension - or if the
   "noncontiguous" dimensions have all extent one (= element access for that
   dimension); a mixture is possible if the left dimensions are contiguous
   and the right ones are elements.  */

bool
PREFIX (is_contiguous) (gfc_descriptor_t *array)
{
  int i;
  ptrdiff_t dim_extent;
  ptrdiff_t extent = 1;
  bool element = false;

  for (i = 0; i < GFC_DESCRIPTOR_RANK(array); i++)
    {
      if (!element && array->dim[i]._stride != extent)
	return false;

      dim_extent = array->dim[i]._ubound - array->dim[i].lower_bound + 1;
      if (dim_extent <= 0)
	return true;  /* Zero-sized array.  */
      else if (dim_extent == 1 && GFC_DESCRIPTOR_RANK(array) == 1)
        element = true;
      else if (element)
	return false;
      extent *= dim_extent;
    }
  return true;
}
