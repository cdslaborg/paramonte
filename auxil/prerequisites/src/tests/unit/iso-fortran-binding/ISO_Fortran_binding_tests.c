/*
   OpenCoarrays is distributed under the OSI-approved BSD 3-clause License:
   OpenCoarrays -- ISO_Fortran_binding standard-compliant interoperability with
   C.
   Copyright (c) 2018, Sourcery, Inc.
   Copyright (c) 2018, Sourcery Institute
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
   2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
   3. Neither the names of the copyright holders nor the names of their
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
 */

#include "../../../iso-fortran-binding/ISO_Fortran_binding.h"
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

int main (void)
{

  CFI_rank_t      rank;
  CFI_attribute_t attribute;
  CFI_type_t      type[10] = {CFI_type_Bool,        CFI_type_short,
                         CFI_type_ucs4_char,   CFI_type_double,
                         CFI_type_float128,    CFI_type_float128_Complex,
                         CFI_type_long_double, CFI_type_long_double_Complex,
                         CFI_type_struct,      CFI_type_other};
  size_t elem_len;
  int    ind;
  size_t base_type;
  size_t base_type_size;
  size_t iso_errno;

  /* Test function establish. */
  /* Fresh descriptor, base address is NULL. */
  /* Loop through type. */
  for (int i = 0; i < 10; i++)
    {
      elem_len = 0;
      if (type[i] == CFI_type_char || type[i] == CFI_type_ucs4_char ||
          type[i] == CFI_type_signed_char || type[i] == CFI_type_struct ||
          type[i] == CFI_type_other)
        {
          base_type      = type[i];
          base_type_size = elem_len;
        }
      else
        {
          base_type      = type[i] & CFI_type_mask;
          base_type_size = (type[i] - base_type) >> CFI_type_kind_shift;
        }
      /* Loop through attribute. */
      for (int j = 1; j <= 3; j++)
        {
          attribute = j;
          /* Loop through rank. */
          for (int k = 0; k <= CFI_MAX_RANK; k++)
            {
              iso_errno = 1;
              rank  = k;
              CFI_CDESC_T (rank) test1;
              /* We do this because C sometimes doesn't make the structures with
               * a null base_addr which leads to weird behaviour inside
               * CFI_establish.
               */
              if (test1.base_addr != NULL)
                {
                  test1.base_addr = NULL;
                  free (test1.base_addr);
                }
              ind = CFI_establish ((CFI_cdesc_t *) &test1, NULL, attribute,
                                   type[i], elem_len, rank, NULL);
              if (ind != CFI_SUCCESS)
                {
                  goto next_attribute1;
                }
              if (attribute != test1.attribute)
                {
                  printf ("CFI_establish: failed to assign attribute.\n");
                  return 1;
                }
              if (type[i] != test1.type)
                {
                  printf ("CFI_establish: failed to assign type.\n");
                  return 1;
                }
              if (rank != test1.rank)
                {
                  printf ("CFI_establish: failed to assign rank.\n");
                  return 1;
                }
              elem_len = base_type_size;
              if (base_type_size == 10)
                {
                  elem_len = 64;
                }
              if (base_type == CFI_type_Complex)
                {
                  elem_len *= 2;
                }
              if (elem_len != test1.elem_len)
                {
                  printf ("CFI_establish: failed to assign element length.\n");
                  return 1;
                }
            }
        next_attribute1:;
        }
    }

  /* Fresh descriptor, base address is not NULL */
  CFI_index_t *extents = NULL;
  /* Loop through type. */
  for (int i = 0; i < 10; i++)
    {
      elem_len = 0;
      if (type[i] == CFI_type_char || type[i] == CFI_type_ucs4_char ||
          type[i] == CFI_type_signed_char || type[i] == CFI_type_struct ||
          type[i] == CFI_type_other)
        {
          base_type      = type[i];
          base_type_size = elem_len;
        }
      else
        {
          base_type      = type[i] & CFI_type_mask;
          base_type_size = (type[i] - base_type) >> CFI_type_kind_shift;
        }
      /* Loop through attribute. */
      for (int j = 1; j <= 3; j++)
        {
          attribute = j;
          /* Loop through rank. */
          for (int k = 0; k <= CFI_MAX_RANK; k++)
            {
              iso_errno = 1;
              rank  = k;
              if (extents != NULL)
                {
                  free (extents);
                }
              extents = malloc (rank * sizeof (CFI_index_t));
              for (int r = 0; r < rank; r++)
                {
                  extents[r] = r + 1;
                }
              CFI_CDESC_T (rank) test2;
              /* We do this because C sometimes doesn't make the structures with
               * a null base_addr which leads to weird behaviour inside
               * CFI_establish.
               */
              if (test2.base_addr != NULL)
                {
                  test2.base_addr = NULL;
                  free (test2.base_addr);
                }
              ind = CFI_establish ((CFI_cdesc_t *) &test2, &ind, attribute,
                                   type[i], elem_len, rank, extents);
              if (ind != CFI_SUCCESS)
                {
                  goto next_attribute2;
                }
              if (attribute != test2.attribute)
                {
                  printf ("CFI_establish: failed to assign attribute.\n");
                  return 1;
                }
              if (type[i] != test2.type)
                {
                  printf ("CFI_establish: failed to assign type.\n");
                  return 1;
                }
              if (rank != test2.rank)
                {
                  printf ("CFI_establish: failed to assign rank.\n");
                  return 1;
                }

              elem_len = base_type_size;
              if (base_type_size == 10)
                {
                  elem_len = 64;
                }
              if (base_type == CFI_type_Complex)
                {
                  elem_len *= 2;
                }
              if (elem_len != test2.elem_len)
                {
                  printf ("CFI_establish: failed to assign element length.\n");
                  return 1;
                }

              for (int r = 0; r < rank; r++)
                {
                  if (extents[r] != test2.dim[r].extent)
                    {
                      printf ("CFI_establish: failed to assign dimension "
                              "extents.\n");
                      return 1;
                    }
                }

              if (attribute == CFI_attribute_pointer)
                {
                  for (int r = 0; r < rank; r++)
                    {
                      if (test2.dim[r].lower_bound != 0)
                        {
                          printf ("CFI_establish: failed to assign dimension "
                                  "lower bounds.\n");
                          return 1;
                        }
                    }
                }
            }
        next_attribute2:;
        }
    }

  /* Fresh descriptor, base address is not NULL */
  CFI_index_t *lower = NULL;
  CFI_index_t *upper = NULL;
  /* Loop through type. */
  for (int i = 0; i < 10; i++)
    {
      elem_len = 0;
      if (type[i] == CFI_type_struct)
        {
          base_type      = type[i];
          base_type_size = 69;
        }
      else if (type[i] == CFI_type_other)
        {
          base_type      = type[i];
          base_type_size = 666;
        }
      else if (type[i] == CFI_type_char || type[i] == CFI_type_ucs4_char ||
               type[i] == CFI_type_signed_char)
        {
          base_type      = type[i] & CFI_type_mask;
          base_type_size = 3;
        }
      else
        {
          base_type      = type[i] & CFI_type_mask;
          base_type_size = (type[i] - base_type) >> CFI_type_kind_shift;
        }

      elem_len = base_type_size;
      if (base_type_size == 10)
        {
          elem_len = 64;
        }
      if (base_type == CFI_type_Complex)
        {
          elem_len *= 2;
        }
      /* Loop through attribute. */
      for (int j = 1; j <= 3; j++)
        {
          attribute = j;
          /* Loop through rank. */
          for (int k = 0; k <= CFI_MAX_RANK; k++)
            {
              iso_errno = 1;
              rank  = k;
              if (extents != NULL)
                {
                  free (extents);
                }
              if (lower != NULL)
                {
                  free (lower);
                }
              if (upper != NULL)
                {
                  free (upper);
                }
              extents = malloc (rank * sizeof (CFI_index_t));
              lower   = malloc (rank * sizeof (CFI_index_t));
              upper   = malloc (rank * sizeof (CFI_index_t));
              for (int r = 0; r < rank; r++)
                {
                  extents[r] = 2;
                  lower[r]   = r;
                  upper[r]   = lower[r] + extents[r];
                }
              CFI_CDESC_T (rank) test3;
              /* We do this because C sometimes doesn't make the structures with
               * a null base_addr which leads to weird behaviour inside
               * CFI_establish.
               */
              if (test3.base_addr != NULL)
                {
                  test3.base_addr = NULL;
                  free (test3.base_addr);
                }
              ind = CFI_establish ((CFI_cdesc_t *) &test3, NULL, attribute,
                                   type[i], elem_len, rank, extents);
              ind =
                  CFI_allocate ((CFI_cdesc_t *) &test3, lower, upper, elem_len);
              if (ind != CFI_SUCCESS)
                {
                  goto next_attribute3;
                }
              for (int r = 0; r < rank; r++)
                {
                  if (lower[r] != test3.dim[r].lower_bound)
                    {
                      printf ("CFI_allocate: failed to reassign dimension "
                              "lower bounds.\n");
                      return 1;
                    }
                  if (upper[r] - test3.dim[r].lower_bound + 1 !=
                      test3.dim[r].extent)
                    {
                      printf ("CFI_allocate: failed to reassign dimension "
                              "extents.\n");
                      return 1;
                    }
                  if (test3.dim[r].sm != test3.elem_len)
                    {
                      printf (
                          "CFI_allocate: failed to assign dimension stride.\n");
                      return 1;
                    }
                }
              if (elem_len != test3.elem_len)
                {
                  printf ("CFI_allocate: failed to reassign element length.\n");
                  return 1;
                }
            }
        next_attribute3:;
        }
    }

  rank  = 1;
  iso_errno = 1;
  CFI_CDESC_T (rank) test4;
  base_type      = type[3] & CFI_type_mask;
  base_type_size = (type[3] - base_type) >> CFI_type_kind_shift;
  attribute      = CFI_attribute_allocatable;
  ind = CFI_establish ((CFI_cdesc_t *) &test4, NULL, attribute, type[3],
                       elem_len, rank, NULL);
  ind = CFI_allocate ((CFI_cdesc_t *) &test4, NULL, NULL, base_type_size);
  if (ind != CFI_INVALID_EXTENT)
    {
      printf ("CFI_allocate: failed to detect invalid extents.\n");
      return 1;
    }

  rank  = 1;
  iso_errno = 1;
  CFI_CDESC_T (rank) test5;
  base_type      = type[3] & CFI_type_mask;
  base_type_size = (type[3] - base_type) >> CFI_type_kind_shift;
  attribute      = CFI_attribute_pointer;
  ind = CFI_establish ((CFI_cdesc_t *) &test5, &ind, attribute, type[3],
                       elem_len, rank, extents);
  ind = CFI_allocate ((CFI_cdesc_t *) &test5, NULL, NULL, base_type_size);
  if (ind != CFI_ERROR_BASE_ADDR_NOT_NULL)
    {
      printf ("CFI_allocate: failed to detect base address is not NULL.\n");
      return 1;
    }

  /* Test CFI_deallocate. */
  rank           = 1;
  iso_errno          = 1;
  base_type      = type[3] & CFI_type_mask;
  base_type_size = (type[3] - base_type) >> CFI_type_kind_shift;
  for (int i = 1; i <= 3; i++)
    {
      attribute = i;
      if (extents != NULL)
        {
          free (extents);
        }
      if (lower != NULL)
        {
          free (lower);
        }
      if (upper != NULL)
        {
          free (upper);
        }
      extents = malloc (rank * sizeof (CFI_index_t));
      lower   = malloc (rank * sizeof (CFI_index_t));
      upper   = malloc (rank * sizeof (CFI_index_t));
      CFI_CDESC_T (rank) test6;
      ind = CFI_establish ((CFI_cdesc_t *) &test6, NULL, attribute, type[i],
                           elem_len, rank, extents);
      ind = CFI_allocate ((CFI_cdesc_t *) &test6, lower, upper, base_type_size);
      if (ind == CFI_SUCCESS)
        {
          ind = CFI_deallocate ((CFI_cdesc_t *) &test6);
          if (ind != CFI_INVALID_ATTRIBUTE && test6.base_addr != NULL)
            {
              printf ("CFI_deallocate: failed to deallocate memory.\n");
              return 1;
            }
        }
    }

  /* Test CFI_is_contiguous. */
  int tmp_ind;
  base_type      = type[3] & CFI_type_mask;
  base_type_size = (type[3] - base_type) >> CFI_type_kind_shift;
  for (int i = 1; i <= 3; i++)
    {
      attribute = i;
      for (int j = 0; j <= 4; j++)
        {
          iso_errno = 1;
          rank  = j;
          if (extents != NULL)
            {
              free (extents);
            }
          if (lower != NULL)
            {
              free (lower);
            }
          if (upper != NULL)
            {
              free (upper);
            }
          extents = malloc (rank * sizeof (CFI_index_t));
          lower   = malloc (rank * sizeof (CFI_index_t));
          upper   = malloc (rank * sizeof (CFI_index_t));
          for (int r = 0; r < rank; r++)
            {
              extents[r] = 2;
              lower[r]   = r;
              upper[r]   = lower[r] + extents[r];
            }
          CFI_CDESC_T (rank) test7;
          ind = CFI_establish ((CFI_cdesc_t *) &test7, NULL, attribute, type[3],
                               elem_len, rank, extents);
          tmp_ind = CFI_allocate ((CFI_cdesc_t *) &test7, lower, upper,
                                  base_type_size);
          if (tmp_ind != CFI_SUCCESS)
            {
              goto next_attribute4;
            }
          ind = CFI_is_contiguous ((CFI_cdesc_t *) &test7);
          if (ind != CFI_INVALID_RANK && rank == 0 &&
              tmp_ind != CFI_INVALID_ATTRIBUTE)
            {
              printf ("CFI_is_contiguous: failed to detect incorrect rank.\n");
              return 1;
            }
          else if (ind == CFI_ERROR_BASE_ADDR_NULL && test7.base_addr != NULL &&
                   tmp_ind != CFI_SUCCESS)
            {
              printf ("CFI_is_contiguous: failed to detect base address is not "
                      "NULL.\n");
              return 1;
            }
        }
      next_attribute4:;
    }

  /* Test CFI_address. */
  CFI_index_t *tr_subscripts;
  CFI_dim_t *  tr_dim;
  /* Loop through type. */
  for (int i = 0; i < 10; i++)
    {
      elem_len = 0;
      if (type[i] == CFI_type_struct)
        {
          base_type      = type[i];
          base_type_size = 69;
        }
      else if (type[i] == CFI_type_other)
        {
          base_type      = type[i];
          base_type_size = 666;
        }
      else if (type[i] == CFI_type_char || type[i] == CFI_type_ucs4_char ||
               type[i] == CFI_type_signed_char)
        {
          base_type      = type[i] & CFI_type_mask;
          base_type_size = 3;
        }
      else
        {
          base_type      = type[i] & CFI_type_mask;
          base_type_size = (type[i] - base_type) >> CFI_type_kind_shift;
        }

      elem_len = base_type_size;
      if (base_type_size == 10)
        {
          elem_len = 64;
        }
      if (base_type == CFI_type_Complex)
        {
          elem_len *= 2;
        }
      /* Loop through attribute. */
      for (int j = 1; j <= 3; j++)
        {
          attribute = j;
          /* Loop through rank. */
          for (int k = 1; k <= CFI_MAX_RANK; k++)
            {
              iso_errno = 1;
              rank  = k;
              CFI_CDESC_T (rank) source;
              if (extents != NULL)
                {
                  free (extents);
                }
              if (lower != NULL)
                {
                  free (lower);
                }
              if (upper != NULL)
                {
                  free (upper);
                }
              extents = malloc (rank * sizeof (CFI_index_t));
              lower   = malloc (rank * sizeof (CFI_index_t));
              upper   = malloc (rank * sizeof (CFI_index_t));
              for (int r = 0; r < rank; r++)
                {
                  extents[r] = rank - r + 1;
                  lower[r]   = rank - r - 3;
                  upper[r]   = lower[r] + extents[r] - 1;
                }
              ind = CFI_establish ((CFI_cdesc_t *) &source, NULL,
                                   CFI_attribute_allocatable, type[i], elem_len,
                                   rank, extents);
              ind = CFI_allocate ((CFI_cdesc_t *) &source, lower, upper,
                                  elem_len);
              if (ind == CFI_SUCCESS)
                {
                  CFI_index_t dif_addr;
                  CFI_index_t n_entries = 1;
                  dif_addr              = (CFI_index_t) (
                      (char *) CFI_address ((CFI_cdesc_t *) &source, upper) -
                      (char *) CFI_address ((CFI_cdesc_t *) &source, lower));
                  for (int r = 0; r < rank; r++)
                    {
                      n_entries = n_entries * (upper[r] - lower[r] + 1);
                    }
                  tr_subscripts = malloc (rank * sizeof (CFI_index_t));
                  tr_dim        = malloc (rank * sizeof (CFI_dim_t));
                  for (int i = 0; i < rank; i++)
                    {
                      CFI_index_t idx  = rank - i - 1;
                      tr_subscripts[i] = upper[idx];
                      tr_dim[i]        = source.dim[idx];
                      /* Normalise the subscripts to start counting the address
                       * from 0. */
                      tr_subscripts[i] -= tr_dim[i].lower_bound;
                    }
                  /* We assume column major order as that is how Fortran stores
                   * arrays. We
                   * calculate the memory address of the specified element via
                   * the canonical
                   * array dimension reduction map and multiplying by the memory
                   * stride. */
                  CFI_index_t index     = tr_subscripts[0] * tr_dim[0].sm;
                  CFI_index_t tmp_index = 1;
                  for (int i = 1; i < rank; i++)
                    {
                      tmp_index *= tr_subscripts[i] * tr_dim[i - 1].extent *
                                   tr_dim[i - 1].sm;
                      index += tmp_index;
                    }
                  free (tr_subscripts);
                  free (tr_dim);
                  if (index - dif_addr != 0)
                    {
                      printf ("CFI_address: difference in address is not being "
                              "properly calculated.\n");
                      return 1;
                    }
                }
              else if (ind == CFI_ERROR_MEM_ALLOCATION)
                {
                  goto next_type;
                }
            }
        }
    next_type:;
    }

  /* Test CFI_setpointer */
  for (int i = 0; i < CFI_MAX_RANK; i++)
    {
      rank           = i;
      iso_errno          = 1;
      base_type      = type[3] & CFI_type_mask;
      base_type_size = (type[3] - base_type) >> CFI_type_kind_shift;
      attribute      = CFI_attribute_other;
      CFI_CDESC_T (rank) test8a, test8b;

      if (extents != NULL)
        {
          free (extents);
        }
      if (lower != NULL)
        {
          free (lower);
        }
      extents = malloc (rank * sizeof (CFI_index_t));
      lower   = malloc (rank * sizeof (CFI_index_t));
      for (int r = 0; r < rank; r++)
        {
          extents[r] = r + 1;
          lower[r]   = r - 2;
        }
      ind = CFI_establish ((CFI_cdesc_t *) &test8a, &ind, attribute, type[3],
                           base_type_size, rank, extents);
      for (int r = 0; r < rank; r++)
        {
          extents[r] = r + 2;
        }
      ind = CFI_establish ((CFI_cdesc_t *) &test8b, &iso_errno, attribute, type[3],
                           base_type_size, rank, extents);
      ind = CFI_setpointer ((CFI_cdesc_t *) &test8a, (CFI_cdesc_t *) &test8b,
                            lower);
      for (int r = 0; r < rank; r++)
        {
          if (test8a.dim[r].lower_bound != lower[r])
            {
              printf ("CFI_setpointer: failed to reassign lower bounds.\n");
              return 1;
            }
          if (test8a.dim[r].extent != test8b.dim[r].extent)
            {
              printf ("CFI_setpointer: failed to reassign extents.\n");
              return 1;
            }
          if (test8a.dim[r].sm != test8b.dim[r].sm)
            {
              printf ("CFI_setpointer: failed to reassign memory strides.\n");
              return 1;
            }
        }
      if (test8a.base_addr != test8b.base_addr)
        {
          printf ("CFI_setpointer: failed to reassign base address.\n");
          return 1;
        }
      if (test8a.version != test8b.version)
        {
          printf ("CFI_setpointer: failed to reassign lower bounds.\n");
          return 1;
        }
      if (test8a.attribute != test8b.attribute)
        {
          printf ("CFI_setpointer: failed to reassign attribute.\n");
          return 1;
        }
    }

  /* NULL source. */
  rank           = 10;
  iso_errno          = 1;
  base_type      = type[3] & CFI_type_mask;
  base_type_size = (type[3] - base_type) >> CFI_type_kind_shift;
  CFI_CDESC_T (rank) test9;

  if (extents != NULL)
    {
      free (extents);
    }
  if (lower != NULL)
    {
      free (lower);
    }
  extents = malloc (rank * sizeof (CFI_index_t));
  lower   = malloc (rank * sizeof (CFI_index_t));
  for (int r = 0; r < rank; r++)
    {
      extents[r] = r + 1;
      lower[r]   = r - 2;
    }
  ind = CFI_establish ((CFI_cdesc_t *) &test9, &ind, attribute, type[3],
                       base_type_size, rank, extents);
  ind = CFI_setpointer ((CFI_cdesc_t *) &test9, NULL, lower);
  if (test9.attribute != CFI_attribute_pointer)
    {
      printf ("CFI_setpointer: failed to set attribute pointer.\n");
      return 1;
    }
  if (test9.base_addr != NULL)
    {
      printf ("CFI_setpointer: failed to set base address to NULL.\n");
      return 1;
    }

  rank      = 3;
  iso_errno     = 1;
  attribute = CFI_attribute_other;
  CFI_CDESC_T (rank) test10a, test10b;
  if (extents != NULL)
    {
      free (extents);
    }
  if (lower != NULL)
    {
      free (lower);
    }
  extents = malloc (rank * sizeof (CFI_index_t));
  lower   = malloc (rank * sizeof (CFI_index_t));
  for (int r = 0; r < rank; r++)
    {
      extents[r] = r + 1;
      lower[r]   = r - 2;
    }
  base_type      = CFI_type_long & CFI_type_mask;
  base_type_size = (CFI_type_long - base_type) >> CFI_type_kind_shift;
  ind = CFI_establish ((CFI_cdesc_t *) &test10a, &ind, attribute, CFI_type_long,
                       base_type_size, rank, extents);
  for (int r = 0; r < rank; r++)
    {
      extents[r] = r + 2;
    }
  base_type      = CFI_type_double & CFI_type_mask;
  base_type_size = (CFI_type_double - base_type) >> CFI_type_kind_shift;
  ind            = CFI_establish ((CFI_cdesc_t *) &test10b, &iso_errno, attribute,
                       CFI_type_double, base_type_size, rank, extents);
  ind = CFI_setpointer ((CFI_cdesc_t *) &test10a, (CFI_cdesc_t *) &test10b,
                        lower);
  if (ind != CFI_INVALID_TYPE)
    {
      printf ("CFI_setpointer: failed to detect invalid type.\n");
      return 1;
    }

  iso_errno          = 1;
  base_type      = CFI_type_other & CFI_type_mask;
  base_type_size = 666;
  ind            = CFI_establish ((CFI_cdesc_t *) &test10a, &ind, attribute,
                       CFI_type_other, base_type_size, rank, extents);
  base_type      = CFI_type_other & CFI_type_mask;
  base_type_size = 69;
  ind            = CFI_establish ((CFI_cdesc_t *) &test10b, &iso_errno, attribute,
                       CFI_type_other, base_type_size, rank, extents);
  ind = CFI_setpointer ((CFI_cdesc_t *) &test10a, (CFI_cdesc_t *) &test10b,
                        lower);
  if (ind != CFI_INVALID_ELEM_LEN)
    {
      printf ("CFI_setpointer: failed to detect invalid element length.\n");
      return 1;
    }

  iso_errno          = 1;
  base_type      = type[3] & CFI_type_mask;
  base_type_size = (CFI_type_long - base_type) >> CFI_type_kind_shift;
  ind = CFI_establish ((CFI_cdesc_t *) &test10a, &ind, attribute, type[3],
                       base_type_size, rank, extents);
  rank++;
  CFI_CDESC_T (rank) test10c;
  if (extents != NULL)
    {
      free (extents);
    }
  if (lower != NULL)
    {
      free (lower);
    }
  extents = malloc (rank * sizeof (CFI_index_t));
  lower   = malloc (rank * sizeof (CFI_index_t));
  for (int r = 0; r < rank; r++)
    {
      extents[r] = r + 1;
      lower[r]   = r - 2;
    }
  base_type      = CFI_type_other & CFI_type_mask;
  base_type_size = (CFI_type_long - base_type) >> CFI_type_kind_shift;
  ind = CFI_establish ((CFI_cdesc_t *) &test10c, &iso_errno, attribute, type[3],
                       base_type_size, rank, extents);
  ind = CFI_setpointer ((CFI_cdesc_t *) &test10a, (CFI_cdesc_t *) &test10c,
                        lower);
  if (ind != CFI_INVALID_RANK)
    {
      printf ("CFI_setpointer: failed to detect invalid rank.\n");
      return 1;
    }

  /* Test CFI_section */
  CFI_index_t *strides;
  /* Loop through type. */
  for (int i = 0; i < 10; i++)
    {
      elem_len = 0;
      if (type[i] == CFI_type_struct)
        {
          base_type      = type[i];
          base_type_size = 69;
        }
      else if (type[i] == CFI_type_other)
        {
          base_type      = type[i];
          base_type_size = 666;
        }
      else if (type[i] == CFI_type_char || type[i] == CFI_type_ucs4_char ||
               type[i] == CFI_type_signed_char)
        {
          base_type      = type[i] & CFI_type_mask;
          base_type_size = 3;
        }
      else
        {
          base_type      = type[i] & CFI_type_mask;
          base_type_size = (type[i] - base_type) >> CFI_type_kind_shift;
        }
      elem_len = base_type_size;
      if (base_type_size == 10)
        {
          elem_len = 64;
        }
      if (base_type == CFI_type_Complex)
        {
          elem_len *= 2;
        }
      /* Loop through rank. */
      for (int k = 1; k <= CFI_MAX_RANK; k++)
        {
          iso_errno = 1;
          rank  = k;
          CFI_CDESC_T (rank) section, source;
          if (extents != NULL)
            {
              free (extents);
            }
          if (lower != NULL)
            {
              free (lower);
            }
          if (upper != NULL)
            {
              free (upper);
            }
          if (strides == NULL)
            {
              free (strides);
            }
          extents = malloc (rank * sizeof (CFI_index_t));
          lower   = malloc (rank * sizeof (CFI_index_t));
          upper   = malloc (rank * sizeof (CFI_index_t));
          strides = malloc (rank * sizeof (CFI_index_t));
          for (int r = 0; r < rank; r++)
            {
              extents[r] = rank - r + 10;
              lower[r]   = rank - r - 5;
              upper[r]   = lower[r] + extents[r] - 1;
            }
          ind = CFI_establish ((CFI_cdesc_t *) &source, NULL,
                               CFI_attribute_allocatable, type[i], elem_len,
                               rank, extents);
          ind = CFI_establish ((CFI_cdesc_t *) &section, NULL,
                               CFI_attribute_other, type[i], elem_len, rank,
                               NULL);
          ind = CFI_allocate ((CFI_cdesc_t *) &source, lower, upper, elem_len);
          if (ind != CFI_SUCCESS)
            {
              goto next_type2;
            }
          /* Lower is within bounds. */
          for (int r = 0; r < rank; r++)
            {
              lower[r]   = rank - r - 3;
              strides[r] = r + 1;
            }
          ind = CFI_section ((CFI_cdesc_t *) &section, (CFI_cdesc_t *) &source,
                             lower, NULL, strides);
          if (ind != CFI_SUCCESS)
            {
              printf ("CFI_section: failed to detect lower bounds are within "
                      "bounds.\n");
              return 1;
            }
          /* Lower is below lower bounds. */
          for (int r = 0; r < rank; r++)
            {
              lower[r]   = rank - r - 6;
              strides[r] = r + 1;
            }
          ind = CFI_section ((CFI_cdesc_t *) &section, (CFI_cdesc_t *) &source,
                             lower, NULL, strides);
          if (ind != CFI_ERROR_OUT_OF_BOUNDS)
            {
              printf ("CFI_section: failed to detect lower bounds are below "
                      "bounds.\n");
              return 1;
            }
          /* Lower is above upper bounds. */
          for (int r = 0; r < rank; r++)
            {
              lower[r]   = upper[r] + 1;
              strides[r] = r + 1;
            }
          ind = CFI_section ((CFI_cdesc_t *) &section, (CFI_cdesc_t *) &source,
                             lower, NULL, strides);
          if (ind != CFI_ERROR_OUT_OF_BOUNDS)
            {
              printf ("CFI_section: failed to detect lower bounds are above "
                      "bounds.\n");
              return 1;
            }
          for (int r = 0; r < rank; r++)
            {
              extents[r] = rank - r + 10;
              lower[r]   = rank - r - 5;
              upper[r]   = lower[r] + extents[r] - 1;
            }
          /* Upper is within bounds. */
          for (int r = 0; r < rank; r++)
            {
              upper[r]   = rank - r - 3;
              strides[r] = r + 1;
            }
          ind = CFI_section ((CFI_cdesc_t *) &section, (CFI_cdesc_t *) &source,
                             NULL, upper, strides);
          if (ind != CFI_SUCCESS)
            {
              printf ("CFI_section: failed to detect upper bounds are within "
                      "bounds.\n");
              return 1;
            }
          /* Upper is below lower bounds. */
          for (int r = 0; r < rank; r++)
            {
              upper[r]   = rank - r - 6;
              strides[r] = r + 1;
            }
          ind = CFI_section ((CFI_cdesc_t *) &section, (CFI_cdesc_t *) &source,
                             NULL, upper, strides);
          if (ind != CFI_ERROR_OUT_OF_BOUNDS)
            {
              printf ("CFI_section: failed to detect upper bounds are below "
                      "bounds.\n");
              return 1;
            }
          /* Upper is above upper bounds. */
          for (int r = 0; r < rank; r++)
            {
              upper[r]   = lower[r] + extents[r];
              strides[r] = r + 1;
            }
          ind = CFI_section ((CFI_cdesc_t *) &section, (CFI_cdesc_t *) &source,
                             NULL, upper, strides);
          if (ind != CFI_ERROR_OUT_OF_BOUNDS)
            {
              printf ("CFI_section: failed to detect lower bounds are above "
                      "bounds.\n");
              return 1;
            }
          for (int r = 0; r < rank; r++)
            {
              extents[r] = rank - r + 10;
              lower[r]   = rank - r - 3;
              upper[r]   = lower[r] + extents[r] - 3;
              strides[r] = r + 1;
            }
          ind = CFI_section ((CFI_cdesc_t *) &section, (CFI_cdesc_t *) &source,
                             lower, upper, strides);
          for (int i = 0; i < rank; i++)
            {
              if (section.dim[i].lower_bound != lower[i])
                {
                  printf ("CFI_section: failed to correctly assign lower "
                          "bounds.\n");
                  return 1;
                }
              if (section.dim[i].extent != upper[i] - lower[i] + 1)
                {
                  printf ("CFI_section: failed to correctly assign extents.\n");
                  return 1;
                }
              if (section.dim[i].sm != strides[i] * section.elem_len)
                {
                  printf ("CFI_section: failed to correctly assign memory "
                          "strides.\n");
                  return 1;
                }
            }
        }
    next_type2:;
    }

  iso_errno = 1;
  rank  = 1;
  CFI_CDESC_T (rank) section, source;
  if (extents != NULL)
    {
      free (extents);
    }
  if (lower != NULL)
    {
      free (lower);
    }
  if (upper != NULL)
    {
      free (upper);
    }
  if (strides != NULL)
    {
      free (strides);
    }
  extents = malloc (rank * sizeof (CFI_index_t));
  lower   = malloc (rank * sizeof (CFI_index_t));
  upper   = malloc (rank * sizeof (CFI_index_t));
  strides = malloc (rank * sizeof (CFI_index_t));
  for (int r = 0; r < rank; r++)
    {
      extents[r] = rank - r + 10;
      lower[r]   = rank - r - 5;
      upper[r]   = lower[r] + extents[r] - 1;
    }
  ind = CFI_establish ((CFI_cdesc_t *) &source, NULL, CFI_attribute_allocatable,
                       type[3], elem_len, rank, extents);
  ind = CFI_establish ((CFI_cdesc_t *) &section, NULL, CFI_attribute_other,
                       type[3], elem_len, rank, NULL);
  ind = CFI_allocate ((CFI_cdesc_t *) &source, lower, upper, elem_len);
  if (ind == CFI_SUCCESS)
    {
      for (int r = 0; r < rank; r++)
        {
          lower[r]   = rank - r - 3;
          strides[r] = r + 1;
          upper[r]   = lower[r] + extents[r] - 3;
        }
      ind = CFI_section ((CFI_cdesc_t *) &section, NULL, lower, upper, strides);
      if (ind != CFI_INVALID_DESCRIPTOR)
        {
          printf ("CFI_section: failed to detect that source is NULL.\n");
          return 1;
        }
      ind = CFI_section (NULL, (CFI_cdesc_t *) &source, lower, upper, strides);
      if (ind != CFI_INVALID_DESCRIPTOR)
        {
          printf ("CFI_section: failed to detect that section is NULL.\n");
          return 1;
        }
      ind =
          CFI_establish ((CFI_cdesc_t *) &section, NULL, CFI_attribute_allocatable,
                         type[3], elem_len, rank, NULL);
      ind = CFI_section ((CFI_cdesc_t *) &section, (CFI_cdesc_t *) &source, lower,
                         upper, strides);
      if (ind != CFI_INVALID_ATTRIBUTE)
        {
          printf ("CFI_section: failed to detect invalid attribute.\n");
          return 1;
        }
      ind = CFI_establish ((CFI_cdesc_t *) &section, NULL, CFI_attribute_other,
                           type[3], elem_len, rank, NULL);
      ind = CFI_deallocate ((CFI_cdesc_t *) &source);
      ind = CFI_section ((CFI_cdesc_t *) &section, (CFI_cdesc_t *) &source, lower,
                         upper, strides);
      if (ind != CFI_ERROR_BASE_ADDR_NULL)
        {
          printf ("CFI_section: failed to detect that the base address is NULL.\n");
          return 1;
        }
    }

  CFI_CDESC_T (0) section2, source2;
  ind = CFI_establish ((CFI_cdesc_t *) &source2, &ind, CFI_attribute_other,
                       type[3], 0, 0, NULL);
  ind = CFI_establish ((CFI_cdesc_t *) &section2, &iso_errno, CFI_attribute_other,
                       type[3], 0, 0, NULL);
  ind = CFI_section ((CFI_cdesc_t *) &section2, (CFI_cdesc_t *) &source2, lower,
                     upper, strides);
  if (ind != CFI_INVALID_RANK)
    {
      printf ("CFI_section: failed to detect invalid rank.\n");
      return 1;
    }

  ind = CFI_establish ((CFI_cdesc_t *) &source, NULL, CFI_attribute_allocatable,
                       type[3], 0, rank, extents);
  ind = CFI_establish ((CFI_cdesc_t *) &section, NULL, CFI_attribute_other,
                       type[6], 0, rank, NULL);
  ind = CFI_allocate ((CFI_cdesc_t *) &source, lower, upper, elem_len);
  if (ind == CFI_SUCCESS)
    {
      for (int r = 0; r < rank; r++)
        {
          lower[r]   = rank - r - 3;
          strides[r] = r + 1;
          upper[r]   = lower[r] + extents[r] - 3;
        }
      ind = CFI_section ((CFI_cdesc_t *) &section, (CFI_cdesc_t *) &source, lower,
                         upper, strides);
      if (ind != CFI_INVALID_ELEM_LEN)
        {
          printf ("CFI_section: failed to detect incompatible element lengths "
                  "between source and section.\n");
          return 1;
        }
      ind = CFI_establish ((CFI_cdesc_t *) &section, NULL, CFI_attribute_other,
                           CFI_type_long, 0, rank, NULL);
      ind = CFI_section ((CFI_cdesc_t *) &section, (CFI_cdesc_t *) &source, lower,
                         upper, strides);
      if (ind != CFI_INVALID_TYPE)
        {
          printf ("CFI_section: failed to detect invalid type.\n");
          return 1;
        }
    }

  for (int i = 1; i < CFI_MAX_RANK; i++)
    {
      iso_errno   = 1;
      rank    = i;
      int ctr = 0;
      CFI_CDESC_T (rank) source;
      if (extents != NULL)
        {
          free (extents);
        }
      if (lower != NULL)
        {
          free (lower);
        }
      if (upper != NULL)
        {
          free (upper);
        }
      if (strides != NULL)
        {
          free (strides);
        }
      extents = malloc (rank * sizeof (CFI_index_t));
      lower   = malloc (rank * sizeof (CFI_index_t));
      upper   = malloc (rank * sizeof (CFI_index_t));
      strides = malloc (rank * sizeof (CFI_index_t));
      for (int r = 0; r < rank; r++)
        {
          extents[r] = rank - r + 10;
          lower[r]   = rank - r - 5;
          upper[r]   = lower[r] + extents[r] - 1;
        }
      ind =
          CFI_establish ((CFI_cdesc_t *) &source, NULL,
                         CFI_attribute_allocatable, type[3], 0, rank, extents);
      ind = CFI_allocate ((CFI_cdesc_t *) &source, lower, upper, elem_len);
      if (ind != CFI_SUCCESS)
        {
          continue;
        }
      for (int r = 0; r < rank; r++)
        {
          lower[r] = rank - r - 3;
          if (r % 2 == 0)
            {
              strides[r] = 0;
              upper[r]   = lower[r];
              ctr++;
            }
          else
            {
              strides[r] = r + 1;
              upper[r]   = lower[r] + extents[r] - 3;
            }
        }
      CFI_CDESC_T (rank - ctr) section;
      ind = CFI_establish ((CFI_cdesc_t *) &section, NULL, CFI_attribute_other,
                           type[3], 0, rank - ctr, NULL);
      ind = CFI_section ((CFI_cdesc_t *) &section, (CFI_cdesc_t *) &source,
                         lower, upper, strides);
      ctr = 0;
      for (int r = 0; r < rank; r++)
        {
          if (strides[r] == 0)
            {
              ctr++;
              continue;
            }
          int idx = r - ctr;
          if (section.dim[idx].lower_bound != lower[r])
            {
              printf ("CFI_section: failed to correctly assign lower bounds in "
                      "rank reduction.\n");
              return 1;
            }
          if (section.dim[idx].extent != upper[r] - lower[r] + 1)
            {
              printf ("CFI_section: failed to correctly assign extents in rank "
                      "reduction.\n");
              return 1;
            }
          if (section.dim[idx].sm != strides[r] * section.elem_len)
            {
              printf ("CFI_section: failed to correctly assign memory strides "
                      "in rank reduction.\n");
              return 1;
            }
          CFI_CDESC_T (rank - ctr - 1) section2;
          ind = CFI_establish ((CFI_cdesc_t *) &section2, NULL,
                               CFI_attribute_other, type[3], 0, rank - ctr - 1,
                               NULL);
          ind = CFI_section ((CFI_cdesc_t *) &section2, (CFI_cdesc_t *) &source,
                             lower, upper, strides);
          if (ind != CFI_SUCCESS && ind != CFI_INVALID_RANK)
            {
              printf ("CFI_section: failed to detect invalid rank.\n");
              return 1;
            }
        }
    }

  /* CFI_section negative strides. */
  iso_errno = 1;
  rank  = 8;
  if (extents != NULL)
    {
      free (extents);
    }
  if (lower != NULL)
    {
      free (lower);
    }
  if (upper != NULL)
    {
      free (upper);
    }
  if (strides != NULL)
    {
      free (strides);
    }
  extents = malloc (rank * sizeof (CFI_index_t));
  lower   = malloc (rank * sizeof (CFI_index_t));
  upper   = malloc (rank * sizeof (CFI_index_t));
  strides = malloc (rank * sizeof (CFI_index_t));
  for (int r = 0; r < rank; r++)
    {
      extents[r] = rank - r + 10;
      lower[r]   = rank - r - 3;
      upper[r]   = lower[r] + extents[r] - 3;
      strides[r] = -(r + 1);
    }
  CFI_CDESC_T (rank) section3, source3;
  ind = CFI_establish ((CFI_cdesc_t *) &source3, NULL,
                       CFI_attribute_allocatable, type[3], 0, rank, extents);
  ind = CFI_establish ((CFI_cdesc_t *) &section3, NULL, CFI_attribute_other,
                       type[3], 0, rank, NULL);
  ind = CFI_allocate ((CFI_cdesc_t *) &source3, lower, upper, elem_len);
  if (ind == CFI_SUCCESS)
    {
      ind = CFI_section ((CFI_cdesc_t *) &section3, (CFI_cdesc_t *) &source3,
                         upper, lower, strides);
      if (ind != CFI_SUCCESS && ind != CFI_INVALID_STRIDE)
        {
          printf ("CFI_section: failed to detect invalid stride.\n");
          return 1;
        }
    }

  /* CFI_select_part */
  typedef struct foo_t
  {
    int w;
    double _Complex p;
    double _Complex y;
    double z;
    double x;
  } foo_t;
  rank = 2;
  CFI_CDESC_T (rank) foo_c, cx, cy;
  int         arr_len = 100;
  foo_t       foo[arr_len][arr_len];
  CFI_index_t extent[] = {arr_len, arr_len};
  /* Establish c descriptor for the structure. */
  ind = CFI_establish ((CFI_cdesc_t *) &foo_c, &foo, CFI_attribute_other,
                       CFI_type_struct, sizeof (foo_t), rank, extent);
  for (int i = 0; i < arr_len; i++)
    {
      for (int j = 0; j < arr_len; j++)
        {
          foo[i][j].x = (double) (i + 1) * 2 - (double) (j + 1) * 11.;
          foo[i][j].y = (double) (i + 1) * 3 - (double) (j + 1) * 13. +
                        ((double) (i + 1) * 5. - (double) (j + 1) * 17) * I;
        }
    }
  /* Establish c descriptor for the x component. */
  ind = CFI_establish ((CFI_cdesc_t *) &cx, NULL, CFI_attribute_other,
                       CFI_type_double, 0, rank, extent);
  ind = CFI_select_part ((CFI_cdesc_t *) &cx, (CFI_cdesc_t *) &foo_c,
                         offsetof (foo_t, x), 0);
  /* Establish c descriptor for the y component. */
  ind = CFI_establish ((CFI_cdesc_t *) &cy, NULL, CFI_attribute_other,
                       CFI_type_double_Complex, 0, rank, extent);
  ind = CFI_select_part ((CFI_cdesc_t *) &cy, (CFI_cdesc_t *) &foo_c,
                         offsetof (foo_t, y), 0);
  CFI_index_t index[2];
  for (int i = 0; i < arr_len; i++)
    {
      index[0] = i + 1;
      for (int j = 0; j < arr_len; j++)
        {
          index[1] = j + 1;
          if (*(double *) (char *) CFI_address ((CFI_cdesc_t *) &cx, index) !=
              foo[i][j].x)
            {
              printf ("CFI_select_part: failed to properly assign memory and "
                      "dimensional information.\n");
              return 1;
            }
          if (*(double *) (char *) CFI_address ((CFI_cdesc_t *) &cy, index) !=
              creal (foo[i][j].y))
            {
              printf ("CFI_select_part: failed to properly assign memory and "
                      "dimensional information.\n");
              return 1;
            }
          if (*(double *) (char *) (CFI_address ((CFI_cdesc_t *) &cy, index) +
                                    cy.elem_len / 2) != cimag (foo[i][j].y))
            {
              printf ("CFI_select_part: failed to properly assign memory and "
                      "dimensional information.\n");
              return 1;
            }
        }
    }

  rank = 1;
  if (extents != NULL)
    {
      free (extents);
    }
  extents    = malloc (rank * sizeof (CFI_index_t));
  extents[0] = arr_len;
  ind        = CFI_establish ((CFI_cdesc_t *) &foo_c, &foo, CFI_attribute_other,
                       CFI_type_struct, sizeof (foo_t), rank, extents);
  ind = CFI_establish ((CFI_cdesc_t *) &cx, NULL, CFI_attribute_allocatable,
                       CFI_type_double, 0, rank, extents);
  ind = CFI_select_part ((CFI_cdesc_t *) &cx, (CFI_cdesc_t *) &foo_c,
                         offsetof (foo_t, x), 0);
  if (ind != CFI_INVALID_ATTRIBUTE)
    {
      printf ("CFI_select_part: failed to detect invalid attribute.\n");
      return 1;
    }

  ind = CFI_establish ((CFI_cdesc_t *) &foo_c, NULL, CFI_attribute_other,
                       CFI_type_struct, sizeof (foo_t), rank, extents);
  ind = CFI_establish ((CFI_cdesc_t *) &cx, NULL, CFI_attribute_other,
                       CFI_type_double, 0, rank, extents);
  ind = CFI_select_part ((CFI_cdesc_t *) &cx, (CFI_cdesc_t *) &foo_c,
                         offsetof (foo_t, x), 0);
  if (ind != CFI_ERROR_BASE_ADDR_NULL)
    {
      printf ("CFI_select_part: failed to detect that base address of source "
              "is NULL.\n");
      return 1;
    }

  ind = CFI_establish ((CFI_cdesc_t *) &foo_c, &foo, CFI_attribute_other,
                       CFI_type_struct, sizeof (foo_t), rank + 1, extents);
  ind = CFI_establish ((CFI_cdesc_t *) &cx, NULL, CFI_attribute_other,
                       CFI_type_double, 0, rank, extents);
  ind = CFI_select_part ((CFI_cdesc_t *) &cx, (CFI_cdesc_t *) &foo_c,
                         offsetof (foo_t, x), 0);
  if (ind != CFI_INVALID_RANK)
    {
      printf ("CFI_select_part: failed to detect invalid rank.\n");
      return 1;
    }

  extents[0] = -1;
  ind        = CFI_establish ((CFI_cdesc_t *) &foo_c, &foo, CFI_attribute_other,
                       CFI_type_struct, sizeof (foo_t), rank, extents);
  extents[0] = arr_len;
  ind        = CFI_establish ((CFI_cdesc_t *) &cx, NULL, CFI_attribute_other,
                       CFI_type_double, 0, rank, extents);
  ind = CFI_select_part ((CFI_cdesc_t *) &cx, (CFI_cdesc_t *) &foo_c,
                         offsetof (foo_t, x), 0);
  if (ind != CFI_INVALID_DESCRIPTOR)
    {
      printf ("CFI_select_part: failed to detect that source is an assumed "
              "size array.\n");
      return 1;
    }

  ind = CFI_establish ((CFI_cdesc_t *) &foo_c, &foo, CFI_attribute_other,
                       CFI_type_struct, sizeof (foo_t), rank, extents);
  ind = CFI_establish ((CFI_cdesc_t *) &cx, NULL, CFI_attribute_other,
                       CFI_type_double, 0, rank, extents);
  ind = CFI_select_part ((CFI_cdesc_t *) &cx, (CFI_cdesc_t *) &foo_c, -1, 0);
  if (ind != CFI_ERROR_OUT_OF_BOUNDS)
    {
      printf (
          "CFI_select_part: failed to detect out of bounds displacement.\n");
      return 1;
    }
  ind = CFI_select_part ((CFI_cdesc_t *) &cx, (CFI_cdesc_t *) &foo_c,
                         foo_c.elem_len, 0);
  if (ind != CFI_ERROR_OUT_OF_BOUNDS)
    {
      printf ("CFI_select_part: failed to detect out of bounds size.\n");
      return 1;
    }

  ind = CFI_establish ((CFI_cdesc_t *) &foo_c, &foo, CFI_attribute_other,
                       CFI_type_struct, sizeof (foo_t), rank, extent);
  ind = CFI_establish ((CFI_cdesc_t *) &cx, NULL, CFI_attribute_other,
                       CFI_type_double, 0, rank, extent);
  ind = CFI_select_part ((CFI_cdesc_t *) &cx, (CFI_cdesc_t *) &foo_c,
                         foo_c.elem_len - 1, 0);
  if (ind != CFI_ERROR_OUT_OF_BOUNDS)
    {
      printf ("CFI_select_part: failed to detect displacement plus element length go beyond the structure bounds.\n");
      return 1;
    }

  return 0;
}
