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

#include "ISO_Fortran_binding.h"
#include <stdio.h>
#include <stdlib.h>

/* Functions. */
int CFI_establish (CFI_cdesc_t *dv, void *base_addr, CFI_attribute_t attribute,
                   CFI_type_t type, size_t elem_len, CFI_rank_t rank,
                   const CFI_index_t extents[])
{
  /* C descriptor must not be NULL. */
  if (dv == NULL)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_establish: C descriptor is "
                       "NULL. (Error No. %d).\n",
               CFI_INVALID_DESCRIPTOR);
      return CFI_INVALID_DESCRIPTOR;
    }

  /* Rank must be between 0 and CFI_MAX_RANK. */
  if (rank < 0 || rank > CFI_MAX_RANK)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_establish: Rank must be "
                       "between 0 and %d, 0 < rank (0 !< %d). (Error No. "
                       "%d).\n",
               CFI_MAX_RANK, rank, CFI_INVALID_RANK);
      return CFI_INVALID_RANK;
    }

  /* C Descriptor must not be an allocated allocatable. */
  if (dv->attribute == CFI_attribute_allocatable && dv->base_addr != NULL)
    {
      fprintf (stderr,
               "ISO_Fortran_binding.c: CFI_establish: If the C Descriptor "
               "represents an allocatable variable (dv->attribute = %d), its "
               "base address must be NULL (dv->base_addr = NULL). (Error No. "
               "%d).\n",
               CFI_attribute_allocatable, CFI_INVALID_DESCRIPTOR);
      return CFI_INVALID_DESCRIPTOR;
    }

  /* If base address is not NULL, the established C Descriptor is for a
   * nonallocatable entity. */
  if (attribute == CFI_attribute_allocatable && base_addr != NULL)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_establish: If base address "
                       "is not NULL (base_addr != NULL), the established C "
                       "Descriptor is for a nonallocatable entity (attribute "
                       "!= %d). (Error No. %d).\n",
               CFI_attribute_allocatable, CFI_INVALID_ATTRIBUTE);
      return CFI_INVALID_ATTRIBUTE;
    }

  dv->base_addr = base_addr;

  /* elem_len is only used if the item is not a type with a kind parameter. */
  if (type == CFI_type_char || type == CFI_type_ucs4_char ||
      type == CFI_type_signed_char || type == CFI_type_struct ||
      type == CFI_type_other)
    {
      dv->elem_len = elem_len;
    }
  else
    {
      /* base_type describes the intrinsic type with kind parameter. */
      size_t base_type = type & CFI_type_mask;
      /* base_type_size is the size in bytes of the variable as given by its
       * kind parameter. */
      size_t base_type_size = (type - base_type) >> CFI_type_kind_shift;
      /* Kind types 10 have a size of 64 bytes. */
      if (base_type_size == 10)
        {
          base_type_size = 64;
        }
      /* Complex numbers are twice the size of their real counterparts. */
      if (base_type == CFI_type_Complex)
        {
          base_type_size *= 2;
        }
      dv->elem_len = base_type_size;
    }

  dv->version   = CFI_VERSION;
  dv->rank      = rank;
  dv->attribute = attribute;
  dv->type      = type;

  /* Extents must not be NULL if rank is greater than zero and base_addr is not
   * NULL */
  if (rank > 0 && base_addr != NULL)
    {
      if (extents == NULL)
        {
          fprintf (stderr, "ISO_Fortran_binding.c: CFI_establish: Extents must "
                           "not be NULL (extents != NULL) if rank (= %d) > 0 "
                           "and base address is not NULL (base_addr != NULL). "
                           "(Error No. %d).\n",
                   rank, CFI_INVALID_EXTENT);
          return CFI_INVALID_EXTENT;
        }
      for (int i = 0; i < rank; i++)
        {
          /* If the C Descriptor is for a pointer then the lower bounds of every
           * dimension are set to zero. */
          if (attribute == CFI_attribute_pointer)
            {
              dv->dim[i].lower_bound = 0;
            }
          else
            {
              dv->dim[i].lower_bound = 1;
            }
          dv->dim[i].extent = extents[i];
          dv->dim[i].sm     = dv->elem_len;
        }
    }

  return CFI_SUCCESS;
}

int CFI_setpointer (CFI_cdesc_t *result, CFI_cdesc_t *source,
                    const CFI_index_t lower_bounds[])
{
  /* Result must not be NULL. */
  if (result == NULL)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_setpointer: Result is NULL. "
                       "(Error No. %d).\n",
               CFI_INVALID_DESCRIPTOR);
      return CFI_INVALID_DESCRIPTOR;
    }

  /* If source is NULL, the result is a C Descriptor that describes a
   * disassociated pointer. */
  if (source == NULL)
    {
      result->base_addr = NULL;
      result->version   = CFI_VERSION;
      result->attribute = CFI_attribute_pointer;
    }
  else
    {
      /* Check that element lengths, ranks and types of source and result are
       * the same. */
      if (result->elem_len != source->elem_len)
        {
          fprintf (stderr, "ISO_Fortran_binding.c: CFI_setpointer: Element "
                           "lengths of result (result->elem_len = %ld) and "
                           "source (source->elem_len = %ld) must be the same. "
                           "(Error No. %d).\n",
                   result->elem_len, source->elem_len, CFI_INVALID_ELEM_LEN);
          return CFI_INVALID_ELEM_LEN;
        }

      if (result->rank != source->rank)
        {
          fprintf (stderr, "ISO_Fortran_binding.c: CFI_setpointer: Ranks of "
                           "result (result->rank = %d) and source "
                           "(source->rank = %d) must be the same. (Error "
                           "No. %d).\n",
                   result->rank, source->rank, CFI_INVALID_RANK);
          return CFI_INVALID_RANK;
        }

      if (result->type != source->type)
        {
          fprintf (stderr, "ISO_Fortran_binding.c: CFI_setpointer: Types of "
                           "result (result->type = %d) and source "
                           "(source->type = %d) must be the same. (Error "
                           "No. %d).\n",
                   result->type, source->type, CFI_INVALID_TYPE);
          return CFI_INVALID_TYPE;
        }

      /* If the source is a disassociated pointer, the result must also describe
       * a disassociated pointer. */
      if (source->base_addr == NULL &&
          source->attribute == CFI_attribute_pointer)
        {
          result->base_addr = NULL;
        }
      else
        {
          result->base_addr = source->base_addr;
        }
      /* Assign components to result. */
      result->version   = source->version;
      result->attribute = source->attribute;

      /* Dimension information. */
      for (int i = 0; i < source->rank; i++)
        {
          if (lower_bounds != NULL)
            {
              result->dim[i].lower_bound = lower_bounds[i];
            }
          else
            {
              result->dim[i].lower_bound = source->dim[i].lower_bound;
            }
          result->dim[i].extent = source->dim[i].extent;
          result->dim[i].sm     = source->dim[i].sm;
        }
    }

  return CFI_SUCCESS;
}

void *CFI_address (const CFI_cdesc_t *dv, const CFI_index_t subscripts[])
{
  /* C Descriptor must not be NULL. */
  if (dv == NULL)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_address: C Descriptor is "
                       "NULL. (Error No. %d).\n",
               CFI_INVALID_DESCRIPTOR);
      return NULL;
    }

  /* Base address of C Descriptor must not be NULL. */
  if (dv->base_addr == NULL)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_address: base address of C "
                       "Descriptor must not be NULL. (Error No. %d).\n",
               CFI_ERROR_BASE_ADDR_NULL);
      return NULL;
    }

  /* Return base address if C descriptor is a scalar. */
  if (dv->rank == 0)
    {
      return dv->base_addr;
    }
  /* Calculate the appropriate base address if dv is not a scalar. */
  else
    {
      /* Base address is the C address of the element of the object specified by
       * subscripts. */
      void *base_addr;

      /* In order to properly account for Fortran's column major order we need
       * to transpose the subscripts, since columns are stored contiguously as
       * opposed to rows like C. */
      CFI_index_t *tr_subscripts;
      CFI_dim_t *  tr_dim;
      tr_subscripts = malloc (dv->rank * sizeof (CFI_index_t));
      tr_dim        = malloc (dv->rank * sizeof (CFI_dim_t));
      for (int i = 0; i < dv->rank; i++)
        {
          CFI_index_t idx  = dv->rank - i - 1;
          tr_subscripts[i] = subscripts[idx];
          tr_dim[i]        = dv->dim[idx];
          /* Normalise the subscripts to start counting the address from 0. */
          tr_subscripts[i] -= tr_dim[i].lower_bound;
        }

      /* We assume column major order as that is how Fortran stores arrays. We
       * calculate the memory address of the specified element via the canonical
       * array dimension reduction map and multiplying by the memory stride. */
      CFI_index_t index = tr_subscripts[0] * tr_dim[0].sm;
      /* Check that the first subscript is within the bounds of the Fortran
       * array. */
      if (subscripts[0] < dv->dim[0].lower_bound ||
          subscripts[0] > dv->dim[0].lower_bound + dv->dim[0].extent - 1)
        {
          fprintf (stderr, "ISO_Fortran_binding.c: CFI_address: subscripts[0], "
                           "is out of bounds. dim->[0].lower_bound <= "
                           "subscripts[0] <= dv->dim[0].lower_bound + "
                           "dv->dim[0].extent - 1 (%ld <= %ld <= %ld). (Error "
                           "No. %d).\n",
                   dv->dim[0].lower_bound, subscripts[0],
                   dv->dim[0].lower_bound + dv->dim[0].extent - 1,
                   CFI_ERROR_OUT_OF_BOUNDS);
          return NULL;
        }

      /* Start calculating the memory offset. We use the transposed subscripts
       * because we assume the array is coming from Fortran and the address is
       * being queried in column-major order. */
      CFI_index_t tmp_index = 1;
      for (int i = 1; i < dv->rank; i++)
        {
          /* Check that the subsequent subscripts are within the bounds of the
           * Fortran array. */
          if (subscripts[i] < dv->dim[i].lower_bound ||
              subscripts[i] > dv->dim[i].lower_bound + dv->dim[i].extent - 1)
            {
              fprintf (stderr, "ISO_Fortran_binding.c: CFI_address: "
                               "subscripts[%d], is out of bounds. "
                               "dv->dim[%d].lower_bound <= subscripts[%d] <= "
                               "dv->dim[%d].lower_bound + dv->dim[%d].extent - "
                               "1 (%ld <= %ld <= %ld). (Error No. %d).\n",
                       i, i, i, i, i, dv->dim[i].lower_bound, subscripts[i],
                       dv->dim[i].extent + dv->dim[i].lower_bound - 1,
                       CFI_ERROR_OUT_OF_BOUNDS);
              return NULL;
            }

          /* Use the canonical dimension reduction mapping to find the memory
           * address of the relevant subscripts. It is assumed the arrays are
           * stored in column-major order like in Fortran, and the provided
           * subscripts are given as if we were operating on a Fortran array. */
          tmp_index *=
              tr_subscripts[i] * tr_dim[i - 1].extent * tr_dim[i - 1].sm;
          index += tmp_index;
        }
      free (tr_subscripts);
      free (tr_dim);

      /* There's no way in C to do general arithmetic on a void pointer so we
       * cast to a char pointer, do the arithmetic and cast back to a
       * void pointer. */
      base_addr = (char *) dv->base_addr + index;

      return base_addr;
    }
}

int CFI_is_contiguous (const CFI_cdesc_t *dv)
{
  /* C descriptor must not be NULL. */
  if (dv == NULL)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_is_contiguous: C descriptor "
                       "is NULL. (Error No. %d).\n",
               CFI_INVALID_DESCRIPTOR);
      return CFI_INVALID_DESCRIPTOR;
    }

  /* Base address must not be NULL. */
  if (dv->base_addr == NULL)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_is_contiguous: Base address "
                       "of C Descriptor is already NULL. (Error No. %d).\n",
               CFI_ERROR_BASE_ADDR_NULL);
      return CFI_ERROR_BASE_ADDR_NULL;
    }

  /* Must be an array. */
  if (dv->rank == 0)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_is_contiguous: C Descriptor "
                       "must describe an array (0 < dv->rank = %d). (Error No. "
                       "%d).\n",
               dv->rank, CFI_INVALID_RANK);
      return CFI_INVALID_RANK;
    }

  /* If an array is not contiguous the memory stride is different to the element
   * length. */
  for (int i = 0; i < dv->rank; i++)
    {
      if (dv->dim[i].sm != dv->elem_len)
        return 0;
    }

  /* Allocatable arrays are always contiguous. */
  if (dv->attribute == CFI_attribute_allocatable)
    {
      return 1;
    }
  else
    {
      return 0;
    }
}

int CFI_allocate (CFI_cdesc_t *dv, const CFI_index_t lower_bounds[],
                  const CFI_index_t upper_bounds[], size_t elem_len)
{
  /* C Descriptor must not be NULL. */
  if (dv == NULL)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_allocate: C Descriptor is "
                       "NULL. (Error No. %d).\n",
               CFI_INVALID_DESCRIPTOR);
      return CFI_INVALID_DESCRIPTOR;
    }

  /* Base address of C Descriptor must be NULL. */
  if (dv->base_addr != NULL)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_allocate: Base address of C "
                       "Descriptor must be NULL. (Error No. %d).\n",
               CFI_ERROR_BASE_ADDR_NOT_NULL);
      return CFI_ERROR_BASE_ADDR_NOT_NULL;
    }

  /* The C Descriptor must be for an allocatable or pointer object. */
  if (dv->attribute == CFI_attribute_other)
    {
      fprintf (stderr,
               "ISO_Fortran_binding.c: CFI_allocate: The object of the C "
               "Descriptor must be a pointer or allocatable variable. "
               "(Error No. %d).\n",
               CFI_INVALID_ATTRIBUTE);
      return CFI_INVALID_ATTRIBUTE;
    }

  /* If the type is a character, the descriptor's element length is replaced
   * by the elem_len argument. */
  if (dv->type == CFI_type_char || dv->type == CFI_type_ucs4_char ||
      dv->type == CFI_type_signed_char)
    {
      dv->elem_len = elem_len;
    }

  /* Dimension information and calculating the array length. */
  size_t arr_len = 1;
  /* If rank is greater than 0, lower_bounds and upper_bounds are used. They're
   * ignored otherwhise. */
  if (dv->rank > 0)
    {
      if (lower_bounds == NULL || upper_bounds == NULL)
        {
          fprintf (stderr, "ISO_Fortran_binding.c: CFI_allocate: If 0 < rank "
                           "(= %d) upper_bounds[] and lower_bounds[], must not "
                           "be NULL. (Error No. %d).\n",
                   dv->rank, CFI_INVALID_EXTENT);
          return CFI_INVALID_EXTENT;
        }
      for (int i = 0; i < dv->rank; i++)
        {
          dv->dim[i].lower_bound = lower_bounds[i];
          dv->dim[i].extent      = upper_bounds[i] - dv->dim[i].lower_bound + 1;
          dv->dim[i].sm          = dv->elem_len;
          arr_len *= dv->dim[i].extent;
        }
    }

  dv->base_addr = calloc (arr_len, dv->elem_len);
  if (dv->base_addr == NULL)
    {
      printf ("ISO_Fortran_binding.c: CFI_allocate: Failure in memory "
              "allocation. (Error no. %d).\n",
              CFI_ERROR_MEM_ALLOCATION);
      return CFI_ERROR_MEM_ALLOCATION;
    }

  return CFI_SUCCESS;
}

int CFI_deallocate (CFI_cdesc_t *dv)
{
  /* C Descriptor must not be NULL */
  if (dv == NULL)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_deallocate: C Descriptor. "
                       "is NULL. (Error No. %d).\n",
               CFI_INVALID_DESCRIPTOR);
      return CFI_INVALID_DESCRIPTOR;
    }

  /* Base address must not be NULL. */
  if (dv->base_addr == NULL)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_deallocate: Base address is "
                       "NULL already. (Error No. %d).\n",
               CFI_ERROR_BASE_ADDR_NULL);
      return CFI_ERROR_BASE_ADDR_NULL;
    }

  /* C Descriptor must be for an allocatable or pointer variable. */
  if (dv->attribute == CFI_attribute_other)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_deallocate: C Descriptor "
                       "must describe a pointer or allocatabale object. (Error "
                       "No. %d).\n",
               CFI_INVALID_ATTRIBUTE);
      return CFI_INVALID_ATTRIBUTE;
    }

  /* Free and nullify memory. */
  free (dv->base_addr);
  dv->base_addr = NULL;

  return CFI_SUCCESS;
}

int CFI_section (CFI_cdesc_t *result, const CFI_cdesc_t *source,
                 const CFI_index_t lower_bounds[],
                 const CFI_index_t upper_bounds[], const CFI_index_t strides[])
{
  /* C Descriptors must not be NULL. */
  if (source == NULL)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_section: Source must not be "
                       "NULL. (Error No. %d).\n",
               CFI_INVALID_DESCRIPTOR);
      return CFI_INVALID_DESCRIPTOR;
    }
  if (result == NULL)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_section: Result must not be "
                       "NULL. (Error No. %d).\n",
               CFI_INVALID_DESCRIPTOR);
      return CFI_INVALID_DESCRIPTOR;
    }

  /* Base address of source must not be NULL. */
  if (source->base_addr == NULL)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_section: Base address of "
                       "source must not be NULL. (Error No. %d).\n",
               CFI_ERROR_BASE_ADDR_NULL);
      return CFI_ERROR_BASE_ADDR_NULL;
    }

  /* Result must not be an allocatable array. */
  if (result->attribute == CFI_attribute_allocatable)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_section: Result must not "
                       "describe an allocatable array. (Error No. %d).\n",
               CFI_INVALID_ATTRIBUTE);
      return CFI_INVALID_ATTRIBUTE;
    }

  /* Source must be some form of array (nonallocatable nonpointer array,
   * allocated allocatable array or an associated pointer array). */
  if (source->rank <= 0)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_section: Source must "
                       "describe an array (0 < source->rank, 0 !< %d). (Error No. "
                       "%d).\n",
               source->rank, CFI_INVALID_RANK);
      return CFI_INVALID_RANK;
    }

  /* Element lengths of source and result must be equal. */
  if (result->elem_len != source->elem_len)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_section: The element "
                       "lengths of source (source->elem_len = %ld) and result "
                       "(result->elem_len = %ld) must be equal. (Error No. "
                       "%d).\n",
               source->elem_len, result->elem_len, CFI_INVALID_ELEM_LEN);
      return CFI_INVALID_ELEM_LEN;
    }

  /* Types must be equal. */
  if (result->type != source->type)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_section: Types of source "
                       "(source->type = %d) and result (result->type = %d) "
                       "must be equal. (Error No. %d).\n",
               source->type, result->type, CFI_INVALID_TYPE);
      return CFI_INVALID_TYPE;
    }

  /* Stride of zero in the i'th dimension means rank reduction in that
   * dimension. */
  int zero_count = 0;
  for (int i = 0; i < source->rank; i++)
    {
      if (strides[i] == 0)
        {
          zero_count++;
        }
    }

  /* Rank of result must be equal the the rank of source minus the number of
   * zeros in strides. */
  if (result->rank != source->rank - zero_count)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_section: Rank of result "
                       "must be equal to the rank of source minus the number "
                       "of zeros in strides (result->rank = source->rank - "
                       "zero_count, %d != %d - %d) (Error No. %d).\n",
               result->rank, source->rank, zero_count, CFI_INVALID_RANK);
      return CFI_INVALID_RANK;
    }

  /* Dimension information. */
  CFI_index_t *lower;
  CFI_index_t *upper;
  CFI_index_t *stride;
  lower  = malloc (source->rank * sizeof (CFI_index_t));
  upper  = malloc (source->rank * sizeof (CFI_index_t));
  stride = malloc (source->rank * sizeof (CFI_index_t));

  /* Lower bounds. */
  if (lower_bounds == NULL)
    {
      for (int i = 0; i < source->rank; i++)
        {
          lower[i] = source->dim[i].lower_bound;
        }
    }
  else
    {
      for (int i = 0; i < source->rank; i++)
        {
          lower[i] = lower_bounds[i];
        }
    }

  /* Upper bounds. */
  if (upper_bounds == NULL)
    {
      if (source->dim[source->rank].extent == -1)
        {
          fprintf (stderr,
                   "ISO_Fortran_binding.c: CFI_section: Source must not "
                   "be an assumed size array if upper_bounds is NULL. (Error "
                   "No. %d).\n",
                   CFI_INVALID_EXTENT);
          return CFI_INVALID_EXTENT;
        }
      for (int i = 0; i < source->rank; i++)
        {
          upper[i] = source->dim[i].lower_bound + source->dim[i].extent - 1;
        }
    }
  else
    {
      for (int i = 0; i < source->rank; i++)
        {
          upper[i] = upper_bounds[i];
        }
    }

  /* Stride */
  if (strides == NULL)
    {
      for (int i = 0; i < source->rank; i++)
        {
          stride[i] = 1;
        }
    }
  else
    {
      for (int i = 0; i < source->rank; i++)
        {
          stride[i] = strides[i];
          /* If stride[i] = then lower[i] and upper[i] must be equal. */
          if (stride[i] == 0 && lower[i] != upper[i])
            {
              fprintf (stderr, "ISO_Fortran_binding.c: CFI_section: If "
                               "strides[%d] = 0, then the lower bounds, "
                               "lower_bounds[%d] = %ld, and upper_bounds[%d] = "
                               "%ld, must be equal. (Error No. %d).\n",
                       i, i, lower_bounds[i], i, upper_bounds[i],
                       CFI_ERROR_OUT_OF_BOUNDS);
              return CFI_ERROR_OUT_OF_BOUNDS;
            }
        }
    }

  /* Check that section upper and lower bounds are within the array bounds. */
  for (int i = 0; i < source->rank; i++)
    {
      if (lower_bounds != NULL &&
          (lower[i] < source->dim[i].lower_bound ||
           lower[i] > source->dim[i].lower_bound + source->dim[i].extent - 1))
        {
          fprintf (stderr, "ISO_Fortran_binding.c: CFI_section: Lower bounds "
                           "must be within the bounds of the fortran array "
                           "(source->dim[%d].lower_bound <= lower_bounds[%d] "
                           "<= source->dim[%d].lower_bound + "
                           "source->dim[%d].extent - 1, %ld <= %ld <= %ld). "
                           "(Error No. %d).\n",
                   i, i, i, i, source->dim[i].lower_bound, lower[i],
                   source->dim[i].lower_bound + source->dim[i].extent - 1,
                   CFI_ERROR_OUT_OF_BOUNDS);
          return CFI_ERROR_OUT_OF_BOUNDS;
        }
      if (upper_bounds != NULL &&
          (upper[i] < source->dim[i].lower_bound ||
           upper[i] > source->dim[i].lower_bound + source->dim[i].extent - 1))
        {
          fprintf (stderr, "ISO_Fortran_binding.c: CFI_section: Upper bounds "
                           "must be within the bounds of the fortran array "
                           "(source->dim[%d].lower_bound <= upper_bounds[%d] "
                           "<= source->dim[%d].lower_bound + "
                           "source->dim[%d].extent - 1, %ld !<= %ld !<= %ld). "
                           "(Error No. %d).\n",
                   i, i, i, i, source->dim[i].lower_bound, upper[i],
                   source->dim[i].lower_bound + source->dim[i].extent - 1,
                   CFI_ERROR_OUT_OF_BOUNDS);
          return CFI_ERROR_OUT_OF_BOUNDS;
        }
      if (upper[i] < lower[i] && stride[i] >= 0)
        {
          fprintf (stderr, "ISO_Fortran_binding.c: CFI_section: If the upper "
                           "bound is smaller than the lower bound for a given "
                           "dimension (upper[%d] < lower[%d], %ld < %ld), then "
                           "the stride for said dimension must be negative "
                           "(stride[%d] < 0, %ld < 0). (Error No. %d)\n",
                   i, i, upper[i], lower[i], i, stride[i], CFI_INVALID_STRIDE);
          return CFI_INVALID_STRIDE;
        }
    }

  /* Update the result to describe the array section. */
  /* Set appropriate memory address. */
  result->base_addr = CFI_address (source, lower);

  /* Set the appropriate dimension information that gives us access to the
   * data. */
  int aux = 0;
  for (int i = 0; i < source->rank; i++)
    {
      if (stride[i] == 0)
        {
          aux++;
          continue;
        }
      int idx                      = i - aux;
      result->dim[idx].lower_bound = lower[i];
      result->dim[idx].extent      = upper[i] - lower[i] + 1;
      result->dim[idx].sm          = stride[i] * source->dim[i].sm;
    }

  free (lower);
  free (upper);
  free (stride);

  return CFI_SUCCESS;
}

int CFI_select_part (CFI_cdesc_t *result, const CFI_cdesc_t *source,
                     size_t displacement, size_t elem_len)
{
  /* C Descriptors must not be NULL. */
  if (source == NULL)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_select_part: Source must "
                       "not be NULL. (Error No. %d).\n",
               CFI_INVALID_DESCRIPTOR);
      return CFI_INVALID_DESCRIPTOR;
    }
  if (result == NULL)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_select_part: Result must "
                       "not be NULL. (Error No. %d).\n",
               CFI_INVALID_DESCRIPTOR);
      return CFI_INVALID_DESCRIPTOR;
    }

  /* Attribute of result will be CFI_attribute_other or CFI_attribute_pointer.
   */
  if (result->attribute == CFI_attribute_allocatable)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_select_part: Result must "
                       "not describe an allocatabale object (result->attribute "
                       "!= %d). (Error No. %d).\n",
               CFI_attribute_allocatable, CFI_INVALID_ATTRIBUTE);
      return CFI_INVALID_ATTRIBUTE;
    }

  /* Base address of source must not be NULL. */
  if (source->base_addr == NULL)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_select_part: Base address "
                       "of source must not be NULL. (Error No. %d).\n",
               CFI_ERROR_BASE_ADDR_NULL);
      return CFI_ERROR_BASE_ADDR_NULL;
    }

  /* Source and result must have the same rank. */
  if (source->rank != result->rank)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_select_part: Source and "
                       "result must have the same rank (source->rank = %d, "
                       "result->rank = %d). (Error No. %d).\n",
               source->rank, result->rank, CFI_INVALID_RANK);
      return CFI_INVALID_RANK;
    }

  /* Nonallocatable nonpointer must not be an assumed size array. */
  if (source->rank > 0 && source->dim[source->rank - 1].extent == -1)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_select_part: Source must "
                       "not describe an assumed size array "
                       "(source->dim[%d].extent != -1). (Error No. %d).\n",
               source->rank - 1, CFI_INVALID_DESCRIPTOR);
      return CFI_INVALID_DESCRIPTOR;
    }

  /* Element length. */
  if (result->type == CFI_type_char || result->type == CFI_type_ucs4_char ||
      result->type == CFI_type_signed_char)
    {
      result->elem_len = elem_len;
    }

  /* Ensure displacement is within the bounds of the element length of source.
   */
  if (displacement < 0 || displacement > source->elem_len - 1)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_select_part: Displacement "
                       "must be within the bounds of source (0 <= displacement "
                       "<= source->elem_len - 1, 0 <= %ld <= %ld). (Error No. "
                       "%d).\n",
               displacement, source->elem_len - 1, CFI_ERROR_OUT_OF_BOUNDS);
      return CFI_ERROR_OUT_OF_BOUNDS;
    }
  /* Ensure displacement and element length of result are less than or equal to
   * the element length of source. */
  if (displacement + result->elem_len > source->elem_len)
    {
      fprintf (stderr, "ISO_Fortran_binding.c: CFI_select_part: Displacement "
                       "plus the element length of result must be less than or "
                       "equal to the element length of source (displacement + "
                       "result->elem_len <= source->elem_len, %ld + %ld = %ld "
                       "<= %ld). (Error No. %d).\n",
               displacement, result->elem_len, displacement + result->elem_len,
               source->elem_len, CFI_ERROR_OUT_OF_BOUNDS);
      return CFI_ERROR_OUT_OF_BOUNDS;
    }
  if (result->rank > 0)
    {
      for (int i = 0; i < result->rank; i++)
        {
          result->dim[i].lower_bound = source->dim[i].lower_bound;
          result->dim[i].extent      = source->dim[i].extent;
          result->dim[i].sm =
              source->dim[i].sm +
              displacement * (source->dim[i].sm / source->elem_len - 1);
        }
    }

  result->base_addr = (char *) source->base_addr + displacement;
  return CFI_SUCCESS;
}
