/* syncimages2 test program

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
*/

#include <stdlib.h>
#include <stdio.h>
#include "libcaf.h"
#include <sys/types.h>
#include <sys/errno.h>
#include <sys/time.h>
#include <stdbool.h>


static struct timezone tz;
static struct timeval start_time, finish_time;

/* Start measuring a time delay */
void start_timer(void)
{
    gettimeofday( &start_time, &tz);
}

/* Retunrn elapsed time in milliseconds */
double elapsed_time(void)
{
    gettimeofday( &finish_time, &tz);
    return(1000.0*(finish_time.tv_sec - start_time.tv_sec) +
           (finish_time.tv_usec - start_time.tv_usec)/1000.0 );
}

/* Return the stopping time in milliseconds */
double stop_time(void)
{
    gettimeofday( &finish_time, &tz);
    return(1000.0*finish_time.tv_sec + finish_time.tv_usec/1000.0);
}

int main(int argc, char **argv)
{
  int info = 0, me,np,n=1,i, *images;
  double *a_d,*d;
  caf_token_t token;
  ptrdiff_t size=sizeof(double);
  char errmsg[255];
  bool check = true;

  /* if(argc == 1) */
  /*   { */
  /*     printf("Please insert message size\n"); */
  /*     return 1; */
  /*   } */

  /* sscanf(argv[1],"%d",&n); */

  /* n = (int)n/sizeof(double); */

  /* size = n*sizeof(double); */

  _gfortran_caf_init (&argc, &argv);

  me = _gfortran_caf_this_image (1);
  np = _gfortran_caf_num_images (1, 1);

  a_d = _gfortran_caf_register(size,CAF_REGTYPE_COARRAY_STATIC,&token,&info,errmsg,255);

  /* start_timer(); */
  images = calloc(np,sizeof(int));

  if(me==1)
    {
      images[0] = -1;
      _gfortran_caf_sync_images(1,images,&info,errmsg,255);
    }
  else
    {
      images[0] = 1;
      _gfortran_caf_sync_images(1,images,&info,errmsg,255);

    }

  printf("proc %d dval: %lf\n",me);

  /* stop_time(); */

  if(info!=0)
    printf("Error\n");

  _gfortran_caf_finalize();

  return 0;
}
