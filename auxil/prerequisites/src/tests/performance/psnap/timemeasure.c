/* PSNAP Test: timemeausure.c

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
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/errno.h>
#include <sys/time.h>

static struct timezone tz;
static struct timeval start_time, finish_time;

/* Start measuring a time delay */
void start_timer(void)
{
    gettimeofday( &start_time, &tz);
}

/* Retunrn elapsed time in microseconds */
int elapsed_time(void)
{
    gettimeofday( &finish_time, &tz);
    return(1000000.0*(finish_time.tv_sec - start_time.tv_sec) +
           (finish_time.tv_usec - start_time.tv_usec) );
}

/* Return the stopping time in microseconds */
double stop_timer(void)
{
    gettimeofday( &finish_time, &tz);
    return(1000000.0*finish_time.tv_sec + finish_time.tv_usec);
}
