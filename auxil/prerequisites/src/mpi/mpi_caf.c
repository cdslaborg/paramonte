/* One-sided MPI implementation of Libcaf
*
* Copyright (c) 2012-2018, Sourcery, Inc.
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of the Sourcery, Inc., nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED. IN NO EVENT SHALL SOURCERY, INC., BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
* THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. */

/****l* mpi/mpi_caf.c
 * NAME
 *   mpi_caf
 * SYNOPSIS
 *   This program implements the LIBCAF_MPI transport layer.
******
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>     /* For memcpy. */
#include <stdarg.h>     /* For variadic arguments. */
#include <float.h>      /* For type conversion of floating point numbers. */
#ifndef ALLOCA_MISSING
#include <alloca.h>     /* Assume functionality provided elsewhere if missing */
#endif
#include <unistd.h>
#include <stdint.h>     /* For int32_t. */
#include <mpi.h>
#include <pthread.h>
#include <signal.h>     /* For raise */

#ifdef HAVE_MPI_EXT_H
#include <mpi-ext.h>
#endif
#ifdef USE_FAILED_IMAGES
  #define WITH_FAILED_IMAGES 1
#endif

#include "libcaf.h"

/* Define GFC_CAF_CHECK to enable run-time checking. */
/* #define GFC_CAF_CHECK  1 */

/* Debug array referencing  */
static char* caf_array_ref_str[] = {
  "CAF_ARR_REF_NONE",
  "CAF_ARR_REF_VECTOR",
  "CAF_ARR_REF_FULL",
  "CAF_ARR_REF_RANGE",
  "CAF_ARR_REF_SINGLE",
  "CAF_ARR_REF_OPEN_END",
  "CAF_ARR_REF_OPEN_START"
};

static char* caf_ref_type_str[] = {
  "CAF_REF_COMPONENT",
  "CAF_REF_ARRAY",
  "CAF_REF_STATIC_ARRAY",
};

#ifndef EXTRA_DEBUG_OUTPUT
#define dprint(...)
#define chk_err(...)
#else
#define dprint(format, ...)                     \
fprintf(stderr, "%d/%d: %s(%d) " format,        \
        caf_this_image, caf_num_images,         \
        __FUNCTION__, __LINE__, ## __VA_ARGS__)
#define chk_err(ierr)                               \
do                                                  \
{                                                   \
  if (ierr != MPI_SUCCESS)                          \
  {                                                 \
    int err_class, err_len;                         \
    char err_str[MPI_MAX_ERROR_STRING];             \
    MPI_Error_class(ierr, &err_class);              \
    MPI_Error_string(ierr, err_str, &err_len);      \
    dprint("MPI-error: err_class=%d ierr=%d [%s]",  \
           err_class, ierr, err_str);               \
  }                                                 \
} while (0)
#endif

#ifdef GCC_GE_7
/* The caf-token of the mpi-library.
 * Objects of this data structure are owned by the library and are treated as a
 * black box by the compiler.  In the coarray-program the tokens are opaque
 * pointers, i.e. black boxes.
 * 
 * For each coarray (allocatable|save|pointer) (scalar|array|event|lock) a
 * token needs to be present. */
typedef struct mpi_caf_token_t
{
  /* The pointer to memory associated to this token's data on the local image.
   * The compiler uses the address for direct access to the memory of the object
   * this token is assocated to, i.e., the memory pointed to be local_memptr is
   * the scalar or array.
   * When the library is responsible for deleting the memory, then this is the
   * one to free. */
  void *memptr;
  /* The MPI window to associated to the object's data.
   * The window is used to access the data on other images. In pre GCC_GE_7
   * installations this was the token. */
  MPI_Win memptr_win;
  /* The pointer to the primary array, i.e., to coarrays that are arrays and
   * not a derived type. */
  gfc_descriptor_t *desc;
} mpi_caf_token_t;

/* For components of derived type coarrays a slave_token is needed when the
 * component has the allocatable or pointer attribute. The token is reduced in
 * size, because the other data is already accessible and has been read from
 * the remote to fullfill the request.
 * 
 *   TYPE t
 *   +------------------+
 *   | comp *           |
 *   | comp_token *     |
 *   +------------------+
 * 
 *   TYPE(t) : o                struct T // the mpi_caf_token to t
 *                              +----------------+
 *                              | ...            |
 *                              +----------------+
 * 
 *   o[2]%.comp                 // using T to get o of [2]
 * 
 *   +-o-on-image-2----+  "copy" of the requierd parts of o[2] on current image
 *   | 0x4711          |  comp * in global_dynamic_window
 *   | 0x2424          |  comp_token * of type slave_token
 *   +-----------------+
 *   now all required data is present on the current image to access the remote
 *   components. This nests without limit. */
typedef struct mpi_caf_slave_token_t
{
  /* The pointer to the memory associated to this slave token's data on the
   * local image.  When the library is responsible for deleting the memory,
   * then this is the one to free.  And this is the only reason why its stored
   * here. */
  void *memptr;
  /* The pointer to the descriptor or NULL for scalars.
   * When referencing a remote component array, then the extensions of the array
   * are needed. Usually the data pointer is at offset zero of the descriptor_t
   * structure, but we don't rely on it. So store the pointer to the base
   * address of the descriptor. The descriptor always is in the window of the
   * master data or the allocated component and is never stored at an address
   * not accessible by a window. */
  gfc_descriptor_t *desc;
} mpi_caf_slave_token_t;

#define TOKEN(X) &(((mpi_caf_token_t *) (X))->memptr_win)
#else
typedef MPI_Win *mpi_caf_token_t;
#define TOKEN(X) ((mpi_caf_token_t) (X))
#endif

/* Forward declaration of prototype. */

static void terminate_internal (int stat_code, int exit_code)
            __attribute__((noreturn));
static void sync_images_internal (int count, int images[], int *stat,
                                  char *errmsg, size_t errmsg_len,
                                  bool internal);
static void error_stop_str (const char *string, size_t len, bool quiet)
            __attribute__((noreturn));

/* Global variables. */
static int caf_this_image;
static int caf_num_images = 0;
static int caf_is_finalized = 0;
static MPI_Win global_dynamic_win;

#if MPI_VERSION >= 3
  MPI_Info mpi_info_same_size;
#endif // MPI_VERSION

/* The size of pointer on this plattform. */
static const size_t stdptr_size = sizeof(void *);

/* Variables needed for syncing images. */

static int *images_full;
MPI_Request *sync_handles;
static int *arrived;
static const int MPI_TAG_CAF_SYNC_IMAGES = 424242;

/* Pending puts */
#if defined(NONBLOCKING_PUT) && !defined(CAF_MPI_LOCK_UNLOCK)
typedef struct win_sync {
  MPI_Win *win;
  int img;
  struct win_sync *next;
} win_sync;

static win_sync *last_elem = NULL;
static win_sync *pending_puts = NULL;
#endif

/* Linked list of static coarrays registered.  Do not expose to public in the
 * header, because it is implementation specific. */
struct caf_allocated_tokens_t
{
  caf_token_t token;
  struct caf_allocated_tokens_t *prev;
} *caf_allocated_tokens = NULL;

#ifdef GCC_GE_7
/* Linked list of slave coarrays registered. */
struct caf_allocated_slave_tokens_t
{
  mpi_caf_slave_token_t *token;
  struct caf_allocated_slave_tokens_t *prev;
} *caf_allocated_slave_tokens = NULL;
#endif

/* Image status variable */
static int img_status = 0;
static MPI_Win *stat_tok;

/* Active messages variables */
char **buff_am;
MPI_Status *s_am;
MPI_Request *req_am;
MPI_Datatype *dts;
char *msgbody;
pthread_mutex_t lock_am;
int done_am = 0;

char err_buffer[MPI_MAX_ERROR_STRING];

/* All CAF runtime calls should use this comm instead of MPI_COMM_WORLD for
 * interoperability purposes. */
MPI_Comm CAF_COMM_WORLD;

static caf_teams_list *teams_list = NULL;
static caf_used_teams_list *used_teams = NULL;

/* Emitted when a theorectically unreachable part is reached. */
const char unreachable[] = "Fatal error: unreachable alternative found.\n";

#ifdef WITH_FAILED_IMAGES
/* The stati of the other images.  image_stati is an array of size
 * caf_num_images at the beginning the status of each image is noted here where
 * the index is the image number minus one. */
int *image_stati;

/* This gives the number of all images that are known to have failed. */
int num_images_failed = 0;

/* This is the number of all images that are known to have stopped. */
int num_images_stopped = 0;

/* The async. request-handle to all participating images. */
MPI_Request alive_request;

/* This dummy is used for the alive request.  Its content is arbitrary and
 * never read.  Its just a memory location where one could put something,
 * which is never done. */
int alive_dummy;

/* The mpi error-handler object associate to CAF_COMM_WORLD. */
MPI_Errhandler failed_CAF_COMM_mpi_err_handler;

/* The monitor comm for detecting failed images. We can not attach the monitor
 * to CAF_COMM_WORLD or the messages send by sync images would be caught by the
 * monitor. */
MPI_Comm alive_comm;

/* Set when entering a sync_images_internal, to prevent the error handler from
 * eating our messages. */
bool no_stopped_images_check_in_errhandler = 0;
#endif

/* For MPI interoperability, allow external initialization
 * (and thus finalization) of MPI. */
bool caf_owns_mpi = false;

/* Foo function pointers for coreduce.
 * The handles when arguments are passed by reference. */
int (*int8_t_by_reference)(void *, void *);
int (*int16_t_by_reference)(void *, void *);
int (*int32_t_by_reference)(void *, void *);
int (*int64_t_by_reference)(void *, void *);
float (*float_by_reference)(void *, void *);
double (*double_by_reference)(void *, void *);
/* Strings are always passed by reference. */
void (*char_by_reference)(void *, int, void *, void *, int, int);
/* The handles when arguments are passed by value. */
int8_t (*int8_t_by_value)(int8_t, int8_t);
int16_t (*int16_t_by_value)(int16_t, int16_t);
int (*int32_t_by_value)(int32_t, int32_t);
int64_t (*int64_t_by_value)(int64_t, int64_t);
float (*float_by_value)(float, float);
double (*double_by_value)(double, double);

/* Define shortcuts for Win_lock and _unlock depending on whether the primitives
 * are available in the MPI implementation.  When they are not available the
 * shortcut is expanded to nothing by the preprocessor else to the API call.
 * This prevents having #ifdef #else #endif constructs strewn all over the code
 * reducing its readability. */
#ifdef CAF_MPI_LOCK_UNLOCK
#define CAF_Win_lock(type, img, win) MPI_Win_lock (type, img, 0, win)
#define CAF_Win_unlock(img, win) MPI_Win_unlock (img, win)
#define CAF_Win_lock_all(win)
#define CAF_Win_unlock_all(win)
#else // CAF_MPI_LOCK_UNLOCK
#define CAF_Win_lock(type, img, win)
#define CAF_Win_unlock(img, win) MPI_Win_flush (img, win)
#if MPI_VERSION >= 3
#define CAF_Win_lock_all(win) MPI_Win_lock_all (MPI_MODE_NOCHECK, win)
#else
#define CAF_Win_lock_all(win)
#endif
#define CAF_Win_unlock_all(win) MPI_Win_unlock_all (win)
#endif // CAF_MPI_LOCK_UNLOCK

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

#if defined(NONBLOCKING_PUT) && !defined(CAF_MPI_LOCK_UNLOCK)
void explicit_flush()
{
  win_sync *w = pending_puts, *t;
  MPI_Win *p;
  int ierr;
  while (w != NULL)
  {
    p = w->win;
    ierr = MPI_Win_flush(w->img,*p); chk_err(ierr);
    t = w;
    w = w->next;
    free(t);
  }
  last_elem = NULL;
  pending_puts = NULL;
}
#endif

#ifdef HELPER
void helperFunction()
{
  int i = 0, flag = 0, msgid = 0, ierr;
  int ndim = 0, position = 0;

  s_am = calloc(caf_num_images, sizeof(MPI_Status));
  req_am = calloc(caf_num_images, sizeof(MPI_Request));
  dts = calloc(caf_num_images, sizeof(MPI_Datatype));

  for (i = 0; i < caf_num_images; i++)
  {
    ierr = MPI_Irecv(buff_am[i], 1000, MPI_PACKED, i, 1, CAF_COMM_WORLD,
                     &req_am[i]); chk_err(ierr);
  }

  while (1)
  {
    pthread_mutex_lock(&lock_am);
    for (i = 0; i < caf_num_images; i++)
    {
      if (!caf_is_finalized)
      {
        ierr = MPI_Test(&req_am[i], &flag, &s_am[i]); chk_err(ierr);
        if (flag == 1)
        {
          position = 0;
          ierr = MPI_Unpack(buff_am[i], 1000, &position, &msgid, 1, MPI_INT,
                            CAF_COMM_WORLD); chk_err(ierr);
          /* msgid=2 was initially assigned to strided transfers,
           * it can be reused
           * Strided transfers Msgid=2
           * You can add you own function */

          if (msgid == 2)
          {
            msgid = 0; position = 0;
          }
          ierr = MPI_Irecv(buff_am[i], 1000, MPI_PACKED, i, 1, CAF_COMM_WORLD,
                           &req_am[i]); chk_err(ierr);
          flag = 0;
        }
      }
      else
      {
        done_am = 1;
        pthread_mutex_unlock(&lock_am);
        return;
      }
    }
    pthread_mutex_unlock(&lock_am);
  }
}
#endif


/* Keep in sync with single.c. */

static void
caf_runtime_error (const char *message, ...)
{
  va_list ap;
  fprintf(stderr, "Fortran runtime error on image %d: ", caf_this_image);
  va_start(ap, message);
  vfprintf(stderr, message, ap);
  va_end(ap);
  fprintf(stderr, "\n");

  /* FIXME: Shutdown the Fortran RTL to flush the buffer.  PR 43849.
   * FIXME: Do some more effort than just to abort. */
  //  MPI_Finalize();

  /* Should be unreachable, but to make sure also call exit. */
  exit(EXIT_FAILURE);
}

/* Forward declaration of the feature unsupported message for failed images
 * functions. */
static void
unsupported_fail_images_message(const char * functionname);

/* Forward declaration of the feature unimplemented message for allocatable
 * components. */
static void
unimplemented_alloc_comps_message(const char * functionname);

static void
locking_atomic_op(MPI_Win win, int *value, int newval,
                  int compare, int image_index, size_t index)
{
  CAF_Win_lock(MPI_LOCK_EXCLUSIVE, image_index - 1, win);
  int ierr = MPI_Compare_and_swap(&newval, &compare,value, MPI_INT,
                                  image_index - 1, index * sizeof(int), win);
  chk_err(ierr);
  CAF_Win_unlock(image_index - 1, win);
}


/* Define a helper to check whether the image at the given index is healthy,
 * i.e., it hasn't failed. */
#ifdef WITH_FAILED_IMAGES
#define check_image_health(image_index, stat)                   \
if (image_stati[image_index - 1] == STAT_FAILED_IMAGE)          \
{                                                               \
  if (stat == NULL) terminate_internal (STAT_FAILED_IMAGE, 0);  \
  *stat = STAT_FAILED_IMAGE;                                    \
  return;                                                       \
}
#else
#define check_image_health(image_index, stat)
#endif

#ifdef WITH_FAILED_IMAGES
/* Handle failed image's errors and try to recover the remaining process to
 * allow the user to detect an image fail and exit gracefully. */
static void
failed_stopped_errorhandler_function (MPI_Comm* pcomm, int* perr, ...)
{
  MPI_Comm comm, shrunk, newcomm;
  int num_failed_in_group, i, err, ierr;
  MPI_Group comm_world_group, failed_group;
  int *ranks_of_failed_in_comm_world, *ranks_failed;
  int ns, srank, crank, rc, flag, drank, newrank;
  bool stopped = false;

  comm = *pcomm;

  MPI_Error_class(*perr, &err);
  if (err != MPIX_ERR_PROC_FAILED && err != MPIX_ERR_REVOKED)
  {
    /* We can handle PROC_FAILED and REVOKED ones only. */
    char errstr[MPI_MAX_ERROR_STRING];
    int errlen;
    MPI_Error_string(err, errstr, &errlen);
    /* We can't use caf_runtime_error here, because that would exit, which
     * means only the one process will stop, but we need to stop MPI
     * completely, which can be done by calling MPI_Abort(). */
    fprintf(stderr,
            "Fortran runtime error on image #%d:\nMPI error: '%s'.\n",
            caf_this_image, errstr);
    MPI_Abort(*pcomm, err);
  }

  dprint("(error = %d)\n", err);

  ierr = MPIX_Comm_failure_ack(comm); chk_err(ierr);
  ierr = MPIX_Comm_failure_get_acked(comm, &failed_group); chk_err(ierr);
  ierr = MPI_Group_size(failed_group, &num_failed_in_group); chk_err(ierr);

  dprint("%d images failed.\n", num_failed_in_group);
  if (num_failed_in_group <= 0)
  {
    *perr = MPI_SUCCESS;
    return;
  }
  if (num_failed_in_group > caf_num_images)
  {
    *perr = MPI_SUCCESS;
    return;
  }

  ierr = MPI_Comm_group(comm, &comm_world_group); chk_err(ierr);
  ranks_of_failed_in_comm_world =
    (int *) alloca(sizeof(int) * num_failed_in_group);
  ranks_failed = (int *) alloca(sizeof(int) * num_failed_in_group);
  for (i = 0; i < num_failed_in_group; ++i)
  {
    ranks_failed[i] = i;
  }
  /* Now translate the ranks of the failed images into communicator world. */
  ierr = MPI_Group_translate_ranks(failed_group, num_failed_in_group,
                                   ranks_failed, comm_world_group,
                                   ranks_of_failed_in_comm_world);
  chk_err(ierr);

  num_images_failed += num_failed_in_group;

  if (!no_stopped_images_check_in_errhandler)
  {
    int buffer, flag;
    MPI_Request req;
    MPI_Status request_status;
    dprint("Checking for stopped images.\n");
    ierr = MPI_Irecv(&buffer, 1, MPI_INT, MPI_ANY_SOURCE,
                     MPI_TAG_CAF_SYNC_IMAGES, CAF_COMM_WORLD, &req);
    chk_err(ierr);
    if (ierr == MPI_SUCCESS)
    {
      ierr = MPI_Test(&req, &flag, &request_status); chk_err(ierr);
      if (flag)
      {
        // Received a result
        if (buffer == STAT_STOPPED_IMAGE)
        {
          dprint("Image #%d found stopped.\n", request_status.MPI_SOURCE);
          stopped = true;
          if (image_stati[request_status.MPI_SOURCE] == 0)
            ++num_images_stopped;
          image_stati[request_status.MPI_SOURCE] = STAT_STOPPED_IMAGE;
        }
      }
      else
      {
        dprint("No stopped images found.\n");
        ierr = MPI_Cancel(&req); chk_err(ierr);
      }
    }
    else
    {
      int err;
      MPI_Error_class(ierr, &err);
      dprint("Error on checking for stopped images %d.\n", err);
    }
  }

  /* TODO: Consider whether removing the failed image from images_full will be
   * necessary. This is more or less politics. */
  for (i = 0; i < num_failed_in_group; ++i)
  {
    if (ranks_of_failed_in_comm_world[i] >= 0
        && ranks_of_failed_in_comm_world[i] < caf_num_images)
    {
      if (image_stati[ranks_of_failed_in_comm_world[i]] == 0)
        image_stati[ranks_of_failed_in_comm_world[i]] = STAT_FAILED_IMAGE;
    }
    else
    {
      dprint("Rank of failed image %d out of range of images 0..%d.\n",
             ranks_of_failed_in_comm_world[i], caf_num_images);
    }
  }

redo:
  dprint("Before shrink. \n");
  ierr = MPIX_Comm_shrink(*pcomm, &shrunk);
  dprint("After shrink, rc = %d.\n", ierr);
  ierr = MPI_Comm_set_errhandler(shrunk, failed_CAF_COMM_mpi_err_handler);
  chk_err(ierr);
  ierr = MPI_Comm_size(shrunk, &ns); chk_err(ierr);
  ierr = MPI_Comm_rank(shrunk, &srank); chk_err(ierr);
  ierr = MPI_Comm_rank(*pcomm, &crank); chk_err(ierr);

  dprint("After getting ranks, ns = %d, srank = %d, crank = %d.\n",
         ns, srank, crank);

  /* Split does the magic: removing spare processes and reordering ranks
   * so that all surviving processes remain at their former place */
  rc = MPI_Comm_split(shrunk, (crank < 0) ? MPI_UNDEFINED : 1, crank, &newcomm);
  ierr = MPI_Comm_rank(newcomm, &newrank); chk_err(ierr);
  dprint("After split, rc = %d, rank = %d.\n", rc, newrank);
  flag = (rc == MPI_SUCCESS);
  /* Split or some of the communications above may have failed if
   * new failures have disrupted the process: we need to
   * make sure we succeeded at all ranks, or retry until it works. */
  flag = MPIX_Comm_agree(newcomm, &flag);
  dprint("After agree, flag = %d.\n", flag);

  ierr = MPI_Comm_rank(newcomm, &drank); chk_err(ierr);
  dprint("After rank, drank = %d.\n", drank);

  ierr = MPI_Comm_free(&shrunk); chk_err(ierr);
  if (MPI_SUCCESS != flag)
  {
    if (MPI_SUCCESS == rc)
    {
      ierr = MPI_Comm_free(&newcomm); chk_err(ierr);
    }
    goto redo;
  }

  {
    int cmpres;
    ierr = MPI_Comm_compare(*pcomm, CAF_COMM_WORLD, &cmpres); chk_err(ierr);
    dprint("Comm_compare(*comm, CAF_COMM_WORLD, res = %d) = %d.\n",
           cmpres, ierr);
    ierr = MPI_Comm_compare(*pcomm, alive_comm, &cmpres); chk_err(ierr);
    dprint("Comm_compare(*comm, alive_comm, res = %d) = %d.\n",
           cmpres, ierr);
    if (cmpres == MPI_CONGRUENT)
    {
      ierr = MPI_Win_detach(*stat_tok, &img_status); chk_err(ierr);
      dprint("detached win img_status.\n");
      ierr = MPI_Win_free(stat_tok); chk_err(ierr);
      dprint("freed win img_status.\n");
      ierr = MPI_Win_create(&img_status, sizeof(int), 1, mpi_info_same_size,
                            newcomm, stat_tok); chk_err(ierr);
      dprint("(re-)created win img_status.\n");
      CAF_Win_lock_all(*stat_tok);
      dprint("Win_lock_all on img_status.\n");
    }
  }
  /* Also free the old communicator before replacing it. */
  ierr = MPI_Comm_free(pcomm); chk_err(ierr);
  *pcomm = newcomm;

  *perr = stopped ? STAT_STOPPED_IMAGE : STAT_FAILED_IMAGE;
}
#endif

void mutex_lock(MPI_Win win, int image_index, size_t index, int *stat,
                int *acquired_lock, char *errmsg, size_t errmsg_len)
{
  const char msg[] = "Already locked";
#if MPI_VERSION >= 3
  int value = 0, compare = 0, newval = caf_this_image, ierr = 0, i = 0;
#ifdef WITH_FAILED_IMAGES
  int flag, check_failure = 100, zero = 0;
#endif

  if (stat != NULL)
    *stat = 0;

#ifdef WITH_FAILED_IMAGES
  ierr = MPI_Test(&alive_request, &flag, MPI_STATUS_IGNORE); chk_err(ierr);
#endif

  locking_atomic_op(win, &value, newval, compare, image_index, index);

  if (value == caf_this_image && image_index == caf_this_image)
    goto stat_error;

  if (acquired_lock != NULL)
  {
    if (value == 0)
      *acquired_lock = 1;
    else
      *acquired_lock = 0;
    return;
  }

  while (value != 0)
  {
    ++i;
#ifdef WITH_FAILED_IMAGES
    if (i == check_failure)
    {
      i = 1;
      ierr = MPI_Test(&alive_request, &flag, MPI_STATUS_IGNORE); chk_err(ierr);
    }
#endif

    locking_atomic_op(win, &value, newval, compare, image_index, index);
#ifdef WITH_FAILED_IMAGES
    if (image_stati[value] == STAT_FAILED_IMAGE)
    {
      CAF_Win_lock(MPI_LOCK_EXCLUSIVE, image_index - 1, win);
      /* MPI_Fetch_and_op(&zero, &newval, MPI_INT, image_index - 1,
       * index * sizeof(int), MPI_REPLACE, win); */
      ierr = MPI_Compare_and_swap(&zero, &value, &newval, MPI_INT,
                                  image_index - 1, index * sizeof(int), win);
      chk_err(ierr);
      CAF_Win_unlock(image_index - 1, win);
      break;
    }
#else
    usleep(caf_this_image * i);
#endif
  }

  if (stat)
    *stat = ierr;
  else if (ierr == STAT_FAILED_IMAGE)
    terminate_internal(ierr, 0);

  return;

stat_error:
  if (errmsg != NULL)
  {
    memset(errmsg,' ',errmsg_len);
    memcpy(errmsg, msg, MIN(errmsg_len,strlen(msg)));
  }

  if (stat != NULL)
    *stat = 99;
  else
    terminate_internal(99, 1);
#else // MPI_VERSION
#warning Locking for MPI-2 is not implemented
  printf("Locking for MPI-2 is not supported, "
         "please update your MPI implementation\n");
#endif // MPI_VERSION
}

void mutex_unlock(MPI_Win win, int image_index, size_t index, int *stat,
                  char* errmsg, size_t errmsg_len)
{
  const char msg[] = "Variable is not locked";
  if (stat != NULL)
    *stat = 0;
#if MPI_VERSION >= 3
  int value = 1, ierr = 0, newval = 0, flag;
#ifdef WITH_FAILED_IMAGES
  ierr = MPI_Test(&alive_request, &flag, MPI_STATUS_IGNORE); chk_err(ierr);
#endif

  CAF_Win_lock(MPI_LOCK_EXCLUSIVE, image_index - 1, win);
  ierr = MPI_Fetch_and_op(&newval, &value, MPI_INT, image_index - 1,
                          index * sizeof(int), MPI_REPLACE, win); chk_err(ierr);
  ierr = CAF_Win_unlock(image_index - 1, win); chk_err(ierr);

  /* Temporarily commented */
  /* if (value == 0)
   *   goto stat_error; */

  if (stat)
    *stat = ierr;
  else if (ierr == STAT_FAILED_IMAGE)
    terminate_internal(ierr, 0);

  return;

stat_error:
  if (errmsg != NULL)
  {
    memset(errmsg,' ',errmsg_len);
    memcpy(errmsg, msg, MIN(errmsg_len,strlen(msg)));
  }
  if (stat != NULL)
    *stat = 99;
  else
    terminate_internal(99, 1);
#else // MPI_VERSION
#warning Locking for MPI-2 is not implemented
  printf("Locking for MPI-2 is not supported, "
         "please update your MPI implementation\n");
#endif // MPI_VERSION
}

/* Initialize coarray program.  This routine assumes that no other
 * MPI initialization happened before. */

void
PREFIX(init) (int *argc, char ***argv)
{
  int flag;
  if (caf_num_images == 0)
  {
    int ierr = 0, i = 0, j = 0, rc, prov_lev = 0;
    int is_init = 0, prior_thread_level = MPI_THREAD_FUNNELED;
    ierr = MPI_Initialized(&is_init); chk_err(ierr);

    if (is_init)
    {
      ierr = MPI_Query_thread(&prior_thread_level); chk_err(ierr);
    }
#ifdef HELPER
    if (is_init)
    {
        prov_lev = prior_thread_level;
        caf_owns_mpi = false;
    }
    else
    {
        ierr = MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE, &prov_lev);
        chk_err(ierr);
        caf_owns_mpi = true;
    }

    if (caf_this_image == 0 && MPI_THREAD_MULTIPLE != prov_lev)
      caf_runtime_error("MPI_THREAD_MULTIPLE is not supported: %d", prov_lev);
#else
    if (is_init)
      caf_owns_mpi = false;
    else
    {
      ierr = MPI_Init_thread(argc, argv, prior_thread_level, &prov_lev);
      chk_err(ierr);
      caf_owns_mpi = true;
      if (caf_this_image == 0 && MPI_THREAD_FUNNELED > prov_lev)
        caf_runtime_error("MPI_THREAD_FUNNELED is not supported: %d %d", MPI_THREAD_FUNNELED, prov_lev);
    }
#endif
    if (unlikely ((ierr != MPI_SUCCESS)))
      caf_runtime_error("Failure when initializing MPI: %d", ierr);

    /* Duplicate MPI_COMM_WORLD so that no CAF internal functions use it.
     * This is critical for MPI-interoperability. */
    rc = MPI_Comm_dup(MPI_COMM_WORLD, &CAF_COMM_WORLD);
#ifdef WITH_FAILED_IMAGES
    flag = (MPI_SUCCESS == rc);
    rc = MPIX_Comm_agree(MPI_COMM_WORLD, &flag);
    if (rc != MPI_SUCCESS)
    {
      dprint("MPIX_Comm_agree(flag = %d) = %d.\n", flag, rc);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD, 10000);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    ierr = MPI_Comm_size(CAF_COMM_WORLD, &caf_num_images); chk_err(ierr);
    ierr = MPI_Comm_rank(CAF_COMM_WORLD, &caf_this_image); chk_err(ierr);

    ++caf_this_image;
    caf_is_finalized = 0;

    /* BEGIN SYNC IMAGE preparation
     * Prepare memory for syncing images. */
    images_full = (int *) calloc(caf_num_images - 1, sizeof(int));
    for (i = 1, j = 0; i <= caf_num_images; ++i)
    {
      if (i != caf_this_image)
        images_full[j++] = i;
    }

    arrived = calloc(caf_num_images, sizeof(int));
    sync_handles = malloc(caf_num_images * sizeof(MPI_Request));
    /* END SYNC IMAGE preparation. */

    stat_tok = malloc(sizeof(MPI_Win));

    teams_list = (caf_teams_list *)calloc(1, sizeof(caf_teams_list));
    teams_list->team_id = -1;
    MPI_Comm *tmp_comm = (MPI_Comm *)calloc(1, sizeof(MPI_Comm));
    *tmp_comm = CAF_COMM_WORLD;
    teams_list->team = tmp_comm;
    teams_list->prev = NULL;
    used_teams = (caf_used_teams_list *)calloc(1, sizeof(caf_used_teams_list));
    used_teams->team_list_elem = teams_list;
    used_teams->prev = NULL;

#ifdef WITH_FAILED_IMAGES
    MPI_Comm_dup(MPI_COMM_WORLD, &alive_comm);
    /* Handling of failed/stopped images is done by setting an error handler
     * on a asynchronous request to each other image.  For a failing image
     * the request will trigger the call of the error handler thus allowing
     * each other image to handle the failed/stopped image. */
    ierr = MPI_Comm_create_errhandler(failed_stopped_errorhandler_function,
                                      &failed_CAF_COMM_mpi_err_handler);
    chk_err(ierr);
    ierr = MPI_Comm_set_errhandler(CAF_COMM_WORLD,
                                   failed_CAF_COMM_mpi_err_handler);
    chk_err(ierr);
    ierr = MPI_Comm_set_errhandler(alive_comm,
                                   failed_CAF_COMM_mpi_err_handler);
    chk_err(ierr);
    ierr = MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
    chk_err(ierr);

    ierr = MPI_Irecv(&alive_dummy, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
                     alive_comm, &alive_request); chk_err(ierr);

    image_stati = (int *) calloc(caf_num_images, sizeof(int));
#endif

#if MPI_VERSION >= 3
    ierr = MPI_Info_create(&mpi_info_same_size); chk_err(ierr);
    ierr = MPI_Info_set(mpi_info_same_size, "same_size", "true"); chk_err(ierr);

    /* Setting img_status */
    ierr = MPI_Win_create(&img_status, sizeof(int), 1, mpi_info_same_size,
                          CAF_COMM_WORLD, stat_tok); chk_err(ierr);
    CAF_Win_lock_all(*stat_tok);
#else
    ierr = MPI_Win_create(&img_status, sizeof(int), 1, MPI_INFO_NULL,
                          CAF_COMM_WORLD, stat_tok); chk_err(ierr);
#endif // MPI_VERSION

    /* Create the dynamic window to allow images to asyncronously attach
     * memory. */
    ierr = MPI_Win_create_dynamic(MPI_INFO_NULL, CAF_COMM_WORLD,
                                  &global_dynamic_win); chk_err(ierr);
    CAF_Win_lock_all(global_dynamic_win);
  }
}


/* Internal finalize of coarray program. */

void
finalize_internal(int status_code)
{
  int ierr;
  dprint("(status_code = %d)\n", status_code);

#ifdef WITH_FAILED_IMAGES
  no_stopped_images_check_in_errhandler = true;
  ierr = MPI_Win_flush_all(*stat_tok); chk_err(ierr);
#endif
  /* For future security enclose setting img_status in a lock. */
  CAF_Win_lock(MPI_LOCK_EXCLUSIVE, caf_this_image - 1, *stat_tok);
  if (status_code == 0)
  {
    img_status = STAT_STOPPED_IMAGE;
#ifdef WITH_FAILED_IMAGES
    image_stati[caf_this_image - 1] = STAT_STOPPED_IMAGE;
#endif
  }
  else
  {
    img_status = status_code;
#ifdef WITH_FAILED_IMAGES
    image_stati[caf_this_image - 1] = status_code;
#endif
  }
  CAF_Win_unlock(caf_this_image - 1, *stat_tok);

  /* Announce to all other images, that this one has changed its execution
   * status. */
  for (int i = 0; i < caf_num_images - 1; ++i)
  {
    ierr = MPI_Send(&img_status, 1, MPI_INT, images_full[i] - 1,
                    MPI_TAG_CAF_SYNC_IMAGES, CAF_COMM_WORLD); chk_err(ierr);
  }

#ifdef WITH_FAILED_IMAGES
  /* Terminate the async request before revoking the comm, or we will get
   * triggered by the errorhandler, which we don't want here anymore. */
  ierr = MPI_Cancel(&alive_request); chk_err(ierr);

  if (status_code == 0)
  {
    /* In finalization do not report stopped or failed images any more. */
    ierr = MPI_Errhandler_set(CAF_COMM_WORLD, MPI_ERRORS_RETURN); chk_err(ierr);
    ierr = MPI_Errhandler_set(alive_comm, MPI_ERRORS_RETURN); chk_err(ierr);
    /* Only add a conventional barrier to prevent images rom quitting too
     * early, when this images is not failing. */
    dprint("Before MPI_Barrier(CAF_COMM_WORLD)\n");
    ierr = MPI_Barrier(CAF_COMM_WORLD); chk_err(ierr);
    dprint("After MPI_Barrier(CAF_COMM_WORLD) = %d\n", ierr);
  }
  else
    return;
#else
  /* Add a conventional barrier to prevent images from quitting too early. */
  if (status_code == 0)
  {
    ierr = MPI_Barrier(CAF_COMM_WORLD); chk_err(ierr);
  }
  else
    /* Without failed images support, but a given status_code, we need to
     * return to the caller, or we will hang in the following instead of
     * terminating the program. */
    return;
#endif

#ifdef GCC_GE_7
  struct caf_allocated_slave_tokens_t
    *cur_stok = caf_allocated_slave_tokens,
    *prev_stok = NULL;
  CAF_Win_unlock_all(global_dynamic_win);
  while (cur_stok)
  {
    prev_stok = cur_stok->prev;
    ierr = MPI_Win_detach(global_dynamic_win, cur_stok); chk_err(ierr);
    if (cur_stok->token->memptr)
    {
      ierr = MPI_Win_detach(global_dynamic_win, cur_stok->token->memptr);
      chk_err(ierr);
      free(cur_stok->token->memptr);
    }
    free(cur_stok->token);
    free(cur_stok);
    cur_stok = prev_stok;
  }
#else
  CAF_Win_unlock_all(global_dynamic_win);
#endif

  dprint("Freed all slave tokens.\n");
  struct caf_allocated_tokens_t
    *cur_tok = caf_allocated_tokens,
    *prev = caf_allocated_tokens;
  MPI_Win *p;

  while (cur_tok)
  {
    prev = cur_tok->prev;
    p = TOKEN(cur_tok->token);
    if (p != NULL)
      CAF_Win_unlock_all(*p);
#ifdef GCC_GE_7
    /* Unregister the window to the descriptors when freeing the token. */
    dprint("MPI_Win_free(p);\n");
    ierr = MPI_Win_free(p); chk_err(ierr);
    free(cur_tok->token);
#else // GCC_GE_7
    ierr = MPI_Win_free(p); chk_err(ierr);
#endif // GCC_GE_7
    free(cur_tok);
    cur_tok = prev;
  }
#if MPI_VERSION >= 3
  ierr = MPI_Info_free(&mpi_info_same_size); chk_err(ierr);
#endif // MPI_VERSION

  /* Free the global dynamic window. */
  ierr = MPI_Win_free(&global_dynamic_win); chk_err(ierr);
#ifdef WITH_FAILED_IMAGES
  if (status_code == 0)
  {
    dprint("before Win_unlock_all.\n");
    CAF_Win_unlock_all(*stat_tok);
    dprint("before Win_free(stat_tok)\n");
    ierr = MPI_Win_free(stat_tok); chk_err(ierr);
    dprint("before Comm_free(CAF_COMM_WORLD)\n");
    ierr = MPI_Comm_free(&CAF_COMM_WORLD); chk_err(ierr);
    ierr = MPI_Comm_free(&alive_comm); chk_err(ierr);
    dprint("after Comm_free(CAF_COMM_WORLD)\n");
  }

  ierr = MPI_Errhandler_free(&failed_CAF_COMM_mpi_err_handler); chk_err(ierr);

  /* Only call Finalize if CAF runtime Initialized MPI. */
  if (caf_owns_mpi)
  {
    ierr = MPI_Finalize(); chk_err(ierr);
  }
#else
  ierr = MPI_Comm_free(&CAF_COMM_WORLD); chk_err(ierr);

  CAF_Win_unlock_all(*stat_tok);
  ierr = MPI_Win_free(stat_tok); chk_err(ierr);

  /* Only call Finalize if CAF runtime Initialized MPI. */
  if (caf_owns_mpi)
  {
    ierr = MPI_Finalize(); chk_err(ierr);
  }
#endif

  pthread_mutex_lock(&lock_am);
  caf_is_finalized = 1;
  pthread_mutex_unlock(&lock_am);
  free(sync_handles);
  dprint("Finalisation done!!!\n");
}


/* Finalize coarray program. */
void
PREFIX(finalize) (void)
{
  finalize_internal(0);
}

/* TODO: This is interface is violating the F2015 standard, but not the gfortran
 * API. Fix it (the fortran API). */
int
PREFIX(this_image) (int distance __attribute__((unused)))
{
  return caf_this_image;
}

/* TODO: This is interface is violating the F2015 standard, but not the gfortran
 * API. Fix it (the fortran API). */
int
PREFIX(num_images) (int distance __attribute__((unused)),
                    int failed __attribute__((unused)))
{
  return caf_num_images;
}

#ifdef GCC_GE_7
/* Register an object with the coarray library creating a token where
 * necessary/requested.
 * See the ABI-documentation of gfortran for the expected behavior.
 * Contrary to this expected behavior is this routine not registering memory
 * in the descriptor, that is already present.  I.e., when the compiler
 * expects the library to allocate the memory for an object in desc, then
 * its data_ptr is NULL. This is still missing here.  At the moment the
 * compiler also does not make use of it, but it is contrary to the
 * documentation. */
void
PREFIX(register) (size_t size, caf_register_t type, caf_token_t *token,
                  gfc_descriptor_t *desc, int *stat, char *errmsg,
                  charlen_t errmsg_len)
{
  void *mem = NULL;
  size_t actual_size;
  int l_var = 0, *init_array = NULL, ierr;

  if (unlikely(caf_is_finalized))
    goto error;

  /* Start GASNET if not already started. */
  if (caf_num_images == 0)
    PREFIX(init) (NULL, NULL);

  if (type == CAF_REGTYPE_LOCK_STATIC || type == CAF_REGTYPE_LOCK_ALLOC ||
      type == CAF_REGTYPE_CRITICAL || type == CAF_REGTYPE_EVENT_STATIC ||
      type == CAF_REGTYPE_EVENT_ALLOC)
  {
    actual_size = size * sizeof(int);
    l_var = 1;
  }
  else
    actual_size = size;

  switch (type)
  {
    case CAF_REGTYPE_COARRAY_ALLOC_REGISTER_ONLY:
    case CAF_REGTYPE_COARRAY_ALLOC_ALLOCATE_ONLY:
      {
        /* Create or allocate a slave token. */
        mpi_caf_slave_token_t *slave_token;
#ifdef EXTRA_DEBUG_OUTPUT
        MPI_Aint mpi_address;
#endif
        CAF_Win_unlock_all(global_dynamic_win);
        if (type == CAF_REGTYPE_COARRAY_ALLOC_REGISTER_ONLY)
        {
          *token = calloc(1, sizeof(mpi_caf_slave_token_t));
          slave_token = (mpi_caf_slave_token_t *)(*token);
          ierr = MPI_Win_attach(global_dynamic_win, *token,
                                sizeof(mpi_caf_slave_token_t)); chk_err(ierr);
#ifdef EXTRA_DEBUG_OUTPUT
          ierr = MPI_Get_address(*token, &mpi_address); chk_err(ierr);
#endif
          dprint("Attach slave token %p (mpi-address: %zd) to "
                 "global_dynamic_window = %d\n",
                 slave_token, mpi_address, global_dynamic_win);

          /* Register the memory for auto freeing. */
          struct caf_allocated_slave_tokens_t *tmp =
            malloc(sizeof(struct caf_allocated_slave_tokens_t));
          tmp->prev  = caf_allocated_slave_tokens;
          tmp->token = *token;
          caf_allocated_slave_tokens = tmp;
        }
        else // (type == CAF_REGTYPE_COARRAY_ALLOC_ALLOCATE_ONLY)
        {
          int ierr;
          slave_token = (mpi_caf_slave_token_t *)(*token);
          mem = malloc(actual_size);
          slave_token->memptr = mem;
          ierr = MPI_Win_attach(global_dynamic_win, mem, actual_size);
          chk_err(ierr); 
#ifdef EXTRA_DEBUG_OUTPUT
          ierr = MPI_Get_address(mem, &mpi_address); chk_err(ierr);
#endif
          dprint("Attach mem %p (mpi-address: %zd) to global_dynamic_window = "
                 "%d on slave_token %p, size %zd, ierr: %d\n",
                 mem, mpi_address, global_dynamic_win, slave_token,
                 actual_size, ierr);
          if (desc != NULL && GFC_DESCRIPTOR_RANK(desc) != 0)
          {
            slave_token->desc = desc;
#ifdef EXTRA_DEBUG_OUTPUT
            ierr = MPI_Get_address(desc, &mpi_address); chk_err(ierr);
#endif
            dprint("Attached descriptor %p (mpi-address: %zd) to "
                   "global_dynamic_window %d at address %p, ierr = %d.\n",
                   desc, mpi_address, global_dynamic_win, &slave_token->desc,
                   ierr);
          }
        }
        CAF_Win_lock_all(global_dynamic_win);
        dprint("Slave token %p on exit: mpi_caf_slave_token_t { desc: %p }\n",
               slave_token, slave_token->desc);
      }
      break;
    default:
      {
        mpi_caf_token_t *mpi_token;
        MPI_Win *p;

        *token = calloc(1, sizeof(mpi_caf_token_t));
        mpi_token = (mpi_caf_token_t *) (*token);
        p = TOKEN(mpi_token);

#if MPI_VERSION >= 3
        ierr = MPI_Win_allocate(actual_size, 1, MPI_INFO_NULL, CAF_COMM_WORLD,
                                &mem, p); chk_err(ierr);
        CAF_Win_lock_all(*p);
#else // MPI_VERSION
        ierr = MPI_Alloc_mem(actual_size, MPI_INFO_NULL, &mem); chk_err(ierr);
        ierr = MPI_Win_create(mem, actual_size, 1, MPI_INFO_NULL,
                              CAF_COMM_WORLD, p); chk_err(ierr);
#endif // MPI_VERSION

#ifndef GCC_GE_8
        if (GFC_DESCRIPTOR_RANK(desc) != 0)
#endif
          mpi_token->desc = desc;

        if (l_var)
        {
          init_array = (int *)calloc(size, sizeof(int));
          CAF_Win_lock(MPI_LOCK_EXCLUSIVE, caf_this_image - 1, *p);
          ierr = MPI_Put(init_array, size, MPI_INT, caf_this_image - 1, 0, size,
                         MPI_INT, *p); chk_err(ierr);
          CAF_Win_unlock(caf_this_image - 1, *p);
          free(init_array);
        }

        struct caf_allocated_tokens_t *tmp =
            malloc(sizeof(struct caf_allocated_tokens_t));
        tmp->prev  = caf_allocated_tokens;
        tmp->token = *token;
        caf_allocated_tokens = tmp;

        if (stat)
          *stat = 0;

        /* The descriptor will be initialized only after the call to
         * register. */
        mpi_token->memptr = mem;
        dprint("Token %p on exit: mpi_caf_token_t "
               "{ (local_)memptr: %p, memptr_win: %d }\n",
               mpi_token, mpi_token->memptr, mpi_token->memptr_win);
      } // default:
      break;
  } // switch

  desc->base_addr = mem;
  return;

error:
  {
    char msg[80];
    strcpy(msg, "Failed to allocate coarray");
    if (caf_is_finalized)
      strcat(msg, " - there are stopped images");

    if (stat)
    {
      *stat = caf_is_finalized ? STAT_STOPPED_IMAGE : 1;
      if (errmsg_len > 0)
      {
        size_t len = (strlen(msg) > (size_t) errmsg_len) ?
                     (size_t) errmsg_len : strlen (msg);
        memcpy(errmsg, msg, len);
        if ((size_t) errmsg_len > len)
          memset(&errmsg[len], ' ', errmsg_len - len);
      }
    }
    else
      caf_runtime_error(msg);
  }
}
#else // GCC_LT_7
void *
PREFIX(register) (size_t size, caf_register_t type, caf_token_t *token,
                  int *stat, char *errmsg, charlen_t errmsg_len)
{
  void *mem;
  size_t actual_size;
  int l_var = 0, *init_array = NULL, ierr;

  if (unlikely(caf_is_finalized))
    goto error;

  /* Start GASNET if not already started. */
  if (caf_num_images == 0)
#ifdef COMPILER_SUPPORTS_CAF_INTRINSICS
    _gfortran_caf_init(NULL, NULL);
#else
    PREFIX(init) (NULL, NULL);
#endif

  /* Token contains only a list of pointers. */
  *token = malloc(sizeof(MPI_Win));
  MPI_Win *p = *token;

  if (type == CAF_REGTYPE_LOCK_STATIC || type == CAF_REGTYPE_LOCK_ALLOC ||
      type == CAF_REGTYPE_CRITICAL || type == CAF_REGTYPE_EVENT_STATIC ||
      type == CAF_REGTYPE_EVENT_ALLOC)
  {
    actual_size = size * sizeof(int);
    l_var = 1;
  }
  else
    actual_size = size;

#if MPI_VERSION >= 3
  ierr = MPI_Win_allocate(actual_size, 1, mpi_info_same_size, CAF_COMM_WORLD,
                          &mem, p); chk_err(ierr);
  CAF_Win_lock_all(*p);
#else // MPI_VERSION
  ierr = MPI_Alloc_mem(actual_size, MPI_INFO_NULL, &mem); chk_err(ierr);
  ierr = MPI_Win_create(mem, actual_size, 1, MPI_INFO_NULL,
                        CAF_COMM_WORLD, p); chk_err(ierr);
#endif // MPI_VERSION

  if (l_var)
  {
    init_array = (int *)calloc(size, sizeof(int));
    CAF_Win_lock(MPI_LOCK_EXCLUSIVE, caf_this_image - 1, *p);
    ierr = MPI_Put(init_array, size, MPI_INT, caf_this_image - 1, 0, size,
                   MPI_INT, *p); chk_err(ierr);
    CAF_Win_unlock(caf_this_image - 1, *p);
    free(init_array);
  }

  PREFIX(sync_all) (NULL, NULL, 0);

  struct caf_allocated_tokens_t *tmp =
         malloc(sizeof(struct caf_allocated_tokens_t));
  tmp->prev  = caf_allocated_tokens;
  tmp->token = *token;
  caf_allocated_tokens = tmp;

  if (stat)
    *stat = 0;
  return mem;

error:
  {
    char msg[80];
    strcpy(msg, "Failed to allocate coarray");
    if (caf_is_finalized)
      strcat(msg, " - there are stopped images");

    if (stat)
    {
      *stat = caf_is_finalized ? STAT_STOPPED_IMAGE : 1;
      if (errmsg_len > 0)
      {
        size_t len = (strlen(msg) > (size_t) errmsg_len) ?
                     (size_t) errmsg_len : strlen (msg);
        memcpy(errmsg, msg, len);
        if (errmsg_len > len)
          memset(&errmsg[len], ' ', errmsg_len - len);
      }
    }
    else
      caf_runtime_error(msg);
  }
  return NULL;
}
#endif // GCC_GE_7


#ifdef GCC_GE_7
void
PREFIX(deregister) (caf_token_t *token, int type, int *stat, char *errmsg,
                    charlen_t errmsg_len)
#else
void
PREFIX(deregister) (caf_token_t *token, int *stat, char *errmsg,
                    charlen_t errmsg_len)
#endif
{
  dprint("deregister(%p)\n", *token);
  int ierr;

  if (unlikely(caf_is_finalized))
  {
    const char msg[] =
      "Failed to deallocate coarray - there are stopped images";
    if (stat)
    {
      *stat = STAT_STOPPED_IMAGE;

      if (errmsg_len > 0)
      {
        size_t len = (sizeof(msg) - 1 > (size_t) errmsg_len) ?
          (size_t) errmsg_len : sizeof (msg) - 1;
        memcpy(errmsg, msg, len);
        if (errmsg_len > len)
          memset(&errmsg[len], ' ', errmsg_len - len);
      }
      return;
    }
    caf_runtime_error(msg);
  }

  if (stat)
    *stat = 0;

#ifdef GCC_GE_7
  if (type != CAF_DEREGTYPE_COARRAY_DEALLOCATE_ONLY)
  {
    /* Sync all images only, when deregistering the token. Just freeing the
     * memory needs no sync. */
#ifdef WITH_FAILED_IMAGES
    ierr = MPI_Barrier(CAF_COMM_WORLD); chk_err(ierr);
#else
    PREFIX(sync_all) (NULL, NULL, 0);
#endif
  }
#endif // GCC_GE_7
  {
    struct caf_allocated_tokens_t
      *cur = caf_allocated_tokens,
      *next = caf_allocated_tokens,
      *prev;
    MPI_Win *p;

    while (cur)
    {
      prev = cur->prev;

      if (cur->token == *token)
      {
        p = TOKEN(*token);
#ifdef GCC_GE_7
        dprint("Found regular token %p for memptr_win: %d.\n",
               *token, ((mpi_caf_token_t *)*token)->memptr_win);
#endif
        CAF_Win_unlock_all(*p);
        ierr = MPI_Win_free(p); chk_err(ierr);

        next->prev = prev ? prev->prev:  NULL;

        if (cur == caf_allocated_tokens)
          caf_allocated_tokens = prev;

        free(cur);
        free(*token);
        return;
      }
      next = cur;
      cur = prev;
    }
  }

#ifdef GCC_GE_7
  /* Feel through: Has to be a component token. */
  {
    struct caf_allocated_slave_tokens_t
      *cur_stok = caf_allocated_slave_tokens,
      *next_stok = caf_allocated_slave_tokens,
      *prev_stok;

    while (cur_stok)
    {
      prev_stok = cur_stok->prev;

      if (cur_stok->token == *token)
      {
        dprint("Found sub token %p.\n", *token);

        mpi_caf_slave_token_t *slave_token = *(mpi_caf_slave_token_t **)token;
        CAF_Win_unlock_all(global_dynamic_win);

        if (slave_token->memptr)
        {
          ierr = MPI_Win_detach(global_dynamic_win, slave_token->memptr);
          chk_err(ierr);
          free(slave_token->memptr);
          slave_token->memptr = NULL;
          if (type == CAF_DEREGTYPE_COARRAY_DEALLOCATE_ONLY)
          {
            CAF_Win_lock_all(global_dynamic_win);
            return; // All done.
          }
        }
        ierr = MPI_Win_detach(global_dynamic_win, slave_token); chk_err(ierr);
        CAF_Win_lock_all(global_dynamic_win);

        next_stok->prev = prev_stok ? prev_stok->prev: NULL;

        if (cur_stok == caf_allocated_slave_tokens)
          caf_allocated_slave_tokens = prev_stok;

        free(cur_stok);
        free(*token);
        return;
      }

      next_stok = cur_stok;
      cur_stok = prev_stok;
    }
  }
#endif // GCC_GE_7
#ifdef EXTRA_DEBUG_OUTPUT
  fprintf(stderr,
          "Fortran runtime warning on image %d: "
          "Could not find token to free %p", caf_this_image, *token);
#endif
}

void
PREFIX(sync_memory) (int *stat __attribute__((unused)),
                     char *errmsg __attribute__((unused)),
                     charlen_t errmsg_len __attribute__((unused)))
{
#if defined(NONBLOCKING_PUT) && !defined(CAF_MPI_LOCK_UNLOCK)
  explicit_flush();
#endif
}


void
PREFIX(sync_all) (int *stat, char *errmsg, charlen_t errmsg_len)
{
  int err = 0, ierr;

  dprint("Entering sync all.\n");
  if (unlikely(caf_is_finalized))
  {
    err = STAT_STOPPED_IMAGE;
  }
  else
  {
#if defined(NONBLOCKING_PUT) && !defined(CAF_MPI_LOCK_UNLOCK)
    explicit_flush();
#endif

#ifdef WITH_FAILED_IMAGES
    ierr = MPI_Barrier(alive_comm); chk_err(ierr);
#else
    ierr = MPI_Barrier(CAF_COMM_WORLD); chk_err(ierr);
#endif
    dprint("MPI_Barrier = %d.\n", err);
    if (ierr == STAT_FAILED_IMAGE)
      err = STAT_FAILED_IMAGE;
    else if (ierr != 0)
      MPI_Error_class(ierr, &err);
  }

  if (stat != NULL)
    *stat = err;
#ifdef WITH_FAILED_IMAGES
  else if (err == STAT_FAILED_IMAGE)
    /* F2015 requests stat to be set for FAILED IMAGES, else error out. */
    terminate_internal(err, 0);
#endif

  if (err != 0 && err != STAT_FAILED_IMAGE)
  {
    char msg[80];
    strcpy(msg, "SYNC ALL failed");
    if (caf_is_finalized)
      strcat(msg, " - there are stopped images");

    if (errmsg_len > 0)
    {
      size_t len = (strlen(msg) > (size_t) errmsg_len) ?
                   (size_t) errmsg_len : strlen (msg);
      memcpy(errmsg, msg, len);
      if (errmsg_len > len)
        memset(&errmsg[len], ' ', errmsg_len - len);
    }
    else if (stat == NULL)
      caf_runtime_error(msg);
  }
  dprint("Leaving sync all.\n");
}

/* Convert kind 4 characters into kind 1 one.
 * Copied from the gcc:libgfortran/caf/single.c. */
static void
assign_char4_from_char1(size_t dst_size, size_t src_size, uint32_t *dst,
                        unsigned char *src)
{
  size_t i, n;
  n = (dst_size > src_size) ? src_size : dst_size;
  for (i = 0; i < n; ++i)
  {
    dst[i] = (int32_t) src[i];
  }
  for (; i < dst_size; ++i)
  {
    dst[i] = (int32_t) ' ';
  }
}


/* Convert kind 1 characters into kind 4 one.
 * Copied from the gcc:libgfortran/caf/single.c. */
static void
assign_char1_from_char4(size_t dst_size, size_t src_size, unsigned char *dst,
                        uint32_t *src)
{
  size_t i, n;
  n = (dst_size > src_size) ? src_size : dst_size;
  for (i = 0; i < n; ++i)
  {
    dst[i] = src[i] > UINT8_MAX ? (unsigned char) '?' : (unsigned char) src[i];
  }
  if (dst_size > n)
    memset(&dst[n], ' ', dst_size - n);
}

/* Convert convertable types.
 * Copied from the gcc:libgfortran/caf/single.c. Can't say much about it. */
static void
convert_type(void *dst, int dst_type, int dst_kind, void *src, int src_type,
             int src_kind, int *stat)
{
#ifdef HAVE_GFC_INTEGER_16
  typedef __int128 int128t;
#else
  typedef int64_t int128t;
#endif

#if defined(GFC_REAL_16_IS_LONG_DOUBLE)
  typedef long double real128t;
  typedef _Complex long double complex128t;
#elif defined(HAVE_GFC_REAL_16)
  typedef _Complex float __attribute__((mode(TC))) __complex128;
  typedef __float128 real128t;
  typedef __complex128 complex128t;
#elif defined(HAVE_GFC_REAL_10)
  typedef long double real128t;
  typedef long double complex128t;
#else
  typedef double real128t;
  typedef _Complex double complex128t;
#endif

  int128t int_val = 0;
  real128t real_val = 0;
  complex128t cmpx_val = 0;

  switch (src_type)
  {
    case BT_INTEGER:
      if (src_kind == 1)
        int_val = *(int8_t*) src;
      else if (src_kind == 2)
        int_val = *(int16_t*) src;
      else if (src_kind == 4)
        int_val = *(int32_t*) src;
      else if (src_kind == 8)
        int_val = *(int64_t*) src;
#ifdef HAVE_GFC_INTEGER_16
      else if (src_kind == 16)
        int_val = *(int128t*) src;
#endif
      else
        goto error;
      break;
    case BT_REAL:
      if (src_kind == 4)
        real_val = *(float*) src;
      else if (src_kind == 8)
        real_val = *(double*) src;
#ifdef HAVE_GFC_REAL_10
      else if (src_kind == 10)
        real_val = *(long double*) src;
#endif
#ifdef HAVE_GFC_REAL_16
      else if (src_kind == 16)
        real_val = *(real128t*) src;
#endif
      else
        goto error;
      break;
    case BT_COMPLEX:
      if (src_kind == 4)
        cmpx_val = *(_Complex float*) src;
      else if (src_kind == 8)
        cmpx_val = *(_Complex double*) src;
#ifdef HAVE_GFC_REAL_10
      else if (src_kind == 10)
        cmpx_val = *(_Complex long double*) src;
#endif
#ifdef HAVE_GFC_REAL_16
      else if (src_kind == 16)
        cmpx_val = *(complex128t*) src;
#endif
      else
        goto error;
      break;
    default:
      goto error;
  }

  switch (dst_type)
  {
    case BT_INTEGER:
      if (src_type == BT_INTEGER)
      {
        if (dst_kind == 1)
          *(int8_t*) dst = (int8_t) int_val;
        else if (dst_kind == 2)
          *(int16_t*) dst = (int16_t) int_val;
        else if (dst_kind == 4)
          *(int32_t*) dst = (int32_t) int_val;
        else if (dst_kind == 8)
          *(int64_t*) dst = (int64_t) int_val;
#ifdef HAVE_GFC_INTEGER_16
        else if (dst_kind == 16)
          *(int128t*) dst = (int128t) int_val;
#endif
        else
          goto error;
      }
      else if (src_type == BT_REAL)
      {
        if (dst_kind == 1)
          *(int8_t*) dst = (int8_t) real_val;
        else if (dst_kind == 2)
          *(int16_t*) dst = (int16_t) real_val;
        else if (dst_kind == 4)
          *(int32_t*) dst = (int32_t) real_val;
        else if (dst_kind == 8)
          *(int64_t*) dst = (int64_t) real_val;
#ifdef HAVE_GFC_INTEGER_16
        else if (dst_kind == 16)
          *(int128t*) dst = (int128t) real_val;
#endif
        else
          goto error;
      }
      else if (src_type == BT_COMPLEX)
      {
        if (dst_kind == 1)
          *(int8_t*) dst = (int8_t) cmpx_val;
        else if (dst_kind == 2)
          *(int16_t*) dst = (int16_t) cmpx_val;
        else if (dst_kind == 4)
          *(int32_t*) dst = (int32_t) cmpx_val;
        else if (dst_kind == 8)
          *(int64_t*) dst = (int64_t) cmpx_val;
#ifdef HAVE_GFC_INTEGER_16
        else if (dst_kind == 16)
          *(int128t*) dst = (int128t) cmpx_val;
#endif
        else
          goto error;
      }
      else
        goto error;
      return;
    case BT_REAL:
      if (src_type == BT_INTEGER)
      {
        if (dst_kind == 4)
          *(float*) dst = (float) int_val;
        else if (dst_kind == 8)
          *(double*) dst = (double) int_val;
#ifdef HAVE_GFC_REAL_10
        else if (dst_kind == 10)
          *(long double*) dst = (long double) int_val;
#endif
#ifdef HAVE_GFC_REAL_16
        else if (dst_kind == 16)
          *(real128t*) dst = (real128t) int_val;
#endif
        else
          goto error;
      }
      else if (src_type == BT_REAL)
      {
        if (dst_kind == 4)
          *(float*) dst = (float) real_val;
        else if (dst_kind == 8)
          *(double*) dst = (double) real_val;
#ifdef HAVE_GFC_REAL_10
        else if (dst_kind == 10)
          *(long double*) dst = (long double) real_val;
#endif
#ifdef HAVE_GFC_REAL_16
        else if (dst_kind == 16)
          *(real128t*) dst = (real128t) real_val;
#endif
        else
          goto error;
      }
      else if (src_type == BT_COMPLEX)
      {
        if (dst_kind == 4)
          *(float*) dst = (float) cmpx_val;
        else if (dst_kind == 8)
          *(double*) dst = (double) cmpx_val;
#ifdef HAVE_GFC_REAL_10
        else if (dst_kind == 10)
          *(long double*) dst = (long double) cmpx_val;
#endif
#ifdef HAVE_GFC_REAL_16
        else if (dst_kind == 16)
          *(real128t*) dst = (real128t) cmpx_val;
#endif
        else
          goto error;
      }
      return;
    case BT_COMPLEX:
      if (src_type == BT_INTEGER)
      {
        if (dst_kind == 4)
          *(_Complex float*) dst = (_Complex float) int_val;
        else if (dst_kind == 8)
          *(_Complex double*) dst = (_Complex double) int_val;
#ifdef HAVE_GFC_REAL_10
        else if (dst_kind == 10)
          *(_Complex long double*) dst = (_Complex long double) int_val;
#endif
#ifdef HAVE_GFC_REAL_16
        else if (dst_kind == 16)
          *(complex128t*) dst = (complex128t) int_val;
#endif
        else
          goto error;
      }
      else if (src_type == BT_REAL)
      {
        if (dst_kind == 4)
          *(_Complex float*) dst = (_Complex float) real_val;
        else if (dst_kind == 8)
          *(_Complex double*) dst = (_Complex double) real_val;
#ifdef HAVE_GFC_REAL_10
        else if (dst_kind == 10)
          *(_Complex long double*) dst = (_Complex long double) real_val;
#endif
#ifdef HAVE_GFC_REAL_16
        else if (dst_kind == 16)
          *(complex128t*) dst = (complex128t) real_val;
#endif
        else
          goto error;
      }
      else if (src_type == BT_COMPLEX)
      {
        if (dst_kind == 4)
          *(_Complex float*) dst = (_Complex float) cmpx_val;
        else if (dst_kind == 8)
          *(_Complex double*) dst = (_Complex double) cmpx_val;
#ifdef HAVE_GFC_REAL_10
        else if (dst_kind == 10)
          *(_Complex long double*) dst = (_Complex long double) cmpx_val;
#endif
#ifdef HAVE_GFC_REAL_16
        else if (dst_kind == 16)
          *(complex128t*) dst = (complex128t) cmpx_val;
#endif
        else
          goto error;
      }
      else
        goto error;
      return;
    default:
      goto error;
  }

error:
  fprintf(stderr,
          "libcaf_mpi RUNTIME ERROR: Cannot convert type %d kind %d "
          "to type %d kind %d\n", src_type, src_kind, dst_type, dst_kind);
  if (stat)
    *stat = 1;
  else
    abort();
}

static void
convert_with_strides(void *dst, int dst_type, int dst_kind,
                     ptrdiff_t byte_dst_stride,
                     void *src, int src_type, int src_kind,
                     ptrdiff_t byte_src_stride, size_t num, int *stat)
{
  /* Compute the step from one item to convert to the next in bytes. The stride
   * is expected to be the one or similar to the array.stride, i.e. *_stride is
   * expected to be >= 1 to progress from one item to the next. */
  for (size_t i = 0; i < num; ++i)
  {
    convert_type(dst, dst_type, dst_kind, src, src_type, src_kind, stat);
    dst += byte_dst_stride;
    src += byte_src_stride;
  }
}

static void
copy_char_to_self(void *src, int src_type, int src_size, int src_kind,
                  void *dst, int dst_type, int dst_size, int dst_kind,
                  size_t size, bool src_is_scalar)
{
#ifdef GFC_CAF_CHECK
  if (dst_type != BT_CHARACTER || src_type != BT_CHARACTER)
    caf_runtime_error("internal error: copy_char_to_self() "
                      "for non-char types called.");
#endif
  const size_t 
      dst_len = dst_size / dst_kind,
      src_len = src_size / src_kind;
  const size_t min_len = (src_len < dst_len) ? src_len : dst_len;
  /* The address of dest passed by the compiler points on the right  memory
   * location. No offset summation is needed. */
  if (dst_kind == src_kind)
  {
    for (size_t c = 0; c < size; ++c)
    {
      memmove(dst, src, min_len * dst_kind);
      /* Fill dest when source is too short. */
      if (dst_len > src_len)
      {
        int32_t * dest_addr = (int32_t *)(dst + dst_kind * src_len);
        const size_t pad_num = dst_len - src_len;
        if (dst_kind == 1)
          memset(dest_addr, ' ', pad_num);
        else if (dst_kind == 4)
        {
          const void * end_addr = &(dest_addr[pad_num]);
          while (dest_addr != end_addr)
            *(dest_addr++) = (int32_t)' ';
        }
        else
          caf_runtime_error(unreachable);
      }
      dst = (void *)((ptrdiff_t)(dst) + dst_size);
      if (!src_is_scalar)
        src = (void *)((ptrdiff_t)(src) + src_size);
    }
  }
  else
  {
    /* Assign using kind-conversion. */
    if (dst_kind == 1 && src_kind == 4)
      for (size_t c = 0; c < size; ++c)
      {
        assign_char1_from_char4(dst_len, src_len, dst, src);
        dst = (void *)((ptrdiff_t)(dst) + dst_size);
        if (!src_is_scalar)
          src = (void *)((ptrdiff_t)(src) + src_size);
      }
    else if (dst_kind == 4 && src_kind == 1)
      for (size_t c = 0; c < size; ++c)
      {
        assign_char4_from_char1(dst_len, src_len, dst, src);
        dst = (void *)((ptrdiff_t)(dst) + dst_size);
        if (!src_is_scalar)
          src = (void *)((ptrdiff_t)(src) + src_size);
      }
    else
      caf_runtime_error("_caf_send(): Unsupported char kinds in same image "
                        "assignment (kind(lhs)= %d, kind(rhs) = %d)",
                        dst_kind, src_kind);
  }
}

static void
copy_to_self(gfc_descriptor_t *src, int src_kind,
              gfc_descriptor_t *dest, int dst_kind, size_t size, int *stat)
{
#ifdef GFC_CAF_CHECK
  if (GFC_DESCRIPTOR_TYPE(dest) == BT_CHARACTER
      || GFC_DESCRIPTOR_TYPE(src) == BT_CHARACTER)
    caf_runtime_error("internal error: copy_to_self() for char types called.");
#endif
  /* The address of dest passed by the compiler points on the right
   * memory location. No offset summation is needed. */
  if (dst_kind == src_kind)
    memmove(dest->base_addr, src->base_addr, size * GFC_DESCRIPTOR_SIZE(dest));
  else
    /* When the rank is 0 then a scalar is copied to a vector and the stride
     * is zero. */
    convert_with_strides(dest->base_addr, GFC_DESCRIPTOR_TYPE(dest), dst_kind,
                         GFC_DESCRIPTOR_SIZE(dest), src->base_addr,
                         GFC_DESCRIPTOR_TYPE(src), src_kind,
                         (GFC_DESCRIPTOR_RANK(src) > 0)
                         ? GFC_DESCRIPTOR_SIZE(src) : 0, size, stat);
}

/* token: The token of the array to be written to. 
 * offset: Difference between the coarray base address and the actual data,
 * used for caf(3)[2] = 8 or caf[4]%a(4)%b = 7. 
 * image_index: Index of the coarray (typically remote,
 * though it can also be on this_image). 
 * data: Pointer to the to-be-transferred data. 
 * size: The number of bytes to be transferred. 
 * asynchronous: Return before the data transfer has been complete */

void selectType(int size, MPI_Datatype *dt)
{
  int t_s;

#define SELTYPE(type) MPI_Type_size(type, &t_s);  \
if (t_s == size)                                  \
{                                                 \
  *dt = type;                                     \
  return;                                         \
}

  SELTYPE(MPI_BYTE)
  SELTYPE(MPI_SHORT)
  SELTYPE(MPI_INT)
  SELTYPE(MPI_DOUBLE)
  SELTYPE(MPI_COMPLEX)
  SELTYPE(MPI_DOUBLE_COMPLEX)

#undef SELTYPE
}

void
PREFIX(sendget) (caf_token_t token_s, size_t offset_s, int image_index_s,
                 gfc_descriptor_t *dest, caf_vector_t *dst_vector,
                 caf_token_t token_g, size_t offset_g, int image_index_g,
                 gfc_descriptor_t *src , caf_vector_t *src_vector,
                 int dst_kind, int src_kind, bool mrt, int *pstat)
{
  int j, ierr = 0;
  size_t i, size;
  ptrdiff_t dimextent;
  const int
    src_rank = GFC_DESCRIPTOR_RANK(src),
    dst_rank = GFC_DESCRIPTOR_RANK(dest);
  const size_t
    src_size = GFC_DESCRIPTOR_SIZE(src),
    dst_size = GFC_DESCRIPTOR_SIZE(dest);
  const int
    src_type = GFC_DESCRIPTOR_TYPE(src),
    dst_type = GFC_DESCRIPTOR_TYPE(dest);
  const bool
    src_contiguous = PREFIX(is_contiguous) (src),
    dst_contiguous = PREFIX(is_contiguous) (dest);
  const bool
    src_same_image = caf_this_image == image_index_g,
    dst_same_image = caf_this_image == image_index_s,
    same_type_and_kind = dst_type == src_type && dst_kind == src_kind;

  MPI_Win *p = TOKEN(token_g);
  ptrdiff_t src_offset = 0, dst_offset = 0;
  void *pad_str = NULL;
  bool free_pad_str = false;
  void *src_t_buff = NULL, *dst_t_buff = NULL;
  bool free_src_t_buff = false, free_dst_t_buff = false;
  const bool
    dest_char_array_is_longer = dst_type == BT_CHARACTER && dst_size > src_size;
  int
    src_remote_image = image_index_g - 1,
    dst_remote_image = image_index_s - 1;

  if (!src_same_image)
  {
    MPI_Group current_team_group, win_group;
    ierr = MPI_Comm_group(CAF_COMM_WORLD, &current_team_group); chk_err(ierr);
    ierr = MPI_Win_get_group(*p, &win_group); chk_err(ierr);
    ierr = MPI_Group_translate_ranks(current_team_group, 1,
                                     (int[]){src_remote_image}, win_group,
                                     &src_remote_image); chk_err(ierr);
    ierr = MPI_Group_free(&current_team_group); chk_err(ierr);
    ierr = MPI_Group_free(&win_group); chk_err(ierr);
  }
  if (!dst_same_image)
  {
    MPI_Group current_team_group, win_group;
    ierr = MPI_Comm_group(CAF_COMM_WORLD, &current_team_group); chk_err(ierr);
    ierr = MPI_Win_get_group(*p, &win_group); chk_err(ierr);
    ierr = MPI_Group_translate_ranks(current_team_group, 1,
                                     (int[]){dst_remote_image}, win_group,
                                     &dst_remote_image); chk_err(ierr);
    ierr = MPI_Group_free(&current_team_group); chk_err(ierr);
    ierr = MPI_Group_free(&win_group); chk_err(ierr);
  }

  /* Ensure stat is always set. */
#ifdef GCC_GE_7
  int * stat = pstat;
  if (stat)
    *stat = 0;
#else
  /* Gcc prior to 7.0 does not have stat here. */
  int * stat = NULL;
#endif

  size = 1;
  for (j = 0; j < dst_rank; ++j)
  {
    dimextent = dest->dim[j]._ubound - dest->dim[j].lower_bound + 1;
    if (dimextent < 0)
      dimextent = 0;
    size *= dimextent;
  }

  if (size == 0)
    return;

  dprint("src_vector = %p, dst_vector = %p, src_image_index = %d, "
         "dst_image_index = %d, offset_src = %zd, offset_dst = %zd.\n",
         src_vector, dst_vector, image_index_g, image_index_s,
         offset_g, offset_s);
  check_image_health(image_index_g, stat);
  check_image_health(image_index_s, stat);

  /* For char arrays: create the padding array, when dst is longer than src. */
  if (dest_char_array_is_longer)
  {
    const size_t pad_num = (dst_size / dst_kind) - (src_size / src_kind);
    const size_t pad_sz = pad_num * dst_kind;
    /* For big arrays alloca() may not be able to get the memory on the stack.
     * Use a regular malloc then. */
    if ((free_pad_str = ((pad_str = alloca(pad_sz)) == NULL)))
    {
      pad_str = malloc(pad_sz);
      if (src_t_buff == NULL)
        caf_runtime_error("Unable to allocate memory "
                          "for internal buffer in sendget().");
    }
    if (dst_kind == 1)
    {
      memset(pad_str, ' ', pad_num);
    }
    else /* dst_kind == 4. */
    {
      for (int32_t *it = (int32_t *) pad_str,
           *itEnd = ((int32_t *) pad_str) + pad_num; it < itEnd; ++it)
      {
        *it = (int32_t) ' ';
      }
    }
  }

  if (src_contiguous && src_vector == NULL)
  {
    if (src_same_image)
    {
      dprint("in caf_this == image_index, size = %zd, "
             "dst_kind = %d, src_kind = %d, dst_size = %zd, src_size = %zd\n",
             size, dst_kind, src_kind, dst_size, src_size);
      src_t_buff = src->base_addr;
      if (same_type_and_kind && !dest_char_array_is_longer)
      {
        dst_t_buff = src_t_buff;
      }
      else
      {
        dprint("allocating %zd bytes for dst_t_buff.\n", dst_size * size);
        if ((free_dst_t_buff = ((dst_t_buff = alloca(dst_size * size)) == NULL)))
        {
          dst_t_buff = malloc(dst_size * size);
          if (dst_t_buff == NULL)
            caf_runtime_error("Unable to allocate memory "
                              "for internal buffer in sendget().");
        }
        if (dst_type == BT_CHARACTER)
        {
          /* The size is encoded in the descriptor's type for char arrays. */
          copy_char_to_self(src->base_addr, src_type, src_size, src_kind,
                            dst_t_buff, dst_type, dst_size, dst_kind,
                            size, src_rank == 0);
        }
        else
        {
          convert_with_strides(dst_t_buff, dst_type, dst_kind, dst_size,
                               src->base_addr, src_type, src_kind,
                               (src_rank > 0) ? src_size : 0, size, stat);
        }
      }
    }
    else
    {
      /* When replication is needed, only access the scalar on the remote. */
      const size_t src_real_size = src_rank > 0 ?
        (src_size * size) : src_size;
      if ((free_dst_t_buff = ((dst_t_buff = alloca(dst_size * size)) == NULL)))
      {
        dst_t_buff = malloc(dst_size * size);
        if (dst_t_buff == NULL)
          caf_runtime_error("Unable to allocate memory "
                            "for internal buffer in sendget().");
      }

      if (dst_kind != src_kind || src_rank == 0 || dest_char_array_is_longer)
      {
        if ((free_src_t_buff = ((src_t_buff = alloca(src_size * size)) == NULL)))
        {
          src_t_buff = malloc(src_size * size);
          if (src_t_buff == NULL)
            caf_runtime_error("Unable to allocate memory "
                              "for internal buffer in sendget().");
        }
      }
      else
        src_t_buff = dst_t_buff;

      CAF_Win_lock(MPI_LOCK_SHARED, src_remote_image, *p);
      if ((same_type_and_kind && dst_rank == src_rank)
          || dst_type == BT_CHARACTER)
      {
        if (!dest_char_array_is_longer
            && (dst_kind == src_kind || dst_type != BT_CHARACTER))
        {
          const size_t trans_size =
            ((dst_size > src_size) ? src_size : dst_size) * size;
          ierr = MPI_Get(dst_t_buff, trans_size, MPI_BYTE, src_remote_image,
                         offset_g, trans_size, MPI_BYTE, *p); chk_err(ierr);
        }
        else
        {
          ierr = MPI_Get(src_t_buff, src_real_size, MPI_BYTE, src_remote_image,
                         offset_g, src_real_size, MPI_BYTE, *p); chk_err(ierr);
          dprint("copy_char_to_self(src_size = %zd, src_kind = %d, "
                 "dst_size = %zd, dst_kind = %d, size = %zd)\n",
                 src_size, src_kind, dst_size, dst_kind, size);
          copy_char_to_self(src_t_buff, src_type, src_size, src_kind,
                            dst_t_buff, dst_type, dst_size, dst_kind,
                            size, src_rank == 0);
          dprint("|%s|\n", (char *)dst_t_buff);
        }
      }
      else
      {
        ierr = MPI_Get(src_t_buff, src_real_size, MPI_BYTE, src_remote_image,
                       offset_g, src_real_size, MPI_BYTE, *p); chk_err(ierr);
        convert_with_strides(dst_t_buff, dst_type, dst_kind, dst_size,
                             src_t_buff, src_type, src_kind,
                             (src_rank > 0) ? src_size: 0, size, stat);
      }
      CAF_Win_unlock(src_remote_image, *p);
    }
  }
#ifdef STRIDED
  else if (!src_same_image && same_type_and_kind && dst_type != BT_CHARACTER)
  {
    /* For strided copy, no type and kind conversion, copy to self or
     * character arrays are supported. */
    MPI_Datatype dt_s, dt_d, base_type_src, base_type_dst;
    int *arr_bl;
    int *arr_dsp_s;

    if ((free_dst_t_buff = ((dst_t_buff = alloca(dst_size * size)) == NULL)))
    {
      dst_t_buff = malloc(dst_size * size);
      if (dst_t_buff == NULL)
        caf_runtime_error("Unable to allocate memory "
                          "for internal buffer in sendget().");
    }

    selectType(src_size, &base_type_src);
    selectType(dst_size, &base_type_dst);

    if (src_rank == 1)
    {
      if (src_vector == NULL)
      {
        ierr = MPI_Type_vector(size, 1, src->dim[0]._stride, base_type_src,
                               &dt_s); chk_err(ierr);
      }
      else
      {
        arr_bl = calloc(size, sizeof(int));
        arr_dsp_s = calloc(size, sizeof(int));

        dprint("Setting up strided vector index.\n");
#define KINDCASE(kind, type)                                            \
case kind:                                                              \
  for (i = 0; i < size; ++i)                                            \
  {                                                                     \
    arr_dsp_s[i] = ((ptrdiff_t)                                         \
      ((type *) src_vector->u.v.vector)[i] - src->dim[0].lower_bound);  \
    arr_bl[i] = 1;                                                      \
  }                                                                     \
  break

        switch (src_vector->u.v.kind)
        {
          KINDCASE(1, int8_t);
          KINDCASE(2, int16_t);
          KINDCASE(4, int32_t);
          KINDCASE(8, int64_t);
#ifdef HAVE_GFC_INTEGER_16
          KINDCASE(16, __int128);
#endif
        default:
          caf_runtime_error(unreachable);
          return;
        }
#undef KINDCASE
        ierr = MPI_Type_indexed(size, arr_bl, arr_dsp_s, base_type_src, &dt_s);
        chk_err(ierr);
        free(arr_bl);
        free(arr_dsp_s);
      }
      ierr = MPI_Type_vector(size, 1, 1, base_type_dst, &dt_d); chk_err(ierr);
    }
    else
    {
      arr_bl = calloc(size, sizeof(int));
      arr_dsp_s = calloc(size, sizeof(int));

      for (i = 0; i < size; ++i)
      {
        arr_bl[i] = 1;
      }

      for (i = 0; i < size; ++i)
      {
        ptrdiff_t array_offset_sr = 0, extent = 1, tot_ext = 1;
        if (src_vector == NULL)
        {
          for (j = 0; j < src_rank - 1; ++j)
          {
            extent = src->dim[j]._ubound - src->dim[j].lower_bound + 1;
            array_offset_sr += ((i / tot_ext) % extent) * src->dim[j]._stride;
            tot_ext *= extent;
          }

          array_offset_sr += (i / tot_ext) * src->dim[src_rank - 1]._stride;
        }
        else
        {
#define KINDCASE(kind, type)                                          \
case kind:                                                            \
  array_offset_sr = ((ptrdiff_t)                                      \
    ((type *) src_vector->u.v.vector)[i] - src->dim[0].lower_bound);  \
  break
          switch (src_vector->u.v.kind)
          {
            KINDCASE(1, int8_t);
            KINDCASE(2, int16_t);
            KINDCASE(4, int32_t);
            KINDCASE(8, int64_t);
#ifdef HAVE_GFC_INTEGER_16
            KINDCASE(16, __int128);
#endif
            default:
              caf_runtime_error(unreachable);
              return;
          }
#undef KINDCASE
        }
        arr_dsp_s[i] = array_offset_sr;
      }

      ierr = MPI_Type_indexed(size, arr_bl, arr_dsp_s, base_type_src, &dt_s);
      chk_err(ierr);
      ierr = MPI_Type_vector(size, 1, 1, base_type_dst, &dt_d); chk_err(ierr);

      free(arr_bl);
      free(arr_dsp_s);
    }

    ierr = MPI_Type_commit(&dt_s); chk_err(ierr);
    ierr = MPI_Type_commit(&dt_d); chk_err(ierr);

    CAF_Win_lock(MPI_LOCK_SHARED, src_remote_image, *p);
    ierr = MPI_Get(dst_t_buff, 1, dt_d, src_remote_image, offset_g, 1,
                   dt_s, *p); chk_err(ierr);
    CAF_Win_unlock(src_remote_image, *p);

#ifdef WITH_FAILED_IMAGES
    check_image_health(image_index_g, stat);

    if (!stat && ierr == STAT_FAILED_IMAGE)
      terminate_internal(ierr, 1);

    if (stat)
      *stat = ierr;
#else
    if (ierr != 0)
    {
      terminate_internal(ierr, 1);
      return;
    }
#endif
    ierr = MPI_Type_free(&dt_s); chk_err(ierr);
    ierr = MPI_Type_free(&dt_d); chk_err(ierr);
  }
#endif // STRIDED
  else
  {
    if ((free_dst_t_buff = ((dst_t_buff = alloca(dst_size * size)) == NULL)))
    {
      dst_t_buff = malloc(dst_size * size);
      if (dst_t_buff == NULL)
        caf_runtime_error("Unable to allocate memory "
                          "for internal buffer in sendget().");
    }

    if (src_same_image)
      src_t_buff = src->base_addr;
    else if (!same_type_and_kind)
    {
      if ((free_src_t_buff = (((src_t_buff = alloca(src_size))) == NULL)))
      {
        src_t_buff = malloc(src_size);
        if (src_t_buff == NULL)
          caf_runtime_error("Unable to allocate memory "
                            "for internal buffer in sendget().");
      }
    }

    if (!src_same_image)
      CAF_Win_lock(MPI_LOCK_SHARED, src_remote_image, *p);
    for (i = 0; i < size; ++i)
    {
      ptrdiff_t array_offset_sr = 0, extent = 1, tot_ext = 1;
      if (src_vector == NULL)
      {
        for (j = 0; j < src_rank - 1; ++j)
        {
          extent = src->dim[j]._ubound - src->dim[j].lower_bound + 1;
          array_offset_sr += ((i / tot_ext) % extent) * src->dim[j]._stride;
          tot_ext *= extent;
        }

        array_offset_sr += (i / tot_ext) * src->dim[src_rank - 1]._stride;
      }
      else
      {
#define KINDCASE(kind, type)                                        \
case kind:                                                          \
  array_offset_sr = ((ptrdiff_t)                                    \
    ((type *)src_vector->u.v.vector)[i] - src->dim[0].lower_bound); \
  break
        switch (src_vector->u.v.kind)
        {
          KINDCASE(1, int8_t);
          KINDCASE(2, int16_t);
          KINDCASE(4, int32_t);
          KINDCASE(8, int64_t);
#ifdef HAVE_GFC_INTEGER_16
          KINDCASE(16, __int128);
#endif
          default:
            caf_runtime_error(unreachable);
            return;
        }
      }
#undef KINDCASE

      src_offset = array_offset_sr * src_size;
      void *dst = (void *)((char *) dst_t_buff + i * dst_size);

      if (!src_same_image)
      {
        // Do the more likely first.
        dprint("kind(dst) = %d, el_sz(dst) = %zd, "
               "kind(src) = %d, el_sz(src) = %zd, lb(dst) = %zd.\n",
               dst_kind, dst_size, src_kind, src_size, src->dim[0].lower_bound);
        if (same_type_and_kind)
        {
          const size_t trans_size = (src_size < dst_size) ? src_size : dst_size;
          ierr = MPI_Get(dst, trans_size, MPI_BYTE, src_remote_image,
                         offset_g + src_offset, trans_size, MPI_BYTE, *p);
          chk_err(ierr);
          if (pad_str)
            memcpy((void *)((char *)dst + src_size), pad_str,
                   dst_size - src_size);
        }
        else if (dst_type == BT_CHARACTER)
        {
          ierr = MPI_Get(src_t_buff, src_size, MPI_BYTE, src_remote_image,
                         offset_g + src_offset, src_size, MPI_BYTE, *p);
          chk_err(ierr);
          copy_char_to_self(src_t_buff, src_type, src_size, src_kind,
                            dst, dst_type, dst_size, dst_kind, 1, true);
        }
        else
        {
          ierr = MPI_Get(src_t_buff, src_size, MPI_BYTE, src_remote_image,
                         offset_g + src_offset, src_size, MPI_BYTE, *p);
          chk_err(ierr);
          convert_type(dst, dst_type, dst_kind, src_t_buff, src_type,
                       src_kind, stat);
        }
      }
      else
      {
        if (!mrt)
        {
          dprint("strided same_image, no temp,  for i = %zd, "
                 "src_offset = %zd, offset = %zd.\n",
                 i, src_offset, offset_g);
          if (same_type_and_kind)
            memmove(dst, src->base_addr + src_offset, src_size);
          else
            convert_type(dst, dst_type, dst_kind,
                         src->base_addr + src_offset, src_type, src_kind, stat);
        }
        else
        {
          dprint("strided same_image, *WITH* temp, for i = %zd.\n", i);
          if (same_type_and_kind)
            memmove(dst, src->base_addr + src_offset, src_size);
          else
            convert_type(dst, dst_type, dst_kind,
                         src->base_addr + src_offset, src_type, src_kind, stat);
        }
      }

#ifndef WITH_FAILED_IMAGES
      if (ierr != 0)
      {
        caf_runtime_error("MPI Error: %d", ierr);
        return;
      }
#endif
    }
    if (!src_same_image)
      CAF_Win_unlock(src_remote_image, *p);
  }

  p = TOKEN(token_s);
  /* Now transfer data to the remote dest. */
  if (dst_contiguous && dst_vector == NULL)
  {
    if (dst_same_image)
      memmove(dest->base_addr, dst_t_buff, dst_size * size);
    else
    {
      CAF_Win_lock(MPI_LOCK_EXCLUSIVE, dst_remote_image, *p);
      const size_t trans_size = size * dst_size;
      ierr = MPI_Put(dst_t_buff, trans_size, MPI_BYTE, dst_remote_image,
                     offset_s, trans_size, MPI_BYTE, *p); chk_err(ierr);
#ifdef CAF_MPI_LOCK_UNLOCK
      MPI_Win_unlock(dst_remote_image, *p);
#elif NONBLOCKING_PUT
      /* Pending puts init */
      if (pending_puts == NULL)
      {
        pending_puts = calloc(1,sizeof(win_sync));
        pending_puts->next = NULL;
        pending_puts->win = token_s;
        pending_puts->img = dst_remote_image;
        last_elem = pending_puts;
        last_elem->next = NULL;
      }
      else
      {
        last_elem->next = calloc(1,sizeof(win_sync));
        last_elem = last_elem->next;
        last_elem->win = token_s;
        last_elem->img = dst_remote_image;
        last_elem->next = NULL;
      }
#else
      ierr = MPI_Win_flush(dst_remote_image, *p); chk_err(ierr);
#endif // CAF_MPI_LOCK_UNLOCK
    }
  }
#ifdef STRIDED
  else if (!dst_same_image && same_type_and_kind && dst_type != BT_CHARACTER)
  {
    /* For strided copy, no type and kind conversion, copy to self or
     * character arrays are supported. */
    MPI_Datatype dt_s, dt_d, base_type_dst;
    int *arr_bl, *arr_dsp_d;

    selectType(dst_size, &base_type_dst);

    if (dst_rank == 1)
    {
      if (dst_vector == NULL)
      {
        ierr = MPI_Type_vector(size, 1, dest->dim[0]._stride, base_type_dst,
                               &dt_d); chk_err(ierr);
      }
      else
      {
        arr_bl = calloc(size, sizeof(int));
        arr_dsp_d = calloc(size, sizeof(int));

        dprint("Setting up strided vector index.\n");
#define KINDCASE(kind, type)                                            \
case kind:                                                              \
  for (i = 0; i < size; ++i)                                            \
  {                                                                     \
    arr_dsp_d[i] = ((ptrdiff_t)                                         \
      ((type *) dst_vector->u.v.vector)[i] - dest->dim[0].lower_bound); \
    arr_bl[i] = 1;                                                      \
  }                                                                     \
  break
        switch (dst_vector->u.v.kind)
        {
          KINDCASE(1, int8_t);
          KINDCASE(2, int16_t);
          KINDCASE(4, int32_t);
          KINDCASE(8, int64_t);
#ifdef HAVE_GFC_INTEGER_16
          KINDCASE(16, __int128);
#endif
        default:
          caf_runtime_error(unreachable);
          return;
        }
#undef KINDCASE
        ierr = MPI_Type_indexed(size, arr_bl, arr_dsp_d, base_type_dst, &dt_d);
        chk_err(ierr);
        free(arr_bl);
        free(arr_dsp_d);
      }
      ierr = MPI_Type_vector(size, 1, 1, base_type_dst, &dt_s); chk_err(ierr);
    }
    else
    {
      arr_bl = calloc(size, sizeof(int));
      arr_dsp_d = calloc(size, sizeof(int));

      for (i = 0; i < size; ++i)
      {
        arr_bl[i] = 1;
      }

      for (i = 0; i < size; ++i)
      {
        ptrdiff_t array_offset_dst = 0, extent = 1, tot_ext = 1;
        if (dst_vector == NULL)
        {
          for (j = 0; j < dst_rank - 1; ++j)
          {
            extent = dest->dim[j]._ubound - dest->dim[j].lower_bound + 1;
            array_offset_dst += ((i / tot_ext) % extent) * dest->dim[j]._stride;
            tot_ext *= extent;
          }

          array_offset_dst += (i / tot_ext) * dest->dim[dst_rank - 1]._stride;
        }
        else
        {
#define KINDCASE(kind, type)                                            \
case kind:                                                              \
  array_offset_dst = ((ptrdiff_t)                                       \
    ((type *) dst_vector->u.v.vector)[i] - dest->dim[0].lower_bound);   \
  break
          switch (dst_vector->u.v.kind)
          {
            KINDCASE(1, int8_t);
            KINDCASE(2, int16_t);
            KINDCASE(4, int32_t);
            KINDCASE(8, int64_t);
#ifdef HAVE_GFC_INTEGER_16
            KINDCASE(16, __int128);
#endif
            default:
              caf_runtime_error(unreachable);
              return;
          }
#undef KINDCASE
        }
        arr_dsp_d[i] = array_offset_dst;
      }

      ierr = MPI_Type_vector(size, 1, 1, base_type_dst, &dt_s); chk_err(ierr);
      ierr = MPI_Type_indexed(size, arr_bl, arr_dsp_d, base_type_dst, &dt_d);
      chk_err(ierr);

      free(arr_bl);
      free(arr_dsp_d);
    }

    ierr = MPI_Type_commit(&dt_s); chk_err(ierr);
    ierr = MPI_Type_commit(&dt_d); chk_err(ierr);

    CAF_Win_lock(MPI_LOCK_EXCLUSIVE, dst_remote_image, *p);
    ierr = MPI_Put(dst_t_buff, 1, dt_s, dst_remote_image, offset_s, 1,
                   dt_d, *p); chk_err(ierr);
    CAF_Win_unlock(dst_remote_image, *p);

#ifdef WITH_FAILED_IMAGES
    check_image_health(image_index_s, stat);

    if (!stat && ierr == STAT_FAILED_IMAGE)
      terminate_internal(ierr, 1);

    if (stat)
      *stat = ierr;
#else
    if (ierr != 0)
    {
      terminate_internal(ierr, 1);
      return;
    }
#endif
    ierr = MPI_Type_free(&dt_s); chk_err(ierr);
    ierr = MPI_Type_free(&dt_d); chk_err(ierr);
  }
#endif // STRIDED
  else
  {
    if (!dst_same_image)
      CAF_Win_lock(MPI_LOCK_EXCLUSIVE, dst_remote_image, *p);
    for (i = 0; i < size; ++i)
    {
      ptrdiff_t array_offset_dst = 0, extent = 1, tot_ext = 1;
      if (dst_vector == NULL)
      {
        for (j = 0; j < dst_rank - 1; ++j)
        {
          extent = dest->dim[j]._ubound - dest->dim[j].lower_bound + 1;
          array_offset_dst += ((i / tot_ext) % extent) * dest->dim[j]._stride;
          tot_ext *= extent;
        }

        array_offset_dst += (i / tot_ext) * dest->dim[dst_rank - 1]._stride;
      }
      else
      {
#define KINDCASE(kind, type)                                         \
case kind:                                                           \
  array_offset_dst = ((ptrdiff_t)                                    \
    ((type *)dst_vector->u.v.vector)[i] - dest->dim[0].lower_bound); \
  break
        switch (dst_vector->u.v.kind)
        {
          KINDCASE(1, int8_t);
          KINDCASE(2, int16_t);
          KINDCASE(4, int32_t);
          KINDCASE(8, int64_t);
#ifdef HAVE_GFC_INTEGER_16
          KINDCASE(16, __int128);
#endif
          default:
            caf_runtime_error(unreachable);
            return;
        }
      }
#undef KINDCASE
      dst_offset = array_offset_dst * dst_size;

      void *sr = (void *)((char *)dst_t_buff + i * dst_size);

      if (!dst_same_image)
      {
        // Do the more likely first.
        ierr = MPI_Put(sr, dst_size, MPI_BYTE, dst_remote_image,
                       offset_s + dst_offset, dst_size, MPI_BYTE, *p);
         chk_err(ierr);
      }
      else
        memmove(dest->base_addr + dst_offset, sr, dst_size);

#ifndef WITH_FAILED_IMAGES
      if (ierr != 0)
      {
        caf_runtime_error("MPI Error: %d", ierr);
        return;
      }
#endif
    } /* for */
    if (!dst_same_image)
      CAF_Win_unlock(dst_remote_image, *p);
  }

  /* Free memory, when not allocated on stack. */
  if (free_src_t_buff)
    free(src_t_buff);
  if (free_dst_t_buff)
    free(dst_t_buff);
  if (free_pad_str)
    free(pad_str);

#ifdef WITH_FAILED_IMAGES
  /* Catch failed images, when failed image support is active. */
  check_image_health(image_index_g , stat);
  check_image_health(image_index_s , stat);
#endif

  if (ierr != MPI_SUCCESS)
  {
    int mpi_error;
    MPI_Error_class(ierr, &mpi_error);
    if (stat)
      *stat = mpi_error;
    else
    {
      int error_len = 2048;
      char error_str[error_len];
      strcpy(error_str, "MPI-error: ");
      MPI_Error_string(mpi_error, &error_str[11], &error_len);
      error_stop_str(error_str, error_len + 11, false);
    }
  }
}


/* Send array data from src to dest on a remote image.
 * The argument mrt means may_require_temporary */

void
PREFIX(send) (caf_token_t token, size_t offset, int image_index,
              gfc_descriptor_t *dest, caf_vector_t *dst_vector,
              gfc_descriptor_t *src, int dst_kind, int src_kind,
              bool mrt, int *pstat)
{
  int j, ierr = 0;
  size_t i, size;
  ptrdiff_t dimextent;
  const int
    src_rank = GFC_DESCRIPTOR_RANK(src),
    dst_rank = GFC_DESCRIPTOR_RANK(dest);
  const size_t
    src_size = GFC_DESCRIPTOR_SIZE(src),
    dst_size = GFC_DESCRIPTOR_SIZE(dest);
  const int
    src_type = GFC_DESCRIPTOR_TYPE(src),
    dst_type = GFC_DESCRIPTOR_TYPE(dest);
  const bool
    src_contiguous = PREFIX(is_contiguous) (src),
    dst_contiguous = PREFIX(is_contiguous) (dest);
  const bool
    same_image = caf_this_image == image_index,
    same_type_and_kind = dst_type == src_type && dst_kind == src_kind;

  MPI_Win *p = TOKEN(token);
  ptrdiff_t dst_offset = 0;
  void *pad_str = NULL, *t_buff = NULL;
  bool free_pad_str = false, free_t_buff = false;
  const bool dest_char_array_is_longer
      = dst_type == BT_CHARACTER && dst_size > src_size && !same_image;
  int remote_image = image_index - 1;
  if (!same_image)
  {
    MPI_Group current_team_group, win_group;
    ierr = MPI_Comm_group(CAF_COMM_WORLD, &current_team_group); chk_err(ierr);
    ierr = MPI_Win_get_group(*p, &win_group); chk_err(ierr);
    ierr = MPI_Group_translate_ranks(current_team_group, 1,
                                     (int[]){remote_image}, win_group,
                                     &remote_image); chk_err(ierr);
    ierr = MPI_Group_free(&current_team_group); chk_err(ierr);
    ierr = MPI_Group_free(&win_group); chk_err(ierr);
  }

  /* Ensure stat is always set. */
#ifdef GCC_GE_7
  int * stat = pstat;
  if (stat)
    *stat = 0;
#else
  /* Gcc prior to 7.0 does not have stat here. */
  int * stat = NULL;
#endif

  size = 1;
  for (j = 0; j < dst_rank; ++j)
  {
    dimextent = dest->dim[j]._ubound - dest->dim[j].lower_bound + 1;
    if (dimextent < 0)
      dimextent = 0;
    size *= dimextent;
  }

  if (size == 0)
    return;

  dprint("dst_vector = %p, image_index = %d, offset = %zd.\n",
         dst_vector, image_index, offset);
  check_image_health(image_index, stat);

  /* For char arrays: create the padding array, when dst is longer than src. */
  if (dest_char_array_is_longer)
  {
    const size_t pad_num = (dst_size / dst_kind) - (src_size / src_kind);
    const size_t pad_sz = pad_num * dst_kind;
    /* For big arrays alloca() may not be able to get the memory on the stack.
     * Use a regular malloc then. */
    if ((free_pad_str = ((pad_str = alloca(pad_sz)) == NULL)))
    {
      pad_str = malloc(pad_sz);
      if (t_buff == NULL)
        caf_runtime_error("Unable to allocate memory "
                          "for internal buffer in send().");
    }
    if (dst_kind == 1)
      memset(pad_str, ' ', pad_num);
    else /* dst_kind == 4. */
    {
      for (int32_t *it = (int32_t *) pad_str,
           *itEnd = ((int32_t *) pad_str) + pad_num; it < itEnd; ++it)
      {
        *it = (int32_t) ' ';
      }
    }
  }

  if (src_contiguous && dst_contiguous && dst_vector == NULL)
  {
    if (same_image)
    {
      dprint("in caf_this == image_index, size = %zd, dst_kind = %d, "
             "src_kind = %d\n", size, dst_kind, src_kind);
      if (dst_type == BT_CHARACTER)
        /* The size is encoded in the descriptor's type for char arrays. */
        copy_char_to_self(src->base_addr, src_type, src_size, src_kind,
                          dest->base_addr, dst_type, dst_size, dst_kind,
                          size, src_rank == 0);
      else
        copy_to_self(src, src_kind, dest, dst_kind, size, stat);
      return;
    }
    else
    {
      CAF_Win_lock(MPI_LOCK_EXCLUSIVE, remote_image, *p);
      if (dst_kind != src_kind || dest_char_array_is_longer || src_rank == 0)
      {
        if ((free_t_buff = ((t_buff = alloca(dst_size * size)) == NULL)))
        {
          t_buff = malloc(dst_size * size);
          if (t_buff == NULL)
            caf_runtime_error("Unable to allocate memory "
                              "for internal buffer in send().");
        }
      }

      if ((same_type_and_kind && dst_rank == src_rank)
          || dst_type == BT_CHARACTER)
        {
          if (dest_char_array_is_longer
              || (dst_kind != src_kind && dst_type == BT_CHARACTER))
          {
            copy_char_to_self(src->base_addr, src_type, src_size,
                              src_kind, t_buff, dst_type, dst_size,
                              dst_kind, size, src_rank == 0);
            ierr = MPI_Put(t_buff, dst_size, MPI_BYTE, remote_image,
                           offset, dst_size, MPI_BYTE, *p); chk_err(ierr);
          }
          else
          {
            const size_t trans_size =
              ((dst_size > src_size) ? src_size : dst_size) * size;
            ierr = MPI_Put(src->base_addr, trans_size, MPI_BYTE, remote_image,
                           offset, trans_size, MPI_BYTE, *p); chk_err(ierr);
          }
        }
      else
      {
        convert_with_strides(t_buff, dst_type, dst_kind, dst_size,
                             src->base_addr, src_type, src_kind,
                             (src_rank > 0) ? src_size: 0, size, stat);
        ierr = MPI_Put(t_buff, dst_size * size, MPI_BYTE, remote_image,
                       offset, dst_size * size, MPI_BYTE, *p); chk_err(ierr);
      }
#ifdef CAF_MPI_LOCK_UNLOCK
      ierr = MPI_Win_unlock(remote_image, *p); chk_err(ierr);
#elif NONBLOCKING_PUT
      /* Pending puts init */
      if (pending_puts == NULL)
      {
        pending_puts = calloc(1,sizeof(win_sync));
        pending_puts->next = NULL;
        pending_puts->win = token;
        pending_puts->img = remote_image;
        last_elem = pending_puts;
        last_elem->next = NULL;
      }
      else
      {
        last_elem->next = calloc(1,sizeof(win_sync));
        last_elem = last_elem->next;
        last_elem->win = token;
        last_elem->img = remote_image;
        last_elem->next = NULL;
      }
#else
      ierr = MPI_Win_flush(remote_image, *p); chk_err(ierr);
#endif // CAF_MPI_LOCK_UNLOCK
    }
  }

#ifdef STRIDED
  else if (!same_image && same_type_and_kind && dst_type != BT_CHARACTER)
  {
    /* For strided copy, no type and kind conversion, copy to self or
     * character arrays are supported. */
    MPI_Datatype dt_s, dt_d, base_type_src, base_type_dst;
    int *arr_bl, *arr_dsp_s, *arr_dsp_d;

    selectType(src_size, &base_type_src);
    selectType(dst_size, &base_type_dst);

    if (dst_rank == 1)
    {
      if (dst_vector == NULL)
      {
        dprint("Setting up mpi datatype vector with "
               "stride %d, size %d and offset %d.\n",
               dest->dim[0]._stride, size, offset);
        ierr = MPI_Type_vector(size, 1, dest->dim[0]._stride, base_type_dst,
                               &dt_d); chk_err(ierr);
      }
      else
      {
        arr_bl = calloc(size, sizeof(int));
        arr_dsp_d = calloc(size, sizeof(int));

        dprint("Setting up strided vector index.\n");
#define KINDCASE(kind, type)                                            \
case kind:                                                              \
  for (i = 0; i < size; ++i)                                            \
  {                                                                     \
    arr_dsp_d[i] = ((ptrdiff_t)                                         \
      ((type *) dst_vector->u.v.vector)[i] - dest->dim[0].lower_bound); \
    arr_bl[i] = 1;                                                      \
  }                                                                     \
  break
        switch (dst_vector->u.v.kind)
        {
          KINDCASE(1, int8_t);
          KINDCASE(2, int16_t);
          KINDCASE(4, int32_t);
          KINDCASE(8, int64_t);
#ifdef HAVE_GFC_INTEGER_16
          KINDCASE(16, __int128);
#endif
        default:
          caf_runtime_error(unreachable);
          return;
        }
#undef KINDCASE
        ierr = MPI_Type_indexed(size, arr_bl, arr_dsp_d, base_type_dst, &dt_d);
        chk_err(ierr);
        free(arr_bl);
        free(arr_dsp_d);
      }
      MPI_Type_vector(size, 1, src->dim[0]._stride, base_type_src, &dt_s);
    }
    else
    {
      arr_bl = calloc(size, sizeof(int));
      arr_dsp_s = calloc(size, sizeof(int));
      arr_dsp_d = calloc(size, sizeof(int));

      for (i = 0; i < size; ++i)
      {
        arr_bl[i] = 1;
      }

      for (i = 0; i < size; ++i)
      {
        ptrdiff_t array_offset_dst = 0, extent = 1, tot_ext = 1;
        if (dst_vector == NULL)
        {
          for (j = 0; j < dst_rank - 1; ++j)
          {
            extent = dest->dim[j]._ubound - dest->dim[j].lower_bound + 1;
            array_offset_dst += ((i / tot_ext) % extent) * dest->dim[j]._stride;
            tot_ext *= extent;
          }

          array_offset_dst += (i / tot_ext) * dest->dim[dst_rank - 1]._stride;
        }
        else
        {
#define KINDCASE(kind, type)                                          \
case kind:                                                            \
  array_offset_dst = ((ptrdiff_t)                                     \
    ((type *) dst_vector->u.v.vector)[i] - dest->dim[0].lower_bound); \
  break
          switch (dst_vector->u.v.kind)
          {
            KINDCASE(1, int8_t);
            KINDCASE(2, int16_t);
            KINDCASE(4, int32_t);
            KINDCASE(8, int64_t);
#ifdef HAVE_GFC_INTEGER_16
            KINDCASE(16, __int128);
#endif
            default:
              caf_runtime_error(unreachable);
              return;
          }
#undef KINDCASE
        }
        arr_dsp_d[i] = array_offset_dst;

        if (src_rank != 0)
        {
          ptrdiff_t array_offset_sr = 0;
          extent = 1;
          tot_ext = 1;
          for (j = 0; j < src_rank - 1; ++j)
          {
            extent = src->dim[j]._ubound - src->dim[j].lower_bound + 1;
            array_offset_sr += ((i / tot_ext) % extent) * src->dim[j]._stride;
            tot_ext *= extent;
          }

          array_offset_sr += (i / tot_ext) * src->dim[src_rank - 1]._stride;
          arr_dsp_s[i] = array_offset_sr;
        }
        else
          arr_dsp_s[i] = 0;
      }

      ierr = MPI_Type_indexed(size, arr_bl, arr_dsp_s, base_type_src, &dt_s);
      chk_err(ierr);
      ierr = MPI_Type_indexed(size, arr_bl, arr_dsp_d, base_type_dst, &dt_d);
      chk_err(ierr);

      free(arr_bl);
      free(arr_dsp_s);
      free(arr_dsp_d);
    }

    ierr = MPI_Type_commit(&dt_s); chk_err(ierr);
    ierr = MPI_Type_commit(&dt_d); chk_err(ierr);

    CAF_Win_lock(MPI_LOCK_EXCLUSIVE, remote_image, *p);
    ierr = MPI_Put(src->base_addr, 1, dt_s, remote_image, offset, 1, dt_d, *p);
    chk_err(ierr);
    CAF_Win_unlock(remote_image, *p);

#ifdef WITH_FAILED_IMAGES
    check_image_health(image_index, stat);

    if (!stat && ierr == STAT_FAILED_IMAGE)
      terminate_internal(ierr, 1);

    if (stat)
      *stat = ierr;
#else
    if (ierr != 0)
    {
      terminate_internal(ierr, 1);
      return;
    }
#endif
    ierr = MPI_Type_free(&dt_s); chk_err(ierr);
    ierr = MPI_Type_free(&dt_d); chk_err(ierr);
  }
#endif // STRIDED
  else
  {
    if (same_image && mrt)
    {
      if ((free_t_buff = (((t_buff = alloca(dst_size * size))) == NULL)))
      {
        t_buff = malloc(dst_size * size);
        if (t_buff == NULL)
          caf_runtime_error("Unable to allocate memory "
                            "for internal buffer in send().");
      }
    }
    else if (!same_type_and_kind && !same_image)
    {
      if ((free_t_buff = (((t_buff = alloca(dst_size))) == NULL)))
      {
        t_buff = malloc(dst_size);
        if (t_buff == NULL)
          caf_runtime_error("Unable to allocate memory "
                            "for internal buffer in send().");
      }
    }

    if (!same_image)
      CAF_Win_lock(MPI_LOCK_EXCLUSIVE, remote_image, *p);
    for (i = 0; i < size; ++i)
    {
      ptrdiff_t array_offset_dst = 0, extent = 1, tot_ext = 1;
      if (!same_image || !mrt)
      {
        /* For same image and may require temp, the dst_offset is
         * computed on storage. */
        if (dst_vector == NULL)
        {
          for (j = 0; j < dst_rank - 1; ++j)
          {
            extent = dest->dim[j]._ubound - dest->dim[j].lower_bound + 1;
            array_offset_dst += ((i / tot_ext) % extent) * dest->dim[j]._stride;
            tot_ext *= extent;
          }
          array_offset_dst += (i / tot_ext) * dest->dim[dst_rank - 1]._stride;
        }
        else
        {
#define KINDCASE(kind, type)                                          \
case kind:                                                            \
  array_offset_dst = ((ptrdiff_t)                                     \
    ((type *)dst_vector->u.v.vector)[i] - dest->dim[0].lower_bound);  \
  break
          switch (dst_vector->u.v.kind)
          {
            KINDCASE(1, int8_t);
            KINDCASE(2, int16_t);
            KINDCASE(4, int32_t);
            KINDCASE(8, int64_t);
  #ifdef HAVE_GFC_INTEGER_16
            KINDCASE(16, __int128);
  #endif
            default:
              caf_runtime_error(unreachable);
              return;
          }
        }
        dst_offset = array_offset_dst * dst_size;
      }

      void *sr;
      if (src_rank != 0)
      {
        ptrdiff_t array_offset_sr = 0;
        extent = 1;
        tot_ext = 1;
        for (j = 0; j < src_rank - 1; ++j)
        {
          extent = src->dim[j]._ubound - src->dim[j].lower_bound + 1;
          array_offset_sr += ((i / tot_ext) % extent) * src->dim[j]._stride;
          tot_ext *= extent;
        }

        array_offset_sr += (i / tot_ext) * src->dim[dst_rank - 1]._stride;
        sr = (void *)((char *)src->base_addr + array_offset_sr * src_size);
      }
      else
        sr = src->base_addr;

      if (!same_image)
      {
        // Do the more likely first.
        dprint("kind(dst) = %d, el_sz(dst) = %zd, "
               "kind(src) = %d, el_sz(src) = %zd, lb(dst) = %zd.\n",
               dst_kind, dst_size, src_kind,
               src_size, dest->dim[0].lower_bound);
        if (same_type_and_kind)
        {
          const size_t trans_size = (src_size < dst_size) ? src_size : dst_size;
          ierr = MPI_Put(sr, trans_size, MPI_BYTE, remote_image,
                         offset + dst_offset, trans_size, MPI_BYTE, *p);
          chk_err(ierr);
          if (pad_str)
          {
            ierr = MPI_Put(pad_str, dst_size - src_size, MPI_BYTE,
                           remote_image, offset + dst_offset + src_size,
                           dst_size - src_size, MPI_BYTE, *p); chk_err(ierr);
          }
        }
        else if (dst_type == BT_CHARACTER)
        {
          copy_char_to_self(sr, src_type, src_size, src_kind,
                            t_buff, dst_type, dst_size, dst_kind, 1, true);
          ierr = MPI_Put(t_buff, dst_size, MPI_BYTE, remote_image,
                         offset + dst_offset, dst_size, MPI_BYTE, *p);
          chk_err(ierr);
        }
        else
        {
          convert_type(t_buff, dst_type, dst_kind,
                       sr, src_type, src_kind, stat);
          ierr = MPI_Put(t_buff, dst_size, MPI_BYTE, remote_image,
                         offset + dst_offset, dst_size, MPI_BYTE, *p);
          chk_err(ierr);
        }
      }
      else
      {
        if (!mrt)
        {
          dprint("strided same_image, no temp, for i = %zd, "
                 "dst_offset = %zd, offset = %zd.\n",
                 i, dst_offset, offset);
          if (same_type_and_kind)
            memmove(dest->base_addr + dst_offset, sr, src_size);
          else
            convert_type(dest->base_addr + dst_offset, dst_type,
                         dst_kind, sr, src_type, src_kind, stat);
        }
        else
        {
          dprint("strided same_image, *WITH* temp, for i = %zd.\n", i);
          if (same_type_and_kind)
            memmove(t_buff + i * dst_size, sr, src_size);
          else
            convert_type(t_buff + i * dst_size, dst_type, dst_kind,
                         sr, src_type, src_kind, stat);
        }
      }

#ifndef WITH_FAILED_IMAGES
      if (ierr != 0)
      {
        caf_runtime_error("MPI Error: %d", ierr);
        return;
      }
#endif
    }
    if (!same_image)
      CAF_Win_unlock(remote_image, *p);

    if (same_image && mrt)
    {
      for (i = 0; i < size; ++i)
      {
        ptrdiff_t array_offset_dst = 0, extent = 1, tot_ext = 1;
        if (dst_vector == NULL)
        {
          for (j = 0; j < dst_rank - 1; ++j)
          {
            extent = dest->dim[j]._ubound - dest->dim[j].lower_bound + 1;
            array_offset_dst += ((i / tot_ext) % extent) * dest->dim[j]._stride;
            tot_ext *= extent;
          }

          array_offset_dst += (i / tot_ext) * dest->dim[dst_rank - 1]._stride;
        }
        else
        {
          switch (dst_vector->u.v.kind)
          {
            // KINDCASE is defined above.
            KINDCASE(1, int8_t);
            KINDCASE(2, int16_t);
            KINDCASE(4, int32_t);
            KINDCASE(8, int64_t);
#ifdef HAVE_GFC_INTEGER_16
            KINDCASE(16, __int128);
#endif
            default:
              caf_runtime_error(unreachable);
              return;
          }
#undef KINDCASE
        }
        dst_offset = array_offset_dst * dst_size;
        memmove(dest->base_addr + dst_offset, t_buff + i * dst_size, dst_size);
      }
    }
  }

  /* Free memory, when not allocated on stack. */
  if (free_t_buff)
    free(t_buff);
  if (free_pad_str)
    free(pad_str);

#ifdef WITH_FAILED_IMAGES
  /* Catch failed images, when failed image support is active. */
  check_image_health(image_index , stat);
#endif

  if (ierr != MPI_SUCCESS)
  {
    int mpi_error;
    MPI_Error_class(ierr, &mpi_error);
    if (stat)
      *stat = mpi_error;
    else
    {
      int error_len = 2048;
      char error_str[error_len];
      strcpy(error_str, "MPI-error: ");
      MPI_Error_string(mpi_error, &error_str[11], &error_len);
      error_stop_str(error_str, error_len + 11, false);
    }
  }
}


/* Get array data from a remote src to a local dest. */

void
PREFIX(get) (caf_token_t token, size_t offset, int image_index,
             gfc_descriptor_t *src, caf_vector_t *src_vector,
             gfc_descriptor_t *dest, int src_kind, int dst_kind,
             bool mrt, int *pstat)
{
  int j, ierr = 0;
  size_t i, size;
  const int
    src_rank = GFC_DESCRIPTOR_RANK(src),
    dst_rank = GFC_DESCRIPTOR_RANK(dest);
  const size_t
    src_size = GFC_DESCRIPTOR_SIZE(src),
    dst_size = GFC_DESCRIPTOR_SIZE(dest);
  const int
    src_type = GFC_DESCRIPTOR_TYPE(src),
    dst_type = GFC_DESCRIPTOR_TYPE(dest);
  const bool
    src_contiguous = PREFIX(is_contiguous) (src),
    dst_contiguous = PREFIX(is_contiguous) (dest);
  const bool
    same_image = caf_this_image == image_index,
    same_type_and_kind = dst_type == src_type && dst_kind == src_kind;

  MPI_Win *p = TOKEN(token);
  ptrdiff_t dimextent, src_offset = 0;
  void *pad_str = NULL, *t_buff = NULL;
  bool free_pad_str = false, free_t_buff = false;
  const bool dest_char_array_is_longer
      = dst_type == BT_CHARACTER && dst_size > src_size && !same_image;
  int remote_image = image_index - 1;
  if (!same_image)
  {
    MPI_Group current_team_group, win_group;
    ierr = MPI_Comm_group(CAF_COMM_WORLD, &current_team_group); chk_err(ierr);
    ierr = MPI_Win_get_group(*p, &win_group); chk_err(ierr);
    ierr = MPI_Group_translate_ranks(current_team_group, 1,
                                     (int[]){remote_image}, win_group,
                                     &remote_image); chk_err(ierr);
    ierr = MPI_Group_free(&current_team_group); chk_err(ierr);
    ierr = MPI_Group_free(&win_group); chk_err(ierr);
  }

  /* Ensure stat is always set. */
#ifdef GCC_GE_7
  int * stat = pstat;
  if (stat)
    *stat = 0;
#else
  /* Gcc prior to 7.0 does not have stat here. */
  int * stat = NULL;
#endif

  size = 1;
  for (j = 0; j < dst_rank; ++j)
  {
    dimextent = dest->dim[j]._ubound - dest->dim[j].lower_bound + 1;
    if (dimextent < 0)
      dimextent = 0;
    size *= dimextent;
  }

  if (size == 0)
    return;

  dprint("src_vector = %p, image_index = %d, offset = %zd.\n",
         src_vector, image_index, offset);
  check_image_health(image_index, stat);

  /* For char arrays: create the padding array, when dst is longer than src. */
  if (dest_char_array_is_longer)
  {
    const size_t pad_num = (dst_size / dst_kind) - (src_size / src_kind);
    const size_t pad_sz = pad_num * dst_kind;
    /* For big arrays alloca() may not be able to get the memory on the stack.
     * Use a regular malloc then. */
    if ((free_pad_str = ((pad_str = alloca(pad_sz)) == NULL)))
    {
      pad_str = malloc(pad_sz);
      if (t_buff == NULL)
        caf_runtime_error("Unable to allocate memory "
                          "for internal buffer in get().");
    }
    if (dst_kind == 1)
      memset(pad_str, ' ', pad_num);
    else /* dst_kind == 4. */
    {
      for (int32_t *it = (int32_t *) pad_str,
           *itEnd = ((int32_t *) pad_str) + pad_num; it < itEnd; ++it)
      {
        *it = (int32_t) ' ';
      }
    }
  }

  if (src_contiguous && dst_contiguous && src_vector == NULL)
  {
    if (same_image)
    {
      dprint("in caf_this == image_index, size = %zd, "
             "dst_kind = %d, src_kind = %d\n",
             size, dst_kind, src_kind);
      if (dst_type == BT_CHARACTER)
        /* The size is encoded in the descriptor's type for char arrays. */
        copy_char_to_self(src->base_addr, src_type, src_size, src_kind,
                          dest->base_addr, dst_type, dst_size, dst_kind,
                          size, src_rank == 0);
      else
        copy_to_self(src, src_kind, dest, dst_kind, size, stat);
      return;
    }
    else
    {
      CAF_Win_lock(MPI_LOCK_SHARED, remote_image, *p);
      if (dst_kind != src_kind || dest_char_array_is_longer || src_rank == 0)
      {
        if ((free_t_buff = ((t_buff = alloca(src_size * size)) == NULL)))
        {
          t_buff = malloc(src_size * size);
          if (t_buff == NULL)
            caf_runtime_error("Unable to allocate memory "
                              "for internal buffer in get().");
        }
      }

      if ((same_type_and_kind && dst_rank == src_rank)
          || dst_type == BT_CHARACTER)
      {
        if (!dest_char_array_is_longer
            && (dst_kind == src_kind || dst_type != BT_CHARACTER))
        {
          const size_t trans_size =
            ((dst_size > src_size) ? src_size : dst_size) * size;
          ierr = MPI_Get(dest->base_addr, trans_size, MPI_BYTE, remote_image,
                         offset, trans_size, MPI_BYTE, *p); chk_err(ierr);
        }
        else
        {
          ierr = MPI_Get(t_buff, src_size, MPI_BYTE, remote_image,
                         offset, src_size, MPI_BYTE, *p); chk_err(ierr);
          copy_char_to_self(t_buff, src_type, src_size, src_kind,
                            dest->base_addr, dst_type, dst_size,
                            dst_kind, size, src_rank == 0);
        }
      }
      else
      {
        ierr = MPI_Get(t_buff, src_size * size, MPI_BYTE, remote_image, offset,
                       src_size * size, MPI_BYTE, *p); chk_err(ierr);
        convert_with_strides(dest->base_addr, dst_type, dst_kind, dst_size,
                             t_buff, src_type, src_kind,
                             (src_rank > 0) ? src_size: 0, size, stat);
      }
      CAF_Win_unlock(remote_image, *p);
    }
  }
#ifdef STRIDED
  else if (!same_image && same_type_and_kind && dst_type != BT_CHARACTER)
  {
    /* For strided copy, no type and kind conversion, copy to self or
     * character arrays are supported. */
    MPI_Datatype dt_s, dt_d, base_type_src, base_type_dst;
    int *arr_bl;
    int *arr_dsp_s, *arr_dsp_d;

    selectType(src_size, &base_type_src);
    selectType(dst_size, &base_type_dst);

    if (src_rank == 1)
    {
      if (src_vector == NULL)
      {
        dprint("Setting up mpi datatype vector with stride %d, "
               "size %d and offset %d.\n",
               src->dim[0]._stride, size, offset);
        ierr = MPI_Type_vector(size, 1, src->dim[0]._stride, base_type_src,
                               &dt_s); chk_err(ierr);
      }
      else
      {
        arr_bl = calloc(size, sizeof(int));
        arr_dsp_s = calloc(size, sizeof(int));

        dprint("Setting up strided vector index.\n",
               caf_this_image, caf_num_images, __FUNCTION__);
#define KINDCASE(kind, type)                                            \
case kind:                                                              \
  for (i = 0; i < size; ++i)                                            \
  {                                                                     \
    arr_dsp_s[i] = ((ptrdiff_t)                                         \
      ((type *) src_vector->u.v.vector)[i] - src->dim[0].lower_bound);  \
    arr_bl[i] = 1;                                                      \
  }                                                                     \
  break
        switch (src_vector->u.v.kind)
        {
          KINDCASE(1, int8_t);
          KINDCASE(2, int16_t);
          KINDCASE(4, int32_t);
          KINDCASE(8, int64_t);
#ifdef HAVE_GFC_INTEGER_16
          KINDCASE(16, __int128);
#endif
        default:
          caf_runtime_error(unreachable);
          return;
        }
#undef KINDCASE
        ierr = MPI_Type_indexed(size, arr_bl, arr_dsp_s, base_type_src, &dt_s);
        chk_err(ierr);
        free(arr_bl);
        free(arr_dsp_s);
      }
      ierr = MPI_Type_vector(size, 1, dest->dim[0]._stride, base_type_dst,
                             &dt_d); chk_err(ierr);
    }
    else
    {
      arr_bl = calloc(size, sizeof(int));
      arr_dsp_s = calloc(size, sizeof(int));
      arr_dsp_d = calloc(size, sizeof(int));

      for (i = 0; i < size; ++i)
      {
        arr_bl[i] = 1;
      }

      for (i = 0; i < size; ++i)
      {
        ptrdiff_t array_offset_sr = 0, extent = 1, tot_ext = 1;
        if (src_vector == NULL)
        {
          for (j = 0; j < src_rank - 1; ++j)
          {
            extent = src->dim[j]._ubound - src->dim[j].lower_bound + 1;
            array_offset_sr += ((i / tot_ext) % extent) * src->dim[j]._stride;
            tot_ext *= extent;
          }

          array_offset_sr += (i / tot_ext) * src->dim[src_rank - 1]._stride;
        }
        else
        {
#define KINDCASE(kind, type)                                          \
case kind:                                                            \
  array_offset_sr = ((ptrdiff_t)                                      \
    ((type *) src_vector->u.v.vector)[i] - src->dim[0].lower_bound);  \
  break
          switch (src_vector->u.v.kind)
          {
            KINDCASE(1, int8_t);
            KINDCASE(2, int16_t);
            KINDCASE(4, int32_t);
            KINDCASE(8, int64_t);
#ifdef HAVE_GFC_INTEGER_16
            KINDCASE(16, __int128);
#endif
            default:
              caf_runtime_error(unreachable);
              return;
          }
#undef KINDCASE
        }
        arr_dsp_s[i] = array_offset_sr;

        if (dst_rank != 0)
        {
          ptrdiff_t array_offset_dst = 0;
          extent = 1;
          tot_ext = 1;
          for (j = 0; j < dst_rank - 1; ++j)
          {
            extent = dest->dim[j]._ubound - dest->dim[j].lower_bound + 1;
            array_offset_dst += ((i / tot_ext) % extent) * dest->dim[j]._stride;
            tot_ext *= extent;
          }

          array_offset_dst += (i / tot_ext) * dest->dim[src_rank - 1]._stride;
          arr_dsp_d[i] = array_offset_dst;
        }
        else
          arr_dsp_d[i] = 0;
      }

      ierr = MPI_Type_indexed(size, arr_bl, arr_dsp_s, base_type_src, &dt_s);
      chk_err(ierr);
      ierr = MPI_Type_indexed(size, arr_bl, arr_dsp_d, base_type_dst, &dt_d);
      chk_err(ierr);

      free(arr_bl);
      free(arr_dsp_s);
      free(arr_dsp_d);
    }

    ierr = MPI_Type_commit(&dt_s); chk_err(ierr);
    ierr = MPI_Type_commit(&dt_d); chk_err(ierr);

    CAF_Win_lock(MPI_LOCK_SHARED, remote_image, *p);
    ierr = MPI_Get(dest->base_addr, 1, dt_d, remote_image, offset, 1, dt_s, *p);
    chk_err(ierr);
    CAF_Win_unlock(remote_image, *p);

#ifdef WITH_FAILED_IMAGES
    check_image_health(image_index, stat);

    if (!stat && ierr == STAT_FAILED_IMAGE)
      terminate_internal(ierr, 1);

    if (stat)
      *stat = ierr;
#else
    if (ierr != 0)
    {
      terminate_internal(ierr, 1);
      return;
    }
#endif
    ierr = MPI_Type_free(&dt_s); chk_err(ierr);
    ierr = MPI_Type_free(&dt_d); chk_err(ierr);
  }
#endif // STRIDED
  else
  {
    if (same_image && mrt)
    {
      if ((free_t_buff = (((t_buff = alloca(src_size * size))) == NULL)))
      {
        t_buff = malloc(src_size * size);
        if (t_buff == NULL)
          caf_runtime_error("Unable to allocate memory "
                            "for internal buffer in get().");
      }
    }
    else if (!same_type_and_kind && !same_image)
    {
      if ((free_t_buff = (((t_buff = alloca(src_size))) == NULL)))
      {
        t_buff = malloc(src_size);
        if (t_buff == NULL)
          caf_runtime_error("Unable to allocate memory "
                            "for internal buffer in get().");
      }
    }

    if (!same_image)
      CAF_Win_lock(MPI_LOCK_SHARED, remote_image, *p);
    for (i = 0; i < size; ++i)
    {
      ptrdiff_t array_offset_sr = 0, extent = 1, tot_ext = 1;
      if (src_vector == NULL)
      {
        for (j = 0; j < src_rank - 1; ++j)
        {
          extent = src->dim[j]._ubound - src->dim[j].lower_bound + 1;
          array_offset_sr += ((i / tot_ext) % extent) * src->dim[j]._stride;
          tot_ext *= extent;
        }

        array_offset_sr += (i / tot_ext) * src->dim[src_rank - 1]._stride;
      }
      else
      {
#define KINDCASE(kind, type)                                        \
case kind:                                                          \
  array_offset_sr = ((ptrdiff_t)                                    \
    ((type *)src_vector->u.v.vector)[i] - src->dim[0].lower_bound); \
  break
        switch (src_vector->u.v.kind)
        {
          KINDCASE(1, int8_t);
          KINDCASE(2, int16_t);
          KINDCASE(4, int32_t);
          KINDCASE(8, int64_t);
#ifdef HAVE_GFC_INTEGER_16
          KINDCASE(16, __int128);
#endif
          default:
            caf_runtime_error(unreachable);
            return;
        }
      }
      src_offset = array_offset_sr * src_size;
#undef KINDCASE

      void *dst;
      if (!same_image || !mrt)
      {
        if (dst_rank != 0)
        {
          ptrdiff_t array_offset_dst = 0;
          extent = 1;
          tot_ext = 1;
          for (j = 0; j < dst_rank - 1; ++j)
          {
            extent = dest->dim[j]._ubound - dest->dim[j].lower_bound + 1;
            array_offset_dst += ((i / tot_ext) % extent) * dest->dim[j]._stride;
            tot_ext *= extent;
          }

          array_offset_dst += (i / tot_ext) * dest->dim[dst_rank - 1]._stride;
          dst = (void *)((char *)dest->base_addr + array_offset_dst * dst_size);
        }
        else
          dst = dest->base_addr;
       }

      if (!same_image)
      {
        // Do the more likely first.
        dprint("kind(dst) = %d, el_sz(dst) = %zd, "
               "kind(src) = %d, el_sz(src) = %zd, lb(dst) = %zd.\n",
               dst_kind, dst_size, src_kind, src_size, src->dim[0].lower_bound);
        if (same_type_and_kind)
        {
          const size_t trans_size = (src_size < dst_size) ? src_size : dst_size;
          ierr = MPI_Get(dst, trans_size, MPI_BYTE, remote_image,
                         offset + src_offset, trans_size, MPI_BYTE, *p);
          chk_err(ierr);
          if (pad_str)
            memcpy((void *)((char *)dst + src_size), pad_str,
                   dst_size - src_size);
        }
        else if (dst_type == BT_CHARACTER)
        {
          ierr = MPI_Get(t_buff, src_size, MPI_BYTE, remote_image,
                         offset + src_offset, src_size, MPI_BYTE, *p);
          chk_err(ierr);
          copy_char_to_self(t_buff, src_type, src_size, src_kind,
                            dst, dst_type, dst_size, dst_kind, 1, true);
        }
        else
        {
          ierr = MPI_Get(t_buff, src_size, MPI_BYTE, remote_image,
                         offset + src_offset, src_size, MPI_BYTE, *p);
          chk_err(ierr);
          convert_type(dst, dst_type, dst_kind, t_buff,
                       src_type, src_kind, stat);
        }
      }
      else
      {
        if (!mrt)
        {
          dprint("strided same_image, no temp, for i = %zd, "
                 "src_offset = %zd, offset = %zd.\n",
                 i, src_offset, offset);
          if (same_type_and_kind)
            memmove(dst, src->base_addr + src_offset, src_size);
          else
            convert_type(dst, dst_type, dst_kind,
                         src->base_addr + src_offset, src_type, src_kind, stat);
        }
        else
        {
          dprint("strided same_image, *WITH* temp, for i = %zd.\n", i);
          if (same_type_and_kind)
            memmove(t_buff + i * dst_size,
                    src->base_addr + src_offset, src_size);
          else
            convert_type(t_buff + i * dst_size, dst_type, dst_kind,
                         src->base_addr + src_offset, src_type, src_kind, stat);
        }
      }

#ifndef WITH_FAILED_IMAGES
      if (ierr != 0)
      {
        caf_runtime_error("MPI Error: %d", ierr);
        return;
      }
#endif
    }
    if (!same_image)
      CAF_Win_unlock(remote_image, *p);

    if (same_image && mrt)
    {
      dprint("Same image temporary move.\n");
      memmove(dest->base_addr, t_buff, size * dst_size);
    }
  }

  /* Free memory, when not allocated on stack. */
  if (free_t_buff)
    free(t_buff);
  if (free_pad_str)
    free(pad_str);

#ifdef WITH_FAILED_IMAGES
  /* Catch failed images, when failed image support is active. */
  check_image_health(image_index , stat);
#endif

  if (ierr != MPI_SUCCESS)
  {
    int mpi_error;
    MPI_Error_class(ierr, &mpi_error);
    if (stat)
      *stat = mpi_error;
    else
    {
      int error_len = 2048 - 11;
      char error_str[error_len + 11];
      strcpy(error_str, "MPI-error: ");
      MPI_Error_string(mpi_error, &error_str[11], &error_len);
      error_stop_str(error_str, error_len + 11, false);
    }
  }
}


#ifdef GCC_GE_7
/* Get a chunk of data from one image to the current one, with type conversion.
 *
 * Copied from the gcc:libgfortran/caf/single.c. Can't say much about it. */
static void
get_data(void *ds, mpi_caf_token_t *token, MPI_Aint offset, int dst_type,
         int src_type, int dst_kind, int src_kind, size_t dst_size,
         size_t src_size, size_t num, int *stat, int image_index)
{
  size_t k;
  int ierr;
  MPI_Win win = (token == NULL) ? global_dynamic_win : token->memptr_win;
#ifdef EXTRA_DEBUG_OUTPUT
  if (token)
    dprint("%p = win(%d): %d -> offset: %zd of size %zd -> %zd, "
           "dst type %d(%d), src type %d(%d)\n",
           ds, win, image_index + 1, offset, src_size, dst_size,
           dst_type, dst_kind, src_type, src_kind);
  else
    dprint("%p = global_win(%d) offset: %zd (%zd) of size %zd -> %zd, "
           "dst type %d(%d), src type %d(%d)\n",
           ds, image_index + 1, offset, offset, src_size, dst_size,
           dst_type, dst_kind, src_type, src_kind);
#endif
  if (dst_type == src_type && dst_kind == src_kind)
  {
    size_t sz = ((dst_size > src_size) ? src_size : dst_size) * num;
    ierr = MPI_Get(ds, sz, MPI_BYTE, image_index, offset, sz, MPI_BYTE, win);
    chk_err(ierr);
    if ((dst_type == BT_CHARACTER || src_type == BT_CHARACTER)
        && dst_size > src_size)
      {
        if (dst_kind == 1)
        {
          memset((void*)(char*) ds + src_size, ' ', dst_size - src_size);
        }
        else /* dst_kind == 4. */
        {
          for (k = src_size / 4; k < dst_size / 4; k++)
            ((int32_t*) ds)[k] = (int32_t) ' ';
        }
      }
  }
  else if (dst_type == BT_CHARACTER && dst_kind == 1)
  {
    /* Get the required amount of memory on the stack. */
    void *srh = alloca(src_size);
    ierr = MPI_Get(srh, src_size, MPI_BYTE, image_index, offset, src_size,
                   MPI_BYTE, win); chk_err(ierr);
    /* Get of the data needs to be finished before converting the data. */
    ierr = MPI_Win_flush(image_index, win); chk_err(ierr);
    assign_char1_from_char4(dst_size, src_size, ds, srh);
  }
  else if (dst_type == BT_CHARACTER)
  {
    /* Get the required amount of memory on the stack. */
    void *srh = alloca(src_size);
    ierr = MPI_Get(srh, src_size, MPI_BYTE, image_index, offset, src_size,
                   MPI_BYTE, win); chk_err(ierr);
    /* Get of the data needs to be finished before converting the data. */
    ierr = MPI_Win_flush(image_index, win); chk_err(ierr);
    assign_char4_from_char1(dst_size, src_size, ds, srh);
  }
  else
  {
    /* Get the required amount of memory on the stack. */
    void *srh = alloca(src_size * num);
    dprint("type/kind convert %zd items: "
           "type %d(%d) -> type %d(%d), local buffer: %p\n",
           num, src_type, src_kind, dst_type, dst_kind, srh);
    ierr = MPI_Get(srh, src_size * num, MPI_BYTE, image_index, offset,
                   src_size * num, MPI_BYTE, win); chk_err(ierr);
    /* Get of the data needs to be finished before converting the data. */
    ierr = MPI_Win_flush(image_index, win); chk_err(ierr);
    dprint("srh[0] = %d, ierr = %d\n", (int)((char *)srh)[0], ierr);
    for (k = 0; k < num; ++k)
    {
      convert_type(ds, dst_type, dst_kind, srh, src_type, src_kind, stat);
      ds += dst_size;
      srh += src_size;
    }
  }
}


/* Compute the number of items referenced.
 *
 * Computes the number of items between lower bound (lb) and upper bound (ub)
 * with respect to the stride taking corner cases into account. */
#define COMPUTE_NUM_ITEMS(num, stride, lb, ub)                  \
do                                                              \
{                                                               \
  ptrdiff_t abs_stride = (stride) > 0 ? (stride) : -(stride);   \
  num = (stride) > 0 ? (ub) + 1 - (lb) : (lb) + 1 - (ub);       \
  if (num <= 0 || abs_stride < 1) return;                       \
  num = (abs_stride > 1) ? (1 + (num - 1) / abs_stride) : num;  \
} while (0)


/* Convenience macro to get the extent of a descriptor in a certain dimension
 *
 * Copied from gcc:libgfortran/libgfortran.h. */
#define GFC_DESCRIPTOR_EXTENT(desc,i) \
((desc)->dim[i]._ubound + 1 - (desc)->dim[i].lower_bound)


#define sizeof_desc_for_rank(rank) \
(sizeof(gfc_descriptor_t) + (rank) * sizeof(descriptor_dimension))

/* Define the descriptor of max rank.
 * 
 *  This typedef is made to allow storing a copy of a remote descriptor on the
 *  stack without having to care about the rank. */
typedef struct gfc_max_dim_descriptor_t
{
  gfc_descriptor_t base;
  descriptor_dimension dim[GFC_MAX_DIMENSIONS];
} gfc_max_dim_descriptor_t;

typedef struct gfc_dim1_descriptor_t
{
  gfc_descriptor_t base;
  descriptor_dimension dim[1];
} gfc_dim1_descriptor_t;

static void
get_for_ref(caf_reference_t *ref, size_t *i, size_t dst_index,
            mpi_caf_token_t *mpi_token, gfc_descriptor_t *dst,
            gfc_descriptor_t *src, void *ds, void *sr,
            ptrdiff_t sr_byte_offset, ptrdiff_t desc_byte_offset,
            int dst_kind, int src_kind, size_t dst_dim, size_t src_dim,
            size_t num, int *stat,
            int global_dynamic_win_rank, int memptr_win_rank,
            bool sr_global, /* access sr through global_dynamic_win */
            bool desc_global /* access desc through global_dynamic_win */
#ifdef GCC_GE_8
            , int src_type)
{
#else
            )
{
  int src_type = -1;
#endif
  ptrdiff_t extent_src = 1, array_offset_src = 0, stride_src;
  size_t next_dst_dim, ref_rank;
  gfc_max_dim_descriptor_t src_desc_data;
  int ierr;

  if (unlikely(ref == NULL))
  {
    /* May be we should issue an error here, because this case should not
     * occur. */
    return;
  }

  dprint("sr_offset = %zd, sr = %p, desc_offset = %zd, src = %p, "
         "sr_glb = %d, desc_glb = %d\n",
         sr_byte_offset, sr, desc_byte_offset, src, sr_global, desc_global);

  if (ref->next == NULL)
  {
    size_t dst_size = GFC_DESCRIPTOR_SIZE(dst);

    switch (ref->type)
    {
      case CAF_REF_COMPONENT:
        dprint("caf_offset = %zd\n", ref->u.c.caf_token_offset);
        if (ref->u.c.caf_token_offset > 0)
        {
          sr_byte_offset += ref->u.c.offset;
          if (sr_global)
          {
            ierr = MPI_Get(&sr, stdptr_size, MPI_BYTE, global_dynamic_win_rank,
                           MPI_Aint_add((MPI_Aint)sr, sr_byte_offset),
                           stdptr_size, MPI_BYTE, global_dynamic_win);
            chk_err(ierr);
            desc_global = true;
          }
          else
          {
            ierr = MPI_Get(&sr, stdptr_size, MPI_BYTE, memptr_win_rank,
                           sr_byte_offset, stdptr_size, MPI_BYTE,
                           mpi_token->memptr_win); chk_err(ierr);
            sr_global = true;
          }
          sr_byte_offset = 0;
        }
        else
          sr_byte_offset += ref->u.c.offset;
        if (sr_global)
        {
          get_data(ds, NULL, MPI_Aint_add((MPI_Aint)sr, sr_byte_offset),
                   GFC_DESCRIPTOR_TYPE(dst),
#ifdef GCC_GE_8
                   (src_type != -1) ? src_type : GFC_DESCRIPTOR_TYPE (dst),
#else
                   GFC_DESCRIPTOR_TYPE(dst),
#endif
                   dst_kind, src_kind, dst_size, ref->item_size, 1, stat,
                   global_dynamic_win_rank);
        }
        else
        {
          get_data(ds, mpi_token, sr_byte_offset, GFC_DESCRIPTOR_TYPE(dst),
#ifdef GCC_GE_8
                   src_type,
#else
                   GFC_DESCRIPTOR_TYPE(src),
#endif
                   dst_kind, src_kind, dst_size, ref->item_size, 1, stat,
                   memptr_win_rank);
        }
        ++(*i);
        return;
      case CAF_REF_STATIC_ARRAY:
        src_type = ref->u.a.static_array_type;
        /* Intentionally fall through. */
      case CAF_REF_ARRAY:
        if (ref->u.a.mode[src_dim] == CAF_ARR_REF_NONE)
        {
          if (sr_global)
          {
            get_data(ds + dst_index * dst_size, NULL,
                     MPI_Aint_add((MPI_Aint)sr, sr_byte_offset),
                     GFC_DESCRIPTOR_TYPE(dst),
#ifdef GCC_GE_8
                     (src_type != -1) ? src_type : GFC_DESCRIPTOR_TYPE (src),
#else
                     (src_type == -1) ? GFC_DESCRIPTOR_TYPE(src) : src_type,
#endif
                     dst_kind, src_kind, dst_size, ref->item_size, num,
                     stat, global_dynamic_win_rank);
          }
          else
          {
            get_data(ds + dst_index * dst_size, mpi_token,
                     sr_byte_offset, GFC_DESCRIPTOR_TYPE(dst),
#ifdef GCC_GE_8
                     (src_type != -1) ? src_type : GFC_DESCRIPTOR_TYPE (src),
#else
                     (src_type == -1) ? GFC_DESCRIPTOR_TYPE(src) : src_type,
#endif
                     dst_kind, src_kind, dst_size, ref->item_size, num,
                     stat, memptr_win_rank);
          }
          *i += num;
          return;
        }
        break;
      default:
        caf_runtime_error(unreachable);
    }
  }

  switch (ref->type)
  {
    case CAF_REF_COMPONENT:
      if (ref->u.c.caf_token_offset > 0)
      {
        sr_byte_offset += ref->u.c.offset;
        desc_byte_offset = sr_byte_offset;
        if (sr_global)
        {
          ierr = MPI_Get(&sr, stdptr_size, MPI_BYTE, global_dynamic_win_rank,
                         MPI_Aint_add((MPI_Aint)sr, sr_byte_offset),
                         stdptr_size, MPI_BYTE, global_dynamic_win);
          chk_err(ierr);
          desc_global = true;
        }
        else
        {
          ierr = MPI_Get(&sr, stdptr_size, MPI_BYTE, memptr_win_rank,
                         sr_byte_offset, stdptr_size, MPI_BYTE,
                         mpi_token->memptr_win); chk_err(ierr);
          sr_global = true;
        }
        sr_byte_offset = 0;
      }
      else
      {
        sr_byte_offset += ref->u.c.offset;
        desc_byte_offset += ref->u.c.offset;
      }
      get_for_ref(ref->next, i, dst_index, mpi_token, dst, NULL, ds,
                  sr, sr_byte_offset, desc_byte_offset, dst_kind, src_kind,
                  dst_dim, 0, 1, stat, global_dynamic_win_rank, memptr_win_rank,
                  sr_global, desc_global
#ifdef GCC_GE_8
                  , src_type
#endif
                  );
      return;
    case CAF_REF_ARRAY:
      if (ref->u.a.mode[src_dim] == CAF_ARR_REF_NONE)
      {
        get_for_ref(ref->next, i, dst_index, mpi_token, dst, src, ds, sr,
                    sr_byte_offset, desc_byte_offset, dst_kind, src_kind,
                    dst_dim, 0, 1, stat, global_dynamic_win_rank, memptr_win_rank,
                    sr_global, desc_global
#ifdef GCC_GE_8
                    , src_type
#endif
                    );
        return;
      }
      /* Only when on the left most index switch the data pointer to the
       * array's data pointer. */
      if (src_dim == 0)
      {
        if (sr_global)
        {
          for (ref_rank = 0; ref->u.a.mode[ref_rank] != CAF_ARR_REF_NONE;
               ++ref_rank) ;
          /* Get the remote descriptor. */
          if (desc_global)
          {
            ierr = MPI_Get(&src_desc_data, sizeof_desc_for_rank(ref_rank),
                           MPI_BYTE, global_dynamic_win_rank,
                           MPI_Aint_add((MPI_Aint)sr, desc_byte_offset),
                           sizeof_desc_for_rank(ref_rank), MPI_BYTE,
                           global_dynamic_win); chk_err(ierr);
          }
          else
          {
            ierr = MPI_Get(&src_desc_data,
                           sizeof_desc_for_rank(ref_rank), MPI_BYTE,
                           memptr_win_rank, desc_byte_offset,
                           sizeof_desc_for_rank(ref_rank),
                           MPI_BYTE, mpi_token->memptr_win); chk_err(ierr);
            desc_global = true;
          }
          src = (gfc_descriptor_t *)&src_desc_data;
        }
        else
          src = mpi_token->desc;
        sr_byte_offset = 0;
        desc_byte_offset = 0;
#ifdef EXTRA_DEBUG_OUTPUT
        dprint("remote desc rank: %zd (ref_rank: %zd)\n",
               GFC_DESCRIPTOR_RANK(src), ref_rank);
        for (int r = 0; r < GFC_DESCRIPTOR_RANK(src); ++r)
        {
          dprint("remote desc dim[%d] = (lb = %zd, ub = %zd, stride = %zd)\n",
                 r, src->dim[r].lower_bound, src->dim[r]._ubound,
                 src->dim[r]._stride);
        }
#endif
      }
      switch (ref->u.a.mode[src_dim])
      {
        case CAF_ARR_REF_VECTOR:
          extent_src = GFC_DESCRIPTOR_EXTENT(src, src_dim);
          array_offset_src = 0;
          for (size_t idx = 0; idx < ref->u.a.dim[src_dim].v.nvec; ++idx)
          {
#define KINDCASE(kind, type)                                        \
case kind:                                                          \
  array_offset_src = (((ptrdiff_t)                                  \
    ((type *)ref->u.a.dim[src_dim].v.vector)[idx])                  \
    - src->dim[src_dim].lower_bound  * src->dim[src_dim]._stride);  \
  break

            switch (ref->u.a.dim[src_dim].v.kind)
            {
              KINDCASE(1, int8_t);
              KINDCASE(2, int16_t);
              KINDCASE(4, int32_t);
              KINDCASE(8, int64_t);
#ifdef HAVE_GFC_INTEGER_16
              KINDCASE(16, __int128);
#endif
              default:
                caf_runtime_error(unreachable);
                return;
            }
#undef KINDCASE

            dprint("vector-index computed to: %zd\n", array_offset_src);
            get_for_ref(ref, i, dst_index, mpi_token, dst, src, ds, sr,
                        sr_byte_offset + array_offset_src * ref->item_size,
                        desc_byte_offset + array_offset_src * ref->item_size,
                        dst_kind, src_kind, dst_dim + 1, src_dim + 1,
                        1, stat, global_dynamic_win_rank, memptr_win_rank,
                        sr_global, desc_global
#ifdef GCC_GE_8
                        , src_type
#endif
                        );
            dst_index += dst->dim[dst_dim]._stride;
          }
          return;
        case CAF_ARR_REF_FULL:
          COMPUTE_NUM_ITEMS(extent_src,
                            ref->u.a.dim[src_dim].s.stride,
                            src->dim[src_dim].lower_bound,
                            src->dim[src_dim]._ubound);
          stride_src =
            src->dim[src_dim]._stride * ref->u.a.dim[src_dim].s.stride;
          array_offset_src = 0;
          for (ptrdiff_t idx = 0; idx < extent_src;
               ++idx, array_offset_src += stride_src)
            {
              get_for_ref(ref, i, dst_index, mpi_token, dst, src, ds, sr,
                          sr_byte_offset + array_offset_src * ref->item_size,
                          desc_byte_offset + array_offset_src * ref->item_size,
                          dst_kind, src_kind, dst_dim + 1, src_dim + 1,
                          1, stat, global_dynamic_win_rank, memptr_win_rank,
                          sr_global, desc_global
#ifdef GCC_GE_8
                          , src_type
#endif
                          );
              dst_index += dst->dim[dst_dim]._stride;
            }
          return;
        case CAF_ARR_REF_RANGE:
          COMPUTE_NUM_ITEMS(extent_src,
                            ref->u.a.dim[src_dim].s.stride,
                            ref->u.a.dim[src_dim].s.start,
                            ref->u.a.dim[src_dim].s.end);
          array_offset_src = 
            (ref->u.a.dim[src_dim].s.start - src->dim[src_dim].lower_bound)
            * src->dim[src_dim]._stride;
          stride_src =
            src->dim[src_dim]._stride * ref->u.a.dim[src_dim].s.stride;
          /* Increase the dst_dim only, when the src_extent is greater than one
           * or src and dst extent are both one. Don't increase when the scalar
           * source is not present in the dst. */
          next_dst_dim = (
            (extent_src > 1) ||
            (GFC_DESCRIPTOR_EXTENT(dst, dst_dim) == 1 && extent_src == 1)
          ) ? (dst_dim + 1) : dst_dim;
          for (ptrdiff_t idx = 0; idx < extent_src; ++idx)
          {
            get_for_ref(ref, i, dst_index, mpi_token, dst, src, ds, sr,
                        sr_byte_offset + array_offset_src * ref->item_size,
                        desc_byte_offset + array_offset_src * ref->item_size,
                        dst_kind, src_kind, next_dst_dim, src_dim + 1,
                        1, stat, global_dynamic_win_rank, memptr_win_rank,
                        sr_global, desc_global
#ifdef GCC_GE_8
                        , src_type
#endif
                        );
            dst_index += dst->dim[dst_dim]._stride;
            array_offset_src += stride_src;
          }
          return;
        case CAF_ARR_REF_SINGLE:
          array_offset_src =
            (ref->u.a.dim[src_dim].s.start - src->dim[src_dim].lower_bound)
            * src->dim[src_dim]._stride;
          get_for_ref(ref, i, dst_index, mpi_token, dst, src, ds, sr,
                      sr_byte_offset + array_offset_src * ref->item_size,
                      desc_byte_offset + array_offset_src * ref->item_size,
                      dst_kind, src_kind, dst_dim, src_dim + 1, 1,
                      stat, global_dynamic_win_rank, memptr_win_rank,
                      sr_global, desc_global
#ifdef GCC_GE_8
                      , src_type
#endif
                      );
          return;
        case CAF_ARR_REF_OPEN_END:
          COMPUTE_NUM_ITEMS(extent_src,
                            ref->u.a.dim[src_dim].s.stride,
                            ref->u.a.dim[src_dim].s.start,
                            src->dim[src_dim]._ubound);
          stride_src =
            src->dim[src_dim]._stride * ref->u.a.dim[src_dim].s.stride;
          array_offset_src = (ref->u.a.dim[src_dim].s.start
                              - src->dim[src_dim].lower_bound)
                             * src->dim[src_dim]._stride;
          for (ptrdiff_t idx = 0; idx < extent_src; ++idx)
          {
            get_for_ref(ref, i, dst_index, mpi_token, dst, src, ds, sr,
                        sr_byte_offset + array_offset_src * ref->item_size,
                        desc_byte_offset + array_offset_src * ref->item_size,
                        dst_kind, src_kind, dst_dim + 1, src_dim + 1,
                        1, stat, global_dynamic_win_rank, memptr_win_rank,
                        sr_global, desc_global
#ifdef GCC_GE_8
                        , src_type
#endif
                        );
            dst_index += dst->dim[dst_dim]._stride;
            array_offset_src += stride_src;
          }
          return;
        case CAF_ARR_REF_OPEN_START:
          COMPUTE_NUM_ITEMS(extent_src,
                            ref->u.a.dim[src_dim].s.stride,
                            src->dim[src_dim].lower_bound,
                            ref->u.a.dim[src_dim].s.end);
          stride_src =
            src->dim[src_dim]._stride * ref->u.a.dim[src_dim].s.stride;
          array_offset_src = 0;
          for (ptrdiff_t idx = 0; idx < extent_src; ++idx)
          {
            get_for_ref(ref, i, dst_index, mpi_token, dst, src, ds, sr,
                        sr_byte_offset + array_offset_src * ref->item_size,
                        desc_byte_offset + array_offset_src * ref->item_size,
                        dst_kind, src_kind, dst_dim + 1, src_dim + 1,
                        1, stat, global_dynamic_win_rank, memptr_win_rank,
                        sr_global, desc_global
#ifdef GCC_GE_8
                        , src_type
#endif
                        );
            dst_index += dst->dim[dst_dim]._stride;
            array_offset_src += stride_src;
          }
          return;
        default:
          caf_runtime_error(unreachable);
      }
      return;
    case CAF_REF_STATIC_ARRAY:
      if (ref->u.a.mode[src_dim] == CAF_ARR_REF_NONE)
      {
        get_for_ref(ref->next, i, dst_index, mpi_token, dst, NULL, ds, sr,
                    sr_byte_offset, desc_byte_offset, dst_kind, src_kind,
                    dst_dim, 0, 1, stat, global_dynamic_win_rank, memptr_win_rank,
                    sr_global, desc_global
#ifdef GCC_GE_8
                    , src_type
#endif
                    );
        return;
      }
      switch (ref->u.a.mode[src_dim])
      {
      case CAF_ARR_REF_VECTOR:
        array_offset_src = 0;
        for (size_t idx = 0; idx < ref->u.a.dim[src_dim].v.nvec; ++idx)
        {
#define KINDCASE(kind, type)                                        \
case kind:                                                          \
  array_offset_src = ((type *)ref->u.a.dim[src_dim].v.vector)[idx]; \
  break

          switch (ref->u.a.dim[src_dim].v.kind)
          {
            KINDCASE(1, int8_t);
            KINDCASE(2, int16_t);
            KINDCASE(4, int32_t);
            KINDCASE(8, int64_t);
#ifdef HAVE_GFC_INTEGER_16
            KINDCASE(16, __int128);
#endif
            default:
              caf_runtime_error(unreachable);
              return;
          }
#undef KINDCASE

          get_for_ref(ref, i, dst_index, mpi_token, dst, NULL, ds, sr,
                      sr_byte_offset + array_offset_src * ref->item_size,
                      desc_byte_offset + array_offset_src * ref->item_size,
                      dst_kind, src_kind, dst_dim + 1, src_dim + 1,
                      1, stat, global_dynamic_win_rank, memptr_win_rank,
                      sr_global, desc_global
#ifdef GCC_GE_8
                      , src_type
#endif
                      );
          dst_index += dst->dim[dst_dim]._stride;
        }
        return;
      case CAF_ARR_REF_FULL:
        for (array_offset_src = 0 ;
             array_offset_src <= ref->u.a.dim[src_dim].s.end;
             array_offset_src += ref->u.a.dim[src_dim].s.stride)
        {
          get_for_ref(ref, i, dst_index, mpi_token, dst, NULL, ds, sr,
                      sr_byte_offset + array_offset_src * ref->item_size,
                      desc_byte_offset + array_offset_src * ref->item_size,
                      dst_kind, src_kind, dst_dim + 1, src_dim + 1,
                      1, stat, global_dynamic_win_rank, memptr_win_rank,
                      sr_global, desc_global
#ifdef GCC_GE_8
                       , src_type
#endif
                       );
          dst_index += dst->dim[dst_dim]._stride;
        }
        return;
      case CAF_ARR_REF_RANGE:
        COMPUTE_NUM_ITEMS(extent_src,
                          ref->u.a.dim[src_dim].s.stride,
                          ref->u.a.dim[src_dim].s.start,
                          ref->u.a.dim[src_dim].s.end);
        array_offset_src = ref->u.a.dim[src_dim].s.start;
        for (ptrdiff_t idx = 0; idx < extent_src; ++idx)
        {
          get_for_ref(ref, i, dst_index, mpi_token, dst, NULL, ds, sr,
                      sr_byte_offset + array_offset_src * ref->item_size,
                      desc_byte_offset + array_offset_src * ref->item_size,
                      dst_kind, src_kind, dst_dim + 1, src_dim + 1,
                      1, stat, global_dynamic_win_rank, memptr_win_rank,
                      sr_global, desc_global
#ifdef GCC_GE_8
                      , src_type
#endif
                      );
          dst_index += dst->dim[dst_dim]._stride;
          array_offset_src += ref->u.a.dim[src_dim].s.stride;
        }
        return;
      case CAF_ARR_REF_SINGLE:
        array_offset_src = ref->u.a.dim[src_dim].s.start;
        get_for_ref(ref, i, dst_index, mpi_token, dst, NULL, ds, sr,
                    sr_byte_offset + array_offset_src * ref->item_size,
                    desc_byte_offset + array_offset_src * ref->item_size,
                    dst_kind, src_kind, dst_dim, src_dim + 1, 1,
                    stat, global_dynamic_win_rank, memptr_win_rank,
                    sr_global, desc_global
#ifdef GCC_GE_8
                    , src_type
#endif
                    );
        return;
        /* The OPEN_* are mapped to a RANGE and therefore can not occur. */
      case CAF_ARR_REF_OPEN_END:
      case CAF_ARR_REF_OPEN_START:
      default:
        caf_runtime_error(unreachable);
      }
      return;
    default:
      caf_runtime_error(unreachable);
  }
}

void
PREFIX(get_by_ref) (caf_token_t token, int image_index,
                    gfc_descriptor_t *dst, caf_reference_t *refs,
                    int dst_kind, int src_kind,
                    bool may_require_tmp __attribute__((unused)),
                    bool dst_reallocatable, int *stat
#ifdef GCC_GE_8
                    , int src_type
#endif
                    )
{
  const char vecrefunknownkind[] =
    "libcaf_mpi::caf_get_by_ref(): unknown kind in vector-ref.\n";
  const char unknownreftype[] =
    "libcaf_mpi::caf_get_by_ref(): unknown reference type.\n";
  const char unknownarrreftype[] =
    "libcaf_mpi::caf_get_by_ref(): unknown array reference type.\n";
  const char rankoutofrange[] =
    "libcaf_mpi::caf_get_by_ref(): rank out of range.\n";
  const char extentoutofrange[] =
    "libcaf_mpi::caf_get_by_ref(): extent out of range.\n";
  const char cannotallocdst[] =
    "libcaf_mpi::caf_get_by_ref(): can not allocate %d bytes of memory.\n";
  const char nonallocextentmismatch[] =
    "libcaf_mpi::caf_get_by_ref(): extent of non-allocatable arrays "
    "mismatch (%lu != %lu).\n";
  const char doublearrayref[] =
    "libcaf_mpi::caf_get_by_ref(): two or more array part references "
    "are not supported.\n";
  size_t size, i, ref_rank, dst_index, src_size;
  int ierr, dst_rank = GFC_DESCRIPTOR_RANK(dst), dst_cur_dim = 0;
  mpi_caf_token_t *mpi_token = (mpi_caf_token_t *) token;
  void *remote_memptr = mpi_token->memptr, *remote_base_memptr = NULL;
  gfc_max_dim_descriptor_t src_desc;
  gfc_descriptor_t *src = (gfc_descriptor_t *)&src_desc;
  caf_reference_t *riter = refs;
  long delta;
  ptrdiff_t data_offset = 0, desc_offset = 0;
  /* Reallocation of dst.data is needed (e.g., array to small). */
  bool realloc_needed;
  /* Reallocation of dst.data is required, because data is not alloced at
   * all. */
  bool realloc_required, extent_mismatch = false;
  /* Set when the first non-scalar array reference is encountered. */
  bool in_array_ref = false, array_extent_fixed = false;
  /* Set when remote data is to be accessed through the 
   * global dynamic window. */
  bool access_data_through_global_win = false;
  /* Set when the remote descriptor is to accessed through the global window. */
  bool access_desc_through_global_win = false;
  caf_array_ref_t array_ref;

  realloc_needed = realloc_required = dst->base_addr == NULL;

  if (stat)
    *stat = 0;

  MPI_Group current_team_group, win_group;
  int global_dynamic_win_rank, memptr_win_rank;
  ierr = MPI_Comm_group(CAF_COMM_WORLD, &current_team_group); chk_err(ierr);
  ierr = MPI_Win_get_group(global_dynamic_win, &win_group); chk_err(ierr);
  ierr = MPI_Group_translate_ranks(current_team_group, 1,
                                   (int[]){image_index - 1}, win_group,
                                   &global_dynamic_win_rank); chk_err(ierr);
  ierr = MPI_Group_free(&win_group); chk_err(ierr);
  ierr = MPI_Win_get_group(mpi_token->memptr_win, &win_group); chk_err(ierr);
  ierr = MPI_Group_translate_ranks(current_team_group, 1,
                                   (int[]){image_index - 1}, win_group,
                                   &memptr_win_rank); chk_err(ierr);
  ierr = MPI_Group_free(&current_team_group); chk_err(ierr);
  ierr = MPI_Group_free(&win_group); chk_err(ierr);

  check_image_health(global_dynamic_win_rank, stat);

  dprint("Entering get_by_ref(may_require_tmp = %d).\n", may_require_tmp);

  /* Compute the size of the result.  In the beginning size just counts the
   * number of elements. */
  size = 1;
  /* Shared lock both windows to prevent bother in the sub-routines. */
  CAF_Win_lock(MPI_LOCK_SHARED, global_dynamic_win_rank, global_dynamic_win);
  CAF_Win_lock(MPI_LOCK_SHARED, memptr_win_rank, mpi_token->memptr_win);
  while (riter)
  {
    dprint("offset = %zd, remote_mem = %p, access_data(global_win) = %d\n",
           data_offset, remote_memptr, access_data_through_global_win);
    switch (riter->type)
    {
      case CAF_REF_COMPONENT:
        if (riter->u.c.caf_token_offset > 0)
        {
          if (access_data_through_global_win)
          {
            data_offset += riter->u.c.offset;
            remote_base_memptr = remote_memptr;
            ierr = MPI_Get(&remote_memptr, stdptr_size, MPI_BYTE, global_dynamic_win_rank,
                           MPI_Aint_add((MPI_Aint)remote_memptr, data_offset),
                           stdptr_size, MPI_BYTE, global_dynamic_win);
            chk_err(ierr);
            /* On the second indirection access also the remote descriptor
             * using the global window. */
            access_desc_through_global_win = true;
          }
          else
          {
            data_offset += riter->u.c.offset;
            ierr = MPI_Get(&remote_memptr, stdptr_size, MPI_BYTE, memptr_win_rank,
                           data_offset, stdptr_size, MPI_BYTE,
                           mpi_token->memptr_win); chk_err(ierr);
            dprint("get(custom_token %d), offset = %zd, res. remote_mem = %p\n",
                   mpi_token->memptr_win, data_offset, remote_memptr);
            /* All future access is through the global dynamic window. */
            access_data_through_global_win = true;
          }
          desc_offset = data_offset;
          data_offset = 0;
        }
        else
        {
          data_offset += riter->u.c.offset;
          desc_offset += riter->u.c.offset;
        }
        break;
      case CAF_REF_ARRAY:
        /* When there has been no CAF_REF_COMP before hand, then the
         * descriptor is stored in the token and the extends are the same on
         * all images, which is taken care of in the else part. */
        if (access_data_through_global_win)
        {
          for (ref_rank = 0; riter->u.a.mode[ref_rank] != CAF_ARR_REF_NONE;
               ++ref_rank) ;
          /* Get the remote descriptor and use the stack to store it. Note,
           * src may be pointing to mpi_token->desc therefore it needs to be
           * reset here. */
          src = (gfc_descriptor_t *)&src_desc;
          if (access_desc_through_global_win)
          {
            dprint("remote desc fetch from %p, offset = %zd\n",
                   remote_base_memptr, desc_offset);
            MPI_Get(src, sizeof_desc_for_rank(ref_rank), MPI_BYTE, global_dynamic_win_rank,
                    MPI_Aint_add((MPI_Aint)remote_base_memptr, desc_offset),
                    sizeof_desc_for_rank(ref_rank), MPI_BYTE,
                    global_dynamic_win);
          }
          else
          {
            dprint("remote desc fetch from win %d, offset = %zd\n",
                   mpi_token->memptr_win, desc_offset);
            MPI_Get(src, sizeof_desc_for_rank(ref_rank), MPI_BYTE, memptr_win_rank,
                    desc_offset, sizeof_desc_for_rank(ref_rank), MPI_BYTE,
                    mpi_token->memptr_win);
            access_desc_through_global_win = true;
          }
        }
        else
          src = mpi_token->desc;

#ifdef EXTRA_DEBUG_OUTPUT
        dprint("remote desc rank: %zd (ref_rank: %zd)\n",
               GFC_DESCRIPTOR_RANK(src), ref_rank);
        for (i = 0; i < GFC_DESCRIPTOR_RANK(src); ++i)
        {
          dprint("remote desc dim[%zd] = (lb = %zd, ub = %zd, stride = %zd)\n",
                 i, src->dim[i].lower_bound, src->dim[i]._ubound,
                 src->dim[i]._stride);
        }
#endif
        for (i = 0; riter->u.a.mode[i] != CAF_ARR_REF_NONE; ++i)
        {
          array_ref = riter->u.a.mode[i];
          dprint("i = %zd, array_ref = %s\n", i, caf_array_ref_str[array_ref]);
          switch (array_ref)
          {
            case CAF_ARR_REF_VECTOR:
              delta = riter->u.a.dim[i].v.nvec;
#define KINDCASE(kind, type)                                            \
case kind:                                                              \
  remote_memptr += (((ptrdiff_t)                                        \
    ((type *)riter->u.a.dim[i].v.vector)[0]) - src->dim[i].lower_bound) \
    * src->dim[i]._stride * riter->item_size;                           \
  break

              switch (riter->u.a.dim[i].v.kind)
              {
                KINDCASE(1, int8_t);
                KINDCASE(2, int16_t);
                KINDCASE(4, int32_t);
                KINDCASE(8, int64_t);
#ifdef HAVE_GFC_INTEGER_16
                KINDCASE(16, __int128);
#endif
                default:
                  caf_runtime_error(vecrefunknownkind, stat, NULL, 0);
                  return;
              }
#undef KINDCASE
              break;
            case CAF_ARR_REF_FULL:
              COMPUTE_NUM_ITEMS(delta,
                                riter->u.a.dim[i].s.stride,
                                src->dim[i].lower_bound,
                                src->dim[i]._ubound);
              /* The memptr stays unchanged when ref'ing the first element
               * in a dimension. */
              break;
            case CAF_ARR_REF_RANGE:
              COMPUTE_NUM_ITEMS(delta,
                                riter->u.a.dim[i].s.stride,
                                riter->u.a.dim[i].s.start,
                                riter->u.a.dim[i].s.end);
              remote_memptr += 
                (riter->u.a.dim[i].s.start - src->dim[i].lower_bound)
                * src->dim[i]._stride * riter->item_size;
              break;
            case CAF_ARR_REF_SINGLE:
              delta = 1;
              remote_memptr +=
                (riter->u.a.dim[i].s.start - src->dim[i].lower_bound)
                * src->dim[i]._stride * riter->item_size;
              break;
            case CAF_ARR_REF_OPEN_END:
              COMPUTE_NUM_ITEMS(delta,
                                riter->u.a.dim[i].s.stride,
                                riter->u.a.dim[i].s.start,
                                src->dim[i]._ubound);
              remote_memptr +=
                (riter->u.a.dim[i].s.start - src->dim[i].lower_bound)
                * src->dim[i]._stride * riter->item_size;
              break;
            case CAF_ARR_REF_OPEN_START:
              COMPUTE_NUM_ITEMS(delta,
                                riter->u.a.dim[i].s.stride,
                                src->dim[i].lower_bound,
                                riter->u.a.dim[i].s.end);
              /* The memptr stays unchanged when ref'ing the first element
               * in a dimension. */
              break;
            default:
              caf_runtime_error(unknownarrreftype, stat, NULL, 0);
              return;
          }
          if (delta <= 0)
            return;
          /* Check the various properties of the destination array.
           * Is an array expected and present? */
          if (delta > 1 && dst_rank == 0)
          {
            /* No, an array is required, but not provided. */
            caf_runtime_error(extentoutofrange, stat, NULL, 0);
            return;
          }
          /* When dst is an array. */
          if (dst_rank > 0)
          {
            /* Check that dst_cur_dim is valid for dst. Can be superceeded
             * only by scalar data. */
            if (dst_cur_dim >= dst_rank && delta != 1)
            {
              caf_runtime_error(rankoutofrange, stat, NULL, 0);
              return;
            }
            /* Do further checks, when the source is not scalar. */
            else if (delta != 1)
            {
              /* Check that the extent is not scalar and we are not in an array
               * ref for the dst side. */
              if (!in_array_ref)
              {
                /* Check that this is the non-scalar extent. */
                if (!array_extent_fixed)
                {
                  /* In an array extent now. */
                  in_array_ref = true;
                  /* Check that we haven't skipped any scalar  dimensions yet
                   * and that the dst is compatible. */
                  if (i > 0 && dst_rank == GFC_DESCRIPTOR_RANK(src))
                  {
                    if (dst_reallocatable)
                    {
                      /* Dst is reallocatable, which means that the bounds are
                       * not set. Set them. */
                      for (dst_cur_dim = 0; dst_cur_dim < (int)i; ++dst_cur_dim)
                      {
                        dst->dim[dst_cur_dim].lower_bound = 1;
                        dst->dim[dst_cur_dim]._ubound = 1;
                        dst->dim[dst_cur_dim]._stride = 1;
                      }
                    }
                    else
                      dst_cur_dim = i;
                  }
                  /* Else press thumbs, that there are enough dimensional refs
                   * to come. Checked below. */
                }
                else
                {
                  caf_runtime_error(doublearrayref, stat, NULL, 0);
                  return;
                }
              }
              /* When the realloc is required, then no extent may have
               * been set. */
              extent_mismatch = realloc_required ||
                GFC_DESCRIPTOR_EXTENT(dst, dst_cur_dim) != delta;
              /* When it already known, that a realloc is needed or the extent
               * does not match the needed one. */
              if (realloc_required || realloc_needed || extent_mismatch)
              {
                /* Check whether dst is reallocatable. */
                if (unlikely(!dst_reallocatable))
                {
                  caf_runtime_error(nonallocextentmismatch, stat,
                                    NULL, 0, delta,
                                    GFC_DESCRIPTOR_EXTENT(dst, dst_cur_dim));
                  return;
                }
                /* Only report an error, when the extent needs to be modified,
                 * which is not allowed. */
                else if (!dst_reallocatable && extent_mismatch)
                {
                  caf_runtime_error(extentoutofrange, stat, NULL, 0);
                  return;
                }
                realloc_needed = true;
              }
              /* Only change the extent when it does not match.  This is to
               * prevent resetting given array bounds. */
              if (extent_mismatch)
              {
                dst->dim[dst_cur_dim].lower_bound = 1;
                dst->dim[dst_cur_dim]._ubound = delta;
                dst->dim[dst_cur_dim]._stride = size;
                if (realloc_required)
                  dst->offset = -1;
              }
            }

            /* Only increase the dim counter, when in an array ref */
            if (in_array_ref && dst_cur_dim < dst_rank)
            {
              /* Mode != CAF_ARR_REF_SINGLE(delta == 1), and no rank
               * reduction */
              if (!(delta == 1 && dst_rank != GFC_DESCRIPTOR_RANK(src)))
                ++dst_cur_dim;
            }
          }
          size *= (ptrdiff_t)delta;
        }
        if (in_array_ref)
        {
          array_extent_fixed = true;
          in_array_ref = false;
        }
        break;
      case CAF_REF_STATIC_ARRAY:
        for (i = 0; riter->u.a.mode[i] != CAF_ARR_REF_NONE; ++i)
        {
          array_ref = riter->u.a.mode[i];
          dprint("i = %zd, array_ref = %s\n", i, caf_array_ref_str[array_ref]);
          switch (array_ref)
          {
            case CAF_ARR_REF_VECTOR:
              delta = riter->u.a.dim[i].v.nvec;
#define KINDCASE(kind, type)                                    \
case kind:                                                      \
  remote_memptr +=                                              \
    ((type *)riter->u.a.dim[i].v.vector)[0] * riter->item_size; \
  break

              switch (riter->u.a.dim[i].v.kind)
              {
                KINDCASE(1, int8_t);
                KINDCASE(2, int16_t);
                KINDCASE(4, int32_t);
                KINDCASE(8, int64_t);
#ifdef HAVE_GFC_INTEGER_16
                KINDCASE(16, __int128);
#endif
                default:
                  caf_runtime_error(vecrefunknownkind, stat, NULL, 0);
                  return;
              }
#undef KINDCASE
              break;
            case CAF_ARR_REF_FULL:
              delta = riter->u.a.dim[i].s.end / riter->u.a.dim[i].s.stride
                      + 1;
              /* The memptr stays unchanged when ref'ing the first element in a
               * dimension. */
              break;
            case CAF_ARR_REF_RANGE:
              COMPUTE_NUM_ITEMS(delta,
                                riter->u.a.dim[i].s.stride,
                                riter->u.a.dim[i].s.start,
                                riter->u.a.dim[i].s.end);
              remote_memptr += riter->u.a.dim[i].s.start
                               * riter->u.a.dim[i].s.stride
                               * riter->item_size;
              break;
            case CAF_ARR_REF_SINGLE:
              delta = 1;
              remote_memptr += riter->u.a.dim[i].s.start
                               * riter->u.a.dim[i].s.stride
                               * riter->item_size;
              break;
            case CAF_ARR_REF_OPEN_END:
              /* This and OPEN_START are mapped to a RANGE and therefore can
               * not occur here. */
            case CAF_ARR_REF_OPEN_START:
            default:
              caf_runtime_error(unknownarrreftype, stat, NULL, 0);
              return;
          }
          if (delta <= 0)
            return;
          /* Check the various properties of the destination array.
           * Is an array expected and present? */
          if (delta > 1 && dst_rank == 0)
          {
            /* No, an array is required, but not provided. */
            caf_runtime_error(extentoutofrange, stat, NULL, 0);
            return;
          }
          /* When dst is an array. */
          if (dst_rank > 0)
          {
            /* Check that dst_cur_dim is valid for dst.  Can be superceeded
             * only by scalar data. */
            if (dst_cur_dim >= dst_rank && delta != 1)
            {
              caf_runtime_error(rankoutofrange, stat, NULL, 0);
              return;
            }
            /* Do further checks, when the source is not scalar. */
            else if (delta != 1)
            {
              /* Check that the extent is not scalar and we are not in an array
               * ref for the dst side. */
              if (!in_array_ref)
              {
                /* Check that this is the non-scalar extent. */
                if (!array_extent_fixed)
                {
                  /* In an array extent now. */
                  in_array_ref = true;
                  /* The dst is not reallocatable, so nothing more to do,
                   * then correct the dim counter. */
                  dst_cur_dim = i;
                }
                else
                {
                  caf_runtime_error(doublearrayref, stat, NULL, 0);
                  return;
                }
              }
              /* When the realloc is required, then no extent may have
               * been set. */
              extent_mismatch = realloc_required ||
                GFC_DESCRIPTOR_EXTENT(dst, dst_cur_dim) != delta;
              /* When it is already known, that a realloc is needed or
               * the extent does not match the needed one. */
              if (realloc_required || realloc_needed || extent_mismatch)
              {
                /* Check whether dst is reallocatable. */
                if (unlikely(!dst_reallocatable))
                {
                  caf_runtime_error(nonallocextentmismatch, stat,
                                    NULL, 0, delta,
                                    GFC_DESCRIPTOR_EXTENT(dst, dst_cur_dim));
                  return;
                }
                /* Only report an error, when the extent needs to be modified,
                 * which is not allowed. */
                else if (!dst_reallocatable && extent_mismatch)
                {
                  caf_runtime_error(extentoutofrange, stat, NULL, 0);
                  return;
                }
                realloc_needed = true;
              }
              /* Only change the extent when it does not match.  This is to
               * prevent resetting given array bounds. */
              if (extent_mismatch)
              {
                dst->dim[dst_cur_dim].lower_bound = 1;
                dst->dim[dst_cur_dim]._ubound = delta;
                dst->dim[dst_cur_dim]._stride = size;
                if (realloc_required)
                  dst->offset = -1;
              }
            }
            /* Only increase the dim counter, when in an array ref */
            if (in_array_ref && dst_cur_dim < dst_rank)
            {
              /* Mode != CAF_ARR_REF_SINGLE(delta == 1), and no rank
               * reduction */
              if (!(delta == 1 && dst_rank != GFC_DESCRIPTOR_RANK(src)))
                ++dst_cur_dim;
            }
          }
          size *= (ptrdiff_t)delta;
        }
        if (in_array_ref)
        {
          array_extent_fixed = true;
          in_array_ref = false;
        }
        break;
      default:
        caf_runtime_error(unknownreftype, stat, NULL, 0);
        return;
    }
    src_size = riter->item_size;
    riter = riter->next;
  }
  if (size == 0 || src_size == 0)
    return;
  /* Postcondition:
   * - size contains the number of elements to store in the destination array,
   * - src_size gives the size in bytes of each item in the destination array.
   */

  if (realloc_needed)
  {
    if (!array_extent_fixed)
    {
      /* This can happen only, when the result is scalar. */
      for (dst_cur_dim = 0; dst_cur_dim < dst_rank; ++dst_cur_dim)
      {
        dst->dim[dst_cur_dim].lower_bound = 1;
        dst->dim[dst_cur_dim]._ubound = 1;
        dst->dim[dst_cur_dim]._stride = 1;
      }
    }
    dst->base_addr = malloc(size * GFC_DESCRIPTOR_SIZE(dst));
    if (unlikely(dst->base_addr == NULL))
    {
      caf_runtime_error(cannotallocdst, stat, size * GFC_DESCRIPTOR_SIZE(dst));
      return;
    }
  }

  /* Reset the token. */
  mpi_token = (mpi_caf_token_t *) token;
  remote_memptr = mpi_token->memptr;
  dst_index = 0;
#ifdef EXTRA_DEBUG_OUTPUT
  dprint("dst_rank: %zd\n", GFC_DESCRIPTOR_RANK(dst));
  for (i = 0; i < GFC_DESCRIPTOR_RANK(dst); ++i)
  {
    dprint("dst_dim[%zd] = (%zd, %zd)\n",
           i, dst->dim[i].lower_bound, dst->dim[i]._ubound);
  }
#endif
  i = 0;
  dprint("get_by_ref() calling get_for_ref.\n");
  get_for_ref(refs, &i, dst_index, mpi_token, dst, mpi_token->desc,
              dst->base_addr, remote_memptr, 0, 0, dst_kind, src_kind, 0, 0,
              1, stat, global_dynamic_win_rank, memptr_win_rank, false, false
#ifdef GCC_GE_8
               , src_type
#endif
               );
  CAF_Win_unlock(global_dynamic_win_rank, global_dynamic_win);
  CAF_Win_unlock(memptr_win_rank, mpi_token->memptr_win);
}

static void
put_data(mpi_caf_token_t *token, MPI_Aint offset, void *sr, int dst_type,
         int src_type, int dst_kind, int src_kind, size_t dst_size,
         size_t src_size, size_t num, int *stat, int image_index)
{
  size_t k;
  int ierr;
  MPI_Win win = (token == NULL) ? global_dynamic_win : token->memptr_win;
#ifdef EXTRA_DEBUG_OUTPUT
  if (token)
    dprint("(win: %d, image: %d, offset: %zd) <- %p, "
           "num: %zd, size %zd -> %zd, dst type %d(%d), src type %d(%d)\n",
           win, image_index + 1, offset, sr, num, src_size, dst_size,
           dst_type, dst_kind, src_type, src_kind);
  else
    dprint("(global_win: %x, image: %d, offset: %zd (%zd)) <- %p, "
           "num: %zd, size %zd -> %zd, dst type %d(%d), src type %d(%d)\n",
           win, image_index + 1, offset, offset, sr, num, src_size,
           dst_size, dst_type, dst_kind, src_type, src_kind);
#endif
  if (dst_type == src_type && dst_kind == src_kind)
  {
    size_t sz = (dst_size > src_size ? src_size : dst_size) * num;
    ierr = MPI_Put(sr, sz, MPI_BYTE, image_index, offset, sz, MPI_BYTE, win);
    chk_err(ierr);
    dprint("sr[] = %d, num = %zd, num bytes = %zd\n",
           (int)((char*)sr)[0], num, sz);
    if ((dst_type == BT_CHARACTER || src_type == BT_CHARACTER)
        && dst_size > src_size)
    {
      const size_t trans_size = dst_size / dst_kind - src_size / src_kind;
      void *pad = alloca(trans_size * dst_kind);
      if (dst_kind == 1)
      {
        memset((void*)(char*) pad, ' ', trans_size);
      }
      else /* dst_kind == 4. */
      {
        for (k = 0; k < trans_size; ++k)
        {
          ((int32_t*) pad)[k] = (int32_t) ' ';
        }
      }
      ierr = MPI_Put(pad, trans_size * dst_kind, MPI_BYTE, image_index,
                     offset + (src_size / src_kind) * dst_kind,
                     trans_size * dst_kind, MPI_BYTE, win); chk_err(ierr);
    }
  }
  else if (dst_type == BT_CHARACTER && dst_kind == 1)
  {
    /* Get the required amount of memory on the stack. */
    void *dsh = alloca(dst_size);
    assign_char1_from_char4(dst_size, src_size, dsh, sr);
    ierr = MPI_Put(dsh, dst_size, MPI_BYTE, image_index, offset, dst_size,
                   MPI_BYTE, win); chk_err(ierr);
  }
  else if (dst_type == BT_CHARACTER)
  {
    /* Get the required amount of memory on the stack. */
    void *dsh = alloca(dst_size);
    assign_char4_from_char1(dst_size, src_size, dsh, sr);
    ierr = MPI_Put(dsh, dst_size, MPI_BYTE, image_index, offset, dst_size,
                   MPI_BYTE, win); chk_err(ierr);
  }
  else
  {
    /* Get the required amount of memory on the stack. */
    void *dsh = alloca(dst_size * num), *dsh_iter = dsh;
    dprint("type/kind convert %zd items: "
           "type %d(%d) -> type %d(%d), local buffer: %p\n",
           num, src_type, src_kind, dst_type, dst_kind, dsh);
    for (k = 0; k < num; ++k)
    {
      convert_type(dsh_iter, dst_type, dst_kind, sr, src_type, src_kind, stat);
      dsh_iter += dst_size;
      sr += src_size;
    }
    // dprint("dsh[0] = %d\n", ((int *)dsh)[0]);
    ierr = MPI_Put(dsh, dst_size * num, MPI_BYTE, image_index, offset,
                   dst_size * num, MPI_BYTE, win); chk_err(ierr);
  }
  ierr = MPI_Win_flush(image_index, win); chk_err(ierr);
}


static void
send_for_ref(caf_reference_t *ref, size_t *i, size_t src_index,
             mpi_caf_token_t *mpi_token, gfc_descriptor_t *dst,
             gfc_descriptor_t *src, void *ds, void *sr,
             ptrdiff_t dst_byte_offset, ptrdiff_t desc_byte_offset,
             int dst_kind, int src_kind, size_t dst_dim, size_t src_dim,
             size_t num, int *stat, int global_dynamic_win_rank, int memptr_win_rank,
             bool ds_global, /* access ds through global_dynamic_win */
             bool desc_global /* access desc through global_dynamic_win */
#ifdef GCC_GE_8
             , int dst_type)
{
#else
             )
{
  int dst_type = -1;
#endif
  ptrdiff_t extent_dst = 1, array_offset_dst = 0, dst_stride, src_stride;
  size_t next_dst_dim, ref_rank;
  gfc_max_dim_descriptor_t dst_desc_data;
  caf_ref_type_t ref_type = ref->type;
  caf_array_ref_t array_ref_src = ref->u.a.mode[src_dim];
  int ierr;

  if (unlikely(ref == NULL))
  {
    /* May be we should issue an error here, because this case should not
     * occur. */
    return;
  }

  dprint("Entering send_for_ref: [i = %zd] src_index = %zd, "
         "dst_offset = %zd, desc_offset = %zd, ds_glb = %d, desc_glb = %d\n",
         *i, src_index, dst_byte_offset, desc_byte_offset,
         ds_global, desc_global);

  if (ref->next == NULL)
  {
    size_t src_size = GFC_DESCRIPTOR_SIZE(src);
    dprint("[next == NULL]: src_size = %zd, ref_type = %s\n",
           src_size, caf_ref_type_str[ref_type]);

    switch (ref_type)
    {
      case CAF_REF_COMPONENT:
        dst_byte_offset += ref->u.c.offset;
        if (ref->u.c.caf_token_offset > 0)
        {
          if (ds_global)
          {
            ierr = MPI_Get(&ds, stdptr_size, MPI_BYTE, global_dynamic_win_rank,
                           MPI_Aint_add((MPI_Aint)ds, dst_byte_offset),
                           stdptr_size, MPI_BYTE, global_dynamic_win);
            chk_err(ierr);
            desc_global = true;
          }
          else
          {
            ierr = MPI_Get(&ds, stdptr_size, MPI_BYTE, memptr_win_rank,
                           dst_byte_offset, stdptr_size, MPI_BYTE,
                           mpi_token->memptr_win); chk_err(ierr);
            ds_global = true;
          }
          dst_byte_offset = 0;
        }

        if (ds_global)
        {
          put_data(NULL, MPI_Aint_add((MPI_Aint)ds, dst_byte_offset), sr,
#ifdef GCC_GE_8
                   dst_type,
#else
                   GFC_DESCRIPTOR_TYPE(src),
#endif
                   GFC_DESCRIPTOR_TYPE(src), dst_kind, src_kind,
                   ref->item_size, src_size, 1, stat, global_dynamic_win_rank);
        }
        else
        {
          put_data(mpi_token, dst_byte_offset, sr,
#ifdef GCC_GE_8
                   dst_type,
#else
                   GFC_DESCRIPTOR_TYPE(dst),
#endif
                   GFC_DESCRIPTOR_TYPE(src), dst_kind, src_kind,
                   ref->item_size, src_size, 1, stat, memptr_win_rank);
        }
        ++(*i);
        return;
      case CAF_REF_STATIC_ARRAY:
        dst_type = ref->u.a.static_array_type;
        /* Intentionally fall through. */
      case CAF_REF_ARRAY:
        if (array_ref_src == CAF_ARR_REF_NONE)
        {
          if (ds_global)
          {
            put_data(NULL, MPI_Aint_add((MPI_Aint)ds, dst_byte_offset),
                     sr + src_index * src_size,
#ifdef GCC_GE_8
                     dst_type, GFC_DESCRIPTOR_TYPE(src),
#else
                     GFC_DESCRIPTOR_TYPE(dst),
                     (dst_type == -1) ? GFC_DESCRIPTOR_TYPE(src) : dst_type,
#endif
                     dst_kind, src_kind, ref->item_size, src_size, num,
                     stat, global_dynamic_win_rank);
          }
          else
          {
            put_data(mpi_token, dst_byte_offset, sr + src_index * src_size,
#ifdef GCC_GE_8
                     dst_type, GFC_DESCRIPTOR_TYPE(src),
#else
                     GFC_DESCRIPTOR_TYPE(dst),
                     (dst_type == -1) ? GFC_DESCRIPTOR_TYPE(src) : dst_type,
#endif
                     dst_kind, src_kind, ref->item_size, src_size, num,
                     stat, memptr_win_rank);
          }
          *i += num;
          return;
        }
        break;
      default:
        caf_runtime_error(unreachable);
    }
  }
  caf_array_ref_t array_ref_dst = ref->u.a.mode[dst_dim];

#if 0
  dprint("image_index = %d, num = %zd, src_dim = %zd, dst_dim = %zd, "
         "ref_type = %s, array_ref_src = %s\n",
         image_index, num, src_dim, dst_dim,
         caf_ref_type_str[ref_type],
         caf_array_ref_str[array_ref_src]);
#endif

  switch (ref_type)
  {
    case CAF_REF_COMPONENT:
      if (ref->u.c.caf_token_offset > 0)
      {
        dst_byte_offset += ref->u.c.offset;
        desc_byte_offset = dst_byte_offset;
        if (ds_global)
        {
          ierr = MPI_Get(&ds, stdptr_size, MPI_BYTE, global_dynamic_win_rank,
                         MPI_Aint_add((MPI_Aint)ds, dst_byte_offset),
                         stdptr_size, MPI_BYTE, global_dynamic_win);
          chk_err(ierr);
          desc_global = true;
        }
        else
        {
          ierr = MPI_Get(&ds, stdptr_size, MPI_BYTE, memptr_win_rank,
                         dst_byte_offset, stdptr_size, MPI_BYTE,
                         mpi_token->memptr_win); chk_err(ierr);
          ds_global = true;
        }
        dst_byte_offset = 0;
      }
      else
      {
        dst_byte_offset += ref->u.c.offset;
        desc_byte_offset += ref->u.c.offset;
      }
      send_for_ref(ref->next, i, src_index, mpi_token, dst, src, ds,
                   sr, dst_byte_offset, desc_byte_offset, dst_kind, src_kind,
                   dst_dim, 0, 1, stat, global_dynamic_win_rank, memptr_win_rank,
                   ds_global, desc_global
#ifdef GCC_GE_8
                   , dst_type
#endif
                   );
      return;
    case CAF_REF_ARRAY:
      if (array_ref_src == CAF_ARR_REF_NONE)
      {
        send_for_ref(ref->next, i, src_index, mpi_token, dst, src, ds, sr,
                     dst_byte_offset, desc_byte_offset, dst_kind, src_kind,
                     dst_dim, 0, 1, stat, global_dynamic_win_rank, memptr_win_rank,
                     ds_global, desc_global
#ifdef GCC_GE_8
                     , dst_type
#endif
                    );
        return;
      }
      /* Only when on the left most index switch the data pointer to
       * the array's data pointer. */
      if (src_dim == 0)
      {
        if (ds_global)
        {
          for (ref_rank = 0; ref->u.a.mode[ref_rank] != CAF_ARR_REF_NONE;
               ++ref_rank) ;
          /* Get the remote descriptor. */
          if (desc_global)
          {
            ierr = MPI_Get(&dst_desc_data, sizeof_desc_for_rank(ref_rank),
                           MPI_BYTE, global_dynamic_win_rank,
                           MPI_Aint_add((MPI_Aint)ds, desc_byte_offset),
                           sizeof_desc_for_rank(ref_rank), MPI_BYTE,
                           global_dynamic_win); chk_err(ierr);
          }
          else
          {
            ierr = MPI_Get(&dst_desc_data, sizeof_desc_for_rank(ref_rank),
                           MPI_BYTE, memptr_win_rank, desc_byte_offset,
                           sizeof_desc_for_rank(ref_rank),
                           MPI_BYTE, mpi_token->memptr_win); chk_err(ierr);
            desc_global = true;
          }
          dst = (gfc_descriptor_t *)&dst_desc_data;
        }
        else
        {
          dst = mpi_token->desc;
        }
        dst_byte_offset = 0;
        desc_byte_offset = 0;
#ifdef EXTRA_DEBUG_OUTPUT
        dprint("remote desc rank: %zd (ref_rank: %zd)\n",
               GFC_DESCRIPTOR_RANK(src), ref_rank);
        for (int r = 0; r < GFC_DESCRIPTOR_RANK(src); ++r)
        {
          dprint("remote desc dim[%d] = (lb = %zd, ub = %zd, stride = %zd)\n",
                 r, src->dim[r].lower_bound, src->dim[r]._ubound,
                 src->dim[r]._stride);
        }
#endif
      }
      dprint("array_ref_dst[%zd] = %s := array_ref_src[%zd] = %s",
             dst_dim, caf_array_ref_str[array_ref_dst],
             src_dim, caf_array_ref_str[array_ref_src]);
      switch (array_ref_dst)
      {
        case CAF_ARR_REF_VECTOR:
          extent_dst = GFC_DESCRIPTOR_EXTENT(dst, dst_dim);
          array_offset_dst = 0;
          for (size_t idx = 0; idx < ref->u.a.dim[dst_dim].v.nvec; ++idx)
          {
#define KINDCASE(kind, type)                        \
case kind:                                          \
  array_offset_dst = (((ptrdiff_t)                  \
    ((type *)ref->u.a.dim[dst_dim].v.vector)[idx])  \
    - dst->dim[dst_dim].lower_bound                 \
    * dst->dim[dst_dim]._stride);                   \
  break

            switch (ref->u.a.dim[dst_dim].v.kind)
            {
              KINDCASE(1, int8_t);
              KINDCASE(2, int16_t);
              KINDCASE(4, int32_t);
              KINDCASE(8, int64_t);
#ifdef HAVE_GFC_INTEGER_16
              KINDCASE(16, __int128);
#endif
              default:
                caf_runtime_error(unreachable);
                return;
            }
#undef KINDCASE

            dprint("vector-index computed to: %zd\n", array_offset_dst);
            send_for_ref(ref, i, src_index, mpi_token, dst, src, ds, sr,
                         dst_byte_offset + array_offset_dst * ref->item_size,
                         desc_byte_offset + array_offset_dst * ref->item_size,
                         dst_kind, src_kind, dst_dim + 1, src_dim + 1,
                         1, stat, global_dynamic_win_rank, memptr_win_rank,
                         ds_global, desc_global
#ifdef GCC_GE_8
                         , dst_type
#endif
                         );
            src_index += dst->dim[dst_dim]._stride;
          }
          return;
        case CAF_ARR_REF_FULL:
          COMPUTE_NUM_ITEMS(extent_dst,
                            ref->u.a.dim[dst_dim].s.stride,
                            dst->dim[dst_dim].lower_bound,
                            dst->dim[dst_dim]._ubound);
          dst_stride = dst->dim[dst_dim]._stride
                       * ref->u.a.dim[dst_dim].s.stride;
          array_offset_dst = 0;
          src_stride = (GFC_DESCRIPTOR_RANK(src) > 0) ?
            src->dim[src_dim]._stride : 0;
          dprint("CAF_ARR_REF_FULL: src_stride = %zd, dst_stride = %zd\n",
                 src_stride, dst_stride);
          for (ptrdiff_t idx = 0; idx < extent_dst;
               ++idx, array_offset_dst += dst_stride)
          {
            send_for_ref(ref, i, src_index, mpi_token, dst, src, ds, sr,
                         dst_byte_offset + array_offset_dst * ref->item_size,
                         desc_byte_offset + array_offset_dst * ref->item_size,
                         dst_kind, src_kind, dst_dim + 1, src_dim + 1,
                         1, stat, global_dynamic_win_rank, memptr_win_rank,
                         ds_global, desc_global
#ifdef GCC_GE_8
                         , dst_type
#endif
                        );
            src_index += src_stride;
          }
          // dprint("CAF_ARR_REF_FULL: return, i = %zd\n", *i);
          return;

        case CAF_ARR_REF_RANGE:
          COMPUTE_NUM_ITEMS(extent_dst,
                            ref->u.a.dim[dst_dim].s.stride,
                            ref->u.a.dim[dst_dim].s.start,
                            ref->u.a.dim[dst_dim].s.end);
          array_offset_dst =
            (ref->u.a.dim[dst_dim].s.start - dst->dim[dst_dim].lower_bound)
            * dst->dim[dst_dim]._stride;
          dst_stride = dst->dim[dst_dim]._stride
                       * ref->u.a.dim[dst_dim].s.stride;
          src_stride = (GFC_DESCRIPTOR_RANK(src) > 0) ?
            src->dim[src_dim]._stride : 0;
          /* Increase the dst_dim only, when the src_extent is greater than one
           * or src and dst extent are both one. Don't increase when the
           * scalar source is not present in the dst. */
          next_dst_dim = (
            (extent_dst > 1) ||
            (GFC_DESCRIPTOR_EXTENT(src, src_dim) == 1 && extent_dst == 1)
          ) ? (dst_dim + 1) : dst_dim;
          for (ptrdiff_t idx = 0; idx < extent_dst; ++idx)
          {
            send_for_ref(ref, i, src_index, mpi_token, dst, src, ds, sr,
                         dst_byte_offset + array_offset_dst * ref->item_size,
                         desc_byte_offset + array_offset_dst * ref->item_size,
                         dst_kind, src_kind, next_dst_dim, src_dim + 1, 
                         1, stat, global_dynamic_win_rank, memptr_win_rank,
                         ds_global, desc_global
#ifdef GCC_GE_8
                         , dst_type
#endif
                         );
            src_index += src_stride;
            array_offset_dst += dst_stride;
          }
          // dprint("CAF_ARR_REF_RANGE: return, i = %zd\n", *i);
          return;

        case CAF_ARR_REF_SINGLE:
          array_offset_dst =
            (ref->u.a.dim[dst_dim].s.start - dst->dim[dst_dim].lower_bound)
            * dst->dim[dst_dim]._stride;
          // FIXME: issue #552
          // next_dst_dim = (
          //   (extent_dst > 1) ||
          //   (GFC_DESCRIPTOR_EXTENT(src, src_dim) == 1 && extent_dst == 1)
          // ) ? (dst_dim + 1) : dst_dim;
          next_dst_dim = dst_dim;
          send_for_ref(ref, i, src_index, mpi_token, dst, src, ds, sr,
                       dst_byte_offset + array_offset_dst * ref->item_size,
                       desc_byte_offset + array_offset_dst * ref->item_size,
                       dst_kind, src_kind, next_dst_dim, src_dim + 1,
                       1, stat, global_dynamic_win_rank, memptr_win_rank,
                       ds_global, desc_global
#ifdef GCC_GE_8
                       , dst_type
#endif
                       );

          // dprint("CAF_ARR_REF_SINGLE: return, i = %zd\n", *i);
          return;
        case CAF_ARR_REF_OPEN_END:
          COMPUTE_NUM_ITEMS(extent_dst,
                            ref->u.a.dim[dst_dim].s.stride,
                            ref->u.a.dim[dst_dim].s.start,
                            dst->dim[dst_dim]._ubound);
          dst_stride = dst->dim[dst_dim]._stride
                       * ref->u.a.dim[dst_dim].s.stride;
          src_stride = (GFC_DESCRIPTOR_RANK(src) > 0) ?
            src->dim[src_dim]._stride : 0;
          array_offset_dst =
            (ref->u.a.dim[dst_dim].s.start - dst->dim[dst_dim].lower_bound)
            * dst->dim[dst_dim]._stride;
          for (ptrdiff_t idx = 0; idx < extent_dst; ++idx)
          {
            send_for_ref(ref, i, src_index, mpi_token, dst, src, ds, sr,
                         dst_byte_offset + array_offset_dst * ref->item_size,
                         desc_byte_offset + array_offset_dst * ref->item_size,
                         dst_kind, src_kind, dst_dim + 1, src_dim + 1,
                         1, stat, global_dynamic_win_rank, memptr_win_rank,
                         ds_global, desc_global
#ifdef GCC_GE_8
                         , dst_type
#endif
                         );
            src_index += src_stride;
            array_offset_dst += dst_stride;
          }
          return;
        case CAF_ARR_REF_OPEN_START:
          COMPUTE_NUM_ITEMS(extent_dst,
                            ref->u.a.dim[dst_dim].s.stride,
                            dst->dim[dst_dim].lower_bound,
                            ref->u.a.dim[dst_dim].s.end);
          dst_stride =
            dst->dim[dst_dim]._stride * ref->u.a.dim[dst_dim].s.stride;
          src_stride = (GFC_DESCRIPTOR_RANK(src) > 0) ?
            src->dim[src_dim]._stride : 0;
          array_offset_dst = 0;
          for (ptrdiff_t idx = 0; idx < extent_dst; ++idx)
          {
            send_for_ref(ref, i, src_index, mpi_token, dst, src, ds, sr,
                         dst_byte_offset + array_offset_dst * ref->item_size,
                         desc_byte_offset + array_offset_dst * ref->item_size,
                         dst_kind, src_kind, dst_dim + 1, src_dim + 1,
                         1, stat, global_dynamic_win_rank, memptr_win_rank,
                         ds_global, desc_global
#ifdef GCC_GE_8
                         , dst_type
#endif
                         );
            src_index += src_stride;
            array_offset_dst += dst_stride;
          }
          return;
        default:
          caf_runtime_error(unreachable);
      }
      return;
    case CAF_REF_STATIC_ARRAY:
      if (array_ref_dst == CAF_ARR_REF_NONE)
      {
        send_for_ref(ref->next, i, src_index, mpi_token, dst, NULL, ds, sr,
                     dst_byte_offset, desc_byte_offset, dst_kind, src_kind,
                     dst_dim, 0, 1, stat, global_dynamic_win_rank, memptr_win_rank,
                     ds_global, desc_global
#ifdef GCC_GE_8
                     , dst_type
#endif
                     );
        return;
      }
      switch (array_ref_dst)
      {
        case CAF_ARR_REF_VECTOR:
          array_offset_dst = 0;
          for (size_t idx = 0; idx < ref->u.a.dim[dst_dim].v.nvec; ++idx)
          {
#define KINDCASE(kind, type)                                        \
case kind:                                                          \
  array_offset_dst = ((type *)ref->u.a.dim[dst_dim].v.vector)[idx]; \
  break

            switch (ref->u.a.dim[dst_dim].v.kind)
            {
              KINDCASE(1, int8_t);
              KINDCASE(2, int16_t);
              KINDCASE(4, int32_t);
              KINDCASE(8, int64_t);
#ifdef HAVE_GFC_INTEGER_16
              KINDCASE(16, __int128);
#endif
              default:
                caf_runtime_error(unreachable);
                return;
            }
#undef KINDCASE

            send_for_ref(ref, i, src_index, mpi_token, dst, NULL, ds, sr,
                         dst_byte_offset + array_offset_dst * ref->item_size,
                         desc_byte_offset + array_offset_dst * ref->item_size,
                         dst_kind, src_kind, dst_dim + 1, src_dim + 1,
                         1, stat, global_dynamic_win_rank, memptr_win_rank,
                         ds_global, desc_global
#ifdef GCC_GE_8
                         , dst_type
#endif
                         );
            src_index += src->dim[src_dim]._stride;
          }
          return;
        case CAF_ARR_REF_FULL:
          src_stride = (GFC_DESCRIPTOR_RANK(src) > 0) ?
            src->dim[src_dim]._stride : 0;
          for (array_offset_dst = 0 ;
               array_offset_dst <= ref->u.a.dim[dst_dim].s.end;
               array_offset_dst += ref->u.a.dim[dst_dim].s.stride)
          {
            send_for_ref(ref, i, src_index, mpi_token, dst, NULL, ds, sr,
                         dst_byte_offset + array_offset_dst * ref->item_size,
                         desc_byte_offset + array_offset_dst * ref->item_size,
                         dst_kind, src_kind, dst_dim + 1, src_dim + 1,
                         1, stat, global_dynamic_win_rank, memptr_win_rank,
                         ds_global, desc_global
#ifdef GCC_GE_8
                         , dst_type
#endif
                         );
            src_index += src_stride;
          }
          return;
        case CAF_ARR_REF_RANGE:
          COMPUTE_NUM_ITEMS(extent_dst,
                            ref->u.a.dim[dst_dim].s.stride,
                            ref->u.a.dim[dst_dim].s.start,
                            ref->u.a.dim[dst_dim].s.end);
          src_stride = (GFC_DESCRIPTOR_RANK (src) > 0) ?
            src->dim[src_dim]._stride : 0;
          array_offset_dst = ref->u.a.dim[dst_dim].s.start;
          for (ptrdiff_t idx = 0; idx < extent_dst; ++idx)
          {
            send_for_ref(ref, i, src_index, mpi_token, dst, NULL, ds, sr,
                         dst_byte_offset + array_offset_dst * ref->item_size,
                         desc_byte_offset + array_offset_dst * ref->item_size,
                         dst_kind, src_kind, dst_dim + 1, src_dim + 1,
                         1, stat, global_dynamic_win_rank, memptr_win_rank,
                         ds_global, desc_global
#ifdef GCC_GE_8
                         , dst_type
#endif
                         );
            src_index += src_stride;
            array_offset_dst += ref->u.a.dim[dst_dim].s.stride;
          }
          return;
        case CAF_ARR_REF_SINGLE:
          array_offset_dst = ref->u.a.dim[dst_dim].s.start;
          send_for_ref(ref, i, src_index, mpi_token, dst, NULL, ds, sr,
                       dst_byte_offset + array_offset_dst * ref->item_size,
                       desc_byte_offset + array_offset_dst * ref->item_size,
                       dst_kind, src_kind, dst_dim, src_dim + 1,
                       1, stat, global_dynamic_win_rank, memptr_win_rank,
                       ds_global, desc_global
#ifdef GCC_GE_8
                       , dst_type
#endif
                       );
          return;
          /* The OPEN_* are mapped to a RANGE and therefore can not occur. */
        case CAF_ARR_REF_OPEN_END:
        case CAF_ARR_REF_OPEN_START:
        default:
          caf_runtime_error(unreachable);
      }
      return;
    default:
      caf_runtime_error(unreachable);
  }
}


void
PREFIX(send_by_ref) (caf_token_t token, int image_index,
                     gfc_descriptor_t *src, caf_reference_t *refs,
                     int dst_kind, int src_kind, bool may_require_tmp,
                     bool dst_reallocatable, int *stat
#ifdef GCC_GE_8
                     , int dst_type
#endif
                     )
{
  const char vecrefunknownkind[] =
    "libcaf_mpi::caf_send_by_ref(): unknown kind in vector-ref.\n";
  const char unknownreftype[] =
    "libcaf_mpi::caf_send_by_ref(): unknown reference type.\n";
  const char unknownarrreftype[] =
    "libcaf_mpi::caf_send_by_ref(): unknown array reference type.\n";
  const char rankoutofrange[] =
    "libcaf_mpi::caf_send_by_ref(): rank out of range.\n";
  const char extentoutofrange[] =
    "libcaf_mpi::caf_send_by_ref(): extent out of range.\n";
  const char cannotallocdst[] =
    "libcaf_mpi::caf_send_by_ref(): can not allocate %d bytes of memory.\n";
  const char unabletoallocdst[] =
    "libcaf_mpi::caf_send_by_ref(): "
    "unable to allocate memory on remote image.\n";
  const char nonallocextentmismatch[] =
    "libcaf_mpi::caf_send_by_ref(): "
    "extent of non-allocatable arrays mismatch (%lu != %lu).\n";

  size_t size, i, ref_rank = 0, src_index, dst_size;
  int dst_rank = -1, src_cur_dim = 0, ierr;
  mpi_caf_token_t *mpi_token = (mpi_caf_token_t *) token;
  void *remote_memptr = mpi_token->memptr, *remote_base_memptr = NULL;
  gfc_max_dim_descriptor_t dst_desc, temp_src;
  gfc_descriptor_t *dst = (gfc_descriptor_t *)&dst_desc;
  caf_reference_t *riter = refs;
  long delta;
  ptrdiff_t data_offset = 0, desc_offset = 0;
  /* Reallocation of data on remote is needed (e.g., array to small).  This is
   * used for error tracking only.  It is not (yet) possible to allocate memory
   * on the remote image. */
  bool realloc_dst = false, extent_mismatch = false;
  /* Set when the first non-scalar array reference is encountered. */
  bool in_array_ref = false;
  /* Set when remote data is to be accessed through the
   * global dynamic window. */
  bool access_data_through_global_win = false;
  /* Set when the remote descriptor is to accessed through the global window. */
  bool access_desc_through_global_win = false;
  bool free_temp_src = false;
  caf_array_ref_t array_ref;

  if (stat)
    *stat = 0;

  MPI_Group current_team_group, win_group;
  int global_dynamic_win_rank, memptr_win_rank;
  ierr = MPI_Comm_group(CAF_COMM_WORLD, &current_team_group); chk_err(ierr);
  ierr = MPI_Win_get_group(global_dynamic_win, &win_group); chk_err(ierr);
  ierr = MPI_Group_translate_ranks(current_team_group, 1,
                                   (int[]){image_index - 1}, win_group,
                                   &global_dynamic_win_rank); chk_err(ierr);
  ierr = MPI_Group_free(&win_group); chk_err(ierr);
  ierr = MPI_Win_get_group(mpi_token->memptr_win, &win_group); chk_err(ierr);
  ierr = MPI_Group_translate_ranks(current_team_group, 1,
                                   (int[]){image_index - 1}, win_group,
                                   &memptr_win_rank); chk_err(ierr);
  ierr = MPI_Group_free(&current_team_group); chk_err(ierr);
  ierr = MPI_Group_free(&win_group); chk_err(ierr);

  check_image_health(global_dynamic_win_rank, stat);

#ifdef GCC_GE_8
  dprint("Entering send_by_ref(may_require_tmp = %d, dst_type = %d)\n",
         may_require_tmp, dst_type);
#else
  dprint("Entering send_by_ref(may_require_tmp = %d)\n", may_require_tmp);
#endif

  /* Compute the size of the result.  In the beginning size just counts the
   * number of elements. */
  size = 1;
  /* Shared lock both windows to prevent bother in the sub-routines. */
  CAF_Win_lock(MPI_LOCK_SHARED, global_dynamic_win_rank, global_dynamic_win);
  CAF_Win_lock(MPI_LOCK_SHARED, memptr_win_rank, mpi_token->memptr_win);
  while (riter)
  {
    dprint("remote_image = %d, offset = %zd, remote_mem = %p\n",
           global_dynamic_win_rank, data_offset, remote_memptr);
    switch (riter->type)
    {
      case CAF_REF_COMPONENT:
        if (riter->u.c.caf_token_offset > 0)
        {
          if (access_data_through_global_win)
          {
            data_offset += riter->u.c.offset;
            remote_base_memptr = remote_memptr;
            ierr = MPI_Get(&remote_memptr, stdptr_size, MPI_BYTE, global_dynamic_win_rank,
                           MPI_Aint_add((MPI_Aint)remote_memptr, data_offset),
                           stdptr_size, MPI_BYTE, global_dynamic_win);
            chk_err(ierr);
            /* On the second indirection access also the remote descriptor
             * using the global window. */
            access_desc_through_global_win = true;
          }
          else
          {
            data_offset += riter->u.c.offset;
            ierr = MPI_Get(&remote_memptr, stdptr_size, MPI_BYTE, memptr_win_rank,
                           data_offset, stdptr_size, MPI_BYTE,
                           mpi_token->memptr_win); chk_err(ierr);
            /* All future access is through the global dynamic window. */
            access_data_through_global_win = true;
          }
          desc_offset = data_offset;
          data_offset = 0;
        }
        else
        {
          data_offset += riter->u.c.offset;
          desc_offset += riter->u.c.offset;
        }
        break;
      case CAF_REF_ARRAY:
        /* When there has been no CAF_REF_COMP before hand, then the descriptor
         * is stored in the token and the extends are the same on all images,
         * which is taken care of in the else part. */
        if (access_data_through_global_win)
        {
          for (ref_rank = 0; 
               riter->u.a.mode[ref_rank] != CAF_ARR_REF_NONE; ++ref_rank) ;
          /* Get the remote descriptor and use the stack to store it
           * Note, dst may be pointing to mpi_token->desc therefore it
           * needs to be reset here. */
          dst = (gfc_descriptor_t *)&dst_desc;
          if (access_desc_through_global_win)
          {
            dprint("remote desc fetch from %p, offset = %zd\n",
                   remote_base_memptr, desc_offset);
            ierr = MPI_Get(dst, sizeof_desc_for_rank(ref_rank), MPI_BYTE,
                           global_dynamic_win_rank,
                           MPI_Aint_add(
                            (MPI_Aint)remote_base_memptr, desc_offset),
                           sizeof_desc_for_rank(ref_rank), MPI_BYTE,
                           global_dynamic_win); chk_err(ierr);
          }
          else
          {
            dprint("remote desc fetch from win %d, offset = %zd\n",
                   mpi_token->memptr_win, desc_offset);
            ierr = MPI_Get(dst, sizeof_desc_for_rank(ref_rank), MPI_BYTE,
                           memptr_win_rank, desc_offset,
                           sizeof_desc_for_rank(ref_rank), MPI_BYTE,
                           mpi_token->memptr_win); chk_err(ierr);
            access_desc_through_global_win = true;
          }
        }
        else
          dst = mpi_token->desc;
#ifdef EXTRA_DEBUG_OUTPUT
        dprint("remote desc rank: %zd (ref_rank: %zd)\n",
               GFC_DESCRIPTOR_RANK(dst), ref_rank);
        for (i = 0; i < GFC_DESCRIPTOR_RANK(dst); ++i)
        {
          dprint("remote desc dim[%zd] = (lb = %zd, ub = %zd, stride = %zd)\n",
                 i, dst->dim[i].lower_bound, dst->dim[i]._ubound,
                 dst->dim[i]._stride);
        }
#endif
        for (i = 0; riter->u.a.mode[i] != CAF_ARR_REF_NONE; ++i)
        {
          array_ref = riter->u.a.mode[i];
          dprint("i = %zd, array_ref = %s\n", i, caf_array_ref_str[array_ref]);
          switch (array_ref)
          {
            case CAF_ARR_REF_VECTOR:
              delta = riter->u.a.dim[i].v.nvec;
#define KINDCASE(kind, type)                                            \
case kind:                                                              \
  remote_memptr += (((ptrdiff_t)                                        \
    ((type *)riter->u.a.dim[i].v.vector)[0]) - src->dim[i].lower_bound) \
    * src->dim[i]._stride * riter->item_size;                           \
  break
              switch (riter->u.a.dim[i].v.kind)
              {
                KINDCASE(1, int8_t);
                KINDCASE(2, int16_t);
                KINDCASE(4, int32_t);
                KINDCASE(8, int64_t);
#ifdef HAVE_GFC_INTEGER_16
                KINDCASE(16, __int128);
#endif
                default:
                  caf_runtime_error(vecrefunknownkind, stat, NULL, 0);
                  return;
              }
#undef KINDCASE
              break;
            case CAF_ARR_REF_FULL:
              COMPUTE_NUM_ITEMS(delta,
                                riter->u.a.dim[i].s.stride,
                                dst->dim[i].lower_bound,
                                dst->dim[i]._ubound);
              /* The memptr stays unchanged when ref'ing the first element in
               * a dimension. */
              break;
            case CAF_ARR_REF_RANGE:
              COMPUTE_NUM_ITEMS(delta,
                                riter->u.a.dim[i].s.stride,
                                riter->u.a.dim[i].s.start,
                                riter->u.a.dim[i].s.end);
              remote_memptr +=
                (riter->u.a.dim[i].s.start - dst->dim[i].lower_bound)
                * dst->dim[i]._stride * riter->item_size;
              break;
            case CAF_ARR_REF_SINGLE:
              delta = 1;
              remote_memptr +=
                (riter->u.a.dim[i].s.start - dst->dim[i].lower_bound)
                * dst->dim[i]._stride * riter->item_size;
              break;
            case CAF_ARR_REF_OPEN_END:
              COMPUTE_NUM_ITEMS(delta,
                                riter->u.a.dim[i].s.stride,
                                riter->u.a.dim[i].s.start,
                                dst->dim[i]._ubound);
              remote_memptr +=
                (riter->u.a.dim[i].s.start - dst->dim[i].lower_bound)
                * dst->dim[i]._stride * riter->item_size;
              break;
            case CAF_ARR_REF_OPEN_START:
              COMPUTE_NUM_ITEMS(delta,
                                riter->u.a.dim[i].s.stride,
                                dst->dim[i].lower_bound,
                                riter->u.a.dim[i].s.end);
              /* The memptr stays unchanged when ref'ing the first element in
               * a dimension. */
              break;
            default:
              caf_runtime_error(unknownarrreftype, stat, NULL, 0);
              return;
          } // switch
          if (delta <= 0)
            return;
          if (dst != NULL)
            dst_rank = GFC_DESCRIPTOR_RANK(dst);
          /* Check the various properties of the destination array.
           * Is an array expected and present? */
          if (delta > 1 && dst_rank == 0)
          {
            /* No, an array is required, but not provided. */
            caf_runtime_error(extentoutofrange, stat, NULL, 0);
            return;
          }
          /* When dst is an array. */
          if (dst_rank > 0)
          {
            /* Check that src_cur_dim is valid for dst.  Can be superceeded
             * only by scalar data. */
            if (src_cur_dim >= dst_rank && delta != 1)
            {
              caf_runtime_error(rankoutofrange, stat, NULL, 0);
              return;
            }
            /* Do further checks, when the source is not scalar. */
            else if (delta != 1)
            {
              in_array_ref = true;
              /* When the realloc is required, then no extent may have
               * been set. */
              extent_mismatch = GFC_DESCRIPTOR_EXTENT(dst, src_cur_dim) < delta;
              /* When it already known, that a realloc is needed or
               * the extent does not match the needed one. */
              if (realloc_dst || extent_mismatch)
              {
                /* Check whether dst is reallocatable. */
                if (unlikely(!dst_reallocatable))
                {
                  caf_runtime_error(nonallocextentmismatch, stat,
                                    NULL, 0, delta,
                                    GFC_DESCRIPTOR_EXTENT(dst, src_cur_dim));
                  return;
                }
                /* Only report an error, when the extent needs to be
                 * modified, which is not allowed. */
                else if (!dst_reallocatable && extent_mismatch)
                {
                  caf_runtime_error(extentoutofrange, stat, NULL, 0);
                  return;
                }
                dprint("extent(dst, %d): %zd != delta: %ld.\n", src_cur_dim,
                       GFC_DESCRIPTOR_EXTENT(dst, src_cur_dim), delta);
                realloc_dst = true;
              }
            }

            if (src_cur_dim < GFC_DESCRIPTOR_RANK(src))
              ++src_cur_dim;
          }
          size *= (ptrdiff_t)delta;
        }
        in_array_ref = false;
        break;
      case CAF_REF_STATIC_ARRAY:
        for (i = 0; riter->u.a.mode[i] != CAF_ARR_REF_NONE; ++i)
        {
          array_ref = riter->u.a.mode[i];
          dprint("i = %zd, array_ref = %s\n", i, caf_array_ref_str[array_ref]);
          switch (array_ref)
          {
            case CAF_ARR_REF_VECTOR:
              delta = riter->u.a.dim[i].v.nvec;
#define KINDCASE(kind, type)                                    \
case kind:                                                      \
  remote_memptr +=                                              \
    ((type *)riter->u.a.dim[i].v.vector)[0] * riter->item_size; \
  break

              switch (riter->u.a.dim[i].v.kind)
              {
                KINDCASE(1, int8_t);
                KINDCASE(2, int16_t);
                KINDCASE(4, int32_t);
                KINDCASE(8, int64_t);
#ifdef HAVE_GFC_INTEGER_16
                KINDCASE(16, __int128);
#endif
                default:
                  caf_runtime_error(vecrefunknownkind, stat, NULL, 0);
                  return;
              }
#undef KINDCASE
              break;
            case CAF_ARR_REF_FULL:
              delta =
                riter->u.a.dim[i].s.end / riter->u.a.dim[i].s.stride + 1;
              /* The memptr stays unchanged when ref'ing the first element
               * in a dimension. */
              break;
            case CAF_ARR_REF_RANGE:
              COMPUTE_NUM_ITEMS(delta,
                                riter->u.a.dim[i].s.stride,
                                riter->u.a.dim[i].s.start,
                                riter->u.a.dim[i].s.end);
              remote_memptr += riter->u.a.dim[i].s.start
                * riter->u.a.dim[i].s.stride * riter->item_size;
              break;
            case CAF_ARR_REF_SINGLE:
              delta = 1;
              remote_memptr += riter->u.a.dim[i].s.start
                               * riter->u.a.dim[i].s.stride
                               * riter->item_size;
              break;
            case CAF_ARR_REF_OPEN_END:
              /* This and OPEN_START are mapped to a RANGE and therefore
               * can not occur here. */
            case CAF_ARR_REF_OPEN_START:
            default:
              caf_runtime_error(unknownarrreftype, stat, NULL, 0);
              return;
          } // switch
          if (delta <= 0)
            return;
          if (dst != NULL)
            dst_rank = GFC_DESCRIPTOR_RANK(dst);
          /* Check the various properties of the destination array.
           * Is an array expected and present? */
          if (delta > 1 && dst_rank == 0)
          {
            /* No, an array is required, but not provided. */
            caf_runtime_error(extentoutofrange, stat, NULL, 0);
            return;
          }
          /* When dst is an array. */
          if (dst_rank > 0)
          {
            /* Check that src_cur_dim is valid for dst.  Can be
             * superceeded only by scalar data. */
            if (src_cur_dim >= dst_rank && delta != 1)
            {
              caf_runtime_error(rankoutofrange, stat, NULL, 0);
              return;
            }
            /* Do further checks, when the source is not scalar. */
            else if (delta != 1)
            {
              in_array_ref = true;
              /* When the realloc is required, then no extent may have
               * been set. */
              extent_mismatch = GFC_DESCRIPTOR_EXTENT(dst, src_cur_dim) < delta;
              /* When it is already known, that a realloc is needed or
               * the extent does not match the needed one. */
              if (realloc_dst || extent_mismatch)
              {
                caf_runtime_error(unabletoallocdst, stat);
                return;
              }
            }
            if (src_cur_dim < GFC_DESCRIPTOR_RANK(src))
              ++src_cur_dim;
          }
          size *= (ptrdiff_t)delta;
        }
        in_array_ref = false;
        break;
      default:
        caf_runtime_error(unknownreftype, stat, NULL, 0);
        return;
    }
    dst_size = riter->item_size;
    riter = riter->next;
  }
  if (size == 0 || dst_size == 0)
    return;
  /* Postcondition:
   * - size contains the number of elements to store in the destination array,
   * - dst_size gives the size in bytes of each item in the destination array.
   */

  if (realloc_dst)
  {
    caf_runtime_error(unabletoallocdst, stat);
    return;
  }

  /* Reset the token. */
  mpi_token = (mpi_caf_token_t *) token;
  remote_memptr = mpi_token->memptr;
  src_index = 0;
#ifdef EXTRA_DEBUG_OUTPUT
  dprint("src_rank: %zd\n", GFC_DESCRIPTOR_RANK(src));
  for (i = 0; i < GFC_DESCRIPTOR_RANK(src); ++i)
  {
    dprint("src_dim[%zd] = (%zd, %zd)\n",
           i, src->dim[i].lower_bound, src->dim[i]._ubound);
  }
#endif
  /* When accessing myself and may_require_tmp is set, then copy the source
   * array. */
  if (caf_this_image == image_index && may_require_tmp)
  {
    dprint("preparing temporary source.\n");
    memcpy(&temp_src, src, sizeof_desc_for_rank(GFC_DESCRIPTOR_RANK(src)));
    size_t cap = 0;
    for (int r = 0; r < GFC_DESCRIPTOR_RANK(src); ++r)
    {
      cap += GFC_DESCRIPTOR_EXTENT(src, r);
    }

    cap *= GFC_DESCRIPTOR_SIZE(src);
    temp_src.base.base_addr = alloca(cap);
    if ((free_temp_src = (temp_src.base.base_addr == NULL)))
    {
      temp_src.base.base_addr = malloc(cap);
      if (temp_src.base.base_addr == NULL)
      {
        caf_runtime_error(cannotallocdst, stat, NULL, cap);
        return;
      }
    }
    memcpy(temp_src.base.base_addr, src->base_addr, cap);
    src = (gfc_descriptor_t *)&temp_src;
  }

  i = 0;
  dprint("calling send_for_ref. num elems: size = %zd, elem size in bytes: "
         "dst_size = %zd\n", size, dst_size);
  send_for_ref(refs, &i, src_index, mpi_token, mpi_token->desc, src,
               remote_memptr, src->base_addr, 0, 0, dst_kind, src_kind, 0, 0,
               1, stat, global_dynamic_win_rank, memptr_win_rank,
               false, false
#ifdef GCC_GE_8
               , dst_type
#endif
               );
  if (free_temp_src)
  {
    free(temp_src.base.base_addr);
  }
  CAF_Win_unlock(global_dynamic_win_rank, global_dynamic_win);
  CAF_Win_unlock(memptr_win_rank, mpi_token->memptr_win);
}


void
PREFIX(sendget_by_ref) (caf_token_t dst_token, int dst_image_index,
                        caf_reference_t *dst_refs, caf_token_t src_token,
                        int src_image_index, caf_reference_t *src_refs,
                        int dst_kind, int src_kind,
                        bool may_require_tmp, int *dst_stat, int *src_stat
#ifdef GCC_GE_8
                        , int dst_type, int src_type
#endif
                        )
{
  const char vecrefunknownkind[] =
    "libcaf_mpi::caf_sendget_by_ref(): unknown kind in vector-ref.\n";
  const char unknownreftype[] =
    "libcaf_mpi::caf_sendget_by_ref(): unknown reference type.\n";
  const char unknownarrreftype[] =
    "libcaf_mpi::caf_sendget_by_ref(): unknown array reference type.\n";
  const char cannotallocdst[] =
    "libcaf_mpi::caf_sendget_by_ref(): can not allocate %d bytes of memory.\n";
  size_t size, i, ref_rank, dst_index, src_index = 0, src_size;
  int dst_rank, ierr;
  mpi_caf_token_t
    *src_mpi_token = (mpi_caf_token_t *) src_token,
    *dst_mpi_token = (mpi_caf_token_t *) dst_token;
  void *remote_memptr = src_mpi_token->memptr, *remote_base_memptr = NULL;
  gfc_max_dim_descriptor_t src_desc;
  gfc_max_dim_descriptor_t temp_src_desc;
  gfc_descriptor_t *src = (gfc_descriptor_t *)&src_desc;
  caf_reference_t *riter = src_refs;
  long delta;
  ptrdiff_t data_offset = 0, desc_offset = 0;
  MPI_Group current_team_group, win_group;
  int global_dst_rank, global_src_rank, memptr_dst_rank, memptr_src_rank;
  /* Set when the first non-scalar array reference is encountered. */
  bool in_array_ref = false;
  bool array_extent_fixed = false;
  /* Set when remote data is to be accessed through the 
   * global dynamic window. */
  bool access_data_through_global_win = false;
  /* Set when the remote descriptor is to accessed through the global window. */
  bool access_desc_through_global_win = false;
  caf_array_ref_t array_ref;
#ifndef GCC_GE_8
  int dst_type = -1, src_type = -1;
#endif

  if (src_stat)
    *src_stat = 0;

  ierr = MPI_Comm_group(CAF_COMM_WORLD, &current_team_group); chk_err(ierr);

  ierr = MPI_Win_get_group(global_dynamic_win, &win_group); chk_err(ierr);
  ierr = MPI_Group_translate_ranks(current_team_group, 1,
                                  (int[]){src_image_index - 1}, win_group,
                                  &global_src_rank); chk_err(ierr);
  ierr = MPI_Group_translate_ranks(current_team_group, 1,
                                  (int[]){dst_image_index - 1}, win_group,
                                  &global_dst_rank); chk_err(ierr);
  ierr = MPI_Group_free(&win_group); chk_err(ierr);

  ierr = MPI_Win_get_group(src_mpi_token->memptr_win, &win_group); chk_err(ierr);
  ierr = MPI_Group_translate_ranks(current_team_group, 1,
                                  (int[]){src_image_index - 1}, win_group,
                                  &memptr_src_rank); chk_err(ierr);
  ierr = MPI_Group_free(&win_group); chk_err(ierr);
  ierr = MPI_Win_get_group(dst_mpi_token->memptr_win, &win_group); chk_err(ierr);
  ierr = MPI_Group_translate_ranks(current_team_group, 1,
                                  (int[]){dst_image_index - 1}, win_group,
                                  &memptr_dst_rank); chk_err(ierr);
  ierr = MPI_Group_free(&win_group); chk_err(ierr);
  ierr = MPI_Group_free(&current_team_group); chk_err(ierr);

  check_image_health(global_src_rank, src_stat);

  dprint("Entering get_by_ref(may_require_tmp = %d, dst_type = %d(%d), "
         "src_type = %d(%d)).\n",
         may_require_tmp, dst_type, dst_kind, src_type, src_kind);

  /* Compute the size of the result.  In the beginning size just counts the
   * number of elements. */
  size = 1;
  /* Shared lock both windows to prevent bother in the sub-routines. */
  CAF_Win_lock(MPI_LOCK_SHARED, global_src_rank, global_dynamic_win);
  CAF_Win_lock(MPI_LOCK_SHARED, memptr_src_rank, src_mpi_token->memptr_win);
  while (riter)
  {
    dprint("offset = %zd, remote_mem = %p\n", data_offset, remote_memptr);
    switch (riter->type)
    {
      case CAF_REF_COMPONENT:
        if (riter->u.c.caf_token_offset > 0)
        {
          if (access_data_through_global_win)
          {
            data_offset += riter->u.c.offset;
            remote_base_memptr = remote_memptr;
            ierr = MPI_Get(&remote_memptr, stdptr_size, MPI_BYTE,
                           global_src_rank,
                           MPI_Aint_add((MPI_Aint)remote_memptr, data_offset),
                           stdptr_size, MPI_BYTE, global_dynamic_win);
            chk_err(ierr);
            /* On the second indirection access also the remote descriptor
             * using the global window. */
            access_desc_through_global_win = true;
          }
          else
          {
            data_offset += riter->u.c.offset;
            ierr = MPI_Get(&remote_memptr, stdptr_size, MPI_BYTE,
                           memptr_src_rank, data_offset, stdptr_size, MPI_BYTE,
                           src_mpi_token->memptr_win); chk_err(ierr);
            /* All future access is through the global dynamic window. */
            access_data_through_global_win = true;
          }
          desc_offset = data_offset;
          data_offset = 0;
        }
        else
        {
          data_offset += riter->u.c.offset;
          desc_offset += riter->u.c.offset;
        }
        break;
      case CAF_REF_ARRAY:
        /* When there has been no CAF_REF_COMP before hand, then the
         * descriptor is stored in the token and the extends are the same on all
         * images, which is taken care of in the else part. */
        if (access_data_through_global_win)
        {
          for (ref_rank = 0; riter->u.a.mode[ref_rank] != CAF_ARR_REF_NONE;
               ++ref_rank) ;
          /* Get the remote descriptor and use the stack to store it. Note,
           * src may be pointing to mpi_token->desc therefore it needs to be
           * reset here. */
          src = (gfc_descriptor_t *)&src_desc;
          if (access_desc_through_global_win)
          {
            dprint("remote desc fetch from %p, offset = %zd\n",
                   remote_base_memptr, desc_offset);
            ierr = MPI_Get(src, sizeof_desc_for_rank(ref_rank), MPI_BYTE,
                           global_src_rank,
                           MPI_Aint_add(
                            (MPI_Aint)remote_base_memptr, desc_offset),
                           sizeof_desc_for_rank(ref_rank), MPI_BYTE,
                           global_dynamic_win); chk_err(ierr);
          }
          else
          {
            dprint("remote desc fetch from win %d, offset = %zd\n",
                   src_mpi_token->memptr_win, desc_offset);
            ierr = MPI_Get(src, sizeof_desc_for_rank(ref_rank), MPI_BYTE,
                           memptr_src_rank, desc_offset,
                           sizeof_desc_for_rank(ref_rank),
                           MPI_BYTE, src_mpi_token->memptr_win); chk_err(ierr);
            access_desc_through_global_win = true;
          }
        }
        else
        {
          src = src_mpi_token->desc;
        }
#ifdef EXTRA_DEBUG_OUTPUT
        dprint("remote desc rank: %zd (ref_rank: %zd)\n",
               GFC_DESCRIPTOR_RANK(src), ref_rank);
        for (i = 0; i < GFC_DESCRIPTOR_RANK(src); ++i)
        {
          dprint("remote desc dim[%zd] = (lb = %zd, ub = %zd, stride = %zd)\n",
                 i, src->dim[i].lower_bound, src->dim[i]._ubound,
                 src->dim[i]._stride);
        }
#endif
        for (i = 0; riter->u.a.mode[i] != CAF_ARR_REF_NONE; ++i)
        {
          array_ref = riter->u.a.mode[i];
          dprint("i = %zd, array_ref = %s\n", i, caf_array_ref_str[array_ref]);
          switch (array_ref)
          {
            case CAF_ARR_REF_VECTOR:
              delta = riter->u.a.dim[i].v.nvec;
#define KINDCASE(kind, type)                                            \
case kind:                                                              \
  remote_memptr += (((ptrdiff_t)                                        \
    ((type *)riter->u.a.dim[i].v.vector)[0]) - src->dim[i].lower_bound) \
    * src->dim[i]._stride * riter->item_size;                           \
  break
              switch (riter->u.a.dim[i].v.kind)
              {
                KINDCASE(1, int8_t);
                KINDCASE(2, int16_t);
                KINDCASE(4, int32_t);
                KINDCASE(8, int64_t);
#ifdef HAVE_GFC_INTEGER_16
                KINDCASE(16, __int128);
#endif
                default:
                  caf_runtime_error(vecrefunknownkind, src_stat, NULL, 0);
                  return;
              }
#undef KINDCASE
              break;
            case CAF_ARR_REF_FULL:
              COMPUTE_NUM_ITEMS(delta,
                                riter->u.a.dim[i].s.stride,
                                src->dim[i].lower_bound,
                                src->dim[i]._ubound);
              /* The memptr stays unchanged when ref'ing the first element
               * in a dimension. */
              break;
            case CAF_ARR_REF_RANGE:
              COMPUTE_NUM_ITEMS(delta,
                                riter->u.a.dim[i].s.stride,
                                riter->u.a.dim[i].s.start,
                                riter->u.a.dim[i].s.end);
              remote_memptr += 
                (riter->u.a.dim[i].s.start - src->dim[i].lower_bound)
                * src->dim[i]._stride * riter->item_size;
              break;
            case CAF_ARR_REF_SINGLE:
              delta = 1;
              remote_memptr +=
                (riter->u.a.dim[i].s.start - src->dim[i].lower_bound)
                * src->dim[i]._stride * riter->item_size;
              break;
            case CAF_ARR_REF_OPEN_END:
              COMPUTE_NUM_ITEMS(delta,
                                riter->u.a.dim[i].s.stride,
                                riter->u.a.dim[i].s.start,
                                src->dim[i]._ubound);
              remote_memptr +=
                (riter->u.a.dim[i].s.start - src->dim[i].lower_bound)
                * src->dim[i]._stride * riter->item_size;
              break;
            case CAF_ARR_REF_OPEN_START:
              COMPUTE_NUM_ITEMS(delta,
                                riter->u.a.dim[i].s.stride,
                                src->dim[i].lower_bound,
                                riter->u.a.dim[i].s.end);
              /* The memptr stays unchanged when ref'ing the first element
               * in a dimension. */
              break;
            default:
              caf_runtime_error(unknownarrreftype, src_stat, NULL, 0);
              return;
          } // switch
          if (delta <= 0)
            return;
          size *= (ptrdiff_t)delta;
        }
        if (in_array_ref)
        {
          array_extent_fixed = true;
          in_array_ref = false;
        }
        break;
      case CAF_REF_STATIC_ARRAY:
        for (i = 0; riter->u.a.mode[i] != CAF_ARR_REF_NONE; ++i)
        {
          array_ref = riter->u.a.mode[i];
          dprint("i = %zd, array_ref = %s\n", i, caf_array_ref_str[array_ref]);
          switch (array_ref)
          {
            case CAF_ARR_REF_VECTOR:
              delta = riter->u.a.dim[i].v.nvec;
#define KINDCASE(kind, type)                                      \
case kind:                                                        \
  remote_memptr +=                                                \
    ((type *)riter->u.a.dim[i].v.vector)[0] * riter->item_size;   \
  break
              switch (riter->u.a.dim[i].v.kind)
              {
                KINDCASE(1, int8_t);
                KINDCASE(2, int16_t);
                KINDCASE(4, int32_t);
                KINDCASE(8, int64_t);
#ifdef HAVE_GFC_INTEGER_16
                KINDCASE(16, __int128);
#endif
                default:
                  caf_runtime_error(vecrefunknownkind, src_stat, NULL, 0);
                  return;
              }
#undef KINDCASE
              break;
            case CAF_ARR_REF_FULL:
              delta = riter->u.a.dim[i].s.end / riter->u.a.dim[i].s.stride + 1;
              /* The memptr stays unchanged when ref'ing the first element
               * in a dimension. */
              break;
            case CAF_ARR_REF_RANGE:
              COMPUTE_NUM_ITEMS(delta,
                                riter->u.a.dim[i].s.stride,
                                riter->u.a.dim[i].s.start,
                                riter->u.a.dim[i].s.end);
              remote_memptr += riter->u.a.dim[i].s.start
                               * riter->u.a.dim[i].s.stride
                               * riter->item_size;
              break;
            case CAF_ARR_REF_SINGLE:
              delta = 1;
              remote_memptr += riter->u.a.dim[i].s.start
                               * riter->u.a.dim[i].s.stride
                               * riter->item_size;
              break;
            case CAF_ARR_REF_OPEN_END:
              /* This and OPEN_START are mapped to a RANGE and therefore
               * can not occur here. */
            case CAF_ARR_REF_OPEN_START:
            default:
              caf_runtime_error(unknownarrreftype, src_stat, NULL, 0);
              return;
          } // switch
          if (delta <= 0)
            return;
          size *= (ptrdiff_t)delta;
        }
        if (in_array_ref)
        {
          array_extent_fixed = true;
          in_array_ref = false;
        }
        break;
      default:
        caf_runtime_error(unknownreftype, src_stat, NULL, 0);
        return;
      } // switch
      src_size = riter->item_size;
      riter = riter->next;
    }
  if (size == 0 || src_size == 0)
    return;
  /* Postcondition:
   * - size contains the number of elements to store in the destination array,
   * - src_size gives the size in bytes of each item in the destination array.
  */

  dst_rank = (size > 1) ? 1 : 0;
  memset(&temp_src_desc, 0, sizeof(gfc_dim1_descriptor_t));
#ifdef GCC_GE_8
  temp_src_desc.base.dtype.elem_len = (dst_type != BT_COMPLEX) ?
    dst_kind : (2 * dst_kind);
  temp_src_desc.base.dtype.rank = 1;
  temp_src_desc.base.dtype.type = dst_type;
#else // GCC_GE_7
  temp_src_desc.base.dtype = GFC_DTYPE_INTEGER_4 | 1;
#endif
  temp_src_desc.base.offset = 0;
  temp_src_desc.dim[0]._ubound = size - 1;
  temp_src_desc.dim[0]._stride = 1;

  temp_src_desc.base.base_addr =
    malloc(size * GFC_DESCRIPTOR_SIZE((gfc_descriptor_t *)&temp_src_desc));
  if (unlikely(temp_src_desc.base.base_addr == NULL))
  {
    caf_runtime_error(
      cannotallocdst, src_stat,
      size * GFC_DESCRIPTOR_SIZE((gfc_descriptor_t *)&temp_src_desc));
    return;
  }

#ifndef GCC_GE_8
  static bool warning_given = false;
  if (!warning_given)
  {
    fprintf(stderr,
            "lib_caf_mpi::sendget_by_ref(): Warning !! sendget_by_ref() is "
            "mostly unfunctional due to a design error. Split up your "
            "statement with coarray refs on both sides of the assignment "
            "when the datatype transfered is non 4-byte-integer compatible, "
            "or use gcc >= 8.\n");
    warning_given = true;
  }
#endif
  /* Reset the token. */
  src_mpi_token = (mpi_caf_token_t *) src_token;
  remote_memptr = src_mpi_token->memptr;
  dst_index = 0;
#ifdef EXTRA_DEBUG_OUTPUT
  dprint("dst_rank: %d\n", dst_rank);
  for (i = 0; i < dst_rank; ++i)
  {
    dprint("temp_src_dim[%zd] = (%zd, %zd)\n",
           i, temp_src_desc.dim[i].lower_bound, temp_src_desc.dim[i]._ubound);
  }
#endif
  i = 0;
  dprint("calling get_for_ref.\n");
  get_for_ref(src_refs, &i, dst_index, src_mpi_token,
              (gfc_descriptor_t *)&temp_src_desc, src_mpi_token->desc,
              temp_src_desc.base.base_addr, remote_memptr, 0, 0, dst_kind,
              src_kind, 0, 0, 1, src_stat, global_src_rank, memptr_src_rank,
              false, false
#ifdef GCC_GE_8
              , src_type
#endif
              );
  CAF_Win_unlock(global_src_rank, global_dynamic_win);
  CAF_Win_unlock(memptr_src_rank, src_mpi_token->memptr_win);
  dprint("calling send_for_ref. num elems: size = %zd, elem size in bytes: "
         "src_size = %zd\n", size, src_size);
  i = 0;

  CAF_Win_lock(MPI_LOCK_EXCLUSIVE, global_dst_rank, global_dynamic_win);
  CAF_Win_lock(MPI_LOCK_EXCLUSIVE, memptr_dst_rank, dst_mpi_token->memptr_win);
  send_for_ref(dst_refs, &i, src_index, dst_mpi_token, dst_mpi_token->desc,
               (gfc_descriptor_t *)&temp_src_desc, dst_mpi_token->memptr,
               temp_src_desc.base.base_addr, 0, 0, dst_kind, src_kind, 0, 0,
               1, dst_stat, global_dst_rank, memptr_dst_rank, false, false
#ifdef GCC_GE_8
               , dst_type
#endif
               );
  CAF_Win_unlock(global_dst_rank, global_dynamic_win);
  CAF_Win_unlock(memptr_dst_rank, src_mpi_token->memptr_win);
}

int
PREFIX(is_present) (caf_token_t token, int image_index, caf_reference_t *refs)
{
  const char unsupportedRefType[] =
    "Unsupported ref-type in caf_is_present().";
  const char unexpectedEndOfRefs[] =
    "Unexpected end of references in caf_is_present.";
  const char remotesInnerRefNA[] =
    "Memory referenced on the remote image is not allocated.";
  const int ptr_size = sizeof(void *);
  const int remote_image = image_index - 1;
  mpi_caf_token_t *mpi_token = (mpi_caf_token_t *)token;
  ptrdiff_t local_offset = 0;
  void *remote_memptr = NULL, *remote_base_memptr = NULL;
  bool carryOn = true, firstDesc = true;
  caf_reference_t *riter = refs, *prev;
  size_t i, ref_rank;
  int ierr;
  gfc_max_dim_descriptor_t src_desc;
  caf_array_ref_t array_ref;

  while (carryOn && riter)
  {
    switch (riter->type)
    {
      case CAF_REF_COMPONENT:
        if (riter->u.c.caf_token_offset)
        {
          CAF_Win_lock(MPI_LOCK_SHARED, remote_image, mpi_token->memptr_win);
          ierr = MPI_Get(&remote_memptr, ptr_size, MPI_BYTE, remote_image,
                         local_offset + riter->u.c.offset, ptr_size,
                         MPI_BYTE, mpi_token->memptr_win); chk_err(ierr);
          CAF_Win_unlock(remote_image, mpi_token->memptr_win);
          dprint("Got first remote address %p from offset %zd\n",
                 remote_memptr, local_offset);
          local_offset = 0;
          carryOn = false;
        }
        else
          local_offset += riter->u.c.offset;
        break;
      case CAF_REF_ARRAY:
        {
          const gfc_descriptor_t *src = 
            (gfc_descriptor_t *)(mpi_token->memptr + local_offset);
          for (i = 0; riter->u.a.mode[i] != CAF_ARR_REF_NONE; ++i)
          {
            array_ref = riter->u.a.mode[i];
            dprint("i = %zd, array_ref = %s\n", i, caf_array_ref_str[array_ref]);
            switch (array_ref)
            {
              case CAF_ARR_REF_FULL:
                /* The local_offset stays unchanged when ref'ing the first
                 * element in a dimension. */
                break;
              case CAF_ARR_REF_SINGLE:
                local_offset +=
                  (riter->u.a.dim[i].s.start - src->dim[i].lower_bound)
                  * src->dim[i]._stride * riter->item_size;
                break;
              case CAF_ARR_REF_VECTOR:
              case CAF_ARR_REF_RANGE:
              case CAF_ARR_REF_OPEN_END:
              case CAF_ARR_REF_OPEN_START:
                /* Intentionally fall through, because these are not
                 * suported here. */
              default:
                caf_runtime_error(unsupportedRefType);
                return false;
            }
          }
        }
        break;
      case CAF_REF_STATIC_ARRAY:
        for (i = 0; riter->u.a.mode[i] != CAF_ARR_REF_NONE; ++i)
        {
          array_ref = riter->u.a.mode[i];
          dprint("i = %zd, array_ref = %s\n", i, caf_array_ref_str[array_ref]);
          switch (array_ref)
          {
            case CAF_ARR_REF_FULL:
              /* The local_offset stays unchanged when ref'ing the first
               * element in a dimension. */
              break;
            case CAF_ARR_REF_SINGLE:
              local_offset += riter->u.a.dim[i].s.start
                              * riter->u.a.dim[i].s.stride * riter->item_size;
              break;
            case CAF_ARR_REF_VECTOR:
            case CAF_ARR_REF_RANGE:
            case CAF_ARR_REF_OPEN_END:
            case CAF_ARR_REF_OPEN_START:
            default:
              caf_runtime_error(unsupportedRefType);
              return false;
          }
        }
        break;
      default:
        caf_runtime_error(unsupportedRefType);
        return false;
    } // switch
    prev = riter;
    riter = riter->next;
  }

  if (carryOn)
  {
    // This can only happen, when riter == NULL.
    caf_runtime_error(unexpectedEndOfRefs);
  }

  CAF_Win_lock(MPI_LOCK_SHARED, remote_image, global_dynamic_win);
  if (remote_memptr != NULL)
    remote_base_memptr = remote_memptr + local_offset;

  dprint("Remote desc address is %p from remote memptr %p and offset %zd\n",
         remote_base_memptr, remote_memptr, local_offset);

  while (riter)
  {
    switch (riter->type)
    {
      case CAF_REF_COMPONENT:
        /* After reffing the first allocatable/pointer component, descriptors
         * need to be picked up from the global_win. */
        firstDesc = firstDesc && riter->u.c.caf_token_offset == 0;
        local_offset += riter->u.c.offset;
        remote_base_memptr = remote_memptr + local_offset;
        ierr = MPI_Get(&remote_memptr, ptr_size, MPI_BYTE, remote_image,
                       (MPI_Aint)remote_base_memptr, ptr_size,
                       MPI_BYTE, global_dynamic_win); chk_err(ierr);
        dprint("Got remote address %p from offset %zd nd base memptr %p\n",
               remote_memptr, local_offset, remote_base_memptr);
        local_offset = 0;
        break;
      case CAF_REF_ARRAY:
        if (remote_base_memptr == NULL)
        {
          /* Refing an unallocated array ends in a full_ref. Check that this
           * is true. Error when not full-refing. */
          for (i = 0; riter->u.a.mode[i] != CAF_ARR_REF_NONE; ++i)
          {
            if (riter->u.a.mode[i] != CAF_ARR_REF_FULL)
              break;
          }
          if (riter->u.a.mode[i] != CAF_ARR_REF_NONE)
            caf_runtime_error(remotesInnerRefNA);
          break;
        }
        if (firstDesc)
        {
          /* The first descriptor is accessible by the mpi_token->memptr_win.
           * Count the dims to fetch. */
          for (ref_rank = 0; riter->u.a.mode[ref_rank] != CAF_ARR_REF_NONE;
               ++ref_rank)
	    ;
          dprint("Getting remote descriptor of rank %zd from win: %d, "
                 "sizeof() %zd\n", ref_rank, mpi_token->memptr_win,
                 sizeof_desc_for_rank(ref_rank));
          ierr = MPI_Get(&src_desc, sizeof_desc_for_rank(ref_rank),
                         MPI_BYTE, remote_image, local_offset,
                         sizeof_desc_for_rank(ref_rank),
                         MPI_BYTE, mpi_token->memptr_win); chk_err(ierr);
          firstDesc = false;
        }
        else
        {
          /* All inner descriptors go by the dynamic window.
           * Count the dims to fetch. */
          for (ref_rank = 0; riter->u.a.mode[ref_rank] != CAF_ARR_REF_NONE;
               ++ref_rank)
	    ;
          dprint("Getting remote descriptor of rank %zd from: %p, "
                 "sizeof() %zd\n", ref_rank, remote_base_memptr,
                 sizeof_desc_for_rank(ref_rank));
          ierr = MPI_Get(&src_desc, sizeof_desc_for_rank(ref_rank), MPI_BYTE,
                         remote_image, (MPI_Aint)remote_base_memptr,
                         sizeof_desc_for_rank(ref_rank), MPI_BYTE,
                         global_dynamic_win); chk_err(ierr);
        }
#ifdef EXTRA_DEBUG_OUTPUT
        {
          gfc_descriptor_t * src = (gfc_descriptor_t *)(&src_desc);
          dprint("remote desc rank: %zd (ref_rank: %zd)\n",
                 GFC_DESCRIPTOR_RANK(src), ref_rank);
          for (i = 0; i < GFC_DESCRIPTOR_RANK(src); ++i)
          {
            dprint("remote desc dim[%zd] = "
                   "(lb = %zd, ub = %zd, stride = %zd)\n",
                   i, src_desc.dim[i].lower_bound, src_desc.dim[i]._ubound,
                   src_desc.dim[i]._stride);
          }
        }
#endif

        for (i = 0; riter->u.a.mode[i] != CAF_ARR_REF_NONE; ++i)
        {
          array_ref = riter->u.a.mode[i];
          dprint("i = %zd, array_ref = %s\n", i, caf_array_ref_str[array_ref]);
          switch (array_ref)
          {
            case CAF_ARR_REF_FULL:
              /* The local_offset stays unchanged when ref'ing the first 
               * element in a dimension. */
              break;
            case CAF_ARR_REF_SINGLE:
              local_offset +=
                (riter->u.a.dim[i].s.start - src_desc.dim[i].lower_bound)
                * src_desc.dim[i]._stride
                * riter->item_size;
              break;
            case CAF_ARR_REF_VECTOR:
            case CAF_ARR_REF_RANGE:
            case CAF_ARR_REF_OPEN_END:
            case CAF_ARR_REF_OPEN_START:
              /* Intentionally fall through, because these are not suported
               * here. */
            default:
              caf_runtime_error(unsupportedRefType);
              CAF_Win_unlock(remote_image, global_dynamic_win);
              return false;
          }
        }
        break;
      case CAF_REF_STATIC_ARRAY:
        for (i = 0; riter->u.a.mode[i] != CAF_ARR_REF_NONE; ++i)
        {
          array_ref = riter->u.a.mode[i];
          dprint("i = %zd, array_ref = %s\n", i, caf_array_ref_str[array_ref]);
          switch (array_ref)
          {
            case CAF_ARR_REF_FULL:
              /* The memptr stays unchanged when ref'ing the first element
               * in a dimension. */
              break;
            case CAF_ARR_REF_SINGLE:
              local_offset += riter->u.a.dim[i].s.start
                              * riter->u.a.dim[i].s.stride * riter->item_size;
              break;
            case CAF_ARR_REF_VECTOR:
            case CAF_ARR_REF_RANGE:
            case CAF_ARR_REF_OPEN_END:
            case CAF_ARR_REF_OPEN_START:
            default:
              caf_runtime_error(unsupportedRefType);
              CAF_Win_unlock(remote_image, global_dynamic_win);
              return false;
          }
        }
        break;
      default:
        caf_runtime_error(unsupportedRefType);
        CAF_Win_unlock(remote_image, global_dynamic_win);
        return false;
    } // switch
    riter = riter->next;
  }
  CAF_Win_unlock(remote_image, global_dynamic_win);

  dprint("Got remote_memptr: %p\n", remote_memptr);
  return remote_memptr != NULL;
}
#endif // GCC_GE_7


/* SYNC IMAGES. Note: SYNC IMAGES(*) is passed as count == -1 while
 * SYNC IMAGES([]) has count == 0. Note further that SYNC IMAGES(*)
 * is not semantically equivalent to SYNC ALL. */

void
PREFIX(sync_images) (int count, int images[], int *stat, char *errmsg,
                     charlen_t errmsg_len)
{
  sync_images_internal(count, images, stat, errmsg, errmsg_len, false);
}

static void
sync_images_internal(int count, int images[], int *stat, char *errmsg,
                     size_t errmsg_len, bool internal)
{
  int ierr = 0, i = 0, j = 0, int_zero = 0, done_count = 0, flag;
  MPI_Status s;

#ifdef WITH_FAILED_IMAGES
  no_stopped_images_check_in_errhandler = true;
#endif
  dprint("Entering\n");
  if (count == 0 || (count == 1 && images[0] == caf_this_image))
  {
    if (stat)
      *stat = 0;
#ifdef WITH_FAILED_IMAGES
    no_stopped_images_check_in_errhandler = false;
#endif
    dprint("Leaving early.\n");
    return;
  }

  /* halt execution if sync images contains duplicate image numbers */
  for (i = 0; i < count; ++i)
  {
    for (j = 0; j < i; ++j)
    {
      if (images[i] == images[j])
      {
        ierr = STAT_DUP_SYNC_IMAGES;
        if (stat)
          *stat = ierr;
        goto sync_images_err_chk;
      }
    }
  }

#ifdef GFC_CAF_CHECK
    for (i = 0; i < count; ++i)
    {
      if (images[i] < 1 || images[i] > caf_num_images)
      {
        fprintf(stderr, "COARRAY ERROR: Invalid image index %d to SYNC IMAGES",
                images[i]);
        terminate_internal(1, 1);
      }
    }
#endif

  if (unlikely(caf_is_finalized))
  {
    ierr = STAT_STOPPED_IMAGE;
  }
  else
  {
    if (count == -1)
    {
      count = caf_num_images - 1;
      images = images_full;
    }

#if defined(NONBLOCKING_PUT) && !defined(CAF_MPI_LOCK_UNLOCK)
    explicit_flush();
#endif

#ifdef WITH_FAILED_IMAGES
    /* Provoke detecting process fails. */
    ierr = MPI_Test(&alive_request, &flag, MPI_STATUS_IGNORE); chk_err(ierr);
#endif
    /* A rather simple way to synchronice:
     * - expect all images to sync with receiving an int,
     * - on the other side, send all processes to sync with an int,
     * - when the int received is STAT_STOPPED_IMAGE the return immediately,
     *   else wait until all images in the current set of images have send
     *   some data, i.e., synced.
     *
     * This approach as best as possible implements the syncing of different
     * sets of images and figuring that an image has stopped.  MPI does not
     * provide any direct means of syncing non-coherent sets of images.
     * The groups/communicators of MPI always need to be consistent, i.e.,
     * have the same members on all images participating.  This is
     * contradictiory to the sync images statement, where syncing, e.g., in a
     * ring pattern is possible.
     *
     * This implementation guarantees, that as long as no image is stopped
     * an image only is allowed to continue, when all its images to sync to
     * also have reached a sync images statement.  This implementation makes
     * no assumption when the image continues or in which order synced
     * images continue. */
    for (i = 0; i < count; ++i)
    {
      /* Need to have the request handlers contigously in the handlers
       * array or waitany below will trip about the handler as illegal. */
      ierr = MPI_Irecv(&arrived[images[i] - 1], 1, MPI_INT, images[i] - 1,
                       MPI_TAG_CAF_SYNC_IMAGES, CAF_COMM_WORLD,
                       &sync_handles[i]); chk_err(ierr);
    }
    for (i = 0; i < count; ++i)
    {
      ierr = MPI_Send(&int_zero, 1, MPI_INT, images[i] - 1,
                      MPI_TAG_CAF_SYNC_IMAGES, CAF_COMM_WORLD); chk_err(ierr);
    }
    done_count = 0;
    while (done_count < count)
    {
      ierr = MPI_Waitany(count, sync_handles, &i, &s);
      if (ierr == MPI_SUCCESS && i != MPI_UNDEFINED)
      {
        ++done_count;
        if (ierr == MPI_SUCCESS && arrived[s.MPI_SOURCE] == STAT_STOPPED_IMAGE)
        {
          /* Possible future extension: Abort pending receives.  At the
           * moment the receives are discarded by the program
           * termination.  For the tested mpi-implementation this is ok. */
          ierr = STAT_STOPPED_IMAGE;
          break;
        }
      }
      else if (ierr != MPI_SUCCESS)
#ifdef WITH_FAILED_IMAGES
      {
        int err;
        MPI_Error_class(ierr, &err);
        if (err == MPIX_ERR_PROC_FAILED)
        {
          dprint("Image failed, provoking error handling.\n");
          ierr = STAT_FAILED_IMAGE;
          /* Provoke detecting process fails. */
          MPI_Test(&alive_request, &flag, MPI_STATUS_IGNORE);
        }
        break;
      }
#else
        break;
#endif // WITH_FAILED_IMAGES
    }
  }

sync_images_err_chk:
#ifdef WITH_FAILED_IMAGES
  no_stopped_images_check_in_errhandler = false;
#endif
  dprint("Leaving\n");
  if (stat)
    *stat = ierr;
#ifdef WITH_FAILED_IMAGES
  else if (ierr == STAT_FAILED_IMAGE)
    terminate_internal(ierr, 0);
#endif

  if (ierr != 0 && ierr != STAT_FAILED_IMAGE)
  {
    char msg[80];
    strcpy(msg, "SYNC IMAGES failed");
    if (caf_is_finalized)
      strcat(msg, " - there are stopped images");

    if (errmsg_len > 0)
    {
      size_t len = (strlen(msg) > errmsg_len) ? errmsg_len : strlen (msg);
      memcpy(errmsg, msg, len);
      if (errmsg_len > len)
        memset(&errmsg[len], ' ', errmsg_len-len);
    }
    else if (!internal && stat == NULL)
      caf_runtime_error(msg);
  }
}


#define GEN_REDUCTION(name, datatype, operator)       \
static void                                           \
name(datatype *invec, datatype *inoutvec, int *len,   \
     MPI_Datatype *datatype __attribute__((unused)))  \
{                                                     \
  for (int i = 0; i < len; ++i)                       \
  {                                                   \
    operator;                                         \
  }                                                   \
}

#define REFERENCE_FUNC(TYPE) TYPE ## _by_reference
#define VALUE_FUNC(TYPE) TYPE ## _by_value

#define GEN_COREDUCE(name, dt)                                \
static void                                                   \
name##_by_reference_adapter(void *invec, void *inoutvec,      \
                            int *len, MPI_Datatype *datatype) \
{                                                             \
  for (int i = 0; i < *len; ++i)                              \
  {                                                           \
    *((dt*)inoutvec) =                                        \
      (dt)(REFERENCE_FUNC(dt)((dt *)invec, (dt *)inoutvec));  \
    invec += sizeof(dt);                                      \
    inoutvec += sizeof(dt);                                   \
  }                                                           \
}                                                             \
static void                                                   \
name##_by_value_adapter(void *invec, void *inoutvec,          \
                        int *len, MPI_Datatype *datatype)     \
{                                                             \
  for (int i = 0; i < *len; ++i)                              \
  {                                                           \
    *((dt*)inoutvec) =                                        \
      (dt)(VALUE_FUNC(dt)(*(dt *)invec, *(dt *)inoutvec));    \
    invec += sizeof(dt);                                      \
    inoutvec += sizeof(dt);                                   \
  }                                                           \
}

GEN_COREDUCE(redux_int8, int8_t)
GEN_COREDUCE(redux_int16, int16_t)
GEN_COREDUCE(redux_int32, int32_t)
GEN_COREDUCE(redux_int64, int64_t)
GEN_COREDUCE(redux_real32, float)
GEN_COREDUCE(redux_real64, double)

static void
redux_char_by_reference_adapter(void *invec, void *inoutvec, int *len,
                                MPI_Datatype *datatype)
{
  MPI_Aint lb, string_len;
  MPI_Type_get_extent(*datatype, &lb, &string_len);
  for (int i = 0; i < *len; i++)
  {
    /* The length of the result is fixed, i.e., no deferred string length is
     * allowed there. */
    REFERENCE_FUNC(char)(
      (char *)inoutvec, string_len, (char *)invec,
      (char *)inoutvec, string_len, string_len
    );
    invec += sizeof(char) * string_len;
    inoutvec += sizeof(char) * string_len;
  }
}

#ifndef MPI_INTEGER1
GEN_REDUCTION(do_sum_int1, int8_t, inoutvec[i] += invec[i])
GEN_REDUCTION(do_min_int1, int8_t,
              inoutvec[i] = (invec[i] >= inoutvec[i] ? inoutvec[i] : invec[i]))
GEN_REDUCTION(do_max_int1, int8_t,
              inoutvec[i] = (invec[i] <= inoutvec[i] ? inoutvec[i] : invec[i]))
#endif

/*
#ifndef MPI_INTEGER2 
GEN_REDUCTION(do_sum_int1, int16_t, inoutvec[i] += invec[i]) 
GEN_REDUCTION(do_min_int1, int16_t,
              inoutvec[i] = (invec[i] >= inoutvec[i] ? inoutvec[i] : invec[i])) 
GEN_REDUCTION(do_max_int1, int16_t, 
             inoutvec[i] = (invec[i] <= inoutvec[i] ? inoutvec[i] : invec[i]))
#endif
*/

#if defined(MPI_INTEGER16) && defined(GFC_INTEGER_16)
GEN_REDUCTION(do_sum_int1, GFC_INTEGER_16, inoutvec[i] += invec[i])
GEN_REDUCTION(do_min_int1, GFC_INTEGER_16,
              inoutvec[i] = (invec[i] >= inoutvec[i] ? inoutvec[i] : invec[i]))
GEN_REDUCTION(do_max_int1, GFC_INTEGER_16,
              inoutvec[i] = (invec[i] <= inoutvec[i] ? inoutvec[i] : invec[i]))
#endif

#if defined(GFC_DTYPE_REAL_10) \
    || (!defined(GFC_DTYPE_REAL_10) && defined(GFC_DTYPE_REAL_16))
GEN_REDUCTION(do_sum_real10, long double, inoutvec[i] += invec[i])
GEN_REDUCTION(do_min_real10, long double,
              inoutvec[i] = (invec[i] >= inoutvec[i] ? inoutvec[i] : invec[i]))
GEN_REDUCTION(do_max_real10, long double,
              inoutvec[i] = (invec[i] <= inoutvec[i] ? inoutvec[i] : invec[i]))
GEN_REDUCTION(do_sum_complex10, _Complex long double, inoutvec[i] += invec[i])
GEN_REDUCTION(do_min_complex10, _Complex long double,
              inoutvec[i] = (invec[i] >= inoutvec[i] ? inoutvec[i] : invec[i]))
GEN_REDUCTION(do_max_complex10, _Complex long double,
              inoutvec[i] = (invec[i] <= inoutvec[i] ? inoutvec[i] : invec[i]))
#endif

#if defined(GFC_DTYPE_REAL_10) && defined(GFC_DTYPE_REAL_16)
GEN_REDUCTION(do_sum_real10, __float128, inoutvec[i] += invec[i])
GEN_REDUCTION(do_min_real10, __float128,
              inoutvec[i] = (invec[i] >= inoutvec[i] ? inoutvec[i] : invec[i]))
GEN_REDUCTION(do_max_real10, __float128,
              inoutvec[i] = (invec[i] <= inoutvec[i] ? inoutvec[i] : invec[i]))
GEN_REDUCTION(do_sum_complex10, _Complex __float128, inoutvec[i] += invec[i])
GEN_REDUCTION(do_mincomplexl10, _Complex __float128,
              inoutvec[i] = (invec[i] >= inoutvec[i] ? inoutvec[i] : invec[i]))
GEN_REDUCTION(do_max_complex10, _Complex __float128,
              inoutvec[i] = (invec[i] <= inoutvec[i] ? inoutvec[i] : invec[i]))
#endif
#undef GEN_REDUCTION


static MPI_Datatype
get_MPI_datatype(gfc_descriptor_t *desc, int char_len)
{
  int ierr;
  /* FIXME: Better check whether the sizes are okay and supported;
   * MPI3 adds more types, e.g. MPI_INTEGER1. */
  switch (GFC_DTYPE_TYPE_SIZE(desc))
  {
#ifdef MPI_INTEGER1
    case GFC_DTYPE_INTEGER_1:
      return MPI_INTEGER1;
#endif
#ifdef MPI_INTEGER2
    case GFC_DTYPE_INTEGER_2:
      return MPI_INTEGER2;
#endif
    case GFC_DTYPE_INTEGER_4:
#ifdef MPI_INTEGER4
      return MPI_INTEGER4;
#else
      return MPI_INTEGER;
#endif
#ifdef MPI_INTEGER8
    case GFC_DTYPE_INTEGER_8:
      return MPI_INTEGER8;
#endif
#if defined(MPI_INTEGER16) && defined(GFC_DTYPE_INTEGER_16)
    case GFC_DTYPE_INTEGER_16:
      return MPI_INTEGER16;
#endif

    case GFC_DTYPE_LOGICAL_4:
      return MPI_INT;

    case GFC_DTYPE_REAL_4:
#ifdef MPI_REAL4
      return MPI_REAL4;
#else
      return MPI_REAL;
#endif
    case GFC_DTYPE_REAL_8:
#ifdef MPI_REAL8
      return MPI_REAL8;
#else
      return MPI_DOUBLE_PRECISION;
#endif

/* Note that we cannot use REAL_16 as we do not know whether it matches REAL(10)
 * or REAL(16), which have both the same bitsize and only make use of less
 * bits. */
    case GFC_DTYPE_COMPLEX_4:
      return MPI_COMPLEX;
    case GFC_DTYPE_COMPLEX_8:
      return MPI_DOUBLE_COMPLEX;
  }
/* gfortran passes character string arguments with a
 * GFC_DTYPE_TYPE_SIZE == GFC_TYPE_CHARACTER + 64*strlen */
  if ((GFC_DTYPE_TYPE_SIZE(desc) - GFC_DTYPE_CHARACTER) % 64 == 0)
  {
    MPI_Datatype string;

    if (char_len == 0)
      char_len = GFC_DESCRIPTOR_SIZE(desc);
    ierr = MPI_Type_contiguous(char_len, MPI_CHARACTER, &string); chk_err(ierr);
    ierr = MPI_Type_commit(&string); chk_err(ierr);
    return string;
  }

  return MPI_BYTE;
  /* caf_runtime_error("Unsupported data type in collective: %zd\n", */
  /*                   GFC_DTYPE_TYPE_SIZE(desc)); */
  /* return 0; */
}


static void
internal_co_reduce(MPI_Op op, gfc_descriptor_t *source, int result_image,
                   int *stat, char *errmsg, int src_len, size_t errmsg_len)
{
  size_t i, size;
  int j, ierr, rank = GFC_DESCRIPTOR_RANK(source);
  ptrdiff_t dimextent;

  MPI_Datatype datatype = get_MPI_datatype(source, src_len);

  size = 1;
  for (j = 0; j < rank; ++j)
  {
    dimextent = source->dim[j]._ubound - source->dim[j].lower_bound + 1;
    if (dimextent < 0)
      dimextent = 0;
    size *= dimextent;
  }

  if (rank == 0 || PREFIX(is_contiguous) (source))
  {
    if (result_image == 0)
    {
      ierr = MPI_Allreduce(MPI_IN_PLACE, source->base_addr, size, datatype, op,
                           CAF_COMM_WORLD); chk_err(ierr);
    }
    else if (result_image == caf_this_image)
    {
      ierr = MPI_Reduce(MPI_IN_PLACE, source->base_addr, size, datatype, op,
                        result_image - 1, CAF_COMM_WORLD); chk_err(ierr);
    }
    else
    {
      ierr = MPI_Reduce(source->base_addr, NULL, size, datatype, op,
                        result_image - 1, CAF_COMM_WORLD); chk_err(ierr);
    }
    if (ierr)
      goto error;
    goto co_reduce_cleanup;
  }

  for (i = 0; i < size; ++i)
  {
    ptrdiff_t array_offset_sr = 0, tot_ext = 1, extent = 1;
    for (j = 0; j < rank - 1; ++j)
    {
      extent = source->dim[j]._ubound - source->dim[j].lower_bound + 1;
      array_offset_sr += ((i / tot_ext) % extent) * source->dim[j]._stride;
      tot_ext *= extent;
    }
    array_offset_sr += (i / tot_ext) * source->dim[rank - 1]._stride;
    void *sr = (void *)((char *)source->base_addr
                        + array_offset_sr * GFC_DESCRIPTOR_SIZE(source));
    if (result_image == 0)
    {
      ierr = MPI_Allreduce(MPI_IN_PLACE, sr, 1, datatype, op, CAF_COMM_WORLD);
      chk_err(ierr);
    }
    else if (result_image == caf_this_image)
    {
      ierr = MPI_Reduce(MPI_IN_PLACE, sr, 1, datatype, op, result_image - 1,
                        CAF_COMM_WORLD); chk_err(ierr);
    }
    else
    {
      ierr = MPI_Reduce(sr, NULL, 1, datatype, op, result_image - 1,
                        CAF_COMM_WORLD); chk_err(ierr);
    }
    if (ierr)
      goto error;
  }

co_reduce_cleanup:
  if (GFC_DESCRIPTOR_TYPE(source) == BT_CHARACTER)
  {
    ierr = MPI_Type_free(&datatype); chk_err(ierr);
  }
  if (stat)
    *stat = 0;
  return;
error:
  /* FIXME: Put this in an extra function and use it elsewhere. */
  if (stat)
  {
    *stat = ierr;
    if (!errmsg)
      return;
  }

  int len = sizeof(err_buffer);
  MPI_Error_string(ierr, err_buffer, &len);
  if (!stat)
  {
    err_buffer[len == sizeof(err_buffer) ? len - 1 : len] = '\0';
    caf_runtime_error("CO_SUM failed with %s\n", err_buffer);
  }
  memcpy(errmsg, err_buffer, (errmsg_len > len) ? len : errmsg_len);
  if (errmsg_len > len)
    memset(&errmsg[len], '\0', errmsg_len - len);
}

void
PREFIX(co_broadcast) (gfc_descriptor_t *a, int source_image, int *stat,
                      char *errmsg, charlen_t errmsg_len)
{
  size_t i, size;
  int j, ierr, rank = GFC_DESCRIPTOR_RANK(a);
  ptrdiff_t dimextent;

  MPI_Datatype datatype = get_MPI_datatype(a, 0);

  size = 1;
  for (j = 0; j < rank; ++j)
  {
    dimextent = a->dim[j]._ubound - a->dim[j].lower_bound + 1;
    if (dimextent < 0)
      dimextent = 0;
    size *= dimextent;
  }

  if (rank == 0)
  {
    if( datatype == MPI_BYTE)
      {
	ierr = MPI_Bcast(a->base_addr, size*GFC_DESCRIPTOR_SIZE(a),
			 datatype, source_image - 1,
			 CAF_COMM_WORLD); chk_err(ierr);
      }
    else if (datatype != MPI_CHARACTER)
    {
      ierr = MPI_Bcast(a->base_addr, size, datatype, source_image - 1,
                       CAF_COMM_WORLD); chk_err(ierr);
    }
    else
    {
      int a_length;
      if (caf_this_image == source_image)
        a_length = strlen(a->base_addr);
      /* Broadcast the string lenth */
      ierr = MPI_Bcast(&a_length, 1, MPI_INT, source_image - 1, CAF_COMM_WORLD);
      chk_err(ierr);
      if (ierr)
        goto error;
      /* Broadcast the string itself */
      ierr = MPI_Bcast(a->base_addr, a_length, datatype, source_image - 1,
                       CAF_COMM_WORLD); chk_err(ierr);
    }

    if (ierr)
      goto error;
    goto co_broadcast_exit;
  }
  else if (datatype == MPI_CHARACTER) /* rank !=0 */
  {
      caf_runtime_error("Co_broadcast of character arrays are "
                        "not yet supported\n");
  }

  for (i = 0; i < size; ++i)
  {
    ptrdiff_t array_offset_sr = 0, tot_ext = 1, extent = 1;
    for (j = 0; j < rank - 1; ++j)
    {
      extent = a->dim[j]._ubound - a->dim[j].lower_bound + 1;
      array_offset_sr += ((i / tot_ext) % extent) * a->dim[j]._stride;
      tot_ext *= extent;
    }
    array_offset_sr += (i / tot_ext) * a->dim[rank - 1]._stride;
    void *sr = (void *)(
      (char *)a->base_addr + array_offset_sr * GFC_DESCRIPTOR_SIZE(a));

    ierr = MPI_Bcast(sr, 1, datatype, source_image - 1, CAF_COMM_WORLD);
    chk_err(ierr);

    if (ierr)
      goto error;
  }

co_broadcast_exit:
  if (stat)
    *stat = 0;
  if (GFC_DESCRIPTOR_TYPE(a) == BT_CHARACTER)
  {
    ierr = MPI_Type_free(&datatype); chk_err(ierr);
  }
  return;

error:
  /* FIXME: Put this in an extra function and use it elsewhere. */
  if (stat)
  {
    *stat = ierr;
    if (!errmsg)
      return;
  }

  int len = sizeof(err_buffer);
  MPI_Error_string(ierr, err_buffer, &len);
  if (!stat)
  {
    err_buffer[len == sizeof(err_buffer) ? len - 1 : len] = '\0';
    caf_runtime_error("CO_SUM failed with %s\n", err_buffer);
  }
  memcpy(errmsg, err_buffer, (errmsg_len > len) ? len : errmsg_len);
  if (errmsg_len > len)
    memset(&errmsg[len], '\0', errmsg_len - len);
}

/* The front-end function for co_reduce functionality.  It sets up the MPI_Op
 * for use in MPI_*Reduce functions. */
void
PREFIX(co_reduce) (gfc_descriptor_t *a, void *(*opr) (void *, void *),
                   int opr_flags, int result_image, int *stat, char *errmsg,
                   int a_len, charlen_t errmsg_len)
{
  MPI_Op op;
  int type_a = GFC_DESCRIPTOR_TYPE(a), ierr;
  /* Integers and logicals can be treated the same. */
  if (type_a == BT_INTEGER || type_a == BT_LOGICAL)
  {
    /* When the ARG_VALUE opr_flag is set, then the user-function expects its
     * arguments to be passed by value. */
    if ((opr_flags & GFC_CAF_ARG_VALUE) > 0)
    {
#define ifTypeGen(type)                                                   \
if (GFC_DESCRIPTOR_SIZE(a) == sizeof(type ## _t))                         \
{                                                                         \
  type ## _t_by_value = (typeof(VALUE_FUNC(type ## _t)))opr;              \
  int ierr = MPI_Op_create(redux_ ## type ## _by_value_adapter, 1, &op);  \
  chk_err(ierr);                                                          \
}
      ifTypeGen(int8)
      else ifTypeGen(int16)
      else ifTypeGen(int32)
      else ifTypeGen(int64)
      else
      {
        caf_runtime_error("CO_REDUCE unsupported integer datatype");
      }
#undef ifTypeGen
    }
    else
    {
      int32_t_by_reference = (typeof(REFERENCE_FUNC(int32_t)))opr;
      ierr = MPI_Op_create(redux_int32_by_reference_adapter, 1, &op);
      chk_err(ierr);
    }
  }
  /* Treat reals/doubles. */
  else if (type_a == BT_REAL)
  {
    /* When the ARG_VALUE opr_flag is set, then the user-function expects its
     * arguments to be passed by value. */
    if (GFC_DESCRIPTOR_SIZE(a) == sizeof(float))
    {
      if ((opr_flags & GFC_CAF_ARG_VALUE) > 0)
      {
        float_by_value = (typeof(VALUE_FUNC(float)))opr;
        ierr = MPI_Op_create(redux_real32_by_value_adapter, 1, &op);
        chk_err(ierr);
      }
      else
      {
        float_by_reference = (typeof(REFERENCE_FUNC(float)))opr;
        ierr = MPI_Op_create(redux_real32_by_reference_adapter, 1, &op);
        chk_err(ierr);
      }
    }
    else
    {
      /* When the ARG_VALUE opr_flag is set, then the user-function expects
       * its arguments to be passed by value. */
      if ((opr_flags & GFC_CAF_ARG_VALUE) > 0)
      {
        double_by_value = (typeof(VALUE_FUNC(double)))opr;
        ierr = MPI_Op_create(redux_real64_by_value_adapter, 1, &op);
        chk_err(ierr);
      }
      else
      {
        double_by_reference = (typeof(REFERENCE_FUNC(double)))opr;
        ierr = MPI_Op_create(redux_real64_by_reference_adapter, 1, &op);
        chk_err(ierr);
      }
    }
  }
  else if (type_a == BT_CHARACTER)
  {
    /* Char array functions always pass by reference. */
    char_by_reference = (typeof(REFERENCE_FUNC(char)))opr;
    ierr = MPI_Op_create(redux_char_by_reference_adapter, 1, &op);
    chk_err(ierr);
  }
  else
  {
    caf_runtime_error("Data type not yet supported for co_reduce\n");
  }

  internal_co_reduce(op, a, result_image, stat, errmsg, a_len, errmsg_len);
}

void
PREFIX(co_sum) (gfc_descriptor_t *a, int result_image, int *stat, char *errmsg,
                charlen_t errmsg_len)
{
  internal_co_reduce(MPI_SUM, a, result_image, stat, errmsg, 0, errmsg_len);
}


void
PREFIX(co_min) (gfc_descriptor_t *a, int result_image, int *stat, char *errmsg,
                int src_len, charlen_t errmsg_len)
{
  internal_co_reduce(MPI_MIN, a, result_image, stat, errmsg, src_len,
                     errmsg_len);
}


void
PREFIX(co_max) (gfc_descriptor_t *a, int result_image, int *stat,
                char *errmsg, int src_len, charlen_t errmsg_len)
{
  internal_co_reduce(MPI_MAX, a, result_image, stat, errmsg, src_len,
                     errmsg_len);
}


/* Locking functions */

void
PREFIX(lock) (caf_token_t token, size_t index, int image_index,
              int *acquired_lock, int *stat, char *errmsg,
              charlen_t errmsg_len)
{
  MPI_Win *p = TOKEN(token);
  mutex_lock(*p, (image_index == 0) ? caf_this_image : image_index,
             index, stat, acquired_lock, errmsg, errmsg_len);
}


void
PREFIX(unlock) (caf_token_t token, size_t index, int image_index,
                int *stat, char *errmsg, charlen_t errmsg_len)
{
  MPI_Win *p = TOKEN(token);
  mutex_unlock(*p, (image_index == 0) ? caf_this_image : image_index,
               index, stat, errmsg, errmsg_len);
}

/* Atomics operations */

void
PREFIX(atomic_define) (caf_token_t token, size_t offset,
                       int image_index, void *value, int *stat,
                       int type __attribute__((unused)), int kind)
{
  MPI_Win *p = TOKEN(token);
  MPI_Datatype dt;
  int ierr = 0,
      image = (image_index != 0) ? image_index - 1 : caf_this_image - 1;

  selectType(kind, &dt);

#if MPI_VERSION >= 3
  CAF_Win_lock(MPI_LOCK_EXCLUSIVE, image, *p);
  ierr = MPI_Accumulate(value, 1, dt, image, offset, 1, dt, MPI_REPLACE, *p);
  chk_err(ierr);
  CAF_Win_unlock(image, *p);
#else // MPI_VERSION
  CAF_Win_lock(MPI_LOCK_EXCLUSIVE, image, *p);
  ierr = MPI_Put(value, 1, dt, image, offset, 1, dt, *p); chk_err(ierr);
  CAF_Win_unlock(image, *p);
#endif // MPI_VERSION

  if (stat)
    *stat = ierr;
  else if (ierr != 0)
    terminate_internal(ierr, 0);

  return;
}

void
PREFIX(atomic_ref) (caf_token_t token, size_t offset,
                    int image_index,
                    void *value, int *stat,
                    int type __attribute__((unused)), int kind)
{
  MPI_Win *p = TOKEN(token);
  MPI_Datatype dt;
  int ierr = 0, 
      image = (image_index != 0) ? image_index - 1 : caf_this_image - 1;

  selectType(kind, &dt);

#if MPI_VERSION >= 3
  CAF_Win_lock(MPI_LOCK_EXCLUSIVE, image, *p);
  ierr = MPI_Fetch_and_op(NULL, value, dt, image, offset, MPI_NO_OP, *p);
  chk_err(ierr);
  CAF_Win_unlock(image, *p);
#else // MPI_VERSION
  CAF_Win_lock(MPI_LOCK_EXCLUSIVE, image, *p);
  ierr = MPI_Get(value, 1, dt, image, offset, 1, dt, *p); chk_err(ierr);
  CAF_Win_unlock(image, *p);
#endif // MPI_VERSION

  if (stat)
    *stat = ierr;
  else if (ierr != 0)
    terminate_internal(ierr, 0);

  return;
}

void
PREFIX(atomic_cas) (caf_token_t token, size_t offset, int image_index,
                    void *old, void *compare, void *new_val, int *stat,
                    int type __attribute__((unused)), int kind)
{
  MPI_Win *p = TOKEN(token);
  MPI_Datatype dt;
  int ierr = 0,
      image = (image_index != 0) ? image_index - 1 : caf_this_image - 1;

  selectType(kind, &dt);

#if MPI_VERSION >= 3
  CAF_Win_lock(MPI_LOCK_EXCLUSIVE, image, *p);
  ierr = MPI_Compare_and_swap(new_val, compare, old, dt, image, offset, *p);
  chk_err(ierr);
  CAF_Win_unlock(image, *p);
#else // MPI_VERSION
#warning atomic_cas for MPI-2 is not yet implemented
  printf("We apologize but atomic_cas for MPI-2 is not yet implemented\n");
  ierr = 1;
#endif // MPI_VERSION

  if (stat)
    *stat = ierr;
  else if (ierr != 0)
    terminate_internal(ierr, 0);

  return;
}

void
PREFIX(atomic_op) (int op, caf_token_t token, size_t offset, int image_index,
                   void *value, void *old, int *stat,
                   int type __attribute__((unused)), int kind)
{
  int ierr = 0;
  MPI_Datatype dt;
  MPI_Win *p = TOKEN(token);
  int image = (image_index != 0) ? image_index - 1 : caf_this_image - 1;

#if MPI_VERSION >= 3
  old = malloc(kind);
  selectType(kind, &dt);

  CAF_Win_lock(MPI_LOCK_EXCLUSIVE, image, *p);
  /* Atomic_add */
  switch(op) {
    case 1:
      ierr = MPI_Fetch_and_op(value, old, dt, image, offset, MPI_SUM, *p);
      chk_err(ierr);
      break;
    case 2:
      ierr = MPI_Fetch_and_op(value, old, dt, image, offset, MPI_BAND, *p);
      chk_err(ierr);
      break;
    case 4:
      ierr = MPI_Fetch_and_op(value, old, dt, image, offset, MPI_BOR, *p);
      chk_err(ierr);
      break;
    case 5:
      ierr = MPI_Fetch_and_op(value, old, dt, image, offset, MPI_BXOR, *p);
      chk_err(ierr);
      break;
    default:
      printf("We apologize but the atomic operation requested for MPI < 3 "
             "is not yet implemented\n");
      break;
    }
  CAF_Win_unlock(image, *p);

  free(old);
#else // MPI_VERSION
  #warning atomic_op for MPI is not yet implemented
  printf("We apologize but atomic_op for MPI < 3 is not yet implemented\n");
#endif // MPI_VERSION
  if (stat)
    *stat = ierr;
  else if (ierr != 0)
    terminate_internal(ierr, 0);

  return;
}

/* Events */

void
PREFIX(event_post) (caf_token_t token, size_t index, int image_index,
                    int *stat, char *errmsg, charlen_t errmsg_len)
{
  int value = 1, ierr = 0, flag;
  MPI_Win *p = TOKEN(token);
  const char msg[] = "Error on event post";
  int image = (image_index == 0) ? caf_this_image - 1 : image_index - 1;

  if (stat != NULL)
    *stat = 0;

#if MPI_VERSION >= 3
  CAF_Win_lock(MPI_LOCK_EXCLUSIVE, image, *p);
  ierr = MPI_Accumulate(&value, 1, MPI_INT, image, index * sizeof(int), 1, 
                        MPI_INT, MPI_SUM, *p); chk_err(ierr);
  CAF_Win_unlock(image, *p);
#else // MPI_VERSION
  #warning Events for MPI-2 are not implemented
  printf("Events for MPI-2 are not supported, "
         "please update your MPI implementation\n");
#endif // MPI_VERSION

  check_image_health(image_index, stat);

  if (!stat && ierr == STAT_FAILED_IMAGE)
    terminate_internal(ierr, 0);

  if (ierr != MPI_SUCCESS)
  {
    if (stat != NULL)
      *stat = ierr;
    if (errmsg != NULL)
    {
      memset(errmsg,' ',errmsg_len);
      memcpy(errmsg, msg, MIN(errmsg_len,strlen(msg)));
    }
  }
}

void
PREFIX(event_wait) (caf_token_t token, size_t index, int until_count,
                    int *stat, char *errmsg, charlen_t errmsg_len)
{
  int ierr = 0, count = 0, i, image = caf_this_image - 1;
  int *var = NULL, flag, old = 0, newval = 0;
  const int spin_loop_max = 20000;
  MPI_Win *p = TOKEN(token);
  const char msg[] = "Error on event wait";

  if (stat != NULL)
    *stat = 0;

  ierr = MPI_Win_get_attr(*p, MPI_WIN_BASE, &var, &flag); chk_err(ierr);

  for (i = 0; i < spin_loop_max; ++i)
  {
    ierr = MPI_Win_sync(*p); chk_err(ierr);
    count = var[index];
    if (count >= until_count)
      break;
  }

  i = 1;
  while (count < until_count)
  {
    ierr = MPI_Win_sync(*p); chk_err(ierr);
    count = var[index];
    usleep(10 * i);
    ++i;
    /* Needed to enforce MPI progress */
    ierr = MPI_Win_flush(image, *p); chk_err(ierr);
  }

  newval = -until_count;

  CAF_Win_lock(MPI_LOCK_EXCLUSIVE, image, *p);
  ierr = MPI_Fetch_and_op(&newval, &old, MPI_INT, image, index * sizeof(int),
                          MPI_SUM, *p); chk_err(ierr);
  CAF_Win_unlock(image, *p);
  check_image_health(image, stat);

  if (!stat && ierr == STAT_FAILED_IMAGE)
    terminate_internal(ierr, 0);

  if (ierr != MPI_SUCCESS)
  {
    if (stat != NULL)
      *stat = ierr;
    if (errmsg != NULL)
    {
      memset(errmsg,' ',errmsg_len);
      memcpy(errmsg, msg, MIN(errmsg_len,strlen(msg)));
    }
  }
}

void
PREFIX(event_query) (caf_token_t token, size_t index,
                     int image_index, int *count, int *stat)
{
  MPI_Win *p = TOKEN(token);
  int ierr = 0,
      image = (image_index == 0) ? caf_this_image - 1 : image_index - 1;

  if (stat != NULL)
    *stat = 0;

#if MPI_VERSION >= 3
  CAF_Win_lock(MPI_LOCK_EXCLUSIVE, image, *p);
  ierr = MPI_Fetch_and_op(NULL, count, MPI_INT, image, index * sizeof(int),
                          MPI_NO_OP, *p); chk_err(ierr);
  CAF_Win_unlock(image, *p);
#else // MPI_VERSION
#warning Events for MPI-2 are not implemented
  printf("Events for MPI-2 are not supported, "
         "please update your MPI implementation\n");
#endif // MPI_VERSION
  if (ierr != MPI_SUCCESS && stat != NULL)
    *stat = ierr;
}


/* Internal function to execute the part that is common to all (error) stop
 * functions. */

static void
terminate_internal(int stat_code, int exit_code)
{
  dprint("terminate_internal (stat_code = %d, exit_code = %d).\n",
         stat_code, exit_code);
  finalize_internal(stat_code);

#ifndef WITH_FAILED_IMAGES
  MPI_Abort(MPI_COMM_WORLD, exit_code);
#endif
  exit(exit_code);
}


#ifdef GCC_GE_8
#undef QUIETARG
#define QUIETARG , bool quiet
#endif

/* STOP function for integer arguments. */

void
PREFIX(stop_numeric) (int stop_code QUIETARG)
{
#ifndef GCC_GE_8
  bool quiet = false;
#endif
  if (!quiet)
    fprintf(stderr, "STOP %d\n", stop_code);

  /* Stopping includes taking down the runtime regularly and returning the
   * stop_code. */
  terminate_internal(STAT_STOPPED_IMAGE, stop_code);
}


/* STOP function for string arguments. */

void
PREFIX(stop_str) (const char *string, charlen_t len QUIETARG)
{
#ifndef GCC_GE_8
  bool quiet = false;
#endif
  if (!quiet)
  {
    fputs("STOP ", stderr);
    while (len--)
      fputc(*(string++), stderr);
    fputs("\n", stderr);
  }
  /* Stopping includes taking down the runtime regularly. */
  terminate_internal(STAT_STOPPED_IMAGE, 0);
}


/* ERROR STOP function for string arguments. */

static void
error_stop_str(const char *string, size_t len, bool quiet)
{
  if (!quiet)
  {
    fputs("ERROR STOP ", stderr);
    while (len--)
      fputc(*(string++), stderr);
    fputs("\n", stderr);
  }
  terminate_internal(STAT_STOPPED_IMAGE, 1);
}


void
PREFIX(error_stop_str) (const char *string, charlen_t len QUIETARG)
{
#ifndef GCC_GE_8
  bool quiet = false;
#endif
  error_stop_str(string, len, quiet);
}


/* ERROR STOP function for numerical arguments. */

void
PREFIX(error_stop) (int error QUIETARG)
{
#ifndef GCC_GE_8
  bool quiet = false;
#endif
  if (!quiet)
    fprintf(stderr, "ERROR STOP %d\n", error);

  terminate_internal(STAT_STOPPED_IMAGE, error);
}


/* FAIL IMAGE statement. */

void
PREFIX(fail_image) (void)
{
  fputs("IMAGE FAILED!\n", stderr);

  raise(SIGKILL);
  /* A failing image is expected to take down the runtime regularly. */
  terminate_internal(STAT_FAILED_IMAGE, 0);
}

int
PREFIX(image_status) (int image)
{
#ifdef GFC_CAF_CHECK
  if (image < 1 || image > caf_num_images)
  {
    char errmsg[60];
    sprintf(errmsg, "Image #%d out of bounds of images 1..%d.",
            image, caf_num_images);
    caf_runtime_error(errmsg);
  }
#endif
#ifdef WITH_FAILED_IMAGES
  if (image_stati[image - 1] == 0)
  {
    int status, ierr;
    /* Check that we are fine before doing anything.
     *
     * Do an MPI-operation to learn about failed/stopped images, that have
     * not been detected yet. */
    ierr = MPI_Test(&alive_request, &status, MPI_STATUSES_IGNORE);
    chk_err(ierr);
    MPI_Error_class(ierr, &status);
    if (ierr == MPI_SUCCESS)
    {
      CAF_Win_lock(MPI_LOCK_SHARED, image - 1, *stat_tok);
      ierr = MPI_Get(&status, 1, MPI_INT, image - 1, 0, 1, MPI_INT, *stat_tok);
      chk_err(ierr);
      dprint("Image status of image #%d is: %d\n", image, status);
      CAF_Win_unlock(image - 1, *stat_tok);
      image_stati[image - 1] = status;
    }
    else if (status == MPIX_ERR_PROC_FAILED)
    {
      image_stati[image - 1] = STAT_FAILED_IMAGE;
    }
    else
    {
      const int strcap = 200;
      char errmsg[strcap];
      int slen, supplied_len;
      sprintf(errmsg, "Image status for image #%d returned mpi error: ",
              image);
      slen = strlen(errmsg);
      supplied_len = strcap - slen;
      MPI_Error_string(status, &errmsg[slen], &supplied_len);
      caf_runtime_error(errmsg);
    }
  }
  return image_stati[image - 1];
#else
  unsupported_fail_images_message("IMAGE_STATUS()");
#endif // WITH_FAILED_IMAGES

  return 0;
}

void
PREFIX(failed_images) (gfc_descriptor_t *array,
                       int team __attribute__((unused)), int * kind)
{
  int local_kind = kind ? *kind : 4; /* GFC_DEFAULT_INTEGER_KIND = 4*/

#ifdef WITH_FAILED_IMAGES
  void *mem = calloc(num_images_failed, local_kind);
  array->base_addr = mem;
  for (int i = 0; i < caf_num_images; ++i)
  {
    if (image_stati[i] == STAT_FAILED_IMAGE)
    {
      switch (local_kind)
      {
        case 1:
          *(int8_t *)mem = i + 1;
          break;
        case 2:
          *(int16_t *)mem = i + 1;
          break;
        case 4:
          *(int32_t *)mem = i + 1;
          break;
        case 8:
          *(int64_t *)mem = i + 1;
          break;
#ifdef HAVE_GFC_INTEGER_16
        case 16:
          *(int128t *)mem = i + 1;
          break;
#endif
        default:
          caf_runtime_error("Unsupported integer kind %1 "
                            "in caf_failed_images.", local_kind);
      }
      mem += local_kind;
    }
  }
  array->dim[0]._ubound = num_images_failed - 1;
#else
  unsupported_fail_images_message("FAILED_IMAGES()");
  array->dim[0]._ubound = -1;
  array->base_addr = NULL;
#endif // WITH_FAILED_IMAGES

#ifdef GCC_GE_8
  array->dtype.type = BT_INTEGER;
  array->dtype.elem_len = local_kind;
#else
  array->dtype = ((BT_INTEGER << GFC_DTYPE_TYPE_SHIFT)
                  | (local_kind << GFC_DTYPE_SIZE_SHIFT));
#endif
  array->dim[0].lower_bound = 0;
  array->dim[0]._stride = 1;
  array->offset = 0;
}

void
PREFIX(stopped_images) (gfc_descriptor_t *array,
                        int team __attribute__((unused)), int * kind)
{
  int local_kind = kind ? *kind : 4; /* GFC_DEFAULT_INTEGER_KIND = 4*/

#ifdef WITH_FAILED_IMAGES
  void *mem = calloc(num_images_stopped, local_kind);
  array->base_addr = mem;
  for (int i = 0; i < caf_num_images; ++i)
  {
    if (image_stati[i])
    {
      switch (local_kind)
      {
        case 1:
          *(int8_t *)mem = i + 1;
          break;
        case 2:
          *(int16_t *)mem = i + 1;
          break;
        case 4:
          *(int32_t *)mem = i + 1;
          break;
        case 8:
          *(int64_t *)mem = i + 1;
          break;
#ifdef HAVE_GFC_INTEGER_16
        case 16:
          *(int128t *)mem = i + 1;
          break;
#endif
        default:
          caf_runtime_error("Unsupported integer kind %1 "
                            "in caf_stopped_images.", local_kind);
      }
      mem += local_kind;
    }
  }
  array->dim[0]._ubound = num_images_stopped - 1;
#else
  unsupported_fail_images_message("STOPPED_IMAGES()");
  array->dim[0]._ubound = -1;
  array->base_addr = NULL;
#endif // WITH_FAILED_IMAGES

#ifdef GCC_GE_8
  array->dtype.type = BT_INTEGER;
  array->dtype.elem_len = local_kind;
#else
  array->dtype = ((BT_INTEGER << GFC_DTYPE_TYPE_SHIFT)
                  | (local_kind << GFC_DTYPE_SIZE_SHIFT));
#endif
  array->dim[0].lower_bound = 0;
  array->dim[0]._stride = 1;
  array->offset = 0;
}

/* Give a descriptive message when failed images support is not available. */
void
unsupported_fail_images_message(const char * functionname)
{
  fprintf(stderr,
          "*** caf_mpi-lib runtime message on image %d:\n"
          "*** The failed images feature '%s' "
          "*** of Fortran 2015 standard is not available in this build."
          "*** You need a compiler with failed images support activated and "
          "*** compile OpenCoarrays with failed images support.\n",
          caf_this_image, functionname);
#ifdef STOP_ON_UNSUPPORTED
  exit(EXIT_FAILURE);
#endif
}

/* Give a descriptive message when support for an allocatable components
 * feature is not available. */
void
unimplemented_alloc_comps_message(const char * functionname)
{
  fprintf(stderr,
          "*** Message from libcaf_mpi runtime function '%s' on image %d:\n"
          "*** Assigning to an allocatable coarray component of a derived type"
          "is not yet supported with GCC 7.\n"
          "*** Either revert to GCC 6 or convert all "
          "puts (type(foo)::x; x%%y[recipient] = z) to "
          "gets (z = x%%y[provider]).\n",
          functionname, caf_this_image );
#ifdef STOP_ON_UNSUPPORTED
  exit(EXIT_FAILURE);
#endif
}

void PREFIX(form_team) (int team_id, caf_team_t *team,
                        int index __attribute__((unused)))
{
  struct caf_teams_list *tmp;
  void * tmp_team;
  MPI_Comm *newcomm;
  MPI_Comm *current_comm = &CAF_COMM_WORLD;
  int ierr;

  newcomm = (MPI_Comm *)calloc(1,sizeof(MPI_Comm));
  ierr = MPI_Comm_split(*current_comm, team_id, caf_this_image, newcomm);
  chk_err(ierr);

  tmp = calloc(1,sizeof(struct caf_teams_list));
  tmp->prev = teams_list;
  teams_list = tmp;
  teams_list->team_id = team_id;
  teams_list->team = newcomm;
  *team = tmp;
}

void PREFIX(change_team) (caf_team_t *team,
                          int coselector __attribute__((unused)))
{
  caf_used_teams_list *tmp_used = NULL;
  caf_teams_list * tmp_list = NULL;
  void *tmp_team;
  MPI_Comm *tmp_comm;

  tmp_list = (struct caf_teams_list *)*team;
  tmp_team = (void *)tmp_list->team;
  tmp_comm = (MPI_Comm *)tmp_team;

  tmp_used = (caf_used_teams_list *)calloc(1,sizeof(caf_used_teams_list));
  tmp_used->prev = used_teams;

  /* We need to look in the teams_list and find the appropriate element. 
   * This is not efficient but can be easily fixed in the future.  
   * Instead of keeping track of the communicator in the compiler  
   * we should keep track of the caf_teams_list element associated with it. */ 

  /*
  tmp_list = teams_list; 

  while (tmp_list) 
  { 
    if (tmp_list->team == tmp_team) 
      break; 
    tmp_list = tmp_list->prev; 
  }
  */

  if (tmp_list == NULL)
    caf_runtime_error("CHANGE TEAM called on a non-existing team");

  tmp_used->team_list_elem = tmp_list;
  used_teams = tmp_used;
  tmp_team = tmp_used->team_list_elem->team;
  tmp_comm = (MPI_Comm *)tmp_team;
  CAF_COMM_WORLD = *tmp_comm;
  int ierr = MPI_Comm_rank(*tmp_comm,&caf_this_image); chk_err(ierr);
  caf_this_image++;
  ierr = MPI_Comm_size(*tmp_comm,&caf_num_images); chk_err(ierr);
  ierr = MPI_Barrier(*tmp_comm); chk_err(ierr);
}

MPI_Fint
PREFIX(get_communicator) (caf_team_t *team)
{
  if (team != NULL) caf_runtime_error("get_communicator does not yet support "
                                      "the optional team argument");

  MPI_Comm* comm_ptr = teams_list->team;
  MPI_Fint ret = MPI_Comm_c2f(*comm_ptr);

  return ret;
  // return  *(int*)comm_ptr;
}

int
PREFIX(team_number) (caf_team_t *team)
{
  if (team != NULL)
    return ((caf_teams_list *)team)->team_id;
  else
    return used_teams->team_list_elem->team_id; /* current team */
}

void PREFIX(end_team) (caf_team_t *team __attribute__((unused)))
{
  caf_used_teams_list *tmp_used = NULL;
  void *tmp_team;
  MPI_Comm *tmp_comm;
  int ierr;

  ierr = MPI_Barrier(CAF_COMM_WORLD); chk_err(ierr);
  if (used_teams->prev == NULL)
    caf_runtime_error("END TEAM called on initial team");

  tmp_used = used_teams;
  used_teams = used_teams->prev;
  free(tmp_used);
  tmp_used = used_teams;
  tmp_team = tmp_used->team_list_elem->team;
  tmp_comm = (MPI_Comm *)tmp_team;
  CAF_COMM_WORLD = *tmp_comm;
  /* CAF_COMM_WORLD = (MPI_Comm)*tmp_used->team_list_elem->team; */
  ierr = MPI_Comm_rank(CAF_COMM_WORLD,&caf_this_image); chk_err(ierr);
  caf_this_image++;
  ierr = MPI_Comm_size(CAF_COMM_WORLD,&caf_num_images); chk_err(ierr);
}

void PREFIX(sync_team) (caf_team_t *team , int unused __attribute__((unused)))
{
  caf_teams_list *tmp_list = NULL;
  caf_used_teams_list *tmp_used = NULL;
  void *tmp_team;
  MPI_Comm *tmp_comm;

  tmp_used = used_teams;
  tmp_list = (struct caf_teams_list *)*team;
  tmp_team = (void *)tmp_list->team;
  tmp_comm = (MPI_Comm *)tmp_team;

  /* if the team is not a child */
  if (tmp_used->team_list_elem != tmp_list->prev)
  /* then search backwards through the team list, first checking if it's the
   * current team, then if it is an ancestor team */
    while (tmp_used)
    {
      if (tmp_used->team_list_elem == tmp_list)
        break;
      tmp_used = tmp_used->prev;
    }

  if (tmp_used == NULL)
    caf_runtime_error("SYNC TEAM called on team different from current, "
                      "or ancestor, or child");

  int ierr = MPI_Barrier(*tmp_comm); chk_err(ierr);
}
