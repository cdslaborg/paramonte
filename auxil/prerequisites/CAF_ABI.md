[This document is formatted with GitHub-Flavored Markdown.              ]:#
[For better viewing, including hyperlinks, read it online at            ]:#
[https://github.com/sourceryinstitute/OpenCoarrays/blob/master/CAF_API.md]:#

OpenCoarrays Application Binary Interface (ABI)
===============================================

[![Download as PDF][pdf img]](http://md2pdf.herokuapp.com/sourceryinstitute/OpenCoarrays/blob/master/CAF_ABI.pdf)

Download this file as a PDF document
[here](http://md2pdf.herokuapp.com/sourceryinstitute/OpenCoarrays/blob/master/CAF_ABI.pdf).

* [To Do](#to-do)
* [Implementation status](#implementation-status)
* [Definitions and types](#definitions-and-types)
* [Provided functions](#provided-functions)

This document describes the OpenCoarrays application binary interface (ABI) through
which a compiler accesses coarray functionality.  As such, the target audience for
this document is compiler developers.  Most application developers need only write
standard-conforming Fortran 2008 or 2015 and compile their code with the OpenCoarrays
`caf` compiler wrapper without knowledge of the ABI.

The actual function names in this document have a PREFIX in the source code to avoid
name clashes.  The prefix can be vendor-specific.

### Warning ###

*This document may be out of date.*

To Do
-----

* [ ] Discuss the current draft
* [ ] Add missing functions of the current gfortran implementation
* [ ] Address the TODO items
* [ ] Extend the functions to match a sensible set
* [ ] Update the implementation status, especially for the ARMCI library

Implementation status
---------------------

The library implementation in this directory should be ABI-compatible
with the wording below, except for some `int errmsg_len` vs. `size_t`
changes that have not yet been implemented.

Definitions and types
---------------------

### 2.1  `caf_token_t` ###

Typedef of type `void *` on the compiler side. Can be any data
type on the library side.

### 2.2  `caf_register_t` ###

Type indicating which kind of coarray variable should be registered.

```c
typedef enum caf_register_t {
  CAF_REGTYPE_COARRAY_STATIC,
  CAF_REGTYPE_COARRAY_ALLOC,
  CAF_REGTYPE_LOCK_STATIC,
  CAF_REGTYPE_LOCK_ALLOC,
  CAF_REGTYPE_CRITICAL,
  CAF_REGTYPE_EVENT_STATIC,
  CAF_REGTYPE_EVENT_ALLOC
  }
caf_register_t;
```

__TODO__:
  Check whether this set is complete and makes sense


### 2.3  `caf_token_t` ###

In terms of the processor, an opaque pointer, which is used to identify a
coarray.  The exact content is implementation-defined by the library.

### 2.4  Stat values ###

```c
#define STAT_UNLOCKED           0
#define STAT_LOCKED             1
#define STAT_LOCKED_OTHER_IMAGE 2
#define STAT_STOPPED_IMAGE      6000
```

__TODO__:
  Define more, allow room for lib-specific values, update for [TS18508].
  Do we need to take care of special vendor choices?

__Note__:
  Some values have to be such that they differ from certain other
  values.


Provided functions
------------------

### 3.1  Initialization function ###

```c
void caf_init (int *argc, char ***argv)
```

This function shall be called at startup of the program before the Fortran main
program.  It takes as arguments the command-line arguments of the program. It is
permitted to pass to NULL pointers as argument; if non-NULL, the library is
permitted to modify the arguments.

| Argument	| `intent`	| description	|
| ------	| ------	| ------	|
| `argc`	| `inout`	| An integer pointer with the number of arguments passed to the program or NULL.	|
| `argv`	| `inout`	| A pointer to an array of strings with the      command-line arguments or NULL.	|

__Note__:
  The function is modeled after the initialization function of the
  Message Passing Interface (MPI) specification.  Due to the way coarray
  registration (3.5) works, it might not be the first call to the libaray. If
  the main program is not written in Fortran and only a library uses coarrays,
  it can happen that this function is never called.  Therefore, it is
  recommended that the library does not rely on the passed arguments and whether
  the call has been done.

__GCC__:
  In gfortran, the function is generated when the Fortran main program is
  compiled with -fcoarray=lib; the call happens before the run-time library
  initialiation such that changes to the command-line arguments will be visible
  when the command-line intrinsics are invoked.


### 3.2  Finalization function ###

```c
void caf_finish (void)
```

This function shall be called at the end of the program to permit a graceful
shutdown.

__Note__:
  It is recommended to add this call at the end of the Fortran main program and
  when invoking `STOP`.  To ensure that the shutdown is also performed for
  programs where this function is not explicitly invoked, for instance
  non-Fortran programs or calls to the system's `exit()` function, the library can
  use a destructor function.  Note that programs can also be terminated using
  the `ERROR STOP` statement, which is handled via its own library call.

__GCC__:
  In gfortran, this function is called at the end of the Fortran main program and
  when before the program stops with a `STOP` command, the respective file has been
  compiled with the `-fcoarray=lib` option.


### 3.3 Querying the image number ###

```c
int caf_this_image (int distance)
```

This function returns the current image number, which is a positive number.

| Argument	| description	|
| ------	| ------	|
| `distance`	| As specified for the `this_image` intrinsic in [TS18508]. Shall be a nonnegative number.	|

__Note__:
  If the Fortran intrinsic `this_image()` is invoked without an argument, which is the only permitted form in Fortran 2008, the processor shall pass 0 as first argument.

__GCC__:
  (No special note.)



### 3.4 Querying the maximal number of images ###

```c
int caf_num_images (int distance, int failed)
```

This function returns the number of images in the current team, if distance is 0
or the number of images in the parent team at the specified distance. If failed
is -1, the function returns the number of all images at the specified
distance; if it is 0, the function returns the number of non-failed images, and
if it is 1, it returns the number of failed images.

| Argument	| description	|
| ------	| ------	|
| `distance`	| the distance from this image to the ancestor. Shall be positive.	|
| `failed`	| shall be -1, 0, or 1	|

__Note__:
  This function follows [TS18508]. If the `num_image` intrinsic has no arguments,
  the processor shall pass `distance = 0` and `failed = -1` to the function.

__GCC__:
  (No special note.)



### 3.5 Registering coarrays ###

```c
void *caf_register (size_t size, caf_register_t type, caf_token_t *token, int *stat, char *errmsg, int errmsg_len)
```

Allocates memory for a coarray and creates a token to identify the coarray. The
function is called for both coarrays with `SAVE` attribute and using an explicit
`ALLOCATE` statement. If an error occurs and `STAT` is a `NULL` pointer, the function
shall abort with printing an error message and starting the error termination.
If no error occurs and `STAT=` is present, it shall be set to zero. Otherwise, it
shall be set to a positive value and, if not-@code{NULL}, @var{ERRMSG} shall be
set to a string describing the failure.  The function returns a pointer to the
requested memory for the local image as a call to `malloc` would do.

For `CAF_REGTYPE_COARRAY_STATIC` and `CAF_REGTYPE_COARRAY_ALLOC`, the passed size is
the byte size requested. For `CAF_REGTYPE_LOCK_STATIC`, `CAF_REGTYPE_LOCK_ALLOC`
and `CAF_REGTYPE_CRITICAL` it is the array size or one for a scalar.

| Argument	| description	|
| ------	| ------	|
| `size`	| For normal coarrays, the byte size of the coarray to be allocated; for lock types, the number of elements.	|
| `type`	| one of the `caf_register_t` types. Possible values: `CAF_REGTYPE_COARRAY_STATIC` - for nonallocatable coarrays `CAF_REGTYPE_COARRAY_ALLOC` - for allocatable coarrays `CAF_REGTYPE_LOCK_STATIC` - for nonallocatable lock variables `CAF_REGTYPE_LOCK_ALLOC` - for allocatable lock variables `CAF_REGTYPE_CRITICAL` - for lock variables used for critical sections	|
| `token`	| `intent(out)` An opaque pointer identifying the coarray.	|
| `stat`	| `intent(out)` For allocatable coarrays, stores the `STAT=`; may be `NULL`	|
| `errmsg`	| intent(out) When an error occurs, this will be set to an error message; may be `NULL`	|
| `errmgs_len`	| the buffer size of errmsg.	|

__TODO__:

  - [ ] Check whether the locking should be handled like that and whether one needs
  more, e.g. for locking types in DT?
  - [ ] Check whether one needs an additional function for to register coarrays
  which are in static memory and used without memory allocation, i.e. just to
  register the address.
  - [ ] Check whether we need an explicit `SYNC ALL` at the beginning of the main
  program or whether we can do without.
  - [ ] Does [TS18508] require more for `SAVE` within teams or within blocks?

__Note__:
  Non-allocatable coarrays have to be registered prior use from remote images.
  In order to guarantee this, they have to be registered before the main
  program. This can be achieved by creating constructor functions.  When using
  `caf_register`, also non-allocatable coarrays the memory is allocated and no
  static memory is used.

  For normal coarrays, the returned pointer is used for accesses on the local
  image. For lock types, the value shall only used for checking the allocation
  status. Note that for critical blocks, the locking is only required on one
  image; in the locking statement, the processor shall always pass always an
  image index of one for critical-section lock variables (`CAF_REGTYPE_CRITICAL`).

__GCC__:
   (no special notes)

__TODO__:
   Change `errmsg_len` to `size_t`



### 3.6  Deregistering coarrays ###

```c
void caf_deregister (const caf_token_t *token, int *stat, char *errmsg, size_t errmsg_len)
```

Called to free the memory of a coarray; the processor calls this function for
automatic and explicit deallocation.  In case of an error, this function shall
fail with an error message, unless the `STAT=` variable is not null.

| Argument	| `intent`	| description	|
| ------	| ------	| -----	|
| `token`	| `inout`	| An opaque pointer identifying the coarray.	|
| `stat`	| `out`	| For allocatable coarrays, stores the `STAT=`; may be `NULL`	|
| `errmsg`	| `out`	| When an error occurs, this will be set to an error message, may be `NULL`	|
| `errmgs_len`	| | the buffersize of `errmsg`.	|

__Note__:
  The implementation is permitted to set the token to `NULL`. However, it is not required to do so.
  For nonalloatable coarrays this function is never called.  If a cleanup is required, it has to be handled via the finish, stop and error stop functions, and via destructors.

__GCC__:
   (no special notes)

__TODO__:
   Change `errmsg_len` to `size_t`


### 3.7  Sending data from a local image to a remote image ###

```c
void caf_send (caf_token_t token, size_t offset, int image_index,
               gfc_descriptor_t *dest, caf_vector_t *dst_vector,
               gfc_descriptor_t *src, int dst_kind, int src_kind)
```

Called to send a scalar, an array section or whole array from a local
to a remote image identified by the `image_index`.

| Argument	| description	|
| ------	| ------	|
| `token`	| `intent(in)`  An opaque pointer identifying the coarray.	|
| `offset`	| By which amount of bytes the actual data is shifted compared to the base address of the coarray.	|
| `image_index`	| The ID of the remote image; must be a positive number.	|
| `dest`	| `intent(in)` Array descriptor for the remote image for the bounds and the size. The `base_addr` shall not be accessed.	|
| `dst_vector`	| `intent(in)`  If not `NULL`, it contains the vector subscript of the destination array; the values are relative to the dimension triplet of the dest argument.	|
| `src`	| `intent(in)` Array descriptor of the local array to be transferred to the remote image	|
| `dst_kind`	| Kind of the destination argument	|
| `src_kind`	| Kind of the source argument	|

__Note__:
  It is permitted to have `image_id` equal the current image; the memory of the
  send-to and the send-from might (partially) overlap in that case. The
  implementation has to take care that it handles this case. Note that the
  assignment of a scalar to an array is permitted. In addition, the library has
  to handle numeric-type conversion and for strings, padding and different
  `character` kinds.

__GCC__:
  Currently, it uses gfortran's private array descriptor. A change to [TS29113]'s
  array descriptor is planned; when that's done, the additional kind arguments
  will be removed.
  Note that the kind arguments permit to distiniguish the `character` kinds and
  `real`/`complex` kinds 10 and 16, which have the same byte size.


__TODO__ `FOR SEND*`:

  - [ ] Wait is missing
  - [ ] Assignment to an address instead of using a token, to handle
    `caf[i]%allocatable%alloc_array(:,:) = ...`
    Or some other means to handle those.
  - [ ] Image index: How to handle references to other TEAMS?

__OTHER TODOs__:

- [ ] 3.x TODO: Handle `GET` and remote-to-remote communication
- [ ] 3.y TODO: Handle `ATOMIC`, `LOCK`, `CRITICAL`
- [ ] 3.z TODO Teams and error recovery

### 3.8  Getting data from a remote image ###

```c
void caf_get_desc (caf_token_t token, size_t offset,
                   int image_index, gfc_descriptor_t *src,
                   caf_vector_t *src_vector, gfc_descriptor_t *dest,
                   int src_kind, int dst_kind)
```

Called to get an array section or whole array from a a remote,
image identified by the `image_index`.

| Argument	| description	|
| ------	| ------	|
| `token`	| `intent(in)` An opaque pointer identifying the coarray.	|
| `offset`	| By which amount of bytes the actual data is shifted compared to the base address of the coarray.	|
| `image_index`	| The ID of the remote image; must be a positive number.	|
| `dest`	| `intent(out)` Array descriptor of the local array to which the data will be transferred	|
| `src`	| `intent(in)` Array descriptor for the remote image for the bounds and the size. The `base_addr` shall not be accessed.
| `src_vector`	| `intent(int)` If not `NULL`, it contains the vector subscript of the destination array; the values are relative to the dimension triplet of the dest argument.	|
| `dst_kind`	| Kind of the destination argument	|
| `src_kind`	| Kind of the source argument	|

__Note__:
  It is permitted to have `image_id` equal the current image; the memory of the
  send-to and the send-from might (partially) overlap in that case. The
  implementation has to take care that it handles this case. Note that the
  library has to handle numeric-type conversion and for strings, padding
  and different `character` kinds.

__GCC__:
  Currently, it uses gfortran's private array descriptor. A change to [TS29113]'s
  array descriptor is planned; when that's done, the additional kind arguments
  will be removed.
  Note that the kind arguments permit to distinguish the `character` kinds and
  `real`/`complex` kinds 10 and 16, which have the same byte size.


### 3.9  Sending data between remote images ###

```c
void caf_sendget (caf_token_t dst_token, size_t dst_offset,
                  int dst_image_index, gfc_descriptor_t *dest,
                  caf_vector_t *dst_vector, caf_token_t src_token,
                  size_t src_offset, int src_image_index,
                  gfc_descriptor_t *src, caf_vector_t *src_vector,
                  int dst_kind, int src_kind)
```

Called to send a scalar, an array section or whole array from a remote image
identified by the `src_image_index` to a remote image identified by the
`dst_image_index`.

| Argument	| description	|
| ------	| ------	|
| `dst_token`	| `intent(in)`  An opaque pointer identifying the destination coarray.	|
| `dst_offset`	| By which amount of bytes the actual data is shifted compared to the base address of the destination coarray.	|
| `dst_image_index`	| The ID of the destination remote image; must be a positive number.	|
| `dest`	| `intent(in)` Array descriptor for the destination remote image for the bounds and the size. The `base_addr` shall not be accessed.	|
| `dst_vector`	| `intent(int)`  If not NULL, it contains the vector subscript of the destination array; the values are relative to the dimension triplet of the dest argument.	|
| `src_token`	| `intent(in)` An opaque pointer identifying the source coarray.	|
| `src_offset`	| By which amount of bytes the actual data is shifted compared to the base address of the source coarray.	|
| `src_image_index`	| The ID of the source remote image; must be a positive number.	|
| `src`	| `intent(in)` Array descriptor of the local array to be transferred to the remote image	|
| `src_vector`	| `intent(in)` Array descriptor of the local array to be transferred to the remote image	|
| `dst_kind`	| Kind of the destination argument	|
| `src_kind`	| Kind of the source argument	|

__Note__:
  It is permitted to have `image_id` equal the current image; the memory of the
  send-to and the send-from might (partially) overlap in that case. The
  implementation has to take care that it handles this case. Note that the
  assignment of a scalar to an array is permitted. In addition, the library has
  to handle numeric-type conversion and for strings, padding and different
  `character` kinds.

__GCC__:
  Currently, it uses gfortran's private array descriptor. A change to [TS29113]'s
  array descriptor is planned; when that's done, the additional kind arguments
  will be removed.
  Note that the kind arguments permit to distinguish the `character` kinds and
  `real`/`complex` kinds 10 and 16, which have the same byte size.



### 3.10  Barriers ###

### 3.10.1  All-Image Barrier ###

```c
void caf_sync_all (int *stat, char *errmsg, size_t errmsg_len)
```

Barrier which waits for all other images, pending asynchronous communication
and other data transfer.

| Argument	| description	|
| ------	| ------	|
| `stat`	| Status variable, if `NULL`, failures are fatal. If non-null, assigned 0 on success, and a stat code (cf. 2.3) in case of an error.	|
| `errmsg`	| If not NULL: Ignored unless stat is present; unmodified when successful, otherwise, an error message is copied into the variable.	|
| `errmsg_len`	| Maximal length of the error string, which is not '\0' terminated. The string should be padded by blanks.	|

__Note__:
  For portability, consider only using 7bit ASCII characters in the error
  message.

__GCC__:
  Implemented in GCC 4.x using an int argument for the length.
  Currently, `size_t` is not implemented.



### 3.10.2  Barrier for Selected Images ###

```c
void sync_images (int count, int images[], int *stat,
                  char *errmsg, size_t errmsg_len)
```

| Argument	| description	|
| ------	| ------	|
| `count`	| Size of the array "images"; has value -1 for `sync images(*)` and value 0 for a zero-sized array.	|
| `image`	| list of images to be synced with.	|
| `stat`	| Status variable, if NULL, failures are fatal. If non-null, assigned 0 on success, and a stat code (cf. 2.3) in case of an error.	|
| `errmsg`	| If not `NULL`: Ignored unless stat is present; unmodified when successful, otherwise, an error message is copied into the variable.	|
| `errmsg_len`	| Maximal length of the error string, which is not `\0` terminated. The string should be padded by blanks.	|

__Note__:
  For portability, consider only using 7bit ASCII characters in the error
  message. Note that the list can contain also the ID of `this_image` or can be
  an empty set. Example use is that image 1 syncs with all others (i.e `sync images(*)`) and the others sync only with that image (`sync image(1)`). Or
  for point-to point communication (`sync image([left_image, right_image]`).

__GCC__:
  Implemented in GCC 4.x using an int argument for the error-string length.
  Currently, `size_t` is not implemented.

### 3.11  Error abort ###

```c
void error_stop_str (const char *string, int32_t str_len);
void error_stop (int32_t exit_error_code)
```

__TODO__

  - [ ] Fix this description by filling-in the missing bits
  - [ ] `STOP` vs `ERROR STOP` handling. Currently, `STOP` calls `finalize` and then the
    normal `STOP` while for `ERROR STOP` directly calls the library
  - [ ] F2008 requires that one prints the raised exceptions with `STOP` and `ERROR
    STOP`. libgfortran's `STOP` and `ERROR STOP` do so - the current implementation
    for `ERROR STOP` does not.


### 3.11  Locking and unlocking ###

#### 3.11.1  Locking a lock variable ####

```c
void caf_lock (caf_token_t token, size_t index, int image_index,
               int *aquired_lock, int *stat, char *errmsg,
               int errmsg_len)
```

Acquire a lock on the given image on a scalar locking variable or for the
given array element for an array-valued variable. If the `acquired_lock`
is `NULL`, the function return after having obtained the lock. If it is
non-null, the result is is assigned the value true (one) when the lock could be
obtained and false (zero) otherwise.  Locking a lock variable which has already
been locked by the same image is an error.

| Argument	| arguments	|
| ------	| ------	|
| `token`	| `intent(in)` An opaque pointer identifying the coarray.	|
| `index`	| Array index; first array index is 0. For scalars, it is always 0.	|
| `image_index`	| The ID of the remote image; must be a positive number.	|
| `aquired_lock`	| `intent(out)` If not NULL, it returns whether lock could be obtained	|
| `stat`	| `intent(out)` For allocatable coarrays, stores the `STAT=`; may be NULL	|
| `errmsg`	| intent(out) When an error occurs, this will be set to an error message; may be NULL	|
| `errmsg_len`	| the buffer size of errmsg.	|

__Note__:
  This function is also called for critical sections; for those, the array index
  is always zero and the image index is one.  Libraries are permitted to use other
  images for critical-section locking variables.

__GCC__:
   (no special notes)

__TODO__:
   Change `errmsg_len` to `size_t`


#### 3.11.2  Unlocking a lock variable ####

```c
void caf_unlock (caf_token_t token, size_t index, int image_index,
                 int *stat, char *errmsg, int errmsg_len)
```

Release a lock on the given image on a scalar locking variable or for the
given array element for an array-valued variable. Unlocking a lock variable
which is unlocked or has been locked by a different image is an error.

| Argument	| description	|
| ------	| ------	|
| `token`	| `intent(in)` An opaque pointer identifying the coarray.	|
| `index`	| Array index; first array index is 0. For scalars, it is always 0.	|
| `image_index`	| The ID of the remote image; must be a positive number.	|
| `stat`	| `intent(out)` For allocatable coarrays, stores the `STAT=`; may be `NULL`	|
| `errmsg`	| `intent(out)` When an error occurs, this will be set to an error message; may be `NULL`	|
| `errmsg_len`	| the buffer size of `errmsg`.	|

__Note__:
  This function is also called for critical sections; for those, the array index
  is always zero and the image index is one.  Libraries are permitted to use other
  images for critical-section locking variables.

__GCC__:
   (no special notes)

__TODO__:
   Change `errmsg_len` to `size_t`

---

[![GitHub forks](https://img.shields.io/github/forks/sourceryinstitute/OpenCoarrays.svg?style=social&label=Fork)](https://github.com/sourceryinstitute/OpenCoarrays/fork)
[![GitHub stars](https://img.shields.io/github/stars/sourceryinstitute/OpenCoarrays.svg?style=social&label=Star)](https://github.com/sourceryinstitute/OpenCoarrays)
[![GitHub watchers](https://img.shields.io/github/watchers/sourceryinstitute/OpenCoarrays.svg?style=social&label=Watch)](https://github.com/sourceryinstitute/OpenCoarrays)
[![Twitter URL](https://img.shields.io/twitter/url/http/shields.io.svg?style=social)](https://twitter.com/intent/tweet?hashtags=HPC,Fortran,PGAS&related=zbeekman,gnutools,HPCwire,HPC_Guru,hpcprogrammer,SciNetHPC,DegenerateConic,jeffdotscience,travisci&text=Stop%20programming%20w%2F%20the%20%23MPI%20docs%20in%20your%20lap%2C%20try%20Coarray%20Fortran%20w%2F%20OpenCoarrays%20%26%20GFortran!&url=https%3A//github.com/sourceryinstitute/OpenCoarrays)


[Hyperlinks]:#

[TS29113]: ftp://ftp.nag.co.uk/sc22wg5/n1901-n1950/n1942.pdf
[TS18508]: http://isotc.iso.org/livelink/livelink?func=ll&objId=17288706&objAction=Open
[To Do]: #to-do
[Implementation status]: #implementation-status
[Definitions and types]: #definitions-and-types
[Provided functions]: #provided-functions
[pdf img]: https://img.shields.io/badge/PDF-CAF_ABI.md-6C2DC7.svg?style=flat-square "Download as PDF"
