/****p* doc/sample-compiler-output.c
 * NAME
 *   sendrcv
 * SYNOPSIS
 *   This C program represents the code an OpenCoarrays-compatible
 *   Fortran compiler might generate from doc/sample-fortran-source.f90.
 *   The code delegates all necessary synchronization and communicaiton
 *   to an OpenCoarrays transport layer.  In this program, image 1 puts
 *   its local elements of an array coarray into the corresponding
 *   elements of image 2.
 *
 * SOURCE
*/
sendrecv ()
{
  struct array2_real(kind=8) d;
  integer(kind=4) i;
  integer(kind=4) me;
  integer(kind=4) np;

  d.data = 0B;
  {
    integer(kind=4) overflow.0;

    overflow.0 = 0;
    if (overflow.0 != 0)
      {
        _gfortran_runtime_error (
         &"Integer overflow when calculating the amount of memory to allocate"[1]{lb: 1 sz: 1}
        );
      }
    else
      {
        if (d.data != 0B)
          {
            _gfortran_runtime_error_at (
              &"At line 9 of file sample-fortran-source.f90"[1]{lb: 1 sz: 1},
              &"Attempting to allocate already allocated variable \'%s\'"[1]{lb: 1 sz: 1}, &"d"[1]{lb: 1 sz: 1}
            );
          }
        else
          {
            d.data = (void * restrict) _gfortran_caf_register (8000, 1, &d.token, 0B, 0B, 0);
          }
      }
    d.dtype = 537;
    d.dim[0].lbound = 1;
    d.dim[0].ubound = 1000;
    d.dim[0].stride = 1;
    d.dim[1].lbound = 1;
    d.offset = -1;
    np = _gfortran_caf_num_images (0, -1);
    me = _gfortran_caf_this_image (0);
    i = 1;
    if (i <= 1000)
      {
        while (1)
          {
            {
              logical(kind=4) D.2368;

              (*(real(kind=8)[0:] * restrict) d.data)[d.offset + (integer(kind=8)) i] = (real(kind=8)) i;
              L.1:;
              D.2368 = i == 1000;
              i = i + 1;
              if (D.2368) goto L.2;
            }
          }
      }
    L.2:;
    __sync_synchronize ();
    _gfortran_caf_sync_all (0B, 0B, 0);
    if (me == 1)
      {
        _gfortran_caf_send (d.token, 0, 3 - (integer(kind=4)) d.dim[1].lbound, &d, 0B, &d, 8, 8);
      }
    L.3:;
    __sync_synchronize ();
    _gfortran_caf_sync_all (0B, 0B, 0);
    if (d.data == 0B)
      {
        _gfortran_runtime_error_at (
         &"At line 24 of file sample-fortran-source.f90"[1]{lb: 1 sz: 1},
         &"Attempt to DEALLOCATE unallocated \'%s\'"[1]{lb: 1 sz: 1}, &"d"[1]{lb: 1 sz: 1}
        );
      }
    else
      {
        _gfortran_caf_deregister (&d.token, 0B, 0B, 0);
      }
    d.data = 0B;
  }
}


main (integer(kind=4) argc, character(kind=1) * * argv)
{
  static integer(kind=4) options.1[9] = {68, 1023, 0, 0, 1, 1, 0, 0, 31};

  _gfortran_caf_init (&argc, &argv);
  _gfortran_set_args (argc, argv);
  _gfortran_set_options (9, &options.1[0]);
  sendrecv ();
  __sync_synchronize ();
  _gfortran_caf_finalize ();
  return 0;
}
/******
