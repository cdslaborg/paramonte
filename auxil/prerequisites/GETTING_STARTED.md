<a name="top"> </a>

[This document is formatted with GitHub-Flavored Markdown.                       ]:#
[For better viewing, including hyperlinks, read it online at                     ]:#
[https://github.com/sourceryinstitute/OpenCoarrays/blob/master/GETTING_STARTED.md]:#

Getting Started
===============

[![Download as PDF][pdf img]](http://md2pdf.herokuapp.com/sourceryinstitute/OpenCoarrays/blob/master/GETTING_STARTED.pdf)

Download this file as a PDF document
[here](http://md2pdf.herokuapp.com/sourceryinstitute/OpenCoarrays/blob/master/GETTING_STARTED.pdf).

* [The caf compiler wrapper]
* [A sample basic workflow]
* [An advanced workflow]

The caf compiler wrapper
--------------------------

The preferred method for compiling a CAF program is by invoking the `caf` bash script
that the OpenCoarrays CMake scripts install in the `bin` subdirectory of the installation
path. This is an experimental script with limited but useful capabilities that will
grow over time.  Please submit bug reports and feature requests via our [Issues] page.

The `caf` script liberates the source code and workflow from explicit dependence on the
underlying compiler and communication library by passing the unmodified source code to 
the compiler with the necessary arguments for building a parallel Fortran 2018 program,
embedding the paths to OpenCoarrays libraries (e.g., `libcaf_mpi.a`) installed
in the `lib` subdirectory of the OpenCoarrays installation path.  

A sample basic workflow
-----------------------

The following program listing, compilation, and execution workflow exemplify
the use of OpenCoarrays with GCC in a Linux bash shell with the `bin`
directory of the chosen installation path in the user's `PATH` environment variable:

```fortran
$ cat tally.f90
      program main
        use iso_c_binding, only : c_int
        use iso_fortran_env, only : error_unit
        implicit none
        integer(c_int) :: tally
        tally = this_image() ! this image's contribution
        call co_sum(tally)
        verify: block
          integer(c_int) :: image
          if (tally/=sum([(image,image=1,num_images())])) then
             write(error_unit,'(a,i5)') "Incorrect tally on image ",this_image()
             error stop
          end if
        end block verify
        ! Wait for all images to pass the test
        sync all
        if (this_image()==1) print *,"Test passed"
      end program
$ caf tally.f90 -o tally
$ cafrun -np 4 ./tally
        Test passed
```

where "4" is the number of images to be launched at program start-up.

An advanced workflow
--------------------

If you prefer to invoke the compiler directly, first run `caf` and `cafrun` with the `--show` flag 
to see the proper linking and file includes.  For example, on a macOS system where OpenCoarrays 
was installed via the [homebrew] package manager, the following results:

```bash
$ caf --show
/usr/local/bin/gfortran -I/usr/local/Cellar/opencoarrays/2.2.0/include/OpenCoarrays-2.2.0_GNU-8.2.0 -fcoarray=lib -Wl,-flat_namespace -Wl,-commons,use_dylibs -L/usr/local/Cellar/libevent/2.1.8/lib -L/usr/local/Cellar/open-mpi/3.1.1/lib ${@} /usr/local/Cellar/opencoarrays/2.2.0/lib/libcaf_mpi.a /usr/local/lib/libmpi_usempif08.dylib /usr/local/lib/libmpi_usempi_ignore_tkr.dylib /usr/local/lib/libmpi_mpifh.dylib /usr/local/lib/libmpi.dylib
$ cafrun --show
/usr/local/bin/mpiexec -n <number_of_images> /path/to/coarray_Fortran_program [arg4 [arg5 [...]]]
```

---

[![GitHub forks](https://img.shields.io/github/forks/sourceryinstitute/OpenCoarrays.svg?style=social&label=Fork)](https://github.com/sourceryinstitute/OpenCoarrays/fork)
[![GitHub stars](https://img.shields.io/github/stars/sourceryinstitute/OpenCoarrays.svg?style=social&label=Star)](https://github.com/sourceryinstitute/OpenCoarrays)
[![GitHub watchers](https://img.shields.io/github/watchers/sourceryinstitute/OpenCoarrays.svg?style=social&label=Watch)](https://github.com/sourceryinstitute/OpenCoarrays)
[![Twitter URL](https://img.shields.io/twitter/url/http/shields.io.svg?style=social)](https://twitter.com/intent/tweet?hashtags=HPC,Fortran,PGAS&related=zbeekman,gnutools,HPCwire,HPC_Guru,hpcprogrammer,SciNetHPC,DegenerateConic,jeffdotscience,travisci&text=Stop%20programming%20w%2F%20the%20%23MPI%20docs%20in%20your%20lap%2C%20try%20Coarray%20Fortran%20w%2F%20OpenCoarrays%20%26%20GFortran!&url=https%3A//github.com/sourceryinstitute/OpenCoarrays)

[Hyperlinks]:#

[The caf compiler wrapper]: #the-caf-compiler-wrapper
[A sample basic workflow]: #a-sample-basic-workflow
[An advanced workflow]:  #an-advanced-workflow

[Sourcery Store]: http://www.sourceryinstitute.org/store
[Issues]: https://github.com/sourceryinstitute/OpenCoarrays/issues
[opencoarrays module]: ./src/extensions/opencoarrays.F90
[GCC]: http://gcc.gnu.org
[TS 18508]: http://isotc.iso.org/livelink/livelink?func=ll&objId=17181227&objAction=Open
[The caf compiler wrapper]: #the-caf-compiler-wrapper
[The cafrun program launcher]: #the-cafrun-program-launcher
[pdf img]: https://img.shields.io/badge/PDF-GETTING_STARTED.md-6C2DC7.svg?style=flat-square "Download as PDF"
[homebrew]: https://brew.sh
