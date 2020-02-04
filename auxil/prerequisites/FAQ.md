Frequently Asked Questions
==========================

* [1. I have installed OpenCoarrays with OpenMPI but I'm having trouble running with more than a few images, why?]
* [2. I have installed OpenCoarrays with a recent version of OpenMPI, but my Coarray Fortran programs won't run when launched with cafrun, why?]
* [3. `install.sh` is trying to download and install GCC/GFortran and its prerequisites, but I want to use GCC version X already present on my system, how can I do this?]
* [4. How can I uninstall OpenCoarrays?]
* [5. OpenCoarrays was built with MPICH as the MPI back end, but I am running into bugs, what should I do?]
* [6. How can I pass additional flags through to the underlying parallel run-time or compiler?]
* [7. When `install.sh` builds the GCC compilers, it takes forever (hours).  How can I speed up the build?]

## 1. I have installed OpenCoarrays with OpenMPI but I'm having trouble running with more than a few images, why? ##

OpenMPI requires oversubscribed jobs (more MPI ranks/coarray
images than logical CPU cores available) to pass the
`--oversubscribe` flag to `mpirun`. To run a Coarray program with
more images than available logical CPU cores, please pass
`--oversubscribe` as an argument to the `cafrun` wrapper
script. `cafrun` will pass any additional flags through to the
underlying run-time launcher (i.e., `mpirun`). Here is an example
invocation of `cafrun` demonstrating this procedure:

```
cafrun -np 32 --oversubscribe ./a.out <arg1> <arg2>
```

## 2. I have installed OpenCoarrays with a recent version of OpenMPI, but my Coarray Fortran programs won't run when launched with `cafrun`, why? ##

Recent versions of OpenMPI require a hostfile to be used,
specifying the number of "slots" (i.e., logical CPU cores)
available on each node or host on the system. If you are running
OpenCoarrays on a production HPC system, your system
administrator should have configured OpenMPI with a default hosts
file or should provide instructions on how to generate your
own. If you are running OpenCoarrays on a shared memory machine,
such as a laptop, desktop, or workstation, then you may have to
provide your own hostfile. On a shared memory machine, this
hostfile should look something like:

```
hostname.local 8
```

where `hostname` is the hostname of the machine you are running on
and `8` is the number of logical CPU cores. On a Macbook Pro (late
2013) with a 4 core Intel core i7 CPU with two hyper threads per
core, 8 is the number of logical CPU cores. When OpenCoarrays is
configured and installed, if OpenMPI is detected a hostfile is
created in the build directory for use when running the test
suite. (This hostfile corresponds to the machine/node used to
configure OpenCoarrays using CMake.) If you installed OpenCoarrays
via the [`install.sh`] script, then it will be located at
`./prerequisites/builds/opencoarrays/<version>/hostfile` or if you
invoked CMake yourself it should be in the top level of your build
directory at `/path/to/build/dir/hostfile`. (You must run CMake
first, to create this file.)

## 3. `install.sh` is trying to download and install GCC/GFortran and its prerequisites, but I want to use GCC version X already present on my system, how can I do this? ##

[`install.sh`] defaults to trying to install the most recent stable
version of GFortran and GCC that has the best functionality with
OpenCoarrays. This is usually the latest stable version of
GCC/GFortran, however, sometimes regressions are present and a
slightly older version is set as the default until most of the serious
regressions are resolved. If you know what you are doing and wish to
use a different Fortran/C compiler, you can explicitly pass arguments
to [`install.sh`] to specify which compilers to use. Please see the
usage information for [`install.sh`] by invoking it with the `--help`
flag. For example, the `-f` or `--with-fortran` flag is passed with
the path to the desired Fortran compiler, the `-c` or `--with-c` flag
is passed with the path to the desired C compiler and the `-C` or
`--with-cxx` flag is passed with the path to the desired C++
compiler. Alternatively the "Developer/quick-start instructions"
(located in the [`INSTALL`] file) may be followed to perform the
installation using CMake directly, without calling [`install.sh`]. The
compilers are specified using the `FC`, `CC` and `CXX` environment
variables, and all prerequisites (GFortran, CMake, and a suitable MPI
implementation) are assumed to be already installed on the system. If
you pass one compiler-specification argument (e.g., `-f`), it is best
to pass all three (e.g., `-f`, `-c`, and `-C`) to ensure consistency.

## 4. How can I uninstall OpenCoarrays? ##

After installing OpenCoarrays, you can enter the build directory
(`prerequisites/builds/opencoarrays/<version>` if installed via
[`install.sh`]) and run `make uninstall`. A script or additional flag
to [`install.sh`] is planned to automate this process, but has yet to
be implemented. When OpenCoarrays is installed a file manifest is
written in this build directory, so you can examine which files get
installed and their locations on your system.

## 5. OpenCoarrays was built with MPICH as the MPI back end, but I am running into bugs, what should I do? ##

There is a bug in MPICH that, as of this writing, is patched, but is
not yet in a stable release. This bug effects the failed images
functionality, which can be disabled during configuration. If you
encounter issues using OpenCoarrays with MPICH, reinstalling
OpenCoarrays _without_ failed image support may resolve the issue. To
do so you will have to perform a manual CMake configuration and
installation as described in the "Developer/quick start" installation
guide documented in the file [`INSTALL`]. Please pass the
`-DCAF_ENABLE_FAILED_IMAGES=FALSE` flag to CMake when configuring
OpenCoarrays. If this doesn't resolve your issue, please file a
[new issue].

## 6. How can I pass additional flags through to the underlying parallel run-time or compiler? ##

The `caf` compiler wrapper script and the `cafrun` Coarray Fortran
program launch script will both pass additional flags through to
the underlying compiler and parallel run-time job launcher,
respectively. Specify the additional flags to be passed to
GFortran or `mpirun` after any flags specific to `caf` or `cafrun`
such as the `-s` flag to show the underlying command, or the `-np <N>`
flag to specify the number of images in the `cafrun` script and before
any files such as Fortran source files or Coarray Fortran executables.

## 7. When `install.sh` builds the GCC compilers, it takes forever (hours).  How can I speed up the build? ##

To increase the odds of success, `install.sh` defaults to a GCC
bootstrap build, which builds a minimal compiler to build the ultimate
compiler (not every version of GCC can build every other version of
GCC).  For a much faster build process that has a somewhat higher
chance of failing, pass the `--disable-bootstrap` or `-z` argument and
use more threads by passing, for example, `--num-threads -4` or `-j 4`
to use four threads.  In combination, these two recommendations can
decrease the GCC build time from several hours to 15 or fewer minutes.


[`install.sh`]: https://github.com/sourceryinstitute/OpenCoarrays/blob/master/install.sh
[`INSTALL']: https://github.com/sourceryinstitute/OpenCoarrays/blob/master/INSTALL
[new issue]: https://github.com/sourceryinstitute/OpenCoarrays/issues/new

[TOC links]: #
[1. I have installed OpenCoarrays with OpenMPI but I'm having trouble running with more than a few images, why?]: #1-i-have-installed-opencoarrays-with-openmpi-but-im-having-trouble-running-with-more-than-a-few-images-why
[2. I have installed OpenCoarrays with a recent version of OpenMPI, but my Coarray Fortran programs won't run when launched with cafrun, why?]: #2-i-have-installed-opencoarrays-with-a-recent-version-of-openmpi-but-my-coarray-fortran-programs-wont-run-when-launched-with-cafrun-why
[3. `install.sh` is trying to download and install GCC/GFortran and its prerequisites, but I want to use GCC version X already present on my system, how can I do this?]: #3-installsh-is-trying-to-download-and-install-gccgfortran-and-its-prerequisites-but-i-want-to-use-gcc-version-x-already-present-on-my-system-how-can-i-do-this
[4. How can I uninstall OpenCoarrays?]: #4-how-can-i-uninstall-opencoarrays
[5. OpenCoarrays was built with MPICH as the MPI back end, but I am running into bugs, what should I do?]: #5-opencoarrays-was-built-with-mpich-as-the-mpi-back-end-but-i-am-running-into-bugs-what-should-i-do
[6. How can I pass additional flags through to the underlying parallel run-time or compiler?]: #6-how-can-i-pass-additional-flags-through-to-the-underlying-parallel-run-time-or-compiler
[7. When `install.sh` builds the GCC compilers, it takes forever (hours).  How can I speed up the build?]: #7-when-installsh-builds-the-gcc-compilers-it-takes-forever-hours--how-can-i-speed-up-the-build
