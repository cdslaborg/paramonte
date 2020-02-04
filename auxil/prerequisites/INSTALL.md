[This document is formatted with GitHub-Flavored Markdown.               ]:#
[For better viewing, including hyperlinks, read it online at             ]:#
[https://github.com/sourceryinstitute/OpenCoarrays/blob/master/INSTALL.md]:#

<a name="top"> </a>

# Installing OpenCoarrays #

[![GitHub release][GH release badge]][latest release]
[![Github All Releases][GH all releases badge]][latest release]
[![Download as PDF][pdf img]][INSTALL.pdf]

Download this file as a PDF document
[here][INSTALL.pdf].

* [Developer Build and Install]
* [End-User Installation]
  * [macOS]
  * [Windows]
  * [Linux]
  * [FreeBSD]
  * [Virtual machine]
  * [Installation Script]
* [Advanced Installation from Source]
  * [Prerequisites]
  * [CMake]
  * [Make]

## Developer Build and Install ##

If you are a GCC developer, a package maintainer building OpenCoarrays
for distribution, or an advanced user who is comfortable building
software from source (using cmake), then we recommend installing
OpenCoarrays directly via CMake. If you do not fit into one of these
categories, we encourage you to skip ahead to review installation
options via Linux or MacOS package management software, or the
[`install.sh`] script. The text below is a condensed version of the
content available in [`INSTALL`]: plain text instructions for installing
OpenCoarrays in a canonical CMake way.

Prerequites for direct CMake installation:

* An MPI 3 implementation (MPICH is preferred, OpenMPI works too)
* A recent version of GCC with GFortran version 6.1 or newer
* CMake version 3.4 or newer

After obtaining the OpenCoarrays source (from git or our [latest release])
the following commands to build and install OpenCoarrays from source
using CMake:

```bash
mkdir opencoarrays-build
cd opencoarrays-build
export FC=/path/to/gfortran
export CC=/path/to/gcc
cmake /path/to/OpenCoarrays/source \
  -DCMAKE_INSTALL_PREFIX=/path/to/desired/installation/location
make
make test # optional; verify build works
make install
```

If you have either of the CMake gui tools installed, `ccmake` or
`cmake-gui` you may explore different configuration options and/or try
to locate/change which MPI version is found by repeating the steps
above and simply replacing `cmake` with `ccmake` or `cmake-gui`.

Please keep in mind that CMake cache variables are sticky and, once
set, can only be changed by using `ccmake`, `cmake-gui`, or explicitly
setting them on the command line: `cmake ../path/to/src -DVAR=VALUE`
If the wrong compiler or MPI implementation is being used and you
cannot determine why, you can try deleting the entire build directly
and re-running CMake.

## End-User Installation ##

Most users will find it easiest and fastest to use package management
software to install OpenCoarrays. Below is the status of OpenCoarrays
in various package managers. If you do not see your favorite package
manager listed or it is badly out of date, please reach out and ask
for it to be included or updated, or contribute it yourself. We are
happy to work with package managers to resolve issues, and adapt our
build system to play nicely so they don't have to maintain patches.

[![Packaging status][repology-badge]][OC-on-repology]

Package management options for
macOS, Windows, and Linux are described first below. Also described
below are options for installing via the Sourcery Institute virtual
machine or via the bash scripts included that are in the OpenCoarrays
source.

[top]

### macOS ###

* [Homebrew]:
  [![homebrew][Homebrew badge]][braumeister link]
  This is the recommended OpenCoarrays installation method on macOS.
  Basic Homebrew installation steps:

  ```
  brew update
  brew install opencoarrays
  ```

  OpenCoarrays also ships with a [`Brewfile`][Brewfile]
  that will make it easier to install opencoarrays using MPICH built
  with the GNU Compiler Collection ([GCC]). To install using the
  [`Brewfile`][Brewfile] with MPICH wrapping GCC, follow these steps:

  ```
  brew tap homebrew/bundle
  brew update
  brew bundle
  ```

* [MacPorts]:
  An unmaintained [OpenCoarrays Portfile] exists for the [MacPorts] package
  manager.  Although the current OpenCoarrays contributors have no plans to
  update the portfile, new contributors are welcome to asssume the port
  maintainer role and to submit a pull request to update this [INSTALL.md] file.

[top]

### Windows ###

Windows users may run the [`install.sh`] script inside the Windows Subsystem
for Linux ([WSL]). The script uses Ubuntu's [APT] package manager to build
[GCC] 5.4.0, [CMake], and [MPICH].  Windows users who desire a newer version
of GCC are welcome to submit a request via our [Issues] page and suggest a
method for updating. Previously attempted upgrade methods are described in
the discussion thread starting with [commit comment 20539810].

[top]

### Linux ###

Access OpenCoarrays on Linux via any of the following package managers
or pre-installed copies:

* The [linuxbrew] package manager installs OpenCoarrays on all Linux distributions.
* Debian-based distributions such as Ubuntu provide an "open-coarrays" [APT package].
* [Arch Linux] provides an [aur package].
* An [HPCLinux] installation script is in the [developer-scripts] subdirectory (available
  via git clone only).
* [EasyBuild] can install OpenCoarrays on Linux distributions
* [Spack], a multiplatform package manager, can also install OpenCoarrays on Linux
   distributions

[linuxbrew] does not require `sudo` privileges and will generally
provide the most up-to-date OpenCoarrays release because linuxbrew
pulls directly from macOS homebrew, which updates automatically.

Note that distributions are often split into two parts: a "binary" package with the
runtime libraries and a "development" package for developing programs yourself.
Be sure to install both packages.

<a name="easybuild"></a>
With [EasyBuild], the following bash commands install OpenCoarrays:

```bash
# Search available specification files (also known as easyconfigs) for OpenCoarrays
eb --search OpenCoarrays

# Automatically download prerequisites (with the --robot flag) and install OpenCoarrays
# with the desired easyconfig, e.g., OpenCoarrays-1.9.0-gompi-2017a.eb
eb OpenCoarrays-1.9.0-gompi-2017a.eb --robot
```

Once installed, use OpenCoarrays by loading the newly created environment
module `OpenCoarrays/1.9.0-gompi-2017a`:

```
module load OpenCoarrays/1.9.0-gompi-2017a
```

<a name="spack"></a>
With [Spack], the following commands install OpenCoarrays in a bash shell:

```bash
# Check build information for OpenCoarrays in the default specification file
spack spec opencoarrays

# To automatically download prerequisites and install OpenCoarrays with the default
# specification.
# (Note: In addition to its own prerequisites, Spack requires gfortran compiler
# to be installed to compile OpenMPI)
spack install opencoarrays

# Or, To install with customisations (e.g., to install OpenCoarrays [version 2.2.0]
# with MPICH [version default] and GCC [version 7.1.0]).
spack install opencoarrays@2.2.0 ^mpich %gcc@7.1.0
```

In the previous example, it was assumed that GCC [version 7.1.0] is
already installed, and is available as a compiler to Spack. Otherwise,
[add a new compiler to Spack]. Once installed, OpenCoarrays can be
used by [loading the environment modules with Spack], e.g.

```
spack module loads --dependencies opencoarrays
```

[top]

### FreeBSD ###

Use the OpenCoarrays FreeBSD, Port to install OpenCoarrays by
executing the following commands as root:

```
pkg install opencoarrays
```

For more information, please review the [FreeBSD ports/packages installation information].

[top]

## Virtual machine ##

Users of macOS, Windows, or Linux have the option to use OpenCoarrays
by installing the Lubuntu Linux virtual machine from the
[Sourcery Institute Store].  The virtual machine boots inside the
open-source [VirtualBox] virtualization package.  In addition to
containing [GCC], [MPICH], and OpenCoarrays, the virtual machine
contains dozens of other open-source software packages that support
modern Fortran software development.  See the
[download and installation instructions] for a partial list of the
included packages.

[top]

## Installation Script ##

If the above package management or virtualization options are
infeasible or unavailable, Linux, macOS, and [WSL] users may also install
OpenCoarrays by downloading and uncompressing our [latest release] and
running our installation script in the top-level OpenCoarrays source
directory (see above for the corresponding [Windows] script):

```
tar xvzf OpenCoarrays-x.y.z.tar.gz
cd OpenCoarrays-x.y.z
./install.sh
```

where `x.y.z` should be replaced with the appropriate version numbers.
A complete installation should result in the creation of the following
directories inside the installation path (.e.g, inside `build` in the
above example):

* `bin`: contains the compiler wrapper (`caf`) and program launcher
  (`cafun`).
* `mod`: contains the `opencoarrays.mod` module file for use with
  non-OpenCoarrays-aware compilers
* `lib`: contains the `libcaf_mpi.a` static library to which codes
  link for CAF support

### Example script invocations ###

Execute `./install.sh --help` or `./install.sh -h` to see a list of flags
that can be passed to the installer.  Below are examples of useful combinations
of flags. Each flag also has a single-character version not shown here.

1. If you don't care about tailoring installation locations,use

   ```
   ./install.sh
   ```

2. Or build faster by multithreading, skipping user queries, and also specify
    an installation destination for all packages that must be built:

   ```
   ./install.sh \
     --prefix-root ${HOME}/software \
     --num-threads 4 \
     --yes-to-all
   ```

3. Specify the compilers to be used (overriding what is in your `PATH`) as follows:

   ```
   ./install.sh --with-fortran <path-to-gcc-bin>/gfortran \
                --with-cxx <path-to-gcc-bin>/g++ \
                --with-c <path-to-gcc-bin>/gcc
   ```

4. Install only a specific prerequisite package (the default version):

   ```
   ./install.sh --package mpich
   ```

5. Install a specific version of a prerequisite:

   ```
   ./install.sh --package cmake --install-version 3.7.0
   ```

6. Download a prerequisite package (e.g., gcc/gfortran/g++ below) but
   don't build or install it:

   ```
   ./install.sh --only-download gcc
   ```

7. Print the default URL, version, or download mechanism that the
   script will use for a given prerequisite package (e.g., mpich
   below) on this system:

   ```
   ./install.sh --print-url mpich
   ./install.sh --print-version mpich
   ./install.sh --print-downloader mpich
   ```

8. Install a prerelease branch (e.g., trunk below) of the GCC repository:

   ```
   ./install.sh --package gcc --install-branch trunk
   ```

[top]

## Advanced Installation from Source ##

### Prerequisites ###

Package managers and the [`install.sh`] attempt to handle the installation
of all OpenCoarrays prerequisites automatically.  Installing with CMake
or the provided, static Makefile burdens the person installing with the
need to ensure that all prerequisites have been built and are in the
expected or specified locations prior to building OpenCoarrays. The
prerquisite package/version dependency tree is as follows:

```text
opencoarrays
├── cmake-3.4.0
└── mpich-3.2
    └── gcc-6.1.0
        ├── flex-2.6.0
        │   └── bison-3.0.4
        │       └── m4-1.4.17
        ├── gmp
        ├── mpc
        └── mpfr
```

[top]

### CMake ###

On most platforms, the [`install.sh`] script ultimately invokes [CMake] after
performing numerous checks, customizations, and installations of any missing
prerequisites. Users wishing to install OpenCoarrays directly with CMake should
have a look at the documentation in the [`./INSTALL`] file if they encounter
issues or need further guidance. A brief summary is also given at the top of
this document [here][Developer Build and Install].

[top]

### Make ###

Unlike the Makefiles that CMake generates automatically for the chosen
platform, static Makefiles require a great deal more maintenance and are
less portable.  Also, the static Makefiles provided in [src] lack several
important capabilities.  In particular, they will not build the tests;
they will not generate the `caf` compiler wrapper that ensures correct linking
and `cafrun` program launcher that ensures support for advanced features such
as Fortran 2015 failed images; they will not build the [opencoarrays] module
that can be used to provide some Fortran 2015 features with non-Fortran-2015
compilers; nor do the static Makefiles provide a `make install` option so you
will need to manually move the resultant library from the build location to your
chosen installation location.

If none of the installation methods mentioned higher in this document are
work on your platform and if CMake is unavailable, build and install the
OpenCoarrays parallel runtime library as follows:

```
tar xvzf opencoarrays.tar.gz
cd opencoarray/src
make
mv mpi/libcaf_mpi.a <insert-install-path>
```

replacing the angular-bracketed text with your desired install path.

For the above steps to succeed, you might need to edit the [make.inc]
file to match your system settings.  For example, you might need to
remove the `-Werror` option from the compiler flags or name a
different compiler.  In order to activate efficient strided-array
transfer support, uncomment the `-DSTRIDED` flag inside the [make.inc]
file.

[top]

---

<div align="center">

[![GitHub forks][fork badge]][fork url]
[![GitHub stars][star badge]][star url]
[![GitHub watchers][watch badge]][star url]
[![Twitter URL][twitter badge]][twitter url]

</div>

[Internal document links]: #

[top]: #top
[Developer Build and Install]: #developer-build-and-install
[End-User Installation]: #end-user-installation
[macOS]: #macos
[Windows]: #windows
[Linux]: #linux
[FreeBSD]: #freebsd
[Virtual machine]: #virtual-machine
[Installation Script]: #installation-script

[Advanced Installation from Source]: #advanced-installation-from-source
[Prerequisites]: #prerequisites
[CMake]: #cmake
[Make]: #make

[Links to source]: #

[`install.sh`]: ./install.sh
[`INSTALL`]: ./INSTALL

[URLs]: #

[OC-on-repology]: https://repology.org/project/opencoarrays/versions
[repology-badge]: https://repology.org/badge/vertical-allrepos/opencoarrays.svg

[FreeBSD ports/packages installation information]: https://www.freebsd.org/doc/en_US.ISO8859-1/books/handbook/ports.html
[GH all releases badge]: https://img.shields.io/github/downloads/sourceryinstitute/OpenCoarrays/total.svg?style=flat-square
[GH release badge]: https://img.shields.io/github/release/sourceryinstitute/OpenCoarrays.svg?style=flat-square
[Homebrew badge]: https://img.shields.io/homebrew/v/opencoarrays.svg?style=flat-square
[linuxbrew]: http://linuxbrew.sh
[braumeister link]: http://braumeister.org/formula/opencoarrays
[APT package]: https://qa.debian.org/popcon.php?package=open-coarrays
[APT]: https://en.wikipedia.org/wiki/APT_(Debian)
[HPCLinux]: http://www.paratools.com/hpclinux/
[Brewfile]: https://github.com/sourceryinstitute/OpenCoarrays/blob/master/Brewfile
[INSTALL.pdf]: https://md2pdf.herokuapp.com/sourceryinstitute/OpenCoarrays/blob/master/INSTALL.pdf
[INSTALL.md]: https://github.com/sourceryinstitute/OpenCoarrays/blob/master/INSTALL.md
[CMake]: https://cmake.org
[Sourcery Institute Store]: http://www.sourceryinstitute.org/store/c1/Featured_Products.html
[VirtualBox]: https://www.virtualbox.org
[download and installation instructions]: http://www.sourceryinstitute.org/uploads/4/9/9/6/49967347/overview.pdf
[yum]: http://yum.baseurl.org
[apt-get]: https://en.wikipedia.org/wiki/Advanced_Packaging_Tool
[Issues]: https://github.com/sourceryinstitute/OpenCoarrays/issues
[make.inc]: ./src/make.inc
[opencoarrays]: ./src/extensions/opencoarrays.F90
[prerequisites]: #prerequisites
[MPICH]: http://www.mpich.org
[MVAPICH]:http://mvapich.cse.ohio-state.edu
[MacPorts]: https://www.macports.org
[GCC]: http://gcc.gnu.org
[TS18508 Additional Parallel Features in Fortran]: https://isotc.iso.org/livelink/livelink/nfetch/-8919044/8919782/8919787/17001078/ISO%2DIECJTC1%2DSC22%2DWG5_N2056_Draft_TS_18508_Additional_Paralle.pdf?nodeid=17181227&vernum=0
[GFortran Binaries]:  https://gcc.gnu.org/wiki/GFortranBinaries#FromSource
[Installing GCC]: https://gcc.gnu.org/install/
[Arch Linux]: https://www.archlinux.org
[aur package]: https://aur.archlinux.org/packages/opencoarrays/
[latest release]: https://github.com/sourceryinstitute/OpenCoarrays/releases/latest
[pdf img]: https://img.shields.io/badge/PDF-INSTALL.md-6C2DC7.svg?style=flat-square
[commit comment 20539810]: https://github.com/sourceryinstitute/OpenCoarrays/commit/26e99919fe732576f7277a0e1b83f43cc7c9d749#commitcomment-20539810
[Homebrew]: https://brew.sh
[dnf]: https://github.com/rpm-software-management/dnf
[port details]: https://www.freshports.org/lang/opencoarrays
[port search]: https://www.freebsd.org/cgi/ports.cgi?query=opencoarrays
[EasyBuild]: https://github.com/easybuilders/easybuild
[Spack]: https://github.com/spack/spack
[add a new compiler to Spack]: http://spack.readthedocs.io/en/latest/tutorial_modules.html#add-a-new-compiler
[loading the environment modules with Spack]: http://spack.readthedocs.io/en/latest/module_file_support.html#cmd-spack-module-loads
[OpenCoarrays Portfile]: https://www.macports.org/ports.php?by=name&substr=opencoarrays
[WSL]: https://blogs.msdn.microsoft.com/commandline/2017/07/10/ubuntu-now-available-from-the-windows-store/
[developer-scripts]: https://github.com/sourceryinstitute/OpenCoarrays/tree/master/developer-scripts
[src]: https://github.com/sourceryinstitute/OpenCoarrays/tree/master/src
[fork badge]: https://img.shields.io/github/forks/sourceryinstitute/OpenCoarrays.svg?style=social&label=Fork
[fork url]: https://github.com/sourceryinstitute/OpenCoarrays/fork
[star badge]: https://img.shields.io/github/stars/sourceryinstitute/OpenCoarrays.svg?style=social&label=Star
[star url]: https://github.com/sourceryinstitute/OpenCoarrays
[watch badge]: https://img.shields.io/github/watchers/sourceryinstitute/OpenCoarrays.svg?style=social&label=Watch
[twitter badge]: https://img.shields.io/twitter/url/http/shields.io.svg?style=social
[twitter url]: https://twitter.com/intent/tweet?hashtags=HPC,Fortran,PGAS&related=zbeekman,gnutools,HPCwire,HPC_Guru,hpcprogrammer,SciNetHPC,DegenerateConic,jeffdotscience,travisci&text=Stop%20programming%20w%2F%20the%20%23MPI%20docs%20in%20your%20lap%2C%20try%20Coarray%20Fortran%20w%2F%20OpenCoarrays%20%26%20GFortran!&url=https%3A//github.com/sourceryinstitute/OpenCoarrays
