<a name="top"> </a>

[This document is formatted with GitHub-Flavored Markdown.              ]:#
[For better viewing, including hyperlinks, read it online at            ]:#
[https://github.com/sourceryinstitute/OpenCoarrays/blob/master/README.md]:#
<div align="center">

[![Sourcery Institute][sourcery-institute logo]][Sourcery, Inc.]

OpenCoarrays
============

[![CI Build Status][build img]](https://travis-ci.org/sourceryinstitute/OpenCoarrays)
[![Release Downloads][download img]][Releases]
[![Gitter](https://img.shields.io/gitter/room/sourceryinstitute/opencoarrays.svg?style=flat-square)](https://gitter.im/sourceryinstitute/opencoarrays)
[![GitHub license][license img]](./LICENSE)
[![GitHub release][release img]](https://github.com/sourceryinstitute/OpenCoarrays/releases/latest)
[![homebrew](https://img.shields.io/homebrew/v/opencoarrays.svg?style=flat-square)](https://formulae.brew.sh/formula/opencoarrays)
[![Download as PDF][pdf img]](https://md2pdf.herokuapp.com/sourceryinstitute/OpenCoarrays/blob/master/README.pdf)
[![Twitter URL][twitter img]][default tweet]

[Overview](#overview) | [Downloads](#downloads) |
[Compatibility](#compatibility) | [Prerequisites](#prerequisites) |
[Installation](#installation) | [Getting Started](#getting-started) |
[Contributing](#contributing) | [Status](#status)
[Support](#support) | [Acknowledgments](#acknowledgments) | [Donate](#donate)

</div>

Overview
--------

[OpenCoarrays] supports [Fortran 2018] compilers by providing a
parallel application binary interface (ABI) that abstracts away the
underlying parallel programming model, which can be the Message
Passing Interface ([MPI]) or [OpenSHMEM].  Parallel Fortran 2018
programs may be written and compiled into object files once and
then linked or re-linked to either MPI or [OpenSHMEM] without modifying
or recompiling the Fortran source.  Not a single line of source code
need change to switch parallel programming models.  The default
programming model is MPI because it provides the broadest capability
for supporting Fortran 2018  features.  However, having the option to
change parallel programming models at link-time may enhance portability
and performance (see [Rouson et al. (2017)]).

OpenCoarrays provides a compiler wrapper (`caf`), parallel runtime
libraries (`libcaf_mpi` and `libcaf_openshmem`), and a parallel
executable file launcher (`cafrun`).  The wrapper and launcher
provide a uniform abstraction for compiling and executing parallel
Fortran 2018 programs without direct reference to the underlying
parallel programming model.

Downloads
---------

Please see our [Releases] page.

Compatibility
-------------

The GNU Compiler Collection ([GCC]) Fortran front end ([gfortran]) has
used OpenCoarrays since the GCC 5.1.0 release .  Discussions are under
way around incorporating OpenCoarrays into other compilers.

Prerequisites
-------------

Building OpenCoarrays requires

* An MPI implementation (default: [MPICH]).
* CMake.
* A Fortran compiler (default: [GCC]).
* _Optional_: An [OpenSHMEM] implementation.

If you use a package manager or the OpenCoarrays installer, any
missing prerequisites will be built for you.


Installation
------------

Please see the [INSTALL.md] file.

Or [try OpenCoarrays online] as a [Jupyter] [notebook kernel]
using [Binder] with no downloads, configuration or installation required.
The default [index.ipynb] notebook is read only, but you can
execute it, copy it to make changes, or create an entirely
new [CAF kernel][notebook kernel] notebook.

Packaged Version
----------------

If you would like to be able to install OpenCoarrays through your
favorite package manager, please ask them to add it, or contribute it
yourself. If you see your favorite package manager has an outdated
version, please ask them to update it, or contribute an update
yourself.

[![Packaging status][repology-badge]][OC-on-repology]

Getting Started
---------------

To start using OpenCoarrays, please see the [GETTING_STARTED.md] file.

Contributing
------------

Please see the [CONTRIBUTING.md] file.

Status
------

A list of open issues can be viewed on the
[issues page](https://github.com/sourceryinstitute/opencoarrays/issues).

Support
-------

Please submit bug reports and feature requests via our [Issues] page.

Acknowledgments
----------------

We gratefully acknowledge support from the following institutions:

* [Arm] for approving compiler engineer contributions of code.
* [National Center for Atmospheric Research] for access to the
  Yellowstone/Caldera supercomputers and for logistics support during
  the initial development of OpenCoarrays.
* [CINECA] for access to Eurora/PLX for the project HyPS- BLAS under
  the ISCRA grant program for 2014.
* [Google] for support of a related [Google Summer of Code] 2014
  project.
* The National Energy Research Scientific Computing Center ([NERSC]),
  which is supported by the Office of Science of the U.S. Department
  of Energy under Contract No. DE-AC02-05CH11231, for access to the
  Hopper and Edison supercomputers under the OpenCoarrays project
  start allocation.
* [Sourcery, Inc.], for financial support for the domain registration,
  web hosting, advanced development, and conference travel.

Donate
------

If you find this software useful, please consider donating
[your time](CONTRIBUTING.md) or
[your money](http://www.sourceryinstitute.org/store/p5/Donation.html)
to aid in development efforts.

---

<div align="center">

[![GitHub forks](https://img.shields.io/github/forks/sourceryinstitute/OpenCoarrays.svg?style=social&label=Fork)](https://github.com/sourceryinstitute/OpenCoarrays/fork)
[![GitHub stars](https://img.shields.io/github/stars/sourceryinstitute/OpenCoarrays.svg?style=social&label=Star)](https://github.com/sourceryinstitute/OpenCoarrays)
[![GitHub watchers](https://img.shields.io/github/watchers/sourceryinstitute/OpenCoarrays.svg?style=social&label=Watch)](https://github.com/sourceryinstitute/OpenCoarrays)
[![Twitter URL][twitter img]][default tweet]

</div>

[Hyperlinks]:#

[Overview]: #overview
[Downloads]: #downloads
[Compatibility]: #compatibility
[Prerequisites]: #prerequisites
[Installation]: #installation
[Contributing]: #contributing
[Acknowledgments]: #acknowledgments

[Fortran 2018]: https://j3-fortran.org/doc/year/18/18-007r1.pdf 
[Arm]: https://www.arm.com

[OpenSHMEM]: http://www.openshmem.org/site/
[sourcery-institute logo]: http://www.sourceryinstitute.org/uploads/4/9/9/6/49967347/sourcery-logo-rgb-hi-rez-1.png
[OpenCoarrays]: http://www.opencoarrays.org
[ABI]: https://gcc.gnu.org/onlinedocs/gfortran/Coarray-Programming.html#Coarray-Programming
[MPI]: https://www.mpi-forum.org/
[GCC]: https://gcc.gnu.org
[gfortran]: https://gcc.gnu.org/wiki/GFortran
[MPICH]: https://www.mpich.org
[Sourcery, Inc.]: http://www.sourceryinstitute.org
[Google]: https://www.google.com
[CINECA]: https://www.cineca.it/en
[NERSC]: https://www.nersc.gov
[National Center for Atmospheric Research]: https://ncar.ucar.edu
[INSTALL.md]: ./INSTALL.md
[GASNet]: https://gasnet.lbl.gov
[CONTRIBUTING.md]: ./CONTRIBUTING.md
[GETTING_STARTED.md]: ./GETTING_STARTED.md
[Google Summer of Code]: https://www.google-melange.com/archive/gsoc/2014/orgs/gcc

[Issues]: https://github.com/sourceryinstitute/OpenCoarrays/issues
[Releases]: https://github.com/sourceryinstitute/OpenCoarrays/releases

[try OpenCoarrays online]: https://bit.ly/CAF-Binder
[notebook kernel]: https://github.com/sourceryinstitute/jupyter-CAF-kernel
[Binder]: https://mybinder.org
[Jupyter]: https://jupyter.org
[index.ipynb]: https://nbviewer.jupyter.org/github/sourceryinstitute/jupyter-CAF-kernel/blob/master/index.ipynb

[OC-on-repology]: https://repology.org/project/opencoarrays/versions
[repology-badge]: https://repology.org/badge/vertical-allrepos/opencoarrays.svg

[build img]: https://img.shields.io/travis/sourceryinstitute/OpenCoarrays.svg?style=flat-square "Build badge"
[CI Master Branch]: https://travis-ci.org/sourceryinstitute/OpenCoarrays?branch=master "View Travis-CI builds"
[download img]: https://img.shields.io/github/downloads/sourceryinstitute/OpenCoarrays/total.svg?style=flat-square "Download count badge"
[license img]: https://img.shields.io/badge/license-BSD--3-blue.svg?style=flat-square "BSD-3 License badge"
[release img]: https://img.shields.io/github/release/sourceryinstitute/OpenCoarrays.svg?style=flat-square "Latest release badge"
[pdf img]: https://img.shields.io/badge/PDF-README.md-6C2DC7.svg?style=flat-square "Download this readme as a PDF"
[twitter img]: https://img.shields.io/twitter/url/http/shields.io.svg?style=social
[Writing Fortran 2018 Today]: https://www.eventbrite.com/e/writing-fortran-2018-today-object-oriented-parallel-programming-tickets-48982176007
[Rouson et al. (2017)]: http://www.opencoarrays.org/uploads/6/9/7/4/69747895/a4-rouson.pdf

[default tweet]: https://twitter.com/intent/tweet?hashtags=HPC,Fortran,PGAS&related=zbeekman,gnutools,HPCwire,HPC_Guru,hpcprogrammer,SciNetHPC,DegenerateConic,jeffdotscience,travisci&text=Stop%20programming%20w%2F%20the%20%23MPI%20docs%20in%20your%20lap%2C%20try%20Coarray%20Fortran%20w%2F%20OpenCoarrays%20%26%20GFortran!&url=https%3A//github.com/sourceryinstitute/OpenCoarrays
