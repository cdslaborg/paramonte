/*
####################################################################################################################################
####################################################################################################################################
####                                                                                                                            ####
####    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           ####
####                                                                                                                            ####
####    Copyright (C) 2012-present, The Computational Data Science Lab                                                          ####
####                                                                                                                            ####
####    This file is part of the ParaMonte library.                                                                             ####
####                                                                                                                            ####
####    LICENSE                                                                                                                 ####
####                                                                                                                            ####
####       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          ####
####                                                                                                                            ####
####################################################################################################################################
####################################################################################################################################
*/
///
/// \file pm_sampling.h
///
/// \defgroup pm_sampling pm_sampling
///
/// \brief
/// This header file contains function interfaces to the ParaMonte library Monte Carlo samplers and integrators.
///
/// \final
///
/// \author
/// \AmirShahmoradi, Monday 00:01 AM, January 1, 2018, Institute for Computational Engineering and Sciences, University of Texas Austin
///
#include <stdint.h>
#ifndef pm_sampling
#define pm_sampling
///
/// \ingroup pm_sampling
///
/// \defgroup runParaDRAM runParaDRAM
///
/// \brief
/// Generate and return a non-zero value (`1`) if the procedure fails to fully accomplish the task of generating
/// a Monte Carlo sample of the specified input mathematical objective function, otherwise, return `0`.
///
/// \details
/// This interface group is the entry point to all **C-style interfaces** to the ParaDRAM samplers of mathematical density functions.<br>
/// Although the procedures of this generic interface return a single scalar of type `int32_t`, the procedures generate
/// massive amounts of information about each simulation which are stored in appropriate external hard drive files.<br>
///
/// See,
/// <ol>
///     <li>    [this generic documentation page](\pmdoc_usage_sampling/paradram/output/)
///             for more information on the generated output files for samplings performed using the [ParaDRAM](@ref runParaDRAM) sampler.
/// </ol>
///
/// \param[in]  getLogFunc  :   The input user-specified procedure pointer to the natural logarithm of the target density function.<br>
///                             On input, the function must take an input vector `state` of size `ndim` of type floating-point of kind \C_RKALL
///                             representing a state (point) from within the domain of the user-specified target density function whose function value must be returned.<br>
///                             On output, the user-specified procedure `getLogFunc()` must return the function value corresponding to the input `state[ndim]`.<br>
///                             The following illustrate the generic interface of input function pointer `getLogFunc(state, ndim)`,
///                             \code{.F90}
///                                 REAL getLogFunc(REAL[] state, int32_t ndim);
///                             \endcode
///                             where `REAL` can be a floating-point type of kind \C_RKALL for the corresponding varying-precision sampler interfaces:<br>
///                             <ol>
///                                 <li>    [runParaDRAMF](@ref runParaDRAMF),
///                                 <li>    [runParaDRAMD](@ref runParaDRAMD),
///                                 <li>    [runParaDRAML](@ref runParaDRAML).
///                             </ol>
/// \param[in]  ndim        :   The input scalar constant of type `int32_t` representing the number of dimensions of the domain of the objective function.
/// \param[in]  input       :   The input scalar pointer of type `char` representing the null-terminated C string
///                             that contains either,
///                             <ol>
///                                 <li>    the path to an external input file containing the namelist group of ParaDRAM sampler specifications
///                                         as outlined in the corresponding page of [ParaMonte library generic documentation website](\pmdoc_usage_sampling/paradram/specifications/).<br>
///                                 <li>    the namelist group of ParaDRAM sampler specifications as the can appear in an external input specification file.<br>
///                             </ol>
///                             While all input simulation specifications are optional, it is highly recommended to pay
///                             attention to the default settings of the domain boundaries and sampler starting point.<br>
///                             (**optional**. It is considered as missing if set to `NULL`.)
///
/// \return
/// `stat`                  :   The output scalar of type `int32_t` that is `0` **if and only if**
///                             the sampler succeeds in sampling the specified density function.<br>
///
/// \interface{runParaDRAM}
/// \code{.c}
///
///     #include "pm_sampling.h"
///
///     stat = runParaDRAMF(getLogFunc, ndim, input)
///     stat = runParaDRAMD(getLogFunc, ndim, input)
///     stat = runParaDRAML(getLogFunc, ndim, input)
///
/// \endcode
///
/// \warning
/// Beware that the definition of extended precision real type `long double` is compiler and platform dependent
/// making the use of `long double` with precompiled ParaMonte libraries problematic and non-functional.<br>
///
/// \warning
/// The condition `0 < ndim` must hold for the corresponding input arguments.<br>
/// \vericon
///
/// \example{mvn}
/// \include{lineno} example/pm_sampling/mvn/main.c
/// \inputfile{mvn, input.nml, example specifications namelist input file}
/// \include{lineno} example/pm_sampling/mvn/input.nml
/// \cmakefile{mvn}
/// \include{lineno} example/pm_sampling/mvn/CMakeLists.txt
/// \cmakerun{mvn}
/// \postproc{mvn}
/// \include{lineno} example/pm_sampling/mvn/main.py
/// \vis{mvn}
/// \image html example/pm_sampling/mvn/mvn.traceplot.png width=700
/// \image html example/pm_sampling/mvn/mvn.scatterplot.png width=700
/// \image html example/pm_sampling/mvn/mvn.proposalAdaptation.png width=700
///
/// \example{himmelblau}
/// \include{lineno} example/pm_sampling/himmelblau/main.c
/// \inputfile{himmelblau, input.nml, example specifications namelist input file}
/// \include{lineno} example/pm_sampling/himmelblau/input.nml
/// \cmakefile{himmelblau}
/// \include{lineno} example/pm_sampling/himmelblau/CMakeLists.txt
/// \cmakerun{himmelblau}
/// \postproc{himmelblau}
/// \include{lineno} example/pm_sampling/himmelblau/main.py
/// \vis{himmelblau}
/// \image html example/pm_sampling/himmelblau/himmelblau.traceplot.png width=700
/// \image html example/pm_sampling/himmelblau/himmelblau.scatterplot.png width=700
/// \image html example/pm_sampling/himmelblau/himmelblau.proposalAdaptation.png width=700
///
/// \test
/// [test_pm_sampling](@ref test_pm_sampling)
///
/// \bug
/// \status \unresolved
/// \source \ifx{2024.0.0 20231017}
/// \desc
/// \ifx This interface yields a segmentation fault error for all `real` types supported when linked with the ParaMonte library built with \ifx.
/// \remedy
/// For now, only \ifort will be used.<br>
///
/// \bug
/// \status \unresolved
/// \source \ifort{2021.11.0 20231010}
/// \desc
/// The `runParaDRAML` interface for `long double` yields a segmentation fault error.
/// \remedy
/// None as of today.<br>
///
/// \todo
/// \pvhigh
/// The current tests for `long double` interface fail with \ifx and \icx compilers.<br>
/// Apparently, there is a mismatch between `long double` and `c_long_double` storage.<br>
/// This issue does not exist with GNU compilers, although the GNU definition of `long double`
/// appears to yield incorrect values in some calculations, e.g., in `isFailedGeomCyclicFit()` of the ParaMonte library.<br>
///
/// \final{runParaDRAM}
///
/// \author
/// \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin
///

/// \ingroup runParaDRAM
int32_t runParaDRAML( long double (*getLogFunc)(long double state[], int32_t ndim)
                    , const int32_t ndim
                    , const char* input
                    );
/// \ingroup runParaDRAM
int32_t runParaDRAMD( double (*getLogFunc)(double state[], int32_t ndim)
                    , const int32_t ndim
                    , const char* input
                    );
/// \ingroup runParaDRAM
int32_t runParaDRAMF( float (*getLogFunc)(float state[], int32_t ndim)
                    , const int32_t ndim
                    , const char* input
                    );
#endif