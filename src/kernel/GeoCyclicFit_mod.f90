!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!
!!!!   MIT License
!!!!
!!!!   ParaMonte: plain powerful parallel Monte Carlo library.
!!!!
!!!!   Copyright (C) 2012-present, The Computational Data Science Lab
!!!!
!!!!   This file is part of the ParaMonte library.
!!!!
!!!!   Permission is hereby granted, free of charge, to any person obtaining a
!!!!   copy of this software and associated documentation files (the "Software"),
!!!!   to deal in the Software without restriction, including without limitation
!!!!   the rights to use, copy, modify, merge, publish, distribute, sublicense,
!!!!   and/or sell copies of the Software, and to permit persons to whom the
!!!!   Software is furnished to do so, subject to the following conditions:
!!!!
!!!!   The above copyright notice and this permission notice shall be
!!!!   included in all copies or substantial portions of the Software.
!!!!
!!!!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!!!!   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
!!!!   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
!!!!   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
!!!!   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
!!!!   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
!!!!   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!!!!
!!!!   ACKNOWLEDGMENT
!!!!
!!!!   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
!!!!   As per the ParaMonte library license agreement terms, if you use any parts of
!!!!   this library for any purposes, kindly acknowledge the use of ParaMonte in your
!!!!   work (education/research/industry/development/...) by citing the ParaMonte
!!!!   library as described on this page:
!!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \brief
!> This module contains the classes and procedures for fitting a Cyclic Geometric distribution.
!> The functions and procedures contained in this module were originally part of the [Statistics_mod](@ref statistics_mod)
!> module. However, they were moved out to bypass the segmentation fault error with internal functions when the library is 
!> compiled and run on Microsoft Windows Subsystem for Linux using the GNU Fortran compiler.
!>
!> \author Amir Shahmoradi

module GeoCyclicFit_mod

    use Optimization_mod, only: PowellMinimum_type
    use Constants_mod, only: IK, RK
    implicit none


    character(len=*), parameter :: MODULE_NAME = "@GeoCyclicFit_mod"

#if defined OS_IS_WSL
    private                     :: getSumDistSq
    integer(IK)                 :: numTrial_WSL         !< This madness bypasses the Microsoft Subsystem for Linux Internal Function call GFortran Segmentation Fault error.
    integer(IK)                 :: maxNumTrial_WSL      !< This madness bypasses the Microsoft Subsystem for Linux Internal Function call GFortran Segmentation Fault error.
    integer(IK) , allocatable   :: SuccessStep_WSL(:)   !< This madness bypasses the Microsoft Subsystem for Linux Internal Function call GFortran Segmentation Fault error.
    real(RK)    , allocatable   :: LogCount_WSL(:)      !< This madness bypasses the Microsoft Subsystem for Linux Internal Function call GFortran Segmentation Fault error.
#endif

    type :: GeoCyclicFit_type
        type(PowellMinimum_type) :: PowellMinimum
    contains
        procedure, nopass :: fit => fitGeoCyclicLogPDF
    end type GeoCyclicFit_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return a fit of the Cyclic Geometric distribution PDF to the input natural logarithm of a sequence of Counts,
    !> the `i`th element of which represents the number of successes after `SuccessStep(i)` tries in a Bernoulli trail.
    !>
    !> \param[in]   maxNumTrial :   The maximum number of trails possible. After `maxNumTrial` tries,
    !                               the Geometric distribution restarts from index `1`.
    !> \param[in]   numTrial    :   The length of the array `SuccessStep` and `LogCount`.
    !>                              Note that `numTrial < maxNumTrial` must hold.
    !> \param[in]   SuccessStep :   A vector of length `(1:numTrial)` of integers that represent
    !>                              the steps at which the Bernoulli successes occur.
    !> \param[in]   LogCount    :   A real-valued vector of length `(1:numTrial)` representing the natural logarithms of the
    !>                              counts of success at the corresponding Bernoulli trials specified by elements of `SuccessStep`.
    !>
    !> \return
    !> `PowellMinimum`          :   An object of class [PowellMinimum_type](@ref optimization_mod::powellminimum_type) containing
    !>                              the best-fit successProb and the normalization constant of the fit in the vector component `xmin`.
    !>
    !> \warning
    !> Any value of SuccessStep must be an integer numbers between `1` and `maxNumTrial`.
    !> The onus is on the user to ensure this condition holds.
    !>
    !> \todo
    !> Update: Amir Shahmoradi, Sunday Nov 29, 2020, 11:19 pm, Dallas, TX
    !> The current implementation of the objective function relies on the definitions of module variables.
    !> Although inefficient and ugly, this was necessary to resolve the viscous Segmentation Fault error
    !> that happens with internal function calls on Windows Subsystem for Linux Ubuntu with GFortran.
    !> Once this error of unknown origin is resolved, the function [getSumDistSq](@ref getsumdistsq)
    !> must be converted back to an internal function within [fitGeoCyclicLogPDF](@ref fitgeocycliclogpdf)
    !> and subsequently, all module variables must be removed.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    function fitGeoCyclicLogPDF(maxNumTrial, numTrial, SuccessStep, LogCount) result(PowellMinimum)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: fitGeoLogPDF
#endif
        use Optimization_mod, only: PowellMinimum_type
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK) , intent(in)    :: maxNumTrial
        integer(IK) , intent(in)    :: numTrial
        integer(IK) , intent(in)    :: SuccessStep(numTrial)
        real(RK)    , intent(in)    :: LogCount(numTrial)
        type(PowellMinimum_type)    :: PowellMinimum

        real(RK)                    :: BestFitSuccessProbNormFac(2) ! vector of the two parameters
        real(RK)    , parameter     :: SUCCESS_PROB_INIT_GUESS = 0.23_RK
        real(RK)    , parameter     :: FISHER_TRANS_SUCCESS_PROB_INIT_GUESS = atanh(2*(SUCCESS_PROB_INIT_GUESS - 0.5_RK))

#if defined OS_IS_WSL
        numTrial_WSL = numTrial
        maxNumTrial_WSL = maxNumTrial
        SuccessStep_WSL = SuccessStep
        LogCount_WSL = LogCount
#endif

        ! do Fisher transformation to make the limits infinity.
        BestFitSuccessProbNormFac = [FISHER_TRANS_SUCCESS_PROB_INIT_GUESS, 0._RK] ! LogCount(1)]

        PowellMinimum = PowellMinimum_type  ( ndim = 2_IK &
                                            , getFuncMD = getSumDistSq &
                                            , StartVec = BestFitSuccessProbNormFac &
                                            )
        if (PowellMinimum%Err%occurred) return
        PowellMinimum%xmin(1) = 0.5_RK * tanh(PowellMinimum%xmin(1)) + 0.5_RK ! reverse Fisher-transform

#if !defined OS_IS_WSL
    contains

        !! doxygen has problems digesting the documentation of Fortran internal functions.
        !!> \brief
        !!>
        !!> \param[in]   ndim                            :   The length of the input vector `successProbFisherTransNormFac`.
        !!> \param[in]   successProbFisherTransNormFac   :   The length of the input vector `successProbFisherTransNormFac`.
        !!>
        !!> \return
        !!> `sumDistSq` : The sum of distances squared.
        !!>
        !!> \remark
        !!> Although `successProbFisherTransNormFac` is a vector on input, it is expected to have a length of one at all times.
        !!> This is solely to fullfile the interface restrictions of [PowellMinimum_type](@ref optimization_mod::powellminimum_type).
        pure function getSumDistSq(ndim,successProbFisherTransNormFac) result(sumDistSq)
            use Statistics_mod, only: getLogProbGeoCyclic
            !use Constants_mod, only: IK, RK
            implicit none
            integer(IK) , intent(in)    :: ndim
            real(RK)    , intent(in)    :: successProbFisherTransNormFac(ndim)
            real(RK)                    :: sumDistSq, successProb
            successProb = 0.5_RK * tanh(successProbFisherTransNormFac(1)) + 0.5_RK ! reverse Fisher-transform
            !sumDistSq = sum( (LogCount - getGeoLogPDF(successProb=successProb,seqLen=numTrial) - successProbFisherTransNormFac(2) )**2 )
            sumDistSq = sum(    ( LogCount &
                                - getLogProbGeoCyclic(successProb=successProb, maxNumTrial=maxNumTrial, numTrial=numTrial, SuccessStep=SuccessStep) &
                                - successProbFisherTransNormFac(2) &
                                )**2 &
                            )
        end function getSumDistSq
#endif
    end function fitGeoCyclicLogPDF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if defined OS_IS_WSL
    !> \brief
    !> Return the sum of the distances squared from the current fit.
    !>
    !> \param[in]   ndim                            :   The length of the input vector `successProbFisherTransNormFac`.
    !> \param[in]   successProbFisherTransNormFac   :   The length of the input vector `successProbFisherTransNormFac`.
    !>
    !> \return
    !> `sumDistSq` : The sum of distances squared.
    !>
    !> \remark
    !> Although `successProbFisherTransNormFac` is a vector on input, it is expected to have a length of one at all times.
    !> This is solely to fullfile the interface restrictions of [PowellMinimum_type](@ref optimization_mod::powellminimum_type).
    pure function getSumDistSq(ndim,successProbFisherTransNormFac) result(sumDistSq)
        use Statistics_mod, only: getLogProbGeoCyclic
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK) , intent(in)    :: ndim
        real(RK)    , intent(in)    :: successProbFisherTransNormFac(ndim)
        real(RK)                    :: sumDistSq, successProb
        successProb = 0.5_RK * tanh(successProbFisherTransNormFac(1)) + 0.5_RK ! reverse Fisher-transform
        !sumDistSq = sum( (LogCount - getGeoLogPDF(successProb=successProb,seqLen=numTrial) - successProbFisherTransNormFac(2) )**2 )
        sumDistSq = sum(    ( LogCount_WSL &
                            - getLogProbGeoCyclic(successProb=successProb, maxNumTrial=maxNumTrial_WSL, numTrial=numTrial_WSL, SuccessStep=SuccessStep_WSL) &
                            - successProbFisherTransNormFac(2) &
                            )**2 &
                        )
    end function getSumDistSq
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    !> \brief
!    !> Return a fit of the Geometric distribution PDF to the input natural logarithm of a sequence of Counts,
!    !> the `i`th element of which represents the number of successes after `SuccessStep(i)` tries in a Bernoulli trail.
!    !>
!    !> \param[in]   numTrial    :   The number of trials. The length of the input vector `LogCount`.
!    !> \param[in]   SuccessStep :   The vector of trials of length `numTrial` at which the first success happens.
!    !> \param[in]   LogCount    :   A vector of real values representing the natural logarithms of the counts
!    !>                              of success at each Bernoulli trial, sequentially, from `1` to `numTrial`.
!    !>
!    !> \return
!    !> `PowellMinimum`  :   An object of class [PowellMinimum_type](@ref optimization_mod::powellminimum_type) containing
!    !>                      the best-fit successProb and the normalization constant of the fit in the vector component `xmin`.
!    !>
!    !> \author
!    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
!    function fitGeoLogPDF_old(numTrial, SuccessStep, LogCount) result(PowellMinimum)
!#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
!        !DEC$ ATTRIBUTES DLLEXPORT :: fitGeoLogPDF_old
!#endif
!        use Optimization_mod, only: PowellMinimum_type
!        use Constants_mod, only: IK, RK, POSINF_RK, NEGINF_RK
!        implicit none
!        integer(IK) , intent(in)    :: numTrial
!        integer(IK) , intent(in)    :: SuccessStep(numTrial)
!        real(RK)    , intent(in)    :: LogCount(numTrial)
!        type(PowellMinimum_type)    :: PowellMinimum
!
!        real(RK)                    :: BestFitSuccessProbNormFac(2) ! vector of the two parameters
!        real(RK)    , parameter     :: SUCCESS_PROB_INIT_GUESS = 0.23_RK
!        real(RK)    , parameter     :: FISHER_TRANS_SUCCESS_PROB_INIT_GUESS = atanh(2*(SUCCESS_PROB_INIT_GUESS - 0.5_RK))
!
!        ! do Fisher transformation to make the limits infinity
!        BestFitSuccessProbNormFac = [FISHER_TRANS_SUCCESS_PROB_INIT_GUESS, LogCount(1)]
!
!        PowellMinimum = PowellMinimum_type  ( ndim = 2_IK &
!                                            , getFuncMD = getSumDistSq &
!                                            , StartVec = BestFitSuccessProbNormFac &
!                                            )
!        if (PowellMinimum%Err%occurred) return
!        PowellMinimum%xmin(1) = 0.5_RK * tanh(PowellMinimum%xmin(1)) + 0.5_RK ! reverse Fisher-transform
!
!    contains
!
!        function getSumDistSq(ndim,successProbFisherTransNormFac) result(sumDistSq)
!            !use Constants_mod, only: IK, RK
!            implicit none
!            integer(IK) , intent(in)    :: ndim
!            real(RK)    , intent(in)    :: successProbFisherTransNormFac(ndim)
!            real(RK)                    :: sumDistSq, successProb
!            successProb = 0.5_RK*tanh(successProbFisherTransNormFac(1)) + 0.5_RK ! reverse Fisher-transform
!            !sumDistSq = sum( (LogCount - getGeoLogPDF(successProb=successProb,seqLen=numTrial) - successProbFisherTransNormFac(2) )**2 )
!            sumDistSq = sum(    ( LogCount &
!                                - numTrial * successProbFisherTransNormFac(2) &
!                                - getLogProbGeo(numTrial = numTrial, SuccessStep = SuccessStep, successProb = successProb) &
!                                )**2 &
!                            )
!        end function getSumDistSq
!
!    end function fitGeoLogPDF_old

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module GeoCyclicFit_mod ! LCOV_EXCL_LINE