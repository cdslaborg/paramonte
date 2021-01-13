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
!!!!       https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \brief
!> This module contains the classes and procedures for setting up the `liveSampleSize` attribute of samplers of class [ParaNest_type](@ref paranest_mod::paranest_type).
!> For more information, see the description of this attribute in the body of the module.
!> \author Amir Shahmoradi

module SpecNest_LiveSampleSize_mod

    use Constants_mod, only: IK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecNest_LiveSampleSize_mod"

    integer(IK)                     :: liveSampleSize ! namelist input

    type                            :: LiveSampleSize_type
        integer(IK)                 :: val
        integer(IK)                 :: def
        integer(IK)                 :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set, checkForSanity, nullifyNameListVar
    end type LiveSampleSize_type

    interface LiveSampleSize_type
        module procedure            :: construct
    end interface LiveSampleSize_type

    private :: construct, set, checkForSanity, nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function construct(nd,methodName) result(self)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: construct
#endif
        use Constants_mod, only: IK, NULL_IK, POSINF_IK
        use String_mod, only: num2str
        implicit none
        integer(IK), intent(in)     :: nd
        character(*), intent(in)    :: methodName
        type(LiveSampleSize_type)   :: self
        self%def  = 2000_IK ! @todo: perhaps an intelligent setting strategy here would be useful.
        self%null = NULL_IK
        self%desc = &
        "liveSampleSize is a positive 32bit integer representing the number of live (active) points that are initially sampled &
        &uniformly from the domain of the objective function. New points will be subsequently added to this sample throughout the &
        &simulation and the lowest-value will be removed such that the size of the live (active) sample of points at any stage &
        &during the simulation remains fixed as specified by the input simulation specification liveSampleSize. &
        &The value of liveSampleSize must be at least larger than nd + 1, where nd is number of dimensions of the domain of the &
        &objective function. It should be preferably an integer on the order of hundreds or thousands. The higher the value of &
        &liveSampleSize is, the better the odd of capturing all modes of the objective function will be. &
        &Essentially, the number of live (active) points acts as a random mesh on the integration domain of objective function. &
        &Therefore, a larger live sample size means a finer mesh for the integration, thus yielding more accurate results. &
        &In general, the higher the number of dimensions of the domain of the objective function, the larger the value of &
        &liveSampleSize should be. The default value is " // num2str(self%def) // "."
    end function construct

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine nullifyNameListVar(self)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(LiveSampleSize_type), intent(in)  :: self
        liveSampleSize = self%null
    end subroutine nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure subroutine set(self,liveSampleSize)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: set
#endif
        use Constants_mod, only: IK
        implicit none
        class(LiveSampleSize_type), intent(inout)  :: self
        integer(IK), intent(in)                         :: liveSampleSize
        self%val = liveSampleSize
        if ( self%val==self%null ) self%val = self%def
    end subroutine set

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure subroutine checkForSanity(self,Err,methodName,nd)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Constants_mod, only: IK
        use Err_mod, only: Err_type
        use String_mod, only: num2str
        implicit none
        class(LiveSampleSize_type), intent(in)  :: self
        character(*), intent(in)                :: methodName
        integer(IK), intent(in)                 :: nd
        type(Err_type), intent(inout)           :: Err
        character(*), parameter                 :: PROCEDURE_NAME = "@checkForSanity()"
        if (self%val < nd + 1) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &The input requested value for liveSampleSize (" // num2str(self%val) // ") &
                        &can not be negative. If you are not sure of the appropriate value for liveSampleSize, drop it &
                        &from the input list. " // methodName // " will automatically assign an appropriate value to it.\n\n"
        end if
    end subroutine checkForSanity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module SpecNest_LiveSampleSize_mod ! LCOV_EXCL_LINE