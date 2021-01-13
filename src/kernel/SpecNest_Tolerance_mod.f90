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
!> This module contains the classes and procedures for setting up the `tolerance` attribute of samplers of class [ParaNest_type](@ref paranest_mod::paranest_type).
!> For more information, see the description of this attribute in the body of the module.
!> \author Amir Shahmoradi

module SpecNest_Tolerance_mod

    use Constants_mod, only: RK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecNest_Tolerance_mod"

    real(RK)                        :: tolerance ! namelist input

    type                            :: Tolerance_type
        real(RK)                    :: val
        real(RK)                    :: def
        real(RK)                    :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set, checkForSanity, nullifyNameListVar
    end type Tolerance_type

    interface Tolerance_type
        module procedure            :: construct
    end interface Tolerance_type

    private :: construct, set, nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function construct(methodName) result(self)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: construct
#endif
        use Constants_mod, only: NULL_RK
        use String_mod, only: num2str
        implicit none
        character(*), intent(in)    :: methodName
        type(Tolerance_type)        :: self
        self%def = 0.5_RK
        self%null = NULL_RK
        self%desc = &
        "tolerance is a 64-bit positive real number representing the threshold below which the integration of the objective function &
        &by the " //methodName// " sampler is assumed to have converged. Typical values are around 0.1-1. Avoid setting this number to &
        &extremely tiny values on the order of 1e-6 or smaller, otherwise neither the integration result is guaranteed to converge nor &
        &nor the sampler will ever stop the simulation. &
        &The default value for tolerance is "//num2str(self%def)//"."
    end function construct

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine nullifyNameListVar(self)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(Tolerance_type), intent(in)   :: self
        tolerance = self%null
    end subroutine nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure subroutine set(self,tolerance)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: set
#endif
        use Constants_mod, only: RK
        implicit none
        class(Tolerance_type), intent(inout)    :: self
        real(RK), intent(in)                    :: tolerance
        self%val = tolerance
        if (self%val==self%null) self%val = self%def
    end subroutine set

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure subroutine checkForSanity(self,Err,methodName)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Constants_mod, only: RK
        use Err_mod, only: Err_type
        use String_mod, only: num2str
        implicit none
        class(Tolerance_type), intent(in)   :: self
        type(Err_type), intent(inout)       :: Err
        character(*), intent(in)            :: methodName
        character(*), parameter             :: PROCEDURE_NAME = "@checkForSanity()"
        if ( self%val <= 0._RK ) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &The input variable tolerance (" // num2str(self%val) // &
                        ") cannot be less than 0. If you are not sure of the appropriate value for tolerance, drop it &
                        &from the input list. " // methodName // " will automatically assign an appropriate value to it.\n\n"
        end if
    end subroutine checkForSanity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module SpecNest_Tolerance_mod ! LCOV_EXCL_LINE