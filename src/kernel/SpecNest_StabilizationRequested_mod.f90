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
!> This module contains the classes and procedures for setting up the `stabilizationRequested` attribute of samplers of class [ParaNest_type](@ref paranest_mod::paranest_type).
!> For more information, see the description of this attribute in the body of the module.
!> \author Amir Shahmoradi

module SpecNest_StabilizationRequested_mod

    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecNest_StabilizationRequested_mod"

    logical                         :: stabilizationRequested ! namelist input

    type                            :: StabilizationRequested_type
        logical                     :: val
        logical                     :: def
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set, nullifyNameListVar
    end type StabilizationRequested_type

    interface StabilizationRequested_type
        module procedure            :: construct
    end interface StabilizationRequested_type

    private :: construct, set, nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function construct(methodName) result(self)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: construct
#endif
        use String_mod, only: log2str
        implicit none
        character(*), intent(in)            :: methodName
        type(StabilizationRequested_type)   :: self
        self%def = .true.
        self%desc = &
        "A logical (boolean) variable. If set to true or .true. or TRUE or .t. from within an input file, or to true from within &
        &MATLAB or True from within Python, then, the volume minimization in the clustering process of the ellipsoidal rejection &
        &sampling will be stabilized to ensure ellipsoidal rejection sampling convergence without clustering failures. &
        &For all other sampling methods, the value of this simulation specification is ignored. This is a low-level simulation &
        &specification whose value should not be modified from the default without knowing the consequences. &
        &The default value is " // log2str(self%def) // "."
    end function construct

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine nullifyNameListVar(self)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(StabilizationRequested_type), intent(in)    :: self
        stabilizationRequested = self%def
    end subroutine nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure subroutine set(self,stabilizationRequested)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: set
#endif
        implicit none
        class(StabilizationRequested_type), intent(inout)   :: self
        logical, intent(in)                                 :: stabilizationRequested
        self%val = stabilizationRequested
    end subroutine set

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module SpecNest_StabilizationRequested_mod ! LCOV_EXCL_LINE