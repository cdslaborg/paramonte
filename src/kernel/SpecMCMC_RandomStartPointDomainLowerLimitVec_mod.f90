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

module SpecMCMC_RandomStartPointDomainLowerLimitVec_mod

    use Constants_mod, only: RK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecMCMC_RandomStartPointDomainLowerLimitVec_mod"

    real(RK), allocatable           :: RandomStartPointDomainLowerLimitVec(:) ! namelist input

    type                            :: RandomStartPointDomainLowerLimitVec_type
        real(RK), allocatable       :: Val(:)
        real(RK)                    :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set, checkForSanity, nullifyNameListVar
    end type RandomStartPointDomainLowerLimitVec_type

    interface RandomStartPointDomainLowerLimitVec_type
        module procedure            :: construct
    end interface RandomStartPointDomainLowerLimitVec_type

    private :: construct, set, checkForSanity, nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function construct(methodName) result(self)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: construct
#endif
        use Constants_mod, only: NULL_RK
        use String_mod, only: num2str
        implicit none
        character(*), intent(in)                        :: methodName
        type(RandomStartPointDomainLowerLimitVec_type)  :: self
        self%null = NULL_RK
        self%desc = &
        "RandomStartPointDomainLowerLimitVec represents the lower boundaries of the cubical domain from which the starting point(s) of &
        &the MCMC chain(s) will be initialized randomly (only if requested via the input variable randomStartPointRequested. &
        &This happens only when some or all of the elements of the input variable StartPoint are missing. &
        &In such cases, every missing value of input StartPoint will be set to the center point between RandomStartPointDomainLowerLimitVec &
        &and RandomStartPointDomainUpperLimit in the corresponding dimension. &
        &If RandomStartPointRequested=TRUE (or True, true, t, all case-insensitive), then the missing &
        &elements of StartPoint will be initialized to values drawn randomly from within the corresponding ranges specified by &
        &the input variable RandomStartPointDomainLowerLimitVec. &
        &As an input variable, RandomStartPointDomainLowerLimitVec is an ndim-dimensional vector of 64-bit real numbers, &
        &where ndim is the number of variables of the objective function. It is also possible to assign only select values of &
        &RandomStartPointDomainLowerLimitVec and leave the rest of the components to be assigned the default value. &
        &This is POSSIBLE ONLY when RandomStartPointDomainLowerLimitVec is defined inside the input file to "//methodName//". &
        &For example, having the following inside the input file, \n\n&
        &    RandomStartPointDomainLowerLimitVec(3:5) = -100\n\n&
        &            will only set the lower limits of the third, fourth, and the fifth dimensions to -100, or,\n\n&
        &    RandomStartPointDomainLowerLimitVec(1) = -100, RandomStartPointDomainLowerLimitVec(2) = -1.e6 \n\n&
        &            will set the lower limit on the first dimension to -100, and 1.e6 on the second dimension, or,\n\n&
        &    RandomStartPointDomainLowerLimitVec = 3*-2.5e100\n\n&
        &            will only set the lower limits on the first, second, and the third dimensions to -2.5*10^100, while the rest of &
                    &the lower limits for the missing dimensions will be automatically set to the default value.\n\n&
        &The default values for all elements of RandomStartPointDomainLowerLimitVec are taken from the corresponding values in the input &
        &variable domainLowerLimitVec."
    end function construct

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine nullifyNameListVar(self,nd)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        use Constants_mod, only: IK
        implicit none
        class(RandomStartPointDomainLowerLimitVec_type), intent(in) :: self
        integer(IK), intent(in)                                     :: nd
        if (allocated(randomStartPointDomainLowerLimitVec)) deallocate(randomStartPointDomainLowerLimitVec)
        allocate(randomStartPointDomainLowerLimitVec(nd))
        randomStartPointDomainLowerLimitVec(:) = self%null
    end subroutine nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine set(self, randomStartPointDomainLowerLimitVec)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: set
#endif
        use Constants_mod, only: IK, RK
        implicit none
        class(RandomStartPointDomainLowerLimitVec_type), intent(inout)  :: self
        real(RK), intent(in), optional                                  :: randomStartPointDomainLowerLimitVec(:)
        if (present(randomStartPointDomainLowerLimitVec)) self%Val = randomStartPointDomainLowerLimitVec
    end subroutine set

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine checkForSanity(self,Err,methodName,SpecBase,randomStartPointRequested)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use SpecBase_mod, only: SpecBase_type
        use Constants_mod, only: RK
        use Err_mod, only: Err_type
        use String_mod, only: num2str
        implicit none
        class(RandomStartPointDomainLowerLimitVec_type), intent(inout)  :: self
        type(Err_type), intent(inout)                                   :: Err
        type(SpecBase_type), intent(in)                                 :: SpecBase
        character(*), intent(in)                                        :: methodName
        logical     , intent(in)                                        :: randomStartPointRequested
        character(*), parameter                                         :: PROCEDURE_NAME = "@checkForSanity()"
        integer                                                         :: i
        do i = 1,size(self%Val(:))

            if (self%Val(i)==self%null) self%Val(i) = SpecBase%DomainLowerLimitVec%Val(i)

            ! check if the domain is set when random start point is requested

            if ( randomStartPointRequested .and. self%Val(i)==SpecBase%DomainLowerLimitVec%def ) then
                Err%occurred = .true.
                Err%msg =   Err%msg // &
                            MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                            &You have requested a random start point by setting randomStartPointRequested to TRUE while the &
                            &element #"//num2str(i)//" of RandomStartPointDomainLowerLimitVec has not been preset to a finite value. &
                            &This information is essential otherwise, how could the sampler draw points randomly from within an unspecified domain?\n\n"
            end if

            ! check if the random start point domain is within the boundaries of the domain of the target.

            if ( self%Val(i)<SpecBase%DomainLowerLimitVec%Val(i) ) then
                Err%occurred = .true.
                Err%msg =   Err%msg // &
                            MODULE_NAME // PROCEDURE_NAME // ": Error occurred. The component " // num2str(i) // &
                            " of the variable RandomStartPointDomainLowerLimitVec (" // num2str(self%Val(i)) // &
                            ") cannot be smaller than the corresponding component of the variable domainLowerLimitVec (" // &
                            num2str(SpecBase%DomainLowerLimitVec%Val(i)) // "). If you don't know &an appropriate value to " // & ! LCOV_EXCL_LINE
                            "set for RandomStartPointDomainLowerLimitVec, drop it from the input list. " // methodName // " will &
                            &automatically assign an appropriate value to it.\n\n"
            end if

        end do
        deallocate(randomStartPointDomainLowerLimitVec)
    end subroutine checkForSanity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module SpecMCMC_RandomStartPointDomainLowerLimitVec_mod ! LCOV_EXCL_LINE