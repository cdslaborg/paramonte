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

module SpecMCMC_StartPointVec_mod

    use Constants_mod, only: RK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecMCMC_StartPointVec_mod"

    real(RK), allocatable           :: startPointVec(:) ! namelist input

    type                            :: StartPointVec_type
        real(RK), allocatable       :: Val(:)
        real(RK)                    :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setStartPointVec, checkForSanity, nullifyNameListVar
    end type StartPointVec_type

    interface StartPointVec_type
        module procedure            :: constructStartPointVec
    end interface StartPointVec_type

    private :: constructStartPointVec, setStartPointVec, checkForSanity, nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function constructStartPointVec() result(StartPointVecObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructStartPointVec
#endif
        use Constants_mod, only: NULL_RK
        use String_mod, only: num2str
        implicit none
        type(StartPointVec_type) :: StartPointVecObj
        StartPointVecObj%null   = NULL_RK
        StartPointVecObj%desc   = &
        "startPointVec is a 64bit real-valued vector of length ndim (the dimension of the domain of the input objective function). &
        &For every element of startPointVec that is not provided as input, the default value will be the center of the domain of &
        &startPointVec as specified by domainLowerLimitVec and domainUpperLimitVec input variables. &
        &If the input variable randomStartPointRequested=TRUE (or true or t, all case-insensitive), then the missing &
        &elements of startPointVec will be initialized to values drawn randomly from within the corresponding &
        &ranges specified by the input variables randomStartPointDomainLowerLimitVec and randomStartPointDomainUpperLimitVec."
    end function constructStartPointVec

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine nullifyNameListVar(StartPointVecObj,nd)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        use Constants_mod, only: IK
        implicit none
        class(StartPointVec_type), intent(in)   :: StartPointVecObj
        integer(IK), intent(in)                 :: nd
        if (allocated(startPointVec)) deallocate(startPointVec)
        allocate(startPointVec(nd), source = StartPointVecObj%null)
    end subroutine nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setStartPointVec (StartPointVecObj,startPointVec &
                                ,randomStartPointDomainLowerLimitVec,randomStartPointDomainUpperLimitVec,randomStartPointRequested &
                                ,domainLowerLimitVec,domainUpperLimitVec)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setStartPointVec
#endif
        use Constants_mod, only: IK, RK
        implicit none
        class(StartPointVec_type), intent(inout)    :: StartPointVecObj
        real(RK), intent(in)                        :: startPointVec(:)
        real(RK), intent(in)                        :: randomStartPointDomainLowerLimitVec(:), randomStartPointDomainUpperLimitVec(:)
        real(RK), intent(in)                        :: domainLowerLimitVec(:), domainUpperLimitVec(:)
        logical, intent(in)                         :: randomStartPointRequested
        real(RK)                                    :: unifrnd
        integer(IK)                                 :: i
        StartPointVecObj%Val = startPointVec
        do i = 1, size(startPointVec)
            if (startPointVec(i)==StartPointVecObj%null) then
                if (randomStartPointRequested) then
                    call random_number(unifrnd)
                    StartPointVecObj%Val(i) = randomStartPointDomainLowerLimitVec(i) + unifrnd * (randomStartPointDomainUpperLimitVec(i)-randomStartPointDomainLowerLimitVec(i))
                else
                    StartPointVecObj%Val(i) = 0.5_RK * ( domainLowerLimitVec(i) + domainUpperLimitVec(i) )
                end if
            end if
        end do
    end subroutine setStartPointVec

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine checkForSanity(StartPointVecObj,Err,methodName,DomainLowerLimitVec,DomainUpperLimitVec)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        use Err_mod, only: Err_type
        implicit none
        class(StartPointVec_type), intent(in)   :: StartPointVecObj
        real(RK), intent(in)                    :: DomainLowerLimitVec(:), DomainUpperLimitVec(:)
        character(*), intent(in)                :: methodName
        type(Err_type), intent(inout)           :: Err
        character(*), parameter                 :: PROCEDURE_NAME = "@checkForSanity()"
        integer(IK)                             :: i
        do i = 1, size(StartPointVecObj%Val)
            if ( StartPointVecObj%Val(i)<DomainLowerLimitVec(i) .or. StartPointVecObj%Val(i)>DomainUpperLimitVec(i) ) then
                Err%occurred = .true.
                Err%msg =   Err%msg // &
                            MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                            &The input requested value for the component " // num2str(i) // " of the vector startPointVec (" // &
                            num2str(StartPointVecObj%Val(i)) // ") must be within the range of the sampling Domain defined &
                            &in the program: (" &
                            // num2str(DomainLowerLimitVec(i)) // "," &
                            // num2str(DomainUpperLimitVec(i)) // "). If you don't &
                            &know an appropriate value for startPointVec, drop it from the input list. " // &
                            methodName // " will automatically assign an appropriate value to it.\n\n"
            end if
        end do
    end subroutine checkForSanity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module SpecMCMC_StartPointVec_mod