!***********************************************************************************************************************************
!***********************************************************************************************************************************
!
!   ParaMonte: plain powerful parallel Monte Carlo library.
!
!   Copyright (C) 2012-present, The Computational Data Science Lab
!
!   This file is part of the ParaMonte library.
!
!   ParaMonte is free software: you can redistribute it and/or modify it
!   under the terms of the GNU Lesser General Public License as published
!   by the Free Software Foundation, version 3 of the License.
!
!   ParaMonte is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with the ParaMonte library. If not, see,
!
!       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
!
!   ACKNOWLEDGMENT
!
!   As per the ParaMonte library license agreement terms,
!   if you use any parts of this library for any purposes,
!   we ask you to acknowledge the use of the ParaMonte library
!   in your work (education/research/industry/development/...)
!   by citing the ParaMonte library as described on this page:
!
!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!
!***********************************************************************************************************************************
!***********************************************************************************************************************************

module SpecBase_SampleSize_mod

    use Constants_mod, only: IK
    implicit none

    integer(IK)                     :: sampleSize ! namelist input

    character(*), parameter         :: MODULE_NAME = "@SpecBase_SampleSize_mod"

    type                            :: SampleSize_type
        integer(IK)                 :: val
        integer(IK)                 :: def
        integer(IK)                 :: null
        character(:), allocatable   :: str
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setSampleSize, checkForSanity, nullifyNameListVar
    end type SampleSize_type

    interface SampleSize_type
        module procedure            :: constructSampleSize
    end interface SampleSize_type

    private :: constructSampleSize, setSampleSize, nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function constructSampleSize(methodName) result(SampleSizeObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructSampleSize
#endif
        use Constants_mod, only: IK, NULL_IK
        use String_mod, only: num2str
        implicit none
        character(*), intent(in)    :: methodName
        type(SampleSize_type)       :: SampleSizeObj
        SampleSizeObj%def  = -1
        SampleSizeObj%null = NULL_IK
        SampleSizeObj%desc = &
        "The variable sampleSize is an integer that dictates the number of (hopefully, independent and identically distributed &
        &[i.i.d.]) samples to be drawn from the user-provided objective function. Three ranges of values are possible:\n\n&
        &    sampleSize < 0:\n\n&
        &            Then, the absolute value of sampleSize dictates the sample size in units of the effective sample size. &
                     &The effective sample is by definition i.i.d., and free from duplicates. The effective sample size &
                     &is determined by " // methodName // " automatically toward the end of the simulation.\n&
        &            For example:\n\n&
        &                    sampleSize = -1 yields the effective i.i.d. sample drawn from the objective function.\n\n&
        &                    sampleSize = -2 yields a (potentially non-i.i.d.) sample twice as big as the effective &
                             &sample.\n\n&
        &    sampleSize > 0:\n\n&
        &            Then, the sample size is assumed to be in units of the number of points to be sampled. &
                     &If sampleSize turns out to be less than effectiveSampleSize, the resulting sample will be i.i.d.. &
                     &If sampleSize turns out to be larger than effectiveSampleSize, the resulting sample will be &
                     &potentially non-i.i.d.. The larger the difference, the more non-i.i.d. the resulting sample will be.\n&
        &            For example:\n\n&
        &                    sampleSize = 1000 yields a 1000-points sample from the objective function.\n\n&
        &    sampleSize = 0:\n\n&
        &            in which case, no sample file will be generated.\n\n&
        &Default value is sampleSize = "// num2str(SampleSizeObj%def) //"."
    end function constructSampleSize

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine nullifyNameListVar(SampleSizeObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(SampleSize_type), intent(in)  :: SampleSizeObj
        sampleSize = SampleSizeObj%null
    end subroutine nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setSampleSize(SampleSizeObj,sampleSize)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setSampleSize
#endif
        use String_mod, only: num2str
        use Constants_mod, only: IK
        implicit none
        class(SampleSize_type), intent(inout)   :: SampleSizeObj
        integer(IK), intent(in)                 :: sampleSize
        SampleSizeObj%val = sampleSize
        if (SampleSizeObj%val==SampleSizeObj%null) then
            SampleSizeObj%val = SampleSizeObj%def
        end if
        SampleSizeObj%str = num2str(SampleSizeObj%val)
    end subroutine setSampleSize

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine checkForSanity(SampleSizeObj,Err,methodName)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Err_mod, only: Err_type
        use String_mod, only: num2str
        implicit none
        class(SampleSize_type), intent(in) :: SampleSizeObj
        character(*), intent(in)           :: methodName
        type(Err_type), intent(inout)      :: Err
        character(*), parameter            :: PROCEDURE_NAME = "@checkForSanity()"
        if ( SampleSizeObj%val<1 ) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &The input value for variable sampleSize must be a positive integer. If you are not sure about the &
                        &appropriate value for this variable, simply drop it from the input. " // methodName // &
                        " will automatically assign an appropriate value to it.\n\n"
        end if
    end subroutine checkForSanity

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SpecBase_SampleSize_mod