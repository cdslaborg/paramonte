!**********************************************************************************************************************************
!**********************************************************************************************************************************
!
!  ParaMonte: plain powerful parallel Monte Carlo library.
!
!  Copyright (C) 2012-present, The Computational Data Science Lab
!
!  This file is part of ParaMonte library. 
!
!  ParaMonte is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, version 3 of the License.
!
!  ParaMonte is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with ParaMonte.  If not, see <https://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!**********************************************************************************************************************************

module SpecMCMC_SampleRefinementMethod_mod

    use Constants_mod, only: IK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecMCMC_SampleRefinementMethod_mod"

    character(*), parameter         :: BATCH_MEANS = "BatchMeans"
    character(*), parameter         :: MAX_CUMSUM_AUTOCORR = "MaxCumSumAutoCorr"
    integer(IK) , parameter         :: LEN_MAX_CUMSUM_AUTOCORR = len(MAX_CUMSUM_AUTOCORR)
    integer(IK) , parameter         :: LEN_BATCH_MEANS = len(BATCH_MEANS)
    integer(IK) , parameter         :: MAX_LEN_SAMPLE_REFINEMENT_METHOD = 63

    character(MAX_LEN_SAMPLE_REFINEMENT_METHOD) :: sampleRefinementMethod ! namelist input

    type                            :: SampleRefinementMethod_type
       !logical                     :: isMaxCumSumAutoCorr
       !logical                     :: isViaCompactChain
       !logical                     :: isBatchMeans
        character(:), allocatable   :: def
        character(:), allocatable   :: val
        character(:), allocatable   :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setSampleRefinementMethod, checkForSanity, nullifyNameListVar
    end type SampleRefinementMethod_type

    interface SampleRefinementMethod_type
        module procedure            :: constructSampleRefinementMethod
    end interface SampleRefinementMethod_type

    private :: constructSampleRefinementMethod, setSampleRefinementMethod, checkForSanity, nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function constructSampleRefinementMethod(methodName) result(SampleRefinementMethodObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructSampleRefinementMethod
#endif
        use Constants_mod, only: NULL_SK, IK
        use Decoration_mod, only: TAB
        use String_mod, only: num2str
        implicit none
        character(*), intent(in)            :: methodName
        type(SampleRefinementMethod_type)   :: SampleRefinementMethodObj

       !SampleRefinementMethodObj%isMaxCumSumAutoCorr = .false.
       !SampleRefinementMethodObj%isViaCompactChain = .false.
       !SampleRefinementMethodObj%isBatchMeans = .false.
        SampleRefinementMethodObj%def = BATCH_MEANS
        SampleRefinementMethodObj%null = repeat(NULL_SK, MAX_LEN_SAMPLE_REFINEMENT_METHOD)

        SampleRefinementMethodObj%desc = &
        "sampleRefinementMethod is a string variable that represents the method of computing the Integrated Autocorrelation Time &
        &(IAC) to be used in "// methodName //" for refining the final output MCMC chain and sample. &
        &The string value must be enclosed by either single or double quotation marks when provided as input. &
        &Options that are currently &
        &supported include:\n\n&
        &    sampleRefinementMethod = '" // BATCH_MEANS // "'\n\n&
        &            This method of computing the Integrated Autocorrelation Time is based on the approach described in &
                    &SCHMEISER, B., 1982, Batch size effects in the analysis of simulation output, Oper. Res. 30 556-568. The &
                    &batch sizes in the BatchMeans method are chosen to be int(N^(2/3)) where N is the length of the MCMC chain. &
                    &As long as the batch size is larger than the IAC of the chain and there are significantly more than 10 &
                    &batches, the BatchMeans method will provide reliable estimates of the IAC.\n\n&
        &Note that in order to obtain i.i.d. samples from a multidimensional chain, "//methodName//" will use the maximum of &
        &IAC among all dimensions of the chain to refine the chain. Also, note that the value specified for sampleRefinementCount &
        &is used only when the variable sampleSize < 0, otherwise, it will be ignored. &
        &The default value is sampleRefinementMethod = '" // SampleRefinementMethodObj%def // "'. &
        &Note that the input values are case-insensitive and white-space characters are ignored."
    end function constructSampleRefinementMethod

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine nullifyNameListVar(SampleRefinementMethodObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(SampleRefinementMethod_type), intent(in) :: SampleRefinementMethodObj
        sampleRefinementMethod = SampleRefinementMethodObj%null
    end subroutine nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setSampleRefinementMethod(SampleRefinementMethodObj,sampleRefinementMethod)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setSampleRefinementMethod
#endif
        use String_mod, only: replaceStr !, getLowerCase
        implicit none
        class(SampleRefinementMethod_type), intent(inout)   :: SampleRefinementMethodObj
        character(*), intent(in)                            :: sampleRefinementMethod
       !character(:), allocatable                           :: sampleRefinementMethodLowerCase
       !integer                                             :: lenSampleRefinementMethod

        SampleRefinementMethodObj%val = trim(adjustl(replaceStr(sampleRefinementMethod," ", "")))
        if ( SampleRefinementMethodObj%val==trim(adjustl(SampleRefinementMethodObj%null)) ) then
            SampleRefinementMethodObj%val = SampleRefinementMethodObj%def
        end if

       !lenSampleRefinementMethod = len(SampleRefinementMethodObj%val)
       !sampleRefinementMethodLowerCase = getLowerCase(SampleRefinementMethodObj%val)
       !if ( sampleRefinementMethodLowerCase == getLowerCase(MAX_CUMSUM_AUTOCORR) ) then
       !    SampleRefinementMethodObj%isMaxCumSumAutoCorr = .true.
       !elseif ( sampleRefinementMethodLowerCase(1:LEN_BATCH_MEANS) == getLowerCase(BATCH_MEANS) ) then
       !    if ( sampleRefinementMethodLowerCase(LEN_BATCH_MEANS+1:lenSampleRefinementMethod) == "compact" ) then
       !        SampleRefinementMethodObj%isViaCompactChain = .true.
       !        SampleRefinementMethodObj%isBatchMeans = .true.
       !    elseif( len_trim(sampleRefinementMethodLowerCase(LEN_BATCH_MEANS+1:lenSampleRefinementMethod)) == 0 ) then
       !        SampleRefinementMethodObj%isBatchMeans = .true.
       !    end if
       !end if

    end subroutine setSampleRefinementMethod

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine checkForSanity(SampleRefinementMethodObj,Err,methodName)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use String_mod, only: getLowerCase
        use String_mod, only: num2str
        use Err_mod, only: Err_type
        implicit none
        class(SampleRefinementMethod_type), intent(in)  :: SampleRefinementMethodObj
        character(*), intent(in)                        :: methodName
        type(Err_type), intent(inout)                   :: Err
        character(*), parameter                         :: PROCEDURE_NAME = "@checkForSanity()"
        character(:), allocatable                       :: sampleRefinementMethodLowerCase
        sampleRefinementMethodLowerCase = getLowerCase(SampleRefinementMethodObj%val)
        if (index(sampleRefinementMethodLowerCase,getLowerCase(BATCH_MEANS))==0 .and. &
            index(sampleRefinementMethodLowerCase,getLowerCase(MAX_CUMSUM_AUTOCORR))==0) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &The input requested method for the computation of the Integrated Autocorrelation Time (" // &
                        SampleRefinementMethodObj%val // ") assigned to the variable sampleRefinementMethod cannot be anything other than " // &
                        BATCH_MEANS  // &
                        ". " // &
                        ! " or " // MAX_CUMSUM_AUTOCORR // ". &
                        "If you are not sure of the appropriate value for SampleRefinementMethod, drop it from the input list. " // &
                        methodName // " will automatically assign an appropriate value to it.\n\n"
        end if
    end subroutine checkForSanity

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SpecMCMC_SampleRefinementMethod_mod