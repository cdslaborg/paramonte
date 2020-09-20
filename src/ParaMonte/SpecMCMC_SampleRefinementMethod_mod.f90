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

module SpecMCMC_SampleRefinementMethod_mod

    use Constants_mod, only: IK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecMCMC_SampleRefinementMethod_mod"

    character(*), parameter         :: BATCH_MEANS_METHOD_NAME = "BatchMeans"
    character(*), parameter         :: CUTOFF_AUTOCORR_METHOD_NAME = "CutOffAutoCorr"
    character(*), parameter         :: MAX_CUMSUM_AUTOCORR_METHOD_NAME = "MaxCumSumAutoCorr"
    integer(IK) , parameter         :: LEN_BATCH_MEANS_METHOD_NAME = len(BATCH_MEANS_METHOD_NAME)
    integer(IK) , parameter         :: LEN_CUTOFF_AUTOCORR_METHOD_NAME = len(CUTOFF_AUTOCORR_METHOD_NAME)
    integer(IK) , parameter         :: LEN_MAX_CUMSUM_AUTOCORR_METHOD_NAME = len(MAX_CUMSUM_AUTOCORR_METHOD_NAME)
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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        SampleRefinementMethodObj%def = BATCH_MEANS_METHOD_NAME
        SampleRefinementMethodObj%null = repeat(NULL_SK, MAX_LEN_SAMPLE_REFINEMENT_METHOD)

        SampleRefinementMethodObj%desc = &
        "sampleRefinementMethod is a string variable that represents the method of computing the Integrated Autocorrelation Time &
        &(IAC) to be used in "// methodName //" for refining the final output MCMC chain and sample. &
        &The string value must be enclosed by either single or double quotation marks when provided as input. &
        &Options that are currently supported include:\n\n&
        &    sampleRefinementMethod = '" // BATCH_MEANS_METHOD_NAME // "'\n\n&
        &            This method of computing the Integrated Autocorrelation Time is based on the approach described in &
                    &SCHMEISER, B., 1982, Batch size effects in the analysis of simulation output, Oper. Res. 30 556-568. The &
                    &batch sizes in the BatchMeans method are chosen to be int(N^(2/3)) where N is the length of the MCMC chain. &
                    &As long as the batch size is larger than the IAC of the chain and there are significantly more than 10 &
                    &batches, the BatchMeans method will provide reliable estimates of the IAC. &
                    &Note that the refinement strategy involves two separate phases of sample decorrelation. At the first stage, &
                    &the Markov chain is decorrelated recursively (for as long as needed) based on the IAC of its compact format, &
                    &where only the the uniquely-visited states are kept in the (compact) chain. Once the Markov chain is refined &
                    &such that its compact format is fully decorrelated, the second phase of the decorrelation begins during which &
                    &the Markov chain is decorrelated based on the IAC of the chain in its verbose (Markov) format. This process &
                    &is repeated recursively for as long as there is any residual autocorrelation in the refined sample.\n\n&
        &    sampleRefinementMethod = '" // BATCH_MEANS_METHOD_NAME // "-compact'\n\n&
        &            This is the same as the first case in the above, except that only the first phase of the sample refinement &
                    &described in the above will be performed, that is, the (verbose) Markov chain is refined only based on the &
                    &IAC computed from the compact format of the Markov chain. This will lead to a larger final refined sample. &
                    &However, the final sample will likely not be fully decorrelated.\n\n&
        &    sampleRefinementMethod = '" // BATCH_MEANS_METHOD_NAME // "-verbose'\n\n&
        &            This is the same as the first case in the above, except that only the second phase of the sample refinement &
                    &described in the above will be performed, that is, the (verbose) Markov chain is refined only based on the &
                    &IAC computed from the verbose format of the Markov chain. While the resulting refined sample will be fully &
                    &decorrelated, the size of the refined sample may be smaller than the default choice in the first case in the &
                    &above.\n\n&
        &Note that in order to obtain i.i.d. samples from a multidimensional chain, "//methodName//" will use the maximum of &
        &IAC among all dimensions of the chain to refine the chain. Also, note that the value specified for sampleRefinementCount &
        &is used only when the variable sampleSize < 0, otherwise, it will be ignored. &
        &The default value is sampleRefinementMethod = '" // SampleRefinementMethodObj%def // "'. &
        &Note that the input values are case-insensitive and white-space characters are ignored."
    end function constructSampleRefinementMethod

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine nullifyNameListVar(SampleRefinementMethodObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(SampleRefinementMethod_type), intent(in) :: SampleRefinementMethodObj
        sampleRefinementMethod = SampleRefinementMethodObj%null
    end subroutine nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
       !if ( sampleRefinementMethodLowerCase == getLowerCase(MAX_CUMSUM_AUTOCORR_METHOD_NAME) ) then
       !    SampleRefinementMethodObj%isMaxCumSumAutoCorr = .true.
       !elseif ( sampleRefinementMethodLowerCase(1:LEN_BATCH_MEANS_METHOD_NAME) == getLowerCase(BATCH_MEANS_METHOD_NAME) ) then
       !    if ( sampleRefinementMethodLowerCase(LEN_BATCH_MEANS_METHOD_NAME+1:lenSampleRefinementMethod) == "compact" ) then
       !        SampleRefinementMethodObj%isViaCompactChain = .true.
       !        SampleRefinementMethodObj%isBatchMeans = .true.
       !    elseif( len_trim(sampleRefinementMethodLowerCase(LEN_BATCH_MEANS_METHOD_NAME+1:lenSampleRefinementMethod)) == 0 ) then
       !        SampleRefinementMethodObj%isBatchMeans = .true.
       !    end if
       !end if

    end subroutine setSampleRefinementMethod

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine checkForSanity(SampleRefinementMethodObj,Err,methodName)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use String_mod, only: getLowerCase, replaceStr
        use String_mod, only: num2str
        use Err_mod, only: Err_type
        implicit none
        class(SampleRefinementMethod_type), intent(in)  :: SampleRefinementMethodObj
        character(*), intent(in)                        :: methodName
        type(Err_type), intent(inout)                   :: Err
        character(*), parameter                         :: PROCEDURE_NAME = "@checkForSanity()"
        character(:), allocatable                       :: sampleRefinementMethodLowerCase
        sampleRefinementMethodLowerCase = getLowerCase(SampleRefinementMethodObj%val)
        if  (index(sampleRefinementMethodLowerCase,getLowerCase(replaceStr(BATCH_MEANS_METHOD_NAME," ","")))==0 &
            .and. &
            (index(sampleRefinementMethodLowerCase,getLowerCase(CUTOFF_AUTOCORR_METHOD_NAME))==0 .and. index(sampleRefinementMethodLowerCase,"cutoff")==0) &
            .and. &
            (index(sampleRefinementMethodLowerCase,getLowerCase(MAX_CUMSUM_AUTOCORR_METHOD_NAME))==0 .and. index(sampleRefinementMethodLowerCase,"cumsum")==0) &
            ) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &The input requested method for the computation of the Integrated Autocorrelation Time (" // &
                        SampleRefinementMethodObj%val // ") assigned to the variable sampleRefinementMethod cannot be anything other than " // &
                        BATCH_MEANS_METHOD_NAME  // &
                        ". " // &
                        ! " or " // MAX_CUMSUM_AUTOCORR_METHOD_NAME // ". &
                        "If you are not sure of the appropriate value for SampleRefinementMethod, drop it from the input list. " // &
                        methodName // " will automatically assign an appropriate value to it.\n\n"
        end if
    end subroutine checkForSanity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module SpecMCMC_SampleRefinementMethod_mod