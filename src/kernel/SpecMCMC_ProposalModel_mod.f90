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

module SpecMCMC_ProposalModel_mod

    use Constants_mod, only: IK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecMCMC_ProposalModel_mod"

    integer(IK), parameter          :: MAX_LEN_PROPOSAL_MODEL = 63_IK

    character(:), allocatable       :: proposalModel ! namelist input

    type                            :: ProposalModel_type
        logical                     :: isUniform
        logical                     :: isNormal
        character(7)                :: uniform
        character(6)                :: normal
        character(:), allocatable   :: val
        character(:), allocatable   :: def
        character(:), allocatable   :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setProposalModel, checkForSanity, nullifyNameListVar
    end type ProposalModel_type

    interface ProposalModel_type
        module procedure            :: constructProposalModel
    end interface ProposalModel_type

    private :: constructProposalModel, setProposalModel, checkForSanity, nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function constructProposalModel() result(ProposalModelObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructProposalModel
#endif
        use Decoration_mod, only: TAB
        use Constants_mod, only: NULL_SK, IK
        use String_mod, only: num2str
        implicit none
        type(ProposalModel_type)    :: ProposalModelObj
        ProposalModelObj%isUniform  = .false.
        ProposalModelObj%isNormal   = .false.
        ProposalModelObj%uniform    = "uniform"
        ProposalModelObj%normal     = "normal"
        ProposalModelObj%def        = ProposalModelObj%normal
        ProposalModelObj%null       = repeat(NULL_SK, MAX_LEN_PROPOSAL_MODEL)
        ProposalModelObj%desc       = &
        "proposalModel is a string variable containing the name of the proposal distribution for the MCMC sampler. &
        &The string value must be enclosed by either single or double quotation marks when provided as input. &
        &One option is currently supported:\n\n" // &
        "    proposalModel = '" // ProposalModelObj%normal // "'\n\n" // &
        "            This is equivalent to the multivariate normal distribution&
                     &, which is the most widely-used proposal model along with MCMC samplers.\n\n&
        &    proposalModel = '" // ProposalModelObj%uniform // "'\n\n" // &
        "            The proposals will be drawn uniformly from within a ndim-dimensional ellipsoid whose covariance matrix &
                     &and scale are initialized by the user and optionally adaptively updated throughout the simulation.\n\n&
        &The default value is '" // ProposalModelObj%def // "'."
    end function constructProposalModel

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine nullifyNameListVar(ProposalModelObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(ProposalModel_type), intent(in) :: ProposalModelObj
        proposalModel = ProposalModelObj%null
    end subroutine nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setProposalModel(ProposalModelObj,proposalModel)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setProposalModel
#endif
        use String_mod, only: getLowerCase
        implicit none
        class(ProposalModel_type), intent(inout)    :: ProposalModelObj
        character(*), intent(in)                    :: proposalModel
        ProposalModelObj%val = getLowerCase( trim(adjustl(proposalModel)) )
        if (ProposalModelObj%val==trim(adjustl(ProposalModelObj%null))) ProposalModelObj%val = trim(adjustl(ProposalModelObj%def))
        if (ProposalModelObj%val==ProposalModelObj%normal) ProposalModelObj%isNormal = .true.
        if (ProposalModelObj%val==ProposalModelObj%uniform) ProposalModelObj%isUniform = .true.
    end subroutine setProposalModel

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine checkForSanity(ProposalModelObj,Err,methodName)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Err_mod, only: Err_type
        use String_mod, only: num2str
        implicit none
        class(ProposalModel_type), intent(in)   :: ProposalModelObj
        character(*), intent(in)                :: methodName
        type(Err_type), intent(inout)           :: Err
        character(*), parameter                 :: PROCEDURE_NAME = "@checkForSanity()"
        if ( .not. (ProposalModelObj%isNormal .or. ProposalModelObj%isUniform) ) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &Invalid requested value for the proposalModel of " // methodName // ". The input requested &
                        &proposal model (" // ProposalModelObj%val // ") is not supported. &
                        &The variable proposalModel cannot be set to anything other than '" // &
                        ProposalModelObj%normal     // "', or '" // &
                        ProposalModelObj%uniform    // "'.\n\n"
        end if
    end subroutine checkForSanity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module SpecMCMC_ProposalModel_mod