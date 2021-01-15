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
!> This module contains the classes and procedures for setting up the `proposalModel` attribute of samplers of class 
!> [ParaNest_type](@ref paranest_mod::paranest_type). For more information, see the description of this attribute in the body of the module.
!> \author Amir Shahmoradi

module SpecNest_ProposalModel_mod

    use Constants_mod, only: IK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecNest_ProposalModel_mod"

    integer(IK), parameter          :: MAX_LEN_PROPOSAL_MODEL = 63_IK

    character(:), allocatable       :: proposalModel ! namelist input

    type                            :: ProposalModel_type
        logical                     :: isRej        !< Rejection
        logical                     :: isRejEll     !< Ellipsoidal Rejection
        logical                     :: isRejBall    !< Spherical Rejection
        logical                     :: isRejCube    !< Cubical Rejection
        character(9)                :: rejection
        character(11)               :: ellipsoidal
        character(9)                :: spherical
        character(7)                :: cubical
        character(4)                :: ball
        character(:), allocatable   :: val
        character(:), allocatable   :: def
        character(:), allocatable   :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set, checkForSanity, nullifyNameListVar
    end type ProposalModel_type

    interface ProposalModel_type
        module procedure            :: construct
    end interface ProposalModel_type

    private :: construct, set, checkForSanity, nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function construct(methodName) result(self)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: construct
#endif
        use Decoration_mod, only: TAB
        use Constants_mod, only: NULL_SK, IK
        use String_mod, only: num2str
        implicit none
        character(*), intent(in)    :: methodName
        type(ProposalModel_type)   :: self
        self%isRej          = .false.
        self%isRejEll       = .false.
        self%isRejBall      = .false.
        self%isRejCube      = .false.
        self%rejection      = "rejection"
        self%ellipsoidal    = "ellipsoidal"
        self%spherical      = "spherical"
        self%ball           = "ball"
        self%def            = self%rejection//"-"//self%ellipsoidal
        self%null           = repeat(NULL_SK, MAX_LEN_PROPOSAL_MODEL)
        self%desc           = &
        "proposalModel is a string variable containing the name of the proposal distribution for the Nest sampler. &
        &The string value must be enclosed by either single or double quotation marks when provided as input. &
        &Options that are currently supported include:\n\n" // &
        "    proposalModel = 'rej' or 'rejection' or 'rej-ell' or 'rejection-ellipsoidal' \n\n" // &
        "            This is equivalent to the constrained rejection sampling &
                     &via bounding ellipsoids around the active point.\n\n&
        &    proposalModel = 'rej-ball'\n\n" // &
        "            The proposals will be drawn uniformly from within a ndim-dimensional ellipsoid whose covariance matrix &
                     &and scale are initialized by the user and optionally adaptively updated throughout the simulation.\n\n&
        &Note that are values are case-INsensitive and all hyphens (dashes, -) and white-space characters are ignored.&
        &The default value is '" // self%def // "'."
    end function construct

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine nullifyNameListVar(self)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(ProposalModel_type), intent(in) :: self
        proposalModel = self%null
    end subroutine nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure subroutine set(self,proposalModel)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: set
#endif
        use String_mod, only: getLowerCase
        implicit none
        class(ProposalModel_type), intent(inout)   :: self
        character(*), intent(in)                    :: proposalModel
        self%val = getLowerCase( trim(adjustl(proposalModel)) )
        if (self%val==self%null) self%val = self%def
        self%isRej      = index(self%val,"rej")
        self%isRejCube  = self%isRej .and. (index(self%val,"cube") > 0 .or. index(self%val,"cubical") > 0)
        self%isRejBall  = self%isRej .and. (index(self%val,"ball") > 0 .or. index(self%val,"spherical") > 0)
        self%isRejEll   = self%isRej .and. (index(self%val,"ell") > 0 .or. .not. (self%isRejCube .or. self%isRejBall))
    end subroutine set

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure subroutine checkForSanity(self,Err,methodName)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Err_mod, only: Err_type
        use String_mod, only: num2str
        implicit none
        class(ProposalModel_type), intent(in)   :: self
        character(*), intent(in)                :: methodName
        type(Err_type), intent(inout)           :: Err
        character(*), parameter                 :: PROCEDURE_NAME = "@checkForSanity()"
        if ( .not. (self%isRejEll .or. self%isRejBall .or. self%isRejCube) ) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &Invalid requested value for the proposalModel of " // methodName // ". The input requested &
                        &proposal model (" // self%val // ") is not supported. &
                        &The variable proposalModel cannot be set to anything other than the values described &
                        &in the description of the simulation specification proposalModel.\n\n"
        end if
    end subroutine checkForSanity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module SpecNest_ProposalModel_mod ! LCOV_EXCL_LINE