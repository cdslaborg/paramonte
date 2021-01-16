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
!> This file implements the main body of the `ParaDRAM_ProposalAbstract_mod` and `ParaDISE_ProposalAbstract_mod` modules.
!>
!> \remark
!> This module requires preprocessing, prior to compilation.
!>
!> \author Amir Shahmoradi

    use Constants_mod, only: IK, RK
    use Err_mod, only: Err_type

    implicit none

#if defined PARADRAM
    character(*), parameter :: MODULE_NAME = "@ParaDRAM_Proposal_mod"
#elif defined PARADISE
    character(*), parameter :: MODULE_NAME = "@ParaDISE_Proposal_mod"
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type, abstract :: ProposalAbstract_type
        type(Err_type) :: Err
    contains
        procedure(getNew_proc)                  , pass      , deferred  :: getNew
#if defined PARADISE
        procedure(getLogProb_proc)              , nopass    , deferred  :: getLogProb
#endif
        procedure(doAdaptation_proc)            , pass      , deferred  :: doAdaptation
       !procedure(writeRestartFile_proc)        , pass      , deferred  :: writeRestartFile
#if defined CAF_ENABLED || defined MPI_ENABLED
        procedure(bcastAdaptation_proc)         , pass      , deferred  :: bcastAdaptation
#endif
    end type ProposalAbstract_type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if defined CAF_ENABLED || defined MPI_ENABLED
    abstract interface
        subroutine bcastAdaptation_proc(self)
            import :: ProposalAbstract_type
            class(ProposalAbstract_type), intent(inout) :: self
        end subroutine bcastAdaptation_proc
    end interface
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    abstract interface
    subroutine getNew_proc  ( self          &
                            , nd            &
                            , counterDRS    &
                            , StateOld      &
                            , StateNew      &
                            )
        use Constants_mod, only: IK, RK
        import :: ProposalAbstract_type
        class(ProposalAbstract_type), intent(inout) :: self
        integer(IK) , intent(in)    :: nd
        integer(IK) , intent(in)    :: counterDRS
        real(RK)    , intent(in)    :: StateOld(nd)
        real(RK)    , intent(out)   :: StateNew(nd)
    end subroutine getNew_proc
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    abstract interface
    function getLogProb_proc( nd                &
                            , counterDRS        &
                            , StateOld          &
                            , StateNew          &
                            ) result (logProb)
        use Constants_mod, only: IK, RK
        import :: ProposalAbstract_type
        integer(IK), intent(in) :: nd
        integer(IK), intent(in) :: counterDRS
        real(RK)   , intent(in) :: StateOld(nd)
        real(RK)   , intent(in) :: StateNew(nd)
        real(RK)                :: logProb
    end function getLogProb_proc
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    abstract interface
    subroutine doAdaptation_proc( self                      &
                                , nd                        &
                                , chainSize                 &
                                , Chain                     &
                                , ChainWeight               &
                                , isFreshRun                &
                                , samplerUpdateIsGreedy     &
                                , meanAccRateSinceStart     &
                                , samplerUpdateSucceeded    &
                                , adaptationMeasure         &
                                )
        use Constants_mod, only: IK, RK
        import :: ProposalAbstract_type
        class(ProposalAbstract_type), intent(inout) :: self
        integer(IK), intent(in)     :: nd
        integer(IK), intent(in)     :: chainSize
        real(RK)   , intent(in)     :: Chain(nd,chainSize)
        integer(IK), intent(in)     :: ChainWeight(chainSize)
        logical    , intent(in)     :: isFreshRun
        logical    , intent(in)     :: samplerUpdateIsGreedy
        real(RK)   , intent(inout)  :: meanAccRateSinceStart
        logical    , intent(out)    :: samplerUpdateSucceeded
        real(RK)   , intent(out)    :: adaptationMeasure
    end subroutine doAdaptation_proc
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !abstract interface
    !subroutine readRestartFileAscii_proc()
    !end subroutine readRestartFileAscii_proc
    !end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !abstract interface
    !subroutine writeRestartFile_proc(self, meanAccRateSinceStart)
    !    use Constants_mod, only: RK
    !    import :: ProposalAbstract_type
    !    class(ProposalAbstract_type), intent(in)    :: self
    !    real(RK), intent(in), optional              :: meanAccRateSinceStart
    !end subroutine writeRestartFile_proc
    !end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


