!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief 
!>  This module contains tests of the module [pm_polation](@ref pm_polation).
!>
!>  \final
!>
!>  \author 
!>  Amir Shahmoradi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_polation

    use pm_polation
    use pm_err, only: err_type
    use pm_test, only: test_type, LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

    real(RK), parameter :: SortedX(*) = [ +3.178053830347950_RK &
                                        , +3.218875824868200_RK &
                                        , +3.295836866004330_RK &
                                        , +3.433987204485150_RK &
                                        , +3.663561646129650_RK &
                                        , +4.007333185232470_RK &
                                        , +4.465908118654580_RK &
                                        , +5.017279836814920_RK &
                                        , +5.631211781821370_RK &
                                        , +6.282266746896010_RK &
                                        , +6.953684210870540_RK &
                                        , +7.635786861395590_RK &
                                        , +8.323365694436080_RK &
                                        , +9.013717030471370_RK &
                                        , +9.705463352014880_RK ]
    real(RK), parameter :: SortedY(*) = [ +1.426588681046810_RK &
                                        , +1.168902956625920_RK &
                                        , +0.882431985413412_RK &
                                        , +0.587833620221075_RK &
                                        , +0.229810034869442_RK &
                                        , -0.174864514871957_RK &
                                        , -0.605947339379271_RK &
                                        , -1.074884950466410_RK &
                                        , -1.510025316819630_RK &
                                        , -1.948747191172130_RK &
                                        , -2.382173017436780_RK &
                                        , -2.805304379208950_RK &
                                        , -3.199663563773090_RK &
                                        , -3.602191628213170_RK &
                                        , -3.989064786493580_RK ]

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)
        call test%run(test_InterpLinear_type_1, SK_"test_InterpLinear_type_1")
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_InterpLinear_type_1() result(assertion)

        use pm_kind, only: RK, IK, LK

        implicit none

        logical(LK)             :: assertion
        real(RK), parameter     :: TOLERANCE = 1.e-10_RK
        real(RK), parameter     :: QueryX(*) =      [ +3.17805383034795_RK &
                                                    , +7.61283103040736_RK &
                                                    , +8.30003171177958_RK &
                                                    , +8.70350676947973_RK &
                                                    , +8.99019232964177_RK &
                                                    , +9.21273749657591_RK &
                                                    , +9.39465993143281_RK &
                                                    , +9.54852542660107_RK &
                                                    , +9.68184287734565_RK &
                                                    , +9.79945948211208_RK ]
        real(RK), parameter     :: QueryY_ref(*) =  [ +1.42658868104681_RK &
                                                    , -2.79106410024333_RK &
                                                    , -3.18628041418956_RK &
                                                    , -3.42131512460047_RK &
                                                    , -3.58847491355328_RK &
                                                    , -3.71349785982781_RK &
                                                    , -3.81524167073542_RK &
                                                    , -3.90129406915411_RK &
                                                    , -3.97585455703407_RK &
                                                    , -4.04163402840252_RK ]
        real(RK)                :: QueryY(size(QueryX))
        type(InterpLinear_type) :: Interp

        Interp = InterpLinear_type(SortedX, SortedY)
        QueryY = Interp%predict(QueryX)
        assertion = all(abs(QueryY - QueryY_ref) <= TOLERANCE)

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "QueryY_ref ", QueryY_ref
            write(test%disp%unit,"(*(g0,:,', '))") "QueryY     ", QueryY
            write(test%disp%unit,"(*(g0,:,', '))") "diff       ", QueryY_ref - QueryY
            write(test%disp%unit,"(*(g0,:,', '))")
            ! LCOV_EXCL_STOP
        end if

    end function test_InterpLinear_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_polation