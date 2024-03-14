program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! All kinds are supported.
    use pm_kind, only: IK ! All kinds are supported.
    use pm_kind, only: CK ! All kinds are supported.
    use pm_kind, only: RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arrayCenter, only: setCentered

    implicit none

    character(:, SK), allocatable   ::   str_SK   ,   strCentered_SK    ! Can be any processor-supported kind.
    character(2, SK), allocatable   :: Array_SK(:), ArrayCentered_SK(:) ! Can be any processor-supported kind.
    logical(LK)     , allocatable   :: Array_LK(:), ArrayCentered_LK(:) ! Can be any processor-supported kind.
    integer(IK)     , allocatable   :: Array_IK(:), ArrayCentered_IK(:) ! Can be any processor-supported kind.
    complex(CK)     , allocatable   :: Array_CK(:), ArrayCentered_CK(:) ! Can be any processor-supported kind.
    real(RK)        , allocatable   :: Array_RK(:), ArrayCentered_RK(:) ! Can be any processor-supported kind.

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand and center an array to a new larger size.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call reset()

    str_SK = SK_"ABCDEF"
    allocate(Array_SK(3:8), source = [SK_"AA", SK_"BB", SK_"CC", SK_"DD", SK_"EE", SK_"FF"])
    allocate(Array_IK(3:8), source = [1_IK, 2_IK, 3_IK, 4_IK, 5_IK, 6_IK])
    allocate(Array_LK(3:8), source = [.true._LK, .true._LK, .true._LK, .true._LK, .true._LK, .true._LK])
    allocate(Array_CK(3:8), source = [(1._CK, -1._CK), (2._CK, -2._CK), (3._CK, -3._CK), (4._CK, -4._CK), (5._CK, -5._CK), (6._CK, -6._CK)])
    allocate(Array_RK(3:8), source = [1._RK, 2._RK, 3._RK, 4._RK, 5._RK, 6._RK])

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Center character scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call reset()
    call disp%show("str_SK")
    call disp%show( str_SK, deliml = SK_"""" )
    call disp%show("len(str_SK)")
    call disp%show( len(str_SK) )
    call disp%show("allocate(character(8,SK) :: strCentered_SK)")
                    allocate(character(8,SK) :: strCentered_SK)
    call disp%show("call setCentered(strCentered_SK, str_SK, fill = SK_'-')")
                    call setCentered(strCentered_SK, str_SK, fill = SK_'-')
    call disp%show("strCentered_SK")
    call disp%show( strCentered_SK , deliml = SK_"""" )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Center character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call reset()
    call disp%show("Array_SK")
    call disp%show( Array_SK, deliml = SK_"""" )
    call disp%show("allocate(character(2,SK) :: ArrayCentered_SK(8))")
                    allocate(character(2,SK) :: ArrayCentered_SK(8))
    call disp%show("call setCentered(ArrayCentered_SK, Array_SK, fill = SK_'%%')")
                    call setCentered(ArrayCentered_SK, Array_SK, fill = SK_'%%')
    call disp%show("ArrayCentered_SK")
    call disp%show( ArrayCentered_SK, deliml = SK_"""" )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Center logical array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call reset()
    call disp%show("Array_LK")
    call disp%show( Array_LK )
    call disp%show("allocate(ArrayCentered_LK(8))")
                    allocate(ArrayCentered_LK(8))
    call disp%show("call setCentered(ArrayCentered_LK, Array_LK, fill = .false._LK)")
                    call setCentered(ArrayCentered_LK, Array_LK, fill = .false._LK)
    call disp%show("ArrayCentered_LK")
    call disp%show( ArrayCentered_LK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Center integer array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call reset()
    call disp%show("Array_IK")
    call disp%show( Array_IK )
    call disp%show("allocate(ArrayCentered_IK(8))")
                    allocate(ArrayCentered_IK(8))
    call disp%show("call setCentered(ArrayCentered_IK, Array_IK, fill = -9999_IK)")
                    call setCentered(ArrayCentered_IK, Array_IK, fill = -9999_IK)
    call disp%show("ArrayCentered_IK")
    call disp%show( ArrayCentered_IK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Center complex array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call reset()
    call disp%show("Array_CK")
    call disp%show( Array_CK )
    call disp%show("allocate(ArrayCentered_CK(8))")
                    allocate(ArrayCentered_CK(8))
    call disp%show("call setCentered(ArrayCentered_CK, Array_CK, fill = (-999._CK,-999._CK))")
                    call setCentered(ArrayCentered_CK, Array_CK, fill = (-999._CK,-999._CK))
    call disp%show("ArrayCentered_CK")
    call disp%show( ArrayCentered_CK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Center real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call reset()
    call disp%show("Array_RK")
    call disp%show( Array_RK )
    call disp%show("allocate(ArrayCentered_RK(8))")
                    allocate(ArrayCentered_RK(8))
    call disp%show("call setCentered(ArrayCentered_RK, Array_RK, fill = 0._RK)")
                    call setCentered(ArrayCentered_RK, Array_RK, fill = 0._RK)
    call disp%show("ArrayCentered_RK")
    call disp%show( ArrayCentered_RK )


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Shrink an array to a new smaller size. Note that fill has no effect because the array shrinks.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    deallocate(Array_SK); allocate(Array_SK(3:8), source = [SK_"AA", SK_"BB", SK_"CC", SK_"DD", SK_"EE", SK_"FF"])
    deallocate(Array_IK); allocate(Array_IK(3:8), source = [1_IK, 2_IK, 3_IK, 4_IK, 5_IK, 6_IK])
    deallocate(Array_LK); allocate(Array_LK(3:8), source = [.true._LK, .true._LK, .true._LK, .true._LK, .true._LK, .true._LK])
    deallocate(Array_CK); allocate(Array_CK(3:8), source = [(1._CK, -1._CK), (2._CK, -2._CK), (3._CK, -3._CK), (4._CK, -4._CK), (5._CK, -5._CK), (6._CK, -6._CK)])
    deallocate(Array_RK); allocate(Array_RK(3:8), source = [1._RK, 2._RK, 3._RK, 4._RK, 5._RK, 6._RK])

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Center character scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call reset()
    str_SK = SK_"ABCDEF"
    call disp%show("str_SK")
    call disp%show( str_SK, deliml = SK_"""" )
    call disp%show("allocate(character(3,SK) :: strCentered_SK)")
                    allocate(character(3,SK) :: strCentered_SK)
    call disp%show("call setCentered(strCentered_SK, str_SK, fill = '~')")
                    call setCentered(strCentered_SK, str_SK, fill = '~')
    call disp%show("strCentered_SK")
    call disp%show( strCentered_SK , deliml = SK_"""" )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Center character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call reset()
    call disp%show("Array_SK")
    call disp%show( Array_SK, deliml = SK_"""" )
    call disp%show("allocate(character(2,SK) :: ArrayCentered_SK(3))")
                    allocate(character(2,SK) :: ArrayCentered_SK(3))
    call disp%show("call setCentered(ArrayCentered_SK, Array_SK, fill = '**')")
                    call setCentered(ArrayCentered_SK, Array_SK, fill = '**')
    call disp%show("ArrayCentered_SK")
    call disp%show( ArrayCentered_SK, deliml = SK_"""" )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Center logical array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call reset()
    call disp%show("Array_LK")
    call disp%show( Array_LK )
    call disp%show("allocate(ArrayCentered_LK(3))")
                    allocate(ArrayCentered_LK(3))
    call disp%show("call setCentered(ArrayCentered_LK, Array_LK, fill = .false._LK)")
                    call setCentered(ArrayCentered_LK, Array_LK, fill = .false._LK)
    call disp%show("ArrayCentered_LK")
    call disp%show( ArrayCentered_LK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Center integer array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call reset()
    call disp%show("Array_IK")
    call disp%show( Array_IK )
    call disp%show("allocate(ArrayCentered_IK(3))")
                    allocate(ArrayCentered_IK(3))
    call disp%show("call setCentered(ArrayCentered_IK, Array_IK, fill = 0_IK)")
                    call setCentered(ArrayCentered_IK, Array_IK, fill = 0_IK)
    call disp%show("ArrayCentered_IK")
    call disp%show( ArrayCentered_IK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Center complex array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call reset()
    call disp%show("Array_CK")
    call disp%show( Array_CK )
    call disp%show("allocate(ArrayCentered_CK(3))")
                    allocate(ArrayCentered_CK(3))
    call disp%show("call setCentered(ArrayCentered_CK, Array_CK, fill = (0._CK, 0._CK))")
                    call setCentered(ArrayCentered_CK, Array_CK, fill = (0._CK, 0._CK))
    call disp%show("ArrayCentered_CK")
    call disp%show( ArrayCentered_CK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Center real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call reset()
    call disp%show("Array_RK")
    call disp%show( Array_RK )
    call disp%show("allocate(ArrayCentered_RK(3))")
                    allocate(ArrayCentered_RK(3))
    call disp%show("call setCentered(ArrayCentered_RK, Array_RK, fill = +0._RK)")
                    call setCentered(ArrayCentered_RK, Array_RK, fill = +0._RK)
    call disp%show("ArrayCentered_RK")
    call disp%show( ArrayCentered_RK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Example of perfect array centering.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call reset()
    str_SK = "ABC"
    call disp%show("str_SK")
    call disp%show( str_SK, deliml = SK_"""" )
    call disp%show("allocate(character(7,SK) :: strCentered_SK)")
                    allocate(character(7,SK) :: strCentered_SK)
    call disp%show("call setCentered(strCentered_SK, str_SK, fill = SK_'-')")
                    call setCentered(strCentered_SK, str_SK, fill = SK_'-')
    call disp%show("strCentered_SK")
    call disp%show( strCentered_SK, deliml = SK_"""" )
    call disp%skip()

    call reset()
    str_SK = SK_"ABC"
    call disp%show("str_SK")
    call disp%show( str_SK, deliml = SK_"""" )
    call disp%show("allocate(character(1,SK) :: strCentered_SK)")
                    allocate(character(1,SK) :: strCentered_SK)
    call disp%show("call setCentered(strCentered_SK, str_SK, fill = SK_'-')")
                    call setCentered(strCentered_SK, str_SK, fill = SK_'-')
    call disp%show("strCentered_SK")
    call disp%show( strCentered_SK, deliml = SK_"""" )
    call disp%skip()

    call reset()
    str_SK = SK_"AB"
    call disp%show("str_SK")
    call disp%show( str_SK, deliml = SK_"""" )
    call disp%show("allocate(character(6,SK) :: strCentered_SK)")
                    allocate(character(6,SK) :: strCentered_SK)
    call disp%show("call setCentered(strCentered_SK, str_SK, fill = SK_'-')")
                    call setCentered(strCentered_SK, str_SK, fill = SK_'-')
    call disp%show("strCentered_SK")
    call disp%show( strCentered_SK, deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Example of imperfect array centering.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call reset()
    str_SK = SK_"ABC"
    call disp%show("str_SK")
    call disp%show( str_SK, deliml = SK_"""" )
    call disp%show("allocate(character(6,SK) :: strCentered_SK)")
                    allocate(character(6,SK) :: strCentered_SK)
    call disp%show("call setCentered(strCentered_SK, str_SK, fill = SK_'-')")
                    call setCentered(strCentered_SK, str_SK, fill = SK_'-')
    call disp%show("strCentered_SK")
    call disp%show( strCentered_SK , deliml = SK_"""" )
    call disp%skip()

    call reset()
    str_SK = SK_"ABC"
    call disp%show("str_SK")
    call disp%show( str_SK, deliml = SK_"""" )
    call disp%show("allocate(character(2,SK) :: strCentered_SK)")
                    allocate(character(2,SK) :: strCentered_SK)
    call disp%show("call setCentered(strCentered_SK, str_SK, fill = SK_'-')")
                    call setCentered(strCentered_SK, str_SK, fill = SK_'-')
    call disp%show("strCentered_SK")
    call disp%show( strCentered_SK , deliml = SK_"""" )
    call disp%skip()

    call reset()
    str_SK = SK_"AB"
    call disp%show("str_SK")
    call disp%show( str_SK, deliml = SK_"""" )
    call disp%show("allocate(character(3,SK) :: strCentered_SK)")
                    allocate(character(3,SK) :: strCentered_SK)
    call disp%show("call setCentered(strCentered_SK, str_SK, fill = SK_'-')")
                    call setCentered(strCentered_SK, str_SK, fill = SK_'-')
    call disp%show("strCentered_SK")
    call disp%show( strCentered_SK , deliml = SK_"""" )
    call disp%skip()

    call reset()
    str_SK = SK_"AB"
    call disp%show("str_SK")
    call disp%show( str_SK, deliml = SK_"""" )
    call disp%show("allocate(character(1,SK) :: strCentered_SK)")
                    allocate(character(1,SK) :: strCentered_SK)
    call disp%show("call setCentered(strCentered_SK, str_SK, fill = SK_'-')")
                    call setCentered(strCentered_SK, str_SK, fill = SK_'-')
    call disp%show("strCentered_SK")
    call disp%show( strCentered_SK , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Centering array contents within a given size, padded by left and right margins.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call reset()
    str_SK = SK_"In God We Trust"
    call disp%show("str_SK")
    call disp%show( str_SK , deliml = SK_"""" )
    call disp%show("allocate(character(45+4+4,SK) :: strCentered_SK)")
                    allocate(character(45+4+4,SK) :: strCentered_SK)
    call disp%show("call setCentered(strCentered_SK, str_SK, fill = SK_' ', lmsize = 4_IK, rmsize = 4_IK, lmfill = SK_'*', rmfill = SK_'*')")
                    call setCentered(strCentered_SK, str_SK, fill = SK_' ', lmsize = 4_IK, rmsize = 4_IK, lmfill = SK_'*', rmfill = SK_'*')
    call disp%show("strCentered_SK")
    call disp%show( strCentered_SK , deliml = SK_"""" )
    call disp%skip()

    call reset()
    str_SK = SK_"ABCDEF"
    call disp%show("str_SK")
    call disp%show( str_SK, deliml = SK_"""" )
    call disp%show("allocate(character(8+3+5,SK) :: strCentered_SK)")
                    allocate(character(8+3+5,SK) :: strCentered_SK)
    call disp%show("call setCentered(strCentered_SK, str_SK, lmsize = 3_IK, rmsize = 5_IK, fill = SK_'-', lmfill = SK_'*', rmfill = SK_'+')")
                    call setCentered(strCentered_SK, str_SK, lmsize = 3_IK, rmsize = 5_IK, fill = SK_'-', lmfill = SK_'*', rmfill = SK_'+')
    call disp%show("strCentered_SK")
    call disp%show( strCentered_SK , deliml = SK_"""" )
    call disp%skip()

    call reset()
    str_SK = SK_"ABCDEF"
    call disp%show("str_SK")
    call disp%show( str_SK, deliml = SK_"""" )
    call disp%show("allocate(character(2+3+5,SK) :: strCentered_SK)")
                    allocate(character(2+3+5,SK) :: strCentered_SK)
    call disp%show("call setCentered(strCentered_SK, str_SK, lmsize = 3_IK, rmsize = 5_IK, fill = SK_'-', lmfill = SK_'*', rmfill = SK_'+')")
                    call setCentered(strCentered_SK, str_SK, lmsize = 3_IK, rmsize = 5_IK, fill = SK_'-', lmfill = SK_'*', rmfill = SK_'+')
    call disp%show("strCentered_SK")
    call disp%show( strCentered_SK , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Call to center() preserves the array lower bound.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call reset()
    deallocate(Array_IK); allocate(Array_IK(-6:-1), source = [1_IK, 2_IK, 3_IK, 4_IK, 5_IK, 6_IK])

    call disp%show("Array_IK")
    call disp%show( Array_IK )
    call disp%show("allocate(ArrayCentered_IK(10 + 3 + 5))")
                    allocate(ArrayCentered_IK(10 + 3 + 5))
    call disp%show("call setCentered(ArrayCentered_IK, Array_IK, lmsize = 3_IK, rmsize = 5_IK, fill = +0_IK, lmfill = -huge(0_IK), rmfill = +huge(0_IK))")
                    call setCentered(ArrayCentered_IK, Array_IK, lmsize = 3_IK, rmsize = 5_IK, fill = +0_IK, lmfill = -huge(0_IK), rmfill = +huge(0_IK))
    call disp%show("ArrayCentered_IK")
    call disp%show( ArrayCentered_IK )

contains

    subroutine reset()
        if (allocated(  strCentered_SK)) deallocate(  strCentered_SK)
        if (allocated(ArrayCentered_SK)) deallocate(ArrayCentered_SK)
        if (allocated(ArrayCentered_IK)) deallocate(ArrayCentered_IK)
        if (allocated(ArrayCentered_RK)) deallocate(ArrayCentered_RK)
        if (allocated(ArrayCentered_CK)) deallocate(ArrayCentered_CK)
        if (allocated(ArrayCentered_LK)) deallocate(ArrayCentered_LK)
    end subroutine

end program example