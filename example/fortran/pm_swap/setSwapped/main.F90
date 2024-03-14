program example

    use pm_kind, only: SK, IK, LK, CK, RK ! all intrinsic types and kinds are supported.
    use pm_io, only: display_type
    use pm_swap, only: setSwapped

    implicit none

    type(display_type) :: disp

    character(:, SK), allocatable   :: a_SK, b_SK, Pair_SK(:)
    integer(IK)                     :: a_IK, b_IK, Pair_IK(2)
    logical(LK)                     :: a_LK, b_LK, Pair_LK(2)
    complex(CK)                     :: a_CK, b_CK, Pair_CK(2)
    real(RK)                        :: a_RK, b_RK, Pair_RK(2)

    a_SK = "Hell"           ; b_SK = "Heaven"
    a_IK = 10_IK            ; b_IK = 5_IK
    a_LK = .true._LK        ; b_LK = .false._LK
    a_CK = (+10._CK, -5._CK); b_CK = (+5._CK, +15._CK)
    a_RK = 10._RK           ; b_RK = 5._RK

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("! Swap scalars.")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip()

    block
        character(:), allocatable :: a, b
        call disp%skip()
        call disp%show("a = 'Heaven'")
                        a = 'Heaven'
        call disp%show("b = 'Hell  '")
                        b = 'Hell  '
        call disp%show("call setSwapped(a, b)")
                        call setSwapped(a, b)
        call disp%show("a")
        call disp%show( a , deliml = SK_"""" )
        call disp%show("b")
        call disp%show( b , deliml = SK_"""" )
        call disp%skip()
    end block

    block
        character(:), allocatable :: a, b
        call disp%skip()
        call disp%show("a = 'p r m n e'")
                        a = 'p r m n e'
        call disp%show("b = ' a a o t '")
                        b = ' a a o t '
        call disp%show("call setSwapped(a, b, inca = 2_IK, incb = 2_IK)")
                        call setSwapped(a, b, inca = 2_IK, incb = 2_IK)
        call disp%show("a")
        call disp%show( a , deliml = SK_"""" )
        call disp%show("b")
        call disp%show( b , deliml = SK_"""" )
        call disp%skip()
    end block

    block
        integer, allocatable :: a, b
        call disp%skip()
        call disp%show("a = 1")
                        a = 1
        call disp%show("b = 2")
                        b = 2
        call disp%show("call setSwapped(a, b)")
                        call setSwapped(a, b)
        call disp%show("a")
        call disp%show( a , deliml = SK_"""" )
        call disp%show("b")
        call disp%show( b , deliml = SK_"""" )
        call disp%skip()
    end block

    block
        logical, allocatable :: a, b
        call disp%skip()
        call disp%show("a = .false.")
                        a = .false.
        call disp%show("b = .true.")
                        b = .true.
        call disp%show("call setSwapped(a, b)")
                        call setSwapped(a, b)
        call disp%show("a")
        call disp%show( a , deliml = SK_"""" )
        call disp%show("b")
        call disp%show( b , deliml = SK_"""" )
        call disp%skip()
    end block

    block
        complex, allocatable :: a, b
        call disp%skip()
        call disp%show("a = 1")
                        a = 1
        call disp%show("b = 2")
                        b = 2
        call disp%show("call setSwapped(a, b)")
                        call setSwapped(a, b)
        call disp%show("a")
        call disp%show( a , deliml = SK_"""" )
        call disp%show("b")
        call disp%show( b , deliml = SK_"""" )
        call disp%skip()
    end block

    block
        real, allocatable :: a, b
        call disp%skip()
        call disp%show("a = 1")
                        a = 1
        call disp%show("b = 2")
                        b = 2
        call disp%show("call setSwapped(a, b)")
                        call setSwapped(a, b)
        call disp%show("a")
        call disp%show( a , deliml = SK_"""" )
        call disp%show("b")
        call disp%show( b , deliml = SK_"""" )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("! Swap vectors.")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip()

    block
        character(6), allocatable :: a(:), b(:)
        call disp%skip()
        call disp%show("a = ['Heaven', 'Hell  ']")
                        a = ['Heaven', 'Hell  ']
        call disp%show("b = ['Hell  ', 'Heaven']")
                        b = ['Hell  ', 'Heaven']
        call disp%show("call setSwapped(a, b)")
                        call setSwapped(a, b)
        call disp%show("a")
        call disp%show( a , deliml = SK_"""" )
        call disp%show("b")
        call disp%show( b , deliml = SK_"""" )
        call disp%skip()
    end block

    block
        integer, allocatable :: a(:), b(:)
        call disp%skip()
        call disp%show("a = [1, 2]")
                        a = [1, 2]
        call disp%show("b = [2, 1]")
                        b = [2, 1]
        call disp%show("call setSwapped(a, b)")
                        call setSwapped(a, b)
        call disp%show("a")
        call disp%show( a , deliml = SK_"""" )
        call disp%show("b")
        call disp%show( b , deliml = SK_"""" )
        call disp%skip()
    end block

    block
        logical, allocatable :: a(:), b(:)
        call disp%skip()
        call disp%show("a = [.false., .true.]")
                        a = [.false., .true.]
        call disp%show("b = [.true., .false.]")
                        b = [.true., .false.]
        call disp%show("call setSwapped(a, b)")
                        call setSwapped(a, b)
        call disp%show("a")
        call disp%show( a , deliml = SK_"""" )
        call disp%show("b")
        call disp%show( b , deliml = SK_"""" )
        call disp%skip()
    end block

    block
        complex, allocatable :: a(:), b(:)
        call disp%skip()
        call disp%show("a = [1, 2]")
                        a = [1, 2]
        call disp%show("b = [2, 1]")
                        b = [2, 1]
        call disp%show("call setSwapped(a, b)")
                        call setSwapped(a, b)
        call disp%show("a")
        call disp%show( a , deliml = SK_"""" )
        call disp%show("b")
        call disp%show( b , deliml = SK_"""" )
        call disp%skip()
    end block

    block
        real, allocatable :: a(:), b(:)
        call disp%skip()
        call disp%show("a = [1, 2]")
                        a = [1, 2]
        call disp%show("b = [2, 1]")
                        b = [2, 1]
        call disp%show("call setSwapped(a, b)")
                        call setSwapped(a, b)
        call disp%show("a")
        call disp%show( a , deliml = SK_"""" )
        call disp%show("b")
        call disp%show( b , deliml = SK_"""" )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Swap strided vectors.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        integer(IK) , parameter :: DUM = huge(DUM)
        integer, allocatable :: a(:), b(:)
        call disp%skip()
        call disp%show("a = [1, 2, 3, 4, 5]")
                        a = [1, 2, 3, 4, 5]
        call disp%show("b = [-1, DUM , -2, DUM , -3, DUM , -4, DUM , -5]")
                        b = [-1, DUM , -2, DUM , -3, DUM , -4, DUM , -5]
        call disp%show("call setSwapped(a, b, inca = 1_IK, incb = 2_IK)")
                        call setSwapped(a, b, inca = 1_IK, incb = 2_IK)
        call disp%show("a")
        call disp%show( a , deliml = SK_"""" )
        call disp%show("b")
        call disp%show( b , deliml = SK_"""" )
        call disp%skip()
    end block

    block
        integer(IK) , parameter :: DUM = huge(DUM)
        integer, allocatable :: a(:), b(:)
        call disp%skip()
        call disp%show("a = [1, 2, 3, 4, 5]")
                        a = [1, 2, 3, 4, 5]
        call disp%show("b = [-1, DUM , -2, DUM , -3, DUM , -4, DUM , -5]")
                        b = [-1, DUM , -2, DUM , -3, DUM , -4, DUM , -5]
        call disp%show("call setSwapped(a, b, inca = 1_IK, incb = -2_IK)")
                        call setSwapped(a, b, inca = 1_IK, incb = -2_IK)
        call disp%show("a")
        call disp%show( a , deliml = SK_"""" )
        call disp%show("b")
        call disp%show( b , deliml = SK_"""" )
        call disp%skip()
    end block

    block
        integer(IK) , parameter :: DUM = huge(DUM)
        integer, allocatable :: a(:), b(:)
        call disp%skip()
        call disp%show("a = [1, 2, 3, 4, 5]")
                        a = [1, 2, 3, 4, 5]
        call disp%show("b = [-1, DUM , -2, DUM , -3, DUM , -4, DUM , -5]")
                        b = [-1, DUM , -2, DUM , -3, DUM , -4, DUM , -5]
        call disp%show("call setSwapped(a(1:1), b, inca = 0_IK, incb = 2_IK)")
                        call setSwapped(a(1:1), b, inca = 0_IK, incb = 2_IK)
        call disp%show("a")
        call disp%show( a , deliml = SK_"""" )
        call disp%show("b")
        call disp%show( b , deliml = SK_"""" )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Swap strided vectors (blass call).")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real(IK) , parameter :: DUM = huge(DUM)
        real, allocatable :: a(:), b(:)
        call disp%skip()
        call disp%show("a = [real :: 1, 2, 3, 4, 5]")
                        a = [real :: 1, 2, 3, 4, 5]
        call disp%show("b = [real :: -1, DUM , -2, DUM , -3, DUM , -4, DUM , -5]")
                        b = [real :: -1, DUM , -2, DUM , -3, DUM , -4, DUM , -5]
        call disp%show("call setSwapped(a, b, inca = 1_IK, incb = 2_IK)")
                        call setSwapped(a, b, inca = 1_IK, incb = 2_IK)
        call disp%show("a")
        call disp%show( a , deliml = SK_"""" )
        call disp%show("b")
        call disp%show( b , deliml = SK_"""" )
        call disp%skip()
    end block

    block
        integer(IK) , parameter :: DUM = huge(DUM)
        real, allocatable :: a(:), b(:)
        call disp%skip()
        call disp%show("a = [real :: 1, 2, 3, 4, 5]")
                        a = [real :: 1, 2, 3, 4, 5]
        call disp%show("b = [real :: -1, DUM , -2, DUM , -3, DUM , -4, DUM , -5]")
                        b = [real :: -1, DUM , -2, DUM , -3, DUM , -4, DUM , -5]
        call disp%show("call setSwapped(a, b, inca = 1_IK, incb = -2_IK)")
                        call setSwapped(a, b, inca = 1_IK, incb = -2_IK)
        call disp%show("a")
        call disp%show( a , deliml = SK_"""" )
        call disp%show("b")
        call disp%show( b , deliml = SK_"""" )
        call disp%skip()
    end block

    block
        real(IK) , parameter :: DUM = huge(DUM)
        real, allocatable :: a(:), b(:)
        call disp%skip()
        call disp%show("a = [real :: 1, 2, 3, 4, 5]")
                        a = [real :: 1, 2, 3, 4, 5]
        call disp%show("b = [real :: -1, DUM , -2, DUM , -3, DUM , -4, DUM , -5]")
                        b = [real :: -1, DUM , -2, DUM , -3, DUM , -4, DUM , -5]
        call disp%show("call setSwapped(a, b(1:1), inca = 1_IK, incb = 0_IK)")
                        call setSwapped(a, b(1:1), inca = 1_IK, incb = 0_IK)
        call disp%show("a")
        call disp%show( a , deliml = SK_"""" )
        call disp%show("b")
        call disp%show( b , deliml = SK_"""" )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Swap higher-rank arrays.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        integer(IK) , parameter :: DUM = huge(DUM)
        integer, allocatable :: a(:,:), b(:,:)
        call disp%skip()
        call disp%show("a = reshape([1, 2, 3, 4, 5, 6], [2,3])")
                        a = reshape([1, 2, 3, 4, 5, 6], [2,3])
        call disp%show("b = reshape([-1 , -2 , -3, DUM, DUM, DUM, -4, -5, -6], [3, 3], order = [2, 1])")
                        b = reshape([-1 , -2 , -3, DUM, DUM, DUM, -4, -5, -6], [3, 3], order = [2, 1])
        call disp%show("a")
        call disp%show( a , deliml = SK_"""" )
        call disp%show("b")
        call disp%show( b , deliml = SK_"""" )
        call disp%show("call setSwapped(a, b(1::2, :))")
                        call setSwapped(a, b(1::2, :))
        call disp%show("a")
        call disp%show( a , deliml = SK_"""" )
        call disp%show("b")
        call disp%show( b , deliml = SK_"""" )
        call disp%skip()
    end block

end program example