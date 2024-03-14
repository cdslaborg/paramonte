program example

    use pm_kind, only: SK, IK, LK, CK, RK
    use pm_io, only: display_type
    use pm_arrayStrip, only: getStripped, left, right, leftRight

    implicit none

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Strip scalar pattern from a string scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        character(:)    , allocatable   :: array, pattern

        call disp%skip()
        call disp%show("array = ''")
                        array = ''
        call disp%show("pattern = ''")
                        pattern = ''
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = 'bbaaa'")
                        array = 'bbaaa'
        call disp%show("pattern = ''")
                        pattern = ''
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = 'bbaaa'")
                        array = 'bbaaa'
        call disp%show("pattern = 'a'")
                        pattern = 'a'
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = 'bbaaa'")
                        array = 'bbaaa'
        call disp%show("pattern = 'b'")
                        pattern = 'b'
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = 'bbaaa'")
                        array = 'bbaaa'
        call disp%show("pattern = 'bbaaa '")
                        pattern = 'bbaaa '
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = 'bbAa'")
                        array = 'bbAa'
        call disp%show("pattern = 'aa'")
                        pattern = 'aa'
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, iseqScalar_SK)")
        call disp%show( getStripped(array, pattern, iseqScalar_SK) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Strip scalar pattern from a string vector.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        character(2), allocatable   :: array(:)
        character(:), allocatable   :: pattern

        call disp%skip()
        call disp%show("array = [character(2) :: ]")
                        array = [character(2) :: ]
        call disp%show("pattern = ''")
                        pattern = ''
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = ['bb', 'a ', 'a ']")
                        array = ['bb', 'a ', 'a ']
        call disp%show("pattern = 'a'")
                        pattern = 'a'
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = ['a ', 'A ', 'bb']")
                        array = ['a ', 'A ', 'bb']
        call disp%show("pattern = 'A '")
                        pattern = 'A '
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%show("getStripped(array, pattern, iseqScalar_SK)")
        call disp%show( getStripped(array, pattern, iseqScalar_SK) )
        call disp%show("getStripped(array, pattern, iseqScalar_SK, leftRight)")
        call disp%show( getStripped(array, pattern, iseqScalar_SK, leftRight) )
        call disp%show("getStripped(array, pattern, iseqScalar_SK, right)")
        call disp%show( getStripped(array, pattern, iseqScalar_SK, right) )
        call disp%show("getStripped(array, pattern, iseqScalar_SK, left)")
        call disp%show( getStripped(array, pattern, iseqScalar_SK, left) )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Strip vector pattern from a string vector.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        character(2), allocatable   :: array(:)
        character(:), allocatable   :: pattern(:)

        call disp%skip()
        call disp%show("array = [character(2) :: ]")
                        array = [character(2) :: ]
        call disp%show("pattern = [character(1) :: ]")
                        pattern = [character(1) :: ]
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = ['bb', 'a ', 'a ']")
                        array = ['bb', 'a ', 'a ']
        call disp%show("pattern = [character(1) :: ]")
                        pattern = [character(1) :: ]
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = ['bb', 'a ', 'a ']")
                        array = ['bb', 'a ', 'a ']
        call disp%show("pattern = ['a', 'a']")
                        pattern = ['a', 'a']
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = ['bb', 'a ', 'a ']")
                        array = ['bb', 'a ', 'a ']
        call disp%show("pattern = 'bbb'")
                        pattern = 'bbb'
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = ['bb', 'a ', 'a ']")
                        array = ['bb', 'a ', 'a ']
        call disp%show("pattern = [array, '  ']")
                        pattern = [array, '  ']
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = ['bb', 'a ', 'A ']")
                        array = ['bb', 'a ', 'A ']
        call disp%show("pattern = ['A ', 'a ']")
                        pattern = ['A ', 'a ']
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%show("getStripped(array, pattern, iseqVector_SK)")
        call disp%show( getStripped(array, pattern, iseqVector_SK) )
        call disp%show("getStripped(array, pattern, iseqVector_SK, leftRight)")
        call disp%show( getStripped(array, pattern, iseqVector_SK, leftRight) )
        call disp%show("getStripped(array, pattern, iseqVector_SK, right)")
        call disp%show( getStripped(array, pattern, iseqVector_SK, right) )
        call disp%show("getStripped(array, pattern, iseqVector_SK, left)")
        call disp%show( getStripped(array, pattern, iseqVector_SK, left) )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Strip scalar pattern an integer array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        integer, allocatable   :: array(:)
        integer, allocatable   :: pattern

        call disp%skip()
        call disp%show("array = [2, 1, 1]")
                        array = [2, 1, 1]
        call disp%show("pattern = 1")
                        pattern = 1
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [2, -1, 1]")
                        array = [2, -1, 1]
        call disp%show("pattern = 1")
                        pattern = 1
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%show("getStripped(array, pattern, iseqScalar_IK)")
        call disp%show( getStripped(array, pattern, iseqScalar_IK) )
        call disp%show("getStripped(array, pattern, iseqScalar_IK, leftRight)")
        call disp%show( getStripped(array, pattern, iseqScalar_IK, leftRight) )
        call disp%show("getStripped(array, pattern, iseqScalar_IK, right)")
        call disp%show( getStripped(array, pattern, iseqScalar_IK, right) )
        call disp%show("getStripped(array, pattern, iseqScalar_IK, left)")
        call disp%show( getStripped(array, pattern, iseqScalar_IK, left) )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Strip vector pattern an integer array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        integer, allocatable   :: array(:)
        integer, allocatable   :: pattern(:)

        call disp%skip()
        call disp%show("array = [integer :: ]")
                        array = [integer :: ]
        call disp%show("pattern = [integer :: ]")
                        pattern = [integer :: ]
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [2, 1, 1]")
                        array = [2, 1, 1]
        call disp%show("pattern = [integer :: ]")
                        pattern = [integer :: ]
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [2, 1, 1]")
                        array = [2, 1, 1]
        call disp%show("pattern = [1, 1]")
                        pattern = [1, 1]
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [2, 1, 1]")
                        array = [2, 1, 1]
        call disp%show("pattern = 3")
                        pattern = 3
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [2, 1, 1]")
                        array = [2, 1, 1]
        call disp%show("pattern = array")
                        pattern = array
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [2, 1, 1]")
                        array = [2, 1, 1]
        call disp%show("pattern = [array, 0]")
                        pattern = [array, 0]
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [2, -1, 1]")
                        array = [2, -1, 1]
        call disp%show("pattern = [+1, -1]")
                        pattern = [+1, -1]
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, iseqVector_IK)")
        call disp%show( getStripped(array, pattern, iseqVector_IK) )
        call disp%show("getStripped(array, pattern, iseqVector_IK, leftRight)")
        call disp%show( getStripped(array, pattern, iseqVector_IK, leftRight) )
        call disp%show("getStripped(array, pattern, iseqVector_IK, right)")
        call disp%show( getStripped(array, pattern, iseqVector_IK, right) )
        call disp%show("getStripped(array, pattern, iseqVector_IK, left)")
        call disp%show( getStripped(array, pattern, iseqVector_IK, left) )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Strip scalar pattern from a logical array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        logical, allocatable   :: array(:)
        logical, allocatable   :: pattern
        logical, allocatable   :: arrayPattern(:)

        call disp%skip()
        call disp%show("array = [.true., .false., .false.]")
                        array = [.true., .false., .false.]
        call disp%show("pattern = .false.")
                        pattern = .false.
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [.false., .true., .true.]")
                        array = [.false., .true., .true.]
        call disp%show("pattern = .false.")
                        pattern = .false.
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%show("getStripped(array, pattern, iseqScalar_LK)")
        call disp%show( getStripped(array, pattern, iseqScalar_LK) )
        call disp%show("getStripped(array, pattern, iseqScalar_LK, leftRight)")
        call disp%show( getStripped(array, pattern, iseqScalar_LK, leftRight) )
        call disp%show("getStripped(array, pattern, iseqScalar_LK, right)")
        call disp%show( getStripped(array, pattern, iseqScalar_LK, right) )
        call disp%show("getStripped(array, pattern, iseqScalar_LK, left)")
        call disp%show( getStripped(array, pattern, iseqScalar_LK, left) )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Strip vector pattern from a logical array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        logical, allocatable   :: array(:)
        logical, allocatable   :: pattern(:)

        call disp%skip()
        call disp%show("array = [logical :: ]")
                        array = [logical :: ]
        call disp%show("pattern = [logical :: ]")
                        pattern = [logical :: ]
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [.true., .false., .false.]")
                        array = [.true., .false., .false.]
        call disp%show("pattern = [logical :: ]")
                        pattern = [logical :: ]
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [.true., .false., .false.]")
                        array = [.true., .false., .false.]
        call disp%show("pattern = [.false., .false.]")
                        pattern = [.false., .false.]
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [.true., .false., .false.]")
                        array = [.true., .false., .false.]
        call disp%show("pattern = [.true., .true.]")
                        pattern = [.true., .true.]
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [.true., .false., .false.]")
                        array = [.true., .false., .false.]
        call disp%show("pattern = array")
                        pattern = array
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [.true., .false., .false.]")
                        array = [.true., .false., .false.]
        call disp%show("pattern = [array, .false.]")
                        pattern = [array, .false.]
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [.true., .true., .false.]")
                        array = [.true., .true., .false.]
        call disp%show("pattern = [.false., .true.]")
                        pattern = [.false., .true.]
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%show("getStripped(array, pattern, iseqVector_LK)")
        call disp%show( getStripped(array, pattern, iseqVector_LK) )
        call disp%show("getStripped(array, pattern, iseqVector_LK, leftRight)")
        call disp%show( getStripped(array, pattern, iseqVector_LK, leftRight) )
        call disp%show("getStripped(array, pattern, iseqVector_LK, right)")
        call disp%show( getStripped(array, pattern, iseqVector_LK, right) )
        call disp%show("getStripped(array, pattern, iseqVector_LK, left)")
        call disp%show( getStripped(array, pattern, iseqVector_LK, left) )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Strip scalar pattern from a complex array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        complex, allocatable   :: array(:)
        complex, allocatable   :: pattern

        call disp%skip()
        call disp%show("array = [2, 1, 1]")
                        array = [2, 1, 1]
        call disp%show("pattern = 1")
                        pattern = 1
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [2, -1, 1]")
                        array = [2, -1, 1]
        call disp%show("pattern = 1")
                        pattern = 1
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%show("getStripped(array, pattern, iseqScalar_CK)")
        call disp%show( getStripped(array, pattern, iseqScalar_CK) )
        call disp%show("getStripped(array, pattern, iseqScalar_CK, leftRight)")
        call disp%show( getStripped(array, pattern, iseqScalar_CK, leftRight) )
        call disp%show("getStripped(array, pattern, iseqScalar_CK, right)")
        call disp%show( getStripped(array, pattern, iseqScalar_CK, right) )
        call disp%show("getStripped(array, pattern, iseqScalar_CK, left)")
        call disp%show( getStripped(array, pattern, iseqScalar_CK, left) )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Strip vector pattern from a complex array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        complex, allocatable   :: array(:)
        complex, allocatable   :: pattern(:)

        call disp%skip()
        call disp%show("array = [complex :: ]")
                        array = [complex :: ]
        call disp%show("pattern = [complex :: ]")
                        pattern = [complex :: ]
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [2, 1, 1]")
                        array = [2, 1, 1]
        call disp%show("pattern = [complex :: ]")
                        pattern = [complex :: ]
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [2, 1, 1]")
                        array = [2, 1, 1]
        call disp%show("pattern = [1, 1]")
                        pattern = [1, 1]
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [2, 1, 1]")
                        array = [2, 1, 1]
        call disp%show("pattern = 3")
                        pattern = 3
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [2, 1, 1]")
                        array = [2, 1, 1]
        call disp%show("pattern = array")
                        pattern = array
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [2, 1, 1]")
                        array = [2, 1, 1]
        call disp%show("pattern = [array, cmplx(0)]")
                        pattern = [array, cmplx(0)]
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [2, -1, 1]")
                        array = [2, -1, 1]
        call disp%show("pattern = [+1, -1]")
                        pattern = [+1, -1]
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%show("getStripped(array, pattern, iseqVector_CK)")
        call disp%show( getStripped(array, pattern, iseqVector_CK) )
        call disp%show("getStripped(array, pattern, iseqVector_CK, leftRight)")
        call disp%show( getStripped(array, pattern, iseqVector_CK, leftRight) )
        call disp%show("getStripped(array, pattern, iseqVector_CK, right)")
        call disp%show( getStripped(array, pattern, iseqVector_CK, right) )
        call disp%show("getStripped(array, pattern, iseqVector_CK, left)")
        call disp%show( getStripped(array, pattern, iseqVector_CK, left) )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Strip scalar pattern from a real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real, allocatable   :: array(:)
        real, allocatable   :: pattern

        call disp%skip()
        call disp%show("array = [2, -1, 1]")
                        array = [2, -1, 1]
        call disp%show("pattern = 1")
                        pattern = 1
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, iseqScalar_RK)")
        call disp%show( getStripped(array, pattern, iseqScalar_RK) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Strip vector pattern from a real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real, allocatable   :: array(:)
        real, allocatable   :: pattern(:)

        call disp%skip()
        call disp%show("array = [real :: ]")
                        array = [real :: ]
        call disp%show("pattern = [real :: ]")
                        pattern = [real :: ]
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [2, 1, 1]")
                        array = [2, 1, 1]
        call disp%show("pattern = [real :: ]")
                        pattern = [real :: ]
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [2, 1, 1]")
                        array = [2, 1, 1]
        call disp%show("pattern = [1, 1]")
                        pattern = [1, 1]
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [2, 1, 1]")
                        array = [2, 1, 1]
        call disp%show("pattern = 3")
                        pattern = 3
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [2, 1, 1]")
                        array = [2, 1, 1]
        call disp%show("pattern = array")
                        pattern = array
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [2, 1, 1]")
                        array = [2, 1, 1]
        call disp%show("pattern = [array, real(0)]")
                        pattern = [array, real(0)]
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [2, -1, 1]")
                        array = [2, -1, 1]
        call disp%show("pattern = [+1, -1]")
                        pattern = [+1, -1]
        call disp%show("getStripped(array, pattern)")
        call disp%show( getStripped(array, pattern) )
        call disp%show("getStripped(array, pattern, leftRight)")
        call disp%show( getStripped(array, pattern, leftRight) )
        call disp%show("getStripped(array, pattern, right)")
        call disp%show( getStripped(array, pattern, right) )
        call disp%show("getStripped(array, pattern, left)")
        call disp%show( getStripped(array, pattern, left) )
        call disp%show("getStripped(array, pattern, iseqVector_RK)")
        call disp%show( getStripped(array, pattern, iseqVector_RK) )
        call disp%show("getStripped(array, pattern, iseqVector_RK, leftRight)")
        call disp%show( getStripped(array, pattern, iseqVector_RK, leftRight) )
        call disp%show("getStripped(array, pattern, iseqVector_RK, right)")
        call disp%show( getStripped(array, pattern, iseqVector_RK, right) )
        call disp%show("getStripped(array, pattern, iseqVector_RK, left)")
        call disp%show( getStripped(array, pattern, iseqVector_RK, left) )
        call disp%skip()

    end block

contains

    pure function iseqScalar_SK(segment, pattern) result(equivalent)
        use pm_strASCII, only: getStrLower
        character(*), intent(in)    :: segment, pattern
        logical(LK)                 :: equivalent
        equivalent = getStrLower(segment) == getStrLower(pattern)
    end

    pure function iseqVector_SK(segment, pattern, sizePattern) result(equivalent)
        use pm_strASCII, only: getStrLower
        integer(IK) , intent(in)    :: sizePattern
        character(*), intent(in)    :: segment(sizePattern), pattern(sizePattern)
        logical(LK)                 :: equivalent
        equivalent = all(getStrLower(segment) == getStrLower(pattern))
    end

    pure function iseqScalar_IK(segment, pattern) result(equivalent)
        integer     , intent(in)    :: segment, pattern
        logical(LK)                 :: equivalent
        equivalent = abs(segment) == abs(pattern)
    end

    pure function iseqVector_IK(segment, pattern, sizePattern) result(equivalent)
        integer(IK) , intent(in)    :: sizePattern
        integer     , intent(in)    :: segment(sizePattern), pattern(sizePattern)
        logical(LK)                 :: equivalent
        equivalent = all(abs(segment) == abs(pattern))
    end

    pure function iseqScalar_LK(segment, pattern) result(equivalent)
        logical     , intent(in)    :: segment, pattern
        logical(LK)                 :: equivalent
        equivalent = segment .eqv. .not. pattern
    end

    pure function iseqVector_LK(segment, pattern, sizePattern) result(equivalent)
        integer(IK) , intent(in)    :: sizePattern
        logical     , intent(in)    :: segment(sizePattern), pattern(sizePattern)
        logical(LK)                 :: equivalent
        equivalent = all(segment .eqv. .not. pattern)
    end

    pure function iseqScalar_CK(segment, pattern) result(equivalent)
        complex     , intent(in)    :: segment, pattern
        logical(LK)                 :: equivalent
        equivalent = abs(segment%re) == abs(pattern%re)
    end

    pure function iseqVector_CK(segment, pattern, sizePattern) result(equivalent)
        integer(IK) , intent(in)    :: sizePattern
        complex     , intent(in)    :: segment(sizePattern), pattern(sizePattern)
        logical(LK)                 :: equivalent
        equivalent = all(abs(segment%re) == abs(pattern%re))
    end

    pure function iseqScalar_RK(segment, pattern) result(equivalent)
        real        , intent(in)    :: segment, pattern
        logical(LK)                 :: equivalent
        equivalent = abs(segment) == abs(pattern)
    end

    pure function iseqVector_RK(segment, pattern, sizePattern) result(equivalent)
        integer(IK) , intent(in)    :: sizePattern
        real        , intent(in)    :: segment(sizePattern), pattern(sizePattern)
        logical(LK)                 :: equivalent
        equivalent = all(abs(segment) == abs(pattern))
    end

end program example