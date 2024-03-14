program example

    use pm_kind, only: SK, IK, LK, CK, RK
    use pm_io, only: display_type
    use pm_arrayStrip, only: getSIL

    implicit none

    integer(IK) :: index

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Strip pattern from the beginning of a string.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        character(:)    , allocatable   :: string, pattern

        call disp%skip()
        call disp%show("string = ''")
                        string = ''
        call disp%show("pattern = ''")
                        pattern = ''
        call disp%show("index = getSIL(string, pattern)")
                        index = getSIL(string, pattern)
        call disp%show("string(index:)")
        call disp%show( string(index:), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("string = 'aaabb'")
                        string = 'aaabb'
        call disp%show("pattern = ''")
                        pattern = ''
        call disp%show("index = getSIL(string, pattern)")
                        index = getSIL(string, pattern)
        call disp%show("string(index:)")
        call disp%show( string(index:), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("string = 'aaabb'")
                        string = 'aaabb'
        call disp%show("pattern = 'a'")
                        pattern = 'a'
        call disp%show("index = getSIL(string, pattern)")
                        index = getSIL(string, pattern)
        call disp%show("string(index:)")
        call disp%show( string(index:), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("string = 'aaabb'")
                        string = 'aaabb'
        call disp%show("pattern = 'b'")
                        pattern = 'b'
        call disp%show("index = getSIL(string, pattern)")
                        index = getSIL(string, pattern)
        call disp%show("string(index:)")
        call disp%show( string(index:), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("string = 'aaabb'")
                        string = 'aaabb'
        call disp%show("pattern = 'aaabb '")
                        pattern = 'aaabb '
        call disp%show("index = getSIL(string, pattern)")
                        index = getSIL(string, pattern)
        call disp%show("string(index:)")
        call disp%show( string(index:), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("string = 'aAbb'")
                        string = 'aAbb'
        call disp%show("pattern = 'aa'")
                        pattern = 'aa'
        call disp%show("index = getSIL(string, pattern)")
                        index = getSIL(string, pattern)
        call disp%show("string(index:)")
        call disp%show( string(index:), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%show("index = getSIL(string, pattern, iseqScalar_SK)")
                        index = getSIL(string, pattern, iseqScalar_SK)
        call disp%show("string(index:)")
        call disp%show( string(index:), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Strip pattern from the beginning of a string array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        character(2), allocatable   :: array(:)
        character(:), allocatable   :: scalarPattern
        character(:), allocatable   :: arrayPattern(:)

        call disp%skip()
        call disp%show("array = [character(2) :: ]")
                        array = [character(2) :: ]
        call disp%show("scalarPattern = ''")
                        scalarPattern = ''
        call disp%show("index = getSIL(array, scalarPattern)")
                        index = getSIL(array, scalarPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [character(2) :: ]")
                        array = [character(2) :: ]
        call disp%show("arrayPattern = [character(1) :: ]")
                        arrayPattern = [character(1) :: ]
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = ['a ', 'a ', 'bb']")
                        array = ['a ', 'a ', 'bb']
        call disp%show("arrayPattern = [character(1) :: ]")
                        arrayPattern = [character(1) :: ]
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = ['a ', 'a ', 'bb']")
                        array = ['a ', 'a ', 'bb']
        call disp%show("scalarPattern = 'a'")
                        scalarPattern = 'a'
        call disp%show("index = getSIL(array, scalarPattern)")
                        index = getSIL(array, scalarPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = ['a ', 'a ', 'bb']")
                        array = ['a ', 'a ', 'bb']
        call disp%show("arrayPattern = ['a', 'a']")
                        arrayPattern = ['a', 'a']
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = ['a ', 'a ', 'bb']")
                        array = ['a ', 'a ', 'bb']
        call disp%show("arrayPattern = 'bbb'")
                        arrayPattern = 'bbb'
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = ['a ', 'a ', 'bb']")
                        array = ['a ', 'a ', 'bb']
        call disp%show("arrayPattern = [array, '  ']")
                        arrayPattern = [array, '  ']
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = ['a ', 'A ', 'bb']")
                        array = ['a ', 'A ', 'bb']
        call disp%show("scalarPattern = 'A '")
                        scalarPattern = 'A '
        call disp%show("index = getSIL(array, scalarPattern)")
                        index = getSIL(array, scalarPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%show("index = getSIL(array, scalarPattern, iseqScalar_SK)")
                        index = getSIL(array, scalarPattern, iseqScalar_SK)
        call disp%show("array(index:)")
        call disp%show( array(index:), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = ['a ', 'A ', 'bb']")
                        array = ['a ', 'A ', 'bb']
        call disp%show("arrayPattern = ['A ', 'a ']")
                        arrayPattern = ['A ', 'a ']
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%show("index = getSIL(array, arrayPattern, iseqVector_SK)")
                        index = getSIL(array, arrayPattern, iseqVector_SK)
        call disp%show("array(index:)")
        call disp%show( array(index:), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Strip pattern from the beginning of a integer array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        integer, allocatable   :: array(:)
        integer, allocatable   :: scalarPattern
        integer, allocatable   :: arrayPattern(:)

        call disp%skip()
        call disp%show("array = [integer :: ]")
                        array = [integer :: ]
        call disp%show("arrayPattern = [integer :: ]")
                        arrayPattern = [integer :: ]
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [1, 1, 2]")
                        array = [1, 1, 2]
        call disp%show("arrayPattern = [integer :: ]")
                        arrayPattern = [integer :: ]
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [1, 1, 2]")
                        array = [1, 1, 2]
        call disp%show("scalarPattern = 1")
                        scalarPattern = 1
        call disp%show("index = getSIL(array, scalarPattern)")
                        index = getSIL(array, scalarPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [1, 1, 2]")
                        array = [1, 1, 2]
        call disp%show("arrayPattern = [1, 1]")
                        arrayPattern = [1, 1]
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [1, 1, 2]")
                        array = [1, 1, 2]
        call disp%show("arrayPattern = 3")
                        arrayPattern = 3
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [1, 1, 2]")
                        array = [1, 1, 2]
        call disp%show("arrayPattern = array")
                        arrayPattern = array
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [1, 1, 2]")
                        array = [1, 1, 2]
        call disp%show("arrayPattern = [array, 0]")
                        arrayPattern = [array, 0]
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [1, -1, 2]")
                        array = [1, -1, 2]
        call disp%show("scalarPattern = 1")
                        scalarPattern = 1
        call disp%show("index = getSIL(array, scalarPattern)")
                        index = getSIL(array, scalarPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%show("index = getSIL(array, scalarPattern, iseqScalar_IK)")
                        index = getSIL(array, scalarPattern, iseqScalar_IK)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [1, -1, 2]")
                        array = [1, -1, 2]
        call disp%show("arrayPattern = [-1, 1]")
                        arrayPattern = [-1, 1]
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%show("index = getSIL(array, arrayPattern, iseqVector_IK)")
                        index = getSIL(array, arrayPattern, iseqVector_IK)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Strip pattern from the beginning of a logical array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        logical, allocatable   :: array(:)
        logical, allocatable   :: scalarPattern
        logical, allocatable   :: arrayPattern(:)

        call disp%skip()
        call disp%show("array = [logical :: ]")
                        array = [logical :: ]
        call disp%show("arrayPattern = [logical :: ]")
                        arrayPattern = [logical :: ]
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [.false., .false., .true.]")
                        array = [.false., .false., .true.]
        call disp%show("arrayPattern = [logical :: ]")
                        arrayPattern = [logical :: ]
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [.false., .false., .true.]")
                        array = [.false., .false., .true.]
        call disp%show("scalarPattern = .false.")
                        scalarPattern = .false.
        call disp%show("index = getSIL(array, scalarPattern)")
                        index = getSIL(array, scalarPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [.false., .false., .true.]")
                        array = [.false., .false., .true.]
        call disp%show("arrayPattern = [.false., .false.]")
                        arrayPattern = [.false., .false.]
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [.false., .false., .true.]")
                        array = [.false., .false., .true.]
        call disp%show("arrayPattern = [.true., .true.]")
                        arrayPattern = [.true., .true.]
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [.false., .false., .true.]")
                        array = [.false., .false., .true.]
        call disp%show("arrayPattern = array")
                        arrayPattern = array
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [.false., .false., .true.]")
                        array = [.false., .false., .true.]
        call disp%show("arrayPattern = [array, .false.]")
                        arrayPattern = [array, .false.]
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [.false., .true., .true.]")
                        array = [.false., .true., .true.]
        call disp%show("scalarPattern = .false.")
                        scalarPattern = .false.
        call disp%show("index = getSIL(array, scalarPattern)")
                        index = getSIL(array, scalarPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%show("index = getSIL(array, scalarPattern, iseqScalar_LK)")
                        index = getSIL(array, scalarPattern, iseqScalar_LK)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [.false., .true., .true.]")
                        array = [.false., .true., .true.]
        call disp%show("arrayPattern = [.true., .false.]")
                        arrayPattern = [.true., .false.]
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%show("index = getSIL(array, arrayPattern, iseqVector_LK)")
                        index = getSIL(array, arrayPattern, iseqVector_LK)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Strip pattern from the beginning of a complex array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        complex, allocatable   :: array(:)
        complex, allocatable   :: scalarPattern
        complex, allocatable   :: arrayPattern(:)

        call disp%skip()
        call disp%show("array = [complex :: ]")
                        array = [complex :: ]
        call disp%show("arrayPattern = [complex :: ]")
                        arrayPattern = [complex :: ]
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [1, 1, 2]")
                        array = [1, 1, 2]
        call disp%show("arrayPattern = [complex :: ]")
                        arrayPattern = [complex :: ]
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [1, 1, 2]")
                        array = [1, 1, 2]
        call disp%show("scalarPattern = 1")
                        scalarPattern = 1
        call disp%show("index = getSIL(array, scalarPattern)")
                        index = getSIL(array, scalarPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [1, 1, 2]")
                        array = [1, 1, 2]
        call disp%show("arrayPattern = [1, 1]")
                        arrayPattern = [1, 1]
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [1, 1, 2]")
                        array = [1, 1, 2]
        call disp%show("arrayPattern = 3")
                        arrayPattern = 3
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [1, 1, 2]")
                        array = [1, 1, 2]
        call disp%show("arrayPattern = array")
                        arrayPattern = array
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [1, 1, 2]")
                        array = [1, 1, 2]
        call disp%show("arrayPattern = [array, cmplx(0)]")
                        arrayPattern = [array, cmplx(0)]
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [1, -1, 2]")
                        array = [1, -1, 2]
        call disp%show("scalarPattern = 1")
                        scalarPattern = 1
        call disp%show("index = getSIL(array, scalarPattern)")
                        index = getSIL(array, scalarPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%show("index = getSIL(array, scalarPattern, iseqScalar_CK)")
                        index = getSIL(array, scalarPattern, iseqScalar_CK)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [1, -1, 2]")
                        array = [1, -1, 2]
        call disp%show("arrayPattern = [-1, 1]")
                        arrayPattern = [-1, 1]
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%show("index = getSIL(array, arrayPattern, iseqVector_CK)")
                        index = getSIL(array, arrayPattern, iseqVector_CK)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Strip pattern from the beginning of a real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real, allocatable   :: array(:)
        real, allocatable   :: scalarPattern
        real, allocatable   :: arrayPattern(:)

        call disp%skip()
        call disp%show("array = [real :: ]")
                        array = [real :: ]
        call disp%show("arrayPattern = [real :: ]")
                        arrayPattern = [real :: ]
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [1, 1, 2]")
                        array = [1, 1, 2]
        call disp%show("arrayPattern = [real :: ]")
                        arrayPattern = [real :: ]
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [1, 1, 2]")
                        array = [1, 1, 2]
        call disp%show("scalarPattern = 1")
                        scalarPattern = 1
        call disp%show("index = getSIL(array, scalarPattern)")
                        index = getSIL(array, scalarPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [1, 1, 2]")
                        array = [1, 1, 2]
        call disp%show("arrayPattern = [1, 1]")
                        arrayPattern = [1, 1]
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [1, 1, 2]")
                        array = [1, 1, 2]
        call disp%show("arrayPattern = 3")
                        arrayPattern = 3
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [1, 1, 2]")
                        array = [1, 1, 2]
        call disp%show("arrayPattern = array")
                        arrayPattern = array
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [1, 1, 2]")
                        array = [1, 1, 2]
        call disp%show("arrayPattern = [array, real(0)]")
                        arrayPattern = [array, real(0)]
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [1, -1, 2]")
                        array = [1, -1, 2]
        call disp%show("scalarPattern = 1")
                        scalarPattern = 1
        call disp%show("index = getSIL(array, scalarPattern)")
                        index = getSIL(array, scalarPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%show("index = getSIL(array, scalarPattern, iseqScalar_RK)")
                        index = getSIL(array, scalarPattern, iseqScalar_RK)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("array = [1, -1, 2]")
                        array = [1, -1, 2]
        call disp%show("arrayPattern = [-1, 1]")
                        arrayPattern = [-1, 1]
        call disp%show("index = getSIL(array, arrayPattern)")
                        index = getSIL(array, arrayPattern)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%show("index = getSIL(array, arrayPattern, iseqVector_RK)")
                        index = getSIL(array, arrayPattern, iseqVector_RK)
        call disp%show("array(index:)")
        call disp%show( array(index:) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

    end block

contains

    pure function iseqScalar_SK(segment, pattern) result(equivalent)
        use pm_strASCII, only: getStrLower
        character(*), intent(in)    :: segment, pattern
        logical(LK)                 :: equivalent
        equivalent = getStrLower(segment) == getStrLower(pattern)
    end

    pure function iseqVector_SK(Segment, pattern, sizePattern) result(equivalent)
        use pm_strASCII, only: getStrLower
        integer(IK) , intent(in)    :: sizePattern
        character(*), intent(in)    :: Segment(sizePattern), pattern(sizePattern)
        logical(LK)                 :: equivalent
        equivalent = all(getStrLower(Segment) == getStrLower(pattern))
    end

    pure function iseqScalar_IK(segment, pattern) result(equivalent)
        integer     , intent(in)    :: segment, pattern
        logical(LK)                 :: equivalent
        equivalent = abs(segment) == abs(pattern)
    end

    pure function iseqVector_IK(Segment, pattern, sizePattern) result(equivalent)
        integer(IK) , intent(in)    :: sizePattern
        integer     , intent(in)    :: Segment(sizePattern), pattern(sizePattern)
        logical(LK)                 :: equivalent
        equivalent = all(abs(Segment) == abs(pattern))
    end

    pure function iseqScalar_LK(segment, pattern) result(equivalent)
        logical     , intent(in)    :: segment, pattern
        logical(LK)                 :: equivalent
        equivalent = segment .eqv. .not. pattern
    end

    pure function iseqVector_LK(Segment, pattern, sizePattern) result(equivalent)
        integer(IK) , intent(in)    :: sizePattern
        logical     , intent(in)    :: Segment(sizePattern), pattern(sizePattern)
        logical(LK)                 :: equivalent
        equivalent = all(Segment .eqv. .not. pattern)
    end

    pure function iseqScalar_CK(segment, pattern) result(equivalent)
        complex     , intent(in)    :: segment, pattern
        logical(LK)                 :: equivalent
        equivalent = abs(segment%re) == abs(pattern%re)
    end

    pure function iseqVector_CK(Segment, pattern, sizePattern) result(equivalent)
        integer(IK) , intent(in)    :: sizePattern
        complex     , intent(in)    :: Segment(sizePattern), pattern(sizePattern)
        logical(LK)                 :: equivalent
        equivalent = all(abs(Segment%re) == abs(pattern%re))
    end

    pure function iseqScalar_RK(segment, pattern) result(equivalent)
        real        , intent(in)    :: segment, pattern
        logical(LK)                 :: equivalent
        equivalent = abs(segment) == abs(pattern)
    end

    pure function iseqVector_RK(Segment, pattern, sizePattern) result(equivalent)
        integer(IK) , intent(in)    :: sizePattern
        real        , intent(in)    :: Segment(sizePattern), pattern(sizePattern)
        logical(LK)                 :: equivalent
        equivalent = all(abs(Segment) == abs(pattern))
    end

end program example