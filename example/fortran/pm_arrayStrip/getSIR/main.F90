program example

    use pm_kind, only: SK, IK, LK, CK, RK
    use pm_io, only: display_type
    use pm_arrayStrip, only: getSIR

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
        call disp%show("index = getSIR(string, pattern)")
                        index = getSIR(string, pattern)
        call disp%show("string(:index)")
        call disp%show( string(:index), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("string = 'bbaaa'")
                        string = 'bbaaa'
        call disp%show("pattern = ''")
                        pattern = ''
        call disp%show("index = getSIR(string, pattern)")
                        index = getSIR(string, pattern)
        call disp%show("string(:index)")
        call disp%show( string(:index), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("string = 'bbaaa'")
                        string = 'bbaaa'
        call disp%show("pattern = 'a'")
                        pattern = 'a'
        call disp%show("index = getSIR(string, pattern)")
                        index = getSIR(string, pattern)
        call disp%show("string(:index)")
        call disp%show( string(:index), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("string = 'bbaaa'")
                        string = 'bbaaa'
        call disp%show("pattern = 'b'")
                        pattern = 'b'
        call disp%show("index = getSIR(string, pattern)")
                        index = getSIR(string, pattern)
        call disp%show("string(:index)")
        call disp%show( string(:index), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("string = 'bbaaa'")
                        string = 'bbaaa'
        call disp%show("pattern = 'bbaaa '")
                        pattern = 'bbaaa '
        call disp%show("index = getSIR(string, pattern)")
                        index = getSIR(string, pattern)
        call disp%show("string(:index)")
        call disp%show( string(:index), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("string = 'bbAa'")
                        string = 'bbAa'
        call disp%show("pattern = 'aa'")
                        pattern = 'aa'
        call disp%show("index = getSIR(string, pattern)")
                        index = getSIR(string, pattern)
        call disp%show("string(:index)")
        call disp%show( string(:index), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%show("index = getSIR(string, pattern, iseqScalar_SK)")
                        index = getSIR(string, pattern, iseqScalar_SK)
        call disp%show("string(:index)")
        call disp%show( string(:index), deliml = """" )
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
        character(2), allocatable   :: Array(:)
        character(:), allocatable   :: scalarPattern
        character(:), allocatable   :: arrayPattern(:)

        call disp%skip()
        call disp%show("Array = [character(2) :: ]")
                        Array = [character(2) :: ]
        call disp%show("scalarPattern = ''")
                        scalarPattern = ''
        call disp%show("index = getSIR(Array, scalarPattern)")
                        index = getSIR(Array, scalarPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [character(2) :: ]")
                        Array = [character(2) :: ]
        call disp%show("arrayPattern = [character(1) :: ]")
                        arrayPattern = [character(1) :: ]
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = ['bb', 'a ', 'a ']")
                        Array = ['bb', 'a ', 'a ']
        call disp%show("arrayPattern = [character(1) :: ]")
                        arrayPattern = [character(1) :: ]
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = ['bb', 'a ', 'a ']")
                        Array = ['bb', 'a ', 'a ']
        call disp%show("scalarPattern = 'a'")
                        scalarPattern = 'a'
        call disp%show("index = getSIR(Array, scalarPattern)")
                        index = getSIR(Array, scalarPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = ['bb', 'a ', 'a ']")
                        Array = ['bb', 'a ', 'a ']
        call disp%show("arrayPattern = ['a', 'a']")
                        arrayPattern = ['a', 'a']
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = ['bb', 'a ', 'a ']")
                        Array = ['bb', 'a ', 'a ']
        call disp%show("arrayPattern = 'bbb'")
                        arrayPattern = 'bbb'
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = ['bb', 'a ', 'a ']")
                        Array = ['bb', 'a ', 'a ']
        call disp%show("arrayPattern = [Array, '  ']")
                        arrayPattern = [Array, '  ']
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = ['a ', 'A ', 'bb']")
                        Array = ['a ', 'A ', 'bb']
        call disp%show("scalarPattern = 'A '")
                        scalarPattern = 'A '
        call disp%show("index = getSIR(Array, scalarPattern)")
                        index = getSIR(Array, scalarPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%show("index = getSIR(Array, scalarPattern, iseqScalar_SK)")
                        index = getSIR(Array, scalarPattern, iseqScalar_SK)
        call disp%show("Array(:index)")
        call disp%show( Array(:index), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = ['bb', 'a ', 'A ']")
                        Array = ['bb', 'a ', 'A ']
        call disp%show("arrayPattern = ['A ', 'a ']")
                        arrayPattern = ['A ', 'a ']
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index), deliml = """" )
        call disp%show("index")
        call disp%show( index )
        call disp%show("index = getSIR(Array, arrayPattern, iseqVector_SK)")
                        index = getSIR(Array, arrayPattern, iseqVector_SK)
        call disp%show("Array(:index)")
        call disp%show( Array(:index), deliml = """" )
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
        integer, allocatable   :: Array(:)
        integer, allocatable   :: scalarPattern
        integer, allocatable   :: arrayPattern(:)

        call disp%skip()
        call disp%show("Array = [integer :: ]")
                        Array = [integer :: ]
        call disp%show("arrayPattern = [integer :: ]")
                        arrayPattern = [integer :: ]
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [2, 1, 1]")
                        Array = [2, 1, 1]
        call disp%show("arrayPattern = [integer :: ]")
                        arrayPattern = [integer :: ]
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [2, 1, 1]")
                        Array = [2, 1, 1]
        call disp%show("scalarPattern = 1")
                        scalarPattern = 1
        call disp%show("index = getSIR(Array, scalarPattern)")
                        index = getSIR(Array, scalarPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [2, 1, 1]")
                        Array = [2, 1, 1]
        call disp%show("arrayPattern = [1, 1]")
                        arrayPattern = [1, 1]
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [2, 1, 1]")
                        Array = [2, 1, 1]
        call disp%show("arrayPattern = 3")
                        arrayPattern = 3
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [2, 1, 1]")
                        Array = [2, 1, 1]
        call disp%show("arrayPattern = Array")
                        arrayPattern = Array
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [2, 1, 1]")
                        Array = [2, 1, 1]
        call disp%show("arrayPattern = [Array, 0]")
                        arrayPattern = [Array, 0]
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [2, -1, 1]")
                        Array = [2, -1, 1]
        call disp%show("scalarPattern = 1")
                        scalarPattern = 1
        call disp%show("index = getSIR(Array, scalarPattern)")
                        index = getSIR(Array, scalarPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%show("index = getSIR(Array, scalarPattern, iseqScalar_IK)")
                        index = getSIR(Array, scalarPattern, iseqScalar_IK)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [2, -1, 1]")
                        Array = [2, -1, 1]
        call disp%show("arrayPattern = [+1, -1]")
                        arrayPattern = [+1, -1]
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%show("index = getSIR(Array, arrayPattern, iseqVector_IK)")
                        index = getSIR(Array, arrayPattern, iseqVector_IK)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
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
        logical, allocatable   :: Array(:)
        logical, allocatable   :: scalarPattern
        logical, allocatable   :: arrayPattern(:)

        call disp%skip()
        call disp%show("Array = [logical :: ]")
                        Array = [logical :: ]
        call disp%show("arrayPattern = [logical :: ]")
                        arrayPattern = [logical :: ]
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [.true., .false., .false.]")
                        Array = [.true., .false., .false.]
        call disp%show("arrayPattern = [logical :: ]")
                        arrayPattern = [logical :: ]
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [.true., .false., .false.]")
                        Array = [.true., .false., .false.]
        call disp%show("scalarPattern = .false.")
                        scalarPattern = .false.
        call disp%show("index = getSIR(Array, scalarPattern)")
                        index = getSIR(Array, scalarPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [.true., .false., .false.]")
                        Array = [.true., .false., .false.]
        call disp%show("arrayPattern = [.false., .false.]")
                        arrayPattern = [.false., .false.]
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [.true., .false., .false.]")
                        Array = [.true., .false., .false.]
        call disp%show("arrayPattern = [.true., .true.]")
                        arrayPattern = [.true., .true.]
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [.true., .false., .false.]")
                        Array = [.true., .false., .false.]
        call disp%show("arrayPattern = Array")
                        arrayPattern = Array
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [.true., .false., .false.]")
                        Array = [.true., .false., .false.]
        call disp%show("arrayPattern = [Array, .false.]")
                        arrayPattern = [Array, .false.]
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [.false., .true., .true.]")
                        Array = [.false., .true., .true.]
        call disp%show("scalarPattern = .false.")
                        scalarPattern = .false.
        call disp%show("index = getSIR(Array, scalarPattern)")
                        index = getSIR(Array, scalarPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%show("index = getSIR(Array, scalarPattern, iseqScalar_LK)")
                        index = getSIR(Array, scalarPattern, iseqScalar_LK)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [.true., .true., .false.]")
                        Array = [.true., .true., .false.]
        call disp%show("arrayPattern = [.false., .true.]")
                        arrayPattern = [.false., .true.]
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%show("index = getSIR(Array, arrayPattern, iseqVector_LK)")
                        index = getSIR(Array, arrayPattern, iseqVector_LK)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
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
        complex, allocatable   :: Array(:)
        complex, allocatable   :: scalarPattern
        complex, allocatable   :: arrayPattern(:)

        call disp%skip()
        call disp%show("Array = [complex :: ]")
                        Array = [complex :: ]
        call disp%show("arrayPattern = [complex :: ]")
                        arrayPattern = [complex :: ]
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [2, 1, 1]")
                        Array = [2, 1, 1]
        call disp%show("arrayPattern = [complex :: ]")
                        arrayPattern = [complex :: ]
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [2, 1, 1]")
                        Array = [2, 1, 1]
        call disp%show("scalarPattern = 1")
                        scalarPattern = 1
        call disp%show("index = getSIR(Array, scalarPattern)")
                        index = getSIR(Array, scalarPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [2, 1, 1]")
                        Array = [2, 1, 1]
        call disp%show("arrayPattern = [1, 1]")
                        arrayPattern = [1, 1]
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [2, 1, 1]")
                        Array = [2, 1, 1]
        call disp%show("arrayPattern = 3")
                        arrayPattern = 3
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [2, 1, 1]")
                        Array = [2, 1, 1]
        call disp%show("arrayPattern = Array")
                        arrayPattern = Array
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [2, 1, 1]")
                        Array = [2, 1, 1]
        call disp%show("arrayPattern = [Array, cmplx(0)]")
                        arrayPattern = [Array, cmplx(0)]
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [2, -1, 1]")
                        Array = [2, -1, 1]
        call disp%show("scalarPattern = 1")
                        scalarPattern = 1
        call disp%show("index = getSIR(Array, scalarPattern)")
                        index = getSIR(Array, scalarPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%show("index = getSIR(Array, scalarPattern, iseqScalar_CK)")
                        index = getSIR(Array, scalarPattern, iseqScalar_CK)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [2, -1, 1]")
                        Array = [2, -1, 1]
        call disp%show("arrayPattern = [+1, -1]")
                        arrayPattern = [+1, -1]
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%show("index = getSIR(Array, arrayPattern, iseqVector_CK)")
                        index = getSIR(Array, arrayPattern, iseqVector_CK)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
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
        real, allocatable   :: Array(:)
        real, allocatable   :: scalarPattern
        real, allocatable   :: arrayPattern(:)

        call disp%skip()
        call disp%show("Array = [real :: ]")
                        Array = [real :: ]
        call disp%show("arrayPattern = [real :: ]")
                        arrayPattern = [real :: ]
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [2, 1, 1]")
                        Array = [2, 1, 1]
        call disp%show("arrayPattern = [real :: ]")
                        arrayPattern = [real :: ]
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [2, 1, 1]")
                        Array = [2, 1, 1]
        call disp%show("scalarPattern = 1")
                        scalarPattern = 1
        call disp%show("index = getSIR(Array, scalarPattern)")
                        index = getSIR(Array, scalarPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [2, 1, 1]")
                        Array = [2, 1, 1]
        call disp%show("arrayPattern = [1, 1]")
                        arrayPattern = [1, 1]
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [2, 1, 1]")
                        Array = [2, 1, 1]
        call disp%show("arrayPattern = 3")
                        arrayPattern = 3
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [2, 1, 1]")
                        Array = [2, 1, 1]
        call disp%show("arrayPattern = Array")
                        arrayPattern = Array
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [2, 1, 1]")
                        Array = [2, 1, 1]
        call disp%show("arrayPattern = [Array, real(0)]")
                        arrayPattern = [Array, real(0)]
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [2, -1, 1]")
                        Array = [2, -1, 1]
        call disp%show("scalarPattern = 1")
                        scalarPattern = 1
        call disp%show("index = getSIR(Array, scalarPattern)")
                        index = getSIR(Array, scalarPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%show("index = getSIR(Array, scalarPattern, iseqScalar_RK)")
                        index = getSIR(Array, scalarPattern, iseqScalar_RK)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%skip()

        call disp%skip()
        call disp%show("Array = [2, -1, 1]")
                        Array = [2, -1, 1]
        call disp%show("arrayPattern = [+1, -1]")
                        arrayPattern = [+1, -1]
        call disp%show("index = getSIR(Array, arrayPattern)")
                        index = getSIR(Array, arrayPattern)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
        call disp%show("index")
        call disp%show( index )
        call disp%show("index = getSIR(Array, arrayPattern, iseqVector_RK)")
                        index = getSIR(Array, arrayPattern, iseqVector_RK)
        call disp%show("Array(:index)")
        call disp%show( Array(:index) )
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