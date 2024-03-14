program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! all processor types and kinds are supported.
    use pm_io, only: display_type
    use pm_arrayReverse, only: getReversed
    use pm_str, only: operator(.alleq.)

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("'' .alleq. ''")
    call disp%show( '' .alleq. '' )
    call disp%skip()

    call disp%skip()
    call disp%show("'paramonte' .alleq. ''")
    call disp%show( 'paramonte' .alleq. '' )
    call disp%skip()

    call disp%skip()
    call disp%show("'paramonte' .alleq. 'paramonte'")
    call disp%show( 'paramonte' .alleq. 'paramonte' )
    call disp%skip()

    call disp%skip()
    call disp%show("'aaa' .alleq. 'AAA'")
    call disp%show( 'aaa' .alleq. 'AAA' )
    call disp%skip()

    call disp%skip()
    call disp%show("'aaa' .alleq. 'aa'")
    call disp%show( 'aaa' .alleq. 'aa' )
    call disp%skip()

    call disp%skip()
    call disp%show("'a' .alleq. 'aaaa'")
    call disp%show( 'a' .alleq. 'aaaa' )
    call disp%skip()

    call disp%skip()
    call disp%show("'aaa' .alleq. 'a'")
    call disp%show( 'aaa' .alleq. 'a' )
    call disp%skip()

    call disp%skip()
    call disp%show("'a' .alleq. 'a'")
    call disp%show( 'a' .alleq. 'a' )
    call disp%skip()

    call disp%skip()
    call disp%show("['aaa', 'AAA'] .alleq. 'a'")
    call disp%show( ['aaa', 'AAA'] .alleq. 'a' )
    call disp%skip()

    call disp%skip()
    call disp%show("['aaa', 'AAA'] .alleq. ['a', 'A']")
    call disp%show( ['aaa', 'AAA'] .alleq. ['a', 'A'] )
    call disp%skip()

    call disp%skip()
    call disp%show("['aaa', 'AAA'] .alleq. getReversed(['a', 'A'])")
    call disp%show( ['aaa', 'AAA'] .alleq. getReversed(['a', 'A']) )
    call disp%skip()

end program example