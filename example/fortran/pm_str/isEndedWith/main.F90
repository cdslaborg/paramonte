program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! all processor types and kinds are supported.
    use pm_io, only: display_type
    use pm_str, only: isEndedWith
    use pm_str, only: css_type

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isEndedWith('', '')")
    call disp%show( isEndedWith('', '') )
    call disp%show("isEndedWith('ParaMonte', '')")
    call disp%show( isEndedWith('ParaMonte', '') )
    call disp%show("isEndedWith('', 'ParaMonte')")
    call disp%show( isEndedWith('', 'ParaMonte') )
    call disp%show("isEndedWith('ParaMonte', 'monte')")
    call disp%show( isEndedWith('ParaMonte', 'monte') )
    call disp%show("isEndedWith('ParaMonte', 'Monte')")
    call disp%show( isEndedWith('ParaMonte', 'Monte') )
    call disp%show("isEndedWith('ParaMonte', 'Monte   ')")
    call disp%show( isEndedWith('ParaMonte', 'Monte   ') )
    call disp%show("isEndedWith('ParaMonte ', 'Monte  ')")
    call disp%show( isEndedWith('ParaMonte ', 'Monte  ') )
    call disp%show("isEndedWith('ParaMonte   ', 'Monte')")
    call disp%show( isEndedWith('ParaMonte   ', 'Monte') )
    call disp%show("isEndedWith('ParaMonte', ['Monte', 'monte'])")
    call disp%show( isEndedWith('ParaMonte', ['Monte', 'monte']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isEndedWith(css_type(''), css_type(''))")
    call disp%show( isEndedWith(css_type(''), css_type('')) )
    call disp%show("isEndedWith(css_type('ParaMonte'), css_type(''))")
    call disp%show( isEndedWith(css_type('ParaMonte'), css_type('')) )
    call disp%show("isEndedWith(css_type(''), css_type('ParaMonte'))")
    call disp%show( isEndedWith(css_type(''), css_type('ParaMonte')) )
    call disp%show("isEndedWith(css_type('ParaMonte'), css_type('monte'))")
    call disp%show( isEndedWith(css_type('ParaMonte'), css_type('monte')) )
    call disp%show("isEndedWith(css_type('ParaMonte'), css_type('Monte'))")
    call disp%show( isEndedWith(css_type('ParaMonte'), css_type('Monte')) )
    call disp%show("isEndedWith(css_type('ParaMonte'), css_type('Monte   '))")
    call disp%show( isEndedWith(css_type('ParaMonte'), css_type('Monte   ')) )
    call disp%show("isEndedWith(css_type('ParaMonte '), css_type('Monte  '))")
    call disp%show( isEndedWith(css_type('ParaMonte '), css_type('Monte  ')) )
    call disp%show("isEndedWith(css_type('ParaMonte   '), css_type('Monte'))")
    call disp%show( isEndedWith(css_type('ParaMonte   '), css_type('Monte')) )
    !The following example fails on Windows with ifort Version 2021.11.1 Build 20231117_000000.
    !forrtl: severe (408): fort: (8): Attempt to fetch from allocatable variable VAL when it is not allocated.
    !Attempts to reproduce this error in an isolated simple case were unsuccessful.
#if !__INTEL_COMPILER
    call disp%show("isEndedWith(css_type('ParaMonte'), [css_type('ParaMonte'), css_type('Monte  '), css_type('Monte  ', trimmed = .true._LK), css_type('monte')])")
    call disp%show( isEndedWith(css_type('ParaMonte'), [css_type('ParaMonte'), css_type('Monte  '), css_type('Monte  ', trimmed = .true._LK), css_type('monte')]) )
#endif
    call disp%skip()

end program example