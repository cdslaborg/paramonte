program example

    use pm_kind, only: LK ! The default `logical` kind.
    use pm_kind, only: SK ! The default `character` kind.
    use pm_kind, only: IK ! All other kinds are also supported.
    use pm_kind, only: CK ! All other kinds are also supported.
    use pm_kind, only: RK ! All other kinds are also supported.
    use pm_arrayCompareLex, only: operator(.llt.)
    use pm_distUnif, only: setUnifRand
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Lexically compare two scalar strings with the `< ` operator.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("SK_'ParaMonte' .llt. SK_'ParaMonte'")
    call disp%show( SK_'ParaMonte' .llt. SK_'ParaMonte' )
    call disp%skip()

    call disp%show("SK_'ParaMonte' .llt. SK_'Paramonte'")
    call disp%show( SK_'ParaMonte' .llt. SK_'Paramonte' )
    call disp%skip()

    call disp%show("SK_'ParaMonte' .llt. SK_'ParaMonteX'")
    call disp%show( SK_'ParaMonte' .llt. SK_'ParaMonteX' )
    call disp%skip()

    call disp%show("SK_'ParaMonteX' .llt. SK_'ParaMonte'")
    call disp%show( SK_'ParaMonteX' .llt. SK_'ParaMonte' )
    call disp%skip()

    call disp%show("SK_'ParaMonte' .llt. SK_'ParaMonte ' ! Compare with the Fortran intrinsic lexical comparison below.")
    call disp%show( SK_'ParaMonte' .llt. SK_'ParaMonte ' )
    call disp%skip()

    call disp%show("SK_'ParaMonte' < SK_'ParaMonte ' ! Compare with the lexical comparison above.")
    call disp%show( SK_'ParaMonte' < SK_'ParaMonte ' )
    call disp%skip()

    call disp%show("SK_'ParaMonte ' .llt. SK_'ParaMonte' ! Compare with the Fortran intrinsic lexical comparison below.")
    call disp%show( SK_'ParaMonte ' .llt. SK_'ParaMonte' )
    call disp%skip()

    call disp%show("SK_'ParaMonte ' < SK_'ParaMonte' ! Compare with the lexical comparison above.")
    call disp%show( SK_'ParaMonte ' < SK_'ParaMonte' )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Lexically compare two vector of logical with the `< ` operator.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("[.false._LK] .llt. [.false._LK]")
    call disp%show( [.false._LK] .llt. [.false._LK] )
    call disp%skip()

    call disp%show("[.false._LK] .llt. [.true._LK]")
    call disp%show( [.false._LK] .llt. [.true._LK] )
    call disp%skip()

    call disp%show("[.true._LK] .llt. [.false._LK]")
    call disp%show( [.true._LK] .llt. [.false._LK] )
    call disp%skip()

    call disp%show("[.true._LK] .llt. [.true._LK]")
    call disp%show( [.true._LK] .llt. [.true._LK] )
    call disp%skip()

    call disp%show("[.true._LK] .llt. [.true._LK, .true._LK]")
    call disp%show( [.true._LK] .llt. [.true._LK, .true._LK] )
    call disp%skip()

    call disp%show("[.true._LK] .llt. [.false._LK, .true._LK]")
    call disp%show( [.true._LK] .llt. [.false._LK, .true._LK] )
    call disp%skip()

    call disp%show("[.true._LK, .true._LK] .llt. [.true._LK, .true._LK]")
    call disp%show( [.true._LK, .true._LK] .llt. [.true._LK, .true._LK] )
    call disp%skip()

    call disp%show("[.true._LK, .false._LK] .llt. [.true._LK, .true._LK]")
    call disp%show( [.true._LK, .false._LK] .llt. [.true._LK, .true._LK] )
    call disp%skip()

    call disp%show("[.true._LK, .true._LK] .llt. [.true._LK, .false._LK]")
    call disp%show( [.true._LK, .true._LK] .llt. [.true._LK, .false._LK] )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Lexically compare two vector of integer with the `< ` operator.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("[+0_IK] .llt. [-1_IK]")
    call disp%show( [+0_IK] .llt. [-1_IK] )
    call disp%skip()

    call disp%show("[+0_IK] .llt. [+0_IK]")
    call disp%show( [+0_IK] .llt. [+0_IK] )
    call disp%skip()

    call disp%show("[+0_IK] .llt. [+1_IK]")
    call disp%show( [+0_IK] .llt. [+1_IK] )
    call disp%skip()

    call disp%show("[+0_IK] .llt. [+0_IK, -1_IK]")
    call disp%show( [+0_IK] .llt. [+0_IK, -1_IK] )
    call disp%skip()

    call disp%show("[+0_IK] .llt. [+0_IK, +0_IK]")
    call disp%show( [+0_IK] .llt. [+0_IK, +0_IK] )
    call disp%skip()

    call disp%show("[+0_IK] .llt. [+0_IK, +1_IK]")
    call disp%show( [+0_IK] .llt. [+0_IK, +1_IK] )
    call disp%skip()

    call disp%show("[+0_IK, +0_IK] .llt. [-1_IK, -1_IK]")
    call disp%show( [+0_IK, +0_IK] .llt. [-1_IK, -1_IK] )
    call disp%skip()

    call disp%show("[+0_IK, +0_IK] .llt. [-1_IK, +0_IK]")
    call disp%show( [+0_IK, +0_IK] .llt. [-1_IK, +0_IK] )
    call disp%skip()

    call disp%show("[+0_IK, +0_IK] .llt. [-1_IK, +1_IK]")
    call disp%show( [+0_IK, +0_IK] .llt. [-1_IK, +1_IK] )
    call disp%skip()

    call disp%show("[+0_IK, +0_IK] .llt. [+0_IK, -1_IK]")
    call disp%show( [+0_IK, +0_IK] .llt. [+0_IK, -1_IK] )
    call disp%skip()

    call disp%show("[+0_IK, +0_IK] .llt. [+0_IK, +0_IK]")
    call disp%show( [+0_IK, +0_IK] .llt. [+0_IK, +0_IK] )
    call disp%skip()

    call disp%show("[+0_IK, +0_IK] .llt. [+0_IK, +1_IK]")
    call disp%show( [+0_IK, +0_IK] .llt. [+0_IK, +1_IK] )
    call disp%skip()

    call disp%show("[+0_IK, +0_IK] .llt. [+1_IK, -1_IK]")
    call disp%show( [+0_IK, +0_IK] .llt. [+1_IK, -1_IK] )
    call disp%skip()

    call disp%show("[+0_IK, +0_IK] .llt. [+1_IK, +0_IK]")
    call disp%show( [+0_IK, +0_IK] .llt. [+1_IK, +0_IK] )
    call disp%skip()

    call disp%show("[+0_IK, +0_IK] .llt. [+1_IK, +1_IK]")
    call disp%show( [+0_IK, +0_IK] .llt. [+1_IK, +1_IK] )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Lexically compare two vector of real with the `< ` operator.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("[+0._RK] .llt. [-1._RK]")
    call disp%show( [+0._RK] .llt. [-1._RK] )
    call disp%skip()

    call disp%show("[+0._RK] .llt. [+0._RK]")
    call disp%show( [+0._RK] .llt. [+0._RK] )
    call disp%skip()

    call disp%show("[+0._RK] .llt. [+1._RK]")
    call disp%show( [+0._RK] .llt. [+1._RK] )
    call disp%skip()

    call disp%show("[+0._RK] .llt. [+0._RK, -1._RK]")
    call disp%show( [+0._RK] .llt. [+0._RK, -1._RK] )
    call disp%skip()

    call disp%show("[+0._RK] .llt. [+0._RK, +0._RK]")
    call disp%show( [+0._RK] .llt. [+0._RK, +0._RK] )
    call disp%skip()

    call disp%show("[+0._RK] .llt. [+0._RK, +1._RK]")
    call disp%show( [+0._RK] .llt. [+0._RK, +1._RK] )
    call disp%skip()

    call disp%show("[+0._RK, +0._RK] .llt. [-1._RK, -1._RK]")
    call disp%show( [+0._RK, +0._RK] .llt. [-1._RK, -1._RK] )
    call disp%skip()

    call disp%show("[+0._RK, +0._RK] .llt. [-1._RK, +0._RK]")
    call disp%show( [+0._RK, +0._RK] .llt. [-1._RK, +0._RK] )
    call disp%skip()

    call disp%show("[+0._RK, +0._RK] .llt. [-1._RK, +1._RK]")
    call disp%show( [+0._RK, +0._RK] .llt. [-1._RK, +1._RK] )
    call disp%skip()

    call disp%show("[+0._RK, +0._RK] .llt. [+0._RK, -1._RK]")
    call disp%show( [+0._RK, +0._RK] .llt. [+0._RK, -1._RK] )
    call disp%skip()

    call disp%show("[+0._RK, +0._RK] .llt. [+0._RK, +0._RK]")
    call disp%show( [+0._RK, +0._RK] .llt. [+0._RK, +0._RK] )
    call disp%skip()

    call disp%show("[+0._RK, +0._RK] .llt. [+0._RK, +1._RK]")
    call disp%show( [+0._RK, +0._RK] .llt. [+0._RK, +1._RK] )
    call disp%skip()

    call disp%show("[+0._RK, +0._RK] .llt. [+1._RK, -1._RK]")
    call disp%show( [+0._RK, +0._RK] .llt. [+1._RK, -1._RK] )
    call disp%skip()

    call disp%show("[+0._RK, +0._RK] .llt. [+1._RK, +0._RK]")
    call disp%show( [+0._RK, +0._RK] .llt. [+1._RK, +0._RK] )
    call disp%skip()

    call disp%show("[+0._RK, +0._RK] .llt. [+1._RK, +1._RK]")
    call disp%show( [+0._RK, +0._RK] .llt. [+1._RK, +1._RK] )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Lexically compare two vector of integer with the `< ` operator.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("[(+0._RK, +0._RK)] .llt. [(-1._RK, -1._RK)]")
    call disp%show( [(+0._RK, +0._RK)] .llt. [(-1._RK, -1._RK)] )
    call disp%skip()

    call disp%show("[(+0._RK, +0._RK)] .llt. [(-1._RK, +1._RK)]")
    call disp%show( [(+0._RK, +0._RK)] .llt. [(-1._RK, +1._RK)] )
    call disp%skip()

    call disp%show("[(+0._RK, +0._RK)] .llt. [(+0._RK, +0._RK)]")
    call disp%show( [(+0._RK, +0._RK)] .llt. [(+0._RK, +0._RK)] )
    call disp%skip()

    call disp%show("[(+0._RK, +0._RK)] .llt. [(+1._RK, +1._RK)]")
    call disp%show( [(+0._RK, +0._RK)] .llt. [(+1._RK, +1._RK)] )
    call disp%skip()

    call disp%show("[(+0._RK, +0._RK)] .llt. [(+1._RK, -1._RK)]")
    call disp%show( [(+0._RK, +0._RK)] .llt. [(+1._RK, -1._RK)] )
    call disp%skip()

    call disp%show("[(+0._RK, +0._RK)] .llt. [(+0._RK, +1._RK), (-1._RK, -1._RK)]")
    call disp%show( [(+0._RK, +0._RK)] .llt. [(+0._RK, +1._RK), (-1._RK, -1._RK)] )
    call disp%skip()

    call disp%show("[(+0._RK, +0._RK), (+0._RK, +0._RK)] .llt. [(+0._RK, +0._RK), (-1._RK, -1._RK)]")
    call disp%show( [(+0._RK, +0._RK), (+0._RK, +0._RK)] .llt. [(+0._RK, +0._RK), (-1._RK, -1._RK)] )
    call disp%skip()

    call disp%show("[(+0._RK, +0._RK), (+0._RK, +0._RK)] .llt. [(+0._RK, +0._RK), (+0._RK, -1._RK)]")
    call disp%show( [(+0._RK, +0._RK), (+0._RK, +0._RK)] .llt. [(+0._RK, +0._RK), (+0._RK, -1._RK)] )
    call disp%skip()

    call disp%show("[(+0._RK, +0._RK), (+0._RK, +0._RK)] .llt. [(+0._RK, +0._RK), (+0._RK, +1._RK)]")
    call disp%show( [(+0._RK, +0._RK), (+0._RK, +0._RK)] .llt. [(+0._RK, +0._RK), (+0._RK, +1._RK)] )
    call disp%skip()

    call disp%show("[(+0._RK, +0._RK), (+0._RK, +0._RK)] .llt. [(+0._RK, +0._RK), (+1._RK, +1._RK)]")
    call disp%show( [(+0._RK, +0._RK), (+0._RK, +0._RK)] .llt. [(+0._RK, +0._RK), (+1._RK, +1._RK)] )
    call disp%skip()

end program example