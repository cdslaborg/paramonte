program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: SKC => SK, IKC => IK32, LKC => LK, RKC => RK32, CKC => CK32 ! all processor types and kinds are supported.
    use pm_matrixCopy, only: setMatCopy, transSymm, transHerm
    use pm_matrixCopy, only: dia, lowDia, uppDia, uppLow
    use pm_matrixCopy, only: rdpack, lfpack, rfpack
    use pm_distUnif, only: getUnifRand, setUnifRand
    use pm_io, only: getFormat
    use pm_io, only: display_type

    implicit none

    character(:, SK), allocatable :: cform
    integer(IK) :: doff

    type(display_type)  :: disp
    disp = display_type(file = "main.out.F90")

    cform = getFormat([cmplx(0., 0., CKC)], ed = SK_'f', ndigit = 1_IK, signed = .true.)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Copy subset to subset in Rectangular Default Format.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        character(2,SKC) :: source(10,10), destin(10,10)
        character(2,SKC), parameter :: EMPTY = SKC_"  "

        call disp%skip()
        call disp%show("call setUnifRand(source, SKC_'AA', SKC_'ZZ')")
                        call setUnifRand(source, SKC_'AA', SKC_'ZZ')
        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:5, 2:6), rdpack, source(2:4, 3:7), rdpack, uppDia)")
                        call setMatCopy(destin(3:5, 2:6), rdpack, source(2:4, 3:7), rdpack, uppDia)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:5, 2:6), rdpack, source(2:4, 3:7), rdpack, dia)")
                        call setMatCopy(destin(3:5, 2:6), rdpack, source(2:4, 3:7), rdpack, dia)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(2:6, 3:5), rdpack, source(2:4, 3:7), rdpack, uppDia, transSymm)")
                        call setMatCopy(destin(2:6, 3:5), rdpack, source(2:4, 3:7), rdpack, uppDia, transSymm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(2:6, 3:5), rdpack, source(2:4, 3:7), rdpack, dia, transSymm)")
                        call setMatCopy(destin(2:6, 3:5), rdpack, source(2:4, 3:7), rdpack, dia, transSymm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(2:6, 3:5), rdpack, source(2:4, 3:7), rdpack, uppDia, transHerm)")
                        call setMatCopy(destin(2:6, 3:5), rdpack, source(2:4, 3:7), rdpack, uppDia, transHerm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(2:6, 3:5), rdpack, source(2:4, 3:7), rdpack, dia, transHerm)")
                        call setMatCopy(destin(2:6, 3:5), rdpack, source(2:4, 3:7), rdpack, dia, transHerm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:5, 2:6), rdpack, source(2:4, 3:7), rdpack, lowDia)")
                        call setMatCopy(destin(3:5, 2:6), rdpack, source(2:4, 3:7), rdpack, lowDia)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:5, 2:6), rdpack, source(2:4, 3:7), rdpack, dia)")
                        call setMatCopy(destin(3:5, 2:6), rdpack, source(2:4, 3:7), rdpack, dia)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(2:6, 3:5), rdpack, source(2:4, 3:7), rdpack, lowDia, transSymm)")
                        call setMatCopy(destin(2:6, 3:5), rdpack, source(2:4, 3:7), rdpack, lowDia, transSymm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(2:6, 3:5), rdpack, source(2:4, 3:7), rdpack, dia, transSymm)")
                        call setMatCopy(destin(2:6, 3:5), rdpack, source(2:4, 3:7), rdpack, dia, transSymm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(2:6, 3:5), rdpack, source(2:4, 3:7), rdpack, lowDia, transHerm)")
                        call setMatCopy(destin(2:6, 3:5), rdpack, source(2:4, 3:7), rdpack, lowDia, transHerm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(2:6, 3:5), rdpack, source(2:4, 3:7), rdpack, dia, transHerm)")
                        call setMatCopy(destin(2:6, 3:5), rdpack, source(2:4, 3:7), rdpack, dia, transHerm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:5, 2:6), rdpack, source(2:4, 3:7), rdpack, uppLow)")
                        call setMatCopy(destin(3:5, 2:6), rdpack, source(2:4, 3:7), rdpack, uppLow)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:5, 2:6), rdpack, source(2:4, 3:7), rdpack, dia)")
                        call setMatCopy(destin(3:5, 2:6), rdpack, source(2:4, 3:7), rdpack, dia)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(2:6, 3:5), rdpack, source(2:4, 3:7), rdpack, uppLow, transSymm)")
                        call setMatCopy(destin(2:6, 3:5), rdpack, source(2:4, 3:7), rdpack, uppLow, transSymm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(2:6, 3:5), rdpack, source(2:4, 3:7), rdpack, dia, transSymm)")
                        call setMatCopy(destin(2:6, 3:5), rdpack, source(2:4, 3:7), rdpack, dia, transSymm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(2:6, 3:5), rdpack, source(2:4, 3:7), rdpack, uppLow, transHerm)")
                        call setMatCopy(destin(2:6, 3:5), rdpack, source(2:4, 3:7), rdpack, uppLow, transHerm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(2:6, 3:5), rdpack, source(2:4, 3:7), rdpack, dia, transHerm)")
                        call setMatCopy(destin(2:6, 3:5), rdpack, source(2:4, 3:7), rdpack, dia, transHerm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 2:4), rdpack, source(2:6, 3:5), rdpack, uppDia)")
                        call setMatCopy(destin(3:7, 2:4), rdpack, source(2:6, 3:5), rdpack, uppDia)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 2:4), rdpack, source(2:6, 3:5), rdpack, dia)")
                        call setMatCopy(destin(3:7, 2:4), rdpack, source(2:6, 3:5), rdpack, dia)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(2:4, 3:7), rdpack, source(2:6, 3:5), rdpack, uppDia, transSymm)")
                        call setMatCopy(destin(2:4, 3:7), rdpack, source(2:6, 3:5), rdpack, uppDia, transSymm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(2:4, 3:7), rdpack, source(2:6, 3:5), rdpack, dia, transSymm)")
                        call setMatCopy(destin(2:4, 3:7), rdpack, source(2:6, 3:5), rdpack, dia, transSymm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(2:4, 3:7), rdpack, source(2:6, 3:5), rdpack, uppDia, transHerm)")
                        call setMatCopy(destin(2:4, 3:7), rdpack, source(2:6, 3:5), rdpack, uppDia, transHerm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(2:4, 3:7), rdpack, source(2:6, 3:5), rdpack, dia, transHerm)")
                        call setMatCopy(destin(2:4, 3:7), rdpack, source(2:6, 3:5), rdpack, dia, transHerm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 2:4), rdpack, source(2:6, 3:5), rdpack, lowDia)")
                        call setMatCopy(destin(3:7, 2:4), rdpack, source(2:6, 3:5), rdpack, lowDia)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 2:4), rdpack, source(2:6, 3:5), rdpack, dia)")
                        call setMatCopy(destin(3:7, 2:4), rdpack, source(2:6, 3:5), rdpack, dia)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(2:4, 3:7), rdpack, source(2:6, 3:5), rdpack, lowDia, transSymm)")
                        call setMatCopy(destin(2:4, 3:7), rdpack, source(2:6, 3:5), rdpack, lowDia, transSymm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(2:4, 3:7), rdpack, source(2:6, 3:5), rdpack, dia, transSymm)")
                        call setMatCopy(destin(2:4, 3:7), rdpack, source(2:6, 3:5), rdpack, dia, transSymm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(2:4, 3:7), rdpack, source(2:6, 3:5), rdpack, lowDia, transHerm)")
                        call setMatCopy(destin(2:4, 3:7), rdpack, source(2:6, 3:5), rdpack, lowDia, transHerm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(2:4, 3:7), rdpack, source(2:6, 3:5), rdpack, dia, transHerm)")
                        call setMatCopy(destin(2:4, 3:7), rdpack, source(2:6, 3:5), rdpack, dia, transHerm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 2:4), rdpack, source(2:6, 3:5), rdpack, uppLow)")
                        call setMatCopy(destin(3:7, 2:4), rdpack, source(2:6, 3:5), rdpack, uppLow)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 2:4), rdpack, source(2:6, 3:5), rdpack, dia)")
                        call setMatCopy(destin(3:7, 2:4), rdpack, source(2:6, 3:5), rdpack, dia)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(2:4, 3:7), rdpack, source(2:6, 3:5), rdpack, uppLow, transSymm)")
                        call setMatCopy(destin(2:4, 3:7), rdpack, source(2:6, 3:5), rdpack, uppLow, transSymm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(2:4, 3:7), rdpack, source(2:6, 3:5), rdpack, dia, transSymm)")
                        call setMatCopy(destin(2:4, 3:7), rdpack, source(2:6, 3:5), rdpack, dia, transSymm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(2:4, 3:7), rdpack, source(2:6, 3:5), rdpack, uppLow, transHerm)")
                        call setMatCopy(destin(2:4, 3:7), rdpack, source(2:6, 3:5), rdpack, uppLow, transHerm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(2:4, 3:7), rdpack, source(2:6, 3:5), rdpack, dia, transHerm)")
                        call setMatCopy(destin(2:4, 3:7), rdpack, source(2:6, 3:5), rdpack, dia, transHerm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, uppDia)")
                        call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, uppDia)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, dia)")
                        call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, dia)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, uppDia, transSymm)")
                        call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, uppDia, transSymm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, dia, transSymm)")
                        call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, dia, transSymm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, uppDia, transHerm)")
                        call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, uppDia, transHerm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, dia, transHerm)")
                        call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, dia, transHerm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, lowDia)")
                        call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, lowDia)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, dia)")
                        call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, dia)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, lowDia, transSymm)")
                        call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, lowDia, transSymm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, dia, transSymm)")
                        call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, dia, transSymm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, lowDia, transHerm)")
                        call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, lowDia, transHerm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, dia, transHerm)")
                        call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, dia, transHerm)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("doff = 0")
                        doff = 0
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, uppLow, doff)")
                        call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, uppLow, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, dia)")
                        call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, dia)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, uppLow, transSymm, doff)")
                        call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, uppLow, transSymm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, dia, transSymm, doff)")
                        call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, dia, transSymm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, uppLow, transHerm, doff)")
                        call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, uppLow, transHerm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, dia, transHerm, doff)")
                        call setMatCopy(destin(3:7, 3:7), rdpack, source(1:5, 1:5), rdpack, dia, transHerm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("doff = -3")
                        doff = -3
        call disp%show("call setMatCopy(destin(3:9, 3:7), rdpack, source(1:7, 1:5), rdpack, uppDia, doff)")
                        call setMatCopy(destin(3:9, 3:7), rdpack, source(1:7, 1:5), rdpack, uppDia, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:9, 3:7), rdpack, source(1:7, 1:5), rdpack, dia, doff)")
                        call setMatCopy(destin(3:9, 3:7), rdpack, source(1:7, 1:5), rdpack, dia, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("doff = -3")
                        doff = -3
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:9), rdpack, source(1:7, 1:5), rdpack, uppDia, transSymm, doff)")
                        call setMatCopy(destin(3:7, 3:9), rdpack, source(1:7, 1:5), rdpack, uppDia, transSymm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:9), rdpack, source(1:7, 1:5), rdpack, dia, transSymm, doff)")
                        call setMatCopy(destin(3:7, 3:9), rdpack, source(1:7, 1:5), rdpack, dia, transSymm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("doff = -3")
                        doff = -3
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:9), rdpack, source(1:7, 1:5), rdpack, uppDia, transHerm, doff)")
                        call setMatCopy(destin(3:7, 3:9), rdpack, source(1:7, 1:5), rdpack, uppDia, transHerm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:9), rdpack, source(1:7, 1:5), rdpack, dia, transHerm, doff)")
                        call setMatCopy(destin(3:7, 3:9), rdpack, source(1:7, 1:5), rdpack, dia, transHerm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("doff = +3")
                        doff = +3
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:9, 3:7), rdpack, source(1:7, 1:5), rdpack, lowDia, doff)")
                        call setMatCopy(destin(3:9, 3:7), rdpack, source(1:7, 1:5), rdpack, lowDia, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:9, 3:7), rdpack, source(1:7, 1:5), rdpack, dia, doff)")
                        call setMatCopy(destin(3:9, 3:7), rdpack, source(1:7, 1:5), rdpack, dia, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("doff = +3")
                        doff = +3
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:9), rdpack, source(1:7, 1:5), rdpack, lowDia, transSymm, doff)")
                        call setMatCopy(destin(3:7, 3:9), rdpack, source(1:7, 1:5), rdpack, lowDia, transSymm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:9), rdpack, source(1:7, 1:5), rdpack, dia, transSymm, doff)")
                        call setMatCopy(destin(3:7, 3:9), rdpack, source(1:7, 1:5), rdpack, dia, transSymm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("doff = +3")
                        doff = +3
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:9), rdpack, source(1:7, 1:5), rdpack, lowDia, transHerm, doff)")
                        call setMatCopy(destin(3:7, 3:9), rdpack, source(1:7, 1:5), rdpack, lowDia, transHerm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:9), rdpack, source(1:7, 1:5), rdpack, dia, transHerm, doff)")
                        call setMatCopy(destin(3:7, 3:9), rdpack, source(1:7, 1:5), rdpack, dia, transHerm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("doff = -4")
                        doff = -4
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:9, 3:7), rdpack, source(1:7, 1:5), rdpack, uppLow, doff)")
                        call setMatCopy(destin(3:9, 3:7), rdpack, source(1:7, 1:5), rdpack, uppLow, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:9, 3:7), rdpack, source(1:7, 1:5), rdpack, dia, doff)")
                        call setMatCopy(destin(3:9, 3:7), rdpack, source(1:7, 1:5), rdpack, dia, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("doff = -4")
                        doff = -4
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:9), rdpack, source(1:7, 1:5), rdpack, uppLow, transSymm, doff)")
                        call setMatCopy(destin(3:7, 3:9), rdpack, source(1:7, 1:5), rdpack, uppLow, transSymm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:9), rdpack, source(1:7, 1:5), rdpack, dia, transSymm, doff)")
                        call setMatCopy(destin(3:7, 3:9), rdpack, source(1:7, 1:5), rdpack, dia, transSymm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("doff = +4")
                        doff = +4
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:9), rdpack, source(1:7, 1:5), rdpack, uppLow, transHerm, doff)")
                        call setMatCopy(destin(3:7, 3:9), rdpack, source(1:7, 1:5), rdpack, uppLow, transHerm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:9), rdpack, source(1:7, 1:5), rdpack, dia, transHerm, doff)")
                        call setMatCopy(destin(3:7, 3:9), rdpack, source(1:7, 1:5), rdpack, dia, transHerm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()



        call disp%skip()
        call disp%show("call setUnifRand(source, SKC_'AA', SKC_'ZZ')")
                        call setUnifRand(source, SKC_'AA', SKC_'ZZ')
        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )

        call disp%show("doff = -3")
                        doff = -3
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:9), rdpack, source(1:5, 1:7), rdpack, uppDia, doff)")
                        call setMatCopy(destin(3:7, 3:9), rdpack, source(1:5, 1:7), rdpack, uppDia, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:9), rdpack, source(1:5, 1:7), rdpack, dia, doff)")
                        call setMatCopy(destin(3:7, 3:9), rdpack, source(1:5, 1:7), rdpack, dia, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("doff = -3")
                        doff = -3
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, uppDia, transSymm, doff)")
                        call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, uppDia, transSymm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, dia, transSymm, doff)")
                        call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, dia, transSymm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("doff = -3")
                        doff = -3
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, uppDia, transHerm, doff)")
                        call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, uppDia, transHerm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, dia, transHerm, doff)")
                        call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, dia, transHerm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("doff = +3")
                        doff = +3
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:9), rdpack, source(1:5, 1:7), rdpack, lowDia, doff)")
                        call setMatCopy(destin(3:7, 3:9), rdpack, source(1:5, 1:7), rdpack, lowDia, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:9), rdpack, source(1:5, 1:7), rdpack, dia, doff)")
                        call setMatCopy(destin(3:7, 3:9), rdpack, source(1:5, 1:7), rdpack, dia, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("doff = +3")
                        doff = +3
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, lowDia, transSymm, doff)")
                        call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, lowDia, transSymm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, dia, transSymm, doff)")
                        call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, dia, transSymm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("doff = +3")
                        doff = +3
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, lowDia, transHerm, doff)")
                        call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, lowDia, transHerm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, dia, transHerm, doff)")
                        call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, dia, transHerm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("doff = -3")
                        doff = -3
        call disp%show("call setMatCopy(destin(3:7, 3:9), rdpack, source(1:5, 1:7), rdpack, uppLow, doff)")
                        call setMatCopy(destin(3:7, 3:9), rdpack, source(1:5, 1:7), rdpack, uppLow, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:9), rdpack, source(1:5, 1:7), rdpack, dia, doff)")
                        call setMatCopy(destin(3:7, 3:9), rdpack, source(1:5, 1:7), rdpack, dia, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("doff = -3")
                        doff = -3
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, uppLow, transSymm, doff)")
                        call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, uppLow, transSymm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, dia, transSymm, doff)")
                        call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, dia, transSymm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("doff = -3")
                        doff = -3
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, uppLow, transHerm, doff)")
                        call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, uppLow, transHerm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, dia, transHerm, doff)")
                        call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, dia, transHerm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("doff = +3")
                        doff = +3
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:9), rdpack, source(1:5, 1:7), rdpack, uppLow, doff)")
                        call setMatCopy(destin(3:7, 3:9), rdpack, source(1:5, 1:7), rdpack, uppLow, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:7, 3:9), rdpack, source(1:5, 1:7), rdpack, dia, doff)")
                        call setMatCopy(destin(3:7, 3:9), rdpack, source(1:5, 1:7), rdpack, dia, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("doff = +3")
                        doff = +3
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, uppLow, transSymm, doff)")
                        call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, uppLow, transSymm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, dia, transSymm, doff)")
                        call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, dia, transSymm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%show("doff = +3")
                        doff = +3
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, uppLow, transHerm, doff)")
                        call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, uppLow, transHerm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, dia, transHerm, doff)")
                        call setMatCopy(destin(3:9, 3:7), rdpack, source(1:5, 1:7), rdpack, dia, transHerm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

    end block



    block

        complex(CKC) :: source(5,5), destin(5,5)
        complex(CKC), parameter :: NONE = (0._CKC, 0._CKC)

        call disp%skip()
        call disp%show("source = cmplx(getUnifRand(1, 9, size(source,1,IK), size(source,2,IK)), getUnifRand(1, 9, size(source,1,IK), size(source,2,IK)), CKC)")
                        source = cmplx(getUnifRand(1, 9, size(source,1,IK), size(source,2,IK)), getUnifRand(1, 9, size(source,1,IK), size(source,2,IK)), CKC)
        call disp%show("source")
        call disp%show( source , format = cform )

        call disp%show("doff = -1")
                        doff = -1
        call disp%show("destin = NONE")
                        destin = NONE
        call disp%show("call setMatCopy(destin(2:3, 2:4), rdpack, source(1:2, 1:3), rdpack, uppDia, doff)")
                        call setMatCopy(destin(2:3, 2:4), rdpack, source(1:2, 1:3), rdpack, uppDia, doff)
        call disp%show("destin")
        call disp%show( destin , format = cform )
        call disp%show("destin = NONE")
                        destin = NONE
        call disp%show("call setMatCopy(destin(2:3, 2:4), rdpack, source(1:2, 1:3), rdpack, dia, doff)")
                        call setMatCopy(destin(2:3, 2:4), rdpack, source(1:2, 1:3), rdpack, dia, doff)
        call disp%show("destin")
        call disp%show( destin , format = cform )
        call disp%skip()

        call disp%show("doff = -1")
                        doff = -1
        call disp%show("destin = NONE")
                        destin = NONE
        call disp%show("call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, uppDia, transSymm, doff)")
                        call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, uppDia, transSymm, doff)
        call disp%show("destin")
        call disp%show( destin , format = cform )
        call disp%show("destin = NONE")
                        destin = NONE
        call disp%show("call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, dia, transSymm, doff)")
                        call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, dia, transSymm, doff)
        call disp%show("destin")
        call disp%show( destin , format = cform )
        call disp%skip()

        call disp%show("doff = -1")
                        doff = -1
        call disp%show("destin = NONE")
                        destin = NONE
        call disp%show("call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, uppDia, transHerm, doff)")
                        call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, uppDia, transHerm, doff)
        call disp%show("destin")
        call disp%show( destin , format = cform )
        call disp%show("destin = NONE")
                        destin = NONE
        call disp%show("call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, dia, transHerm, doff)")
                        call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, dia, transHerm, doff)
        call disp%show("destin")
        call disp%show( destin , format = cform )
        call disp%skip()

        call disp%show("destin = NONE")
                        destin = NONE
        call disp%show("doff = +1")
                        doff = +1
        call disp%show("call setMatCopy(destin(2:3, 2:4), rdpack, source(1:2, 1:3), rdpack, lowDia, doff)")
                        call setMatCopy(destin(2:3, 2:4), rdpack, source(1:2, 1:3), rdpack, lowDia, doff)
        call disp%show("destin")
        call disp%show( destin , format = cform )
        call disp%skip()

        call disp%show("doff = +1")
                        doff = +1
        call disp%show("destin = NONE")
                        destin = NONE
        call disp%show("call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, lowDia, transSymm, doff)")
                        call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, lowDia, transSymm, doff)
        call disp%show("destin")
        call disp%show( destin , format = cform )
        call disp%show("destin = NONE")
                        destin = NONE
        call disp%show("call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, dia, transSymm, doff)")
                        call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, dia, transSymm, doff)
        call disp%show("destin")
        call disp%show( destin , format = cform )
        call disp%skip()

        call disp%show("doff = +1")
                        doff = +1
        call disp%show("destin = NONE")
                        destin = NONE
        call disp%show("call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, lowDia, transHerm, doff)")
                        call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, lowDia, transHerm, doff)
        call disp%show("destin")
        call disp%show( destin , format = cform )
        call disp%show("destin = NONE")
                        destin = NONE
        call disp%show("call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, dia, transHerm, doff)")
                        call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, dia, transHerm, doff)
        call disp%show("destin")
        call disp%show( destin , format = cform )
        call disp%skip()

        call disp%show("destin = NONE")
                        destin = NONE
        call disp%show("doff = -1")
                        doff = -1
        call disp%show("call setMatCopy(destin(2:3, 2:4), rdpack, source(1:2, 1:3), rdpack, uppLow, doff)")
                        call setMatCopy(destin(2:3, 2:4), rdpack, source(1:2, 1:3), rdpack, uppLow, doff)
        call disp%show("destin")
        call disp%show( destin , format = cform )
        call disp%skip()

        call disp%show("doff = -1")
                        doff = -1
        call disp%show("destin = NONE")
                        destin = NONE
        call disp%show("call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, uppLow, transSymm, doff)")
                        call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, uppLow, transSymm, doff)
        call disp%show("destin")
        call disp%show( destin , format = cform )
        call disp%show("destin = NONE")
                        destin = NONE
        call disp%show("call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, dia, transSymm, doff)")
                        call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, dia, transSymm, doff)
        call disp%show("destin")
        call disp%show( destin , format = cform )
        call disp%skip()

        call disp%show("doff = -1")
                        doff = -1
        call disp%show("destin = NONE")
                        destin = NONE
        call disp%show("call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, uppLow, transHerm, doff)")
                        call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, uppLow, transHerm, doff)
        call disp%show("destin")
        call disp%show( destin , format = cform )
        call disp%show("destin = NONE")
                        destin = NONE
        call disp%show("call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, dia, transHerm, doff)")
                        call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, dia, transHerm, doff)
        call disp%show("destin")
        call disp%show( destin , format = cform )
        call disp%skip()

        call disp%show("doff = +1")
                        doff = +1
        call disp%show("destin = NONE")
                        destin = NONE
        call disp%show("call setMatCopy(destin(2:3, 2:4), rdpack, source(1:2, 1:3), rdpack, uppLow, doff)")
                        call setMatCopy(destin(2:3, 2:4), rdpack, source(1:2, 1:3), rdpack, uppLow, doff)
        call disp%show("destin")
        call disp%show( destin , format = cform )
        call disp%show("destin = NONE")
                        destin = NONE
        call disp%show("call setMatCopy(destin(2:3, 2:4), rdpack, source(1:2, 1:3), rdpack, dia, doff)")
                        call setMatCopy(destin(2:3, 2:4), rdpack, source(1:2, 1:3), rdpack, dia, doff)
        call disp%show("destin")
        call disp%show( destin , format = cform )
        call disp%skip()

        call disp%show("doff = +1")
                        doff = +1
        call disp%show("destin = NONE")
                        destin = NONE
        call disp%show("call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, uppLow, transSymm, doff)")
                        call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, uppLow, transSymm, doff)
        call disp%show("destin")
        call disp%show( destin , format = cform )
        call disp%show("destin = NONE")
                        destin = NONE
        call disp%show("call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, dia, transSymm, doff)")
                        call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, dia, transSymm, doff)
        call disp%show("destin")
        call disp%show( destin , format = cform )
        call disp%skip()

        call disp%show("doff = +1")
                        doff = +1
        call disp%show("destin = NONE")
                        destin = NONE
        call disp%show("call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, uppLow, transHerm, doff)")
                        call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, uppLow, transHerm, doff)
        call disp%show("destin")
        call disp%show( destin , format = cform )
        call disp%show("destin = NONE")
                        destin = NONE
        call disp%show("call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, dia, transHerm, doff)")
                        call setMatCopy(destin(2:4, 2:3), rdpack, source(1:2, 1:3), rdpack, dia, transHerm, doff)
        call disp%show("destin")
        call disp%show( destin , format = cform )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Copy the diagonal elements in Linear Full Packing format.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        integer(IK), parameter :: ndim = 10
        character(2,SKC), parameter :: EMPTY = SKC_"  "
        character(2,SKC) :: source((ndim + 1) * ndim / 2), desnew((ndim + 1) * ndim / 2), destin(ndim, ndim)

        call disp%skip()
        call disp%show("call setUnifRand(source, SKC_'AA', SKC_'ZZ')")
                        call setUnifRand(source, SKC_'AA', SKC_'ZZ')
        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )

        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:9,3:7), rdpack, source(2:6), lfpack, dia)")
                        call setMatCopy(destin(3:9,3:7), rdpack, source(2:6), lfpack, dia)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("desnew = EMPTY")
                        desnew = EMPTY
        call disp%show("call setMatCopy(desnew(2:6), lfpack, destin(3:9,3:7), rdpack, dia)")
                        call setMatCopy(desnew(2:6), lfpack, destin(3:9,3:7), rdpack, dia)
        call disp%show("desnew")
        call disp%show( desnew , deliml = SK_"""" )
        call disp%skip()

        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("doff = -3")
                        doff = -3
        call disp%show("call setMatCopy(destin(3:9,3:7), rdpack, source(2:5), lfpack, dia, doff)")
                        call setMatCopy(destin(3:9,3:7), rdpack, source(2:5), lfpack, dia, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("desnew = EMPTY")
                        desnew = EMPTY
        call disp%show("call setMatCopy(desnew(2:5), lfpack, destin(3:9,3:7), rdpack, dia, doff)")
                        call setMatCopy(desnew(2:5), lfpack, destin(3:9,3:7), rdpack, dia, doff)
        call disp%show("desnew")
        call disp%show( desnew , deliml = SK_"""" )
        call disp%skip()


        call disp%skip()
        call disp%show("call setUnifRand(source, SKC_'AA', SKC_'ZZ')")
                        call setUnifRand(source, SKC_'AA', SKC_'ZZ')
        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )

        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("doff = +3;")
                        doff = +3;
        call disp%show("call setMatCopy(destin(3:7, 3:9), rdpack, source(2:5), lfpack, dia, doff)")
                        call setMatCopy(destin(3:7, 3:9), rdpack, source(2:5), lfpack, dia, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("desnew = EMPTY")
                        desnew = EMPTY
        call disp%show("call setMatCopy(desnew(2:5), lfpack, destin(3:7, 3:9), rdpack, dia, doff)")
                        call setMatCopy(desnew(2:5), lfpack, destin(3:7, 3:9), rdpack, dia, doff)
        call disp%show("desnew")
        call disp%show( desnew , deliml = SK_"""" )
        call disp%skip()


        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("doff = +3;")
                        doff = +3;
        call disp%show("call setMatCopy(destin(3:9, 3:7), rdpack, source(2:3), lfpack, dia, doff)")
                        call setMatCopy(destin(3:9, 3:7), rdpack, source(2:3), lfpack, dia, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("desnew = EMPTY")
                        desnew = EMPTY
        call disp%show("call setMatCopy(desnew(2:3), lfpack, destin(3:9,3:7), rdpack, dia, doff)")
                        call setMatCopy(desnew(2:3), lfpack, destin(3:9,3:7), rdpack, dia, doff)
        call disp%show("desnew")
        call disp%show( desnew , deliml = SK_"""" )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Copy upper/lower triangle in Linear Full Packing format.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        integer(IK), parameter :: ndim = 10
        character(2,SKC), parameter :: EMPTY = SKC_"  "
        character(2,SKC) :: source((ndim + 1) * ndim / 2), desnew((ndim + 1) * ndim / 2), destin(ndim, ndim)

        call disp%skip()
        call disp%show("call setUnifRand(source, SKC_'AA', SKC_'ZZ')")
                        call setUnifRand(source, SKC_'AA', SKC_'ZZ')
        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("doff = -3")
                        doff = -3
        call disp%show("call setMatCopy(destin(3:9,3:7), rdpack, source, lfpack, uppDia, doff)")
                        call setMatCopy(destin(3:9,3:7), rdpack, source, lfpack, uppDia, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("desnew = EMPTY")
                        desnew = EMPTY
        call disp%show("call setMatCopy(desnew, lfpack, destin(3:9,3:7), rdpack, uppDia, doff)")
                        call setMatCopy(desnew, lfpack, destin(3:9,3:7), rdpack, uppDia, doff)
        call disp%show("desnew")
        call disp%show( desnew , deliml = SK_"""" )
        call disp%skip()

        call disp%skip()
        call disp%show("call setUnifRand(source, SKC_'AA', SKC_'ZZ')")
                        call setUnifRand(source, SKC_'AA', SKC_'ZZ')
        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("doff = -3")
                        doff = -3
        call disp%show("call setMatCopy(destin(3:9,3:7), rdpack, source, lfpack, uppDia, transSymm, doff)")
                        call setMatCopy(destin(3:9,3:7), rdpack, source, lfpack, uppDia, transSymm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("desnew = EMPTY")
                        desnew = EMPTY
        call disp%show("call setMatCopy(desnew, lfpack, destin(3:9,3:7), rdpack, lowDia, transSymm, -doff)")
                        call setMatCopy(desnew, lfpack, destin(3:9,3:7), rdpack, lowDia, transSymm, -doff)
        call disp%show("desnew")
        call disp%show( desnew , deliml = SK_"""" )
        call disp%skip()


        call disp%skip()
        call disp%show("call setUnifRand(source, SKC_'AA', SKC_'ZZ')")
                        call setUnifRand(source, SKC_'AA', SKC_'ZZ')
        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )

        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("doff = +3;")
                        doff = +3;
        call disp%show("call setMatCopy(destin(3:7, 3:9), rdpack, source, lfpack, lowDia, doff)")
                        call setMatCopy(destin(3:7, 3:9), rdpack, source, lfpack, lowDia, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(3:9, 3:7), rdpack, source, lfpack, lowDia, transSymm, doff)")
                        call setMatCopy(destin(3:9, 3:7), rdpack, source, lfpack, lowDia, transSymm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("desnew = EMPTY")
                        desnew = EMPTY
        call disp%show("call setMatCopy(desnew, lfpack, destin(3:9,3:7), rdpack, lowDia, transSymm, doff)")
                        call setMatCopy(desnew, lfpack, destin(3:9,3:7), rdpack, lowDia, transSymm, doff)
        call disp%show("desnew")
        call disp%show( desnew , deliml = SK_"""" )
        call disp%skip()


        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("doff = +3;")
                        doff = +3;
        call disp%show("call setMatCopy(destin(3:9, 3:7), rdpack, source, lfpack, lowDia, transHerm, doff)")
                        call setMatCopy(destin(3:9, 3:7), rdpack, source, lfpack, lowDia, transHerm, doff)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("desnew = EMPTY")
                        desnew = EMPTY
        call disp%show("call setMatCopy(desnew, lfpack, destin(3:9,3:7), rdpack, lowDia, transHerm, doff)")
                        call setMatCopy(desnew, lfpack, destin(3:9,3:7), rdpack, lowDia, transHerm, doff)
        call disp%show("desnew")
        call disp%show( desnew , deliml = SK_"""" )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Copy upper/lower triangle in Rectangular Full Packing format.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        character(2,SKC) :: source(10,10), destin(10,10), desnew(10,10)
        character(2,SKC), parameter :: EMPTY = SKC_"  "

        call disp%skip()
        call disp%show("call setUnifRand(source, SKC_'AA', SKC_'ZZ')")
                        call setUnifRand(source, SKC_'AA', SKC_'ZZ')
        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(4:8, 3:5), rfpack, source(2:6, 3:7), rdpack, uppDia)")
                        call setMatCopy(destin(4:8, 3:5), rfpack, source(2:6, 3:7), rdpack, uppDia)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("desnew = EMPTY")
                        desnew = EMPTY
        call disp%show("call setMatCopy(desnew(3:5, 4:8), rfpack, source(2:6, 3:7), rdpack, uppDia) ! transpose")
                        call setMatCopy(desnew(3:5, 4:8), rfpack, source(2:6, 3:7), rdpack, uppDia) ! transpose
        call disp%show("desnew")
        call disp%show( desnew , deliml = SK_"""" )
        call disp%show("desnew = EMPTY")
                        desnew = EMPTY
        call disp%show("call setMatCopy(desnew(2:6, 3:7), rdpack, destin(4:8, 3:5), rfpack, uppDia)")
                        call setMatCopy(desnew(2:6, 3:7), rdpack, destin(4:8, 3:5), rfpack, uppDia)
        call disp%show("desnew")
        call disp%show( desnew , deliml = SK_"""" )
        call disp%skip()

        call disp%skip()
        call disp%show("call setUnifRand(source, SKC_'AA', SKC_'ZZ')")
                        call setUnifRand(source, SKC_'AA', SKC_'ZZ')
        call disp%show("source")
        call disp%show( source , deliml = SK_"""" )
        call disp%show("destin = EMPTY")
                        destin = EMPTY
        call disp%show("call setMatCopy(destin(4:8, 3:5), rfpack, source(2:6, 3:7), rdpack, lowDia)")
                        call setMatCopy(destin(4:8, 3:5), rfpack, source(2:6, 3:7), rdpack, lowDia)
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%show("desnew = EMPTY")
                        desnew = EMPTY
        call disp%show("call setMatCopy(desnew(2:6, 3:7), rdpack, destin(4:8, 3:5), rfpack, lowDia)")
                        call setMatCopy(desnew(2:6, 3:7), rdpack, destin(4:8, 3:5), rfpack, lowDia)
        call disp%show("desnew")
        call disp%show( desnew , deliml = SK_"""" )
        call disp%skip()

    end block

end program example