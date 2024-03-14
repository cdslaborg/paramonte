program example

    use pm_kind, only: SK, IK, LK, CK, RK
    use pm_matrixSubset, only: dia, uppDia, lowDia, uppLow, uppLowDia
    use pm_distUnif, only: setUnifRand
    use pm_matrixTrans, only: setMatTrans
    use pm_io, only: display_type

    implicit none

    type(display_type)  :: disp
    disp = display_type(file = "main.out.F90")

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        character(2) :: matA(5,10), matB(10,5)

        call disp%skip()
        call disp%show("call setUnifRand(matA, 'AA', 'ZZ')")
                        call setUnifRand(matA, 'AA', 'ZZ')
        call disp%show("matA")
        call disp%show( matA , deliml = SK_"""" )
        call disp%show("call setMatTrans(matA, matB)")
                        call setMatTrans(matA, matB)
        call disp%show("matB")
        call disp%show( matB , deliml = SK_"""" )
        call disp%skip()

    end block

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        character(2) :: matA(10,10)

        call disp%skip()
        call disp%show("call setUnifRand(matA, 'AA', 'ZZ')")
                        call setUnifRand(matA, 'AA', 'ZZ')
        call disp%show("matA")
        call disp%show( matA , deliml = SK_"""" )
        call disp%show("call setMatTrans(matA)")
                        call setMatTrans(matA)
        call disp%show("matA")
        call disp%show( matA , deliml = SK_"""" )
        call disp%skip()

    end block

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end program example