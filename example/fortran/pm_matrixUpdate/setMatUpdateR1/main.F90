program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKC => RK32 ! all processor type kinds are supported.
    use pm_kind, only: CKC => CK32 ! all processor type kinds are supported.
    use pm_io, only: display_type
    use pm_matrixUpdate, only: transHerm
    use pm_matrixUpdate, only: setMatUpdateR1

    implicit none

    type(display_type) :: disp

    real(RKC)   , parameter     :: dummyr = -huge(0._RKC)
    complex(CKC), parameter     :: dumm_cmplx_value = cmplx(-huge(0._CKC), -huge(0._CKC), CKC)

    disp = display_type(file = "main.out.F90")

    block

        real(RKC), allocatable :: RefA(:,:), mat(:,:), vecA(:), vecB(:)

        mat = reshape( [ 1._RKC, 2._RKC, 3._RKC &
                        , 2._RKC, 2._RKC, 4._RKC &
                        , 3._RKC, 2._RKC, 2._RKC &
                        , 4._RKC, 2._RKC, 1._RKC &
                        ], shape = [4, 3], order = [2, 1])
        RefA = reshape( [ 4._RKC, 8._RKC,12._RKC &
                        , 4._RKC, 6._RKC,10._RKC &
                        , 4._RKC, 4._RKC, 5._RKC &
                        , 8._RKC,10._RKC,13._RKC &
                        ], shape = [4, 3], order = [2, 1])
        call disp%skip()
        call disp%show("dummyr ! some dummy value to illustrate functionality of `incA` and `incB` arguments.")
        call disp%show( dummyr )
        call disp%show("mat")
        call disp%show( mat )
        call disp%show("vecA = [3._RKC, 2._RKC, 1._RKC, 4._RKC]")
                        vecA = [3._RKC, 2._RKC, 1._RKC, 4._RKC]
        call disp%show("vecB = [1._RKC, dummyr, 2._RKC, dummyr, 3._RKC]")
                        vecB = [1._RKC, dummyr, 2._RKC, dummyr, 3._RKC]
        call disp%show("call setMatUpdateR1(mat, vecA, vecB, incA = 1_IK, incB = 2_IK, roff = 0_IK)")
                        call setMatUpdateR1(mat, vecA, vecB, incA = 1_IK, incB = 2_IK, roff = 0_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%show("mat - RefA")
        call disp%show( mat - RefA )
        call disp%skip()

        mat = reshape( [ 1._RKC, 2._RKC, 3._RKC &
                        , 2._RKC, 2._RKC, 4._RKC &
                        , 3._RKC, 2._RKC, 2._RKC &
                        , 4._RKC, 2._RKC, 1._RKC &
                        , dummyr, dummyr, dummyr &
                        , dummyr, dummyr, dummyr &
                        , dummyr, dummyr, dummyr &
                        , dummyr, dummyr, dummyr &
                        , dummyr, dummyr, dummyr &
                        , dummyr, dummyr, dummyr &
                        ], shape = [10, 3], order = [2, 1])
        call disp%skip()
        call disp%show("dummyr ! some dummy value to illustrate functionality of `incA` and `incB` arguments.")
        call disp%show( dummyr )
        call disp%show("mat")
        call disp%show( mat )
        call disp%show("vecA = [3._RKC, 2._RKC, 1._RKC, 4._RKC]")
                        vecA = [3._RKC, 2._RKC, 1._RKC, 4._RKC]
        call disp%show("vecB = [3._RKC, dummyr, 2._RKC, dummyr, 1._RKC]")
                        vecB = [3._RKC, dummyr, 2._RKC, dummyr, 1._RKC]
        call disp%show("call setMatUpdateR1(mat, vecA, vecB, incA = 1_IK, incB = -2_IK, roff = 0_IK)")
                        call setMatUpdateR1(mat, vecA, vecB, incA = 1_IK, incB = -2_IK, roff = 0_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

        mat = reshape( [ dummyr, dummyr, dummyr &
                        , dummyr, dummyr, dummyr &
                        , dummyr, dummyr, dummyr &
                        , 1._RKC, 2._RKC, 3._RKC &
                        , 2._RKC, 2._RKC, 4._RKC &
                        , 3._RKC, 2._RKC, 2._RKC &
                        , 4._RKC, 2._RKC, 1._RKC &
                        , dummyr, dummyr, dummyr &
                        , dummyr, dummyr, dummyr &
                        , dummyr, dummyr, dummyr &
                        ], shape = [10, 3], order = [2, 1])
        call disp%skip()
        call disp%show("dummyr ! some dummy value to illustrate functionality of `incA` and `incB` arguments.")
        call disp%show( dummyr )
        call disp%show("mat")
        call disp%show( mat )
        call disp%show("vecA = [3._RKC, dummyr, dummyr, 2._RKC, dummyr, dummyr, 1._RKC, dummyr, dummyr, 4._RKC]")
                        vecA = [3._RKC, dummyr, dummyr, 2._RKC, dummyr, dummyr, 1._RKC, dummyr, dummyr, 4._RKC]
        call disp%show("vecB = [1._RKC, 2._RKC, 3._RKC]")
                        vecB = [1._RKC, 2._RKC, 3._RKC]
        call disp%show("call setMatUpdateR1(mat, vecA, vecB, incA = 3_IK, incB = 1_IK, roff = 3_IK)")
                        call setMatUpdateR1(mat, vecA, vecB, incA = 3_IK, incB = 1_IK, roff = 3_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

        mat = reshape( [ 1._RKC, 2._RKC, 3._RKC &
                        , 2._RKC, 2._RKC, 4._RKC &
                        , 3._RKC, 2._RKC, 2._RKC &
                        , 4._RKC, 2._RKC, 1._RKC &
                        ], shape = [4, 3], order = [2, 1])
        RefA = reshape( [ 7._RKC,14._RKC,21._RKC &
                        , 6._RKC,10._RKC,16._RKC &
                        , 5._RKC, 6._RKC, 8._RKC &
                        ,12._RKC,18._RKC,25._RKC &
                        ], shape = [4, 3], order = [2, 1])
        call disp%skip()
        call disp%show("dummyr ! some dummy value to illustrate functionality of `incA` and `incB` arguments.")
        call disp%show( dummyr )
        call disp%show("mat")
        call disp%show( mat )
        call disp%show("vecA = [3._RKC, 2._RKC, 1._RKC, 4._RKC]")
                        vecA = [3._RKC, 2._RKC, 1._RKC, 4._RKC]
        call disp%show("vecB = [1._RKC, dummyr, 2._RKC, dummyr, 3._RKC]")
                        vecB = [1._RKC, dummyr, 2._RKC, dummyr, 3._RKC]
        call disp%show("call setMatUpdateR1(mat, vecA, vecB, alpha = 2._RKC, incA = 1_IK, incB = 2_IK, roff = 0_IK)")
                        call setMatUpdateR1(mat, vecA, vecB, alpha = 2._RKC, incA = 1_IK, incB = 2_IK, roff = 0_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%show("mat - RefA")
        call disp%show( mat - RefA )
        call disp%skip()

    end block

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Complex matrix update of rank one.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    block

        complex(CKC), allocatable   :: RefA(:,:), mat(:,:), vecA(:), vecB(:)

        mat = reshape( [ (1._CKC, 2._CKC), (3._CKC, 5._CKC), (2._CKC, 0._CKC) &
                        , (2._CKC, 3._CKC), (7._CKC, 9._CKC), (4._CKC, 8._CKC) &
                        , (7._CKC, 4._CKC), (1._CKC, 4._CKC), (6._CKC, 0._CKC) &
                        , (8._CKC, 2._CKC), (2._CKC, 5._CKC), (8._CKC, 0._CKC) &
                        , (9._CKC, 1._CKC), (3._CKC, 6._CKC), (1._CKC, 0._CKC) &
                        , dumm_cmplx_value, dumm_cmplx_value, dumm_cmplx_value &
                        , dumm_cmplx_value, dumm_cmplx_value, dumm_cmplx_value &
                        , dumm_cmplx_value, dumm_cmplx_value, dumm_cmplx_value &
                        , dumm_cmplx_value, dumm_cmplx_value, dumm_cmplx_value &
                        , dumm_cmplx_value, dumm_cmplx_value, dumm_cmplx_value &
                        ], shape = [10, 3], order = [2, 1])
        RefA = reshape( [ (-2._CKC, 6._CKC), ( 7._CKC,13._CKC), ( 5._CKC, 1._CKC) &
                        , ( 6._CKC,11._CKC), (23._CKC, 9._CKC), ( 8._CKC, 4._CKC) &
                        , ( 6._CKC, 7._CKC), ( 5._CKC, 8._CKC), ( 8._CKC, 0._CKC) &
                        , ( 3._CKC,12._CKC), (14._CKC,21._CKC), (15._CKC, 1._CKC) &
                        , (11._CKC, 5._CKC), (11._CKC, 6._CKC), ( 3._CKC,-2._CKC) &
                        ,  dumm_cmplx_value,  dumm_cmplx_value,  dumm_cmplx_value &
                        ,  dumm_cmplx_value,  dumm_cmplx_value,  dumm_cmplx_value &
                        ,  dumm_cmplx_value,  dumm_cmplx_value,  dumm_cmplx_value &
                        ,  dumm_cmplx_value,  dumm_cmplx_value,  dumm_cmplx_value &
                        ,  dumm_cmplx_value,  dumm_cmplx_value,  dumm_cmplx_value &
                        ], shape = [10, 3], order = [2, 1])

        call disp%skip()
        call disp%show("dummyr ! some dummy value to illustrate functionality of `incA` and `incB` arguments.")
        call disp%show( dummyr )
        call disp%show("mat")
        call disp%show( mat )
        call disp%show("vecA = [(1._CKC, 2._CKC), (4._CKC, 0._CKC), (1._CKC, 1._CKC), (3._CKC, 4._CKC), (2._CKC, 0._CKC)]")
                        vecA = [(1._CKC, 2._CKC), (4._CKC, 0._CKC), (1._CKC, 1._CKC), (3._CKC, 4._CKC), (2._CKC, 0._CKC)]
        call disp%show("vecB = [(1._CKC, 2._CKC), (4._CKC, 0._CKC), (1._CKC, -1._CKC)]")
                        vecB = [(1._CKC, 2._CKC), (4._CKC, 0._CKC), (1._CKC, -1._CKC)]
        call disp%show("call setMatUpdateR1(mat, vecA, vecB, incA = 1_IK, incB = 1_IK, roff = 0_IK)")
                        call setMatUpdateR1(mat, vecA, vecB, incA = 1_IK, incB = 1_IK, roff = 0_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

        call disp%skip
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!Complex matrix update of rank one with complex conjugate transpose of `Y`.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip


        mat = reshape( [ (1._CKC, 2._CKC), (3._CKC, 5._CKC), (2._CKC, 0._CKC) &
                        , (2._CKC, 3._CKC), (7._CKC, 9._CKC), (4._CKC, 8._CKC) &
                        , (7._CKC, 4._CKC), (1._CKC, 4._CKC), (6._CKC, 0._CKC) &
                        , (8._CKC, 2._CKC), (2._CKC, 5._CKC), (8._CKC, 0._CKC) &
                        , (9._CKC, 1._CKC), (3._CKC, 6._CKC), (1._CKC, 0._CKC) &
                        , dumm_cmplx_value, dumm_cmplx_value, dumm_cmplx_value &
                        , dumm_cmplx_value, dumm_cmplx_value, dumm_cmplx_value &
                        , dumm_cmplx_value, dumm_cmplx_value, dumm_cmplx_value &
                        , dumm_cmplx_value, dumm_cmplx_value, dumm_cmplx_value &
                        , dumm_cmplx_value, dumm_cmplx_value, dumm_cmplx_value &
                        ], shape = [10, 3], order = [2, 1])
        RefA = reshape( [ ( 6._CKC, 2._CKC), ( 7._CKC,13._CKC), (1._CKC,  3._CKC) &
                        , ( 6._CKC,-5._CKC), (23._CKC, 9._CKC), (8._CKC, 12._CKC) &
                        , (10._CKC, 3._CKC), ( 5._CKC, 8._CKC), (6._CKC,  2._CKC) &
                        , (19._CKC, 0._CKC), (14._CKC,21._CKC), (7._CKC,  7._CKC) &
                        , (11._CKC,-3._CKC), (11._CKC, 6._CKC), (3._CKC,  2._CKC) &
                        ,  dumm_cmplx_value,  dumm_cmplx_value,  dumm_cmplx_value &
                        ,  dumm_cmplx_value,  dumm_cmplx_value,  dumm_cmplx_value &
                        ,  dumm_cmplx_value,  dumm_cmplx_value,  dumm_cmplx_value &
                        ,  dumm_cmplx_value,  dumm_cmplx_value,  dumm_cmplx_value &
                        ,  dumm_cmplx_value,  dumm_cmplx_value,  dumm_cmplx_value &
                        ], shape = [10, 3], order = [2, 1])

        call disp%skip()
        call disp%show("dummyr ! some dummy value to illustrate functionality of `incA` and `incB` arguments.")
        call disp%show( dummyr )
        call disp%show("mat")
        call disp%show( mat )
        call disp%show("vecA = [(1._CKC, 2._CKC), (4._CKC, 0._CKC), (1._CKC, 1._CKC), (3._CKC, 4._CKC), (2._CKC, 0._CKC)]")
                        vecA = [(1._CKC, 2._CKC), (4._CKC, 0._CKC), (1._CKC, 1._CKC), (3._CKC, 4._CKC), (2._CKC, 0._CKC)]
        call disp%show("vecB = [(1._CKC, 2._CKC), (4._CKC, 0._CKC), (1._CKC, -1._CKC)]")
                        vecB = [(1._CKC, 2._CKC), (4._CKC, 0._CKC), (1._CKC, -1._CKC)]
        call disp%show("call setMatUpdateR1(mat, vecA, vecB, operationB = transHerm, incA = 1_IK, incB = 1_IK, roff = 0_IK)")
                        call setMatUpdateR1(mat, vecA, vecB, operationB = transHerm, incA = 1_IK, incB = 1_IK, roff = 0_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

    end block

end program example