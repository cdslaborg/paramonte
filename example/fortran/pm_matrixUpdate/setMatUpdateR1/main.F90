program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKG => RKS ! all processor type kinds are supported.
    use pm_kind, only: CKG => CKS ! all processor type kinds are supported.
    use pm_io, only: display_type
    use pm_matrixUpdate, only: transHerm
    use pm_matrixUpdate, only: setMatUpdateR1

    implicit none

    type(display_type) :: disp

    real(RKG)   , parameter     :: dummyr = -huge(0._RKG)
    complex(CKG), parameter     :: dumm_cmplx_value = cmplx(-huge(0._CKG), -huge(0._CKG), CKG)

    disp = display_type(file = "main.out.F90")

    block

        real(RKG), allocatable :: RefA(:,:), mat(:,:), vecA(:), vecB(:)

        mat = reshape( [ 1._RKG, 2._RKG, 3._RKG &
                        , 2._RKG, 2._RKG, 4._RKG &
                        , 3._RKG, 2._RKG, 2._RKG &
                        , 4._RKG, 2._RKG, 1._RKG &
                        ], shape = [4, 3], order = [2, 1])
        RefA = reshape( [ 4._RKG, 8._RKG,12._RKG &
                        , 4._RKG, 6._RKG,10._RKG &
                        , 4._RKG, 4._RKG, 5._RKG &
                        , 8._RKG,10._RKG,13._RKG &
                        ], shape = [4, 3], order = [2, 1])
        call disp%skip()
        call disp%show("dummyr ! some dummy value to illustrate functionality of `incA` and `incB` arguments.")
        call disp%show( dummyr )
        call disp%show("mat")
        call disp%show( mat )
        call disp%show("vecA = [3._RKG, 2._RKG, 1._RKG, 4._RKG]")
                        vecA = [3._RKG, 2._RKG, 1._RKG, 4._RKG]
        call disp%show("vecB = [1._RKG, dummyr, 2._RKG, dummyr, 3._RKG]")
                        vecB = [1._RKG, dummyr, 2._RKG, dummyr, 3._RKG]
        call disp%show("call setMatUpdateR1(mat, vecA, vecB, incA = 1_IK, incB = 2_IK, roff = 0_IK)")
                        call setMatUpdateR1(mat, vecA, vecB, incA = 1_IK, incB = 2_IK, roff = 0_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%show("mat - RefA")
        call disp%show( mat - RefA )
        call disp%skip()

        mat = reshape( [ 1._RKG, 2._RKG, 3._RKG &
                        , 2._RKG, 2._RKG, 4._RKG &
                        , 3._RKG, 2._RKG, 2._RKG &
                        , 4._RKG, 2._RKG, 1._RKG &
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
        call disp%show("vecA = [3._RKG, 2._RKG, 1._RKG, 4._RKG]")
                        vecA = [3._RKG, 2._RKG, 1._RKG, 4._RKG]
        call disp%show("vecB = [3._RKG, dummyr, 2._RKG, dummyr, 1._RKG]")
                        vecB = [3._RKG, dummyr, 2._RKG, dummyr, 1._RKG]
        call disp%show("call setMatUpdateR1(mat, vecA, vecB, incA = 1_IK, incB = -2_IK, roff = 0_IK)")
                        call setMatUpdateR1(mat, vecA, vecB, incA = 1_IK, incB = -2_IK, roff = 0_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

        mat = reshape( [ dummyr, dummyr, dummyr &
                        , dummyr, dummyr, dummyr &
                        , dummyr, dummyr, dummyr &
                        , 1._RKG, 2._RKG, 3._RKG &
                        , 2._RKG, 2._RKG, 4._RKG &
                        , 3._RKG, 2._RKG, 2._RKG &
                        , 4._RKG, 2._RKG, 1._RKG &
                        , dummyr, dummyr, dummyr &
                        , dummyr, dummyr, dummyr &
                        , dummyr, dummyr, dummyr &
                        ], shape = [10, 3], order = [2, 1])
        call disp%skip()
        call disp%show("dummyr ! some dummy value to illustrate functionality of `incA` and `incB` arguments.")
        call disp%show( dummyr )
        call disp%show("mat")
        call disp%show( mat )
        call disp%show("vecA = [3._RKG, dummyr, dummyr, 2._RKG, dummyr, dummyr, 1._RKG, dummyr, dummyr, 4._RKG]")
                        vecA = [3._RKG, dummyr, dummyr, 2._RKG, dummyr, dummyr, 1._RKG, dummyr, dummyr, 4._RKG]
        call disp%show("vecB = [1._RKG, 2._RKG, 3._RKG]")
                        vecB = [1._RKG, 2._RKG, 3._RKG]
        call disp%show("call setMatUpdateR1(mat, vecA, vecB, incA = 3_IK, incB = 1_IK, roff = 3_IK)")
                        call setMatUpdateR1(mat, vecA, vecB, incA = 3_IK, incB = 1_IK, roff = 3_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

        mat = reshape( [ 1._RKG, 2._RKG, 3._RKG &
                        , 2._RKG, 2._RKG, 4._RKG &
                        , 3._RKG, 2._RKG, 2._RKG &
                        , 4._RKG, 2._RKG, 1._RKG &
                        ], shape = [4, 3], order = [2, 1])
        RefA = reshape( [ 7._RKG,14._RKG,21._RKG &
                        , 6._RKG,10._RKG,16._RKG &
                        , 5._RKG, 6._RKG, 8._RKG &
                        ,12._RKG,18._RKG,25._RKG &
                        ], shape = [4, 3], order = [2, 1])
        call disp%skip()
        call disp%show("dummyr ! some dummy value to illustrate functionality of `incA` and `incB` arguments.")
        call disp%show( dummyr )
        call disp%show("mat")
        call disp%show( mat )
        call disp%show("vecA = [3._RKG, 2._RKG, 1._RKG, 4._RKG]")
                        vecA = [3._RKG, 2._RKG, 1._RKG, 4._RKG]
        call disp%show("vecB = [1._RKG, dummyr, 2._RKG, dummyr, 3._RKG]")
                        vecB = [1._RKG, dummyr, 2._RKG, dummyr, 3._RKG]
        call disp%show("call setMatUpdateR1(mat, vecA, vecB, alpha = 2._RKG, incA = 1_IK, incB = 2_IK, roff = 0_IK)")
                        call setMatUpdateR1(mat, vecA, vecB, alpha = 2._RKG, incA = 1_IK, incB = 2_IK, roff = 0_IK)
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

        complex(CKG), allocatable   :: RefA(:,:), mat(:,:), vecA(:), vecB(:)

        mat = reshape( [ (1._CKG, 2._CKG), (3._CKG, 5._CKG), (2._CKG, 0._CKG) &
                        , (2._CKG, 3._CKG), (7._CKG, 9._CKG), (4._CKG, 8._CKG) &
                        , (7._CKG, 4._CKG), (1._CKG, 4._CKG), (6._CKG, 0._CKG) &
                        , (8._CKG, 2._CKG), (2._CKG, 5._CKG), (8._CKG, 0._CKG) &
                        , (9._CKG, 1._CKG), (3._CKG, 6._CKG), (1._CKG, 0._CKG) &
                        , dumm_cmplx_value, dumm_cmplx_value, dumm_cmplx_value &
                        , dumm_cmplx_value, dumm_cmplx_value, dumm_cmplx_value &
                        , dumm_cmplx_value, dumm_cmplx_value, dumm_cmplx_value &
                        , dumm_cmplx_value, dumm_cmplx_value, dumm_cmplx_value &
                        , dumm_cmplx_value, dumm_cmplx_value, dumm_cmplx_value &
                        ], shape = [10, 3], order = [2, 1])
        RefA = reshape( [ (-2._CKG, 6._CKG), ( 7._CKG,13._CKG), ( 5._CKG, 1._CKG) &
                        , ( 6._CKG,11._CKG), (23._CKG, 9._CKG), ( 8._CKG, 4._CKG) &
                        , ( 6._CKG, 7._CKG), ( 5._CKG, 8._CKG), ( 8._CKG, 0._CKG) &
                        , ( 3._CKG,12._CKG), (14._CKG,21._CKG), (15._CKG, 1._CKG) &
                        , (11._CKG, 5._CKG), (11._CKG, 6._CKG), ( 3._CKG,-2._CKG) &
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
        call disp%show("vecA = [(1._CKG, 2._CKG), (4._CKG, 0._CKG), (1._CKG, 1._CKG), (3._CKG, 4._CKG), (2._CKG, 0._CKG)]")
                        vecA = [(1._CKG, 2._CKG), (4._CKG, 0._CKG), (1._CKG, 1._CKG), (3._CKG, 4._CKG), (2._CKG, 0._CKG)]
        call disp%show("vecB = [(1._CKG, 2._CKG), (4._CKG, 0._CKG), (1._CKG, -1._CKG)]")
                        vecB = [(1._CKG, 2._CKG), (4._CKG, 0._CKG), (1._CKG, -1._CKG)]
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


        mat = reshape( [ (1._CKG, 2._CKG), (3._CKG, 5._CKG), (2._CKG, 0._CKG) &
                        , (2._CKG, 3._CKG), (7._CKG, 9._CKG), (4._CKG, 8._CKG) &
                        , (7._CKG, 4._CKG), (1._CKG, 4._CKG), (6._CKG, 0._CKG) &
                        , (8._CKG, 2._CKG), (2._CKG, 5._CKG), (8._CKG, 0._CKG) &
                        , (9._CKG, 1._CKG), (3._CKG, 6._CKG), (1._CKG, 0._CKG) &
                        , dumm_cmplx_value, dumm_cmplx_value, dumm_cmplx_value &
                        , dumm_cmplx_value, dumm_cmplx_value, dumm_cmplx_value &
                        , dumm_cmplx_value, dumm_cmplx_value, dumm_cmplx_value &
                        , dumm_cmplx_value, dumm_cmplx_value, dumm_cmplx_value &
                        , dumm_cmplx_value, dumm_cmplx_value, dumm_cmplx_value &
                        ], shape = [10, 3], order = [2, 1])
        RefA = reshape( [ ( 6._CKG, 2._CKG), ( 7._CKG,13._CKG), (1._CKG,  3._CKG) &
                        , ( 6._CKG,-5._CKG), (23._CKG, 9._CKG), (8._CKG, 12._CKG) &
                        , (10._CKG, 3._CKG), ( 5._CKG, 8._CKG), (6._CKG,  2._CKG) &
                        , (19._CKG, 0._CKG), (14._CKG,21._CKG), (7._CKG,  7._CKG) &
                        , (11._CKG,-3._CKG), (11._CKG, 6._CKG), (3._CKG,  2._CKG) &
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
        call disp%show("vecA = [(1._CKG, 2._CKG), (4._CKG, 0._CKG), (1._CKG, 1._CKG), (3._CKG, 4._CKG), (2._CKG, 0._CKG)]")
                        vecA = [(1._CKG, 2._CKG), (4._CKG, 0._CKG), (1._CKG, 1._CKG), (3._CKG, 4._CKG), (2._CKG, 0._CKG)]
        call disp%show("vecB = [(1._CKG, 2._CKG), (4._CKG, 0._CKG), (1._CKG, -1._CKG)]")
                        vecB = [(1._CKG, 2._CKG), (4._CKG, 0._CKG), (1._CKG, -1._CKG)]
        call disp%show("call setMatUpdateR1(mat, vecA, vecB, operationB = transHerm, incA = 1_IK, incB = 1_IK, roff = 0_IK)")
                        call setMatUpdateR1(mat, vecA, vecB, operationB = transHerm, incA = 1_IK, incB = 1_IK, roff = 0_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

    end block

end program example