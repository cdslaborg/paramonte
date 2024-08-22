program example

    use pm_kind, only: IK, LK, SK, RKS, RKD, RKH
    use pm_arrayResize, only: setResized
    use pm_matrixCopy, only: setMatCopy
    use pm_matrixLUP, only: setMatLUP
    use pm_io, only: display_type

    implicit none

    integer(IK), allocatable :: rperm(:), rperm_ref(:)
    integer(IK) :: info

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the LU-Pivoted decomposition of a real square matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    block
        use pm_kind, only: TKG => RKS
        real(TKG), allocatable :: mat_lup(:,:), lup_ref(:,:)
        mat_lup = reshape(  [ 1.0_TKG, +1.2_TKG, 1.4_TKG, 1.6_TKG, 1.8_TKG, 2.0_TKG, 2.2_TKG, 2.4_TKG, 2.6_TKG &
                            , 1.2_TKG, +1.0_TKG, 1.2_TKG, 1.4_TKG, 1.6_TKG, 1.8_TKG, 2.0_TKG, 2.2_TKG, 2.4_TKG &
                            , 1.4_TKG, +1.2_TKG, 1.0_TKG, 1.2_TKG, 1.4_TKG, 1.6_TKG, 1.8_TKG, 2.0_TKG, 2.2_TKG &
                            , 1.6_TKG, +1.4_TKG, 1.2_TKG, 1.0_TKG, 1.2_TKG, 1.4_TKG, 1.6_TKG, 1.8_TKG, 2.0_TKG &
                            , 1.8_TKG, +1.6_TKG, 1.4_TKG, 1.2_TKG, 1.0_TKG, 1.2_TKG, 1.4_TKG, 1.6_TKG, 1.8_TKG &
                            , 2.0_TKG, +1.8_TKG, 1.6_TKG, 1.4_TKG, 1.2_TKG, 1.0_TKG, 1.2_TKG, 1.4_TKG, 1.6_TKG &
                            , 2.2_TKG, +2.0_TKG, 1.8_TKG, 1.6_TKG, 1.4_TKG, 1.2_TKG, 1.0_TKG, 1.2_TKG, 1.4_TKG &
                            , 2.4_TKG, +2.2_TKG, 2.0_TKG, 1.8_TKG, 1.6_TKG, 1.4_TKG, 1.2_TKG, 1.0_TKG, 1.2_TKG &
                            , 2.6_TKG, +2.4_TKG, 2.2_TKG, 2.0_TKG, 1.8_TKG, 1.6_TKG, 1.4_TKG, 1.2_TKG, 1.0_TKG &
                            ], shape = [9, 9], order = [2, 1])
        lup_ref = reshape(  [ 2.6_TKG,  2.4_TKG, 2.2_TKG, 2.0_TKG, 1.8_TKG, 1.6_TKG, 1.4_TKG, 1.2_TKG, 1.0_TKG &
                            , 0.4_TKG,  0.3_TKG, 0.6_TKG, 0.8_TKG, 1.1_TKG, 1.4_TKG, 1.7_TKG, 1.9_TKG, 2.2_TKG &
                            , 0.5_TKG, -0.4_TKG, 0.4_TKG, 0.8_TKG, 1.2_TKG, 1.6_TKG, 2.0_TKG, 2.4_TKG, 2.8_TKG &
                            , 0.5_TKG, -0.3_TKG, 0.0_TKG, 0.4_TKG, 0.8_TKG, 1.2_TKG, 1.6_TKG, 2.0_TKG, 2.4_TKG &
                            , 0.6_TKG, -0.3_TKG, 0.0_TKG, 0.0_TKG, 0.4_TKG, 0.8_TKG, 1.2_TKG, 1.6_TKG, 2.0_TKG &
                            , 0.7_TKG, -0.2_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.4_TKG, 0.8_TKG, 1.2_TKG, 1.6_TKG &
                            , 0.8_TKG, -0.2_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.4_TKG, 0.8_TKG, 1.2_TKG &
                            , 0.8_TKG, -0.1_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.4_TKG, 0.8_TKG &
                            , 0.9_TKG, -0.1_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.4_TKG &
                            ], shape = [9, 9], order = [2, 1])
        rperm_ref = [9, 9, 9, 9, 9, 9, 9, 9, 9]
        call disp%skip
        call disp%show("mat_lup")
        call disp%show( mat_lup )
        call disp%show("call setResized(rperm, size(mat_lup, 1, IK))")
                        call setResized(rperm, size(mat_lup, 1, IK))
        call disp%show("call setMatLUP(mat_lup, rperm, info); if (info /= 0) error stop")
                        call setMatLUP(mat_lup, rperm, info); if (info /= 0) error stop
        call disp%show("mat_lup")
        call disp%show( mat_lup )
        call disp%show("call setResized(rperm, size(mat_lup, 1, IK)) ! reconstruct the original matrix.")
                        call setResized(rperm, size(mat_lup, 1, IK))
        call disp%show("lup_ref ! reference matrix rounded to 1 significant digit.")
        call disp%show( lup_ref )
        call disp%show("lup_ref - mat_lup, format = SK_'(*(f0.1,:,"", ""))'")
        call disp%show( lup_ref - mat_lup, format = SK_"(*(f0.1,:,"", ""))")
        call disp%show("rperm_ref")
        call disp%show( rperm_ref )
        call disp%show("rperm")
        call disp%show( rperm )
        call disp%skip
    end block

    block
        use pm_kind, only: TKG => RKD
        real(TKG), allocatable :: mat_lup(:,:), lup_ref(:,:)
        mat_lup = reshape(  [ 1.0_TKG, +1.2_TKG, 1.4_TKG, 1.6_TKG, 1.8_TKG, 2.0_TKG, 2.2_TKG, 2.4_TKG, 2.6_TKG &
                            , 1.2_TKG, +1.0_TKG, 1.2_TKG, 1.4_TKG, 1.6_TKG, 1.8_TKG, 2.0_TKG, 2.2_TKG, 2.4_TKG &
                            , 1.4_TKG, +1.2_TKG, 1.0_TKG, 1.2_TKG, 1.4_TKG, 1.6_TKG, 1.8_TKG, 2.0_TKG, 2.2_TKG &
                            , 1.6_TKG, +1.4_TKG, 1.2_TKG, 1.0_TKG, 1.2_TKG, 1.4_TKG, 1.6_TKG, 1.8_TKG, 2.0_TKG &
                            , 1.8_TKG, +1.6_TKG, 1.4_TKG, 1.2_TKG, 1.0_TKG, 1.2_TKG, 1.4_TKG, 1.6_TKG, 1.8_TKG &
                            , 2.0_TKG, +1.8_TKG, 1.6_TKG, 1.4_TKG, 1.2_TKG, 1.0_TKG, 1.2_TKG, 1.4_TKG, 1.6_TKG &
                            , 2.2_TKG, +2.0_TKG, 1.8_TKG, 1.6_TKG, 1.4_TKG, 1.2_TKG, 1.0_TKG, 1.2_TKG, 1.4_TKG &
                            , 2.4_TKG, +2.2_TKG, 2.0_TKG, 1.8_TKG, 1.6_TKG, 1.4_TKG, 1.2_TKG, 1.0_TKG, 1.2_TKG &
                            , 2.6_TKG, +2.4_TKG, 2.2_TKG, 2.0_TKG, 1.8_TKG, 1.6_TKG, 1.4_TKG, 1.2_TKG, 1.0_TKG &
                            ], shape = [9, 9], order = [2, 1])
        lup_ref = reshape(  [ 2.6_TKG,  2.4_TKG, 2.2_TKG, 2.0_TKG, 1.8_TKG, 1.6_TKG, 1.4_TKG, 1.2_TKG, 1.0_TKG &
                            , 0.4_TKG,  0.3_TKG, 0.6_TKG, 0.8_TKG, 1.1_TKG, 1.4_TKG, 1.7_TKG, 1.9_TKG, 2.2_TKG &
                            , 0.5_TKG, -0.4_TKG, 0.4_TKG, 0.8_TKG, 1.2_TKG, 1.6_TKG, 2.0_TKG, 2.4_TKG, 2.8_TKG &
                            , 0.5_TKG, -0.3_TKG, 0.0_TKG, 0.4_TKG, 0.8_TKG, 1.2_TKG, 1.6_TKG, 2.0_TKG, 2.4_TKG &
                            , 0.6_TKG, -0.3_TKG, 0.0_TKG, 0.0_TKG, 0.4_TKG, 0.8_TKG, 1.2_TKG, 1.6_TKG, 2.0_TKG &
                            , 0.7_TKG, -0.2_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.4_TKG, 0.8_TKG, 1.2_TKG, 1.6_TKG &
                            , 0.8_TKG, -0.2_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.4_TKG, 0.8_TKG, 1.2_TKG &
                            , 0.8_TKG, -0.1_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.4_TKG, 0.8_TKG &
                            , 0.9_TKG, -0.1_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.4_TKG &
                            ], shape = [9, 9], order = [2, 1])
        rperm_ref = [9, 9, 9, 9, 9, 9, 9, 9, 9]
        call disp%skip
        call disp%show("mat_lup")
        call disp%show( mat_lup )
        call disp%show("call setResized(rperm, size(mat_lup, 1, IK))")
                        call setResized(rperm, size(mat_lup, 1, IK))
        call disp%show("call setMatLUP(mat_lup, rperm, info); if (info /= 0) error stop")
                        call setMatLUP(mat_lup, rperm, info); if (info /= 0) error stop
        call disp%show("mat_lup")
        call disp%show( mat_lup )
        call disp%show("call setResized(rperm, size(mat_lup, 1, IK)) ! reconstruct the original matrix.")
                        call setResized(rperm, size(mat_lup, 1, IK))
        call disp%show("lup_ref ! reference matrix rounded to 1 significant digit.")
        call disp%show( lup_ref )
        call disp%show("lup_ref - mat_lup, format = SK_'(*(f0.1,:,"", ""))'")
        call disp%show( lup_ref - mat_lup, format = SK_"(*(f0.1,:,"", ""))")
        call disp%show("rperm_ref")
        call disp%show( rperm_ref )
        call disp%show("rperm")
        call disp%show( rperm )
        call disp%skip
    end block

    block
        use pm_kind, only: TKG => RKH
        real(TKG), allocatable :: mat_lup(:,:), lup_ref(:,:)
        mat_lup = reshape(  [ 1.0_TKG, +1.2_TKG, 1.4_TKG, 1.6_TKG, 1.8_TKG, 2.0_TKG, 2.2_TKG, 2.4_TKG, 2.6_TKG &
                            , 1.2_TKG, +1.0_TKG, 1.2_TKG, 1.4_TKG, 1.6_TKG, 1.8_TKG, 2.0_TKG, 2.2_TKG, 2.4_TKG &
                            , 1.4_TKG, +1.2_TKG, 1.0_TKG, 1.2_TKG, 1.4_TKG, 1.6_TKG, 1.8_TKG, 2.0_TKG, 2.2_TKG &
                            , 1.6_TKG, +1.4_TKG, 1.2_TKG, 1.0_TKG, 1.2_TKG, 1.4_TKG, 1.6_TKG, 1.8_TKG, 2.0_TKG &
                            , 1.8_TKG, +1.6_TKG, 1.4_TKG, 1.2_TKG, 1.0_TKG, 1.2_TKG, 1.4_TKG, 1.6_TKG, 1.8_TKG &
                            , 2.0_TKG, +1.8_TKG, 1.6_TKG, 1.4_TKG, 1.2_TKG, 1.0_TKG, 1.2_TKG, 1.4_TKG, 1.6_TKG &
                            , 2.2_TKG, +2.0_TKG, 1.8_TKG, 1.6_TKG, 1.4_TKG, 1.2_TKG, 1.0_TKG, 1.2_TKG, 1.4_TKG &
                            , 2.4_TKG, +2.2_TKG, 2.0_TKG, 1.8_TKG, 1.6_TKG, 1.4_TKG, 1.2_TKG, 1.0_TKG, 1.2_TKG &
                            , 2.6_TKG, +2.4_TKG, 2.2_TKG, 2.0_TKG, 1.8_TKG, 1.6_TKG, 1.4_TKG, 1.2_TKG, 1.0_TKG &
                            ], shape = [9, 9], order = [2, 1])
        lup_ref = reshape(  [ 2.6_TKG,  2.4_TKG, 2.2_TKG, 2.0_TKG, 1.8_TKG, 1.6_TKG, 1.4_TKG, 1.2_TKG, 1.0_TKG &
                            , 0.4_TKG,  0.3_TKG, 0.6_TKG, 0.8_TKG, 1.1_TKG, 1.4_TKG, 1.7_TKG, 1.9_TKG, 2.2_TKG &
                            , 0.5_TKG, -0.4_TKG, 0.4_TKG, 0.8_TKG, 1.2_TKG, 1.6_TKG, 2.0_TKG, 2.4_TKG, 2.8_TKG &
                            , 0.5_TKG, -0.3_TKG, 0.0_TKG, 0.4_TKG, 0.8_TKG, 1.2_TKG, 1.6_TKG, 2.0_TKG, 2.4_TKG &
                            , 0.6_TKG, -0.3_TKG, 0.0_TKG, 0.0_TKG, 0.4_TKG, 0.8_TKG, 1.2_TKG, 1.6_TKG, 2.0_TKG &
                            , 0.7_TKG, -0.2_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.4_TKG, 0.8_TKG, 1.2_TKG, 1.6_TKG &
                            , 0.8_TKG, -0.2_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.4_TKG, 0.8_TKG, 1.2_TKG &
                            , 0.8_TKG, -0.1_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.4_TKG, 0.8_TKG &
                            , 0.9_TKG, -0.1_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.0_TKG, 0.4_TKG &
                            ], shape = [9, 9], order = [2, 1])
        rperm_ref = [9, 9, 9, 9, 9, 9, 9, 9, 9]
        call disp%skip
        call disp%show("mat_lup")
        call disp%show( mat_lup )
        call disp%show("call setResized(rperm, size(mat_lup, 1, IK))")
                        call setResized(rperm, size(mat_lup, 1, IK))
        call disp%show("call setMatLUP(mat_lup, rperm, info); if (info /= 0) error stop")
                        call setMatLUP(mat_lup, rperm, info); if (info /= 0) error stop
        call disp%show("mat_lup")
        call disp%show( mat_lup )
        call disp%show("call setResized(rperm, size(mat_lup, 1, IK)) ! reconstruct the original matrix.")
                        call setResized(rperm, size(mat_lup, 1, IK))
        call disp%show("lup_ref ! reference matrix rounded to 1 significant digit.")
        call disp%show( lup_ref )
        call disp%show("lup_ref - mat_lup, format = SK_'(*(f0.1,:,"", ""))'")
        call disp%show( lup_ref - mat_lup, format = SK_"(*(f0.1,:,"", ""))")
        call disp%show("rperm_ref")
        call disp%show( rperm_ref )
        call disp%show("rperm")
        call disp%show( rperm )
        call disp%skip
    end block

#if 1 || LAPACK_ENABLED
    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the LU-Pivoted decomposition of a complex square matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    block
        use pm_kind, only: TKG => CKS
        complex(TKG), allocatable :: mat_lup(:,:), lup_ref(:,:)
         mat_lup = reshape( [ (2.0, 1.0), (2.4,-1.0), (2.8,-1.0), (3.2,-1.0), (3.6,-1.0), (4.0,-1.0), (4.4,-1.0), (4.8,-1.0), (5.2,-1.0) &
                            , (2.4, 1.0), (2.0, 1.0), (2.4,-1.0), (2.8,-1.0), (3.2,-1.0), (3.6,-1.0), (4.0,-1.0), (4.4,-1.0), (4.8,-1.0) &
                            , (2.8, 1.0), (2.4, 1.0), (2.0, 1.0), (2.4,-1.0), (2.8,-1.0), (3.2,-1.0), (3.6,-1.0), (4.0,-1.0), (4.4,-1.0) &
                            , (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0), (2.4,-1.0), (2.8,-1.0), (3.2,-1.0), (3.6,-1.0), (4.0,-1.0) &
                            , (3.6, 1.0), (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0), (2.4,-1.0), (2.8,-1.0), (3.2,-1.0), (3.6,-1.0) &
                            , (4.0, 1.0), (3.6, 1.0), (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0), (2.4,-1.0), (2.8,-1.0), (3.2,-1.0) &
                            , (4.4, 1.0), (4.0, 1.0), (3.6, 1.0), (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0), (2.4,-1.0), (2.8,-1.0) &
                            , (4.8, 1.0), (4.4, 1.0), (4.0, 1.0), (3.6, 1.0), (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0), (2.4,-1.0) &
                            , (5.2, 1.0), (4.8, 1.0), (4.4, 1.0), (4.0, 1.0), (3.6, 1.0), (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0) &
                            ], shape = [9, 9], order = [2, 1])
        lup_ref = reshape(  [ (5.2, 1.0), (4.8, 1.0), (4.4, 1.0), (4.0, 1.0), (+3.6, 1.0), (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0) &
                            , (0.4, 0.1), (0.6,-2.0), (1.1,-1.9), (1.7,-1.9), (+2.3,-1.8), (2.8,-1.8), (3.4,-1.7), (3.9,-1.7), (4.5,-1.6) &
                            , (0.5, 0.1), (0.0,-0.1), (0.6,-1.9), (1.2,-1.8), (+1.8,-1.7), (2.5,-1.6), (3.1,-1.5), (3.7,-1.4), (4.3,-1.3) &
                            , (0.6, 0.1), (0.0,-0.1),(-0.1,-0.1), (0.7,-1.9), (+1.3,-1.7), (2.0,-1.6), (2.7,-1.5), (3.4,-1.4), (4.0,-1.2) &
                            , (0.6, 0.1), (0.0,-0.1),(-0.1,-0.1),(-0.1, 0.0), (+0.7,-1.9), (1.5,-1.7), (2.2,-1.6), (2.9,-1.5), (3.7,-1.3) &
                            , (0.7, 0.1), (0.0,-0.1), (0.0, 0.0),(-0.1, 0.0), (-0.1, 0.0), (0.8,-1.9), (1.6,-1.8), (2.4,-1.6), (3.2,-1.5) &
                            , (0.8, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (+0.0, 0.0), (0.0, 0.0), (0.8,-1.9), (1.7,-1.8), (2.5,-1.8) &
                            , (0.9, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (+0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.8,-2.0), (1.7,-1.9) &
                            , (0.9, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (+0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.8,-2.0) &
                            ], shape = [9, 9], order = [2, 1])
        rperm_ref = [9, 9, 9, 9, 9, 9, 9, 9, 9]
        call disp%skip
        call disp%show("mat_lup")
        call disp%show( mat_lup )
        call disp%show("call setResized(rperm, size(mat_lup, 1, IK))")
                        call setResized(rperm, size(mat_lup, 1, IK))
        call disp%show("call setMatLUP(mat_lup, rperm, info); if (info /= 0) error stop")
                        call setMatLUP(mat_lup, rperm, info); if (info /= 0) error stop
        call disp%show("mat_lup")
        call disp%show( mat_lup )
        call disp%show("call setResized(rperm, size(mat_lup, 1, IK)) ! reconstruct the original matrix.")
                        call setResized(rperm, size(mat_lup, 1, IK))
        call disp%show("lup_ref ! reference matrix rounded to 1 significant digit.")
        call disp%show( lup_ref )
        call disp%show("lup_ref - mat_lup, format = SK_'(*(""("",f0.1,"","",f0.1,"")"",:,"", ""))'")
        call disp%show( lup_ref - mat_lup, format = SK_"(*(""("",f0.1,"","",f0.1,"")"",:,"", ""))")
        call disp%show("rperm_ref")
        call disp%show( rperm_ref )
        call disp%show("rperm")
        call disp%show( rperm )
        call disp%skip
    end block

    block
        use pm_kind, only: TKG => CKD
        complex(TKG), allocatable :: mat_lup(:,:), lup_ref(:,:)
         mat_lup = reshape( [ (2.0, 1.0), (2.4,-1.0), (2.8,-1.0), (3.2,-1.0), (3.6,-1.0), (4.0,-1.0), (4.4,-1.0), (4.8,-1.0), (5.2,-1.0) &
                            , (2.4, 1.0), (2.0, 1.0), (2.4,-1.0), (2.8,-1.0), (3.2,-1.0), (3.6,-1.0), (4.0,-1.0), (4.4,-1.0), (4.8,-1.0) &
                            , (2.8, 1.0), (2.4, 1.0), (2.0, 1.0), (2.4,-1.0), (2.8,-1.0), (3.2,-1.0), (3.6,-1.0), (4.0,-1.0), (4.4,-1.0) &
                            , (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0), (2.4,-1.0), (2.8,-1.0), (3.2,-1.0), (3.6,-1.0), (4.0,-1.0) &
                            , (3.6, 1.0), (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0), (2.4,-1.0), (2.8,-1.0), (3.2,-1.0), (3.6,-1.0) &
                            , (4.0, 1.0), (3.6, 1.0), (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0), (2.4,-1.0), (2.8,-1.0), (3.2,-1.0) &
                            , (4.4, 1.0), (4.0, 1.0), (3.6, 1.0), (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0), (2.4,-1.0), (2.8,-1.0) &
                            , (4.8, 1.0), (4.4, 1.0), (4.0, 1.0), (3.6, 1.0), (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0), (2.4,-1.0) &
                            , (5.2, 1.0), (4.8, 1.0), (4.4, 1.0), (4.0, 1.0), (3.6, 1.0), (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0) &
                            ], shape = [9, 9], order = [2, 1])
        lup_ref = reshape(  [ (5.2, 1.0), (4.8, 1.0), (4.4, 1.0), (4.0, 1.0), (+3.6, 1.0), (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0) &
                            , (0.4, 0.1), (0.6,-2.0), (1.1,-1.9), (1.7,-1.9), (+2.3,-1.8), (2.8,-1.8), (3.4,-1.7), (3.9,-1.7), (4.5,-1.6) &
                            , (0.5, 0.1), (0.0,-0.1), (0.6,-1.9), (1.2,-1.8), (+1.8,-1.7), (2.5,-1.6), (3.1,-1.5), (3.7,-1.4), (4.3,-1.3) &
                            , (0.6, 0.1), (0.0,-0.1),(-0.1,-0.1), (0.7,-1.9), (+1.3,-1.7), (2.0,-1.6), (2.7,-1.5), (3.4,-1.4), (4.0,-1.2) &
                            , (0.6, 0.1), (0.0,-0.1),(-0.1,-0.1),(-0.1, 0.0), (+0.7,-1.9), (1.5,-1.7), (2.2,-1.6), (2.9,-1.5), (3.7,-1.3) &
                            , (0.7, 0.1), (0.0,-0.1), (0.0, 0.0),(-0.1, 0.0), (-0.1, 0.0), (0.8,-1.9), (1.6,-1.8), (2.4,-1.6), (3.2,-1.5) &
                            , (0.8, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (+0.0, 0.0), (0.0, 0.0), (0.8,-1.9), (1.7,-1.8), (2.5,-1.8) &
                            , (0.9, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (+0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.8,-2.0), (1.7,-1.9) &
                            , (0.9, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (+0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.8,-2.0) &
                            ], shape = [9, 9], order = [2, 1])
        rperm_ref = [9, 9, 9, 9, 9, 9, 9, 9, 9]
        call disp%skip
        call disp%show("mat_lup")
        call disp%show( mat_lup )
        call disp%show("call setResized(rperm, size(mat_lup, 1, IK))")
                        call setResized(rperm, size(mat_lup, 1, IK))
        call disp%show("call setMatLUP(mat_lup, rperm, info); if (info /= 0) error stop")
                        call setMatLUP(mat_lup, rperm, info); if (info /= 0) error stop
        call disp%show("mat_lup")
        call disp%show( mat_lup )
        call disp%show("call setResized(rperm, size(mat_lup, 1, IK)) ! reconstruct the original matrix.")
                        call setResized(rperm, size(mat_lup, 1, IK))
        call disp%show("lup_ref ! reference matrix rounded to 1 significant digit.")
        call disp%show( lup_ref )
        call disp%show("lup_ref - mat_lup, format = SK_'(*(""("",f0.1,"","",f0.1,"")"",:,"", ""))'")
        call disp%show( lup_ref - mat_lup, format = SK_"(*(""("",f0.1,"","",f0.1,"")"",:,"", ""))")
        call disp%show("rperm_ref")
        call disp%show( rperm_ref )
        call disp%show("rperm")
        call disp%show( rperm )
        call disp%skip
    end block

    block
        use pm_kind, only: TKG => CKH
        complex(TKG), allocatable :: mat_lup(:,:), lup_ref(:,:)
         mat_lup = reshape( [ (2.0, 1.0), (2.4,-1.0), (2.8,-1.0), (3.2,-1.0), (3.6,-1.0), (4.0,-1.0), (4.4,-1.0), (4.8,-1.0), (5.2,-1.0) &
                            , (2.4, 1.0), (2.0, 1.0), (2.4,-1.0), (2.8,-1.0), (3.2,-1.0), (3.6,-1.0), (4.0,-1.0), (4.4,-1.0), (4.8,-1.0) &
                            , (2.8, 1.0), (2.4, 1.0), (2.0, 1.0), (2.4,-1.0), (2.8,-1.0), (3.2,-1.0), (3.6,-1.0), (4.0,-1.0), (4.4,-1.0) &
                            , (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0), (2.4,-1.0), (2.8,-1.0), (3.2,-1.0), (3.6,-1.0), (4.0,-1.0) &
                            , (3.6, 1.0), (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0), (2.4,-1.0), (2.8,-1.0), (3.2,-1.0), (3.6,-1.0) &
                            , (4.0, 1.0), (3.6, 1.0), (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0), (2.4,-1.0), (2.8,-1.0), (3.2,-1.0) &
                            , (4.4, 1.0), (4.0, 1.0), (3.6, 1.0), (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0), (2.4,-1.0), (2.8,-1.0) &
                            , (4.8, 1.0), (4.4, 1.0), (4.0, 1.0), (3.6, 1.0), (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0), (2.4,-1.0) &
                            , (5.2, 1.0), (4.8, 1.0), (4.4, 1.0), (4.0, 1.0), (3.6, 1.0), (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0) &
                            ], shape = [9, 9], order = [2, 1])
        lup_ref = reshape(  [ (5.2, 1.0), (4.8, 1.0), (4.4, 1.0), (4.0, 1.0), (+3.6, 1.0), (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0) &
                            , (0.4, 0.1), (0.6,-2.0), (1.1,-1.9), (1.7,-1.9), (+2.3,-1.8), (2.8,-1.8), (3.4,-1.7), (3.9,-1.7), (4.5,-1.6) &
                            , (0.5, 0.1), (0.0,-0.1), (0.6,-1.9), (1.2,-1.8), (+1.8,-1.7), (2.5,-1.6), (3.1,-1.5), (3.7,-1.4), (4.3,-1.3) &
                            , (0.6, 0.1), (0.0,-0.1),(-0.1,-0.1), (0.7,-1.9), (+1.3,-1.7), (2.0,-1.6), (2.7,-1.5), (3.4,-1.4), (4.0,-1.2) &
                            , (0.6, 0.1), (0.0,-0.1),(-0.1,-0.1),(-0.1, 0.0), (+0.7,-1.9), (1.5,-1.7), (2.2,-1.6), (2.9,-1.5), (3.7,-1.3) &
                            , (0.7, 0.1), (0.0,-0.1), (0.0, 0.0),(-0.1, 0.0), (-0.1, 0.0), (0.8,-1.9), (1.6,-1.8), (2.4,-1.6), (3.2,-1.5) &
                            , (0.8, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (+0.0, 0.0), (0.0, 0.0), (0.8,-1.9), (1.7,-1.8), (2.5,-1.8) &
                            , (0.9, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (+0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.8,-2.0), (1.7,-1.9) &
                            , (0.9, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (+0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.8,-2.0) &
                            ], shape = [9, 9], order = [2, 1])
        rperm_ref = [9, 9, 9, 9, 9, 9, 9, 9, 9]
        call disp%skip
        call disp%show("mat_lup")
        call disp%show( mat_lup )
        call disp%show("call setResized(rperm, size(mat_lup, 1, IK))")
                        call setResized(rperm, size(mat_lup, 1, IK))
        call disp%show("call setMatLUP(mat_lup, rperm, info); if (info /= 0) error stop")
                        call setMatLUP(mat_lup, rperm, info); if (info /= 0) error stop
        call disp%show("mat_lup")
        call disp%show( mat_lup )
        call disp%show("call setResized(rperm, size(mat_lup, 1, IK)) ! reconstruct the original matrix.")
                        call setResized(rperm, size(mat_lup, 1, IK))
        call disp%show("lup_ref ! reference matrix rounded to 1 significant digit.")
        call disp%show( lup_ref )
        call disp%show("lup_ref - mat_lup, format = SK_'(*(""("",f0.1,"","",f0.1,"")"",:,"", ""))'")
        call disp%show( lup_ref - mat_lup, format = SK_"(*(""("",f0.1,"","",f0.1,"")"",:,"", ""))")
        call disp%show("rperm_ref")
        call disp%show( rperm_ref )
        call disp%show("rperm")
        call disp%show( rperm )
        call disp%skip
    end block
#endif

end program example