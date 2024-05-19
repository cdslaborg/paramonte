program example

    use pm_kind, only: SK
    use pm_kind, only: IK, RK ! all real kinds are supported.
    use pm_io, only: display_type
    use pm_distNorm, only: getZigNorm
    use pm_arrayRebind, only: setRebound

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: RKG => RKH
        real(RKG), allocatable :: zig(:,:)
        integer(IK) :: nlayer
        real(RKG) :: abserr
        call disp%skip()
        call disp%show("nlayer = 256")
                        nlayer = 256
        call disp%show("call setRebound(zig, [1_IK, 0_IK], [2_IK, nlayer])")
                        call setRebound(zig, [1_IK, 0_IK], [2_IK, nlayer])
        call disp%show("zig(:,:) = getZigNorm(nlayer, abserr, abstol = epsilon(0._RKG)) ! rectangle rightmost corners and the corresponding function values.")
                        zig(:,:) = getZigNorm(nlayer, abserr, abstol = epsilon(0._RKG))
        call disp%show("[nlayer, shape(zig, IK)]")
        call disp%show( [nlayer, shape(zig, IK)] )
        call disp%show("abserr")
        call disp%show( abserr )
        call disp%show("zig(1, 1) ! The upper rightmost corner of the lowest rectangle (which has a tail) in the ziggurat set.")
        call disp%show( zig(1, 1) )
        call disp%show("transpose(zig)")
        call disp%show( transpose(zig) )!, format = "(sp,*(ES0.37E1,:,', '))")
        call disp%show("reshape(zig(1, 1 : nlayer - 1) * (zig(2, 2 : nlayer) - zig(2, 1 : nlayer - 1)), [nlayer - 1, 1]) ! area of each partition.")
        call disp%show( reshape(zig(1, 1 : nlayer - 1) * (zig(2, 2 : nlayer) - zig(2, 1 : nlayer - 1)), [nlayer - 1, 1]) )
        call disp%skip()
    end block

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !Output an example array for illustration.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_io, only: getErrTableWrite, trans
        use pm_kind, only: RKG => RKH
        real(RKG) :: abserr
        if (0 /= getErrTableWrite(SK_"getZigNorm.RK.txt", getZigNorm(64_IK, abserr), trans)) error stop "Failed to write the table to the file."
    end block

end program example