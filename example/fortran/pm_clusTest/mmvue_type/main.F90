program example

    use pm_kind, only: SK, IK, LK
    use pm_clusTest, only: mmvue_type, rngf_type
    use pm_io, only: getErrTableWrite, trans
    use pm_io, only: display_type
    use pm_val2str, only: getStr

    implicit none

    integer(IK) :: itry
    type(mmvue_type) :: mmvue
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    do itry = 0, 9
        call disp%skip()
        call disp%show("mmvue = mmvue_type(rng = rngf_type())")
                        mmvue = mmvue_type(rng = rngf_type())
        call disp%show("[mmvue%ndim, mmvue%nell, mmvue%nsam, mmvue%nsim]")
        call disp%show( [mmvue%ndim, mmvue%nell, mmvue%nsam, mmvue%nsim] )
        call disp%show("if (0 /= getErrTableWrite(SK_'mmvue_type.'//getStr(itry)//SK_'.txt', reshape([transpose(mmvue%sample), real(mmvue%membership, kind(mmvue%sample))], [mmvue%nsam, mmvue%ndim + 1_IK]))) error stop 'Failed to write the random vectors file.'")
                        if (0 /= getErrTableWrite(SK_'mmvue_type.'//getStr(itry)//SK_'.txt', reshape([transpose(mmvue%sample), real(mmvue%membership, kind(mmvue%sample))], [mmvue%nsam, mmvue%ndim + 1_IK]))) error stop 'Failed to write the random vectors file.'
        call disp%skip()
    end do

end program example