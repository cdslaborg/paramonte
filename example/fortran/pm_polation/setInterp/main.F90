program example

    use pm_kind, only: SK, IK, LK
    use pm_polation, only: setInterp, neimean, neinear, neinext, neiprev, piwilin, monopol
    use pm_arrayResize, only: setResized
    use pm_arraySpace, only: getLinSpace
    use pm_io, only: getErrTableWrite
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    integer(IK) :: itry, ntry = 5
    integer(IK) :: dim, isam, ndim, nsam
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%show("!1D interpolation.")
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: RKC => RKS ! all processor kinds are supported.
        use pm_mathConst, only: PI_RKH => PI
        real(RKC), parameter :: PI = PI_RKH
        real(RKC), allocatable :: crdx(:), func(:), queryx(:), interp(:), relerr(:)
        integer(IK) :: ncrd, nquery
        call disp%skip()
        call disp%show("ncrd = 9; nquery = 33")
                        ncrd = 9; nquery = 33
        call disp%show("crdx = getLinSpace(0._RKC, 2 * PI, ncrd)")
                        crdx = getLinSpace(0._RKC, 2 * PI, ncrd)
        call disp%show("crdx")
        call disp%show( crdx )
        call disp%show("func = sin(crdx)")
                        func = sin(crdx)
        call disp%show("func")
        call disp%show( func )
        call disp%show("queryx = getLinSpace(0._RKC, 2 * PI, nquery)")
                        queryx = getLinSpace(0._RKC, 2 * PI, nquery)
        call disp%show("queryx")
        call disp%show( queryx )
        call disp%show("call setResized(interp, size(queryx, 1, IK))")
                        call setResized(interp, size(queryx, 1, IK))
        call disp%show("call setResized(relerr, size(queryx, 1, IK))")
                        call setResized(relerr, size(queryx, 1, IK))
        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!1D mean neighbors interpolation.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()
        call disp%show("if (0 /= getErrTableWrite(SK_'setInterp.neimean.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'interpolation node outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setInterp.neimean.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'interpolation node outputting failed.'
        call disp%show("call setInterp(neimean, crdx, func, queryx, interp)")
                        call setInterp(neimean, crdx, func, queryx, interp)
        call disp%show("interp")
        call disp%show( interp )
        call disp%show("if (0 /= getErrTableWrite(SK_'setInterp.neimean.interp.txt', reshape([queryx, interp], [nquery, 2_IK]), header = SK_'queryx,interp')) error stop 'interpolation outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setInterp.neimean.interp.txt', reshape([queryx, interp], [nquery, 2_IK]), header = SK_'queryx,interp')) error stop 'interpolation outputting failed.'
        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!1D nearest neighbor interpolation.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()
        call disp%show("if (0 /= getErrTableWrite(SK_'setInterp.neinear.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'interpolation node outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setInterp.neinear.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'interpolation node outputting failed.'
        call disp%show("call setInterp(neinear, crdx, func, queryx, interp)")
                        call setInterp(neinear, crdx, func, queryx, interp)
        call disp%show("interp")
        call disp%show( interp )
        call disp%show("if (0 /= getErrTableWrite(SK_'setInterp.neinear.interp.txt', reshape([queryx, interp], [nquery, 2_IK]), header = SK_'queryx,interp')) error stop 'interpolation outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setInterp.neinear.interp.txt', reshape([queryx, interp], [nquery, 2_IK]), header = SK_'queryx,interp')) error stop 'interpolation outputting failed.'
        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!1D next neighbor interpolation.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()
        call disp%show("if (0 /= getErrTableWrite(SK_'setInterp.neinext.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'interpolation node outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setInterp.neinext.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'interpolation node outputting failed.'
        call disp%show("call setInterp(neinext, crdx, func, queryx, interp)")
                        call setInterp(neinext, crdx, func, queryx, interp)
        call disp%show("interp")
        call disp%show( interp )
        call disp%show("if (0 /= getErrTableWrite(SK_'setInterp.neinext.interp.txt', reshape([queryx, interp], [nquery, 2_IK]), header = SK_'queryx,interp')) error stop 'interpolation outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setInterp.neinext.interp.txt', reshape([queryx, interp], [nquery, 2_IK]), header = SK_'queryx,interp')) error stop 'interpolation outputting failed.'
        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!1D previous neighbor interpolation.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()
        call disp%show("if (0 /= getErrTableWrite(SK_'setInterp.neiprev.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'interpolation node outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setInterp.neiprev.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'interpolation node outputting failed.'
        call disp%show("call setInterp(neiprev, crdx, func, queryx, interp)")
                        call setInterp(neiprev, crdx, func, queryx, interp)
        call disp%show("interp")
        call disp%show( interp )
        call disp%show("if (0 /= getErrTableWrite(SK_'setInterp.neiprev.interp.txt', reshape([queryx, interp], [nquery, 2_IK]), header = SK_'queryx,interp')) error stop 'interpolation outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setInterp.neiprev.interp.txt', reshape([queryx, interp], [nquery, 2_IK]), header = SK_'queryx,interp')) error stop 'interpolation outputting failed.'
        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!1D piecewise linear interpolation.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()
        call disp%show("if (0 /= getErrTableWrite(SK_'setInterp.piwilin.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'interpolation node outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setInterp.piwilin.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'interpolation node outputting failed.'
        call disp%show("call setInterp(piwilin, crdx, func, queryx, interp)")
                        call setInterp(piwilin, crdx, func, queryx, interp)
        call disp%show("interp")
        call disp%show( interp )
        call disp%show("if (0 /= getErrTableWrite(SK_'setInterp.piwilin.interp.txt', reshape([queryx, interp], [nquery, 2_IK]), header = SK_'queryx,interp')) error stop 'interpolation outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setInterp.piwilin.interp.txt', reshape([queryx, interp], [nquery, 2_IK]), header = SK_'queryx,interp')) error stop 'interpolation outputting failed.'
        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!1D single polynomial interpolation.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()
        call disp%show("if (0 /= getErrTableWrite(SK_'setInterp.monopol.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'interpolation node outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setInterp.monopol.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'interpolation node outputting failed.'
        call disp%show("call setInterp(monopol, crdx, func, queryx, interp, relerr)")
                        call setInterp(monopol, crdx, func, queryx, interp, relerr)
        call disp%show("interp")
        call disp%show( interp )
        call disp%show("relerr")
        call disp%show( relerr )
        call disp%show("if (0 /= getErrTableWrite(SK_'setInterp.monopol.interp.txt', reshape([queryx, interp], [nquery, 2_IK]), header = SK_'queryx,interp')) error stop 'interpolation outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setInterp.monopol.interp.txt', reshape([queryx, interp], [nquery, 2_IK]), header = SK_'queryx,interp')) error stop 'interpolation outputting failed.'
        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!1D single polynomial interpolation showing Runge effect.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()
        block
            use pm_distNorm, only: getNormLogPDF
            call disp%show("ncrd = 9; nquery = 50 * ncrd")
                            ncrd = 9; nquery = 50 * ncrd
            call disp%show("crdx = getLinSpace(-5., +5., ncrd)")
                            crdx = getLinSpace(-5., +5., ncrd)
            call disp%show("crdx")
            call disp%show( crdx )
            call disp%show("func = exp(getNormLogPDF(crdx))")
                            func = exp(getNormLogPDF(crdx))
            call disp%show("func")
            call disp%show( func )
            call disp%show("if (0 /= getErrTableWrite(SK_'setInterp.rungeEffect.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'interpolation node outputting failed.'")
                            if (0 /= getErrTableWrite(SK_'setInterp.rungeEffect.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'interpolation node outputting failed.'
            call disp%show("queryx = getLinSpace(minval(crdx, 1), maxval(crdx, 1), nquery)")
                            queryx = getLinSpace(minval(crdx, 1), maxval(crdx, 1), nquery)
            call disp%show("queryx")
            call disp%show( queryx )
            call disp%show("call setResized(interp, size(queryx, 1, IK))")
                            call setResized(interp, size(queryx, 1, IK))
            call disp%show("call setResized(relerr, size(queryx, 1, IK))")
                            call setResized(relerr, size(queryx, 1, IK))
            call disp%show("call setInterp(monopol, crdx, func, queryx, interp, relerr)")
                            call setInterp(monopol, crdx, func, queryx, interp, relerr)
            call disp%show("interp")
            call disp%show( interp )
            call disp%show("relerr")
            call disp%show( relerr )
            call disp%show("if (0 /= getErrTableWrite(SK_'setInterp.rungeEffect.interp.txt', reshape([queryx, interp], [nquery, 2_IK]), header = SK_'queryx,interp')) error stop 'interpolation outputting failed.'")
                            if (0 /= getErrTableWrite(SK_'setInterp.rungeEffect.interp.txt', reshape([queryx, interp], [nquery, 2_IK]), header = SK_'queryx,interp')) error stop 'interpolation outputting failed.'
            call disp%skip()
        end block
    end block

end program example