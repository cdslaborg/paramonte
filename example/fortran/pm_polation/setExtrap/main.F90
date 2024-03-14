program example

    use pm_kind, only: SK, IK, LK
    use pm_polation, only: setExtrap, neimean, neinear, neinext, neiprev, piwilin, monopol
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
    call disp%show("!1D extrapolation.")
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: RKC => RKS ! all processor kinds are supported.
        use pm_mathConst, only: PI_RKH => PI
        real(RKC), parameter :: PI = PI_RKH
        real(RKC), allocatable :: crdx(:), func(:), queryx(:), extrap(:), relerr(:)
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
        call disp%show("queryx = getLinSpace(-PI/2, 5 * PI/2, nquery)")
                        queryx = getLinSpace(-PI/2, 5 * PI/2, nquery)
        call disp%show("queryx")
        call disp%show( queryx )
        call disp%show("call setResized(extrap, size(queryx, 1, IK))")
                        call setResized(extrap, size(queryx, 1, IK))
        call disp%show("call setResized(relerr, size(queryx, 1, IK))")
                        call setResized(relerr, size(queryx, 1, IK))
        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!1D mean neighbors extrapolation.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()
        call disp%show("if (0 /= getErrTableWrite(SK_'setExtrap.neimean.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'extrapolation node outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setExtrap.neimean.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'extrapolation node outputting failed.'
        call disp%show("call setExtrap(neimean, crdx, func, queryx, extrap)")
                        call setExtrap(neimean, crdx, func, queryx, extrap)
        call disp%show("extrap")
        call disp%show( extrap )
        call disp%show("if (0 /= getErrTableWrite(SK_'setExtrap.neimean.extrap.txt', reshape([queryx, extrap], [nquery, 2_IK]), header = SK_'queryx,extrap')) error stop 'extrapolation outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setExtrap.neimean.extrap.txt', reshape([queryx, extrap], [nquery, 2_IK]), header = SK_'queryx,extrap')) error stop 'extrapolation outputting failed.'
        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!1D nearest neighbor extrapolation.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()
        call disp%show("if (0 /= getErrTableWrite(SK_'setExtrap.neinear.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'extrapolation node outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setExtrap.neinear.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'extrapolation node outputting failed.'
        call disp%show("call setExtrap(neinear, crdx, func, queryx, extrap)")
                        call setExtrap(neinear, crdx, func, queryx, extrap)
        call disp%show("extrap")
        call disp%show( extrap )
        call disp%show("if (0 /= getErrTableWrite(SK_'setExtrap.neinear.extrap.txt', reshape([queryx, extrap], [nquery, 2_IK]), header = SK_'queryx,extrap')) error stop 'extrapolation outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setExtrap.neinear.extrap.txt', reshape([queryx, extrap], [nquery, 2_IK]), header = SK_'queryx,extrap')) error stop 'extrapolation outputting failed.'
        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!1D next neighbor extrapolation.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()
        call disp%show("if (0 /= getErrTableWrite(SK_'setExtrap.neinext.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'extrapolation node outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setExtrap.neinext.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'extrapolation node outputting failed.'
        call disp%show("call setExtrap(neinext, crdx, func, queryx, extrap)")
                        call setExtrap(neinext, crdx, func, queryx, extrap)
        call disp%show("extrap")
        call disp%show( extrap )
        call disp%show("if (0 /= getErrTableWrite(SK_'setExtrap.neinext.extrap.txt', reshape([queryx, extrap], [nquery, 2_IK]), header = SK_'queryx,extrap')) error stop 'extrapolation outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setExtrap.neinext.extrap.txt', reshape([queryx, extrap], [nquery, 2_IK]), header = SK_'queryx,extrap')) error stop 'extrapolation outputting failed.'
        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!1D previous neighbor extrapolation.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()
        call disp%show("if (0 /= getErrTableWrite(SK_'setExtrap.neiprev.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'extrapolation node outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setExtrap.neiprev.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'extrapolation node outputting failed.'
        call disp%show("call setExtrap(neiprev, crdx, func, queryx, extrap)")
                        call setExtrap(neiprev, crdx, func, queryx, extrap)
        call disp%show("extrap")
        call disp%show( extrap )
        call disp%show("if (0 /= getErrTableWrite(SK_'setExtrap.neiprev.extrap.txt', reshape([queryx, extrap], [nquery, 2_IK]), header = SK_'queryx,extrap')) error stop 'extrapolation outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setExtrap.neiprev.extrap.txt', reshape([queryx, extrap], [nquery, 2_IK]), header = SK_'queryx,extrap')) error stop 'extrapolation outputting failed.'
        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!1D piecewise linear extrapolation.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()
        call disp%show("if (0 /= getErrTableWrite(SK_'setExtrap.piwilin.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'extrapolation node outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setExtrap.piwilin.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'extrapolation node outputting failed.'
        call disp%show("call setExtrap(piwilin, crdx, func, queryx, extrap)")
                        call setExtrap(piwilin, crdx, func, queryx, extrap)
        call disp%show("extrap")
        call disp%show( extrap )
        call disp%show("if (0 /= getErrTableWrite(SK_'setExtrap.piwilin.extrap.txt', reshape([queryx, extrap], [nquery, 2_IK]), header = SK_'queryx,extrap')) error stop 'extrapolation outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setExtrap.piwilin.extrap.txt', reshape([queryx, extrap], [nquery, 2_IK]), header = SK_'queryx,extrap')) error stop 'extrapolation outputting failed.'
        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!1D single polynomial extrapolation.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()
        call disp%show("if (0 /= getErrTableWrite(SK_'setExtrap.monopol.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'extrapolation node outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setExtrap.monopol.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'extrapolation node outputting failed.'
        call disp%show("call setExtrap(monopol, crdx, func, queryx, extrap, relerr)")
                        call setExtrap(monopol, crdx, func, queryx, extrap, relerr)
        call disp%show("extrap")
        call disp%show( extrap )
        call disp%show("relerr")
        call disp%show( relerr )
        call disp%show("if (0 /= getErrTableWrite(SK_'setExtrap.monopol.extrap.txt', reshape([queryx, extrap], [nquery, 2_IK]), header = SK_'queryx,extrap')) error stop 'extrapolation outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setExtrap.monopol.extrap.txt', reshape([queryx, extrap], [nquery, 2_IK]), header = SK_'queryx,extrap')) error stop 'extrapolation outputting failed.'
        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!1D single polynomial extrapolation showing Runge effect.")
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
            call disp%show("if (0 /= getErrTableWrite(SK_'setExtrap.rungeEffect.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'extrapolation node outputting failed.'")
                            if (0 /= getErrTableWrite(SK_'setExtrap.rungeEffect.crd.txt', reshape([crdx, func], [ncrd, 2_IK]), header = SK_'crdx,func')) error stop 'extrapolation node outputting failed.'
            call disp%show("queryx = getLinSpace(1.1 * minval(crdx, 1), 1.1 * maxval(crdx, 1), nquery)")
                            queryx = getLinSpace(1.1 * minval(crdx, 1), 1.1 * maxval(crdx, 1), nquery)
            call disp%show("queryx")
            call disp%show( queryx )
            call disp%show("call setResized(extrap, size(queryx, 1, IK))")
                            call setResized(extrap, size(queryx, 1, IK))
            call disp%show("call setResized(relerr, size(queryx, 1, IK))")
                            call setResized(relerr, size(queryx, 1, IK))
            call disp%show("call setExtrap(monopol, crdx, func, queryx, extrap, relerr)")
                            call setExtrap(monopol, crdx, func, queryx, extrap, relerr)
            call disp%show("extrap")
            call disp%show( extrap )
            call disp%show("relerr")
            call disp%show( relerr )
            call disp%show("if (0 /= getErrTableWrite(SK_'setExtrap.rungeEffect.extrap.txt', reshape([queryx, extrap], [nquery, 2_IK]), header = SK_'queryx,extrap')) error stop 'extrapolation outputting failed.'")
                            if (0 /= getErrTableWrite(SK_'setExtrap.rungeEffect.extrap.txt', reshape([queryx, extrap], [nquery, 2_IK]), header = SK_'queryx,extrap')) error stop 'extrapolation outputting failed.'
            call disp%skip()
        end block
    end block

end program example