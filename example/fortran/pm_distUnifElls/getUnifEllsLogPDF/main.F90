program example

    use pm_kind, only: SK
    use pm_kind, only: IK, LK
    use pm_io, only: display_type
    use pm_arrayResize, only: setResized
    use pm_distUnif, only: getUnifRand, rngf_type
    use pm_matrixChol, only: getMatChol, uppDia, lowDia
    use pm_distUnifElls, only: getUnifEllsLogPDF
    use pm_distChol, only: getCholRand, uppDia
    use pm_matrixInv, only: getMatInv, choUpp
    use pm_io, only: getErrTableWrite, trans

    implicit none

    integer(IK) :: ndim, iell, nell, nsim

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: TKG => RKD
        real(TKG), allocatable :: logPDF, mean(:,:), chol(:, :, :), invGram(:, :, :)

        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("! Random vectors from Multiple Multivariate Uniform Balls with particular means.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()

        call disp%skip()
        call disp%show("ndim = 2; nell = getUnifRand(5, 10); nsim = 2000")
                        ndim = 2; nell = getUnifRand(5, 10); nsim = 2000
        call disp%show("[ndim, nell, nsim]")
        call disp%show( [ndim, nell, nsim] )
        call disp%show("mean = getUnifRand(-5._TKG, +5._TKG, ndim, nell)")
                        mean = getUnifRand(-5._TKG, +5._TKG, ndim, nell)
        call disp%show("mean")
        call disp%show( mean )
        call disp%show("logPDF = getUnifEllsLogPDF(rngf_type(), mean)")
                        logPDF = getUnifEllsLogPDF(rngf_type(), mean)
        call disp%show("logPDF")
        call disp%show( logPDF )
        call disp%skip()

        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("! Random vectors from Multiple Multivariate Uniform Ellipsoids with particular means.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()

        call disp%show("call setResized(chol, [ndim, ndim, nell])")
                        call setResized(chol, [ndim, ndim, nell])
        call disp%show("call setResized(invGram, [ndim, ndim, nell])")
                        call setResized(invGram, [ndim, ndim, nell])
        call disp%show("do iell = 1, nell; chol(:, :, iell) = getCholRand(0._TKG, ndim, uppDia); invGram(:, :, iell) = getMatInv(chol(:, :, iell), choUpp); end do")
                        do iell = 1, nell; chol(:, :, iell) = getCholRand(0._TKG, ndim, uppDia); invGram(:, :, iell) = getMatInv(chol(:, :, iell), choUpp); end do
        call disp%show("logPDF = getUnifEllsLogPDF(rngf_type(), mean, chol, uppDia, invGram)")
                        logPDF = getUnifEllsLogPDF(rngf_type(), mean, chol, uppDia, invGram)
        call disp%show("logPDF")
        call disp%show( logPDF )
        call disp%skip()

    end block

end program example