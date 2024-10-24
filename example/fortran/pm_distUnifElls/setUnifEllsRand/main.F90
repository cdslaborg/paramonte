program example

    use pm_kind, only: SK
    use pm_kind, only: IK, LK
    use pm_io, only: display_type
    use pm_arrayResize, only: setResized
    use pm_distUnif, only: getUnifRand, rngf_type
    use pm_matrixChol, only: getMatChol, uppDia, lowDia
    use pm_distUnifElls, only: setUnifEllsRand
    use pm_distChol, only: getCholRand, uppDia
    use pm_matrixInv, only: getMatInv, choUpp
    use pm_io, only: getErrTableWrite, trans

    implicit none

    integer(IK) :: ndim, iell, nell, nsam

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: TKG => RKD
        real(TKG), allocatable :: mean(:,:), invGram(:, :, :), chol(:, :, :), rand(:, :), mahalSq(:, :), invmul(:)
        integer(IK), allocatable :: membership(:)

        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("! Random vectors from Multiple Multivariate Uniform Balls with particular means.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()

        call disp%skip()
        call disp%show("ndim = 2; nell = getUnifRand(5, 20); nsam = nell * 1000")
                        ndim = 2; nell = getUnifRand(5, 20); nsam = nell * 1000
        call disp%show("[ndim, nell, nsam]")
        call disp%show( [ndim, nell, nsam] )
        call disp%show("mean = getUnifRand(-5._TKG, +5._TKG, ndim, nell)")
                        mean = getUnifRand(-5._TKG, +5._TKG, ndim, nell)
        call disp%show("mean")
        call disp%show( mean )
        call disp%show("call setResized(mahalSq, [nell, nsam])")
                        call setResized(mahalSq, [nell, nsam])
        call disp%show("call setResized(rand, [ndim, nsam])")
                        call setResized(rand, [ndim, nsam])
        call disp%show("call setResized(membership, nsam)")
                        call setResized(membership, nsam)
        call disp%show("call setResized(invmul, nsam)")
                        call setResized(invmul, nsam)
        call disp%show("call setUnifEllsRand(rngf_type(), rand, mahalSq, invmul, membership, mean)")
                        call setUnifEllsRand(rngf_type(), rand, mahalSq, invmul, membership, mean)
        call disp%show("if (0 /= getErrTableWrite(SK_'setUnifEllsRandMean.txt', reshape([transpose(rand), real(membership, TKG)], [nsam, ndim + 1_IK]))) error stop 'Failed to write the random vectors file.'")
                        if (0 /= getErrTableWrite(SK_'setUnifEllsRandMean.txt', reshape([transpose(rand), real(membership, TKG)], [nsam, ndim + 1_IK]))) error stop 'Failed to write the random vectors file.'
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
        call disp%show("chol")
        call disp%show( chol )
        call disp%show("call setUnifEllsRand(rngf_type(), rand, mahalSq, invmul, membership, mean, chol, uppDia, invGram)")
                        call setUnifEllsRand(rngf_type(), rand, mahalSq, invmul, membership, mean, chol, uppDia, invGram)
        call disp%show("if (0 /= getErrTableWrite(SK_'setUnifEllsRandMeanChol.txt', reshape([transpose(rand), real(membership, TKG)], [nsam, ndim + 1_IK]))) error stop 'Failed to write the random vectors file.'")
                        if (0 /= getErrTableWrite(SK_'setUnifEllsRandMeanChol.txt', reshape([transpose(rand), real(membership, TKG)], [nsam, ndim + 1_IK]))) error stop 'Failed to write the random vectors file.'
        call disp%skip()

    end block

end program example