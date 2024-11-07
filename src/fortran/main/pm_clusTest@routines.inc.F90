!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief
!>  This file contains implementations of procedures [pm_clusTest](@ref pm_clusTest).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%
#if     mmvue_typer_ENABLED
        !%%%%%%%%%%%%%%%%%%

        integer(IK) :: iell

        if (present(range)) then
            CHECK_ASSERTION(__LINE__, 0 < range%ndim(1) .and. range%ndim(1) <= range%ndim(2), SK_"@mmvue_typer(): The condition `0 < range%ndim(1) .and. range%ndim(1) <= range%ndim(2)` must hold. range%ndim = "//getStr(range%ndim))
            CHECK_ASSERTION(__LINE__, 0 < range%nell(1) .and. range%nell(1) <= range%nell(2), SK_"@mmvue_typer(): The condition `0 < range%nell(1) .and. range%nell(1) <= range%nell(2)` must hold. range%nell = "//getStr(range%nell))
            CHECK_ASSERTION(__LINE__, 0 < range%nsam(1) .and. range%nsam(1) <= range%nsam(2), SK_"@mmvue_typer(): The condition `0 < range%nsam(1) .and. range%nsam(1) <= range%nsam(2)` must hold. range%nsam = "//getStr(range%nsam))
            CHECK_ASSERTION(__LINE__, 0 < range%std(1) .and. range%std(1) <= range%std(2), SK_"@mmvue_typer(): The condition `0 < range%std(1) .and. range%std(1) <= range%std(2)` must hold. range%std = "//getStr(range%std))
            self%range = range
        end if

        if (self%range%ndim(1) /= self%range%ndim(2)) then; call setUnifRand(rng, self%ndim, self%range%ndim(1), self%range%ndim(2)); else; self%ndim = self%range%ndim(1); end if
        if (self%range%nell(1) /= self%range%nell(2)) then; call setUnifRand(rng, self%nell, self%range%nell(1), self%range%nell(2)); else; self%nell = self%range%nell(1); end if
        if (self%range%nsam(1) /= self%range%nsam(2)) then; call setUnifRand(rng, self%nsam, self%range%nsam(1), self%range%nsam(2)); else; self%nsam = self%range%nsam(1); end if

        call setResized(self%mean, [self%ndim, self%nell])
        call setUnifRand(rng, self%mean, self%range%mean(1), self%range%mean(2))

        call setResized(self%std, [self%ndim, self%nell])
        call setUnifRand(rng, self%std, self%range%std(1), self%range%std(2))

        if (present(nsim)) self%nsim = nsim

        !!!!
        !!!! Generate random Gramian matrices, compute the corresponding lower Cholesky factors, and invervse Gramians.
        !!!!

        call setResized(self%choLowGramUpp, [self%ndim, self%ndim + 1, self%nell])
        call setResized(self%invGram, [self%ndim, self%ndim, self%nell])
        do iell = 1, self%nell
            call setCovRand(rng, self%choLowGramUpp(:, 2 : self%ndim + 1, iell), scale = self%std(1 : self%ndim, iell))
            call setMatChol(self%choLowGramUpp(:, 2 : self%ndim + 1, iell), subset = uppDia, info = self%err%stat, chol = self%choLowGramUpp(:, 1 : self%ndim, iell), operation = transHerm)
            if (self%err%stat /= 0_IK) error stop MODULE_NAME//SK_"@mmuve_typer: setMatChol() failed."
            call setMatInv(self%invGram(:, 1 : self%ndim, iell), self%choLowGramUpp(:, 1 : self%ndim, iell), auxil = choLow)
        end do

        !!!!
        !!!! Generate random sample from the random ellipsoids.
        !!!!

        call setResized(self%invmul, self%nsam)
        call setResized(self%membership, self%nsam)
        call setResized(self%logVolNormed, self%nell)
        call setResized(self%cumPropVolNormed, self%nell)
        call setResized(self%mahalSq, [self%nell, self%nsam])
        call setResized(self%sample, [self%ndim, self%nsam])
        call setUnifEllsRand( rng & ! LCOV_EXCL_LINE
                            , rand = self%sample & ! LCOV_EXCL_LINE
                            , mahalSq = self%mahalSq & ! LCOV_EXCL_LINE
                            , invmul = self%invmul & ! LCOV_EXCL_LINE
                            , membership = self%membership & ! LCOV_EXCL_LINE
                            , mean = self%mean & ! LCOV_EXCL_LINE
                            , chol = self%choLowGramUpp & ! LCOV_EXCL_LINE
                            , subset = lowDia & ! LCOV_EXCL_LINE
                            , invGram = self%invGram & ! LCOV_EXCL_LINE
                            )
        self%size = count(self%mahalSq <= 1._RKG, dim = 2, kind = IK)
        call setLogVolUnitBall(self%logVolUnitBall, self%ndim)

        do  iell = 1, self%nell
            self%logVolNormed(iell) = getMatMulTraceLog(self%choLowGramUpp(:, 1 : self%ndim, iell))
        end do
        call setCumPropExp(self%cumPropVolNormed, self%logVolNormed, maxval(self%logVolNormed), control = sequence)

        !!!!
        !!!! Compute the effective sum of the volumes of all ellipsoids by taking the overlaps into account.
        !!!!

        self%logSumVolNormedEff = -getUnifEllsLogPDF(rng, self%mean, chol = self%choLowGramUpp, subset = lowDia, invGram = self%invGram, nsim = self%nsim, normed = .false._LK)

        !!!!
        !!!! Compute the effective volume of a single sample.
        !!!!

        do  iell = 1, self%nell
            self%logDensity = log(real(self%size(iell), RKG)) - self%logSumVolNormedEff
        end do

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif