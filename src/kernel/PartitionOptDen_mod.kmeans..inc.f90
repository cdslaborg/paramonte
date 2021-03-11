
                        if (Partition%nt > 1) then

                            itmin = 0
                            potential = POSINF_RK
                            do it = 1, Partition%nt
                                call runKmeans  ( nd = nd & ! LCOV_EXCL_LINE
                                                , np = ParLev(is)%np & ! LCOV_EXCL_LINE
                                                , nc = ParLev(is)%nc & ! LCOV_EXCL_LINE
                                                , minSize = Partition%minSize & ! LCOV_EXCL_LINE
                                                , relTol = Partition%kmeansRelTol & ! LCOV_EXCL_LINE
                                                , nfailMax = Partition%maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                                , niterMax = Partition%maxAllowedKmeansRecursion & ! LCOV_EXCL_LINE
#if defined STAN_ENABLED
                                                , Point = NormedPoint(1:nd,ParLev(is)%nps:ParLev(is)%npe) & ! LCOV_EXCL_LINE
#else
                                                , Point = Point(1:nd,ParLev(is)%nps:ParLev(is)%npe) & ! LCOV_EXCL_LINE
#endif
                                                , Membership = KmeansTry(it)%Membership(ParLev(is)%nps:ParLev(is)%npe) & ! LCOV_EXCL_LINE
                                                , MinDistanceSq = MinDistanceSq(ParLev(is)%nps:ParLev(is)%npe) & ! LCOV_EXCL_LINE
                                                , Center = KmeansTry(it)%Center(1:nd,1:ParLev(is)%nc) & ! LCOV_EXCL_LINE
                                                , Size = KmeansTry(it)%Size(1:ParLev(is)%nc) & ! LCOV_EXCL_LINE
                                                , potential = KmeansTry(it)%potential & ! LCOV_EXCL_LINE
                                                , Err = KmeansTry(it)%Err & ! LCOV_EXCL_LINE
                                                , niter = niter & ! LCOV_EXCL_LINE
                                                , nfail = nfail & ! LCOV_EXCL_LINE
                                                )
                                if ( KmeansTry(it)%potential < potential .and. .not. any(KmeansTry(it)%Size(1:ParLev(is)%nc) < Partition%minSize) ) then
                                    potential = KmeansTry(it)%potential
                                    itmin = it
                                end if
                            end do

                            if (itmin > 0) then
                                if (KmeansTry(itmin)%Err%occurred) then
                                    error stop ! xxx must be fixed
                                    neopt = 0 ! LCOV_EXCL_LINE
                                    return ! LCOV_EXCL_LINE
                                else
                                    ParLev(is)%Size(1:ParLev(is)%nc) = KmeansTry(itmin)%Size(1:ParLev(is)%nc)
                                    ParLev(is)%Center(1:nd,1:ParLev(is)%nc) = KmeansTry(itmin)%Center(1:nd,1:ParLev(is)%nc)
                                    ParLev(is)%Membership(ParLev(is)%nps:ParLev(is)%npe) = KmeansTry(itmin)%Membership(ParLev(is)%nps:ParLev(is)%npe)
                                end if
                            else
                                ! LCOV_EXCL_START
                                Partition%Err%msg = "All Kmeans clustering tries failed."
                                Partition%Err%occurred = .true.
                                Partition%Err%stat = 1
                                error stop ! xxx this must be fixed.
                                return
                                ! LCOV_EXCL_STOP
                            end if

                        else

                            call runKmeans  ( nd = nd & ! LCOV_EXCL_LINE
                                            , np = ParLev(is)%np & ! LCOV_EXCL_LINE
                                            , nc = ParLev(is)%nc & ! LCOV_EXCL_LINE
                                            , minSize = Partition%minSize & ! LCOV_EXCL_LINE
                                            , relTol = Partition%kmeansRelTol & ! LCOV_EXCL_LINE
                                            , nfailMax = Partition%maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                            , niterMax = Partition%maxAllowedKmeansRecursion & ! LCOV_EXCL_LINE
#if defined STAN_ENABLED
                                            , Point = NormedPoint(1:nd,ParLev(is)%nps:ParLev(is)%npe) & ! LCOV_EXCL_LINE
#else
                                            , Point = Point(1:nd,ParLev(is)%nps:ParLev(is)%npe) & ! LCOV_EXCL_LINE
#endif
                                            , Membership = ParLev(is)%Membership(ParLev(is)%nps:ParLev(is)%npe) & ! LCOV_EXCL_LINE
                                            , MinDistanceSq = MinDistanceSq(ParLev(is)%nps:ParLev(is)%npe) & ! LCOV_EXCL_LINE
                                            , Center = ParLev(is)%Center(1:nd,1:ParLev(is)%nc) & ! LCOV_EXCL_LINE
                                            , Size = ParLev(is)%Size(1:ParLev(is)%nc) & ! LCOV_EXCL_LINE
                                            , potential = potential & ! LCOV_EXCL_LINE
                                            , Err = Partition%Err & ! LCOV_EXCL_LINE
                                            , niter = niter & ! LCOV_EXCL_LINE
                                            , nfail = nfail & ! LCOV_EXCL_LINE
                                            )

                            if (Partition%Err%occurred) then
                                error stop ! xxx must be fixed
                                neopt = 0 ! LCOV_EXCL_LINE
                                return ! LCOV_EXCL_LINE
                            end if

                        end if
