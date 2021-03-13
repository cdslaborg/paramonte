!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!
!!!!   MIT License
!!!!
!!!!   ParaMonte: plain powerful parallel Monte Carlo library.
!!!!
!!!!   Copyright (C) 2012-present, The Computational Data Science Lab
!!!!
!!!!   This file is part of the ParaMonte library.
!!!!
!!!!   Permission is hereby granted, free of charge, to any person obtaining a
!!!!   copy of this software and associated documentation files (the "Software"),
!!!!   to deal in the Software without restriction, including without limitation
!!!!   the rights to use, copy, modify, merge, publish, distribute, sublicense,
!!!!   and/or sell copies of the Software, and to permit persons to whom the
!!!!   Software is furnished to do so, subject to the following conditions:
!!!!
!!!!   The above copyright notice and this permission notice shall be
!!!!   included in all copies or substantial portions of the Software.
!!!!
!!!!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!!!!   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
!!!!   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
!!!!   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
!!!!   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
!!!!   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
!!!!   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!!!!
!!!!   ACKNOWLEDGMENT
!!!!
!!!!   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
!!!!   As per the ParaMonte library license agreement terms, if you use any parts of
!!!!   this library for any purposes, kindly acknowledge the use of ParaMonte in your
!!!!   work (education/research/industry/development/...) by citing the ParaMonte
!!!!   library as described on this page:
!!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \brief
!> This file is part of PartitionOptDen_mod.f90 source file.
!> \author
!> Amir Shahmoradi
!> Tuesday 2:13 am, March 10, 2021, Dallas, TX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        if (Partition%nt > 1) then

                            itmin = 0
                            potential = POSINF_RK
                            do kt = 1, Partition%nt
                                call runKmeans  ( nd = nd & ! LCOV_EXCL_LINE
                                                , np = Partition%Try(it)%np & ! LCOV_EXCL_LINE
                                                , nc = Partition%Try(it)%nc & ! LCOV_EXCL_LINE
                                                , minSize = Partition%minSize & ! LCOV_EXCL_LINE
                                                , relTol = Partition%kmeansRelTol & ! LCOV_EXCL_LINE
                                                , nfailMax = Partition%maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                                , niterMax = Partition%maxAllowedKmeansRecursion & ! LCOV_EXCL_LINE
#if defined STAN_ENABLED
                                                , Point = NormedPoint(1:nd,Partition%Try(it)%nps:Partition%Try(it)%npe) & ! LCOV_EXCL_LINE
#else
                                                , Point = Point(1:nd,Partition%Try(it)%nps:Partition%Try(it)%npe) & ! LCOV_EXCL_LINE
#endif
                                                , Membership = KmeansTry(kt)%Membership(Partition%Try(it)%nps:Partition%Try(it)%npe) & ! LCOV_EXCL_LINE
                                                , MinDistanceSq = MinDistanceSq(Partition%Try(it)%nps:Partition%Try(it)%npe) & ! LCOV_EXCL_LINE
                                                , Center = KmeansTry(kt)%Center(1:nd,1:Partition%Try(it)%nc) & ! LCOV_EXCL_LINE
                                                , Size = KmeansTry(kt)%Size(1:Partition%Try(it)%nc) & ! LCOV_EXCL_LINE
                                                , potential = KmeansTry(kt)%potential & ! LCOV_EXCL_LINE
                                                , Err = KmeansTry(kt)%Err & ! LCOV_EXCL_LINE
                                                , niter = niter & ! LCOV_EXCL_LINE
                                                , nfail = nfail & ! LCOV_EXCL_LINE
                                                )
                                if (KmeansTry(kt)%Err%stat/= 2 .and. KmeansTry(kt)%potential < potential) then
                                    potential = KmeansTry(kt)%potential
                                    itmin = kt
                                end if
                            end do

                            if (itmin > 0) then
                                if (KmeansTry(itmin)%Err%occurred) then
                                    error stop ! xxx must be fixed
                                    neopt = 0 ! LCOV_EXCL_LINE
                                    return ! LCOV_EXCL_LINE
                                else
                                    Partition%Try(it)%Size(1:Partition%Try(it)%nc) = KmeansTry(itmin)%Size(1:Partition%Try(it)%nc)
                                    Partition%Try(it)%Center(1:nd,1:Partition%Try(it)%nc) = KmeansTry(itmin)%Center(1:nd,1:Partition%Try(it)%nc)
                                    Partition%Try(it)%Membership(Partition%Try(it)%nps:Partition%Try(it)%npe) = KmeansTry(itmin)%Membership(Partition%Try(it)%nps:Partition%Try(it)%npe)
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
                                            , np = Partition%Try(it)%np & ! LCOV_EXCL_LINE
                                            , nc = Partition%Try(it)%nc & ! LCOV_EXCL_LINE
                                            , minSize = Partition%minSize & ! LCOV_EXCL_LINE
                                            , relTol = Partition%kmeansRelTol & ! LCOV_EXCL_LINE
                                            , nfailMax = Partition%maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                            , niterMax = Partition%maxAllowedKmeansRecursion & ! LCOV_EXCL_LINE
#if defined STAN_ENABLED
                                            , Point = NormedPoint(1:nd,Partition%Try(it)%nps:Partition%Try(it)%npe) & ! LCOV_EXCL_LINE
#else
                                            , Point = Point(1:nd,Partition%Try(it)%nps:Partition%Try(it)%npe) & ! LCOV_EXCL_LINE
#endif
                                            , Membership = Partition%Try(it)%Membership(Partition%Try(it)%nps:Partition%Try(it)%npe) & ! LCOV_EXCL_LINE
                                            , MinDistanceSq = MinDistanceSq(Partition%Try(it)%nps:Partition%Try(it)%npe) & ! LCOV_EXCL_LINE
                                            , Center = Partition%Try(it)%Center(1:nd,1:Partition%Try(it)%nc) & ! LCOV_EXCL_LINE
                                            , Size = Partition%Try(it)%Size(1:Partition%Try(it)%nc) & ! LCOV_EXCL_LINE
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
