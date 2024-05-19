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
!>  This file contains the implementation details of the routines for testing the procedures under the generic interface [pm_sampleWeight::setECDF](@ref pm_sampleWeight::setECDF).
!>
!>  \author
!>  \FatemehBagheri, Saturday 4:40 PM, August 21, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        use pm_arraySort, only: setSorted, isSorted

        real(RK)    , parameter     :: TOLERANCE = epsilon(1._RK) * 100_IK
        real(RK)    , parameter     :: SampleOrg(*) =   [ 0.706046088019609031832846377421_RK, 0.031832846377421276922984960890_RK, 0.276922984960890706046088019609_RK &
                                                        , 0.046171390631154097131781235848_RK, 0.097131781235848823457828327293_RK, 0.823457828327293046171390631154_RK &
                                                        , 0.694828622975817317099480060861_RK, 0.317099480060861950222048838355_RK, 0.950222048838355694828622975817_RK &
                                                        , 0.034446080502909438744359656398_RK, 0.438744359656398381558457093008_RK, 0.381558457093008034446080502909_RK &
                                                        , 0.765516788149002795199901137063_RK, 0.795199901137063186872604554379_RK, 0.186872604554379765516788149002_RK &
                                                        ]
        real(RK)    , parameter     :: ECDF_ref(*) =    [ 0.666666666666666666666666666666666659E-1_RK &
                                                        , 0.133333333333333333333333333333333332_RK &
                                                        , 0.200000000000000000000000000000000010_RK &
                                                        , 0.266666666666666666666666666666666663_RK &
                                                        , 0.333333333333333333333333333333333317_RK &
                                                        , 0.399999999999999999999999999999999971_RK &
                                                        , 0.466666666666666666666666666666666625_RK &
                                                        , 0.533333333333333333333333333333333327_RK &
                                                        , 0.599999999999999999999999999999999981_RK &
                                                        , 0.666666666666666666666666666666666635_RK &
                                                        , 0.733333333333333333333333333333333288_RK &
                                                        , 0.799999999999999999999999999999999942_RK &
                                                        , 0.866666666666666666666666666666666596_RK &
                                                        , 0.933333333333333333333333333333333250_RK &
                                                        , 0.999999999999999999999999999999999904_RK &
                                                        ]
        real(RK)    , parameter     :: Lower_ref(*) =   [ 0.000000000000000000000000000000000000E+0_RK &
                                                        , 0.000000000000000000000000000000000000E+0_RK &
                                                        , 0.000000000000000000000000000000000000E+0_RK &
                                                        , 0.000000000000000000000000000000000000E+0_RK &
                                                        , 0.000000000000000000000000000000000000E+0_RK &
                                                        , 4.933969647183537603134775151883088650E-2_RK &
                                                        , 0.116006363138502042698014418185497540E+0_RK &
                                                        , 0.182673029805168709364681084852164242E+0_RK &
                                                        , 0.249339696471835376031347751518830896E+0_RK &
                                                        , 0.316006363138502042698014418185497550E+0_RK &
                                                        , 0.382673029805168709364681084852164204E+0_RK &
                                                        , 0.449339696471835376031347751518830858E+0_RK &
                                                        , 0.516006363138502042698014418185497560E+0_RK &
                                                        , 0.582673029805168709364681084852164213E+0_RK &
                                                        , 0.649339696471835376031347751518830867E+0_RK &
                                                        ]
        real(RK)    , parameter     :: Upper_ref(*) =   [ 0.417326970194831290635318915147835738_RK &
                                                        , 0.483993636861497957301985581814502440_RK &
                                                        , 0.550660303528164623968652248481169094_RK &
                                                        , 0.617326970194831290635318915147835748_RK &
                                                        , 0.683993636861497957301985581814502402_RK &
                                                        , 0.750660303528164623968652248481169056_RK &
                                                        , 0.817326970194831290635318915147835710_RK &
                                                        , 0.883993636861497957301985581814502363_RK &
                                                        , 0.950660303528164623968652248481169017_RK &
                                                        , 1.000000000000000000000000000000000000_RK &
                                                        , 1.000000000000000000000000000000000000_RK &
                                                        , 1.000000000000000000000000000000000000_RK &
                                                        , 1.000000000000000000000000000000000000_RK &
                                                        , 1.000000000000000000000000000000000000_RK &
                                                        , 1.000000000000000000000000000000000000_RK &
                                                        ]
        integer(IK) , parameter     :: np = size(SampleOrg,kind=IK)
        real(RK)                    :: diff(np)
        real(RK)                    :: ECDF(np)
        real(RK)                    :: Lower(np)
        real(RK)                    :: Upper(np)
#if setECDF_D1_IK_ENABLED
        integer(IKG), parameter     :: Sample_ref(np) = int(SampleOrg, kind = IKG) * huge(0_IKG)
        integer(IKG)                :: Sample(np)
#elif setECDF_D1_RK_ENABLED
        real(RK)    , parameter     :: Sample_ref(np) = real(SampleOrg, kind = RK) * huge(0._RK)
        real(RK)                    :: Sample(np)
#else
#error "Unrecognized interface."
#endif

        Sample = Sample_ref

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! The following tests are ordered. Do not change the order of the tests.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setECDF(Sample, ECDF)
        diff = abs(ECDF - ECDF_ref)

        assertion = all(diff < TOLERANCE)
        call report(.false._LK, .false._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `ECDF` for a random input sample.")

        assertion = assertion .and. isSorted(Sample)
        call report(.false._LK, .false._LK)
        call test%assert(assertion, SK_"setECDF() must return a sorted sample given a random input sample.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setECDF(Sample, ECDF)
        diff = abs(ECDF - ECDF_ref)

        assertion = all(diff < TOLERANCE)
        call report(.false._LK, .false._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `ECDF` for a sorted input sample.")

        assertion = assertion .and. isSorted(Sample)
        call report(.false._LK, .false._LK)
        call test%assert(assertion, SK_"setECDF() must return a sorted Sample a sorted input sample.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setECDF(Sample, ECDF, sorted = .false._LK)
        diff = abs(ECDF - ECDF_ref)

        assertion = all(diff < TOLERANCE)
        call report(.false._LK, .false._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `ECDF` for a sorted input sample with `sorted = .false.`.")

        assertion = assertion .and. isSorted(Sample)
        call report(.false._LK, .false._LK)
        call test%assert(assertion, SK_"setECDF() must return a sorted Sample a sorted input sample with `sorted = .false.`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setECDF(Sample, ECDF, sorted = .true._LK)
        diff = abs(ECDF - ECDF_ref)

        assertion = all(diff < TOLERANCE)
        call report(.false._LK, .false._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `ECDF` for a sorted input sample with `sorted = .true.`.")

        assertion = assertion .and. isSorted(Sample)
        call report(.false._LK, .false._LK)
        call test%assert(assertion, SK_"setECDF() must return a sorted Sample a sorted input sample with `sorted = .true.`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Sample = Sample_ref
        call setECDF(Sample, ECDF, sorted = .false._LK)
        diff = abs(ECDF - ECDF_ref)

        assertion = all(diff < TOLERANCE)
        call report(.false._LK, .false._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `ECDF` for a random input sample with `sorted = .false.`.")

        assertion = assertion .and. isSorted(Sample)
        call report(.false._LK, .false._LK)
        call test%assert(assertion, SK_"setECDF() must return a sorted Sample a random input sample with `sorted = .false.`.")


        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Test with the presence of `Lower`
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setECDF(Sample, ECDF, Lower = Lower)

        diff = abs(Lower - Lower_ref)
        assertion = all(diff < TOLERANCE)
        call report(.true._LK, .false._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `Lower` for a random input sample with `Lower`.")

        diff = abs(ECDF - ECDF_ref)
        assertion = all(diff < TOLERANCE)
        call report(.true._LK, .false._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `ECDF` for a random input sample with `Lower`.")

        assertion = assertion .and. isSorted(Sample)
        call report(.true._LK, .false._LK)
        call test%assert(assertion, SK_"setECDF() must return a sorted sample given a random input sample with `Lower`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setECDF(Sample, ECDF, Lower = Lower)

        diff = abs(Lower - Lower_ref)
        assertion = all(diff < TOLERANCE)
        call report(.true._LK, .false._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `Lower` for a sorted input sample with `Lower`.")

        diff = abs(ECDF - ECDF_ref)
        assertion = all(diff < TOLERANCE)
        call report(.true._LK, .false._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `ECDF` for a sorted input sample with `Lower`.")

        assertion = assertion .and. isSorted(Sample)
        call report(.true._LK, .false._LK)
        call test%assert(assertion, SK_"setECDF() must return a sorted Sample a sorted input sample with `Lower`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setECDF(Sample, ECDF, Lower = Lower, sorted = .false._LK)

        diff = abs(Lower - Lower_ref)
        assertion = all(diff < TOLERANCE)
        call report(.true._LK, .false._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `Lower` for a sorted input sample with `sorted = .false.`.")

        diff = abs(ECDF - ECDF_ref)
        assertion = all(diff < TOLERANCE)
        call report(.true._LK, .false._LK)
        call test%assert(assertion, SK_"setECDF() with `Lower` must return must return the correct `ECDF` for a sorted input sample with `sorted = .false.`.")

        assertion = assertion .and. isSorted(Sample)
        call report(.true._LK, .false._LK)
        call test%assert(assertion, SK_"setECDF() with `Lower` must return a sorted Sample a sorted input sample with `sorted = .false.`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setECDF(Sample, ECDF, Lower = Lower, sorted = .true._LK)
        diff = abs(ECDF - ECDF_ref)

        diff = abs(Lower - Lower_ref)
        assertion = all(diff < TOLERANCE)
        call report(.true._LK, .false._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `Lower` for a sorted input sample with `sorted = .true.`.")

        assertion = all(diff < TOLERANCE)
        call report(.true._LK, .false._LK)
        call test%assert(assertion, SK_"setECDF() with `Lower` must return must return the correct `ECDF` for a sorted input sample with `sorted = .true.`.")

        assertion = assertion .and. isSorted(Sample)
        call report(.true._LK, .false._LK)
        call test%assert(assertion, SK_"setECDF() with `Lower` must return a sorted Sample a sorted input sample with `sorted = .true.`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Sample = Sample_ref
        call setECDF(Sample, ECDF, Lower = Lower, sorted = .false._LK)
        diff = abs(ECDF - ECDF_ref)

        diff = abs(Lower - Lower_ref)
        assertion = all(diff < TOLERANCE)
        call report(.true._LK, .false._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `Lower` for a random input sample with `sorted = .false.`.")

        assertion = all(diff < TOLERANCE)
        call report(.true._LK, .false._LK)
        call test%assert(assertion, SK_"setECDF() with `Lower` must return must return the correct `ECDF` for a random input sample with `sorted = .false.`.")

        assertion = assertion .and. isSorted(Sample)
        call report(.true._LK, .false._LK)
        call test%assert(assertion, SK_"setECDF() with `Lower` must return a sorted Sample a random input sample with `sorted = .false.`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Test with the presence of `Upper`
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setECDF(Sample, ECDF, Upper = Upper)

        diff = abs(Upper - Upper_ref)
        assertion = all(diff < TOLERANCE)
        call report(.false._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `Upper` for a random input sample with `Upper`.")

        call setECDF(Sample, ECDF, Upper = Upper, alpha = 0.05_RK)

        diff = abs(Upper - Upper_ref)
        assertion = all(diff < TOLERANCE)
        call report(.false._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `Upper` for a random input sample with `Upper` with `alpha`.")

        diff = abs(ECDF - ECDF_ref)
        assertion = all(diff < TOLERANCE)
        call report(.false._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `ECDF` for a random input sample with `Upper`.")

        assertion = assertion .and. isSorted(Sample)
        call report(.false._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return a sorted sample given a random input sample with `Upper`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setECDF(Sample, ECDF, Upper = Upper)

        diff = abs(Upper - Upper_ref)
        assertion = all(diff < TOLERANCE)
        call report(.false._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `Upper` for a sorted input sample with `Upper`.")

        call setECDF(Sample, ECDF, Upper = Upper, alpha = 0.05_RK)

        diff = abs(Upper - Upper_ref)
        assertion = all(diff < TOLERANCE)
        call report(.false._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `Upper` for a sorted input sample with `Upper` with `alpha`.")

        diff = abs(ECDF - ECDF_ref)
        assertion = all(diff < TOLERANCE)
        call report(.false._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `ECDF` for a sorted input sample with `Upper` with `alpha`.")

        assertion = assertion .and. isSorted(Sample)
        call report(.false._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return a sorted Sample a sorted input sample with `Upper` with `alpha`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setECDF(Sample, ECDF, Upper = Upper, sorted = .false._LK)

        diff = abs(Upper - Upper_ref)
        assertion = all(diff < TOLERANCE)
        call report(.false._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `Upper` for a sorted input sample with `sorted = .false.`.")

        call setECDF(Sample, ECDF, Upper = Upper, alpha = 0.05_RK, sorted = .false._LK)

        diff = abs(Upper - Upper_ref)
        assertion = all(diff < TOLERANCE)
        call report(.false._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `Upper` for a sorted input sample with `sorted = .false.` with `alpha`.")

        diff = abs(ECDF - ECDF_ref)
        assertion = all(diff < TOLERANCE)
        call report(.false._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() with `Upper` must return must return the correct `ECDF` for a sorted input sample with `sorted = .false.`.")

        assertion = assertion .and. isSorted(Sample)
        call report(.false._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() with `Upper` must return a sorted Sample a sorted input sample with `sorted = .false.`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setECDF(Sample, ECDF, Upper = Upper, sorted = .true._LK)

        diff = abs(Upper - Upper_ref)
        assertion = all(diff < TOLERANCE)
        call report(.false._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `Upper` for a sorted input sample with `sorted = .true.`.")

        call setECDF(Sample, ECDF, Upper = Upper, alpha = 0.05_RK, sorted = .true._LK)

        diff = abs(Upper - Upper_ref)
        assertion = all(diff < TOLERANCE)
        call report(.false._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `Upper` for a sorted input sample with `sorted = .true.` with `alpha`.")

        diff = abs(ECDF - ECDF_ref)
        assertion = all(diff < TOLERANCE)
        call report(.false._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `ECDF` for a sorted input sample with `sorted = .true.` with `alpha`.")

        assertion = all(diff < TOLERANCE)
        call report(.false._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() with `Upper` must return must return the correct `ECDF` for a sorted input sample with `sorted = .true.`.")

        assertion = assertion .and. isSorted(Sample)
        call report(.false._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() with `Upper` must return a sorted Sample a sorted input sample with `sorted = .true.`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Sample = Sample_ref
        call setECDF(Sample, ECDF, Upper = Upper, sorted = .false._LK)
        diff = abs(ECDF - ECDF_ref)

        diff = abs(Upper - Upper_ref)
        assertion = all(diff < TOLERANCE)
        call report(.false._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `Upper` for a random input sample with `sorted = .false.`.")

        Sample = Sample_ref
        call setECDF(Sample, ECDF, Upper = Upper, alpha = 0.05_RK, sorted = .false._LK)
        diff = abs(ECDF - ECDF_ref)

        diff = abs(Upper - Upper_ref)
        assertion = all(diff < TOLERANCE)
        call report(.false._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `Upper` for a random input sample with `sorted = .false.` with `alpha`.")

        assertion = all(diff < TOLERANCE)
        call report(.false._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() with `Upper` must return must return the correct `ECDF` for a random input sample with `sorted = .false.`.")

        assertion = assertion .and. isSorted(Sample)
        call report(.false._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() with `Upper` must return a sorted Sample a random input sample with `sorted = .false.`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Test with the presence of `Lower` and `Upper`
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setECDF(Sample, ECDF, Lower = Lower, Upper = Upper)

        diff = abs(Upper - Upper_ref)
        assertion = all(diff < TOLERANCE)
        call report(.true._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `Upper` for a random input sample with `Lower` and `Upper`.")

        call setECDF(Sample, ECDF, Lower = Lower, Upper = Upper, alpha = 0.05_RK)

        diff = abs(Upper - Upper_ref)
        assertion = all(diff < TOLERANCE)
        call report(.true._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `Upper` for a random input sample with `Lower` and `Upper` with `alpha`.")

        diff = abs(ECDF - ECDF_ref)
        assertion = all(diff < TOLERANCE)
        call report(.true._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `ECDF` for a random input sample with `Lower` and `Upper`.")

        assertion = assertion .and. isSorted(Sample)
        call report(.true._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return a sorted sample given a random input sample with `Lower` and `Upper`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setECDF(Sample, ECDF, Lower = Lower, Upper = Upper)

        diff = abs(Upper - Upper_ref)
        assertion = all(diff < TOLERANCE)
        call report(.true._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `Upper` for a sorted input sample with `Lower` and `Upper`.")

        call setECDF(Sample, ECDF, Lower = Lower, Upper = Upper, alpha = 0.05_RK)

        diff = abs(Upper - Upper_ref)
        assertion = all(diff < TOLERANCE)
        call report(.true._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `Upper` for a sorted input sample with `Lower` and `Upper` with `alpha`.")

        diff = abs(ECDF - ECDF_ref)
        assertion = all(diff < TOLERANCE)
        call report(.true._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `ECDF` for a sorted input sample with `Lower` and `Upper`.")

        assertion = assertion .and. isSorted(Sample)
        call report(.true._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return a sorted Sample a sorted input sample with `Lower` and `Upper`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setECDF(Sample, ECDF, Lower = Lower, Upper = Upper, sorted = .false._LK)

        diff = abs(Upper - Upper_ref)
        assertion = all(diff < TOLERANCE)
        call report(.true._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `Upper` for a sorted input sample with `sorted = .false.`.")

        call setECDF(Sample, ECDF, Lower = Lower, Upper = Upper, alpha = 0.05_RK, sorted = .false._LK)

        diff = abs(Upper - Upper_ref)
        assertion = all(diff < TOLERANCE)
        call report(.true._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `Upper` for a sorted input sample with `sorted = .false.` with `alpha`.")

        diff = abs(ECDF - ECDF_ref)
        assertion = all(diff < TOLERANCE)
        call report(.true._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() with `Lower` and `Upper` must return must return the correct `ECDF` for a sorted input sample with `sorted = .false.`.")

        assertion = assertion .and. isSorted(Sample)
        call report(.true._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() with `Lower` and `Upper` must return a sorted Sample a sorted input sample with `sorted = .false.`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setECDF(Sample, ECDF, Lower = Lower, Upper = Upper, sorted = .true._LK)

        diff = abs(Upper - Upper_ref)
        assertion = all(diff < TOLERANCE)
        call report(.true._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `Upper` for a sorted input sample with `sorted = .true.`.")

        call setECDF(Sample, ECDF, Lower = Lower, Upper = Upper, alpha = 0.05_RK, sorted = .true._LK)

        diff = abs(Upper - Upper_ref)
        assertion = all(diff < TOLERANCE)
        call report(.true._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `Upper` for a sorted input sample with `sorted = .true.` with `alpha`.")

        diff = abs(ECDF - ECDF_ref)
        assertion = all(diff < TOLERANCE)
        call report(.true._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() with `Lower` and `Upper` must return must return the correct `ECDF` for a sorted input sample with `sorted = .true.`.")

        assertion = assertion .and. isSorted(Sample)
        call report(.true._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() with `Lower` and `Upper` must return a sorted Sample a sorted input sample with `sorted = .true.`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Sample = Sample_ref
        call setECDF(Sample, ECDF, Lower = Lower, Upper = Upper, sorted = .false._LK)
        diff = abs(ECDF - ECDF_ref)

        diff = abs(Upper - Upper_ref)
        assertion = all(diff < TOLERANCE)
        call report(.true._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `Upper` for a random input sample with `sorted = .false.`.")

        Sample = Sample_ref
        call setECDF(Sample, ECDF, Lower = Lower, Upper = Upper, alpha = 0.05_RK, sorted = .false._LK)
        diff = abs(ECDF - ECDF_ref)

        diff = abs(Upper - Upper_ref)
        assertion = all(diff < TOLERANCE)
        call report(.true._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() must return must return the correct `Upper` for a random input sample with `sorted = .false.` with `alpha`.")

        assertion = all(diff < TOLERANCE)
        call report(.true._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() with `Lower` and `Upper` must return must return the correct `ECDF` for a random input sample with `sorted = .false.`.")

        assertion = assertion .and. isSorted(Sample)
        call report(.true._LK, .true._LK)
        call test%assert(assertion, SK_"setECDF() with `Lower` and `Upper` must return a sorted Sample a random input sample with `sorted = .false.`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(lowerEnabled, upperEnabled)
            use pm_strASCII, only: LF

            logical(LK), intent(in) :: lowerEnabled
            logical(LK), intent(in) :: upperEnabled

            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "Sample_ref             ", Sample_ref, LF
                write(test%disp%unit,"(*(g0,:,', '))") "Sample                 ", Sample, LF
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "ECDF_ref               ", ECDF_ref, LF
                write(test%disp%unit,"(*(g0,:,', '))") "ECDF                   ", ECDF, LF
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "isSorted(Sample_ref)   ", isSorted(Sample_ref)
                write(test%disp%unit,"(*(g0,:,', '))") "isSorted(Sample)       ", isSorted(Sample)
                write(test%disp%unit,"(*(g0,:,', '))")
                if (lowerEnabled) then
                write(test%disp%unit,"(*(g0,:,', '))") "Lower_ref              ", Lower_ref, LF
                write(test%disp%unit,"(*(g0,:,', '))") "Lower                  ", Lower, LF
                write(test%disp%unit,"(*(g0,:,', '))")
                end if
                if (upperEnabled) then
                write(test%disp%unit,"(*(g0,:,', '))") "Upper_ref              ", Upper_ref, LF
                write(test%disp%unit,"(*(g0,:,', '))") "Upper                  ", Upper, LF
                write(test%disp%unit,"(*(g0,:,', '))")
                end if
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "diff                   ", diff
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "TOLERANCE              ", TOLERANCE
                ! LCOV_EXCL_STOP
            end if
        end subroutine
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%