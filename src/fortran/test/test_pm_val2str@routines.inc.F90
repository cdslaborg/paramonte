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
!>  This module contains implementations of the tests of the procedures under the generic interface [getStr](@ref pm_val2str::getStr).
!>
!>  \fintest
!>
!>  \author
!>  \AmirShahmoradi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !use pm_distUnif, only: setUnifRand

        character(:, SK), allocatable   :: str
#if     getStr_SK_ENABLED
        character(3, SK)    :: val(3,3)
#elif   getStr_RK_ENABLED
        real(RK)            :: val(3,3)
#elif   getStr_IK_ENABLED
        integer(IKC)        :: val(3,3)
#elif   getStr_CK_ENABLED
        complex(CK)         :: val(3,3)
#elif   getStr_LK_ENABLED
        logical(LK)         :: val(3,3)
#else
#error "Unrecognized Interface."
#endif
        !call setUnifRand(val)

        assertion = .true._LK

        call runTestWith()
        call runTestWith(signed = .true._LK)
        call runTestWith(signed = .false._LK)
        call runTestWith(signed = .true._LK, format = SK_"(*(g0,:,','))")
        call runTestWith(signed = .false._LK, format = SK_"(*(g0,:,','))")
        call runTestWith(len = 10000_IK, signed = .true._LK)
        call runTestWith(len = 10000_IK, signed = .false._LK)
        call runTestWith(len = 10000_IK, signed = .true._LK, format = SK_"(*(g0,:,','))")
        call runTestWith(len = 10000_IK, signed = .false._LK, format = SK_"(*(g0,:,','))")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTestWith(format, len, signed)

            character(*, SK), intent(in), optional :: format
            logical(LK)     , intent(in), optional :: signed
            integer(IK)     , intent(in), optional :: len

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         getStr_SK_ENABLED
            val = SK_"ParaMonte"
#elif       getStr_RK_ENABLED
            val = 1._RK
#elif       getStr_IK_ENABLED
            val = 1_IKC
#elif       getStr_CK_ENABLED
            val = 1._CK
#elif       getStr_LK_ENABLED
            val = .false._LK
#else
#error      "Unrecognized Interface."
#endif

            str = getStr(val(1,1), format, len, signed)
            call test%assert(assertion, SK_"getStr() must properly convert an input scalar value to string with present(format), present(len), present(signed) = "//getStr([present(format), present(len), present(signed)]))
            call checkLen(len)

            str = getStr(val(:,1), format, len, signed)
            call test%assert(assertion, SK_"getStr() must properly convert an input vector value to string with present(format), present(len), present(signed) = "//getStr([present(format), present(len), present(signed)]))
            call checkLen(len)

            str = getStr(val(:,:), format, len, signed)
            call test%assert(assertion, SK_"getStr() must properly convert an input matrix value to string with present(format), present(len), present(signed) = "//getStr([present(format), present(len), present(signed)]))
            call checkLen(len)

        end subroutine
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine checkLen(length)
            integer(IK)     , intent(in), optional :: length
            if (present(length)) then
                assertion = assertion .and. length == len(str, kind = IK)
                call test%assert(assertion, SK_"getStr() must properly must properly set the length of the output string: len, len(str) = "//getStr([length,len(str,IK)]))
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
