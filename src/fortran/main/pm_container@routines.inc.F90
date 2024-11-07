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
!>  This file contains the implementation details of the routines under the generic interfaces of [pm_container](@ref pm_container).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     isless_ENABLED
        itis = allocated(con1%val) .and. allocated(con2%val)
        if (itis) itis = con1%val < con2%val
#elif   ismore_ENABLED
        itis = allocated(con1%val) .and. allocated(con2%val)
        if (itis) itis = con1%val > con2%val
#elif   isleq_ENABLED
        itis = allocated(con1%val) .and. allocated(con2%val)
        if (itis) itis = con1%val <= con2%val
#elif   ismeq_ENABLED
        itis = allocated(con1%val) .and. allocated(con2%val)
        if (itis) itis = con1%val >= con2%val
#elif   isneq_ENABLED
        itis = allocated(con1%val) .and. allocated(con2%val)
        if (itis) itis = con1%val /= con2%val
#elif   iseq_ENABLED
        itis = allocated(con1%val) .and. allocated(con2%val)
        if (itis) itis = con1%val == con2%val
#elif   assign_ENABLED
        if (allocated(source%val)) destin%val = source%val
#elif   constructCon_ENABLED && (IK_ENABLED || LK_ENABLED || CK_ENABLED || RK_ENABLED)
        container%val = val
#elif   constructCon_ENABLED && PK_ENABLED
        allocate(container%val, source = val)
#elif   constructCon_ENABLED && SK_ENABLED
        if (present(trimmed)) then
            if (trimmed) then
                container%val = val
                return
            end if
        end if
        container%val = trim(val)
#elif   getVal_ENABLED
        if (allocated(con%val)) val = con%val
#else
#error  "Unrecognized interface."
#endif