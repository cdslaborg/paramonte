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
!>  This include file contains implementations of the procedures in module [pm_kind](@ref pm_kind).
!>
!>  \finmain
!>
!>  \author
!>  \FatemehBagheri, Sunday 11:23 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%
#if     modeli_typer_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        model%digits        = digits(mold)
        model%huge          = huge(mold)
        model%kind          = kind(mold)
        model%range         = range(mold)
        model%radix         = radix(mold)
        model%storage_size  = storage_size(mold)

        !%%%%%%%%%%%%%%%%%%%%%%%
#elif   modelb_typer_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        model%modeli_type = modeli_type(mold)
        model%bit_size = bit_size(mold)

        !%%%%%%%%%%%%%%%%%%%%%%%
#elif   modelr_typer_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        model%digits        = digits(mold)
        model%epsilon       = epsilon(mold)
        model%huge          = huge(mold)
        model%kind          = kind(mold)
        model%maxexponent   = maxexponent(mold)
        model%minexponent   = minexponent(mold)
        model%precision     = precision(mold)
        model%range         = range(mold)
        model%radix         = radix(mold)
        model%storage_size  = storage_size(mold)
        model%tiny          = tiny(mold)
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif