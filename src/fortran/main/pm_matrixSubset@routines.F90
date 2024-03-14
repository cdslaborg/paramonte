
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
!>  This file contains procedure implementations of [pm_matrixSubset](@ref pm_matrixSubset).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_matrixSubset) routines ! LCOV_EXCL_LINE

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getSubSymmXXD; end procedure
    module procedure getSubSymmUXX; end procedure
    module procedure getSubSymmXLX; end procedure
    module procedure getSubSymmUXD; end procedure
    module procedure getSubSymmXLD; end procedure
    module procedure getSubSymmULX; end procedure
    module procedure getSubSymmULD; end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getSubCompXXD; end procedure
    module procedure getSubCompUXX; end procedure
    module procedure getSubCompXLX; end procedure
    module procedure getSubCompUXD; end procedure
    module procedure getSubCompXLD; end procedure
    module procedure getSubCompULX; end procedure
    module procedure getSubCompULD; end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getSubUnion_XXX_XXX; end procedure
    module procedure getSubUnion_XXX_UXX; end procedure
    module procedure getSubUnion_XXX_XLX; end procedure
    module procedure getSubUnion_XXX_XXD; end procedure
    module procedure getSubUnion_XXX_UXD; end procedure
    module procedure getSubUnion_XXX_XLD; end procedure
    module procedure getSubUnion_XXX_ULX; end procedure
    module procedure getSubUnion_XXX_ULD; end procedure

    module procedure getSubUnion_UXX_XXX; end procedure
    module procedure getSubUnion_UXX_UXX; end procedure
    module procedure getSubUnion_UXX_XLX; end procedure
    module procedure getSubUnion_UXX_XXD; end procedure
    module procedure getSubUnion_UXX_UXD; end procedure
    module procedure getSubUnion_UXX_XLD; end procedure
    module procedure getSubUnion_UXX_ULX; end procedure
    module procedure getSubUnion_UXX_ULD; end procedure

    module procedure getSubUnion_XLX_XXX; end procedure
    module procedure getSubUnion_XLX_UXX; end procedure
    module procedure getSubUnion_XLX_XLX; end procedure
    module procedure getSubUnion_XLX_XXD; end procedure
    module procedure getSubUnion_XLX_UXD; end procedure
    module procedure getSubUnion_XLX_XLD; end procedure
    module procedure getSubUnion_XLX_ULX; end procedure
    module procedure getSubUnion_XLX_ULD; end procedure

    module procedure getSubUnion_XXD_XXX; end procedure
    module procedure getSubUnion_XXD_UXX; end procedure
    module procedure getSubUnion_XXD_XLX; end procedure
    module procedure getSubUnion_XXD_XXD; end procedure
    module procedure getSubUnion_XXD_UXD; end procedure
    module procedure getSubUnion_XXD_XLD; end procedure
    module procedure getSubUnion_XXD_ULX; end procedure
    module procedure getSubUnion_XXD_ULD; end procedure

    module procedure getSubUnion_UXD_XXX; end procedure
    module procedure getSubUnion_UXD_UXX; end procedure
    module procedure getSubUnion_UXD_XLX; end procedure
    module procedure getSubUnion_UXD_XXD; end procedure
    module procedure getSubUnion_UXD_UXD; end procedure
    module procedure getSubUnion_UXD_XLD; end procedure
    module procedure getSubUnion_UXD_ULX; end procedure
    module procedure getSubUnion_UXD_ULD; end procedure

    module procedure getSubUnion_XLD_XXX; end procedure
    module procedure getSubUnion_XLD_UXX; end procedure
    module procedure getSubUnion_XLD_XLX; end procedure
    module procedure getSubUnion_XLD_XXD; end procedure
    module procedure getSubUnion_XLD_UXD; end procedure
    module procedure getSubUnion_XLD_XLD; end procedure
    module procedure getSubUnion_XLD_ULX; end procedure
    module procedure getSubUnion_XLD_ULD; end procedure

    module procedure getSubUnion_ULX_XXX; end procedure
    module procedure getSubUnion_ULX_UXX; end procedure
    module procedure getSubUnion_ULX_XLX; end procedure
    module procedure getSubUnion_ULX_XXD; end procedure
    module procedure getSubUnion_ULX_UXD; end procedure
    module procedure getSubUnion_ULX_XLD; end procedure
    module procedure getSubUnion_ULX_ULX; end procedure
    module procedure getSubUnion_ULX_ULD; end procedure

    module procedure getSubUnion_ULD_XXX; end procedure
    module procedure getSubUnion_ULD_UXX; end procedure
    module procedure getSubUnion_ULD_XLX; end procedure
    module procedure getSubUnion_ULD_XXD; end procedure
    module procedure getSubUnion_ULD_UXD; end procedure
    module procedure getSubUnion_ULD_XLD; end procedure
    module procedure getSubUnion_ULD_ULX; end procedure
    module procedure getSubUnion_ULD_ULD; end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines