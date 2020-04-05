!**********************************************************************************************************************************
!**********************************************************************************************************************************
!
!  ParaMonte: plain powerful parallel Monte Carlo library.
!
!  Copyright (C) 2012-present, The Computational Data Science Lab
!
!  This file is part of the ParaMonte library. 
!
!  ParaMonte is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, either version 3 of the License.
!
!  ParaMonte is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with ParaMonte.  If not, see <https://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!**********************************************************************************************************************************

module JaggedArray_mod

    use Constants_mod, only: IK, RK
    implicit none

    type :: IntVec_type
        integer(IK)     , allocatable   :: Vector(:)
    end type IntVec_type

    type :: RealVec_type
        real(RK)        , allocatable   :: Vector(:)
    end type RealVec_type

    type :: CharVec_type
        character(:)    , allocatable   :: record
    end type CharVec_type

end module JaggedArray_mod