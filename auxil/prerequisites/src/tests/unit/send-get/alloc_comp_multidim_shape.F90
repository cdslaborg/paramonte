! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! * Redistributions of source code must retain the above copyright notice, this
!   list of conditions and the following disclaimer.
!
! * Redistributions in binary form must reproduce the above copyright notice,
!   this list of conditions and the following disclaimer in the documentation
!   and/or other materials provided with the distribution.
!
! * Neither the name of the copyright holder nor the names of its
!   contributors may be used to endorse or promote products derived from
!   this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

! Comments preceded by "!!" are formatted for the FOR_codim1 docoumentation generator
program alloc_comp_multidim_shape
  !! summary: Test shape of multidimensional allocatable array coarray components of a
  !!          derived type shape(object%comp(:,:)[1])
  !! author: _codim1amian Rouson , 2018
  !! date: 2018-03-08
  !!
  !! [OpenCoarrays issue #511](https://github.com/sourceryinstitute/opencoarrays/issues/511)

  implicit none

  ! TO_DO: add tests for other types and kinds, including integer, logical, character, and derived types

  type reals
    ! Array coarray component with rank + corank = 2
    real, allocatable :: dim01_codim01(:)[:]
    ! Array coarray components with rank + corank = 3
    real, allocatable :: dim01_codim02(:)[:,:]
    real, allocatable :: dim02_codim01(:,:)[:]
    ! Array coarray components with rank + corank = 4
    real, allocatable :: dim01_codim03(:)[:,:,:]
    real, allocatable :: dim02_codim02(:,:)[:,:]
    real, allocatable :: dim03_codim01(:,:,:)[:]
    ! Array coarray components with rank + corank = 5
    real, allocatable :: dim01_codim04(:)[:,:,:,:]
    real, allocatable :: dim02_codim03(:,:)[:,:,:]
    real, allocatable :: dim03_codim02(:,:,:)[:,:]
    real, allocatable :: dim04_codim01(:,:,:,:)[:]
    ! Array coarray components with rank + corank = 6
    real, allocatable :: dim01_codim05(:)[:,:,:,:,:]
    real, allocatable :: dim02_codim04(:,:)[:,:,:,:]
    real, allocatable :: dim03_codim03(:,:,:)[:,:,:]
    real, allocatable :: dim04_codim02(:,:,:,:)[:,:]
    real, allocatable :: dim05_codim01(:,:,:,:,:)[:]
    ! Array coarray components with rank + corank = 7
    real, allocatable :: dim01_codim06(:)[:,:,:,:,:,:]
    real, allocatable :: dim02_codim05(:,:)[:,:,:,:,:]
    real, allocatable :: dim03_codim04(:,:,:)[:,:,:,:]
    real, allocatable :: dim04_codim03(:,:,:,:)[:,:,:]
    real, allocatable :: dim05_codim02(:,:,:,:,:)[:,:]
    real, allocatable :: dim06_codim01(:,:,:,:,:,:)[:]
#ifdef GCC_GE_8
    ! Array coarray components with rank + corank = 8
    real, allocatable :: dim01_codim07(:)[:,:,:,:,:,:,:]
    real, allocatable :: dim02_codim06(:,:)[:,:,:,:,:,:]
    real, allocatable :: dim03_codim05(:,:,:)[:,:,:,:,:]
    real, allocatable :: dim04_codim04(:,:,:,:)[:,:,:,:]
    real, allocatable :: dim05_codim03(:,:,:,:,:)[:,:,:]
    real, allocatable :: dim06_codim02(:,:,:,:,:,:)[:,:]
    real, allocatable :: dim07_codim01(:,:,:,:,:,:,:)[:]
    ! Array coarray components with rank + corank = 9
    real, allocatable :: dim01_codim08(:)[:,:,:,:,:,:,:,:]
    real, allocatable :: dim02_codim07(:,:)[:,:,:,:,:,:,:]
    real, allocatable :: dim03_codim06(:,:,:)[:,:,:,:,:,:]
    real, allocatable :: dim04_codim05(:,:,:,:)[:,:,:,:,:]
    real, allocatable :: dim05_codim04(:,:,:,:,:)[:,:,:,:]
    real, allocatable :: dim06_codim03(:,:,:,:,:,:)[:,:,:]
    real, allocatable :: dim07_codim02(:,:,:,:,:,:,:)[:,:]
    real, allocatable :: dim08_codim01(:,:,:,:,:,:,:,:)[:]
    ! Array coarray components with rank + corank = 10
    real, allocatable :: dim01_codim09(:)[:,:,:,:,:,:,:,:,:]
    real, allocatable :: dim02_codim08(:,:)[:,:,:,:,:,:,:,:]
    real, allocatable :: dim03_codim07(:,:,:)[:,:,:,:,:,:,:]
    real, allocatable :: dim04_codim06(:,:,:,:)[:,:,:,:,:,:]
    real, allocatable :: dim05_codim05(:,:,:,:,:)[:,:,:,:,:]
    real, allocatable :: dim06_codim04(:,:,:,:,:,:)[:,:,:,:]
    real, allocatable :: dim07_codim03(:,:,:,:,:,:,:)[:,:,:]
    real, allocatable :: dim08_codim02(:,:,:,:,:,:,:,:)[:,:]
    real, allocatable :: dim09_codim01(:,:,:,:,:,:,:,:,:)[:]
    ! Array coarray components with rank + corank = 11
    real, allocatable :: dim01_codim10(:)[:,:,:,:,:,:,:,:,:,:]
    real, allocatable :: dim02_codim09(:,:)[:,:,:,:,:,:,:,:,:]
    real, allocatable :: dim03_codim08(:,:,:)[:,:,:,:,:,:,:,:]
    real, allocatable :: dim04_codim07(:,:,:,:)[:,:,:,:,:,:,:]
    real, allocatable :: dim05_codim06(:,:,:,:,:)[:,:,:,:,:,:]
    real, allocatable :: dim06_codim05(:,:,:,:,:,:)[:,:,:,:,:]
    real, allocatable :: dim07_codim04(:,:,:,:,:,:,:)[:,:,:,:]
    real, allocatable :: dim08_codim03(:,:,:,:,:,:,:,:)[:,:,:]
    real, allocatable :: dim09_codim02(:,:,:,:,:,:,:,:,:)[:,:]
    real, allocatable :: dim10_codim01(:,:,:,:,:,:,:,:,:,:)[:]
    ! Array coarray components with rank + corank = 12
    real, allocatable :: dim01_codim11(:)[:,:,:,:,:,:,:,:,:,:,:]
    real, allocatable :: dim02_codim10(:,:)[:,:,:,:,:,:,:,:,:,:]
    real, allocatable :: dim03_codim09(:,:,:)[:,:,:,:,:,:,:,:,:]
    real, allocatable :: dim04_codim08(:,:,:,:)[:,:,:,:,:,:,:,:]
    real, allocatable :: dim05_codim07(:,:,:,:,:)[:,:,:,:,:,:,:]
    real, allocatable :: dim06_codim06(:,:,:,:,:,:)[:,:,:,:,:,:]
    real, allocatable :: dim07_codim05(:,:,:,:,:,:,:)[:,:,:,:,:]
    real, allocatable :: dim08_codim04(:,:,:,:,:,:,:,:)[:,:,:,:]
    real, allocatable :: dim09_codim03(:,:,:,:,:,:,:,:,:)[:,:,:]
    real, allocatable :: dim10_codim02(:,:,:,:,:,:,:,:,:,:)[:,:]
    real, allocatable :: dim11_codim01(:,:,:,:,:,:,:,:,:,:,:)[:]
    ! Array coarray components with rank + corank = 13
    real, allocatable :: dim01_codim12(:)[:,:,:,:,:,:,:,:,:,:,:,:]
    real, allocatable :: dim02_codim11(:,:)[:,:,:,:,:,:,:,:,:,:,:]
    real, allocatable :: dim03_codim10(:,:,:)[:,:,:,:,:,:,:,:,:,:]
    real, allocatable :: dim04_codim09(:,:,:,:)[:,:,:,:,:,:,:,:,:]
    real, allocatable :: dim05_codim08(:,:,:,:,:)[:,:,:,:,:,:,:,:]
    real, allocatable :: dim06_codim07(:,:,:,:,:,:)[:,:,:,:,:,:,:]
    real, allocatable :: dim07_codim06(:,:,:,:,:,:,:)[:,:,:,:,:,:]
    real, allocatable :: dim08_codim05(:,:,:,:,:,:,:,:)[:,:,:,:,:]
    real, allocatable :: dim09_codim04(:,:,:,:,:,:,:,:,:)[:,:,:,:]
    real, allocatable :: dim10_codim03(:,:,:,:,:,:,:,:,:,:)[:,:,:]
    real, allocatable :: dim11_codim02(:,:,:,:,:,:,:,:,:,:,:)[:,:]
    real, allocatable :: dim12_codim01(:,:,:,:,:,:,:,:,:,:,:,:)[:]
    ! Array coarray components with rank + corank = 14
    real, allocatable :: dim01_codim13(:)[:,:,:,:,:,:,:,:,:,:,:,:,:]
    real, allocatable :: dim02_codim12(:,:)[:,:,:,:,:,:,:,:,:,:,:,:]
    real, allocatable :: dim03_codim11(:,:,:)[:,:,:,:,:,:,:,:,:,:,:]
    real, allocatable :: dim04_codim10(:,:,:,:)[:,:,:,:,:,:,:,:,:,:]
    real, allocatable :: dim05_codim09(:,:,:,:,:)[:,:,:,:,:,:,:,:,:]
    real, allocatable :: dim06_codim08(:,:,:,:,:,:)[:,:,:,:,:,:,:,:]
    real, allocatable :: dim07_codim07(:,:,:,:,:,:,:)[:,:,:,:,:,:,:]
    real, allocatable :: dim08_codim06(:,:,:,:,:,:,:,:)[:,:,:,:,:,:]
    real, allocatable :: dim09_codim05(:,:,:,:,:,:,:,:,:)[:,:,:,:,:]
    real, allocatable :: dim10_codim04(:,:,:,:,:,:,:,:,:,:)[:,:,:,:]
    real, allocatable :: dim11_codim03(:,:,:,:,:,:,:,:,:,:,:)[:,:,:]
    real, allocatable :: dim12_codim02(:,:,:,:,:,:,:,:,:,:,:,:)[:,:]
    real, allocatable :: dim13_codim01(:,:,:,:,:,:,:,:,:,:,:,:,:)[:]
    ! Array coarray components with rank + corank = 15
    real, allocatable :: dim01_codim14(:)[:,:,:,:,:,:,:,:,:,:,:,:,:,:]
    real, allocatable :: dim02_codim13(:,:)[:,:,:,:,:,:,:,:,:,:,:,:,:]
    real, allocatable :: dim03_codim12(:,:,:)[:,:,:,:,:,:,:,:,:,:,:,:]
    real, allocatable :: dim04_codim11(:,:,:,:)[:,:,:,:,:,:,:,:,:,:,:]
    real, allocatable :: dim05_codim10(:,:,:,:,:)[:,:,:,:,:,:,:,:,:,:]
    real, allocatable :: dim06_codim09(:,:,:,:,:,:)[:,:,:,:,:,:,:,:,:]
    real, allocatable :: dim07_codim08(:,:,:,:,:,:,:)[:,:,:,:,:,:,:,:]
    real, allocatable :: dim08_codim07(:,:,:,:,:,:,:,:)[:,:,:,:,:,:,:]
    real, allocatable :: dim09_codim06(:,:,:,:,:,:,:,:,:)[:,:,:,:,:,:]
    real, allocatable :: dim10_codim05(:,:,:,:,:,:,:,:,:,:)[:,:,:,:,:]
    real, allocatable :: dim11_codim04(:,:,:,:,:,:,:,:,:,:,:)[:,:,:,:]
    real, allocatable :: dim12_codim03(:,:,:,:,:,:,:,:,:,:,:,:)[:,:,:]
    real, allocatable :: dim13_codim02(:,:,:,:,:,:,:,:,:,:,:,:,:)[:,:]
    real, allocatable :: dim14_codim01(:,:,:,:,:,:,:,:,:,:,:,:,:,:)[:]
#endif
  end type

  type(reals), save :: object

  logical :: error_printed=.false.

  integer, parameter :: shape01D(*)=[2]
  integer, parameter :: shape02D(*)=[2,2]
  integer, parameter :: shape03D(*)=[2,2,2]
  integer, parameter :: shape04D(*)=[2,2,2,2]
  integer, parameter :: shape05D(*)=[2,2,2,2,2]
  integer, parameter :: shape06D(*)=[2,2,2,2,2,2]
  integer, parameter :: shape07D(*)=[2,2,2,2,2,2,2]
  integer, parameter :: shape08D(*)=[2,2,2,2,2,2,2,2]
  integer, parameter :: shape09D(*)=[2,2,2,2,2,2,2,2,2]
  integer, parameter :: shape10D(*)=[2,2,2,2,2,2,2,2,2,2]
  integer, parameter :: shape11D(*)=[2,2,2,2,2,2,2,2,2,2,2]
  integer, parameter :: shape12D(*)=[2,2,2,2,2,2,2,2,2,2,2,2]
  integer, parameter :: shape13D(*)=[2,2,2,2,2,2,2,2,2,2,2,2,2]
  integer, parameter :: shape14D(*)=[2,2,2,2,2,2,2,2,2,2,2,2,2,2]

  associate(me => this_image(), np => num_images())

  ! Array coarray component with rank + corank = 2
  allocate(object%dim01_codim01(2)[*])
  if (any( shape(object%dim01_codim01(:)[1]) /= [2] )) call print_log('misshapen dim01_codim01')

  ! Array coarray components with rank + corank = 3
  allocate(object%dim01_codim02(2)[1,*])
  allocate(object%dim02_codim01(2,2)[*])
  if (any( shape(object%dim01_codim02(:)[1,1]) /= shape01D )) call print_log('misshapen dim01_codim02')
  if (any( shape(object%dim02_codim01(:,:)[1]) /= shape02D )) call print_log('misshapen dim02_codim01')

  ! Array coarray components with rank + corank = 4
  allocate(object%dim01_codim03(2)[1,1,*])
  allocate(object%dim02_codim02(2,2)[1,*])
  allocate(object%dim03_codim01(2,2,2)[*])
  if (any( shape(object%dim01_codim03(:)[1,1,1]) /= shape01D )) call print_log('misshapen dim01_codim03')
  if (any( shape(object%dim02_codim02(:,:)[1,1]) /= shape02D )) call print_log('misshapen dim02_codim02')
  if (any( shape(object%dim03_codim01(:,:,:)[1]) /= shape03D )) call print_log('misshapen dim03_codim01')

  ! Array coarray components with rank + corank = 5
  allocate(object%dim01_codim04(2)[1,1,1,*])
  allocate(object%dim02_codim03(2,2)[1,1,*])
  allocate(object%dim03_codim02(2,2,2)[1,*])
  allocate(object%dim04_codim01(2,2,2,2)[*])
  if (any( shape(object%dim01_codim04(:)[1,1,1,1]) /= shape01D )) call print_log('misshapen dim01_codim04')
  if (any( shape(object%dim02_codim03(:,:)[1,1,1]) /= shape02D )) call print_log('misshapen dim02_codim03')
  if (any( shape(object%dim03_codim02(:,:,:)[1,1]) /= shape03D )) call print_log('misshapen dim03_codim02')
  if (any( shape(object%dim04_codim01(:,:,:,:)[1]) /= shape04D )) call print_log('misshapen dim04_codim01')

  ! Array coarray components with rank + corank = 6
  allocate(object%dim01_codim05(2)[1,1,1,1,*])
  allocate(object%dim02_codim04(2,2)[1,1,1,*])
  allocate(object%dim03_codim03(2,2,2)[1,1,*])
  allocate(object%dim04_codim02(2,2,2,2)[1,*])
  allocate(object%dim05_codim01(2,2,2,2,2)[*])
  if (any( shape(object%dim01_codim05(:)[1,1,1,1,1]) /= shape01D )) call print_log('misshapen dim01_codim05')
  if (any( shape(object%dim02_codim04(:,:)[1,1,1,1]) /= shape02D )) call print_log('misshapen dim02_codim04')
  if (any( shape(object%dim03_codim03(:,:,:)[1,1,1]) /= shape03D )) call print_log('misshapen dim03_codim03')
  if (any( shape(object%dim04_codim02(:,:,:,:)[1,1]) /= shape04D )) call print_log('misshapen dim04_codim02')
  if (any( shape(object%dim05_codim01(:,:,:,:,:)[1]) /= shape05D )) call print_log('misshapen dim05_codim01')

  ! Array coarray components with rank + corank = 7
  allocate(object%dim01_codim06(2)[1,1,1,1,1,*])
  allocate(object%dim02_codim05(2,2)[1,1,1,1,*])
  allocate(object%dim03_codim04(2,2,2)[1,1,1,*])
  allocate(object%dim04_codim03(2,2,2,2)[1,1,*])
  allocate(object%dim05_codim02(2,2,2,2,2)[1,*])
  allocate(object%dim06_codim01(2,2,2,2,2,2)[*])
  if (any( shape(object%dim01_codim06(:)[1,1,1,1,1,1]) /= shape01D )) call print_log('misshapen dim01_codim06')
  if (any( shape(object%dim02_codim05(:,:)[1,1,1,1,1]) /= shape02D )) call print_log('misshapen dim02_codim05')
  if (any( shape(object%dim03_codim04(:,:,:)[1,1,1,1]) /= shape03D )) call print_log('misshapen dim03_codim04')
  if (any( shape(object%dim04_codim03(:,:,:,:)[1,1,1]) /= shape04D )) call print_log('misshapen dim04_codim03')
  if (any( shape(object%dim05_codim02(:,:,:,:,:)[1,1]) /= shape05D )) call print_log('misshapen dim05_codim02')
  if (any( shape(object%dim06_codim01(:,:,:,:,:,:)[1]) /= shape06D )) call print_log('misshapen dim06_codim01')

#ifdef GCC_GE_8

  ! Array coarray components with rank + corank = 8
  allocate(object%dim01_codim07(2)[1,1,1,1,1,1,*])
  allocate(object%dim02_codim06(2,2)[1,1,1,1,1,*])
  allocate(object%dim03_codim05(2,2,2)[1,1,1,1,*])
  allocate(object%dim04_codim04(2,2,2,2)[1,1,1,*])
  allocate(object%dim05_codim03(2,2,2,2,2)[1,1,*])
  allocate(object%dim06_codim02(2,2,2,2,2,2)[1,*])
  allocate(object%dim07_codim01(2,2,2,2,2,2,2)[*])
  if (any( shape(object%dim01_codim07(:)[1,1,1,1,1,1,1]) /= shape01D )) call print_log('misshapen dim01_codim07')
  if (any( shape(object%dim02_codim06(:,:)[1,1,1,1,1,1]) /= shape02D )) call print_log('misshapen dim02_codim06')
  if (any( shape(object%dim03_codim05(:,:,:)[1,1,1,1,1]) /= shape03D )) call print_log('misshapen dim03_codim05')
  if (any( shape(object%dim04_codim04(:,:,:,:)[1,1,1,1]) /= shape04D )) call print_log('misshapen dim04_codim04')
  if (any( shape(object%dim05_codim03(:,:,:,:,:)[1,1,1]) /= shape05D )) call print_log('misshapen dim05_codim03')
  if (any( shape(object%dim06_codim02(:,:,:,:,:,:)[1,1]) /= shape06D )) call print_log('misshapen dim06_codim02')
  if (any( shape(object%dim07_codim01(:,:,:,:,:,:,:)[1]) /= shape07D )) call print_log('misshapen dim07_codim01')

  ! Array coarray components with rank + corank = 9
  allocate(object%dim01_codim08(2)[1,1,1,1,1,1,1,*])
  allocate(object%dim02_codim07(2,2)[1,1,1,1,1,1,*])
  allocate(object%dim03_codim06(2,2,2)[1,1,1,1,1,*])
  allocate(object%dim04_codim05(2,2,2,2)[1,1,1,1,*])
  allocate(object%dim05_codim04(2,2,2,2,2)[1,1,1,*])
  allocate(object%dim06_codim03(2,2,2,2,2,2)[1,1,*])
  allocate(object%dim07_codim02(2,2,2,2,2,2,2)[1,*])
  allocate(object%dim08_codim01(2,2,2,2,2,2,2,2)[*])
  if (any(shape(object%dim01_codim08(:)[1,1,1,1,1,1,1,1]) /= shape01D )) call print_log('misshapen dim01_codim08')
  if (any(shape(object%dim02_codim07(:,:)[1,1,1,1,1,1,1]) /= shape02D )) call print_log('misshapen dim02_codim07')
  if (any(shape(object%dim03_codim06(:,:,:)[1,1,1,1,1,1]) /= shape03D )) call print_log('misshapen dim03_codim06')
  if (any(shape(object%dim04_codim05(:,:,:,:)[1,1,1,1,1]) /= shape04D )) call print_log('misshapen dim04_codim05')
  if (any(shape(object%dim05_codim04(:,:,:,:,:)[1,1,1,1]) /= shape05D )) call print_log('misshapen dim05_codim04')
  if (any(shape(object%dim06_codim03(:,:,:,:,:,:)[1,1,1]) /= shape06D )) call print_log('misshapen dim06_codim03')
  if (any(shape(object%dim07_codim02(:,:,:,:,:,:,:)[1,1]) /= shape07D )) call print_log('misshapen dim07_codim02')
  if (any(shape(object%dim08_codim01(:,:,:,:,:,:,:,:)[1]) /= shape08D )) call print_log('misshapen dim08_codim01')

  ! Array coarray components with rank + corank = 10
  allocate(object%dim01_codim09(2)[1,1,1,1,1,1,1,1,*])
  allocate(object%dim02_codim08(2,2)[1,1,1,1,1,1,1,*])
  allocate(object%dim03_codim07(2,2,2)[1,1,1,1,1,1,*])
  allocate(object%dim04_codim06(2,2,2,2)[1,1,1,1,1,*])
  allocate(object%dim05_codim05(2,2,2,2,2)[1,1,1,1,*])
  allocate(object%dim06_codim04(2,2,2,2,2,2)[1,1,1,*])
  allocate(object%dim07_codim03(2,2,2,2,2,2,2)[1,1,*])
  allocate(object%dim08_codim02(2,2,2,2,2,2,2,2)[1,*])
  allocate(object%dim09_codim01(2,2,2,2,2,2,2,2,2)[*])
  if (any( shape(object%dim01_codim09(:)[1,1,1,1,1,1,1,1,1]) /= shape01D )) call print_log('misshapen dim01_codim09')
  if (any( shape(object%dim02_codim08(:,:)[1,1,1,1,1,1,1,1]) /= shape02D )) call print_log('misshapen dim02_codim08')
  if (any( shape(object%dim03_codim07(:,:,:)[1,1,1,1,1,1,1]) /= shape03D )) call print_log('misshapen dim03_codim07')
  if (any( shape(object%dim04_codim06(:,:,:,:)[1,1,1,1,1,1]) /= shape04D )) call print_log('misshapen dim04_codim06')
  if (any( shape(object%dim05_codim05(:,:,:,:,:)[1,1,1,1,1]) /= shape05D )) call print_log('misshapen dim05_codim05')
  if (any( shape(object%dim06_codim04(:,:,:,:,:,:)[1,1,1,1]) /= shape06D )) call print_log('misshapen dim06_codim04')
  if (any( shape(object%dim07_codim03(:,:,:,:,:,:,:)[1,1,1]) /= shape07D )) call print_log('misshapen dim07_codim03')
  if (any( shape(object%dim08_codim02(:,:,:,:,:,:,:,:)[1,1]) /= shape08D )) call print_log('misshapen dim08_codim02')
  if (any( shape(object%dim09_codim01(:,:,:,:,:,:,:,:,:)[1]) /= shape09D )) call print_log('misshapen dim09_codim01')

  ! Array coarray components with rank + corank = 11
  allocate(object%dim01_codim10(2)[1,1,1,1,1,1,1,1,1,*])
  allocate(object%dim02_codim09(2,2)[1,1,1,1,1,1,1,1,*])
  allocate(object%dim03_codim08(2,2,2)[1,1,1,1,1,1,1,*])
  allocate(object%dim04_codim07(2,2,2,2)[1,1,1,1,1,1,*])
  allocate(object%dim05_codim06(2,2,2,2,2)[1,1,1,1,1,*])
  allocate(object%dim06_codim05(2,2,2,2,2,2)[1,1,1,1,*])
  allocate(object%dim07_codim04(2,2,2,2,2,2,2)[1,1,1,*])
  allocate(object%dim08_codim03(2,2,2,2,2,2,2,2)[1,1,*])
  allocate(object%dim09_codim02(2,2,2,2,2,2,2,2,2)[1,*])
  allocate(object%dim10_codim01(2,2,2,2,2,2,2,2,2,2)[*])
  if (any( shape(object%dim01_codim10(:)[1,1,1,1,1,1,1,1,1,1]) /= shape01D )) call print_log('misshapen dim01_codim10')
  if (any( shape(object%dim02_codim09(:,:)[1,1,1,1,1,1,1,1,1]) /= shape02D )) call print_log('misshapen dim02_codim09')
  if (any( shape(object%dim03_codim08(:,:,:)[1,1,1,1,1,1,1,1]) /= shape03D )) call print_log('misshapen dim03_codim08')
  if (any( shape(object%dim04_codim07(:,:,:,:)[1,1,1,1,1,1,1]) /= shape04D )) call print_log('misshapen dim04_codim07')
  if (any( shape(object%dim05_codim06(:,:,:,:,:)[1,1,1,1,1,1]) /= shape05D )) call print_log('misshapen dim05_codim06')
  if (any( shape(object%dim06_codim05(:,:,:,:,:,:)[1,1,1,1,1]) /= shape06D )) call print_log('misshapen dim06_codim05')
  if (any( shape(object%dim07_codim04(:,:,:,:,:,:,:)[1,1,1,1]) /= shape07D )) call print_log('misshapen dim07_codim04')
  if (any( shape(object%dim08_codim03(:,:,:,:,:,:,:,:)[1,1,1]) /= shape08D )) call print_log('misshapen dim08_codim03')
  if (any( shape(object%dim09_codim02(:,:,:,:,:,:,:,:,:)[1,1]) /= shape09D )) call print_log('misshapen dim09_codim02')
  if (any( shape(object%dim10_codim01(:,:,:,:,:,:,:,:,:,:)[1]) /= shape10D )) call print_log('misshapen dim10_codim01')

  ! Array coarray components with rank + corank = 12
  allocate(object%dim01_codim11(2)[1,1,1,1,1,1,1,1,1,1,*])
  allocate(object%dim02_codim10(2,2)[1,1,1,1,1,1,1,1,1,*])
  allocate(object%dim03_codim09(2,2,2)[1,1,1,1,1,1,1,1,*])
  allocate(object%dim04_codim08(2,2,2,2)[1,1,1,1,1,1,1,*])
  allocate(object%dim05_codim07(2,2,2,2,2)[1,1,1,1,1,1,*])
  allocate(object%dim06_codim06(2,2,2,2,2,2)[1,1,1,1,1,*])
  allocate(object%dim07_codim05(2,2,2,2,2,2,2)[1,1,1,1,*])
  allocate(object%dim08_codim04(2,2,2,2,2,2,2,2)[1,1,1,*])
  allocate(object%dim09_codim03(2,2,2,2,2,2,2,2,2)[1,1,*])
  allocate(object%dim10_codim02(2,2,2,2,2,2,2,2,2,2)[1,*])
  allocate(object%dim11_codim01(2,2,2,2,2,2,2,2,2,2,2)[*])
  if (any( shape(object%dim01_codim11(:)[1,1,1,1,1,1,1,1,1,1,1]) /= shape01D )) call print_log('misshapen dim01_codim11')
  if (any( shape(object%dim02_codim10(:,:)[1,1,1,1,1,1,1,1,1,1]) /= shape02D )) call print_log('misshapen dim02_codim10')
  if (any( shape(object%dim03_codim09(:,:,:)[1,1,1,1,1,1,1,1,1]) /= shape03D )) call print_log('misshapen dim03_codim09')
  if (any( shape(object%dim04_codim08(:,:,:,:)[1,1,1,1,1,1,1,1]) /= shape04D )) call print_log('misshapen dim04_codim08')
  if (any( shape(object%dim05_codim07(:,:,:,:,:)[1,1,1,1,1,1,1]) /= shape05D )) call print_log('misshapen dim05_codim07')
  if (any( shape(object%dim06_codim06(:,:,:,:,:,:)[1,1,1,1,1,1]) /= shape06D )) call print_log('misshapen dim06_codim06')
  if (any( shape(object%dim07_codim05(:,:,:,:,:,:,:)[1,1,1,1,1]) /= shape07D )) call print_log('misshapen dim07_codim05')
  if (any( shape(object%dim08_codim04(:,:,:,:,:,:,:,:)[1,1,1,1]) /= shape08D )) call print_log('misshapen dim08_codim04')
  if (any( shape(object%dim09_codim03(:,:,:,:,:,:,:,:,:)[1,1,1]) /= shape09D )) call print_log('misshapen dim09_codim03')
  if (any( shape(object%dim10_codim02(:,:,:,:,:,:,:,:,:,:)[1,1]) /= shape10D )) call print_log('misshapen dim10_codim02')
  if (any( shape(object%dim11_codim01(:,:,:,:,:,:,:,:,:,:,:)[1]) /= shape11D )) call print_log('misshapen dim11_codim01')

  ! Array coarray components with rank + corank = 13
  allocate(object%dim01_codim12(2)[1,1,1,1,1,1,1,1,1,1,1,*])
  allocate(object%dim02_codim11(2,2)[1,1,1,1,1,1,1,1,1,1,*])
  allocate(object%dim03_codim10(2,2,2)[1,1,1,1,1,1,1,1,1,*])
  allocate(object%dim04_codim09(2,2,2,2)[1,1,1,1,1,1,1,1,*])
  allocate(object%dim05_codim08(2,2,2,2,2)[1,1,1,1,1,1,1,*])
  allocate(object%dim06_codim07(2,2,2,2,2,2)[1,1,1,1,1,1,*])
  allocate(object%dim07_codim06(2,2,2,2,2,2,2)[1,1,1,1,1,*])
  allocate(object%dim08_codim05(2,2,2,2,2,2,2,2)[1,1,1,1,*])
  allocate(object%dim09_codim04(2,2,2,2,2,2,2,2,2)[1,1,1,*])
  allocate(object%dim10_codim03(2,2,2,2,2,2,2,2,2,2)[1,1,*])
  allocate(object%dim11_codim02(2,2,2,2,2,2,2,2,2,2,2)[1,*])
  allocate(object%dim12_codim01(2,2,2,2,2,2,2,2,2,2,2,2)[*])
  if (any( shape(object%dim01_codim12(:)[1,1,1,1,1,1,1,1,1,1,1,1]) /= shape01D )) call print_log('misshapen dim01_codim12')
  if (any( shape(object%dim02_codim11(:,:)[1,1,1,1,1,1,1,1,1,1,1]) /= shape02D )) call print_log('misshapen dim02_codim11')
  if (any( shape(object%dim03_codim10(:,:,:)[1,1,1,1,1,1,1,1,1,1]) /= shape03D )) call print_log('misshapen dim03_codim10')
  if (any( shape(object%dim04_codim09(:,:,:,:)[1,1,1,1,1,1,1,1,1]) /= shape04D )) call print_log('misshapen dim04_codim09')
  if (any( shape(object%dim05_codim08(:,:,:,:,:)[1,1,1,1,1,1,1,1]) /= shape05D )) call print_log('misshapen dim05_codim08')
  if (any( shape(object%dim06_codim07(:,:,:,:,:,:)[1,1,1,1,1,1,1]) /= shape06D )) call print_log('misshapen dim06_codim07')
  if (any( shape(object%dim07_codim06(:,:,:,:,:,:,:)[1,1,1,1,1,1]) /= shape07D )) call print_log('misshapen dim07_codim06')
  if (any( shape(object%dim08_codim05(:,:,:,:,:,:,:,:)[1,1,1,1,1]) /= shape08D )) call print_log('misshapen dim08_codim05')
  if (any( shape(object%dim09_codim04(:,:,:,:,:,:,:,:,:)[1,1,1,1]) /= shape09D )) call print_log('misshapen dim09_codim04')
  if (any( shape(object%dim10_codim03(:,:,:,:,:,:,:,:,:,:)[1,1,1]) /= shape10D )) call print_log('misshapen dim10_codim03')
  if (any( shape(object%dim11_codim02(:,:,:,:,:,:,:,:,:,:,:)[1,1]) /= shape11D )) call print_log('misshapen dim11_codim02')
  if (any( shape(object%dim12_codim01(:,:,:,:,:,:,:,:,:,:,:,:)[1]) /= shape12D )) call print_log('misshapen dim12_codim01')

  ! Array coarray components with rank + corank = 14
  allocate(object%dim01_codim13(2)[1,1,1,1,1,1,1,1,1,1,1,1,*])
  allocate(object%dim02_codim12(2,2)[1,1,1,1,1,1,1,1,1,1,1,*])
  allocate(object%dim03_codim11(2,2,2)[1,1,1,1,1,1,1,1,1,1,*])
  allocate(object%dim04_codim10(2,2,2,2)[1,1,1,1,1,1,1,1,1,*])
  allocate(object%dim05_codim09(2,2,2,2,2)[1,1,1,1,1,1,1,1,*])
  allocate(object%dim06_codim08(2,2,2,2,2,2)[1,1,1,1,1,1,1,*])
  allocate(object%dim07_codim07(2,2,2,2,2,2,2)[1,1,1,1,1,1,*])
  allocate(object%dim08_codim06(2,2,2,2,2,2,2,2)[1,1,1,1,1,*])
  allocate(object%dim09_codim05(2,2,2,2,2,2,2,2,2)[1,1,1,1,*])
  allocate(object%dim10_codim04(2,2,2,2,2,2,2,2,2,2)[1,1,1,*])
  allocate(object%dim11_codim03(2,2,2,2,2,2,2,2,2,2,2)[1,1,*])
  allocate(object%dim12_codim02(2,2,2,2,2,2,2,2,2,2,2,2)[1,*])
  allocate(object%dim13_codim01(2,2,2,2,2,2,2,2,2,2,2,2,2)[*])
  if (any( shape(object%dim01_codim13(:)[1,1,1,1,1,1,1,1,1,1,1,1,1]) /= shape01D )) call print_log('misshapen dim01_codim13')
  if (any( shape(object%dim02_codim12(:,:)[1,1,1,1,1,1,1,1,1,1,1,1]) /= shape02D )) call print_log('misshapen dim02_codim12')
  if (any( shape(object%dim03_codim11(:,:,:)[1,1,1,1,1,1,1,1,1,1,1]) /= shape03D )) call print_log('misshapen dim03_codim11')
  if (any( shape(object%dim04_codim10(:,:,:,:)[1,1,1,1,1,1,1,1,1,1]) /= shape04D )) call print_log('misshapen dim04_codim10')
  if (any( shape(object%dim05_codim09(:,:,:,:,:)[1,1,1,1,1,1,1,1,1]) /= shape05D )) call print_log('misshapen dim05_codim09')
  if (any( shape(object%dim06_codim08(:,:,:,:,:,:)[1,1,1,1,1,1,1,1]) /= shape06D )) call print_log('misshapen dim06_codim08')
  if (any( shape(object%dim07_codim07(:,:,:,:,:,:,:)[1,1,1,1,1,1,1]) /= shape07D )) call print_log('misshapen dim07_codim07')
  if (any( shape(object%dim08_codim06(:,:,:,:,:,:,:,:)[1,1,1,1,1,1]) /= shape08D )) call print_log('misshapen dim08_codim06')
  if (any( shape(object%dim09_codim05(:,:,:,:,:,:,:,:,:)[1,1,1,1,1]) /= shape09D )) call print_log('misshapen dim09_codim05')
  if (any( shape(object%dim10_codim04(:,:,:,:,:,:,:,:,:,:)[1,1,1,1]) /= shape10D )) call print_log('misshapen dim10_codim04')
  if (any( shape(object%dim11_codim03(:,:,:,:,:,:,:,:,:,:,:)[1,1,1]) /= shape11D )) call print_log('misshapen dim11_codim03')
  if (any( shape(object%dim12_codim02(:,:,:,:,:,:,:,:,:,:,:,:)[1,1]) /= shape12D )) call print_log('misshapen dim12_codim02')
  if (any( shape(object%dim13_codim01(:,:,:,:,:,:,:,:,:,:,:,:,:)[1]) /= shape13D )) call print_log('misshapen dim13_codim01')

  ! Array coarray components with rank + corank = 15
  allocate(object%dim01_codim14(2)[1,1,1,1,1,1,1,1,1,1,1,1,1,*])
  allocate(object%dim02_codim13(2,2)[1,1,1,1,1,1,1,1,1,1,1,1,*])
  allocate(object%dim03_codim12(2,2,2)[1,1,1,1,1,1,1,1,1,1,1,*])
  allocate(object%dim04_codim11(2,2,2,2)[1,1,1,1,1,1,1,1,1,1,*])
  allocate(object%dim05_codim10(2,2,2,2,2)[1,1,1,1,1,1,1,1,1,*])
  allocate(object%dim06_codim09(2,2,2,2,2,2)[1,1,1,1,1,1,1,1,*])
  allocate(object%dim07_codim08(2,2,2,2,2,2,2)[1,1,1,1,1,1,1,*])
  allocate(object%dim08_codim07(2,2,2,2,2,2,2,2)[1,1,1,1,1,1,*])
  allocate(object%dim09_codim06(2,2,2,2,2,2,2,2,2)[1,1,1,1,1,*])
  allocate(object%dim10_codim05(2,2,2,2,2,2,2,2,2,2)[1,1,1,1,*])
  allocate(object%dim11_codim04(2,2,2,2,2,2,2,2,2,2,2)[1,1,1,*])
  allocate(object%dim12_codim03(2,2,2,2,2,2,2,2,2,2,2,2)[1,1,*])
  allocate(object%dim13_codim02(2,2,2,2,2,2,2,2,2,2,2,2,2)[1,*])
  allocate(object%dim14_codim01(2,2,2,2,2,2,2,2,2,2,2,2,2,2)[*])
  if (any( shape(object%dim01_codim14(:)[1,1,1,1,1,1,1,1,1,1,1,1,1,1]) /= shape01D )) call print_log('misshapen dim01_codim14')
  if (any( shape(object%dim02_codim13(:,:)[1,1,1,1,1,1,1,1,1,1,1,1,1]) /= shape02D )) call print_log('misshapen dim02_codim13')
  if (any( shape(object%dim03_codim12(:,:,:)[1,1,1,1,1,1,1,1,1,1,1,1]) /= shape03D )) call print_log('misshapen dim03_codim12')
  if (any( shape(object%dim04_codim11(:,:,:,:)[1,1,1,1,1,1,1,1,1,1,1]) /= shape04D )) call print_log('misshapen dim04_codim11')
  if (any( shape(object%dim05_codim10(:,:,:,:,:)[1,1,1,1,1,1,1,1,1,1]) /= shape05D )) call print_log('misshapen dim05_codim10')
  if (any( shape(object%dim06_codim09(:,:,:,:,:,:)[1,1,1,1,1,1,1,1,1]) /= shape06D )) call print_log('misshapen dim06_codim09')
  if (any( shape(object%dim07_codim08(:,:,:,:,:,:,:)[1,1,1,1,1,1,1,1]) /= shape07D )) call print_log('misshapen dim07_codim08')
  if (any( shape(object%dim08_codim07(:,:,:,:,:,:,:,:)[1,1,1,1,1,1,1]) /= shape08D )) call print_log('misshapen dim08_codim07')
  if (any( shape(object%dim09_codim06(:,:,:,:,:,:,:,:,:)[1,1,1,1,1,1]) /= shape09D )) call print_log('misshapen dim09_codim06')
  if (any( shape(object%dim10_codim05(:,:,:,:,:,:,:,:,:,:)[1,1,1,1,1]) /= shape10D )) call print_log('misshapen dim10_codim05')
  if (any( shape(object%dim11_codim04(:,:,:,:,:,:,:,:,:,:,:)[1,1,1,1]) /= shape11D )) call print_log('misshapen dim01_codim04')
  if (any( shape(object%dim12_codim03(:,:,:,:,:,:,:,:,:,:,:,:)[1,1,1]) /= shape12D )) call print_log('misshapen dim12_codim03')
  if (any( shape(object%dim13_codim02(:,:,:,:,:,:,:,:,:,:,:,:,:)[1,1]) /= shape13D )) call print_log('misshapen dim13_codim02')
  if (any( shape(object%dim14_codim01(:,:,:,:,:,:,:,:,:,:,:,:,:,:)[1]) /= shape14D )) call print_log('misshapen dim14_codim01')

#endif

    check_global_success: block

      logical :: no_error_printed

      no_error_printed = .not. error_printed
      call co_all(no_error_printed,result_image=1)

      if (me==1) then
         if (no_error_printed) then
          print *,"Test passed."
        else
          error stop "Errors encountered in alloc_comp_multidim_shape"
        end if
      end if

    end block check_global_success

  end associate

contains

  subroutine print_log(error_message)
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: error_message
    write(error_unit,*) error_message
    error_printed=.true.
  end subroutine

  pure function both(lhs,rhs) RESULT(lhs_and_rhs)
    logical, intent(in) :: lhs,rhs
    logical :: lhs_and_rhs
    lhs_and_rhs = lhs .and. rhs
  end function

  subroutine co_all(boolean,result_image)
    logical, intent(inout) :: boolean
    integer, intent(in) :: result_image
    call co_reduce(boolean,both,result_image=result_image)
  end subroutine

end program alloc_comp_multidim_shape
