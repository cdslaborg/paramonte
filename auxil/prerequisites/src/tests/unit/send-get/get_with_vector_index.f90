program get_with_vector_index
  use iso_fortran_env
  implicit none
  integer, parameter :: nloc=8, nhl=2, ivsize=nloc+2*nhl
  real    :: xv(ivsize)[*]
  real, allocatable :: expected(:,:)
  integer :: rmt_idx(2,2), loc_idx(2,2)
  integer, allocatable :: xchg(:)
  integer :: nrcv, me, np, nxch, i, ip, iv
  character(len=120) :: fmt

  me = this_image()
  np = num_images()
  
  if (np==1) then
    !    allocate(xchg(0))
    xchg = [ integer :: ]
    
  else if (me == 1) then
    xchg = [me+1]
  else if (me == np) then
    xchg = [me-1]
  else
    xchg = [me-1, me+1]
  end if
  nxch = size(xchg)
  nrcv = nxch * nhl 

  allocate(expected(nxch, nhl))
  xv(1:nloc) = [(i,i=(me-1)*nloc+1,me*nloc)]
  iv = nloc + 1
  do ip=1, nxch
    loc_idx(1:nhl,ip) = [ (i,i=iv,iv+nhl-1) ]
    if (xchg(ip) == me-1) then
      rmt_idx(1:nhl,ip) = [ (i,i=nloc-nhl+1,nloc) ]
    else
      rmt_idx(1:nhl,ip) = [ (i,i=1,nhl) ]
    end if
    iv = iv + nhl
  end do

  sync images(xchg)
  iv = nloc + 1
  do ip=1, nxch
    !    xv(loc_idx(1:nhl,ip)) = xv(rmt_idx(1:nhl,ip))[xchg(ip)]
    write(*,*) me,': Reading from :',xchg(ip),' :',rmt_idx(1:nhl,ip),' :',loc_idx(1:nhl,ip)
          
    xv(iv:iv+nhl-1) = xv(rmt_idx(1:nhl,ip))[xchg(ip)]
    iv = iv + nhl
  end do
  if (me == 1) then
    expected(1,:) = nloc + rmt_idx(1:nhl, 1)
  else if (me == np) then
    expected(1,:) = ((np - 2) * nloc) + rmt_idx(1:nhl, 1)
  else
    expected(1,:) = ((me - 2) * nloc) + rmt_idx(1:nhl, 1)
    expected(2,:) = ( me      * nloc) + rmt_idx(1:nhl, 2)
  end if

 
  sync all 
  do ip=1, nxch
    write(fmt,*) '( i0,a,',nloc,'(f5.0,1x),a,',nrcv,'(f5.0,1x) )'
    write(*,fmt) me,': From halo exchange :',xv(1:nloc),':',xv(nloc+1:nloc+nrcv)
    if (any(xv(loc_idx(1:nhl, ip)) /= expected(ip,:))) error stop 'Test failed.'
  end do

  sync all
  if (me == 1) print *, 'Test passed.'
end program get_with_vector_index

