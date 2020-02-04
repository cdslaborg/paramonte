program whole_array_get
  implicit none

  integer,allocatable :: x1(:)[:],y1(:)
  integer,allocatable :: x2(:,:)[:],y2(:,:)
  integer,allocatable :: x3(:,:,:)[:],y3(:,:,:)
  integer,parameter   :: n = 10
  integer :: me,np,i,j,k

  me = this_image()
  np = num_images()

  allocate(x1(n)[*],y1(n))

  x1 = me
  y1 = 0

  sync all

  if(me == 1) then
     y1 = x1(:)[me+1]
     if(any(y1 /= 2)) then
        error stop 'Test 1 fails'
     end if
  end if

  deallocate(x1)
  allocate(x2(1:n,0:n-1)[*],y2(n,n))

  x2 = me
  y1 = 0; y2 = 0

  sync all
  if(me == 1) then
     y2 = x2(:,:)[np]
     if(any(y2 /= np)) then
        error stop 'Test 2 fails'
     end if
  end if

  sync all

  x2(:,n/2) = x2(:,n/2) + n/2

  sync all

  if(me == 1) then
     y1 = x2(:,n/2)[me+1]
     if(any(y1 /= 2+n/2)) then
        error stop 'Test 3 fails'
     end if
  end if

  deallocate(y1,x2,y2)
  allocate(x3(0:n-1,1:n,-1:n-2)[*],y3(n,n,n))

  x3 = me; y3 = 0

  sync all

  if(me == 1) then
     y3 = x3(:,:,:)[me+1]
     if(any(y3 /= me+1)) then
        error stop 'Test 4 fails'
     end if
  endif

  sync all

  x3(:,:,n/2) = me + n/2
  y3 = 0

  sync all

  if(me == 1) then
     y3(:,n/2,:) = x3(:,:,n/2)[me+1]
     if(any(y3(:,n/2,:) /= me+1+n/2)) then
        error stop 'Test 5 fails'
     end if
  endif

  if(me == 1) write(*,*) 'Test passed.'

end program whole_array_get
