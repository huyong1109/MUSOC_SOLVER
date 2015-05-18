module interfext

 use grids

!luc 2012/06/08 (ncar)
 use ufileformat

 interface init
   module procedure init_2D
   module procedure init_2D_nomodel
 end interface

 interface getfield
   module procedure getfield_2D
 end interface

 type field2D
   character(len=256) :: initfname, key

   ! description of the model grid

   type(grid2D) :: model

   ! description of the data grid

   type(grid2D) :: g2D
   logical :: initialised
   real, allocatable :: array(:,:,:)
   real, allocatable :: day(:)
   real, allocatable :: field1(:,:),field2(:,:)

   ! scale factor
   real :: scale
   integer :: index

 end type field2D

 ! fortran unit for output

 integer :: stdout = 6

contains

 subroutine init_2D(f,x0,dx,y0,dy,imax,jmax,initfname,key,mask,scale)
  use ufileformat
  use initfile
  implicit none
  type(field2D), intent(out) :: f
  real, intent(in) :: x0,dx,y0,dy
  integer, intent(in) :: imax,jmax
  character(len=*), intent(in) :: initfname,key
  real, optional, intent(in) :: scale
  integer(2), optional, intent(in) :: mask(:,:)

  character(len=256) :: str,path
  real, allocatable :: d(:)
  real :: valex, dummy
  integer :: dummyi,dummyb
  integer :: i,j,k

  f%model%typex = reg
  f%model%x0=x0
  f%model%dx=dx
  f%model%imax=imax

  f%model%typey = reg
  f%model%y0=y0
  f%model%dy=dy
  f%model%jmax=jmax
 
  f%initfname = initfname
  f%key = key

  if (present(scale)) then
    f%scale = scale
  else
    f%scale = 1.
  end if
!  write(99,*) 'scales',f%scale

  allocate(f%model%mask(imax,jmax))
  if (present(mask)) then
    f%model%mask = mask
  else
    f%model%mask = 1
  end if

!  write(99,*) 'jmbtest',initfname
  call getInitValue(initfname,trim(key)//'path',path); 

  f%g2D%typex = reg
  call getInitValue(initfname,trim(key)//'gridX',str); 
  call uload(trim(path)//str,f%g2D%x0,f%g2D%dx,dummy,dummy,f%g2D%imax,f%g2D%jmax,dummyi)

  f%g2D%typey = reg
  call getInitValue(initfname,trim(key)//'gridY',str); 
  call uload(trim(path)//str,f%g2D%y0,dummy,f%g2D%dy,dummy,f%g2D%imax,f%g2D%jmax,dummyi)

  allocate(f%g2D%mask(f%g2D%imax,f%g2D%jmax),f%field1(imax,jmax),f%field2(imax,jmax))
  f%g2D%mask = 1

  call getInitValue(initfname,trim(key)//'mjd',str); 
  call uload(trim(path)//str,d,valex)
  allocate(f%day(size(d)))
  f%day = d
  deallocate(d)

!  write(99,*) 'jmbtestvalue',initfname
  call getInitValue(initfname,trim(key)//'value',str); 
  call uload(trim(path)//str,f%array,valex)

!  write(99,*) 'jmbtestvaluesuccess',initfname
  dummyi = count(f%array.eq.valex)
  dummyb = count(f%array.ne.valex)
!  write(99,*) 'jmbtestvaluesuccesso whats',dummyi,dummyb


  if (dummyi.ne.0) then
    write(stdout,*) 'init_2D: WARNING ',trim(f%key),' contains exclusion points '

    where (f%array(:,:,1).eq.valex)
      f%g2D%mask = 0
    elsewhere
      f%g2D%mask = 1
    end where
  end if
!  write(99,*) 'jmbtestvaluesuccesso whatsva',valex
! DOES NOT WORK???
!JMB  f%array = f%scale * f%array
! ?????
  if (present(scale)) then
    do k=1,size(f%array,3)
    do j=1,size(f%array,2)
    do i=1,size(f%array,1)
    f%array(i,j,k) = scale*f%array(i,j,k)
    end do
    end do
    end do
  end if
!  write(99,*) 'jmbtestvaluesuccesso whats scaled'
  f%initialised = .false.
  f%index = 1
!  write(99,*) 'jmbtestvaluesuccesso whatsexit'
 end subroutine init_2D








 subroutine init_2D_nomodel(f,initfname,key,mask,scale)
  use ufileformat
  use initfile
  implicit none
  type(field2D), intent(out) :: f
  character(len=*), intent(in) :: initfname,key
  real, optional, intent(in) :: scale
  integer(2), optional, intent(in) :: mask(:,:)

  character(len=256) :: str,path
  real, allocatable :: temp1(:),temp2(:),d(:)
  real :: valex, dummy
  integer :: dummyi,dummyb
  integer :: i,j,k
  real :: rotate_lon

  f%initfname = initfname
  f%key = key

  if (present(scale)) then
    f%scale = scale
  else
    f%scale = 1.
  end if

  call getInitValue(initfname,trim(key)//'path',path);

  !COORDINATES
  !data is on a regular degenerated grid
  f%g2D%typex = reg
  f%g2D%typey = reg
  !read non-degenerated X,Y files
     call getInitValue(initfname,trim(key)//'gridX',str);
     call uload(trim(path)//str,temp1,f%g2D%valex) !ncep -> lon,lat = single array
     if (presentInitValue(initfname,trim(key)//'rotate_lon')) then
        call getInitValue(initfname,trim(key)//'rotate_lon',rotate_lon)
        temp1=temp1+rotate_lon
     endif
     call getInitValue(initfname,trim(key)//'gridY',str);
     call uload(trim(path)//str,temp2,f%g2D%valex) !roms -> lon,lat = matrix
  !compute degenerated values
  f%g2D%imax=size(temp1,1)
  f%g2D%jmax=size(temp2,1)
  f%g2D%dx=temp1(2)-temp1(1)
  f%g2D%x0=temp1(1)-f%g2D%dx
  f%g2D%dy=temp2(2)-temp2(1)
  f%g2D%y0=temp2(1)-f%g2D%dy
  deallocate(temp1,temp2)

  !MASK
  allocate(f%g2D%mask(f%g2D%imax,f%g2D%jmax))
                  !                         ,f%field1(imax,jmax),f%field2(imax,jmax))
  f%g2D%mask = 1

  !TIME
  call getInitValue(initfname,trim(key)//'mjd',str);
  call uload(trim(path)//str,d,valex)
  allocate(f%day(size(d)))
  f%day = d
  deallocate(d)

  !VALUES
  call getInitValue(initfname,trim(key)//'value',str);
  call uload(trim(path)//str,f%array,valex)

  dummyi = count(f%array.eq.valex)
  dummyb = count(f%array.ne.valex)

  if (dummyi.ne.0) then
    write(stdout,*) 'init_2D: WARNING ',trim(f%key),' contains exclusion points '

    where (f%array(:,:,1).eq.valex)
      f%g2D%mask = 0
    elsewhere
      f%g2D%mask = 1
    end where
  end if

  if (present(scale)) then
    do k=1,size(f%array,3)
    do j=1,size(f%array,2)
    do i=1,size(f%array,1)
    f%array(i,j,k) = scale*f%array(i,j,k)
    end do
    end do
    end do
  end if

  f%initialised = .false.
  f%index = 1

 end subroutine init_2D_nomodel










 subroutine getfield_2D(f,mjd,field,interp_in_space)
  use grids
  use ufileformat
  implicit none

  type(field2D), intent(inout) :: f
  real, intent(in) :: mjd
  real,intent(out) :: field(f%model%imax,f%model%jmax)
  logical, intent(in), optional :: interp_in_space

  logical :: out,space
  integer :: i,j,k
  real :: alpha

  if (present(interp_in_space)) then
    if (interp_in_space) then
       space=.true.
    else
       space=.false.
    end if
  else
    space=.true.
  end if


  if (f%day(f%index).le.mjd.and.mjd.le.f%day(f%index+1).and.f%initialised) then
    ! all data are already interpolated in space

  elseif (f%day(f%index+1).le.mjd.and.mjd.le.f%day(f%index+2).and.f%initialised) then
!!$omp single
    ! cronological reading, next slice

    f%index = f%index+1
    write(stdout,*) 'getfield_2D: Interpolate fields ',trim(f%key),' with index ',f%index,1+f%index

!    toleranceX = huge(toleranceX)
!    toleranceY = huge(toleranceY)
!    call interp_field2D(f%g2D,f%array(:,:,f%index+1),f%model,f%field2)

    if (space) then
       f%field1=f%field2
       call regreg(f%array(1,1,f%index+1),f%g2D%x0,f%g2D%dx,0.,f%g2D%y0,0.,f%g2D%dy,f%g2D%imax,f%g2D%jmax, &
                f%field2,f%model%x0,f%model%dx,0.,f%model%y0,0.,f%model%dy,f%model%imax,f%model%jmax)
    end if
    !if we don't interpolate in space, just remember index and (index+1)
    !write(stdout,*) 'getfield_2D: Interpolate fields ',trim(f%key),' with index ',f%index,f%index+1


!!$omp end single
  else  ! unordered reading or not f%initialised
!!$omp single
    if (mjd.lt.f%day(1)) then
      i = 1
      write(stdout,*) 'getfield: WARNING day ',mjd, &
           ' not covered by data ranging from ',f%day(1),' to ',f%day(size(f%day)), ' for field ', trim(f%key)
    elseif (mjd.gt.f%day(size(f%day))) then
      i = size(f%day)-1
      write(stdout,*) 'getfield: WARNING day ',mjd, &
           ' not covered by data ranging from ',f%day(1),' to ',f%day(size(f%day)), ' for field ', trim(f%key)
    else
      do i=1,size(f%day)-1
        if (f%day(i).le.mjd.and.mjd.le.f%day(i+1)) exit     
      end do
    end if

    if (.not.f%initialised.or.i.ne.f%index) then
      if (space) then
        toleranceX = huge(toleranceX)
        toleranceY = huge(toleranceY)
        call interp_field2D(f%g2D,f%array(:,:,i),f%model,f%field1)
        call interp_field2D(f%g2D,f%array(:,:,i+1),f%model,f%field2)
      end if
      write(stdout,*) 'getfield_2D: Interpolate fields ',trim(f%key),' with index ',i,i+1 
    end if

    f%index = i
    f%initialised = .true.
!!$omp end single
  end if

  alpha = (mjd-f%day(f%index+1))/(f%day(f%index)-f%day(f%index+1))

  ! interpolate to the nearest if mjd outside of data range

  if (alpha.lt.0) alpha = 0
  if (alpha.gt.1) alpha = 1

!!$omp do private(i)
  if (space) then
    do j=1,size(field,2)
      do i=1,size(field,1)
        field(i,j) = alpha * f%field1(i,j) + (1-alpha) * f%field2(i,j)
      end do
    end do
  else
    do j=1,size(field,2)
      do i=1,size(field,1)
        field(i,j) = alpha * f%array(i,j,f%index) + (1-alpha) * f%array(i,j,f%index+1)
      end do
    end do
  end if
!!$omp end do

  !write(stdout,*) 'getfield_2D: alpha ',alpha

 end subroutine getfield_2D

end module interfext





