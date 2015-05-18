
!#include the fortran preprocessor definitions
#include "ppdef.h"

module grids

real, save :: toleranceX = 0.05, &   ! degrees
              toleranceY = 0.05, &   ! degrees
              toleranceZ = 5.        ! meters
 

integer, parameter ::   &
  reg             = 1, &
  parall          = 2, &
  irreg           = 3


! description for a 2D-grid

type grid2D
  integer typex,typey
  integer imax,jmax
  real x0,dx,y0,dy, valex
  real, dimension(:,:), allocatable :: x,y
  integer, dimension(:,:), allocatable :: mask
end type

! description for a 3D-grid

type grid3D
  integer typex, typey, typez
  integer imax,jmax,kmax
  real x0,dx,y0,dy,z0,dz, valex
  real, dimension(:,:,:), allocatable :: x,y,z
  integer, dimension(:,:,:),allocatable :: mask
end type


interface interp
  module procedure interp_1D, interp_2D, interp_3D
end interface

interface cinterp
  module procedure cinterp_1D, cinterp_2D, cinterp_3D
end interface

interface coord
  module procedure coord_3D, coord_2D
end interface

interface transCoord
  module procedure transCoord_3D
end interface

contains

!_______________________________________________________
!

subroutine cinterp_1D(z,mask,zi,sk,sc,n)
  implicit none
  real, intent(in)     :: z(:), zi
  integer, intent(in)  :: mask(:)
  integer, intent(out) :: sk(2)
  real, intent(out)    :: sc(2)
  integer, intent(out) :: n

  real :: zp(size(z)+1), rk, surf, depth
  integer :: k,kp,kmax

!  write(99,*) 'z ',z(:)
!  write(99,*) 'mask ',mask(:)
!  write(99,*) 'zi ',zi

  n = 0
  kmax = size(z)

  if (count(mask.eq.1).eq.0) then
! no sea points
    n = 0
    return
  end if

  k = -1

  do kp = 1,kmax-1
    if (z(kp).le.zi.and.zi.lt.z(kp+1)) then
      k = kp
      exit
    end if
  end do

!  write(99,*) 'linear k ',k

  if (k.ne.-1) then
    if (mask(k)+mask(k+1).eq.2) then
! linear
      n = 2
      rk = k+(zi-z(k))/(z(k+1)-z(k)) 
      sk(1) = k;   sc(1) = 1-rk+k;
      sk(2) = k+1; sc(2) = rk-k;
      return
    end if
  end if

!  k = minloc(abs(z-zi),mask.eq.1)

  k = transfer(minloc(abs(z-zi),mask.eq.1),k)

  if (abs(z(k)-zi).lt.toleranceZ) then
! nearest
      n = 1
      sk(1) = k;   sc(1) = 1;
  else
! nearest point masked
      n = 0
  end if

!!$ 
!!$  zp(kmax+1) = 0
!!$  
!!$  do k=kmax,1,-1
!!$      zp(k) = 2*z(k) - zp(k+1) 
!!$  end do
!!$
!!$  k = -1
!!$
!!$  do kp = 1,kmax
!!$        if (zp(kp).le.zi.and.zi.lt.zp(kp+1).and.mask(kp).eq.1) then
!!$          k = kp
!!$          exit
!!$		end if
!!$  end do
!!$
!!$  write(99,*) 'nearest kp ',kp
!!$
!!$  if (k.eq.-1) then
!!$    surf  = 0
!!$    depth = minval(zp,mask.eq.1)
!!$
!!$    if (surf.le.zi.and.zi.lt.surf+toleranceZ) then
!!$          k=kmax
!!$    elseif (depth-toleranceZ.le.zi.and.zi.lt.depth) then
!!$          k=1; 
!!$    else
!!$          n = 0
!!$  ! Out of the grid
!!$          return
!!$    end if
!!$  end if
!!$
!!$  if (mask(k).eq.1) then
!!$! nearest
!!$      n = 1
!!$      sk(1) = k;   sc(1) = 1;
!!$  else
!!$! nearest point masked
!!$      n = 0
!!$  end if
end subroutine

!_______________________________________________________
!

subroutine cinterp_2D(xy,xi,yi,si,sj,sc,n)
implicit none
type(grid2D), intent(in) :: xy
real,    intent(in)  :: xi, yi
integer, intent(out) :: si(4),sj(4), n
real,    intent(out) :: sc(4)

logical              :: out
real                 :: ri,rj,a,b,alphax,alphay
integer              :: i,j,ni,nj
integer :: istat
n = 0

if (xy%typex.eq.reg.and.xy%typey.eq.reg) then

      out = xi.lt.xy%x0+xy%dx-toleranceX.or.xi.gt.xy%x0+xy%dx*xy%imax+toleranceX.or. &
          yi.lt.xy%y0+xy%dy-toleranceY.or.yi.gt.xy%y0+xy%dy*xy%jmax+toleranceY

      if (out) return

      ri = (xi - xy%x0)/xy%dx;  
      rj = (yi - xy%y0)/xy%dy;  
      i = ri
      j = rj

      if (i.lt.1)       then; i=1;         ri = i;       end if
      if (i.ge.xy%imax) then; i=xy%imax-1; ri = xy%imax; end if
      if (j.lt.1)       then; j=1;         rj = j;       end if
      if (j.ge.xy%jmax) then; j=xy%jmax-1; rj = xy%jmax; end if

      if (xy%mask(i,j)   + xy%mask(i+1,j)    + &
          xy%mask(i,j+1) + xy%mask(i+1,j+1).eq.4) then
        n = 4
        a = ri-i
        b = rj-j

        si(1)=i;   sj(1)=j;   sc(1)=(1-b)*(1-a)
        si(2)=i+1; sj(2)=j;   sc(2)=(1-b)*a
        si(3)=i;   sj(3)=j+1; sc(3)=b*(1-a)
        si(4)=i+1; sj(4)=j+1; sc(4)=b*a
      else
        i = (xi - xy%x0)/xy%dx+.5
        j = (yi - xy%y0)/xy%dy+.5

        if (i.lt.1)       i=1;       
        if (i.ge.xy%imax) i=xy%imax; 
        if (j.lt.1)       j=1; 
        if (j.ge.xy%jmax) j=xy%jmax;

        if (xy%mask(i,j).eq.1) then 
          n = 1
          si(1) = i; sj(1) = j; sc(1) = 1
        else
          n = 0
		end if
      end if

elseif (xy%typex.eq.parall.and.xy%typey.eq.parall) then

      out = xi.lt.xy%x(1,1).or.xi.gt.xy%x(xy%imax,xy%jmax).or. &
          yi.lt.xy%y(1,1).or.yi.gt.xy%y(xy%imax,xy%jmax)
      if (out) return

      !find lower left corner
      do i=1,((xy%imax)-1)
        if (((xy%x(i,1)).le.xi).and.(xi.le.(xy%x(i+1,1)))) exit
      end do
      do j=1,((xy%jmax)-1)
        if (((xy%y(1,j)).le.yi).and.(yi.le.(xy%y(1,j+1)))) exit
      end do
      !coefficients in x and y directions
      alphax=(xi-(xy%x(i,1)))/((xy%x(i+1,1))-(xy%x(i,1)))
      alphay=(yi-(xy%y(1,j)))/((xy%y(1,j+1))-(xy%y(1,j)))
      !square interpolation indexes and coefficients
      si(1)=i   ; sj(1)=j   ; sc(1)=(1-alphax)*(1-alphay)
      si(2)=i+1 ; sj(2)=j   ; sc(2)= alphax   *(1-alphay)
      si(3)=i+1 ; sj(3)=j+1 ; sc(3)= alphax   * alphay
      si(4)=i   ; sj(4)=j+1 ; sc(4)=(1-alphax)* alphay
      n=4
else

  write(stderr,*) 'cinterp_2D: not yet implemented'
  call flush(stderr,istat)
  stop
end if

end subroutine


!_______________________________________________________
!

subroutine cinterp_3D(xy,xi,yi,zi,si,sj,sk,sc,n)
implicit none
type(grid3D), intent(in) :: xy
real, intent(in) :: xi, yi, zi
integer, intent(out) :: si(8),sj(8),sk(8),n
real,    intent(out) :: sc(8)
logical              :: out

real :: ri,rj,rk,a,b,c, surf,depth
integer :: i,j,k,kp, ni,nj,nk, n1,n2,n3,n4
integer :: istat

n = 0

if (xy%typex.eq.reg.and.xy%typey.eq.reg.and.xy%typez.eq.reg) then
    out = xi.lt.xy%x0+xy%dx-toleranceX.or.xi.gt.xy%x0+xy%dx*xy%imax+toleranceX.or. &
          yi.lt.xy%y0+xy%dy-toleranceY.or.yi.gt.xy%y0+xy%dy*xy%jmax+toleranceY.or. &
          zi.lt.xy%z0+xy%dz-toleranceZ.or.zi.gt.xy%z0+xy%dz*xy%kmax+toleranceZ

    if (out) return

    ri = (xi - xy%x0)/xy%dx;  
    rj = (yi - xy%y0)/xy%dy;  
    rk = (zi - xy%z0)/xy%dz;  

    i = ri
    j = rj
    k = rk

    if (i.lt.1)       then; i=1;         ri = i;       end if
    if (i.ge.xy%imax) then; i=xy%imax-1; ri = xy%imax; end if
    if (j.lt.1)       then; j=1;         rj = j;       end if
    if (j.ge.xy%jmax) then; j=xy%jmax-1; rj = xy%jmax; end if
    if (k.lt.1)       then; k=1;         rk = k;       end if
    if (k.ge.xy%kmax) then; k=xy%kmax-1; rk = xy%kmax; end if

    if (xy%mask(i,j,  k)   + xy%mask(i+1,j,  k)    + &
        xy%mask(i,j+1,k)   + xy%mask(i+1,j+1,k)    + &
        xy%mask(i,j,  k+1) + xy%mask(i+1,j,  k+1)  + &
        xy%mask(i,j+1,k+1) + xy%mask(i+1,j+1,k+1).eq.8) then
! trilinear
        n= 8
        a = ri-i
        b = rj-j
        c = rk-k

        si(1)=i;   sj(1)=j;   sk(1)=k;   sc(1)=(1-a)*(1-b)*(1-c)
        si(2)=i+1; sj(2)=j;   sk(2)=k;   sc(2)=   a *(1-b)*(1-c)
        si(3)=i;   sj(3)=j+1; sk(3)=k;   sc(3)=(1-a)*   b *(1-c)
        si(4)=i+1; sj(4)=j+1; sk(4)=k;   sc(4)=   a *   b *(1-c)
        si(5)=i;   sj(5)=j;   sk(5)=k+1; sc(5)=(1-a)*(1-b)*c
        si(6)=i+1; sj(6)=j;   sk(6)=k+1; sc(6)=   a *(1-b)*c
        si(7)=i;   sj(7)=j+1; sk(7)=k+1; sc(7)=(1-a)*   b *c
        si(8)=i+1; sj(8)=j+1; sk(8)=k+1; sc(8)=   a *   b *c
    else
        i = (xi - xy%x0)/xy%dx+.5
        j = (yi - xy%y0)/xy%dy+.5
        k = (zi - xy%z0)/xy%dz+.5

        if (i.lt.1)       i=1;       
        if (i.gt.xy%imax) i=xy%imax; 
        if (j.lt.1)       j=1; 
        if (j.gt.xy%jmax) j=xy%jmax;
        if (k.lt.1)       k=1; 
        if (k.gt.xy%kmax) k=xy%kmax;

        if (xy%mask(i,j,k).eq.1) then 
	  ! nearest
          n = 1
          si(1) = i; sj(1) = j; sk(1) = k; sc(1) = 1
        else
          n = 0
		end if
	end if
elseif (xy%typex.eq.reg.and.xy%typey.eq.reg.and.xy%typez.eq.irreg) then
      out = xi.lt.xy%x0+xy%dx-toleranceX.or.xi.gt.xy%x0+xy%dx*xy%imax+toleranceX.or. &
            yi.lt.xy%y0+xy%dy-toleranceY.or.yi.gt.xy%y0+xy%dy*xy%jmax+toleranceY

      if (out) return

      ri = (xi - xy%x0)/xy%dx;  
      rj = (yi - xy%y0)/xy%dy;  
      i = ri
      j = rj

      if (i.lt.1)       then; i=1;         ri = i;       end if
      if (i.ge.xy%imax) then; i=xy%imax-1; ri = xy%imax; end if
      if (j.lt.1)       then; j=1;         rj = j;       end if
      if (j.ge.xy%jmax) then; j=xy%jmax-1; rj = xy%jmax; end if

      call cinterp_1D(xy%z(i,  j,  :),xy%mask(i,  j,  :),zi,sk(1:2),sc(1:2),n1)
      call cinterp_1D(xy%z(i,  j+1,:),xy%mask(i,  j+1,:),zi,sk(3:4),sc(3:4),n2)
      call cinterp_1D(xy%z(i+1,j,  :),xy%mask(i+1,j,  :),zi,sk(5:6),sc(5:6),n3)
      call cinterp_1D(xy%z(i+1,j+1,:),xy%mask(i+1,j+1,:),zi,sk(7:8),sc(7:8),n4)

!      write(stderr,*) 'n1,n2,n3,n4 ',n1,n2,n3,n4

      if (n1+n2+n3+n4.eq.8) then
      ! trilinear
        n= 8
        a = ri-i
        b = rj-j

        si(1)=i;   sj(1)=j;   sc(1)=(1-a)*(1-b)*sc(1)
        si(2)=i;   sj(2)=j;   sc(2)=(1-a)*(1-b)*sc(2)
        si(3)=i;   sj(3)=j+1; sc(3)=(1-a)*   b *sc(3)
        si(4)=i;   sj(4)=j+1; sc(4)=(1-a)*   b *sc(4)
        si(5)=i+1; sj(5)=j;   sc(5)=   a *(1-b)*sc(5)
        si(6)=i+1; sj(6)=j;   sc(6)=   a *(1-b)*sc(6)
        si(7)=i+1; sj(7)=j+1; sc(7)=   a *   b *sc(7)
        si(8)=i+1; sj(8)=j+1; sc(8)=   a *   b *sc(8)
 
      elseif (n1.eq.1.and.n2.eq.1.and.n3.eq.1.and.n4.eq.1) then
      ! horiziontal: bilinear
      ! vertical: nearest

        n= 4
        a = ri-i
        b = rj-j

        si(1)=i;   sj(1)=j;   sk(1)=sk(1);  sc(1)=(1-a)*(1-b)*sc(1)
        si(2)=i;   sj(2)=j+1; sk(2)=sk(3);  sc(2)=(1-a)*   b *sc(3)
        si(3)=i+1; sj(3)=j;   sk(3)=sk(5);  sc(3)=   a *(1-b)*sc(5)
        si(4)=i+1; sj(4)=j+1; sk(4)=sk(7);  sc(4)=   a *   b *sc(7)

      else
       ! horizontal: nearest 
       ! vertical: nearest or linear

        i = (xi - xy%x0)/xy%dx+.5
        j = (yi - xy%y0)/xy%dy+.5

        if (i.lt.1)       i=1;       
        if (i.ge.xy%imax) i=xy%imax; 
        if (j.lt.1)       j=1; 
        if (j.ge.xy%jmax) j=xy%jmax;

        si(1)=i;   sj(1)=j;   
        si(2)=i;   sj(2)=j;   

        call cinterp_1D(xy%z(i,  j,  :),xy%mask(i,  j,  :),zi,sk(1:2),sc(1:2),n)
      end if
else
  write(stderr,*) 'cinterp_3D: not jet implemented'
  call flush(stderr,istat)
  stop
end if

end subroutine


!_______________________________________________________
!

subroutine interp_1D(z,mask,v,zi,vi,out)
implicit none

  real, intent(in)     :: z(:)
  integer, intent(in)  :: mask(size(z))
  real, intent(in) :: v(size(z)), zi
  real, intent(out) :: vi
  logical, intent(out) :: out

  real :: sc(2)
  integer :: i,sk(2), n

  call cinterp(z,mask,zi,sk,sc,n)
  out = n.eq.0
  vi = 0
  do i=1,n
    vi = vi+v(sk(i)) * sc(i)
  end do
end subroutine

!_______________________________________________________
!

subroutine interp_2D(xy,v,xi,yi,vi,out)
implicit none
type(grid2D), intent(in) :: xy
real, intent(in) :: v(xy%imax,xy%jmax), xi, yi
real, intent(out) :: vi
logical, intent(out) :: out

  real :: sc(4)
  integer :: i,si(4),sj(4), n
  call cinterp(xy,xi,yi,si,sj,sc,n)
  out = n.eq.0
  vi = 0
  do i=1,n
    vi = vi+v(si(i),sj(i)) * sc(i)
  end do

end subroutine


!_______________________________________________________
!


subroutine interp_3D(xy,v,xi,yi,zi,vi,out)
implicit none
type(grid3D), intent(in) :: xy
real, intent(in) :: v(xy%imax,xy%jmax,xy%kmax), xi, yi, zi
real, intent(out) :: vi
logical, intent(out) :: out

  real :: sc(8)
  integer :: i,si(8),sj(8),sk(8), n

  call cinterp(xy,xi,yi,zi,si,sj,sk,sc,n)
  out = n.eq.0
  vi = 0
  do i=1,n
    vi = vi+v(si(i),sj(i),sk(i)) * sc(i)
  end do
end subroutine

!_______________________________________________________
!


subroutine transCoord_3D(grid1,grid2,x1,y1,z1,x2,y2,z2,out)
implicit none
type(grid3D), intent(in) :: grid1, grid2
real, intent(in)  :: x1,y1,z1
real, intent(out) :: x2,y2,z2
logical, intent(out) :: out

integer :: i,si(8),sj(8),sk(8),n
real    :: sc(8),x,y,z


call cinterp(grid1,x1,y1,z1,si,sj,sk,sc,n)
out = n.eq.0

if (out) return

x2 = 0
y2 = 0
z2 = 0

do i=1,n
    call coord(grid2,si(i),sj(i),sk(i),x,y,z,out)
    x2 = x2 + x*sc(i)
    y2 = y2 + y*sc(i)
    z2 = z2 + z*sc(i)
end do

end subroutine



!_______________________________________________________
!

subroutine coord_2D(grid,i,j,x,y,out)
implicit none
type(grid2D), intent(in) :: grid
integer, intent(in) :: i,j
real, intent(out) :: x,y
logical, intent(out) :: out
integer :: istat

out = i.lt.1.or.i.gt.grid%imax.or. &
      j.lt.1.or.j.gt.grid%jmax

if (out) return

if (grid%typex.eq.reg) then
  x = grid%x0+i*grid%dx
else if (grid%typex.eq.parall.or.grid%typex.eq.irreg) then
  x = grid%x(i,j)
else
  write(stderr,*) 'coord: not jet implemented grid%typex',grid%typex
  call flush(stderr,istat)
  stop
end if

if (grid%typey.eq.reg) then
  y = grid%y0+j*grid%dy
else if (grid%typey.eq.parall.or.grid%typey.eq.irreg) then
  y = grid%y(i,j)
else
  write(stderr,*) 'coord: not jet implemented grid%typey',grid%typey
  call flush(stderr,istat)
  stop
end if

end subroutine

!_______________________________________________________
!


subroutine coord_3D(grid,i,j,k,x,y,z,out)
implicit none
type(grid3D), intent(in) :: grid
integer, intent(in) :: i,j,k
real, intent(out) :: x,y,z
logical, intent(out) :: out

integer :: istat
out = i.lt.1.or.i.gt.grid%imax.or. &
      j.lt.1.or.j.gt.grid%jmax.or. &
      k.lt.1.or.k.gt.grid%kmax

if (out) return

if (grid%typex.eq.reg) then
  x = grid%x0+i*grid%dx
else if (grid%typex.eq.irreg) then
  x = grid%x(i,j,k)
else
  write(stderr,*) 'coord: not jet implemented grid%typex',grid%typex
  call flush(stderr,istat)
  stop
end if

if (grid%typey.eq.reg) then
  y = grid%y0+j*grid%dy
else if (grid%typey.eq.irreg) then
  y = grid%y(i,j,k)
else
  write(stderr,*) 'coord: not jet implemented grid%typey',grid%typey
  call flush(stderr,istat)
  stop
end if

if (grid%typez.eq.reg) then
  z = grid%z0+k*grid%dz
else if (grid%typez.eq.irreg) then
  z = grid%z(i,j,k)
else
  write(stderr,*) 'coord: not jet implemented grid%typez',grid%typez
  call flush(stderr,istat)
  stop
end if


end subroutine

!_______________________________________________________
!

subroutine interp_field2D(grid1,field1,grid2,field2,outgrid)
implicit none
type(grid2D), intent(in) :: grid1, grid2
real, intent(in) :: field1(:,:)
logical, intent(out), optional :: outgrid(grid2%imax,grid2%jmax)
real, intent(out) :: field2(grid2%imax,grid2%jmax)

logical :: out
integer :: i2,j2
real :: x2,y2

do j2=1,grid2%jmax
do i2=1,grid2%imax
if (grid2%mask(i2,j2).eq.1) then
  call coord(grid2,i2,j2,x2,y2,out)
  if (.not.out) then
    call interp(grid1,field1,x2,y2,field2(i2,j2),out)
  end if

  if (present(outgrid)) outgrid(i2,j2) = out
else
  field2(i2,j2)=0
end if
end do
end do

end subroutine

!_______________________________________________________
!

function interp_field3D(grid1,field1,grid2,outgrid) result(field2)
implicit none
type(grid3D), intent(in) :: grid1, grid2
real, intent(in) :: field1(:,:,:)
logical, intent(out), optional :: outgrid(grid2%imax,grid2%jmax,grid2%kmax)
real :: field2(grid2%imax,grid2%jmax,grid2%kmax)

logical :: out
integer :: i2,j2,k2
real :: x2,y2,z2

do k2=1,grid2%kmax
do j2=1,grid2%jmax
do i2=1,grid2%imax
if (grid2%mask(i2,j2,k2).eq.1) then
  call coord(grid2,i2,j2,k2,x2,y2,z2,out)
  if (.not.out) then
    call interp(grid1,field1,x2,y2,z2,field2(i2,j2,k2),out)
  end if

  if (present(outgrid)) outgrid(i2,j2,k2) = out
else
  field2(i2,j2,k2)=0
end if
end do
end do
end do

end function


end module

