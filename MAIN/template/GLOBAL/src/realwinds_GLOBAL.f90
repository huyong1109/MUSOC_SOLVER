MODULE REALWINDS_GLOBAL
! reading value from init subroutine then private use in subroutine
  TYPE :: REALWINDS_TYPE
    character(len=5) :: WINDNAME,WINDTYPE
    character(len=4) :: YEARS,HOURS
    character(len=40) :: DIRNAME
    integer :: nxdeg,nydeg
    real,allocatable :: xdeg(:),ydeg(:)
  END TYPE REALWINDS_TYPE

  INTEGER,PARAMETER,PRIVATE :: MXWIND=12
  TYPE(REALWINDS_TYPE) :: WIND(MXWIND)

  INTEGER,PRIVATE :: WINDID(MXWIND)=0
  PRIVATE :: ncin_windfield,nc_check

! domain size for realwind data
  INTEGER,PARAMETER,PRIVATE :: No=144,Na=73
  REAL,PRIVATE :: temp(No,Na),lon(No),lat(Na)

! permanent format for data access
  CHARACTER*5,PARAMETER,PRIVATE :: tar(2)=(/"E-WSS","N-SSS"/)
  CHARACTER*4,PARAMETER,PRIVATE :: tat(2)=(/"ewss","nsss"/)

CONTAINS

! ======================================================================

SUBROUTINE REALWIND_INITIAL(wn,wt,yr,hr,dn,nxdeg,xdeg,nydeg,ydeg,id)
  character(len=*) :: wn,wt,yr,hr,dn
  real :: x0 
  integer :: nxdeg,nydeg
  real :: xdeg(nxdeg),ydeg(nydeg)
  integer,intent(out) :: id

  do id = 1,MXWIND
     if ( WINDID(id) == 0 ) then
       WINDID=id
       exit
     endif
  enddo
  allocate(WIND(id)%xdeg(nxdeg))
  allocate(WIND(id)%ydeg(nydeg))
  WIND(id)%WINDNAME = trim(wn)
  WIND(id)%WINDTYPE = trim(wt)
  WIND(id)%YEARS = trim(yr)
  WIND(id)%HOURS = trim(hr)
  WIND(id)%DIRNAME = trim(dn)
  WIND(id)%nydeg = nxdeg
  WIND(id)%nydeg = nydeg
  WIND(id)%xdeg = xdeg
  WIND(id)%ydeg = ydeg
END SUBROUTINE REALWIND_INITIAL

SUBROUTINE REALWIND_RELEASE(id)
  WINDID(id)=0
  WIND(id)%WINDNAME = trim("")
  WIND(id)%WINDTYPE = trim("")
  WIND(id)%YEARS = trim("")
  WIND(id)%HOURS = trim("")
  WIND(id)%DIRNAME = trim("")
  WIND(id)%nxdeg = 0
  WIND(id)%nydeg = 0
  deallocate(WIND(id)%xdeg)
  deallocate(WIND(id)%ydeg)
END SUBROUTINE REALWIND_RELEASE 

! ======================================================================

SUBROUTINE REALWIND_INFLOW(id,ift1,I0,J0,ipl,I1,J1)
!  id = id for specified wind 
!  ift1= time you select
!  ipl = output from interpolation
! variables declared for argument 
  integer :: id
  integer :: ivy,ift1
  character(len=5) :: var
  integer :: I0,J0
  real :: ipl(I0,J0)
  integer :: I1,J1

  !in subroutine variables declaration
  integer :: ip
  character(len=60) :: fn
  real :: XLON,YLAT

  do i=1,2
     if (tar(i) .eq. WIND(id)%WINDTYPE) ip=i
  enddo

  write(fn,15)trim(WIND(id)%DIRNAME),  &
              trim(WIND(id)%WINDNAME), &
              trim(WIND(id)%WINDTYPE), &
              trim(WIND(id)%YEARS),    &
              trim(WIND(id)%HOURS)
  15 format(A,'/',A,'-',A5,'-',I4,'-',A,'.nc')

  write(*,*)trim(fn)
  call ncin_windfield(fn,tat(ip),ift1,No,Na,lon,lat,temp)
  DO JJ = 2, J1
  DO II = 2, I1
     !XLON = WIND(id)%x0deg +(II-1.5)*XINC
     XLON = WIND(id)%xdeg(II)
     YLAT = WIND(id)%ydeg(JJ)
     IF (XLON.GT.360.) XLON=XLON-360.
     DO i=1,NO
        if (abs(lon(i)-XLON) .le. 2.5) ien=i
     ENDDO
     DO i=1,NA
         if (abs(lat(i)-YLAT) .le. 2.5) jen=i
     ENDDO 
     iwn=ien-1
     jwn=jen-1
     XNE=abs(lon(iwn)-XLON)/2.5
     XNW=1-XNE
     XSS=abs(lat(jwn)-YLAT)/2.5
     XNN=1-XSS
     DNW=temp(iwn,jwn)
     DNE=temp(ien,jwn)
     DSW=temp(iwn,jen)
     DSE=temp(ien,jen)
     DN=DNW*XNW+DNE*XNE
     DS=DSW*XNW+DSE*XNE
     ipl(II,JJ)=DN*XNN+DS*XSS
!convert from (Newtons/m-m)*21600 to dynes/cm-cm
     ipl(II,JJ)=ipl(II,JJ)/2160.
! 21600 is the number of seconds in 6 hours and relates to the accumulation
! of data snapshots
  ENDDO
  ENDDO
END SUBROUTINE REALWIND_INFLOW 

! ======================================================================

subroutine ncin_windfield(fn,tat,tfr,n1,n2,x1,x2,field)
use netcdf
  character(len=*) :: fn
  character(len=4) :: tat
  integer :: tfr,n1,n2
  real :: x1(n1),x2(n2),field(n1,n2)
  integer :: ncid,varid
  real :: scl,ado
  ! how data to be read
  integer :: istart(3),icount(3),istrid(3)

  istart=(/1,1,1/)
  icount=(/n1,n2,1/)
  istrid=(/1,1,1/)
  istart(3)=tfr
  istrid(3)=tfr+1

  ! open the nc file, and get the fileid=ncid
  call nc_check( nf90_open(fn , nf_nowrite , ncid) )
  ! get the tat data from the fileid=ncid
  call nc_check( nf90_inq_varid(ncid, tat, varid))
  call nc_check( nf90_get_att(ncid, varid,'scale_factor',scl) )
  call nc_check( nf90_get_att(ncid, varid,'add_offset',ado) )
  call nc_check( nf90_get_var(ncid, varid, field, &
                 start=istart, count=icount, stride=istrid) )
  field = field*scl + ado
  ! get the longitude data from the fileid=ncid
  call nc_check( nf90_inq_varid(ncid, 'longitude', varid) )
  call nc_check( nf90_get_var(ncid, varid, x1) )
  ! get the latitude data from the fileid=ncid
  call nc_check( nf90_inq_varid(ncid, 'latitude', varid) )
  call nc_check( nf90_get_var(ncid, varid, x2) )
  ! close the fileid=ncid
  call nc_check( nf90_close(ncid) )

end subroutine ncin_windfield

subroutine nc_check(status)
use netcdf
integer, intent ( in) :: status

   if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "netcdf error" 
    end if
end subroutine nc_check


END MODULE REALWINDS_GLOBAL
