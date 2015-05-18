MODULE NCIO
USE NETCDF
USE TIMCOM_GENERAL, ONLY: VAR_I0J0K1_HDF5
IMPLICIT NONE
INTEGER,PARAMETER,PRIVATE :: INT2 = SELECTED_INT_KIND(4)
INTEGER,PARAMETER,PRIVATE :: INT2BITS = 16
CHARACTER(LEN=*),PARAMETER,PRIVATE :: LON_NAME="Longitude", &
                                      LAT_NAME="Latitude",   &
                                      DEP_NAME="Depth",      &
                                      DAY_NAME="Time"
CHARACTER(LEN=*),PARAMETER,PRIVATE :: LON_UNITS="degrees_east", &
                                      LAT_UNITS="degrees_north",  &
                                      DEP_UNITS="cm",       &
                                      DAY_UNITS="days"
CHARACTER(LEN=*),PARAMETER,PRIVATE :: UNITS="units",             &
                                      SCA_NAME="scale_factor", &
                                      ADO_NAME="add_offset", &
                                      MIS_NAME="missing_value"

PRIVATE :: nc_check

CONTAINS

SUBROUTINE NC_FILECREATE(fn,NCID)
character(len=*) :: fn
integer :: NCID
  if (VAR_I0J0K1_HDF5.eq..false.) CALL nc_check( nf90_create(fn, NF90_CLOBBER, NCID) )
  if (VAR_I0J0K1_HDF5.eq..true.)  CALL nc_check( nf90_create(fn, NF90_HDF5, NCID) )
  call nc_check( nf90_put_att(ncid,NF90_GLOBAL,'Title','DieCAST/TIMCOM model') )
END SUBROUTINE NC_FILECREATE 

! ================================== !

SUBROUTINE NC_4DST_DEF_DIM(NCID,NLON,NLAT,NDEP,DIMIDS,VARIDS)
INTEGER :: NCID
INTEGER :: NLON,NLAT,NDEP
INTEGER :: DIMIDS(4),VARIDS(4)
!Temporory variables
INTEGER :: LON_DIMID,LAT_DIMID,DEP_DIMID,DAY_DIMID,LON_VARID,LAT_VARID,DEP_VARID,DAY_VARID

!Define the dimensions
   CALL nc_check( NF90_DEF_DIM(NCID,LON_NAME,NLON,LON_DIMID) )
   CALL nc_check( NF90_DEF_DIM(NCID,LAT_NAME,NLAT,LAT_DIMID) )
   CALL nc_check( NF90_DEF_DIM(NCID,DEP_NAME,NDEP,DEP_DIMID) )
   CALL nc_check( NF90_DEF_DIM(NCID,DAY_NAME,1,DAY_DIMID) )

!Define the coordinate variables
   if (VAR_I0J0K1_HDF5.eq..true.) then
     CALL nc_check( NF90_DEF_VAR(NCID,LON_NAME,NF90_REAL,LON_DIMID,LON_VARID,shuffle=.true.,deflate_level=9) )
     CALL nc_check( NF90_DEF_VAR(NCID,LAT_NAME,NF90_REAL,LAT_DIMID,LAT_VARID,shuffle=.true.,deflate_level=9) )
     CALL nc_check( NF90_DEF_VAR(NCID,DEP_NAME,NF90_REAL,DEP_DIMID,DEP_VARID,shuffle=.true.,deflate_level=9) )
     CALL nc_check( NF90_DEF_VAR(NCID,DAY_NAME,NF90_REAL,DAY_DIMID,DAY_VARID,shuffle=.true.,deflate_level=9) )
   else
     CALL nc_check( NF90_DEF_VAR(NCID,LON_NAME,NF90_REAL,LON_DIMID,LON_VARID) )
     CALL nc_check( NF90_DEF_VAR(NCID,LAT_NAME,NF90_REAL,LAT_DIMID,LAT_VARID) )
     CALL nc_check( NF90_DEF_VAR(NCID,DEP_NAME,NF90_REAL,DEP_DIMID,DEP_VARID) )
     CALL nc_check( NF90_DEF_VAR(NCID,DAY_NAME,NF90_REAL,DAY_DIMID,DAY_VARID) )
   endif
!Assign units attributes to coordinate var data
   CALL nc_check( NF90_PUT_ATT(NCID,LON_VARID,UNITS,LON_UNITS) )
   CALL nc_check( NF90_PUT_ATT(NCID,LAT_VARID,UNITS,LAT_UNITS) )
   CALL nc_check( NF90_PUT_ATT(NCID,DEP_VARID,UNITS,DEP_UNITS) )
   CALL nc_check( NF90_PUT_ATT(NCID,DAY_VARID,UNITS,DAY_UNITS) )

   DIMIDS = (/ LON_DIMID, LAT_DIMID, DEP_DIMID, DAY_DIMID /)
   VARIDS = (/ LON_VARID, LAT_VARID, DEP_VARID, DAY_VARID /)

END SUBROUTINE NC_4DST_DEF_DIM

! ================================== !

SUBROUTINE NC_4DST_DEF_FIELD(NCID,DIMIDS,UNAME,U,UUNITS,SCA,ADO,VARID)
INTEGER :: NCID,DIMIDS(4)
CHARACTER(LEN=*) :: UNAME,UUNITS
REAL :: U(:,:,:)
REAL :: SCA,ADO
INTEGER :: VARID
REAL :: UMAX,UMIN
INTRINSIC :: MAXVAL,MINVAL
!Define the netCDF variables
  if (VAR_I0J0K1_HDF5.eq..false.) CALL nc_check( NF90_DEF_VAR(NCID,UNAME,NF90_SHORT,DIMIDS(1:4),VARID) )
  if (VAR_I0J0K1_HDF5.eq..true.)  CALL nc_check( NF90_DEF_VAR(NCID,UNAME,NF90_SHORT,DIMIDS(1:4),VARID,shuffle=.true.,deflate_level=9) )
  UMAX=MAXVAL(U)
  UMIN=MINVAL(U)
  SCA = ((UMAX - UMIN)*1.01) / REAL(2**INT2BITS - 1)
  ADO = UMIN + 2**(INT2BITS - 1) * SCA
!Assign units attributes to the field
   CALL nc_check( NF90_PUT_ATT(NCID,VARID,UNITS,UUNITS) )
   CALL nc_check( NF90_PUT_ATT(NCID,VARID,SCA_NAME,SCA) )
   CALL nc_check( NF90_PUT_ATT(NCID,VARID,ADO_NAME,ADO) )
   call nc_check( nf90_put_att(ncid,varid,MIS_NAME,32767) )
END SUBROUTINE NC_4DST_DEF_FIELD

! ================================== !

SUBROUTINE NC_4DST_PUT_DIM(NCID,VARIDS,NLON,LONS,NLAT,LATS,NDEP,DEPS,DAYS)
INTEGER :: NCID,VARIDS(4)
INTEGER :: NLON,NLAT,NDEP
REAL :: LONS(NLON),LATS(NLAT),DEPS(NDEP),DAYS
!Write the coordinate variable data
   CALL nc_check( NF90_PUT_VAR(NCID, VARIDS(1), LONS) )
   CALL nc_check( NF90_PUT_VAR(NCID, VARIDS(2), LATS) )
   CALL nc_check( NF90_PUT_VAR(NCID, VARIDS(3), DEPS) )
   CALL nc_check( NF90_PUT_VAR(NCID, VARIDS(4), DAYS) )
END SUBROUTINE NC_4DST_PUT_DIM

! ================================== !

SUBROUTINE NC_4DST_PUT_FIELD(NCID,VARID,NLON,NLAT,NDEP,U,SCA,ADO,IN)
  INTEGER :: NCID
  INTEGER :: NLON,NLAT,NDEP
  REAL :: U(NLON,NLAT,NDEP)
  INTEGER(KIND=INT2) :: IN(NLON,NLAT,NDEP)
  REAL :: SCA,ADO
  INTEGER :: VARID
  INTEGER :: I,J,K
  integer(kind=int2) :: mytemp4(NLON,NLAT,NDEP,1)

  mytemp4=32767
  where (IN.eq.1)
    mytemp4(:,:,:,1)=int((u-ado)/sca,int2)
  end where
  CALL nc_check( NF90_PUT_VAR(NCID,VARID,mytemp4) )

END SUBROUTINE NC_4DST_PUT_FIELD

! ================================== !

SUBROUTINE NC_4DST_DEF_END(NCID)
  INTEGER :: NCID
!End define mode
  CALL nc_check( NF90_ENDDEF(NCID) )
END SUBROUTINE NC_4DST_DEF_END

! #################################################################### !

SUBROUTINE NC_FILEOPEN(fn,NCID)
character(len=*) :: fn
integer :: NCID
  CALL nc_check( nf90_open(fn, NF90_NOWRITE, NCID) )
END SUBROUTINE NC_FILEOPEN

! ================================== !

SUBROUTINE NC_4DST_INQUIRE(NCID,ndims,nvars,ngatts,unlimdimid)
integer,intent(in) :: NCID
integer,intent(out) :: ndims, nvars,ngatts,unlimdimid
  CALL nc_check( nf90_inquire(NCID, ndims, nvars, ngatts, unlimdimid) )
  ! In this case we know that there are 3 netCDF dimensions,
  ! , no global attributes, and no unlimited dimension.
  if (ndims /= 4 .or. ngatts /= 0 .or. unlimdimid /= -1) stop 2
END SUBROUTINE NC_4DST_INQUIRE


SUBROUTINE NC_4DST_INQ_DIM(NCID,NLON,NLAT,NDEP,DIMIDS,VARIDS)
INTEGER :: NCID
INTEGER :: NLON,NLAT,NDEP,NDAY
INTEGER :: DIMIDS(4),VARIDS(4)
INTEGER :: LON_DIMID,LAT_DIMID,DEP_DIMID,DAY_DIMID,LON_VARID,LAT_VARID,DEP_VARID,DAY_VARID
  ! Get the dimids of the coordinate variables
  CALL nc_check( nf90_inq_dimid(NCID, LON_NAME, LON_DIMID) )
  CALL nc_check( nf90_inq_dimid(NCID, LAT_NAME, LAT_DIMID) )
  CALL nc_check( nf90_inq_dimid(NCID, DEP_NAME, DEP_DIMID) )
  CALL nc_check( nf90_inq_dimid(NCID, DAY_NAME, DAY_DIMID) )

  ! Get the length of the coordinate variables
  CALL nc_check( nf90_Inquire_Dimension(NCID, LON_DIMID, len = NLON) )
  CALL nc_check( nf90_Inquire_Dimension(NCID, LAT_DIMID, len = NLAT) )
  CALL nc_check( nf90_Inquire_Dimension(NCID, DEP_DIMID, len = NDEP) )
  CALL nc_check( nf90_Inquire_Dimension(NCID, DAY_DIMID, len = NDAY) )

  ! Get the varids of the latitude and longitude coordinate variables.
  CALL nc_check( nf90_inq_varid(NCID, LON_NAME, LON_VARID) )
  CALL nc_check( nf90_inq_varid(NCID, LAT_NAME, LAT_VARID) )
  CALL nc_check( nf90_inq_varid(NCID, DEP_NAME, DEP_VARID) )
  CALL nc_check( nf90_inq_varid(NCID, DAY_NAME, DAY_VARID) )

  DIMIDS = (/ LON_DIMID, LAT_DIMID, DEP_DIMID, DAY_DIMID /)
  VARIDS = (/ LON_VARID, LAT_VARID, DEP_VARID, DAY_VARID /)
END SUBROUTINE NC_4DST_INQ_DIM

! ================================== !

SUBROUTINE NC_4DST_GET_DIM(NCID,VARIDS,NLON,NLAT,NDEP,LONS,LATS,DEPS,DAYS)
INTEGER :: NCID,VARIDS(4)
INTEGER :: NLON,NLAT,NDEP
REAL :: LONS(NLON),LATS(NLAT),DEPS(NDEP),DAYS
  ! Read the latitude and longitude data.
  CALL nc_check( nf90_get_var(NCID, VARIDS(1), LONS) )
  CALL nc_check( nf90_get_var(NCID, VARIDS(2), LATS) )
  CALL nc_check( nf90_get_var(NCID, VARIDS(3), DEPS) )
  CALL nc_check( nf90_get_var(NCID, VARIDS(4), DAYS) )
END SUBROUTINE NC_4DST_GET_DIM

! ================================== !

SUBROUTINE NC_4DST_GET_FIELD(NCID,UNAME,NLON,NLAT,NDEP,U,UUNITS)
INTEGER :: NCID
CHARACTER(LEN=*) :: UNAME,UUNITS
INTEGER :: NLON,NLAT,NDEP
REAL :: U(NLON,NLAT,NDEP)
INTEGER :: UVARID
REAL :: SCA,ADO
INTEGER :: I,J,K
INTEGER(KIND=INT2) :: ITEMP(NLON,NLAT,NDEP)

  ! Get the varids of the field variables.
  CALL nc_check( nf90_inq_varid(NCID, UNAME, UVARID) )
  CALL nc_check( nf90_get_att(NCID, UVARID, UNITS, UUNITS) )
  CALL nc_check( nf90_get_att(NCID, UVARID, SCA_NAME, SCA) )
  CALL nc_check( nf90_get_att(NCID, UVARID, ADO_NAME, ADO) )
  ! Get the data
  CALL nc_check( nf90_get_var(NCID, UVARID, ITEMP) )
  DO K=1,NDEP
  DO J=1,NLAT
  DO I=1,NLON
     U(I,J,K) = REAL(ITEMP(I,J,K))*SCA + ADO
  ENDDO
  ENDDO
  ENDDO
  
END SUBROUTINE NC_4DST_GET_FIELD

! ================================== !

SUBROUTINE NC_FILECLOSE(NCID) 
  INTEGER :: NCID
!Close the file
  CALL nc_check( NF90_CLOSE(NCID) )
END SUBROUTINE NC_FILECLOSE

! ================================== !

subroutine nc_check(status)
integer, intent ( in) :: status

   if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop 2
    end if
end subroutine nc_check

END MODULE NCIO
