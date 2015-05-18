MODULE STW  !  SALINITY, TEMPERATURE, WIND

! NZDAT:   DEPTH ARRAY FOR LEVITUS
! SEASON:  FILE NAME WITH RELATIVE PATH
! TXC,TYC: FILE NAME WITH RELATIVE PATH

INTEGER(2), DIMENSION(33)     ::NZDAT
CHARACTER(128), DIMENSION(4)  ::SEASON
CHARACTER(128), DIMENSION(12) ::TXC,TYC
CHARACTER(10)::HELFMT,LEVFMT

NAMELIST /LEV_FMT/   LEVFMT
NAMELIST /LEV_DEPTH/ NZDAT
NAMELIST /TS_FILE/ SEASON

NAMELIST /HEL_FMT/ HELFMT
NAMELIST /TX_FILE/ TXC
NAMELIST /TY_FILE/ TYC

END MODULE STW

!######################################################################

MODULE BATHY

!  BATHY DATA INFORMATION
!  FORMAT: 'BIG_ENDIAN' OR 'LITTLE_ENDIAN' BY NOW
!  RESOLUTION: 1,2 OR 5 (IN MINUTE)
!  BATHY FILE NAME INCLUDING FULL PATH

CHARACTER(16) ::BATHY_FMT
INTEGER       ::BATHY_RES
CHARACTER(128)::BATHY_FILE

NAMELIST /BATHYMETRY/ BATHY_FILE,BATHY_RES,BATHY_FMT

END MODULE BATHY

