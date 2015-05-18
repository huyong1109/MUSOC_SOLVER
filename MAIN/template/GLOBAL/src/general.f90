MODULE TIMCOM_GENERAL
use interfext, interfinit => init, interfgetfield => getfield
!MEMBER NUMBER
  INTEGER::member

!PARAMETER
  REAL,PARAMETER::G=980.

!TITLE
  CHARACTER(10)::CASE_NAME
  CHARACTER(63)::DSCRIB
  CHARACTER(256)::WORKDIR
  CHARACTER(256)::ATMOSFILE
  character(256)::icfile
  INTEGER::RUNS
  INTEGER::ISAV
  INTEGER::SCRNOUT,TOPTS
  INTEGER::CPSAV=-5    ! CPSAV > 0, every CPSAV timestep ; CPSAV < 0, every ABS(CPSAV) days
  INTEGER::t0year,t0month,t0day,t0hour,t0minut
  REAL::t0mjd
  INTEGER::ics_top,ics_bottom
  NAMELIST /CASE_INFO/ CASE_NAME,DSCRIB

!RUNNING SIGNAL
  INTEGER,PARAMETER :: &
  TIMCOM_INIT=0,TIMCOM_RUN=1,TIMCOM_FINAL=2, &
  TIMCOM_COND=0,TIMCOM_FORCE=1

!GLOBAL VARIABLES FOR ATMOSPHERIC FORCING FIELDS
type(field2D) :: SSTfield,&
     shortwavedown,longwaveup,sensibleup,latentup,Qup, &
     MomentfluxX,MomentfluxY,Precip,EvapField,PresMSL, &
     WindU10m,WindV10m,AirTemperature2m,DewTemperature2m, &
     RelHum2m,CloudCoverage
integer :: windtype, heattype, raintype, salttype, humtype, cloudtype
logical :: bulk

!VARIABLES OUTPUT SIGNAL
  LOGICAL::VAR_I0J0K1_OUTPUT=.TRUE.,&
           VAR_I0J0K1_HDF5=.FALSE.,&
           UVAR_OUTPUT=.TRUE., &
           VVAR_OUTPUT=.TRUE., &
           SVAR_OUTPUT=.TRUE., &
           TVAR_OUTPUT=.TRUE., &
           PVAR_OUTPUT=.TRUE.
  LOGICAL::VAR_I0J0K0_OUTPUT=.TRUE.,&
           WVAR_OUTPUT=.TRUE.
  LOGICAL::VAR_FLUXES_OUTPUT=.FALSE.

!LOG FILE
  INTEGER,PARAMETER :: FNO_LOG=14

!ITF GENERAL
  INTEGER :: ITF_GENERAL=0
  REAL :: DAODT_GENERAL=0.

!RESTART TERMOUT TIMESTEP
  INTEGER :: ITFTO=0
  INTEGER :: NITFTO=10

!RESTART_SIGNAL RESTARTSIG
  INTEGER :: RESTARTSIG = 0
  INTEGER :: icSIG = 0
  INTEGER :: smoothSIG = 0
  INTEGER :: stabilizeSIG = 0

!CHECKPOINT_SIGNAL CHKPTSIG
  INTEGER :: CHKPTSIG = 0

!STOP_SIGNAL
  INTEGER,PARAMETER :: &
  TIMCOM_FAILURE=-2,   &
  TIMCOM_CONTINUE=-1,  &
  TIMCOM_FINISH=0,     &
  TIMCOM_STOP=1,       &
  TIMCOM_REASSIGN=2
  INTEGER :: STOPSIG = -1

!MPASS
  INTEGER :: MPASS = 1
  REAL    :: FLT1  = 0.
  REAL    :: FLT2  = 1.
  NAMELIST /MULTI_TIMESTEP/ MPASS,FLT1,FLT2
CONTAINS

SUBROUTINE TIMCOM_GETARGS
  CHARACTER(128) :: ARG,tempstr
  INTEGER,EXTERNAL :: IARGC
  INTEGER :: istat,i,j

  do i = 1,iargc()
    CALL getarg(i,arg); read(arg,*,iostat=istat) tempstr
    if (istat /= 0 ) then
       write(*,*)"something went wrong with the command line arguments,"
       write(*,*)"I'd better stop !"
       STOPSIG=TIMCOM_FAILURE
       call TIMCOM_STOPSIG_SIGNAL
    endif
    j=index(tempstr,"restart")
       if (j .gt. 0) then
         RESTARTSIG=1
       endif
    j=index(tempstr,"member")
       if (j .gt. 0) then
         read(tempstr(j+7:len(tempstr)),*) member
       endif
    j=index(tempstr,"ic")
       if (j .gt. 0) then
         icSIG=1
         read(tempstr(j+3:len(tempstr)),*) icfile
       endif
    j=index(tempstr,"smooth")
       if (j .gt. 0) then
       !  read(tempstr(j+7:len(tempstr)),*) smoothSIG
         smoothSIG=1
       endif
    j=index(tempstr,"stabilize")
       if (j .gt. 0) then
         stabilizeSIG=1
       endif
     !examine for other command line arguments here
  end do
END SUBROUTINE TIMCOM_GETARGS

SUBROUTINE TIMCOM_OPENLOG
CHARACTER(LEN=160) :: FILENAME
!-----------------------
! OPEN General log files
!-----------------------
! Run history data

  !FILENAME=TRIM(WORKDIR)//'run.log'
  write(FILENAME,"(A,i3.3,A)") 'member',member,'.log'
  OPEN(FNO_LOG,file=FILENAME)

END SUBROUTINE TIMCOM_OPENLOG

!######################################

SUBROUTINE TIMCOM_STOPSIG_READ
  LOGICAL :: alive
  INTEGER :: istat
  INQUIRE(FILE="stoprun",EXIST=alive)
  IF ( alive ) THEN
    OPEN(42,file='stoprun')
    READ(42,"(I7)",IOSTAT=istat) STOPSIG 
    CLOSE(42,STATUS="DELETE")
    IF ( istat /= 0 ) THEN
       WRITE(*,*)"I got wrong stoprun signal!"
       WRITE(*,*)"stoprun signal must be"
       WRITE(*,*)" 0 : finish the job"
       WRITE(*,*)" 1 : checkpoint then stop the job"
       WRITE(*,*)" 2 or above : re-assign the max timestep"
       WRITE(*,*)""
       WRITE(*,*)"I'll continue to run the job"
       STOPSIG = TIMCOM_CONTINUE
    ENDIF
    IF ( STOPSIG < TIMCOM_FINISH ) THEN
       STOPSIG = TIMCOM_CONTINUE
    ENDIF
  ENDIF
END SUBROUTINE TIMCOM_STOPSIG_READ

!######################################

SUBROUTINE TIMCOM_CLOSELOG
  LOGICAL :: lopen
  INQUIRE(UNIT=FNO_LOG,OPENED=lopen)
  IF ( lopen ) THEN
    CLOSE(FNO_LOG)
  ENDIF
END SUBROUTINE TIMCOM_CLOSELOG

!######################################

SUBROUTINE TIMCOM_STOPSIG_SIGNAL
  SELECT CASE (STOPSIG)
  CASE (-1)
    WRITE(*,*)"running the job continully"
  CASE (0)
    CALL TIMCOM_CLOSELOG
    WRITE(*,*)"job successful complete"
    STOP 0 
  CASE (1)
    CALL TIMCOM_CLOSELOG
    WRITE(*,*)"checkpoint, then stop the job!"
    STOP 1 
  CASE (-2)
    CALL TIMCOM_CLOSELOG
    WRITE(*,*)"It's blow-up, job failure!"
    STOP 2 
  CASE DEFAULT
    WRITE(*,*)"re-assign the timestep"
  END SELECT
END SUBROUTINE TIMCOM_STOPSIG_SIGNAL

END MODULE TIMCOM_GENERAL


pure function dnum (string)
  character(len=*), intent(in) :: string
  real(KIND(0.D0))             :: dnum
  read(string, *) dnum
end function dnum

