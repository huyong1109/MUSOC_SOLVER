! *******************PROGRAM TIMCOM************************************
! TaIwanese Multi-scale Community Ocean Models (TIMCOM)   January, 2009
! with BIR (domain decomposition) for periodic solver
! *********************************************************************
PROGRAM MAIN 
! Controls calculation, performs diagnostics, and saves graphics data.
! Graphics and computer animation is done by postprocessors

USE date
USE TIMCOM_GENERAL
! Yu-Chiao 20192011
!USE SCA_MASTR

! **********************************************************************

! -------------------------------
! GET ARGUMENT FOR RESTART SIGNAL
! -------------------------------
  CALL TIMCOM_GETARGS

! -------------
! OPEN LOG FILE
! ------------- 
  CALL TIMCOM_OPENLOG
  write(FNO_LOG,*) 'Simulation member : ',member

! --------------
! INITIALIZATION
! --------------
  CALL INIT

! -------------
! CHECK RESTART
! -------------
  CALL RESTART_CHECK

! ------------------------------
! MAIN TIME INTEGRATION LOOP 100
! ------------------------------
time: DO WHILE(.TRUE.)

      ITF_GENERAL=ITF_GENERAL+1
      ITFTO=ITFTO+1

      IF ( (ITF_GENERAL .LT. TOPTS) .OR. (MOD(ITF_GENERAL,TOPTS) .EQ. 0) .OR. (ITFTO .LE. NITFTO) ) THEN
        CALL SHOW_STEPS_TO_DEV(6)
        CALL SHOW_STEPS_TO_DEV(FNO_LOG)
      ENDIF

    ! reset the checkpoint signal to zero(0 for no, 1 for yes)
      CHKPTSIG=0

      CALL TIMCOM

    ! check checkpoint
      CALL CHECKPOINT_CHECK

! ===============
! job stop signal
! ===============
      !#New Way# CALL TIMCOM_STOPSIG_SIGNAL
      ! STOPSIG  = TIMCOM_FAILURE  =-2 : it's blow up, job failure
      ! STOPSIG  = TIMCOM_CONTINUE =-1 : continue running the job
      ! STOPSIG  = TIMCOM_FINISH   = 0 : finish job normally
      ! STOPSIG  = TIMCOM_STOP     = 1 : checkpoint and stop the job
      ! STOPSIG >= TIMCOM_REASSIGN = 2 : re-assign the timestep(MXIT value)
      IF ( STOPSIG == TIMCOM_FAILURE ) THEN
         EXIT  ! job failure
      ELSEIF ( STOPSIG == TIMCOM_FINISH ) THEN
         EXIT  ! job successful complete
      ELSEIF ( STOPSIG == TIMCOM_STOP ) THEN
         EXIT  ! job checkpoint then stop
      ENDIF

    ! check every 5 step   
      IF (MOD(ITF_GENERAL,5).EQ.0) THEN
        CALL TIMCOM_STOPSIG_READ
      ENDIF

ENDDO time 

!---------------
! CLOSE LOG FILE
!---------------
  CALL TIMCOM_CLOSELOG

!------------------------
! RETURN SIGNAL TO SYSTEM
!------------------------
  CALL TIMCOM_STOPSIG_SIGNAL

CONTAINS

SUBROUTINE SHOW_STEPS_TO_DEV(FNO)
  INTEGER :: FNO
  WRITE(FNO,*)""
  WRITE(FNO,*)"==========================================="
  WRITE(FNO,*)" PROGRAM EXECUTION STEP:",ITF_GENERAL
  WRITE(FNO,*)"==========================================="
  WRITE(FNO,*)""
END SUBROUTINE SHOW_STEPS_TO_DEV

SUBROUTINE RESTART_CHECK
  character(20)::tempostring
  IF ( RESTARTSIG == 1 ) THEN
    write(tempostring,'(A,i3.3)')"restart_timestep",member
    OPEN(21,FILE=tempostring,FORM="FORMATTED",STATUS="OLD",IOSTAT=ISTAT)
    IF ( ISTAT /= 0 ) THEN
      WRITE(*,*)""
      WRITE(*,*)"**************************************************************************************************"
      WRITE(*,*)"ERROR: timestep file OPEN failed , I can't restart the run!"
      WRITE(*,*)"       You might need to reference to restart.complete (for example: TEMP_GLOBAL/restart.complete)"
      WRITE(*,*)"       and generate the file by your self."
      WRITE(*,*)"**************************************************************************************************"
      WRITE(*,*)""
      STOP
    ELSE
      READ(21,*)ITF_GENERAL
      CLOSE(21)
    ENDIF
  ENDIF
END SUBROUTINE RESTART_CHECK

SUBROUTINE CHECKPOINT_CHECK
  character(20)::tempostring
  IF ( CHKPTSIG == 1 ) THEN
    write(tempostring,'(A,i3.3)')"restart_timestep",member
    OPEN(21,FILE=tempostring,FORM="FORMATTED",STATUS="REPLACE")
    WRITE(21,"(I10)")ITF_GENERAL
    CLOSE(21)
  ENDIF
END SUBROUTINE CHECKPOINT_CHECK

END PROGRAM MAIN 
