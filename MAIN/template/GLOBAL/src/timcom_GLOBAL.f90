SUBROUTINE TIMCOM_GLOBAL 
USE TIMCOM_GENERAL 
! Controls calculation, performs diagnostics, and saves graphics data.
! Graphics and computer animation is done by postprocessors

! *** User-defined resolution parameters ***
USE OCN_PARA_GLOBAL
! Horizontal resolution parameters are I0,I1,I2,I3, J0,J1,J2,J3
! Vertical resolution parameters are K0,K1,K2
! NB0,NB1 are dimension parameters for SEVP elliptic solver
 
! Double precision SEVP elliptic solver arrays
USE WINDID_GLOBAL
USE XYFS_GLOBAL
USE ZCO_GLOBAL
USE GRID_VAR_GLOBAL, ONLY: RINV,RINV1,DUM0,DUM1,DUM2,XX,H,X,S,AL,AB,AC,OAC,AR,AT,SRC,CL,CB,CC,CR,CT,IE,VAR_ALLOCATE_SEVP
USE GRID_VAR_GLOBAL, ONLY: RHO,Mcrho,McTP,McS,McTPTP,McSTP,P0,PX,PY,U1,U2,V1,V2,S1,S2,T1,T2,P,ULF,VLF,SLF,TLF,EV,HV,F,TANPHI,SUMIN,VAR_ALLOCATE_B
USE GRID_VAR_GLOBAL, ONLY: U,V,W,VAR_ALLOCATE_CGRID
USE GRID_VAR_GLOBAL, ONLY: Z,ODZ,ODZW,VAR_ALLOCATE_ZFS
USE GRID_VAR_GLOBAL, ONLY: PBAR,PVAR,XBAR,UCLI,VCLI,RMSV,SBAR,TBAR,VAR_ALLOCATE_TAVG_MAIN
USE GRID_VAR_GLOBAL, ONLY: VGLO,UAVG,VAVG,SAVG,TAVG,VAR_ALLOCATE_GLO2NAB_MAIN
USE GRID_VAR_GLOBAL, ONLY: A,XDEG,YV,YVDEG,YDEG,CS,OCS,DX,ODX,DY,ODY,CSV,OCSV,DXV,ODXV,DYV,ODYV,VAR_ALLOCATE_METRIC
USE GRID_VAR_GLOBAL, ONLY: KB,IU0,IV0,IN,IU,IV,IW,VAR_ALLOCATE_BATHY
USE CLIMAT_GLOBAL
USE GRID_VAR_GLOBAL, ONLY: TNUDGE,QAVG,WAVG,SCLI,TCLI,QDAVG,SDAVG,SSURF,TSURF,SSSP,SNSP,TSSP,TNSP,VAR_ALLOCATE_CLIMAT
USE CONTROL_GLOBAL
USE SCA_GLOBAL

USE INIT_VAR_GLOBAL
USE INPUT_GLOBAL
USE REALWINDS_GLOBAL
USE PMOVE_GLOBAL
USE DIR_GLOBAL
USE OUTPUT_GLOBAL
USE RESTART_GLOBAL 

USE GRID_VAR_GLOBAL, ONLY: U2,V2,S2,P,T2
use date

CHARACTER(LEN=128) :: cfmt 
! Yu-Chiao 20111001
CHARACTER(LEN=1) :: EOS_NUM 

! **********************************************************************
!EOS_NUM=CHAR(FL_EOS_OP)
WRITE(EOS_NUM,'(I1.1)') FL_EOS_OP
!PRINT*,'test',EOS_NUM,FL_EOS_OP
!PAUSE

MPASS=DAODT/DAODT_GENERAL
! ------------------------------
! MAIN TIME INTEGRATION LOOP 100
! ------------------------------
GLOBAL: DO WHILE(.TRUE.)

1000     ITF=ITF+1
         DAYS=ITF/DAODT
         N=int(DAYS/360.+.001)
         TMP=DAYS/360.
         ! Luc: commented the following 3 lines, we now use real days in years
         ! During year, N360 goes from 0 to 360
         ! N360=int((TMP-N)*360.+.001)
         ! NYR=N+1
         NEXP=8
         DTS=DT/NEXP
         modelmjd=t0mjd+days
         call gregd2(modelmjd,year,month,day,hour,minut)
         N360=mod(int(modelmjd-t0mjd),365) ! days from the start of simulation
         NYR=int(modelmjd-t0mjd)/365

         IF ( (ITF_GENERAL .LT. TOPTS) .OR. (MOD(ITF_GENERAL,TOPTS) .EQ. 0) .OR. (ITFTO .LE. NITFTO) ) THEN
           WRITE(*,'(A,i10,A,f10.3,A,f10.3)')"GLOBAL timestep :: ", itf,", ",DAYS," days has been run -- mjd=",modelmjd
#ifdef _DEBUG_
           print *, "Before FS_GLOBAL:" 
           print *, "Variables:  maxval     minval     maxloc     minloc     sum"
           print *, "GLOBAL u2 ", maxval(u2),minval(u2),maxloc(u2),minloc(u2),sum(u2)
           print *, "GLOBAL v2 ", maxval(v2),minval(v2),maxloc(v2),minloc(v2),sum(v2)
           print *, "GLOBAL p  ", maxval(P),minval(P),maxloc(P),minloc(P),sum(P)
           print *, "GLOBAL s2 ", maxval(s2),minval(s2),maxloc(s2),minloc(s2),sum(s2)
           print *, "GLOBAL t2 ", maxval(t2),minval(t2),maxloc(t2),minloc(t2),sum(t2)
           print *, "GLOBAL u  ", maxval(u),minval(u),maxloc(u),minloc(u),sum(u)
           print *, "GLOBAL rho", "rho", maxval(rho),minval(rho),maxloc(rho),minloc(rho),sum(rho)
#endif
         ENDIF

! ==================================================
! FS marches one time step each time it is called
! ==================================================
         CALL FS_GLOBAL

#ifdef _DEBUG_
         IF ( (ITF_GENERAL .LT. TOPTS) .OR. (MOD(ITF_GENERAL,TOPTS) .EQ. 0) .OR. (ITFTO .LE. NITFTO) ) THEN
           print *, "After FS_GLOBAL:"
           print *, "Variables:  maxval     minval     maxloc     minloc     sum"
           print *, "GLOBAL u2", maxval(u2),minval(u2),maxloc(u2),minloc(u2),sum(u2)
           print *, "GLOBAL v2", maxval(v2),minval(v2),maxloc(v2),minloc(v2),sum(v2)
           print *, "GLOBAL p ", maxval(P),minval(P),maxloc(P),minloc(P),sum(P)
           print *, "GLOBAL s2", maxval(s2),minval(s2),maxloc(s2),minloc(s2),sum(s2)
           print *, "GLOBAL t2", maxval(t2),minval(t2),maxloc(t2),minloc(t2),sum(t2)
         ENDIF
#endif

! =========
! realwinds
! =========
         CALL REALWIND(TIMCOM_RUN)

! =======
! tracers
! =======
         CALL PMOVE(TIMCOM_RUN)

! ======================
! diagnostic for blow up
! ======================
         IF ( DIAG_BLOWUP() == .true. ) THEN
           WRITE(*,*)"TIMCOM_GLOBAL: DIAG_BLOWUP occured!"
           WRITE(*,*)"I must be stop it!"
           CALL PMOVE(TIMCOM_FINAL)
           STOPSIG=TIMCOM_FAILURE
           EXIT
         ENDIF


! ========
! average
! ========

         CALL AVERAGE
! ===========
! diagnostics
! ===========
! special
!!           CALL DIAG_CROSS


! ==================
! reach max timestep
! ==================
         IF (ITF.GE.MXIT) THEN
           STOPSIG=TIMCOM_FINISH   ! finish job normally
         ENDIF

! ============================================
! OUTPUT DATA WHILE REACH THE OUTPUT CONDITION
! ============================================
         IF ( STOPSIG == TIMCOM_CONTINUE ) THEN
           CALL SNAPSHOT_OUTPUT(TIMCOM_COND)
         ELSEIF ( STOPSIG == TIMCOM_FINISH ) THEN
           CALL SNAPSHOT_OUTPUT(TIMCOM_FORCE)
         ENDIF

! =============================================
! restart output or checkpoint then stop signal
! =============================================
         IF ( STOPSIG == TIMCOM_CONTINUE ) THEN
           CALL WRITE_RESTART(TIMCOM_COND)
         ELSEIF ( STOPSIG == TIMCOM_STOP ) THEN
           CALL PMOVE(TIMCOM_FINAL)
           CALL WRITE_RESTART(TIMCOM_FORCE)
           EXIT
         ENDIF

! =============
! finish signal
! =============
         IF ( STOPSIG == TIMCOM_FINISH ) THEN ! finish job normally!
           CALL PMOVE(TIMCOM_FINAL)
           CALL WRITE_RESTART(TIMCOM_FORCE)
           EXIT
         ENDIF

! ================
! re-assign signal
! ================
         IF ( STOPSIG >= TIMCOM_REASSIGN ) THEN ! re-assign the MXIT value
           MXIT=STOPSIG
           STOPSIG=TIMCOM_CONTINUE
         ENDIF

! ===================================================
! check computation complete of this general timestep
! ===================================================
         IF (MOD(ITF,MPASS).NE.0) THEN
           CYCLE
         ELSE
           EXIT
         ENDIF

ENDDO GLOBAL

CONTAINS

INCLUDE './diagnostic.f90'
INCLUDE './average.f90'
INCLUDE './realwind.f90'
INCLUDE './pmove.f90'
INCLUDE './output.f90'
INCLUDE './write_restart.f90'

END SUBROUTINE TIMCOM_GLOBAL 
