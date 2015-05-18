MODULE INPUT_GLOBAL
#include "namelist.in_GLOBAL"
USE TIMCOM_GENERAL 
USE OCN_PARA_GLOBAL
USE GRID_VAR_GLOBAL, ONLY: VBK,HBK,INFX,INFY
USE CONTROL_GLOBAL
USE SCA_GLOBAL
USE GRID_VAR_GLOBAL, ONLY: RHO,Mcrho,McTP,McS,McTPTP,McSTP,P0,PX,PY,U1,U2,V1,V2,S1,S2,T1,T2,P,ULF,VLF,SLF,TLF,EV,HV,F,TANPHI,SUMIN
USE GRID_VAR_GLOBAL, ONLY: A,XDEG,YV,YVDEG,YDEG,CS,OCS,DX,ODX,DY,ODY,CSV,OCSV,DXV,ODXV,DYV,ODYV
USE XYFS_GLOBAL
USE GRID_VAR_GLOBAL, ONLY: Z,ODZ,ODZW
USE GRID_VAR_GLOBAL, ONLY: KB,IU0,IV0,IN,IU,IV,IW
USE GRID_VAR_GLOBAL, ONLY: SCR
USE CLIMAT_GLOBAL
USE GRID_VAR_GLOBAL, ONLY: TNUDGE,QAVG,WAVG,SCLI,TCLI,QDAVG,SDAVG,SSURF,TSURF,SSSP,SNSP,TSSP,TNSP
USE GRID_VAR_GLOBAL, ONLY:TAUX,TAUY,VAR_ALLOCATE_TAUXY
USE GRID_VAR_GLOBAL, ONLY: SXY,TXY,SXYCLI,TXYCLI
USE GRID_VAR_GLOBAL, ONLY: RINV,RINV1,DUM0,DUM1,DUM2,XX,H,X,S,AL,AB,AC,OAC,AR,AT,SRC,CL,CB,CC,CR,CT,IE
USE DIR_GLOBAL, ONLY: IDIR
USE OUTPUT_GLOBAL
!DTRAC ONLY
USE GRID_VAR_GLOBAL, ONLY:U,V
use initfile
use ufileformat

IMPLICIT NONE

CHARACTER(LEN=256),PRIVATE :: FNAME

INTEGER,PRIVATE,PARAMETER :: FN_PARAM=20, &
                             FN_TURBMIX=21,  &
                             FN_RUNDATA=22,  &
                             FN_ZKB=23,      &
                             FN_winds=24,    &
                             FN_initial=25,  &
                             FN_EVP=26,      &
                             FN_boundaries=27,&
                             FN_eos=29

 CONTAINS

SUBROUTINE READ_INPUT_PARAMETERS
  INTEGER :: I_I,J_J,K_K
  CHARACTER(22)::tempfname
  character(256)::tempo_str
  integer :: i,j,k,tempo_int
  real :: tempo_real
  real, allocatable ::tempo_array2(:,:), tempo_array3(:,:,:)

  write(tempfname,"(A,i3.3,2A)") "namelist",member,".inp_","GLOBAL"

  call getInitValue(tempfname,'COORD_SYS',COORD_SYS)
  call getInitValue(tempfname,'X0DEG',X0DEG)
  call getInitValue(tempfname,'Y0DEG',Y0DEG)
  call getInitValue(tempfname,'DXMNUT',DXMNUT)
  call getInitValue(tempfname,'DXDEG',DXDEG)
  call getInitValue(tempfname,'DYDEG',DYDEG)
  call getInitValue(tempfname,'DYDX',DYDX)
  call getInitValue(tempfname,'DX0STD',DX0STD)

  call getInitValue(tempfname,'DAODT',DAODT)

  if (presentinitvalue(tempfname,'DM0_CONST')) then
    call getInitValue(tempfname,'DM0_CONST',tempo_real)
    allocate(dm0(i0,j0,k0-1))
    DM0=tempo_real
  elseif (presentinitvalue(tempfname,'DM0_FILE')) then
    call getInitValue(tempfname,'DM0_FILE',tempo_str)
    call uload(trim(tempo_str),dm0,tempo_real)
  else
    write(*,*) "GLOBAL init file : couldn't find DM0 (constant or filename)"
  endif
  if (presentinitvalue(tempfname,'DE0_CONST')) then
    call getInitValue(tempfname,'DE0_CONST',tempo_real)
    allocate(de0(i0,j0,k0-1))
    DE0=tempo_real
  elseif (presentinitvalue(tempfname,'DE0_FILE')) then
    call getInitValue(tempfname,'DE0_FILE',tempo_str)
    call uload(trim(tempo_str),DE0,tempo_real)
  else
    write(*,*) "GLOBAL init file : couldn't find DE0 (constant or filename)"
  endif
  if (presentinitvalue(tempfname,'VBK_CONST')) then
    call getInitValue(tempfname,'VBK_CONST',tempo_real)
    allocate(vbk(2:i0-1,2:j0-1,k0-2))
    VBK=tempo_real
  elseif (presentinitvalue(tempfname,'VBK_FILE')) then
    call getInitValue(tempfname,'VBK_FILE',tempo_str)
    call uload(trim(tempo_str),VBK,tempo_real)
  else
    write(*,*) "GLOBAL init file : couldn't find VBK (constant or filename)"
  endif
  call getInitValue(tempfname,'DMZ0',DMZ0)
  call getInitValue(tempfname,'TL',TL)
  call getInitValue(tempfname,'RZMX',RZMX)
  call getInitValue(tempfname,'DRAG',DRAG)

  call getInitValue(tempfname,'N0',N0)
  call getInitValue(tempfname,'NBIR',NBIR)
  call getInitValue(tempfname,'MXMASK',MXMASK)
  call getInitValue(tempfname,'TOLERANCE',TOLERANCE)
  call getInitValue(tempfname,'FL_EVP_STP',FL_EVP_STP)

  call getInitValue(tempfname,'UINFLOW',UINFLOW)
  call getInitValue(tempfname,'VINFLOW',VINFLOW)

  call getInitValue(tempfname,'FL_TURB_H',FL_TURB_H)
  call getInitValue(tempfname,'FL_TURB_V',FL_TURB_V)
  call getInitValue(tempfname,'FL_TURB_V_ADD',FL_TURB_V_ADD)
  call getInitValue(tempfname,'ITF_PHYS',ITF_PHYS)

  call getInitValue(tempfname,'FL_TRACER_ON',FL_TRACER_ON)

  call getInitValue(tempfname,'FL_SURF_EP',FL_SURF_EP)
  call getInitValue(tempfname,'FL_EVAP1',FL_EVAP1)

  call getInitValue(tempfname,'FL_FNEWFOLD_CALC',FL_FNEWFOLD_CALC)
  call getInitValue(tempfname,'FL_NUDGE_HMEAN_T',FL_NUDGE_HMEAN_T)
  call getInitValue(tempfname,'FL_NUDGE_HMEAN_S',FL_NUDGE_HMEAN_S)
  call getInitValue(tempfname,'NUDGING_HMEAN',NUDGING_HMEAN)
  call getInitValue(tempfname,'FL_NUDGE_OP_T',FL_NUDGE_OP_T)
  call getInitValue(tempfname,'FL_NUDGE_OP_S',FL_NUDGE_OP_S)
  call getInitValue(tempfname,'NUDGE_W',NUDGE_W)
  call getInitValue(tempfname,'FL_ARBR_P0',FL_ARBR_P0)
  if (presentinitvalue(tempfname,'TAU_CONST')) then
    call getInitValue(tempfname,'TAU_CONST',tempo_real)
    allocate(tau2D(i0,j0))
    TAU2D=tempo_real
  elseif (presentinitvalue(tempfname,'TAU_FILE')) then
    call getInitValue(tempfname,'TAU_FILE',tempo_str)
    call uload(trim(tempo_str),tempo_array2,tempo_real)
    allocate(tau2D(i0,j0))
    tau2D=tempo_array2
    deallocate(tempo_array2)
  else
    write(*,*) "GLOBAL init file : couldn't find TAU (constant or filename)"
  endif
  if (presentinitvalue(tempfname,'TAUN_CONST')) then
    call getInitValue(tempfname,'TAUN_CONST',tempo_real)
    allocate(tau3D(i0,j0,k0-1))
    TAU3D=tempo_real
  elseif (presentinitvalue(tempfname,'TAUN_FILE')) then
    call getInitValue(tempfname,'TAUN_FILE',tempo_str)
    call uload(trim(tempo_str),tempo_array3,tempo_real)
    allocate(tau3D(i0,j0,k0-1))
    tau3D=tempo_array3
    deallocate(tempo_array3)
  else
    write(*,*) "GLOBAL init file : couldn't find TAUN (constant or filename)"
  endif

  call getInitValue(tempfname,'FL_SPL_ON',FL_SPL_ON)
  call getInitValue(tempfname,'FL_SWAMP_RE',FL_SWAMP_RE)

  call getInitValue(tempfname,'FLTW',FLTW)
  call getInitValue(tempfname,'WRAF',WRAF)
  call getInitValue(tempfname,'FL_BI_FIL',FL_BI_FIL)

  call getInitValue(tempfname,'FL_INI_VIS',FL_INI_VIS)
  call getInitValue(tempfname,'FL_BK_VIS',FL_BK_VIS)

  call getInitValue(tempfname,'FL_PMOVE_ON',FL_PMOVE_ON)
  call getInitValue(tempfname,'FL_PMOVE_NEWDAY',FL_PMOVE_NEWDAY)

  call getInitValue(tempfname,'RIVERFILE',RIVERFILE)

  if (smoothSIG.eq.1) then
    call getInitValue(tempfname,'ICS_TOP',ICS_TOP)
    call getInitValue(tempfname,'ICS_BOTTOM',ICS_BOTTOM)
  endif

! Compute some variables from reading variables
  DAODT_GENERAL=DAODT
  do i=1,i0
  do j=1,j0
    tau2D(i,j)=1./(tau2d(i,j)*DAODT)
    do k=1,k0-1
      tau3D(i,j,k)=1./(tau3D(i,j,k)*DAODT)
    end do
  end do
  end do
  ORZMX=1./RZMX 
  OFLTW=1.-2.*FLTW
  FLT1=MPASS/(.5*DAODT)
  FLT2=1.-FLT1
  IF ( RUNS .LT. 0 ) THEN
     MXIT=ABS(RUNS)*INT(DAODT+.0001)
  ELSE
     MXIT=RUNS
  ENDIF

  ! check if ITF reaches max timestep
  IF (ITF.GE.MXIT) THEN
    WRITE(*,*)"INPUT_GLOBAL: ITF ",ITF, ">= MAX TIMESTEPS", MXIT
    WRITE(*,*)"I'd be stop better!"
    STOPSIG=TIMCOM_FINISH
    CALL TIMCOM_STOPSIG_SIGNAL
  END IF

!++++++ need to be calculated +++++
  ITFDAY=int(DAODT+.0001)

  DT=172800./DAODT
  ODT=1./DT

!  PRN=DE0/DM0

  IF ( SCRNOUT .LT. 0 ) THEN
    TOPTS=ABS(SCRNOUT)*ITFDAY
  ELSEIF ( SCRNOUT .GT. 0 ) THEN
    TOPTS=SCRNOUT
  ELSE
    TOPTS=0
  ENDIF
!+++++++++++++++++++++++++++++++++

  FNAME=TRIM(WORKDIR)//TRIM(IDIR)//'GLOBAL'
  OPEN(FN_PARAM,CONVERT='BIG_ENDIAN',FILE=TRIM(FNAME),form='UNFORMATTED')
  READ(FN_PARAM)I_I,J_J,K_K 
  ! LATERAL_BOUND_FLAGS
  READ(FN_PARAM) LOPENN,LOPENS,LOPENE,LOPENW,FL_INFL_N,FL_INFL_S,FL_INFL_E,FL_INFL_W
  ! TSC_FLAGS except FL_TRACER_ON
  READ(FN_PARAM) FL_RL_TS,TMAX,FL_EOS_OP
  ! WIND_FLAGS
  READ(FN_PARAM) FL_WD_ON,FL_RL_WD
  CLOSE(FN_PARAM)

END SUBROUTINE READ_INPUT_PARAMETERS

SUBROUTINE READ_INPUT_FILES
  INTEGER :: I,J,K,N,M
  real :: d_r_var1,d_r_var2(j0),d_r_var3(j0,k1),d_r_var4(j0,k1,12)
  !d_r_var1 = dummy real variable
  character :: d_c_var1*5,d_c_var2*63
  real :: TMP

  print *, "Read GLOBAL's input files"
  CALL VAR_ALLOCATE_TAUXY ! this also allocate heat flux variables (qdot etc)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! read input files in 4 steps:                                 !
  ! 1. read file ZKB for KB, IN, DEPTH, ...                      !
  ! 2. read file initial for S and T                             !
  ! 3. read file winds for wind shear stress TAUX, TAUY          !
  ! 4. read file boundaries for velocities or surface conditions !
  ! 5. read file TURBMIX for VBK, HBK                            !
  ! 6. read file RUNDATA for F, EV, HV, ...                      !
  ! 7. read file EVP for AL, AR, ...                             !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! 1. read file ZKB for KB, IN, DEPTH, ...
  FNAME=TRIM(WORKDIR)//TRIM(IDIR)//'ZKB'
  OPEN(FN_ZKB,CONVERT='BIG_ENDIAN',FILE=TRIM(FNAME),form='unformatted')
  READ(FN_ZKB) Z,ODZ,ODZW
  READ(FN_ZKB) KB,IN,IU,IV,IW,IU0,IV0
  ! Depth
  READ(FN_ZKB) ((SCR(I,J,1),I=1,I0),J=1,J0)
  DO J=2,J0-1
    DO I=2,I0-1
      ! Convert to kilometers
      SCR(I,J,1)=1.E-5*SCR(I,J,1)
    END DO
  END DO
  CLOSE(FN_ZKB)

#ifdef _DEBUG_  
  print *, "GLOBAL zkb "
  print *, "  iu  ", MAXVAL(iu),MINVAL(iu),MAXLOC(iu),SUM(iu(1:i0-1,1:j0-1,1:k1))
  print *, "  iu0 ", MAXVAL(iu0),MINVAL(iu0),MAXLOC(iu0),SUM(iu0(1:i0-1,1:j0-1))
  print *, "  iv  ", MAXVAL(iv),MINVAL(iv),MAXLOC(iv),SUM(iv(1:i0-1,1:j0-1,1:k1))
  print *, "  iv0 ", MAXVAL(iv0),MINVAL(iv0),MAXLOC(iv0),SUM(iv0(1:i0-1,1:j0-1))
  print *, "  z   ", maxval(z),minval(z),maxloc(z),minloc(z),sum(z) !ok
  print *, "  odz ", maxval(odz),minval(odz),maxloc(odz),minloc(odz),sum(odz) !ok
  print *, "  odzw", maxval(odzw),minval(odzw),maxloc(odzw),minloc(odzw),sum(odzw) !ok
  print *, "  kb  ", maxval(kb),minval(kb),maxloc(kb),minloc(kb),sum(kb)
  print *, "  in  ", maxval(in),minval(in),maxloc(in),minloc(in),sum(in)  !ok
  print *, "  iu  ", maxval(iu),minval(iu),maxloc(iu),minloc(iu),sum(iu)
  print *, "  iv  ", maxval(iv),minval(iv),maxloc(iv),minloc(iv),sum(iv)  !ok
  print *, "  iw  ", maxval(iw),minval(iw),maxloc(iw),minloc(iw),sum(iw)  !ok
  print *, "  iu0 ", maxval(iu0),minval(iu0),maxloc(iu0),minloc(iu0),sum(iu0) !ok
  print *, "  iv0 ", maxval(iv0),minval(iv0),maxloc(iv0),minloc(iv0),sum(iv0) !ok
  print *, "  SCR ", maxval(SCR(:,:,1)),minval(SCR(:,:,1)),maxloc(SCR(:,:,1)),minloc(SCR(:,:,1)),sum(SCR(:,:,1))
#endif

! Yu-Chiao 20110919
    FNAME=TRIM(WORKDIR)//TRIM(IDIR)//'eos'
    OPEN(FN_eos,CONVERT='BIG_ENDIAN',FILE=TRIM(FNAME),FORM='UNFORMATTED')

  ! 2. read file initial for S and T
! Yu-Chiao 10152011
#ifdef FLAG_TS_T
      !E.G. ILE2D
      !E.G. DTRAC CASE ONLY CONSIDERS T (T IS ALSO AN IDEALIZED CASE)
      REAL_T: SELECT CASE(FL_RL_TS)
      CASE(0) REAL_T
        FNAME=TRIM(WORKDIR)//TRIM(IDIR)//'initial'
        OPEN(FN_initial,CONVERT='BIG_ENDIAN',FILE=TRIM(FNAME),FORM='UNFORMATTED')
        READ(FN_initial) TCLI ! I.E. DTRAC
        CLOSE(FN_initial)
      CASE DEFAULT REAL_T
          ! NEED FURTHER TREATMENT
      END SELECT REAL_T
#endif
#ifdef FLAG_TS_T
#ifdef FLAG_TS_S
      REAL_TS: SELECT CASE(FL_RL_TS)
      CASE(0) REAL_TS
          !IDEALIZED T&S CAN BE FURTHER IMPLEMENTED
      CASE DEFAULT REAL_TS
        FNAME=TRIM(WORKDIR)//TRIM(IDIR)//'initial'
        OPEN(FN_initial,CONVERT='BIG_ENDIAN',FILE=TRIM(FNAME),FORM='UNFORMATTED')
        READ(FN_initial) SCLI,TCLI
! Yu-Chiao 20110919 
          IF (FL_EOS_OP.EQ.4) THEN
            READ(FN_eos) Mcrho,McTP,McS
          ELSE IF (FL_EOS_OP.EQ.5) THEN
            READ(FN_eos) Mcrho,McTP,McS,McTPTP
          ELSE IF (FL_EOS_OP.EQ.6) THEN
            READ(FN_eos) Mcrho,McTP,McS,McTPTP,McSTP
          END IF
!          PRINT*,'Mcrho in MAIN',SUM(Mcrho)
!          PAUSE
      END SELECT REAL_TS
#ifdef _DEBUG_    
      print *, "GLOBAL annual"
      print *, "  scli", maxval(scli),minval(scli),maxloc(scli),minloc(scli),sum(scli)
      print *, "  tcli", maxval(tcli),minval(tcli),maxloc(tcli),minloc(tcli),sum(tcli)
      CLOSE(FN_initial)
#endif
#endif
#endif


  ! 3. read file winds for wind shear stress TAUX, TAUY
  !write(*,*) "windtype=",windtype
  if (windtype.eq.-1) then
     IF ((FL_WD_ON == 1) .AND. (FL_RL_WD == 1)) THEN
       FNAME=TRIM(WORKDIR)//TRIM(IDIR)//'winds'
       OPEN(FN_winds,CONVERT='BIG_ENDIAN',FILE=TRIM(FNAME),FORM='UNFORMATTED')
       DO N=1,12
         ! Read wind stress (dynes/cm-cm)
         READ(FN_winds) ((TAUX(I,J,N),I=1,I2),J=1,J2),((TAUY(I,J,N),I=1,I2),J=1,J2)
#ifdef _DEBUG_
         print *, "local taux",maxval(taux(:,:,n)),minval(taux(:,:,n)),maxloc(taux(:,:,n)),minloc(taux(:,:,n)),sum(taux(:,:,n))
         print *, "local tauy",maxval(tauy(:,:,n)),minval(tauy(:,:,n)),maxloc(tauy(:,:,n)),minloc(tauy(:,:,n)),sum(tauy(:,:,n))
#endif
       END DO
       CLOSE(FN_winds)
     END IF
  !elseif (windtype.eq.0)
     ! we have the wind speed U10M and V10M ,  we will convert and interpolate in space and time later
  !elseif (windtype.eq.1)
     ! taux and tauy are already there in Momentflux{X,Y}, we will interpolate in space and time later
  !else
     ! write(*,*) "We got a serious problem, windtype is not set correctly. Timcom should exit now."
  end if

 
  ! 4. read file boundaries for velocities or surface conditions
  FNAME=TRIM(WORKDIR)//TRIM(IDIR)//'boundaries'
  OPEN(FN_boundaries,CONVERT='BIG_ENDIAN',FILE=TRIM(FNAME),FORM='UNFORMATTED')
  IF(TRIM(CASE_NAME) == 'ILE2D') THEN
  ! specified inflows
  READ(FN_boundaries) V,V1,V2,VLF
  ELSEIF(TRIM(CASE_NAME) == 'DTRAC') THEN
  IF (RESTARTSIG.EQ.0) READ(FN_boundaries)   &
     ((U(I1,J,K),J=1,J0),K=1,K1),((U2(I0,J,K),J=1,J0),K=1,K1),  &
     ((U(1,J,K) ,J=1,J0),K=1,K1),((U2(1,J,K) ,J=1,J0),K=1,K1),  &
     ((V(I,J1,K),I=1,I0),K=1,K1),((V2(I,J0,K),I=1,I0),K=1,K1),  &
     ((V(I,1,K) ,I=1,I0),K=1,K1),((V2(I,1,K) ,I=1,I0),K=1,K1) 
  ENDIF
  
  ! Surface climatology and specified inflows
#ifdef FLAG_TS_T
#ifdef FLAG_TS_S
  IF (FL_RL_TS>=1) THEN
    READ(FN_boundaries) SSURF,TSURF,SSSP,SNSP,TSSP,TNSP
    READ(FN_boundaries) SXYCLI,TXYCLI
#ifdef _DEBUG_  
    print *, "GLOBAL bound"
    print *, "  ssurf", maxval(ssurf),minval(ssurf),maxloc(ssurf),minloc(ssurf),sum(ssurf)
    print *, "  tsurf", maxval(tsurf),minval(tsurf),maxloc(tsurf),minloc(tsurf),sum(tsurf)
    print *, "  sssp ", maxval(sssp),minval(sssp),maxloc(sssp),minloc(sssp),sum(sssp)
    print *, "  tssp ", maxval(tssp),minval(tssp),maxloc(tssp),minloc(tssp),sum(tssp)
#endif
  ENDIF 
#endif
#endif
  CLOSE(FN_boundaries)
  
  ! 5. read file TURBMIX for VBK, HBK
  ! STILL NEED SOME MODIFICATION
  !IF (LTURB == 1) THEN
  IF (FL_TURB_V_ADD == 1) THEN
    FNAME=TRIM(WORKDIR)//TRIM(IDIR)//'TURBMIX'
    OPEN(FN_TURBMIX,CONVERT='BIG_ENDIAN',FILE=TRIM(FNAME),FORM='UNFORMATTED')
    READ(FN_TURBMIX) HBK,HBK
    CLOSE(FN_TURBMIX)

#ifdef _DEBUG_
    !print *, "GLOBAL vbk", maxval(VBK),minval(VBK),maxloc(VBK),minloc(VBK),sum(VBK)
    print *, "GLOBAL hbk", maxval(hbk),minval(hbk),maxloc(hbk),minloc(hbk),sum(hbk)
#endif

  END IF

  ! 6. read file RUNDATA for F, EV, HV, ...
  FNAME=TRIM(WORKDIR)//TRIM(IDIR)//'RUNDATA'
  OPEN(FN_RUNDATA,CONVERT='BIG_ENDIAN',FILE=TRIM(FNAME),form='UNFORMATTED')
  ! read variables of BGLO in grid_var.f90
  ! read variables of METRICGLO in grid_var.f90
  READ(FN_RUNDATA) F,TANPHI,EV,HV
  READ(FN_RUNDATA) XDEG,YVDEG,YDEG,CS,CSV,OCS,OCSV,DX,DXV,ODX,ODXV,DY,DYV,ODY,ODYV
  CLOSE(FN_RUNDATA)

#ifdef _DEBUG_ 
  print *, "GLOBAL rundata"
  print *, "  f     ", maxval(f),minval(f),maxloc(f),minloc(f),sum(f)
  print *, "  tanphi", maxval(tanphi),minval(tanphi),maxloc(tanphi),minloc(tanphi),sum(tanphi)
  print *, "  ev    ", maxval(ev),minval(ev),maxloc(ev),minloc(ev),sum(ev)
  print *, "  hv    ", maxval(hv),minval(hv),maxloc(hv),minloc(hv),sum(hv)
  print *, "  yvdeg ", maxval(yvdeg),minval(yvdeg),maxloc(yvdeg),minloc(yvdeg),sum(yvdeg)
  print *, "  ydeg  ", maxval(ydeg),minval(ydeg),maxloc(ydeg),minloc(ydeg),sum(ydeg)
  print *, "  cs    ", maxval(cs),minval(cs),maxloc(cs),minloc(cs),sum(cs)
  print *, "  csv   ", maxval(csv),minval(csv),maxloc(csv),minloc(csv),sum(csv)
  print *, "  ocs   ", maxval(ocs),minval(ocs),maxloc(ocs),minloc(ocs),sum(ocs)
  print *, "  ocsv  ", maxval(ocsv),minval(ocsv),maxloc(ocsv),minloc(ocsv),sum(ocsv)
  print *, "  dx    ", maxval(dx),minval(dx),maxloc(dx),minloc(dx),sum(dx)
  print *, "  dxv   ", maxval(dxv),minval(dxv),maxloc(dxv),minloc(dxv),sum(dxv)
  print *, "  odx   ", maxval(odx),minval(odx),maxloc(odx),minloc(odx),sum(odx)
  print *, "  odxv  ", maxval(odxv),minval(odxv),maxloc(odxv),minloc(odxv),sum(odxv)
  print *, "  dy    ", maxval(dy),minval(dy),maxloc(dy),minloc(dy),sum(dy)
  print *, "  dyv   ", maxval(dyv),minval(dyv),maxloc(dyv),minloc(dyv),sum(dyv)
  print *, "  ody   ", maxval(ody),minval(ody),maxloc(ody),minloc(ody),sum(ody)
  print *, "  odyv  ", maxval(odyv),minval(odyv),maxloc(odyv),minloc(odyv),sum(odyv)
#endif

  ! 7. read file EVP for AL, AR, ...
  ! EVP ELLIPTIC SOLVER DATA  !
  FNAME=TRIM(WORKDIR)//TRIM(IDIR)//'EVP'
  OPEN(FN_EVP,CONVERT='BIG_ENDIAN',FILE=TRIM(FNAME),FORM='unformatted')
  READ(FN_EVP) AL,AR,AB,AT,AC,OAC,RINV,RINV1,IE
  CLOSE(FN_EVP)

#ifdef _DEBUG_
  print *, "GLOBAL evp"
  print *, "  al   ", maxval(al),minval(al),maxloc(al),minloc(al),sum(al)
  print *, "  ar   ", maxval(ar),minval(ar),maxloc(ar),minloc(ar),sum(ar)
  print *, "  ab   ", maxval(ab),minval(ab),maxloc(ab),minloc(ab),sum(ab)
  print *, "  at   ", maxval(at),minval(at),maxloc(at),minloc(at),sum(at)
  print *, "  ac   ", maxval(ac),minval(ac),maxloc(ac),minloc(ac),sum(ac)
  print *, "  rinv ", maxval(rinv),minval(rinv),maxloc(rinv),minloc(rinv),sum(rinv)
  print *, "  rinv1", maxval(rinv1),minval(rinv1),maxloc(rinv1),minloc(rinv1),sum(rinv1)
  print *, "  ie   ", maxval(ie),minval(ie),maxloc(ie),minloc(ie),sum(ie)  
#endif
 
END SUBROUTINE READ_INPUT_FILES

END MODULE INPUT_GLOBAL
