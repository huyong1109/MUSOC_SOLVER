SUBROUTINE PREP_TIMCOM
#include "namelist.in_GLOBAL"
  USE INIT_VAR_GLOBAL
  USE GRID_VAR_GLOBAL, ONLY: RHO,P0,PX,PY,U1,U2,V1,V2,S1,S2,T1,T2,P,ULF,VLF,SLF,TLF,EV,HV,F,TANPHI,SUMIN,VAR_ALLOCATE_B
  USE GRID_VAR_GLOBAL, ONLY: VNOR,VSUD,UEST,UWST,VAR_ALLOCATE_B_PRE
  USE GRID_VAR_GLOBAL, ONLY: KB,IU0,IV0,IN,IU,IV,IW
  USE GRID_VAR_GLOBAL, ONLY: U,V,W,VAR_ALLOCATE_CGRID
  USE GRID_VAR_GLOBAL, ONLY: SSURF,TSURF,SSSP,SNSP,TSSP,TNSP,VAR_ALLOCATE_CLIMAT2
  USE GRID_VAR_GLOBAL, ONLY: A,XDEG,YV,YVDEG,YDEG,CS,OCS,DX,ODX,DY,ODY,CSV,OCSV,DXV,ODXV,DYV,ODYV,VAR_ALLOCATE_METRIC
  USE GRID_VAR_GLOBAL, ONLY: Y,VAR_ALLOCATE_METRIC_PRE
  USE GRID_VAR_GLOBAL, ONLY: Z,ODZ,ODZW,VAR_ALLOCATE_ZFS
  USE GRID_VAR_GLOBAL, ONLY: RINV,RINV1,DUM0,DUM1,DUM2,XX,H,X,S,AL,AB,AC,AR,AT,SRC,CL,CB,CC,CR,CT,IE,VAR_ALLOCATE_SEVP
  USE GRID_VAR_GLOBAL, ONLY: DEPTH
  USE GRID_VAR_GLOBAL, ONLY: VBK,HBK,VAR_ALLOCATE_WINDMX_MAIN
  
  IMPLICIT NONE

  INTEGER :: I,J,K,NSUM,NN,IT,ILO,IHI,JL,JH,JHP,JHM,NB,IM1,M,NM,MVI,INC,IL,IR,IR2,JG,NG,IG,IP1,IA,IB,N1,L
  CHARACTER(160)::TEMPDIR
  REAL::TINC,ZZ,TDIF,SMIN,SMAX,TMIN,TMP,AVG1,AVG2,AREAXY,TEMP,TEM,SVIN,SV,DDDY,&
        DMX,DMY,DZX,DZY,DZYB,DZYT,X1DEG,Y1DEG,SXYT
  
  INTEGER :: FNO_RUNDATA,FNO_EVP,FNO_TURBMIX

  CALL PREP

  CONTAINS

  !*******************************************************************
  SUBROUTINE PREP
  
    INTEGER :: I,J,K,N
    REAL    :: Y1DEG,PI_180,TMP,OMEGA2      !ZWM DXDEG,DYDEG removed, declared in MODULE XYFS_GLOBAL
  
    TEMPDIR=TRIM(WORKDIR)//'/INPUT_GLOBAL/RUNDATA'
    CALL OPEN_BIGENDIAN_OLD_BIN_FILE(TRIM(TEMPDIR),FNO_RUNDATA)
  
    IF (LEVP==1) THEN
      TEMPDIR=TRIM(WORKDIR)//'/INPUT_GLOBAL/EVP'
      CALL OPEN_BIGENDIAN_NEW_BIN_FILE(TRIM(TEMPDIR),FNO_EVP)
    END IF
  
    IE(NB0)=J1
  
    CALL INITFS 

    CLOSE(FNO_RUNDATA)
    CLOSE(FNO_EVP)
  
  
  END SUBROUTINE PREP

  !*******************************************************************
  ! INITIALIZATION *
  !*******************************************************************
  SUBROUTINE INITFS
  ! INITFS initializes all controlling arrays and derived scalars
  !(scalar control parameters are set in BLOCK DATA at start of this file)
  ! ----------------------------------------------------------------------
    USE COUPLING
    INTEGER :: STATUS(4),N

! Free surface variable
    REAL :: GAMMA
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !THERE ARE TWO MAJOR FUNCTIONS IN THIS SUBROUTINE !
    ! 1) SET "VERICAL" EDDY VISCOSITY & DIFFUSIVITY   !
    ! 2) COMPLETE PREPROCESS FOR ELLITIC SOLVER       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL VAR_ALLOCATE_WINDMX_MAIN
    PRINT *, "INITFS is CALLED!"
    MVI=0.25*DAODT !NUMBER OF TIME STEPS BETWEEN MOVIE FRAME STORAGE
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 1) SET "VERICAL" EDDY VISCOSITY & DIFFUSIVITY     !
    ! INTRODUCE TWO FLAGS:                              !
    ! a) FL_INI_VIS = 1: CALCULATE INITIAL VISCOSITY    !
    !    FL_INI_VIS = 0: NO NEED TO CALCULATE INI. VIS. !
    ! b) FL_BK_VIS  = 1: CALCULATE BACKGROUND VISCOSITY !
    !    FL_BK_VIS  = 0: NO NEED TO CALCULATE BK. VIS.  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
    ! CHECK 2D OR 3D CASE. ONLY 3D CASE NEEDS THIS PROCEDURE (NO NEED FOR ILE2D)
    IF (DIMEN==2) THEN  ! NO NEED FOR ILE2D
      PRINT *,"Case is 2D"
    ELSE IF (DIMEN==3) THEN  
      PRINT *,"Case is 3D"
      ! SET INITIAL VERTICAL DIFFUSIVITIES
      PRN=DE0_CONST/DM0_CONST
      IF (FL_INI_VIS == 1) THEN
        !DTRAC, NPB, & TAI
        DO K=1,K2
          ! SET VERTICAL DIFFUSIVITIES TO TL AND DE0/DM0*TL CM**2/SEC
          TMP=TL*DMZ0*ODZW(K+1) !TL IS GIVEN IN THE NANELIST
          TEMP=PRN*TMP
          DO J=1,J2
            DO I=1,I2
              EV(I,J,K)=TMP*IW(I+1,J+1,K+1)
              HV(I,J,K)=TEMP*IW(I+1,J+1,K+1)
            END DO
          END DO
        END DO
      END IF

      ! SET BACKGROUND VERICAL EDDY VISCOSITY & DIFFUSIVITY
      IF (FL_BK_VIS == 1) THEN
        !NPB, TAI & GLOBAL
          DO K=1,K2
            ! augment molecular viscosity & diffusivity by parameterized synoptic wind
            ! events, and by breaking (u,v,T,S, near surface only) and non-breaking
            ! (u,v only, at all depths) internal waves.
            ! synoptic wind forced mixing having 20m e-folding scale.
            !     TMP=EXP(-.0005*Z(2*K+1))
            ! increased synoptic wind forced mixing scale to 50m during year 30.
            TMP=EXP(-.0001*Z(2*K+1))
            ! bigger and deeper augmentation of momentum mixing due to pressure-xfers
            ! associated with internal waves that have no counterpart in T,S xfers.
            VBK(K)=20.+10.*TMP+50.*EXP(-.65E-10*(500.E2-Z(2*K+1))**2)
            HBK(K)=.001+10.*TMP !   HBK(K)=.001+TMP
          END DO
        
        PRINT *, 'SUM(VBK,HBK',SUM(VBK),SUM(HBK)
        PRINT *,CASE_NAME,DSCRIB,DT,TLZ,G,DM0_CONST,DE0_CONST,DMZ0,SUM(F),SUM(TANPHI),RZMX,TAU_CONST
        PRINT *,TAUN_CONST,FLTW,WRAF,DAODT,ODT
        PRINT *,PRN,ORZMX,OFLTW,SUM(Y),SUM(YV),SUM(YVDEG),SUM(YDEG),SUM(CS),SUM(CSV)
        PRINT *,SUM(OCS),SUM(OCSV),SUM(DX),SUM(DXV),SUM(ODX),SUM(ODXV),SUM(DY),SUM(DYV)
        PRINT *,SUM(ODY),SUM(ODYV),RUNS,ISAV
        PRINT *,LOPENN,LOPENS,LOPENE,LOPENW,MVI
    
        
        TEMPDIR=TRIM(WORKDIR)//'/INPUT_GLOBAL/TURBMIX'
        CALL OPEN_BIGENDIAN_NEW_BIN_FILE(TRIM(TEMPDIR),FNO_TURBMIX)
        WRITE(FNO_TURBMIX) VBK,HBK
        CLOSE(FNO_TURBMIX)
        
      END IF
    END IF

    WRITE(FNO_RUNDATA) F,TANPHI,EV,HV
    WRITE(FNO_RUNDATA) XDEG,YVDEG,YDEG,CS,CSV,OCS,OCSV,DX,DXV,ODX,ODXV,DY,DYV,ODY,ODYV
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 2) PREPROCESSOR FOR ELLITIC SOLVER !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (LEVP.NE.1) RETURN
    !IF (DIMEN == 2) NN = J2/NB0
    !IF (DIMEN == 3) NN = 6
    NN = ANINT(REAL(J2)/REAL(NB0))
    
    DO N=1,NB0
      IE(N)=MIN(1+NN*N,J1)
    END DO
    
    WRITE(*,251) IE
    251  FORMAT('IE'/(20I4))
    
    IF (DIMEN == 2) THEN
      IF (IE(NB0).GT.IE(NB1).AND.J1.LT.IE(NB1)+NN+1) GO TO 260
    ELSE IF (DIMEN == 3) THEN
      IF (IE(NB0).GT.IE(NB1).AND.J1.LT.IE(NB1)+NN+3) GO TO 260
    END IF
    
    WRITE(*,252) IE(NB1),J1
    252  FORMAT('IE(NB1),J1=',2I5,' not properly defined. '                          &
    'Fix above DO 250 loop and restart.')
    stop 1 
    
    260  IE(NB0)=J1
  
    IF (DIMEN==2) THEN
      DZX=ODX(1)**2/ODZ(1)
      DZY=ODY(1)**2/ODZ(1)
    
      DO J=1,J2
        DO I=1,I2
          AL(I,J)=DZX
          AR(I,J)=DZX
          AB(I,J)=DZY
          AT(I,J)=DZY
        END DO
      END DO
            
    ELSE IF (DIMEN==3) THEN
      ! first, do top layer (all water to regularize EVP domain)
      TMP=1./ODZ(1)
      DO J=1,J2
        DZX=ODX(J+1)**2*TMP
        DZYB=CSV(J)*OCS(J+1)*ODYV(J)*ODY(J+1)*TMP
        DZYT=CSV(J+1)*OCS(J+1)*ODYV(J+1)*ODY(J+1)*TMP
        DO I=1,I2
          AL(I,J)=DZX
          AR(I,J)=DZX
          AB(I,J)=DZYB
          AT(I,J)=DZYT
        END DO
      END DO
      
      ! remaining layers
      DO K=2,K1
        TMP=1./ODZ(K)
        DO J=1,J2
          DZX=ODX(J+1)**2*TMP
          DZYB=CSV(J)*OCS(J+1)*ODYV(J)*ODY(J+1)*TMP
          DZYT=CSV(J+1)*OCS(J+1)*ODYV(J+1)*ODY(J+1)*TMP
          DO I=1,I2
            AL(I,J)=AL(I,J)+DZX*IU(I,J+1,K)
            AR(I,J)=AR(I,J)+DZX*IU(I+1,J+1,K)
            AB(I,J)=AB(I,J)+DZYB*IV(I+1,J,K)
            AT(I,J)=AT(I,J)+DZYT*IV(I+1,J+1,K)
          END DO
        END DO
      END DO
    END IF
  
    IF (NBIR==1) THEN
      DO J=1,J2
        AL(1,J)=0.
        AR(I2,J)=0.
      END DO
    END IF
  
    DO I=1,I2
      AB(I,1)=0.
      AT(I,J2)=0.
    END DO
  
    !pin down pressure
    IF (NBIR==1) THEN
      AB(I2/2,1)=AB(I2/2,2)
    END IF
  
! Yu-Chiao Free Surface
 GAMMA=0.
#ifdef FLAG_FS
! PRINT*,'hahaha we have free surface now!!!',DT,G,1./(.5*DT**2*G)
! PAUSE 
 GAMMA=1./(.5*DT**2*G)
#endif

    DO J=1,J2
      DO I=1,I2
        AC(I,J)=-AL(I,J)-AR(I,J)-AB(I,J)-AT(I,J)-GAMMA
      END DO
    END DO

    !!!! test ideal case  hyedit
    AL = -0.25
    AB = -0.25
    AC =  1.0
    AR = -0.25
    AT = -0.25

    ! Neumann  boundary on NS
    DO I=1,I2
        J = J2 ! north 
        AC(i,j)=AC(i,j)+AT(i,j)
        AT(i,j) = 0.
        J = 1  ! south
        AC(i,j)=AC(i,j)+AB(i,j)
        AB(i,j) = 0.
    END DO



    WRITE(*,295)
    295  FORMAT('ENTER EVP PREPROCESSOR')
  
    IF (FL_EVP_STP == 0) THEN !LINE #316
      !CALL PRE(AL,AB,AC,AR,AT,RINV,RINV1,H,IE,I2,I0,I2,NB0)
      !CALL PRE_PERIOD(AL,AB,AC,AR,AT,RINV,RINV1,H,IE,I2,I0,I2,NB0)
      CALL PRE1(AL,AB,AC,AR,AT,RINV,RINV1,H,IE,I0,I2,NB0)
    ELSE 
      !CASE_NAME=='GLOBAL'
      ! we assume uniform BIR strip widths
      I=2*I2
      IF (IBIR*NBIR.EQ.I) GO TO 300
      WRITE(*,298) I,NBIR,I2
      298  FORMAT('2*I2=',I4,                                                &
      'which must be evenly divisible by user-specified NBIR.'/       &
      'NOTE: it is possible to avoid this limitation,'/               &
      'but not considered worthwhile.'/                               &
      'Assigned NBIR value (',I2,') must be changed. STOP.')
  
      300  INC=IBIR/2
      
      DO N=1,NBIR !LINE #332
        WRITE(*,340) N
        340  FORMAT('processing BIR strip #:',I3)
        IL=(N-1)*INC+1
        
        IF (NBIR==1 .OR. N/=NBIR) THEN !LINE #337
          CALL PRE(AL(IL,1),AB(IL,1),AC(IL,1),AR(IL,1),AT(IL,1),RINV(1,1,N),RINV1(1,1,N),H,IE,I2,IBIR+2,IBIR,NB0)
        ELSE
          ! last BIR strip
          IR=IL+INC-1
          IR2=IL+IBIR-1
          
          DO J=1,J2
            DO I=IL,IR
              CL(I-IL+1,J)=AL(I,J)
              CR(I-IL+1,J)=AR(I,J)
              CB(I-IL+1,J)=AB(I,J)
              CT(I-IL+1,J)=AT(I,J)
              CC(I-IL+1,J)=AC(I,J)
            END DO
          
            DO I=IR+1,IR2
              CL(I-IL+1,J)=AL(I-IR,J)
              CR(I-IL+1,J)=AR(I-IR,J)
              CB(I-IL+1,J)=AB(I-IR,J)
              CT(I-IL+1,J)=AT(I-IR,J)
              CC(I-IL+1,J)=AC(I-IR,J)
            END DO
          END DO
          CALL PRE(CL,CB,CC,CR,CT,RINV(1,1,N),RINV1(1,1,N),H,IE,IBIR,IBIR+2,IBIR,NB0)      
        END IF !IF STARTED AT LINE #337, IF (NBIR==1 .OR. N/=NBIR)
      END DO  !DO STARTED AT LINE #332    
    END IF !IF STARTED AT LINE #316, IF (FL_EVP_STP == 0)
    
    WRITE(FNO_EVP) AL,AR,AB,AT,AC,RINV,RINV1,IE
  
  END SUBROUTINE INITFS
! ******************************
! ELLIPTIC SOLVER PREPROCESSOR *
! ******************************
! ----------------------------------------------------------------------
      SUBROUTINE PRE1(AX,AY,BB,CX,CY,RINV,RINV1,H,IE,I0,I2,NBLK)
! ----------------------------------------------------------------------
      REAL*8 RINV,RINV1,H
      REAL:: AX,AY,BB,CX,CY
      INTEGER:: I0, I2,IE,NBLK,N
      DIMENSION AX(I2,1),AY(I2,1),BB(I2,1),CX(I2,1),CY(I2,1),&
        RINV(I2,I2,1),RINV1(I2,I2,1),H(I0,1),IE(1)
! Periodic
      I1=I0-1
      JL=1
      NB=0
 100  NB=NB+1
      WRITE(*,111) NB,NBLK
 111  FORMAT ('Processing block #',I3,' out of',I3,' total evp solver blocks')
      JH=IE(NB)
      JHP=JH+1
      JHM=JH-2
      JG=JL+1
      DO 250 NG=1,I2
      IG=NG+1
      DO 210 J=JL,JHP
      DO 210 I=1,I0
 210  H(I,J)=0.
      H(IG,JG)=1.
      H(1,JG)=H(I1,JG)
      H(I0,JG)=H(2,JG)
      IF (NB.EQ.1) GO TO 220
      DO 218 N=1,I2
 218  H(N+1,JL)=RINV1(NG,N,NB-1)*CY(IG-1,JG-2)
 220  DO 225 J=JL,JHM
      DO 225 I=1,I2
! Periodic
      H(1,J+1)=H(I1,J+1)
      H(I0,J+1)=H(2,J+1)
 225  H(I+1,J+2)=-(AX(I,J)*H(I,J+1)+AY(I,J)*H(I+1,J)+BB(I,J)*H(I+1,J+1)+ &
          CX(I,J)*H(I+2,J+1))/CY(I,J)
      J=JH-1
      DO 230 I=1,I2
      H(1,J+1)=H(I1,J+1)
      H(I0,J+1)=H(2,J+1)
 230  RINV(NG,I,NB)=AX(I,J)*H(I,J+1)+AY(I,J)*H(I+1,J)+BB(I,J)* &
          H(I+1,J+1)+CX(I,J)*H(I+2,J+1)
      IF (NB.EQ.NBLK) GO TO 250
      J=IE(NB)
      DO 240 N=1,I2
 240  RINV(NG,N,NBLK)=H(N+1,J)
! Periodic
      H(1,JG)=0.
      H(I0,JG)=0.
 250  H(IG,JG)=0.
      CALL MATINV(RINV(1,1,NB),I2)
      IF (NB.EQ.NBLK) RETURN
      DO 260 I=1,I2
      DO 260 J=1,I2
      RINV1(I,J,NB)=0.
      DO 260 K=1,I2
 260  RINV1(I,J,NB)=RINV1(I,J,NB)-RINV(I,K,NB)*RINV(K,J,NBLK)
      JL=JH
      GO TO 100
      END subroutine

  !*******************************************************************
  ! ELLIPTIC SOLVER PREPROCESSOR *
  !*******************************************************************
  SUBROUTINE PRE_PERIOD(AX,AY,BB,CX,CY,RINV,RINV1,H,IE,I2,M0,M2,NBLK)
    ! ----------------------------------------------------------------------
    ! M2 is the working BIR strip width (I2 is EVENLY divisible by M2)
    REAL(8) RINV,RINV1,H
    REAL:: AX,AY,BB,CX,CY
    INTEGER:: I2,M2,M0,IE,NBLK,N
    DIMENSION AX(I2,*),AY(I2,*),BB(I2,*),CX(I2,*),CY(I2,*),   &
    RINV(M2,M2,*),RINV1(M2,M2,*),H(M0,*),IE(*)
    I1=I0-1 ! periodic
    JL=1
    NB=0
    100  NB=NB+1
    WRITE(*,111) NB,NBLK
    111  FORMAT ('Processing block #',I3,' out of',I3,' total evp solver blocks')
    JH=IE(NB)
    JHP=JH+1
    JHM=JH-2
    JG=JL+1
    DO 250 NG=1,M2
    IG=NG+1
    DO 210 J=JL,JHP
    DO 210 I=1,M0
    210  H(I,J)=0.
    H(IG,JG)=1.
    H(1,JG)=H(I1,JG)
    H(I0,JG)=H(2,JG)
    IF (NB.EQ.1) GO TO 220
    DO N=1,M2
    H(N+1,JL)=RINV1(NG,N,NB-1)*CY(IG-1,JG-2)
    END DO
    220  DO 225 J=JL,JHM
    DO 225 I=1,M2
      H(1,J+1)=H(I1,J+1)
      H(I0,J+1)=H(2,J+1)
    225  H(I+1,J+2)=-(AX(I,J)*H(I,J+1)+AY(I,J)*H(I+1,J)+BB(I,J)*H(I+1,J+1)+   &
    CX(I,J)*H(I+2,J+1))/CY(I,J)
    J=JH-1
    DO 230 I=1,M2
    H(1,J+1)=H(I1,J+1)
    H(I0,J+1)=H(2,J+1)
    230  RINV(NG,I,NB)=AX(I,J)*H(I,J+1)+AY(I,J)*H(I+1,J)+BB(I,J)*             &
    H(I+1,J+1)+CX(I,J)*H(I+2,J+1)
    
    IF (NB.EQ.NBLK) GO TO 250
    J=IE(NB)
    DO N=1,M2
      RINV(NG,N,NBLK)=H(N+1,J)
    END DO
      H(1,JG)=0.
      H(I0,JG)=0.
    250  H(IG,JG)=0.
    CALL MATINV(RINV(1,1,NB),M2)
    
    IF (NB.EQ.NBLK) RETURN
    DO I=1,M2
      DO J=1,M2
        RINV1(I,J,NB)=0.
        DO K=1,M2
          RINV1(I,J,NB)=RINV1(I,J,NB)-RINV(I,K,NB)*RINV(K,J,NBLK)
        END DO
      END DO
    END DO
    JL=JH
    GO TO 100
  
  END SUBROUTINE PRE_PERIOD
  SUBROUTINE PRE(AX,AY,BB,CX,CY,RINV,RINV1,H,IE,I2,M0,M2,NBLK)
    ! ----------------------------------------------------------------------
    ! M2 is the working BIR strip width (I2 is EVENLY divisible by M2)
    REAL(8) RINV,RINV1,H
    REAL:: AX,AY,BB,CX,CY
    INTEGER:: I2,M2,M0,IE,NBLK,N
    DIMENSION AX(I2,*),AY(I2,*),BB(I2,*),CX(I2,*),CY(I2,*),   &
    RINV(M2,M2,*),RINV1(M2,M2,*),H(M0,*),IE(*)
    JL=1
    NB=0
    100  NB=NB+1
    WRITE(*,111) NB,NBLK
    111  FORMAT ('Processing block #',I3,' out of',I3,' total evp solver blocks')
    JH=IE(NB)
    JHP=JH+1
    JHM=JH-2
    JG=JL+1
    DO 250 NG=1,M2
    IG=NG+1
    DO 210 J=JL,JHP
    DO 210 I=1,M0
    210  H(I,J)=0.
    H(IG,JG)=1.
    IF (NB.EQ.1) GO TO 220
    DO N=1,M2
    H(N+1,JL)=RINV1(NG,N,NB-1)*CY(IG-1,JG-2)
    END DO
    220  DO 225 J=JL,JHM
    DO 225 I=1,M2
    225  H(I+1,J+2)=-(AX(I,J)*H(I,J+1)+AY(I,J)*H(I+1,J)+BB(I,J)*H(I+1,J+1)+   &
    CX(I,J)*H(I+2,J+1))/CY(I,J)
    J=JH-1
    DO 230 I=1,M2
    230  RINV(NG,I,NB)=AX(I,J)*H(I,J+1)+AY(I,J)*H(I+1,J)+BB(I,J)*             &
    H(I+1,J+1)+CX(I,J)*H(I+2,J+1)
    
    IF (NB.EQ.NBLK) GO TO 250
    J=IE(NB)
    DO N=1,M2
      RINV(NG,N,NBLK)=H(N+1,J)
    END DO
    250  H(IG,JG)=0.
    CALL MATINV(RINV(1,1,NB),M2)
    
    IF (NB.EQ.NBLK) RETURN
    DO I=1,M2
      DO J=1,M2
        RINV1(I,J,NB)=0.
        DO K=1,M2
          RINV1(I,J,NB)=RINV1(I,J,NB)-RINV(I,K,NB)*RINV(K,J,NBLK)
        END DO
      END DO
    END DO
    JL=JH
    GO TO 100
  
  END SUBROUTINE PRE
  !*******************************************************************
  
  SUBROUTINE MATINV(B,N)
    IMPLICIT REAL*8 (A-H,O-Z)
    INTEGER:: N
    DIMENSION B(N,*),B1(I2),B2(I2)
    
    N1=N-1
    DO I=1,N1
      B1(1)=1./B(I,I)
      B(I,I)=1.0
      DO J=1,N
        B(I,J)=B(I,J)*B1(1)
      END DO
      IP1=I+1
      DO IA=IP1,N
        B1(IA)=B(IA,I)
      END DO
      DO IA=IP1,N
        B(IA,I)=0.
      END DO
      DO J=1,N
        B2(J)=B(I,J)
      END DO
      DO IA=IP1,N
        DO J=1,N
          B(IA,J)=B(IA,J)-B1(IA)*B2(J)
        END DO
      END DO
    END DO
    B1(1)=1./B(N,N)
    B(N,N)=1.
    DO J=1,N
      B(N,J)=B(N,J)*B1(1)
    END DO
    DO I=2,N
      DO IB=1,I
        B1(IB)=B(IB,I)
      END DO
      IM1=I-1
      DO IB=1,IM1
        B(IB,I)=0.
      END DO
      DO J=1,N
        B2(J)=B(I,J)
      END DO
      IM1=I-1
      DO IB=1,IM1
        DO J=1,N
          B(IB,J)=B(IB,J)-B1(IB)*B2(J)
        END DO
      END DO
    END DO
  
  END SUBROUTINE MATINV
  !*******************************************************************
  
  !*******************************************************************
  FUNCTION r(t,s,p)
    REAL:: r,t,s,p
    REAL:: p02,p01
    p02=1747.4508988+t*(11.51588-0.046331033*t)              &
    -s*(3.85429655+0.01353985*t)
    p01=p/10d0+5884.81703666+t*(39.803732+t*(-0.3191477      &
    !   p01=       5884.81703666+t*(39.803732+t*(-0.3191477
    +t*0.0004291133))+2.6126277*s
    r=p01/(p02+0.7028423*p01)
    RETURN
  END FUNCTION r
  !*******************************************************************

END SUBROUTINE PREP_TIMCOM

