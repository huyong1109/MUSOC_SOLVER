SUBROUTINE METGEN
  USE TIMCOM_GENERAL
  USE OCN_PARA            !  I0,J0,K0
  USE OCN_PARA_GLOBAL      !   
  USE XYFS_GLOBAL          !  VARS CONTROLING XY COORDINATE: NAMELIST /GRID_XY/
  USE ZCO_GLOBAL           ! VAR CONTROLING Z COORDINATE
  USE TPS_NAMELIST_GLOBAL      !  NAMELIST INVOLVING ALL CONTROLING VARS
  USE GRID_VAR_GLOBAL, ONLY: K0,K1,Z,ODZ,ODZW,VAR_ALLOCATE_ZFS 
  USE GRID_VAR_GLOBAL, ONLY: RHO,P0,PX,PY,U1,U2,V1,V2,S1,S2,T1,T2,P,ULF,VLF,SLF,TLF,EV,HV,F,TANPHI,SUMIN,VAR_ALLOCATE_B
  USE GRID_VAR_GLOBAL, ONLY: VNOR,VSUD,UEST,UWST,VAR_ALLOCATE_B_PRE
  USE GRID_VAR_GLOBAL, ONLY: A,XDEG,YV,YVDEG,YDEG,CS,OCS,DX,ODX,DY,ODY,CSV,OCSV,DXV,ODXV,DYV,ODYV,VAR_ALLOCATE_METRIC
  USE GRID_VAR_GLOBAL, ONLY: Y,VAR_ALLOCATE_METRIC_PRE
  
  IMPLICIT NONE
  
  !=======================
  !--DECLARING TEMP VARS
  !=======================
  
  INTEGER ::IOST
  INTEGER ::J,K
  REAL	::DXRAD,DYRAD,Y1DEG,PI_180,TMP,TEMP,DEG_LAT,OMEGA2,PHI
  INTEGER	::I,N,IMAX,IMIN,G_SIZE(3)   !grid size
  CHARACTER(160)::TEMPDIR
  INTEGER :: OPENSTATUS
  INTEGER ::FNO_TPS,FNO_KBVIEW,FNO_YZGRID

  ! initialize the default namelist variables 
  TLX=0.;TLY=0.

! Read the namelist.tps_GLOBAL
  CALL OPEN_TPS_NAMELIST_FILE(FNO_TPS)
  REWIND(FNO_TPS)
  READ(FNO_TPS,NML=CASE_INFO)
  REWIND(FNO_TPS)
  READ(FNO_TPS,NML=WORK_DIR)
  REWIND(FNO_TPS)
  READ(FNO_TPS,NML=DIMENSIONS)
  REWIND(FNO_TPS)
  READ(FNO_TPS,NML=GLOBAL)
  REWIND(FNO_TPS)
  READ(FNO_TPS,NML=GRID_XY)
  REWIND(FNO_TPS)
  READ(FNO_TPS,NML=TEMPORAL)
  REWIND(FNO_TPS)
  READ(FNO_TPS,NML=PHYSICS)
  REWIND(FNO_TPS)
  READ(FNO_TPS,NML=SOLVER)
  REWIND(FNO_TPS)
  READ(FNO_TPS,NML=Z_OPTION)
  SELECT CASE(ZOP)
  CASE(1)
    REWIND(FNO_TPS)
    READ(FNO_TPS,NML=Z_LINEAR)
  END SELECT
  CALL CLOSE_TPS_NAMELIST_FILE

  I0=I0_GLOBAL
  J0=J0_GLOBAL
  K0=K0_GLOBAL

  CALL VAR_SPECIFY_OCN_PARA
  CALL VAR_ALLOCATE_B
  CALL VAR_ALLOCATE_B_PRE
  CALL VAR_ALLOCATE_ZFS
  CALL VAR_ALLOCATE_METRIC
  CALL VAR_ALLOCATE_METRIC_PRE

  !=======================
  !--CREATE Z  COORDINATE
  !=======================
  SELECT CASE(ZOP)
  CASE(1)
    WRITE(*,*), 'Z_LINEAR_C=',C,'D=',D
    CALL SUB_Z_LINEAR(TLZ,  C,D,Z, K01)
  CASE(2)
    CALL SUB_Z_STRETCH(TLZ,  ZTOP,Z,K01)
  CASE(3)
    CALL SUB_Z_USRDEF(Z,K01)
  CASE DEFAULT
    PRINT *, 'WRONG Z COORDINATE OPTIONS'
    STOP
  END SELECT

  !==================
  ! vertical metrics
  !==================
  DO K=1,K1
    ODZ(K)=1./(Z(2*K+1)-Z(2*K-1))	!C
  END DO
  DO K=2,K1
    ODZW(K)=1./(Z(2*K)-Z(2*K-2))	!C
  END DO

  !======================
  ! output data : KBVIEW
  !======================
  TEMPDIR=TRIM(WORKDIR)//'/TEMP_GLOBAL/KBVIEW'
  CALL OPEN_BIGENDIAN_NEW_TXT_FILE(TRIM(TEMPDIR),FNO_KBVIEW)
  WRITE(*,9) (.01*Z(K),K=1,K0+K1,2)
  WRITE(FNO_KBVIEW,9) (.01*Z(K),K=1,K0+K1,2)
  WRITE(FNO_KBVIEW,10) (.01*Z(K),K=2,K0+K1,2)
9     FORMAT('z-levels(m)'/(10F8.2))
10    FORMAT('z-levels(m)-midpoint'/(10F8.2))

  CLOSE(FNO_KBVIEW)

!==============
! read-in PARA
!==============
!@  REWIND(FNO_TPS)
!@  IF (DIMEN==3) THEN	! check 2D or 3D case	ZWM
!@    PRINT *,"METGEN Case is 3D"
!@    READ(FNO_TPS,NML=PARA_3D)	!ZWM	note that nml must be in order. if para_2d then para_3d, then there'll be error
!@  ELSE IF (DIMEN==2) THEN	! check 2D or 3D case	ZWM
!@    PRINT *,"METGEN Case is 2D."
!@    READ(FNO_TPS,NML=PARA_2D)	!ZWM
!@  ELSE
!@    PRINT *,"Default case, nothing done."
!@    PAUSE
!@  END IF

  !================
  !setup TLX & TLY
  !================
  IF ( abs(TLY-0.) < 1.E-5 ) TLY=TLX


  !=======================
  !--CREATE XY COORDINATE
  !=======================
  PI_180=3.141592654/180.		!P but simply an unchanged unconstant
  OMEGA2=3.141592654/21600.

  COORD_SYSTEM: IF (TRIM(COORD_SYS)=="SPHERE") THEN !for global & multi case

    DXDEG=DXMNUT/60.
    XDEG(1)=X0DEG-0.5*DXDEG
    DO I=2,I0
      XDEG(I)=XDEG(I-1)+DXDEG
    END DO
    YVDEG(1)=Y0DEG
    YDEG(1)=Y0DEG
    YV(1)=0.
    DY(1)=0.
    ODY(1)=0.
  
    SELECT CASE(DX0STD)
    CASE(0)
      ! southernmost DX
        DX0=R0*DXDEG*PI_180*COS(YVDEG(1)*PI_180)	  
    CASE(1)
      ! Equatorial DX
        DX0=R0*DXDEG*PI_180
    CASE(2)
    CASE DEFAULT
        PRINT *, "WRONG DX0STD VALUE"
        PRINT *, "DX0STD must be"
        PRINT *, "  0 : southernmost DX"
        PRINT *, "  1 : equatorial DX"
        PRINT *, "  2 : user-defined DX"
    END SELECT

    DO J=2,J0
    ! initial guess
      YVDEG(J)=YVDEG(J-1)+DXDEG*DYDX*COS(YVDEG(J-1)*PI_180)
      DX(J)=DX0*COS(.5*(YVDEG(J-1)+YVDEG(J))*PI_180)
      N=0
      DO WHILE(.TRUE.)
        N=N+1
        TMP=YVDEG(J)
        YVDEG(J)=YVDEG(J-1)+DYDX*DX(J)/(R0*PI_180)
        YDEG(J)=.5*(YVDEG(J-1)+YVDEG(J))
        CS(J)=COS(YDEG(J)*PI_180)
        DXRAD=DXDEG*PI_180*CS(J)
        DX(J)=DXRAD*R0
        DYRAD=(YVDEG(J)-YVDEG(J-1))*PI_180
        TEMP=DYRAD/DXRAD
        IF (ABS(TMP-YVDEG(J)).LT..0001) EXIT
      ENDDO
      DYDEG=YVDEG(J)-YVDEG(J-1)
      DY(J)=DYDEG*R0*PI_180
      DYV(J-1)=(YDEG(J)-YDEG(J-1))*R0*PI_180
      ODY(J)=1./DY(J)
      ODYV(J-1)=1./DYV(J-1)
      YV(J)=YV(J-1)+DY(J)
      Y(J)=.5*(YV(J)+YV(J-1))
      OCS(J)=1./CS(J)
      ODX(J)=1./DX(J)
      CSV(J-1)=COS(PI_180*YVDEG(J-1))
      OCSV(J-1)=1./CSV(J-1)
      DXV(J-1)=DX0*CSV(J-1)
      ODXV(J-1)=1./DXV(J-1)
    END DO

    Y1DEG=YVDEG(J1)
    DO J=2,J1
    ! spherical curvature parameter
      TANPHI(J-1)=TAN(PI_180*YDEG(J))/R0
      F(J-1)=OMEGA2*SIN(PI_180*YDEG(J))
    END DO

  ELSE IF (TRIM(COORD_SYS)=="CART") THEN COORD_SYSTEM  ! for DTRAC case
    TLY=400.0*100000  
    YDEG(1)=0.
    DO J=2,J0
      DX(J)=TLX/I2	!ALL C
      DY(J)=TLY/J2
      DYV(J-1)=DX(J)
      ODY(J)=1./DY(J)
      ODYV(J-1)=1./DYV(J-1)
      YV(J)=YV(J-1)+DY(J)
      Y(J)=.5*(YV(J)+YV(J-1))
      OCS(J)=1.
      ODX(J)=1./DX(J)
      CSV(J-1)=1.
      OCSV(J-1)=1.
      DXV(J-1)=DX(J)
      ODXV(J-1)=1./DXV(J-1)
      YDEG(J)=DY(J)/2.+YDEG(J-1) 
    END DO
    ODX(1)=ODX(2)	!C
    DO I=1,I0
      XDEG(I)=DX(2)*REAL(I-1)  
    END DO
    DO J=2,J1
      F(J-1)=OMEGA2*SIN(PI_180*Y0DEG)
    END DO

  ELSE IF (TRIM(COORD_SYS)=="OTHER") THEN COORD_SYSTEM  ! for ILE2D

  ! DX is in cm at this point
    DX=TLX/I2		!C
    DY=TLY/J2		!C
    ODX=1./DX
    ODY=1./DY
    DO I=1,I0
      XDEG(I)=DX(1)*REAL(I-1)   !lat/lon is simply index
    END DO
    DO J=1,J0
      YDEG(J)=DY(1)*REAL(J-1)
    END DO

  END IF COORD_SYSTEM

  !output YZGRID data
  TEMPDIR=TRIM(WORKDIR)//'/TEMP_GLOBAL/YZGRID'
  CALL OPEN_BIGENDIAN_NEW_BIN_FILE(TRIM(TEMPDIR),FNO_YZGRID)
  WRITE(FNO_YZGRID) F,TLZ,Z,ODZ,ODZW,DXMNUT,Y,YV,YVDEG,YDEG,CS,CSV,OCS,OCSV,DX,DXV,ODX,ODXV,DY,DYV,ODY,ODYV
  CLOSE(FNO_YZGRID)
  
  CONTAINS

  SUBROUTINE SUB_Z_LINEAR(DEPTH,C,D,Z,ZK01)
    IMPLICIT NONE
    REAL	::TMP,DZZ,A,B,ZZ,C,D,DEPTH
    INTEGER	::K,ZK01    !  ZK01=K01
    INTENT(IN)::DEPTH,C,D,ZK01
    REAL,INTENT(OUT),DIMENSION(ZK01)::Z
    
    PRINT *,'DEPTH:',DEPTH,'C',C,'D',D,'ZK01',ZK01
    tmp=.01*DEPTH
    dzz=tmp/(ZK01-1)
    A=tmp*(1.-D)/(1.-exp(C*tmp))
    B=-A
    zz=0.
    
    DO k=1,ZK01
      Z(k)=100.*(A+B*exp(C*zz)+D*zz)
!      PRINT *, Z(K)
      ZZ=ZZ+DZZ
    END DO
  END SUBROUTINE SUB_Z_LINEAR

  !----------------------------------------------------------------
  !   creating stretched z coordinate
  !   CALL   SUB_Z_STRETCH(DLZ,  ZTOP,Z,K01)
  SUBROUTINE SUB_Z_STRETCH(DEPTH,ZTOP,Z,ZK01)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ZK01
    REAL,INTENT(IN)::ZTOP,DEPTH
    REAL,INTENT(OUT),DIMENSION(ZK01)::Z
    INTEGER	::N,K,K1,K01
    REAL	::ZB0,GX0,GX,DGX,DZ,ZB,ERROR,ZBOT,SLOPE
    
    ZBOT=DEPTH
    Z(1)=0.
    ! SET TOP PRESSURE LOCATION (NEAR MIDDLE OF TOP LAYER)
    Z(2)=ZTOP
    ! INITIAL ZB
    ZB0=(ZK01-1)*ZTOP
    GX0=1.
    GX=1.1
    DGX=GX-GX0
    N=0
    DO WHILE(.TRUE.)
      N=N+1
      DZ=ZTOP
      DO K=3,ZK01
        DZ=GX*DZ
        Z(K)=Z(K-1)+DZ
      END DO
      ZB=Z(ZK01)
      ERROR=ZB-ZBOT
      SLOPE=(ZB-ZB0)/DGX
      DGX=-ERROR/SLOPE
      IF (ABS(ERROR).LT..00002*ZBOT) EXIT
      GX0=GX
      GX=GX0+DGX
      ZB0=ZB
    END DO
    
    WRITE(*,51) GX,ZB,ZBOT
51  FORMAT(//'NEWTON-RAPHSON CONVERGED GRID EXPANSION FACTOR = ',F8.5/&
             'EXACT DEPTH = ',1PE12.5,', MODEL BOTTOM DEPTH = ',1PE12.5)
  END SUBROUTINE SUB_Z_STRETCH

  !----------------------------------------------------------------
  !     read in user defined z coordinate
  !     CALL SUBROUTINE SUB_Z_USRDEF
  SUBROUTINE SUB_Z_USRDEF(Z,ZK01)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ZK01
    REAL,INTENT(INOUT)::Z(ZK01)
    INTEGER	::FNO,K
    
    CALL OPEN_OLD_TXT_FILE("Z_CUSTOM_GLOBAL.TXT",FNO)
    DO K=1,ZK01
      READ(FNO,*) Z(K)
    END DO
    CLOSE(FNO)
  END SUBROUTINE SUB_Z_USRDEF

END SUBROUTINE METGEN
