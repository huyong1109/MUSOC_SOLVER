MODULE COUPLING
  !coupling related subroutine
  USE OCN_PARA_GLOBAL
  IMPLICIT NONE
  SAVE
  REAL,ALLOCATABLE :: TSP(:,:,:),SSP(:,:,:) !NEW FOR MULTI-DOMAIN

CONTAINS
  
  !==============================================================
  SUBROUTINE DE_INI_VAR_COUPLE
    INTEGER :: STATUS(2)
    if (allocated(TSP)) then  
      deallocate (SSP,STAT=status(1)); deallocate (TSP,STAT=status(2))
      if (status(1)/=0 .or. status(2)/=0) STOP "Cannot deallocate memory"
    end if
  END SUBROUTINE DE_INI_VAR_COUPLE
  
  !==============================================================
  SUBROUTINE VAR_ALLOCATE_CLIMAT_COUPLE
    INTEGER :: STATUS(4)
    allocate (TSP(J0,K1,12),STAT=status(1)); allocate (SSP(J0,K1,12),STAT=status(2))
    if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot allocate memory"
    SSP=0.;TSP=0.
  END SUBROUTINE VAR_ALLOCATE_CLIMAT_COUPLE
  
  !==============================================================
  SUBROUTINE COUPLE_BOUND_COND(POSITION)
    !BOUNDARY CONDITIONS FOR COUPLING
    USE INIT_VAR_GLOBAL
    USE GRID_VAR_GLOBAL, ONLY: RHO,P0,PX,PY,U1,U2,V1,V2,S1,S2,T1,T2,P,ULF,VLF,SLF,TLF,EV,HV,F,TANPHI,SUMIN
    USE GRID_VAR_GLOBAL, ONLY: VNOR,VSUD,UEST,UWST
    USE GRID_VAR_GLOBAL, ONLY: KB,IU0,IV0,IN,IU,IV,IW
    USE GRID_VAR_GLOBAL, ONLY: U,V,W
    USE GRID_VAR_GLOBAL, ONLY: SSURF,TSURF,SSSP,SNSP,TSSP,TNSP
    USE CONTROL_GLOBAL
    USE XYFS_GLOBAL  ! R0,DYDX,Y0DEG,WESTDEG
    USE GRID_VAR_GLOBAL, ONLY: A,XDEG,YV,YVDEG,YDEG,CS,OCS,DX,ODX,DY,ODY,CSV,OCSV,DXV,ODXV,DYV,ODYV
    USE GRID_VAR_GLOBAL, ONLY: Y
    USE GRID_VAR_GLOBAL, ONLY: Z,ODZ,ODZW
    USE SCA_GLOBAL    ! Derived scalars
    USE GRID_VAR_GLOBAL, ONLY: RINV,RINV1,DUM0,DUM1,DUM2,XX,H,X,S,AL,AB,AC,AR,AT,SRC,CL,CB,CC,CR,CT,IE
    USE GRID_VAR_GLOBAL, ONLY: DEPTH
    USE TIMCOM_GENERAL
    USE OPFREQ	!ZWM
    
    IMPLICIT NONE
    REAL::TMIN,TMAXX,TMP1,TMP2,SUMINN,TMP
    INTEGER, INTENT(IN) :: POSITION !select position when BC is called
    INTEGER :: I,J,K,JLO,NM,JJ,N,ILO
    
    IF (TRIM(CASE_NAME)=='NPB') THEN
      IF (POSITION==1) THEN
        !hardwired: artificial closed boundary and shelf west of Kamchitka
        DO I=2,K1
          TMP=.01*Z(2*I+1)
          DO J=J0-59,J1
            DEPTH(I,J)=MIN(DEPTH(I,J),TMP)
          END DO
        END DO
        !MUST have depth=0 at nothern and southern boundary of overlap zone!!!
        !Else, one must open up lateral boundary NE, SE corners of TAI grid.
        DO J=120,126
          DO I=2,10
            DEPTH(I,J)=(1.-EXP(-.5*(I-2.)))*DEPTH(I,J)
          END DO
        END DO
        JLO=(488-3)/2+128
        DO J=JLO,JLO+6
          DO I=2,10
            DEPTH(I,J)=(1.-EXP(-.5*(I-2.)))*DEPTH(I,J)
          END DO
        END DO
      ELSE IF (POSITION==2) THEN
        !hardwired: deepen shelf to open connection to Sea of Okhotsk
        DO 40 I=2,K1
        K=KB(I,J0-59)
        J=J0-59
        35   J=J-1
        K=K+1
        KB(I,J)=MIN(KB(I,J),K)
        DEPTH(I,J)=Z(2*KB(I,J)+1)
        IF (K.EQ.K1) GO TO 40
        GO TO 35
        40   CONTINUE
      END IF
    
    ELSE IF (TRIM(CASE_NAME)=='TAI') THEN
      IF (POSITION==1) THEN    
        !extend NPB artificial shelf along northern boundary
        DO J=J0-30,J0-1
          TMP=.01*Z(2*KB(I1,J)+1)
          DO I=345,I1
            DEPTH(I,J)=MIN(DEPTH(I,J),TMP)
          END DO
        END DO
        DO J=2,J0-30
          TMP=.01*Z(2*KB(I1,J)+1)
          DO I=I3,I1
            DEPTH(I,J)=TMP
          END DO
        END DO
        !Reset KB for summation
        DO J=2,J1
          KB(I1,J)=0
        END DO
      END IF

   ELSE IF (TRIM(CASE_NAME)=='NAB') THEN
      IF (POSITION==1) THEN    
! GOM grid overlap zones must be closed latitudinally!!!
      DO J=35,39
        DO I=2,10
          DEPTH(I,J)=(1.-EXP(-.5*(I-2.)))*DEPTH(I,J)
        END DO
      END DO
! Eastern shelf is in IBE domain.Open using FLAG_LBS_N
! western shelf
  ! luc commented this out due to artificial current in the NE corner of IBE, 07/01/2014
  !    JJ=MIN(K1/2,J2/8)
  !    n=int(K1/JJ)
  !    DO I=1,JJ
  !      TMP=.01*Z(2*n*I+1)
  !      DO J=267,J1
  !        DEPTH(I,J) =MIN(DEPTH(I,J),TMP)
  !      END DO
  !    END DO

! Engineer cross-GSL shelf (shortcircuits GSL)
! close off GSL completely
      DEPTH(2,219)=0.
      TMP=.01*Z(5)
      DO J=211,218
        DO I=2,31
          DEPTH(I,J)=MIN(DEPTH(I,J),TMP)
        END DO
      END DO
      DEPTH(7,209)=.01*Z(7)
      DEPTH(8,209)=.01*Z(7)
      DEPTH(9,209)=.01*Z(9)
      DEPTH(10,209)=.01*Z(9)
      DO J=219,267
        ILO=31
        DO WHILE(DEPTH(ILO-1,J).GT..01*Z(5))
          ILO=ILO-1
        END DO
        N=3
        DO I=ILO,ILO+7
        N=N+1
          DEPTH(I,J)=MIN(DEPTH(I,J),.01*Z(N))
        END DO
      END DO
      ELSEIF (POSITION==2) THEN     
! close off area west of Greenland (questionable climate data)
      DO J=381,J1
        KB(25,J)=MIN(4,KB(25,J))
        KB(26,J)=MIN(3,KB(26,J))
        KB(27,J)=MIN(2,KB(27,J))
        KB(31:41,J)=0
      END DO
      END IF      
   ELSE IF (TRIM(CASE_NAME)=='IBE') THEN
      IF (POSITION==1) THEN    
! narrow sharp eastern shelf
  ! luc commented this out due to artificial current in the NE corner of IBE, 07/01/2014
  !	JJ=MIN(K1/2,J2/8)
  !    n=int(K1/JJ)
  !    DO I=1,JJ
  !      TMP=.01*Z(2*n*I+1)
  !      DO J=521,J1
  !        DEPTH(I0-I,J) =MIN(DEPTH(I0-I,J),TMP)
  !      END DO
  !    END DO
! With a 2:1 grid increment ratio, we match THREE overlap zones, because:
! a) we must not specify any cross-boundary flow where there is land
!   in either domain during the grid coupling process
! b) we sometimes need to upwind data from west of the overlap zones
!    and the depth matching is needed to get meaningful upstream averages
! This accuracy sacrifice in the high resolution grid near the interface
! gives simple and robust grid coupling based on sound CFD practice.
      DO J=2,J1
        TMP=.01*Z(2*KB(2,J)+1)
        DO I=2,4
          DEPTH(I,J)=TMP
        END DO
      END DO
      END IF
   ELSE IF (TRIM(CASE_NAME)=='GOM') THEN
      IF (POSITION==1) THEN
      DEPTH(136:138,162)=.01*Z(5)
      DEPTH(136,167)=.01*Z(5)
      DEPTH(136,165)=.01*Z(5)
      DEPTH(136,164)=.01*Z(5)
      DEPTH(135,164)=.01*Z(5)
! Set southernmost line of depths in 1/8 deg GOM region to zero to save
! storage in 1/16 deg grid comparison. This clips one zone from southern
! Carribean region and closes off a shallow narrow channel that cuts
! southward into Mexico from the western Caribbean. Effects are expected to
! be extremely small except locally. This allows a compact, highly focused
! high resolution nested GOM grid.
      DEPTH(2:120,49)=0.
! With a 2:1 grid increment ratio, we match THREE overlap zones, because:
! a) we must not specify any cross-boundary flow where there is land
!    in either domain during the grid coupling process
! b) we sometimes need to upwind data from west of the overlap zones
!    and the depth matching is needed to get meaningful upstream averages
! This accuracy sacrifice in the high resolution grid near the interface
! gives simple and robust grid coupling based on sound CFD practice.
! compare with eastern domain depths
      DO J=2,J1
        TMP=.01*Z(2*KB(I1,J)+1)
        DO I=I3,I1
          DEPTH(I,J)=TMP
        END DO
      END DO
      ENDIF
   ELSE IF (TRIM(CASE_NAME)=='GIB') THEN
      IF (POSITION==1) THEN
! use only if GIB depths are read from unit 21
      DO J=2,J1
        TMP=.01*Z(2*KB(I1,J)+1)
        DO I=I0-3,I0
          DEPTH(I,J)=TMP
        END DO
      END DO
      ENDIF
   ELSE IF (TRIM(CASE_NAME)=='GIB') THEN
      IF (POSITION==1) THEN
! use only if VIS depths are read from unit 22
      DO J=2,J1
        TMP=.01*Z(2*KB(2,J)+1)
        DO I=1,3
          DEPTH(I,J)=TMP
        END DO
      END DO
      ENDIF      
   END IF
  END SUBROUTINE COUPLE_BOUND_COND
END MODULE COUPLING

