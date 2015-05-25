SUBROUTINE INDATA
#include "namelist.in_GLOBAL"
  USE NETCDF
  USE TPS_NAMELIST_GLOBAL
  USE TPSDATA 
  USE GRID_VAR_GLOBAL, ONLY: A,XDEG,YV,YVDEG,YDEG,CS,OCS,DX,ODX,DY,ODY,CSV,OCSV,DXV,ODXV,DYV,ODYV,VAR_ALLOCATE_METRIC
  USE GRID_VAR_GLOBAL, ONLY: Y,VAR_ALLOCATE_METRIC_PRE
  USE GRID_VAR_GLOBAL, ONLY: KB,IU0,IV0,IN,IU,IV,IW,VAR_ALLOCATE_BATHY
  USE GRID_VAR_GLOBAL, ONLY: K0,K1,Z,ODZ,ODZW,VAR_ALLOCATE_ZFS
  USE GRID_VAR_GLOBAL, ONLY: DEPTH
  USE GRID_VAR_GLOBAL, ONLY:S1,T1,F,VAR_ALLOCATE_B
  USE OCN_PARA_GLOBAL !,ONLY:I0,I1,I2,J2
  USE XYFS_GLOBAL
  USE CONTROL_GLOBAL,ONLY:TLZ
  USE COUPLING
! Yu-Chiao 20110916
  USE DENSITY_EOS_GLOBAL  

  IMPLICIT NONE
  
  INTEGER::I,J,K,N,IT,NSUM,NN,ILO,IHI,L,NMONTH
  CHARACTER(LEN=10)::FMT
  REAL::TX(360,180),TY(360,180),TDIF
  REAL(8)::RHO(360,180,33)
  REAL,ALLOCATABLE :: TAUX(:,:),TAUY(:,:)
  CHARACTER(160)::TEMPDIR
  REAL::ZZ,DZZ,TMP,TMIN,TINC
  INTEGER(2),ALLOCATABLE :: INN(:,:),METERS(:,:)
  
  INTEGER :: OPENSTATUS
  INTEGER :: ISTAT
  
  REAL::XINC,NWIDTH,XLON,YLAT,YLON,TEMP,TMP1,TMP2,DSW,DSE,DNW,DNE,DN,DS,XX_TMP,YY_TMP
  INTEGER::II,JJ,JP,IP,M
  INTEGER::NX,NY,NSEC,IMAX,IMIN

  REAL(8)::TMPD,SLTD,PISD,p01,p02,TEST8,ccrho,ccTP,ccS,ccTPTP,ccSTP
  REAL(8),ALLOCATABLE :: D1(:,:,:)
  REAL,ALLOCATABLE :: TAVE(:,:,:),SAVE(:,:,:)         
  INTEGER :: NM_TMP

  INTEGER :: FNO_TPS,FNO_DATA
  INTEGER :: FNO_depth,FNO_ZKB,FNO_annu,FNO_bathy,FNO_initial,FNO_winds
  INTEGER :: FNO_KB,FNO_T0,FNO_IN,FNO_eos

  CALL VAR_ALLOCATE_BATHY

  ALLOCATE(TAUX(I2,J2),STAT=ISTAT) ; CALL ALLOC_STATUS("TAUX") ; TAUX=0.
  ALLOCATE(TAUY(I2,J2),STAT=ISTAT) ; CALL ALLOC_STATUS("TAUY") ; TAUY=0.
  ALLOCATE(INN(I0,J0),STAT=ISTAT)  ; CALL ALLOC_STATUS("INN")  ; INN=0.
  ALLOCATE(D1(I0,J0,K1),STAT=ISTAT); CALL ALLOC_STATUS("D1")
! Yu-Chiao 20110919
  ALLOCATE(TAVE(I0,J0,K1),STAT=ISTAT);CALL ALLOC_STATUS("TAVE");TAVE=0.
  ALLOCATE(SAVE(I0,J0,K1),STAT=ISTAT);CALL ALLOC_STATUS("SAVE");SAVE=0.

  CALL OPEN_TPS_NAMELIST_FILE(FNO_TPS)
  CALL OPEN_TPSDATA_NAMELIST_FILE(FNO_DATA)

  !===================================================
  !  read in another variables user defined for indata
  !---------------------------------------------------
  REWIND(FNO_TPS)
  READ(FNO_TPS,NML=LATERAL_BOUND_FLAGS)
  REWIND(FNO_TPS)
  READ(FNO_TPS,NML=INFLOW_BOUND_COND)

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! indata contains five major functions:          !
  ! 1) create bathymetry data                      !
  ! 2) calculate mask arrays, e.g. KB, IN, IU, ... !
  ! 3) give initial Temperature and Salinity fileds!
  ! 4) input wind data                             !
  ! 5) deal with inflow/outflow boundaries         !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 1) create bathymetry data                                                    !
  ! INTRODUCE SEVERAL FLAGS:                                                     !
  ! a) FL_RL_BATHY = 0: IDEAL BATHYMETRY WHICH CAN BE SPECIFIED BY THE CODE      !
  !    FL_RL_BATHY = 1: REAL  BATHYMETRY WHICH CAN BE SPECIFIED BY THE INPUT DATA!
  ! b) FL_ITPL = 0: USER SPECIFIED DATA, NO NEED OF INTERPOLATION                !
  !    FL_ITPL = 1: INTERPOLATE BATHYMETRY DATA TO GRID POINTS                   !
  ! c) FL_LBS = 0: NO LATITUDINAL BOUNDARY SHOLLOWING TREATMENT                  !
  !    FL_LBS = 1: DO LATITUDINAL BOUNDARY SHOLLOWING TREATMENT                  !
  !    IN GLOBAL CASE, LBS TREATMENT IS RECOMMENDED                              !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  READ(FNO_TPS,NML=BATHY_FLAGS)
  DIMENSION: IF (DIMEN == 2) THEN !LINE #94

    SELECT CASE(FL_RL_BATHY)  
    CASE(0)
        !I.E. ILE2D
        !NOTE THAT XX_TMP(YY_TMP) IS IN KM AND DX(1)(DY(1)) IS IN CM
        DO J=2,J1
          YY_TMP=-200.+DY(1)*(J+J2/4-1.5)/100000
          DO I=2,I1
            XX_TMP=-200.+DX(1)*(I-1.5)/100000
            TMP=SQRT(XX_TMP**2+YY_TMP**2)
            DEPTH(I,J)=.01*TLZ !THE UNIT OF DEPTH IS METER (TLZ IS IN CM)
            IF (TMP.LT.25.) DEPTH(I,J)=0.
          END DO
        END DO
    CASE(1)
        !(FL_RL_BATHY == 1) CAN BE IMPLEMENTED BY USER
    END SELECT
    
  ELSE DIMENSION !(DIMEN == 3)

    PRINT *,"INDATA case is 3D." 

    FL_RL_BATHY_3D:SELECT CASE(FL_RL_BATHY)
    CASE(0) FL_RL_BATHY_3D
      ! I.E. DTRAC(DOME) OR OTHER IDEALIZED BATHYMETRY 
      ! THAT CAN BE SPECIFIED THROUGH MATH. FORMULATION
        ZZ=600. !ZZ IS THE SHOLLOWEST DEPTH AT THE NORTH BOUND. (IN METER)
        DZZ=.01*TLX/I2/100. !GENERAL FORM
        DO J=J1,2,-1
          DO I=2,I1
            INN(I,J)=1
            DEPTH(I,J)=MIN(TLZ/100.,ZZ+(J1-J)*DZZ) !REVISED ON 2010/10/06
          END DO
        END DO
    CASE(1) FL_RL_BATHY_3D
      ! E.G. GLOBAL, NPB+TAI, MEDINA, REAL CASES
      ! BATHYMETRY IS SPECIFIED BY REAL DATA
  
      ! DIRECTLY READ USER SPECIFIED DATA (NO NEED OF INTERPOLATION)
      ! OR READ Itopo5 AND THEN INTERPOLATE BATHYMETRY DATA TO GRID POINTS
        FL_ITPL_3D: SELECT CASE(FL_ITPL)
        CASE(0) FL_ITPL_3D
          ! E.G. NPBTAI CASE DIRECTLY READS BATHY. DATA FROM INPUT FILE. NO NEED TO PERFORM INTERPOLATION
            CALL READ_DEPTH_ONLY !USER DEFINED SUBROUTINE
        CASE DEFAULT FL_ITPL_3D
          ! E.G. GLOBAL CASE READS INPUT DATA AND THEN DO INTERPOLATION
          !---------- READ BATHY INFO INCLUDING FILE, RES, AND FMT ----------!
            READ(FNO_DATA,NML=BATHYMETRY)    
          !---------- END OF READ BATHY. INFO ----------!

          ! SELECT RESOLUTION, E.G. BAYHY_RES = 5
            SELECT CASE(BATHY_RES)
            CASE(1) 
                NX=360*60
                NY=180*60
            CASE(2)
                NX=360*30
                NY=180*30
            CASE(5)
                NX=360*12
                NY=180*12
            END SELECT
            ALLOCATE(METERS(NX,NY))
      
          !---------- OPEN FILE, E.G. Itopo5, TO INPUT BATHY DATA ----------!
            print *,"Read METERS from bathy"
            CALL OPEN_OLD_BIN_FILE(TRIM(BATHY_FILE),TRIM(BATHY_FMT),FNO_bathy)
            READ(FNO_bathy) METERS ! Read FULL etopo5 FILE INTO MEMORY
            CLOSE(FNO_bathy)
          !---------- END OF READ BATHY. DATA ----------!
  
          !---------- START OF INTERPOLATION ----------!
            DO I=2,I1
              DEPTH(I,2)=X0DEG+(DXMNUT*(I-1.5))/60.
            END DO

            XINC=.2*DXMNUT
            DO JJ=1,J0
            ! YLAT is measured from north pole.  First point is at pole.
            ! XLON is measured units 1/12 deg east, from Grenwich, with first point
            ! at 1/12 deg east (METERS(1,J) is at 1/12 deg east; METERS(I,1) is
            ! at North Pole, or YDEG=90)
              YLAT=12*(90-YDEG(JJ))+1
            ! Average depths in lat-long window
              NWIDTH=DXMNUT/5.
              NWIDTH=NWIDTH/2
              J=YLAT
              JP=J+1
              XLON=12*X0DEG-.5*XINC
              DO II=2,I1
                XLON=XLON+XINC
                I=XLON
                TMP1=XLON-I
                TMP2=1.-TMP1
                IP=MOD(I,4320)+1
                I=MOD(I-1,4320)+1
                DSW=METERS(I,JP)
                DSE=METERS(IP,JP)
                DNW=METERS(I,J)
                DNE=METERS(IP,J)
                DN=TMP1*DNE+TMP2*DNW
                DS=TMP1*DSE+TMP2*DSW
                DEPTH(II,JJ)=(YLAT-J)*DS+(1-YLAT+J)*DN ! interpolation
                DEPTH(II,JJ)=-DEPTH(II,JJ)
                DEPTH(II,JJ)=MAX(0.,DEPTH(II,JJ))
                INN(II,JJ)=1
                IF (DEPTH(II,JJ).LT..01*Z(2)) INN(II,JJ)=0
              END DO
            END DO
          !---------- END OF INTERPOLATION ----------!    
        END SELECT FL_ITPL_3D
    END SELECT FL_RL_BATHY_3D
  END IF DIMENSION !(2D/3D CASE)
  
!---------- DO LATITUDINAL BOUNDARY SHALLOWING TREATMENT ----------!
! LBS TREATMENT IS RECOMMENDED FOR GLOBAL AND NPB CASES (NO NEED FOR TAI)
! ORIGINALLY, SETTINGS IN THE GLOBAL CASE AND NPB CASE ARE SLIGHTLY DIFFERENT
! WE USE THE SAME TREATMENT

! Preprocessor test

#ifdef FLAG_test1 
PRINT*,'FLAG test is sucessful'
#endif


#ifdef FLAG_LBS_S
!=======================Jay Modified===========================
      JJ=MIN(K1/2,J2/8)
      n=int(K1/JJ)
      DO J=2,JJ+1
        TMP=.01*Z(2*n*(J-1)+1)
        DO I=2,I1
          DEPTH(I,J) =MIN(DEPTH(I,J),TMP)
        END DO
      END DO
!=====================================================
#endif
#ifdef FLAG_LBS_N
!=======================Jay Modified===========================
      JJ=MIN(K1/2,J2/8)
      n=int(K1/JJ)
      DO J=2,JJ+1
        TMP=.01*Z(2*n*(J-1)+1)
        DO I=2,I1
          DEPTH(I,J0-J+1)=MIN(DEPTH(I,J0-J+1),TMP)
        END DO
      END DO
!=====================================================
#endif

!---------- END OF LATITUDINAL BOUNDARY SHALLOWING TREATMENT ----------!

!---------- COVERT WATER TO LAND ----------!
  TMP=.01*Z(5)   
  TMP1=.01*Z(3) 
  DO J=2,J1
    DO I=2,I1
      IF (DEPTH(I,J).LT.TMP1) DEPTH(I,J)=0.
      !convert land to water
      IF (DEPTH(I,J).GT.TMP1) DEPTH(I,J)=MAX(DEPTH(I,J),TMP)
    END DO
  END DO
!---------- END OF COVERSION ----------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SPECIAL TREATMENT FOR NPB, TAI and MEDINA (CASE DEPENDENT)                          !
! FOR NPB: hardwired: artificial closed boundary and shelf west of Kamchitka!
! FOR TAI: extend NPB artificial shelf along northern boundary              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (TRIM(CASE_NAME) == 'TAI') THEN
    OPEN(20,file='KBI1J')
    READ(20,25) (KB(I1,J),J=2,J1)
    CLOSE(20)
    CALL COUPLE_BOUND_COND(1)
  END IF
  IF (TRIM(CASE_NAME) == 'NPB'.OR.(TRIM(CASE_NAME) == 'NAB')) THEN
    CALL COUPLE_BOUND_COND(1)
  END IF          
  IF (TRIM(CASE_NAME) == 'IBE') THEN
    OPEN(20,file='KB2J')
    READ(20,25) (KB(2,J),J=2,J1)
    CLOSE(20)
    CALL COUPLE_BOUND_COND(1)
  END IF
  IF (TRIM(CASE_NAME) == 'GOM') THEN
    OPEN(20,file='KBI1J')
    READ(20,25) (KB(I1,J),J=2,J1)
    CLOSE(20)
    CALL COUPLE_BOUND_COND(1)
  END IF
  IF (TRIM(CASE_NAME) == 'GIB') THEN
    OPEN(20,file='KBI1J')
    READ(20,25) (KB(I1,J),J=2,J1)
    CLOSE(20)
    CALL COUPLE_BOUND_COND(1)
  END IF
 25  FORMAT(50I3)
 
!---------- CALCULATE TOTAL WET POINTS ----------!
  M=0
  N=0
  TMP=0.
  TEMP=.01*TLZ
  DO J=2,J1
    DO I=2,I1
      M=M+INN(I,J)
      TMP=MAX(TMP,DEPTH(I,J))
      IF (DEPTH(I,J).GT.TEMP) N=N+1
    END DO
  END DO
  WRITE(*,103) M,N,TEMP,TMP
  103 FORMAT(I6,' total wet points,',I6,' deeper than ',F5.0,' m, deepest= ',F6.0,' m')
!---------- END OF WET-POINT CALCULATION ----------!

! Reset KB for summation
  DO J=2,J1
    KB(2,J)=0
    KB(I1,J)=0
  END DO

! CONVERT TO CM
  DO J=2,J1
    DO I=2,I1
      DEPTH(I,J)=MIN(TLZ,100.*DEPTH(I,J))
    END DO
  END DO

!in GLOBAL case, adapt some points in the Gulf of Guinea
!luc november 2013
  IF (TRIM(CASE_NAME) == 'GLOBAL') THEN
    !!strong modification for the gulf of guinea, following core2 blowup there
    !i=47
    !depth(1:3,46:53)=depth(2,46)
    !depth(I2:I0,46:53)=depth(2,46)
    !depth(2:4,54)=0
    !softer modification
    depth(2,50:51)=depth(2,49)
    depth(2:4,54)=0
  END IF
	
!---------- CHECK PERIODIC B.C. FOR DEPTH ----------!2010/10/13
  IF (LOPENW==1) THEN
    DO J=2,J1
      DEPTH(1,J)=DEPTH(I1,J)
    END DO
  END IF
  IF (LOPENE==1) THEN
    DO J=2,J1
      DEPTH(I0,J)=DEPTH(2,J)
    END DO
  END IF
  IF (LOPENS==1) THEN
    DO I=2,I1
      DEPTH(I,1)=DEPTH(I,J1)
    END DO
  END IF
  IF (LOPENN==1) THEN
    DO I=2,I1
      DEPTH(I,J0)=DEPTH(I,2)
    END DO
  END IF
!---------- END OF PERIODIC B.C. FOR DEPTH ----------!  
  print *,"write DEPTH to depth file"
  TEMPDIR=TRIM(WORKDIR)//'/TEMP_GLOBAL/depth'
  CALL OPEN_BIGENDIAN_NEW_BIN_FILE(TRIM(TEMPDIR),FNO_depth)
  WRITE(FNO_depth) DEPTH
  CLOSE(FNO_depth)
!---------- END OF BATHYMETRY DATA ----------!
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2) calculate mask arrays, e.g. KB, IN, IU, ... !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! Determine "logical grid" depth array KB(I,J)
! KB(I,J)=number of layers at pressure point (I,J)
  DO K=1,K1
    ZZ=.5*(Z(2*K-1)+Z(2*K+1))
    DO J=2,J1
      DO I=2,I1
        ! If DEPTH is deeper than mid-layer (ZZ), increment KB by 1
        L=DEPTH(I,J)/ZZ        
        KB(I,J)=KB(I,J)+MIN(1,L)
      END DO        
    END DO
  END DO

  ! hardwired: deepen shelf to open connection to Sea of Okhotsk
  IF (TRIM(CASE_NAME) == 'NPB') THEN
    CALL COUPLE_BOUND_COND(2) ! THIS SUBROUTINE CAN WORK WHEN KB IS KNOWN
    ! KB(I,J)=MIN(KB(I,J),K); ! DEPTH(I,J)=Z(2*KB(I,J)+1)
    TEMPDIR=TRIM(WORKDIR)//'/TEMP_GLOBAL/depth'
    print *,"write DEPTH to depth file"
    CALL OPEN_BIGENDIAN_NEW_BIN_FILE(TRIM(TEMPDIR),FNO_depth)
    WRITE(FNO_depth) DEPTH
    CLOSE(FNO_depth)
  END IF
  IF (TRIM(CASE_NAME) == 'NAB') THEN
    CALL COUPLE_BOUND_COND(2)
  ENDIF
! ==========================
! Derive land/sea mask array
! ==========================
  DO J=2,J1
    DO I=2,I1
      IT=KB(I,J)
      IF (IT/=0) THEN
        DO K=1,IT
          IN(I,J,K)=1
        END DO
      ELSE
        KB(I,J)=1  
      END IF
    END DO
  END DO
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !---------- CASE DEPENDENT TREATMENTS ----------!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! THE FOLLOWING ARE SPECIAL TREATMENTS ILE2D, NPBTAI, AND GLOBAL MODEL, RESPECTIVELY
  IF (K1.LE.3) THEN
    !ILE2D CASE RESETS KB FOR USE AT BOTTOM DRAG LOCATIONS
    DO J=2,J1
      DO I=2,I1
        K=1  
        DO WHILE (IN(I,J,K).EQ.1.AND.K.LT.K1)
          K=K+1
        END DO  
        KB(I,J)=MAX(1,K-1)
      END DO
    END DO
  END IF
  
  ! Eliminate isolated pools outside connected southern hemisphere ocean
  ! Use IN=2 to denote water that is connected to OR within ocean 
  REALBATHY_3D: IF ((DIMEN == 3) .AND. (FL_RL_BATHY == 1)) THEN !LINE #352       
    IF (TRIM(CASE_NAME) == 'NPB'.OR.TRIM(CASE_NAME) == 'NAB'.OR.TRIM(CASE_NAME) == 'IBE') THEN
      IN(I0/2,J0/2,1)=2 !FOR NPB
    ELSE IF (TRIM(CASE_NAME) == 'TAI'.OR.TRIM(CASE_NAME) == 'GOM') THEN
      IN(I1,J0/2,1)=2   !FOR TAI
    ELSE IF (TRIM(CASE_NAME) == 'MED') THEN
      IN(I0-30,30,1)=2   !FOR MED
    ELSE !CASE_NAME == 'GLOBAL'
        IN(2*I0/3,J0/2,1)=2 ! start with single point in central equatorial Pacific
        IN(2,J1,1)=2        ! we need a second point in the Arctic ocean due to periodic logical grid boundary at 0 deg longitude    
    END IF
      
    N=0  
    DO WHILE(.TRUE.)
      N=N+1
      NSUM=0
      DO J=2,J1
        DO I=2,I1
          NN=IN(I+1,J,1)
          IF (IN(I+1,J,1).EQ.1.AND.IN(I,J,1).EQ.2) IN(I+1,J,1)=2
          NSUM=NSUM+IN(I+1,J,1)-NN
          NN=IN(I-1,J,1)
          IF (IN(I-1,J,1).EQ.1.AND.IN(I,J,1).EQ.2) IN(I-1,J,1)=2
          NSUM=NSUM+IN(I-1,J,1)-NN
          NN=IN(I,J+1,1)
          IF (IN(I,J+1,1).EQ.1.AND.IN(I,J,1).EQ.2) IN(I,J+1,1)=2
          NSUM=NSUM+IN(I,J+1,1)-NN
          NN=IN(I,J-1,1)
          IF (IN(I,J-1,1).EQ.1.AND.IN(I,J,1).EQ.2) IN(I,J-1,1)=2
          NSUM=NSUM+IN(I,J-1,1)-NN
        END DO
        DO I=I1,2,-1
          NN=IN(I+1,J,1)
          IF (IN(I+1,J,1).EQ.1.AND.IN(I,J,1).EQ.2) IN(I+1,J,1)=2
          NSUM=NSUM+IN(I+1,J,1)-NN
          NN=IN(I-1,J,1)
          IF (IN(I-1,J,1).EQ.1.AND.IN(I,J,1).EQ.2) IN(I-1,J,1)=2
          NSUM=NSUM+IN(I-1,J,1)-NN
          NN=IN(I,J+1,1)
          IF (IN(I,J+1,1).EQ.1.AND.IN(I,J,1).EQ.2) IN(I,J+1,1)=2
          NSUM=NSUM+IN(I,J+1,1)-NN
          NN=IN(I,J-1,1)
          IF (IN(I,J-1,1).EQ.1.AND.IN(I,J,1).EQ.2) IN(I,J-1,1)=2
          NSUM=NSUM+IN(I,J-1,1)-NN
         END DO
      END DO
      WRITE(*,156) NSUM,N
156   FORMAT(I6,' changes during iteration ',I6)
      IF (NSUM.NE.0) CYCLE
      EXIT
    ENDDO
   
    DO J=2,J1
      DO I=2,I1
        N=IN(I,J,1)
        IF (IN(I,J,1).NE.2) KB(I,J)=0
      END DO
    END DO
    DO K=1,K1
      DO J=2,J1
        DO I=2,I1
          IN(I,J,K)=0
        END DO
      END DO
    END DO
    
    ! redefine IN with exterraneous pools eliminated
    DO J=2,J1
      DO I=2,I1
        IT=KB(I,J)
        IF (IT.NE.0) THEN
          DO K=1,IT
            IN(I,J,K)=1
          ENDDO
        ELSE
          KB(I,J)=1
        ENDIF
      ENDDO
    ENDDO
    !---------- END OF ELIMINATION ----------!
      
    !---------- CHECK PERIODIC B.C. FOR KB AND IN  ----------!
    !CALL BOUND_COND(1)
    IF (LOPENW==1) THEN ! PERIODIC
      DO J=1,J0
        KB(1,J)=KB(I1,J)
        DO K=1,K1
          IN(1,J,K)=IN(I1,J,K)
        END DO
      END DO
    END IF
    
    IF (LOPENE==1) THEN ! PERIODIC
      DO J=1,J0
        KB(I0,J)=KB(2,J)
        DO K=1,K1
          IN(I0,J,K)=IN(2,J,K)
        END DO
      END DO
    END IF
  IF (TRIM(CASE_NAME) == 'NAB') THEN
      OPEN(20,file='../GOM/KBI1J')
      OPEN(21,file='../IBE/KB2J')
      WRITE(20,25) ((IN(2,J,1)*KB(2,J),N=1,2),J=40,J1)
      CLOSE(20)
      WRITE(21,25) ((IN(I1,J,1)*KB(I1,J),N=1,2),J=2,J1)
      CLOSE(21)
  ELSEIF(TRIM(CASE_NAME) == 'MED') THEN
      OPEN(20,file='../GIB/KBI1J')
      WRITE(20,25) ((IN(2,J,1)*KB(2,J),N=1,3),J=32,J1)
      CLOSE(20)
  ENDIF
  END IF REALBATHY_3D
    
  CALL IU_V_W_CAL
  
  print *,"Write data to ZKB"
  TEMPDIR=TRIM(WORKDIR)//'/INPUT_GLOBAL/ZKB'
  CALL OPEN_BIGENDIAN_NEW_BIN_FILE(TRIM(TEMPDIR),FNO_ZKB)
  WRITE(FNO_ZKB) Z,ODZ,ODZW  !REMOVE TLZ
  WRITE(FNO_ZKB) KB,IN,IU,IV,IW,IU0,IV0
  WRITE(FNO_ZKB) DEPTH
  CLOSE(FNO_ZKB)
  !---------- END OF MASK-ARRAY CALCULATION ----------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3) give initial Temperature and Salinity fileds                              !
! INTRODUCE SEVERAL FLAGS:                                                     !
!  FL_RL_TS = 0: IDEAL CASE                                                  !
!  FL_RL_TS = 1: USE REAL DATA                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  READ(FNO_TPS,NML=TSC_FLAGS)
!Yu-Chiao 10152011
#ifdef FLAG_TS_T
      !E.G. ILE2D
      !E.G. DTRAC CASE ONLY CONSIDERS T (T IS ALSO AN IDEALIZED CASE)
      REAL_T: SELECT CASE(FL_RL_TS)
      CASE(0) REAL_T
        ! IN TWO STEPS (I). INITIALIZE T FIELD & (II). SPECIFY T AT INFLOW
  
        ! STEP(I). INITIALIZE T FIELD
        ! linear stratification
        ! DOME initial stratification was 2 kg/m3 over 3600m depth
        ! this is equivalent to 10 deg with our e.o.s.
        ! Our DOME stratification for 2160m/3600m of DOME depth is thus
        ! 3/5 X 10 deg
        ! GENERAL FORM: TLZ cm/360000 cm 
  
          !TMAX=12. !TS !READ FROM NML FILE
          TMIN=TMAX-10.*TLZ/360000. !GENERAL FORM 2010/10/07
          TINC=-(TMAX-TMIN)/K1  
          TMP=TMAX
          DO K=1,K1
            DO J=1,J0
              DO I=1,I0
                T1(I,J,K)=TMP
              END DO
            END DO
            TMP=TMP+TINC
            N=.01*Z(2*K)            
            WRITE(*,301) K,N,T1(1,1,K)
          END DO 
301   FORMAT('k,depth(m), temperature:',2I5,F6.2)
      
        ! STEP(II). SPECIFY T AT INFLOW
        ! DOME density current problem patterned after Ezer and Mellor
        ! domain width=1100km
        ! 600m depth sill, total depth is 2200m
        ! warmest layer 1 inflow matches interior layer 1 i.c.'s
      
          ILO=8./11.*I2
          IHI=9./11.*I2
          !N=60000./TLZ*K1
          !TINC=(-10.*TLZ/360000.)/(N-1)
          N=6.*K1/22.
          TINC=(-10.*22./36.)/(N-1)
          DO K=1,N
            ZZ=.01*Z(2*K)
            TDIF=4.*EXP((-ZZ)/600.)
            TMP=TDIF/((IHI-ILO))**2
            DO I=ILO,IHI
              T1(I,J0,K)=TMAX-TMP*(IHI-I)**2
            END DO
            DO I=2,ILO
              T1(I,J0,K)=T1(ILO,J0,K)
            END DO
            DO I=IHI+1,I1
              T1(I,J0,K)=T1(IHI,J0,K)
            END DO
            TMAX=TMAX+TINC      !ZWM TMAX do change
          END DO
 
          print *,"Write T1 to initial" 
          TEMPDIR=TRIM(WORKDIR)//'/INPUT_GLOBAL/initial'
          CALL OPEN_BIGENDIAN_NEW_BIN_FILE(TRIM(TEMPDIR),FNO_initial)
          WRITE(FNO_initial) T1
          CLOSE(FNO_initial)
          
      CASE DEFAULT REAL_T
          ! NEED FURTHER TREATMENT
      END SELECT REAL_T
#endif
! Yu-Chiao 10152011
#ifdef FLAG_TS_T
#ifdef FLAG_TS_S
      REAL_TS: SELECT CASE(FL_RL_TS)
      CASE(0) REAL_TS
          !IDEALIZED T&S CAN BE FURTHER IMPLEMENTED
      CASE DEFAULT REAL_TS 
        ! ==================================================
        ! Prepare seasonal S,T climatology
        ! ==================================================
          READ(FNO_DATA,NML=LEVITUS)
          !READ(FNO_DATA,NML=LEV_FMT)
          !READ(FNO_DATA,NML=LEV_DEPTH)
          !READ(FNO_DATA,NML=TS_FILE)
  
        ! =================
        ! Output data files
        ! =================
        ! Annual cycle Levitus climatology on model grid
        ! S,T at surface and lateral boundaries
          print *,"Write S1, T1 to annualevitus"
          TEMPDIR=TRIM(WORKDIR)//'/TEMP_GLOBAL/annualevitus'
          CALL OPEN_BIGENDIAN_NEW_BIN_FILE(TRIM(TEMPDIR),FNO_annu)
          TAVE=0.;SAVE=0.
          DO JJ=1,4
            CALL TSLEV(SEASON(JJ),FL_EOS_OP)
            DO II=1,3
              NMONTH=II+(JJ-1)*3
              WRITE(*,40) NMONTH
              CALL EOS_COMPUTATION_STABILITY1(T1,S1,Z,D1,FL_EOS_OP,0,0)
              WRITE(*,41)NMONTH,minval(T1),maxval(T1)
!              PRINT*,'EOS test,D1',SUM(D1),MAXVAL(D1),MINVAL(D1),FL_EOS_OP
            ! CHECK PERIODIC B.C. !  CALL BOUND_COND(4)
              IF (LOPENW==1) THEN
                DO K=1,K1           
                DO J=2,J1
                  T1(1,J,K)=T1(I1,J,K)
                  S1(1,J,K)=S1(I1,J,K)
                END DO
                END DO
              END IF
              IF (LOPENE==1) THEN
                DO K=1,K1           
                DO J=2,J1
                  T1(I0,J,K)=T1(2,J,K)
                  S1(I0,J,K)=S1(2,J,K)
                END DO
                END DO
              END IF
              DO K=1,K1
                DO J=1,J0
                  DO I=1,I0
                    SAVE(I,J,K)=SAVE(I,J,K)+S1(I,J,K)/12.
                    TAVE(I,J,K)=TAVE(I,J,K)+T1(I,J,K)/12.
                  END DO
                END DO
              END DO
              WRITE(FNO_annu) S1,T1
            END DO
          END DO
          CALL EOS_COMPUTATION_STABILITY1(TAVE,SAVE,Z,D1,FL_EOS_OP,0,1)

          CLOSE (FNO_annu)

          !@RESTART_CONDITION: IF (LRSTRT==0) THEN
          ! FOR THE GLOBAL CASE, SMOOTH PROCEDURE WAS NEGLATED!
          ! STARTING AND ENDDED POINTS FOR UNSTABLE-PT CHECKING
          
          print *, "write S1, T1 data to initial"
          TEMPDIR=TRIM(WORKDIR)//'/INPUT_GLOBAL/initial'
          CALL OPEN_BIGENDIAN_OLD_BIN_FILE(TRIM(TEMPDIR),FNO_initial)
          WRITE(FNO_initial)SAVE,TAVE
          CLOSE (FNO_initial)
40    FORMAT('MONTH',I3)
41    FORMAT('Month',I3,' Tmin,Tmax=',2F6.2)

      END SELECT REAL_TS
#endif
#endif

  !---------- END OF TEMPERATURE AND SALINITY FIELD ----------!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 4) input wind data                                                           !
  ! INTRODUCE SEVERAL FLAGS:                                                     !
  ! a) FL_WD_ON = 0: NO WIND EFFECTS                                             !
  !    FL_WD_ON = 1: CONSIDER WIND EFFECTS                                       !
  ! b) FL_RL_WD = 0: IDEAL WIND                                                  !
  !    FL_RL_WD = 1: REAL WIND DATA                                              !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REWIND(FNO_TPS)
  READ(FNO_TPS,NML=WIND_FLAGS)
  WIND_ON_OR_OFF: SELECT CASE(FL_WD_ON)
  CASE(1) WIND_ON_OR_OFF
      REALWIND: SELECT CASE(FL_RL_WD)
      CASE(0) REALWIND
          !IDEALIZED WIND CONDITION CAN BE IMPLEMENTED FURTHER
      CASE DEFAULT REALWIND
        ! ===================== !
        ! Prepare MONTHly winds !
        ! ===================== !
          READ(FNO_DATA,NML=HELLERMAN)
          !READ(FNO_DATA,NML=TX_FILE)
          !READ(FNO_DATA,NML=TY_FILE)
      
        ! ===================== !
        ! Output data files     !
        ! ===================== !
          print *,"Write data to winds"
          TEMPDIR=TRIM(WORKDIR)//'/INPUT_GLOBAL/winds'
          CALL OPEN_BIGENDIAN_NEW_BIN_FILE(TRIM(TEMPDIR),FNO_winds)  
          DO J=1,4
            DO I=1,3
              NMONTH=I+(J-1)*3
              WRITE(*,'(I3)') NMONTH
              CALL TAUXY(TXC(NMONTH),TYC(NMONTH),NMONTH)  !  READ
              print *, "write data to winds"
            ! Write top level Levitus climatology and Hellerman winds on model grid
              WRITE(FNO_winds) TAUX,TAUY
            END DO
          END DO
          CLOSE(FNO_winds)
      END SELECT REALWIND
  END SELECT WIND_ON_OR_OFF
  !---------- END OF WIND EFFECT FIELD ----------!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 5) deal with inflow/outflow boundaries         !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 2010/10/08
  ! E.G. IN DTRAC CASE, BOUNDS IS USED TO OUTPUT THE FOLLOWING B.C.

  CALL BOUNDS

  !---------- END OF INFLOW/OUTFLOW B.C. ----------!

  CALL CLOSE_TPS_NAMELIST_FILE
  CALL CLOSE_TPSDATA_NAMELIST_FILE

  CALL NC_OUTPUT_KB
  CALL NC_OUTPUT_T0
  CALL NC_OUTPUT_IN

  CONTAINS
  
! ----------------------------------------------------------------------
SUBROUTINE READ_DEPTH_ONLY
!USER DEFINED SUBROUTINE
  USE INIT_VAR_GLOBAL
  USE GRID_VAR_GLOBAL, ONLY: RHO,P0,PX,PY,U1,U2,V1,V2,S1,S2,T1,T2,P,ULF,VLF,SLF,TLF,EV,HV,F,TANPHI,SUMIN
  USE GRID_VAR_GLOBAL, ONLY: VNOR,VSUD,UEST,UWST
  USE GRID_VAR_GLOBAL, ONLY: KB,IU0,IV0,IN,IU,IV,IW
  USE GRID_VAR_GLOBAL, ONLY: U,V,W
  USE GRID_VAR_GLOBAL, ONLY: SSURF,TSURF,SSSP,SNSP,TSSP,TNSP
!  USE CONTROL_GLOBAL
  USE GRID_VAR_GLOBAL, ONLY: A,XDEG,YV,YVDEG,YDEG,CS,OCS,DX,ODX,DY,ODY,CSV,OCSV,DXV,ODXV,DYV,ODYV
  USE GRID_VAR_GLOBAL, ONLY: Y
  USE GRID_VAR_GLOBAL, ONLY: Z,ODZ,ODZW
  USE SCA_GLOBAL    ! Derived scalars
  USE GRID_VAR_GLOBAL, ONLY: RINV,RINV1,DUM0,DUM1,DUM2,XX,H,X,S,AL,AB,AC,AR,AT,SRC,CL,CB,CC,CR,CT,IE
  USE GRID_VAR_GLOBAL, ONLY: DEPTH
  USE TIMCOM_GENERAL
  USE OPFREQ    !ZWM
  USE CONTROL_GLOBAL, ONLY: AROUT
!  USE GRID_VAR_GLOBAL, ONLY:AROUT
  USE GRID_VAR_GLOBAL, ONLY: VBK,HBK
  USE GRID_VAR_GLOBAL, ONLY: SXYCLI,TXYCLI

  IMPLICIT NONE
  
  CHARACTER(160)::TEMPDIR
  REAL::TMIN,TMAXX,TMP1,TMP2,SUMINN,TMP,SUMS,SUMT,TEMP,RLON,RLAT
!  INTEGER, INTENT(IN) :: POSITION !select position when BC is called
!  INTEGER, INTENT(IN), OPTIONAL :: NM !required in some cases
  INTEGER :: POSITION !select position when BC is called
  INTEGER :: NMONTH !required in some cases

  INTEGER :: I,J,K,JLO,IMIN,JMIN,KMIN,IMAX,JMAX,KMAX,STATUS(2)
  !REAL, ALLOCATABLE :: SAVG(:),TAVG(:),TTAVG(:,:,:),SSAVG(:,:,:)
  
  REAL:: TTAVG(J2,K1,12),SSAVG(J2,K1,12),SAVG(K1),TAVG(K1),AVG1,AVG2,AREAXY
  INTEGER :: FNO_depth

  IF (TRIM(CASE_NAME)=='NPB') THEN
    DEPTH=0.
    CALL OPEN_BIGENDIAN_OLD_TXT_FILE('depth_NPB',FNO_depth)
    ! WHY FROM 1 TO 200
    DO I=1,200
      DO J=1,J0
    READ(FNO_depth,11) RLON,RLAT,DEPTH(I,J)
      END DO
    END DO
    DO I=1,I0
      DO J=1,J0
        READ(FNO_depth,11) RLON,RLAT,DEPTH(I,J)
        DEPTH(I,J)=MAX(0.,-DEPTH(I,J))
      END DO
    END DO
    CLOSE (FNO_depth)
  ELSE IF (TRIM(CASE_NAME)=='TAI') THEN
    DEPTH=0.
    CALL OPEN_BIGENDIAN_OLD_TXT_FILE('depth_TAI',FNO_depth)
    DO I=2,I1
      DO J=1,J1
        READ(FNO_depth,11) RLON,RLAT,DEPTH(I,J)
        DEPTH(I,J)=MAX(0.,-DEPTH(I,J))
      END DO
    END DO
    CLOSE (FNO_depth)
  ELSE IF (TRIM(CASE_NAME)=='NAB') THEN
    DEPTH=0.
    CALL OPEN_BIGENDIAN_OLD_BIN_FILE('depth_NAB',FNO_depth)
    READ(FNO_depth)DEPTH
    CLOSE (FNO_depth)
  ELSE IF (TRIM(CASE_NAME)=='MED') THEN
    DEPTH=0.
    CALL OPEN_BIGENDIAN_OLD_BIN_FILE('depth_MED',FNO_depth)
    READ(FNO_depth)DEPTH
    CLOSE (FNO_depth)
  ELSE IF (TRIM(CASE_NAME)=='IBE') THEN
    DEPTH=0.
    CALL OPEN_BIGENDIAN_OLD_BIN_FILE('depth_IBE',FNO_depth)
    READ(FNO_depth)DEPTH
    CLOSE (FNO_depth)
  ELSE IF (TRIM(CASE_NAME)=='GOM') THEN
    DEPTH=0.
    CALL OPEN_BIGENDIAN_OLD_BIN_FILE('depth_GOM',FNO_depth)
    READ(FNO_depth)DEPTH
    CLOSE (FNO_depth)
  ELSE IF (TRIM(CASE_NAME)=='GIB') THEN
    DEPTH=0.
    CALL OPEN_BIGENDIAN_OLD_BIN_FILE('depth_GIB',FNO_depth)
    READ(FNO_depth)DEPTH
    CLOSE (FNO_depth)
  END IF
  11   FORMAT(2F10.4,F10.2)  
END SUBROUTINE READ_DEPTH_ONLY
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
SUBROUTINE IU_V_W_CAL
!calculate IU/V/W
  USE GRID_VAR_GLOBAL, ONLY: KB,IU0,IV0,IN,IU,IV,IW
  USE OCN_PARA_GLOBAL
  USE TIMCOM_GENERAL

  INTEGER :: I,J,K

! other mask arrays
  DO K=1,K1
    DO J=1,J0
      DO I=1,I1
        IU(I,J,K)=IN(I,J,K)*IN(I+1,J,K)
      END DO
    END DO
    DO J=1,J1
      DO I=1,I0
        IV(I,J,K)=IN(I,J,K)*IN(I,J+1,K)
      END DO
    END DO
  END DO
  DO K=1,K2
    DO J=1,J0
      DO I=1,I0
        IW(I,J,K+1)=IN(I,J,K)*IN(I,J,K+1)
      END DO
    END DO
  END DO

! CHECK PERIODIC B.C. FOR IU
  IF (LOPENE==1) THEN   ! PERIODIC
    DO K=1,K1
      DO J=1,J0
        IU(I0,J,K)=IU(2,J,K)
      END DO
    END DO
  END IF

  IF ((LOPENW==1).AND.(LOPENE==1)) THEN
    DO J=2,J1
      DO I=1,I0
        IU0(I,J)=IU(I,J,1)
        IU(I,J,1)=1
      END DO
    END DO
  ELSE
    DO J=2,J1
      DO I=2,I2
        IU0(I,J)=IU(I,J,1)
        IU(I,J,1)=1
      END DO
    END DO
  ENDIF

  DO J=2,J2
    DO I=2,I1
      IV0(I,J)=IV(I,J,1)
      IV(I,J,1)=1
    END DO
  END DO

END SUBROUTINE IU_V_W_CAL
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
SUBROUTINE TSLEV(FULLNAME,FL_EOS_OP)
IMPLICIT NONE

  INTEGER::I,J,II,JJ,IP,JP,K,LT,LB,N    !  TEMP VARS
  CHARACTER(128)::FULLNAME    !  BATHYMATRY DATA FILE WITH FULL PATH TO BE READ IN
  REAL:: ZLEV,SUMM    !  TEMP VARS
  REAL::RHOK,RHOKM,XINC
  REAL::YLAT,XLON,XPLUS,TMP1,TMP2,DNW,DNE,DSW,DSE,DN,DS,TMAX,TMIN
  REAL::Z1,Z2,ZC,C1,C2

! Yu-Chiao 20110927
  REAL::TLEV(360,180,33),SLEV(360,180,33)
  REAL(8)::RHO1(360,180,33)
  REAL,ALLOCATABLE :: T(:,:,:),S(:,:,:)
  INTEGER :: STATUS(4),FNO_fix,FNO
! Yu-Chiao 20110915
  REAL(8)::TMPD,SLTD,PISD,TEST8,ccrho,ccTP,ccS,ccTPTP,ccSTP  
  INTEGER,INTENT(IN)::FL_EOS_OP
  REAL::ZZLEV(33)

  allocate (T(I0,J0,33),STAT=status(1)); allocate (S(I0,J0,33),STAT=status(2))
  if (status(1)/=0 .or. status(2)/=0) STOP "Cannot allocate memory"

  T=0.;S=0.

  TEMPDIR=TRIM(WORKDIR)//'/TEMP_GLOBAL/fixedlev'
  CALL OPEN_BIGENDIAN_OLD_BIN_FILE(TRIM(FULLNAME),FNO_fix)
  CALL OPEN_BIGENDIAN_NEW_BIN_FILE(TRIM(TEMPDIR),FNO)

! ===================================================
! Read modelfied (with landfills) Levitus climatology
! and fix negative stratification areas
! ===================================================

! Yu-Chiao 20110927
  DO K=1,33
  READ(FNO_fix)  ((TLEV(I,J,K),I=1,360),J=1,180),((SLEV(I,J,K),I=1,360),J=1,180)
  END DO

! NZDAT are Levitus levels
! Yu-Chiao 20110927
!!  ZLEV = FLOAT(NZDAT(1))
! Yu-Chiao 20110915
! EOS modification, controlled by the FL_EOS_OP

  TEST8=0.d0
! Yu-Chiao 20110927
  DO K=1,33
    DO J=1,180
      DO I=1,360
        IF ((TLEV(I,J,K) .LT. -1.8)) TLEV(I,J,K) = 0.0
        IF ((SLEV(I,J,K) .LT. -1.8)) SLEV(I,J,K) = 0.0
        PISD=FLOAT(NZDAT(K))
        TMPD=TLEV(I,J,K) 
        SLTD=SLEV(I,J,K)
        ! Yu-Chiao 20110929
        ! field temperature to potential temperature
        TLEV(I,J,K)=theta(SLTD,TMPD,PISD,TEST8)
      END DO
    END DO
  END DO
!  PRINT*,'in TSLEV,NZDAT',SIZE(NZDAT)
!  PAUSE
 DO K=1,SIZE(NZDAT)
   ZZLEV(K)=REAL(NZDAT(K))
 END DO 

 CALL EOS_COMPUTATION_STABILITY1(TLEV,SLEV,ZZLEV,RHO1,FL_EOS_OP,1,0)

! Yu-Chiao 20110928
  DO K=1,33  
    WRITE(FNO) ((TLEV(I,J,K),I=1,360),J=1,180),((SLEV(I,J,K),I=1,360),J=1,180)
  END DO


! =================================================
! Interpolate Levitus data to model horizontal grid
! =================================================

  CLOSE(FNO_fix)
  REWIND FNO
  WRITE(*,24)
24    FORMAT('interpolate Levitus data to model horizontal grid')
  DO K=1,33
    READ(FNO)                                                          &
    ((TLEV(I,J,1),I=1,360),J=1,180),((SLEV(I,J,1),I=1,360),J=1,180)
    WRITE(*,25) K
25    FORMAT('level',I4)
    XINC=DXMNUT/60.
    DO JJ=1,J0
    ! XLON is MODEL data point, measured in degrees east of Grenwich
    ! YLAT is MODEL data point, measured in degrees north of South Pole
    ! YDEG is measured from -90 at South Pole to 90 at North Pole
    ! TLEV(1,1) is at 0.5 deg E, -89.5 deg lat
    ! TLEV(360,180,) is at 359.5 deg E, 89.5 deg lat
    ! allowed model latitudes: -89.5.LT.YDEG.LT.89.5
    ! allowed model longitudes: 0.LE.XLON.LE.720
      YLAT=90.+YDEG(JJ)
      J=YLAT+0.5
      JP=J+1
      DO II=1,I0
        XLON=X0DEG+(II-1.5)*XINC
        IF (XLON.GT.360.) XLON=XLON-360.
      ! XPLUS is model data point distance east of 0.5 deg WEST(!)
      ! NOTE that XPLUS is ALWAYS positive due to XLON restrictions
        XPLUS=XLON+0.5
        I=XPLUS
        IP=I+1
      ! TMP1 is weighting of east point
        TMP1=XPLUS-I
        TMP2=1.-TMP1
      !
      ! Special treatments when 359.5 < XLON < 360 and 0 < XLON < 0.5
        IF (IP.EQ.361) IP=1
        IF (I.EQ.0) I=359
      !
        DNW=TLEV(I,JP,1)
        DNE=TLEV(IP,JP,1)
        DSW=TLEV(I,J,1)
        DSE=TLEV(IP,J,1)
        DN=TMP1*DNE+TMP2*DNW
        DS=TMP1*DSE+TMP2*DSW
        T(II,JJ,K)=(.5+YLAT-J)*DN+(.5-YLAT+J)*DS
        DNW=SLEV(I,JP,1)
        DNE=SLEV(IP,JP,1)
        DSW=SLEV(I,J,1)
        DSE=SLEV(IP,J,1)
        DN=TMP1*DNE+TMP2*DNW
        DS=TMP1*DSE+TMP2*DSW
        S(II,JJ,K)=(.5+YLAT-J)*DN+(.5-YLAT+J)*DS
      END DO
    END DO
  END DO

  DO K=1,33
    DO J=2,J1
      DO I=2,I1
         S(I,J,K)=MIN(MAX(S(I,J,K),20.),41.)
         T(I,J,K)=MIN(MAX(-1.8,T(I,J,K)),40.)
      END DO
    END DO
  END DO

! ===============================================
! Interpolate from Levitus levels to model levels
! ===============================================
  Z2=100.*NZDAT(1)
  K=1
  ZC=Z(2)
  N=1

  DO WHILE(.TRUE.)
    N=N+1
    Z1=Z2
    Z2=100.*NZDAT(N)

    IF (Z2.LT.ZC) CYCLE
322   CONTINUE
    C1=(Z2-ZC)/(Z2-Z1)
    C2=(ZC-Z1)/(Z2-Z1)

    WRITE(*,323) K,N
323   FORMAT('FOR DIECAST MODEL LEVEL',I3,', LOWER INPUT LEVEL=',I3)
    DO J=1,J0
      DO I=1,I0
        T1(I,J,K)=C1*T(I,J,N-1)+C2*T(I,J,N)
        S1(I,J,K)=C1*S(I,J,N-1)+C2*S(I,J,N)
      END DO
    END DO
    K=K+1
    IF (K.EQ.K0) EXIT
    ZC=Z(2*K)
    IF (ZC.GT.Z2) CYCLE
    GO TO 322

  END DO
  CLOSE(FNO)

END SUBROUTINE TSLEV
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
SUBROUTINE TAUXY(TXC,TYC,M)

  IMPLICIT NONE
  INTEGER:: M   !  MONTH
  INTEGER:: I,J,II,JJ,IP,JP
  REAL:: X,XINC,YLAT,XLON,XPLUS,TMP1,TMP2,DNW,DNE,DSW,DSE,DS3,DN,DS
  CHARACTER(10)::FMT
  CHARACTER(128)::TXC,TYC
  INTEGER :: FNO_TXC,FNO_TYC


  CALL OPEN_BIGENDIAN_OLD_TXT_FILE(TRIM(TXC),FNO_TXC)
  CALL OPEN_BIGENDIAN_OLD_TXT_FILE(TRIM(TYC),FNO_TYC)

! Hellerman
  READ(FNO_TXC,400) I,J,X,X,X,X,FMT
  READ(FNO_TYC,400) I,J,X,X,X,X,FMT
400   FORMAT(2I4/4(3X,F10.7),A10)

! ================
! Read Hellerman winds
! ================
  READ(FNO_TXC,415) TX
  READ(FNO_TYC,415) TY
415   FORMAT(5E15.7)
  WRITE(*,416)
416   FORMAT('read wind stress')

  CLOSE(FNO_TYC)
  CLOSE(FNO_TYC)


! ===================================================
! Interpolate Hellerman data to model horizontal grid
! ===================================================
! Increment X0DEG by one grid interval, and start at second latitude
! (no ghost zones in winds)
  XINC=DXMNUT/60.
  X0DEG=X0DEG+XINC
  DO JJ=2,J1
  ! ============================================================
  ! TX(I,1) is at latitude -89.5; TX(I,180) is at latitude 89.5
  ! TX(1,J) is at longitude 0.5; TX(360,J) is at longitude 359.5
  ! XLON is measured in degrees east from 0.5 (XLON=0)
  ! YDEG is measured in degrees north from -89.5 (YDEG=0)
  ! ============================================================
  ! XLON is measured in degrees east from 0 at Grenwich
  ! YDEG is measured from -90 at South Pole to 90 at North Pole
    YLAT=90.+YDEG(JJ)
    J=YLAT+0.5
    JP=J+1
    DO II=2,I1
      XLON=X0DEG+(II-1.5)*XINC
      IF (XLON.GT.360.) XLON=XLON-360.
    ! XPLUS is model data point distance east of 0.5 deg WEST(!)
    ! NOTE that XPLUS is ALWAYS positive due to XLON restrictions
      XPLUS=XLON+0.5
      I=XPLUS
      IP=I+1
    ! TMP1 is weighting of east point
      TMP1=XPLUS-I
      TMP2=1.-TMP1
    !
    ! Special treatments when 359.5 < XLON < 360 and 0 < XLON < 0.5
      IF (IP.EQ.361) IP=1
      IF (I.EQ.0) I=359
    !
      DNW=TX(I,JP)
      DNE=TX(IP,JP)
      DSW=TX(I,J)
      DSE=TX(IP,J)
      DN=TMP1*DNE+TMP2*DNW
      DS=TMP1*DSE+TMP2*DSW
      TAUX(II-1,JJ-1)=(.5+YLAT-J)*DN+(.5-YLAT+J)*DS
      DNW=TY(I,JP)
      DNE=TY(IP,JP)
      DSW=TY(I,J)
      DSE=TY(IP,J)
      DN=TMP1*DNE+TMP2*DNW
      DS=TMP1*DSE+TMP2*DSW
      TAUY(II-1,JJ-1)=(.5+YLAT-J)*DN+(.5-YLAT+J)*DS
    END DO
  END DO


END SUBROUTINE TAUXY
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
REAL(8) function r(t,s,p)

  IMPLICIT NONE
  REAL(8)::T,S,P,P02,P01

  p02=1747.4508988+t*(11.51588-0.046331033*t)&
     -s*(3.85429655+0.01353985*t)
! Yu-Chiao 20110915 why divide 10d0 ???
!  p01=p/10d0+5884.81703666+t*(39.803732+t*(-0.3191477&
!     +t*0.0004291133))+2.6126277*s
  p01=p+5884.81703666+t*(39.803732+t*(-0.3191477&
     +t*0.0004291133))+2.6126277*s
  r=p01/(p02+0.7028423*p01)
  return
END FUNCTION R
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
REAL(8) function theta(sal,ts,ps,pr)

  IMPLICIT NONE
  REAL(8)::SAL,TS,PS,PR,DELP,HAFP,DELT1,DELT2,DELT3,DELT4,T1,T2,T3,T4

  delp=pr-ps
  hafp=ps+.5*delp
  delt1=delp*gamma0(sal,ts,ps)
  t1=ts+.5*delt1
  delt2=delp*gamma0(sal,t1,hafp)
  t2=t1+.2928932*(delt2-delt1)
  delt3=delp*gamma0(sal,t2,hafp)
  t3=t2+1.707107*(delt3-0.5857864*delt2-0.1213203*delt1)
  delt4=delp*gamma0(sal,t3,pr)
  t4=t3+0.16666667*(delt4-6.828427*delt3+&
     4.828427*delt2+delt1)
  theta=t4
  return
END FUNCTION theta
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
REAL(8) function gamma0(ss,tt,p)

  IMPLICIT NONE
  REAL(8)::SS,TT,P,XX

  xx=ss-35.
  gamma0=0.35803e-4+0.18932e-5*xx+p*(0.18741e-7-&
         0.11351e-9*xx-0.46206e-12*p)&
        +tt*(0.85258e-5-0.42393e-7*xx+p*(-0.67795e-9+&
         0.27759e-11*xx+0.18676e-13*p)&
        +tt*(-0.68360e-7+p*(0.87330e-11-0.21687e-15*p)&
        +tt*(0.66228e-9-0.54481e-13*p)))
  return
END function gamma0
! ----------------------------------------------------------------------

SUBROUTINE ALLOC_STATUS(STR)
CHARACTER(LEN=*) :: STR

  IF (ISTAT/=0) THEN
    PRINT *,"Array:",TRIM(STR)
    STOP "Cannot allocate memory"
  ENDIF
END SUBROUTINE ALLOC_STATUS

SUBROUTINE NC_OUTPUT_KB
  INTEGER :: NCID
  INTEGER :: DIMIDS(2),LON_DIMID,LAT_DIMID
  INTEGER :: LON_VARID,LAT_VARID,KB_VARID
  REAL :: USCA,UADO
  CHARACTER(LEN=*),PARAMETER:: UNITS="units",SCA_NAME="scale_factor",ADO_NAME="add_offset"
  TEMPDIR=TRIM(WORKDIR)//'/OUTPUT_GLOBAL/KB.nc'
  CALL nc_check(nf90_create(TEMPDIR,NF90_CLOBBER,NCID))
  CALL nc_check(NF90_DEF_DIM(NCID,"longtitude",I0,LON_DIMID))
  CALL nc_check(NF90_DEF_DIM(NCID,"latitude",J0,LAT_DIMID))
  CALL nc_check( NF90_DEF_VAR(NCID,"longtitude",NF90_REAL,LON_DIMID,LON_VARID) )
  CALL nc_check( NF90_DEF_VAR(NCID,"latitude",NF90_REAL,LAT_DIMID,LAT_VARID) )
  CALL nc_check( NF90_PUT_ATT(NCID,LON_VARID,UNITS,"degrees_east") )
  CALL nc_check( NF90_PUT_ATT(NCID,LAT_VARID,UNITS,"degrees_north") )
  DIMIDS=(/LON_DIMID,LAT_DIMID/)
  CALL nc_check(NF90_DEF_VAR(NCID,"KB",NF90_SHORT,DIMIDS(1:2),KB_VARID))
  CALL nc_check( NF90_PUT_ATT(NCID,KB_VARID,UNITS,"KB") )
  CALL nc_check( NF90_ENDDEF(NCID) )
  CALL nc_check( NF90_PUT_VAR(NCID, LON_VARID, XDEG) )
  CALL nc_check( NF90_PUT_VAR(NCID, LAT_VARID, YDEG) )
  CALL nc_check( NF90_PUT_VAR(NCID, KB_VARID, KB) )
  CALL nc_check( NF90_CLOSE(NCID) )
END SUBROUTINE NC_OUTPUT_KB

SUBROUTINE NC_OUTPUT_T0
  INTEGER,PARAMETER :: INT2 = SELECTED_INT_KIND(4)
  INTEGER,PARAMETER :: INT2BITS = 16
  INTEGER ::NCID
  INTEGER ::DIMIDS(3),LON_DIMID,LAT_DIMID,DEP_DIMID
  INTEGER ::LON_VARID,LAT_VARID,DEP_VARID,T0_VARID
  REAL :: SCA,ADO
  REAL :: T0MAX,T0MIN
  INTRINSIC :: MAXVAL,MINVAL
  CHARACTER(LEN=*),PARAMETER:: UNITS="units",SCA_NAME="scale_factor",ADO_NAME="add_offset"
  REAL :: RTEMP
  INTEGER(KIND=INT2) :: ITEMP(I0,J0,K1)
  INTEGER :: I,J,K
!========================================
  T0MAX=MAXVAL(T1)
  T0MIN=MINVAL(T1)
  SCA=(T0MAX-T0MIN)/REAL(2**INT2BITS-1)
  ADO=T0MIN+2**(INT2BITS-1)*SCA
  DO K=1,K1
  DO J=1,J0
  DO I=1,I0
     RTEMP = (T1(I,J,K) - ADO)/SCA
     ITEMP(I,J,K) = INT(RTEMP,INT2)
  ENDDO
  ENDDO
  ENDDO
!=======================================
  TEMPDIR=TRIM(WORKDIR)//'/OUTPUT_GLOBAL/T0.nc'
  CALL nc_check(nf90_create(TEMPDIR,NF90_CLOBBER,NCID))
  CALL nc_check(NF90_DEF_DIM(NCID,"longtitude",I0,LON_DIMID))
  CALL nc_check(NF90_DEF_DIM(NCID,"latitude",J0,LAT_DIMID))
  CALL nc_check(NF90_DEF_DIM(NCID,"depth",K1,DEP_DIMID))
  CALL nc_check( NF90_DEF_VAR(NCID,"longtitude",NF90_REAL,LON_DIMID,LON_VARID) )
  CALL nc_check( NF90_DEF_VAR(NCID,"latitude",NF90_REAL,LAT_DIMID,LAT_VARID) )
  CALL nc_check( NF90_DEF_VAR(NCID,"depth",NF90_REAL,DEP_DIMID,DEP_VARID) )
  CALL nc_check( NF90_PUT_ATT(NCID,LON_VARID,UNITS,"degrees_east") )
  CALL nc_check( NF90_PUT_ATT(NCID,LAT_VARID,UNITS,"degrees_north") )
  CALL nc_check( NF90_PUT_ATT(NCID,DEP_VARID,UNITS,"meters") )
  DIMIDS=(/LON_DIMID,LAT_DIMID,DEP_DIMID/)
  CALL nc_check(NF90_DEF_VAR(NCID,"T1",NF90_SHORT,DIMIDS(1:3),T0_VARID))
  CALL nc_check( NF90_PUT_ATT(NCID,T0_VARID,UNITS,"T1") )
  CALL nc_check( NF90_PUT_ATT(NCID,T0_VARID,SCA_NAME,SCA) )
  CALL nc_check( NF90_PUT_ATT(NCID,T0_VARID,ADO_NAME,ADO) )
  CALL nc_check( NF90_ENDDEF(NCID) )
  CALL nc_check( NF90_PUT_VAR(NCID, LON_VARID, XDEG) )
  CALL nc_check( NF90_PUT_VAR(NCID, LAT_VARID, YDEG) )
  CALL nc_check( NF90_PUT_VAR(NCID, DEP_VARID, Z(2:K0+K1:2)) )
  CALL nc_check( NF90_PUT_VAR(NCID, T0_VARID, ITEMP))
  CALL nc_check( NF90_CLOSE(NCID) )
END SUBROUTINE NC_OUTPUT_T0

SUBROUTINE NC_OUTPUT_IN
  INTEGER :: NCID
  INTEGER :: DIMIDS(3),LON_DIMID,LAT_DIMID,DEP_DIMID
  INTEGER :: LON_VARID,LAT_VARID,DEP_VARID,IN_VARID
  REAL :: USCA,UADO
  CHARACTER(LEN=*),PARAMETER:: UNITS="units",SCA_NAME="scale_factor",ADO_NAME="add_offset"
  TEMPDIR=TRIM(WORKDIR)//'/OUTPUT_GLOBAL/IN.nc'
  CALL nc_check(nf90_create(TEMPDIR,NF90_CLOBBER,NCID))
  CALL nc_check(NF90_DEF_DIM(NCID,"longtitude",I0,LON_DIMID))
  CALL nc_check(NF90_DEF_DIM(NCID,"latitude",J0,LAT_DIMID))
  CALL nc_check(NF90_DEF_DIM(NCID,"depth",K1,DEP_DIMID))
  CALL nc_check( NF90_DEF_VAR(NCID,"longtitude",NF90_REAL,LON_DIMID,LON_VARID) )
  CALL nc_check( NF90_DEF_VAR(NCID,"latitude",NF90_REAL,LAT_DIMID,LAT_VARID) )
  CALL nc_check( NF90_DEF_VAR(NCID,"depth",NF90_REAL,DEP_DIMID,DEP_VARID) )
  CALL nc_check( NF90_PUT_ATT(NCID,LON_VARID,UNITS,"degrees_east") )
  CALL nc_check( NF90_PUT_ATT(NCID,LAT_VARID,UNITS,"degrees_north") )
  CALL nc_check( NF90_PUT_ATT(NCID,DEP_VARID,UNITS,"meters") )
  DIMIDS=(/LON_DIMID,LAT_DIMID,DEP_DIMID/)
  CALL nc_check(NF90_DEF_VAR(NCID,"IN",NF90_SHORT,DIMIDS,IN_VARID))
  CALL nc_check( NF90_PUT_ATT(NCID,IN_VARID,UNITS,"IN") )
  CALL nc_check( NF90_ENDDEF(NCID) )
  CALL nc_check( NF90_PUT_VAR(NCID, LON_VARID, XDEG) )
  CALL nc_check( NF90_PUT_VAR(NCID, LAT_VARID, YDEG) )
  CALL nc_check( NF90_PUT_VAR(NCID, DEP_VARID, Z(2:K0+K1:2)) )
  CALL nc_check( NF90_PUT_VAR(NCID, IN_VARID, IN) )
  CALL nc_check( NF90_CLOSE(NCID) )
END SUBROUTINE NC_OUTPUT_IN

subroutine nc_check(status)
integer, intent ( in) :: status

   if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop 2
    end if
end subroutine nc_check

END SUBROUTINE INDATA
