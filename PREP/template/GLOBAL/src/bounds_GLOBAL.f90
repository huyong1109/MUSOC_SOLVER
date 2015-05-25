SUBROUTINE BOUNDS
#include "namelist.in_GLOBAL"
  USE INIT_VAR_GLOBAL
  USE GRID_VAR_GLOBAL, ONLY: RHO,P0,PX,PY,U1,U2,V1,V2,S1,S2,T1,T2,P,ULF,VLF,SLF,TLF,EV,HV,F,TANPHI,SUMIN,VAR_ALLOCATE_B
  USE GRID_VAR_GLOBAL, ONLY: SSURF,TSURF,SSSP,SNSP,TSSP,TNSP,VAR_ALLOCATE_CLIMAT2
  USE GRID_VAR_GLOBAL, ONLY: VNOR,VSUD,UEST,UWST
  USE GRID_VAR_GLOBAL, ONLY: KB,IU0,IV0,IN,IU,IV,IW
  USE GRID_VAR_GLOBAL, ONLY: U,V,W,VAR_ALLOCATE_CGRID
  USE GRID_VAR_GLOBAL, ONLY: A,XDEG,YV,YVDEG,YDEG,CS,OCS,DX,ODX,DY,ODY,CSV,OCSV,DXV,ODXV,DYV,ODYV
  USE GRID_VAR_GLOBAL, ONLY: Z,F,ODZ,ODZW,IE,VAR_ALLOCATE_SEVP
  USE CONTROL_GLOBAL, ONLY:AROUT
!  USE GRID_VAR_GLOBAL, ONLY:AROUT
  USE GRID_VAR_GLOBAL, ONLY: SXYCLI,TXYCLI,VAR_ALLOCATE_XYMEANS_MAIN
  USE COUPLING
! ----------------------------------------------------------------------
  IMPLICIT NONE
  REAL::RHO_TEMP,DMIN,DMAX,DDDX,PRE_INF,PRE_INF2
  REAL(8)::TMPD,SLTD,PISD,p01,p02 !check now changed to real(8)
  REAL(8),ALLOCATABLE ::D1(:,:,:) !check now changed to real(8)
  CHARACTER(160)::TEMPDIR
  REAL :: VTMP,SVIN,DDDY,SV,AVG1,AVG2,AREAXY,TMP,TEM
  INTEGER :: FNO_bounds,FNO_annu,FNO_inflow
  INTEGER :: FL_SURF_ITP, FL_ANN_ITP, FL_SPG_ITP

  INTEGER :: status(4)
  INTEGER :: ILO,IHI,NM
  INTEGER :: I,J,K,N


  CALL VAR_ALLOCATE_CGRID
  CALL VAR_ALLOCATE_SEVP
  CALL VAR_ALLOCATE_CLIMAT2
  CALL VAR_ALLOCATE_XYMEANS_MAIN
  allocate (D1(I0,J0,K1),STAT=status(3))
  if (status(3)/=0) STOP "Cannot allocate memory"

  PRINT *, "BOUNDS is CALLED!"

  DT=172800./DAODT
  ODT=1./DT
  TEMPDIR=TRIM(WORKDIR)//'/INPUT_GLOBAL/boundaries'
  CALL OPEN_BIGENDIAN_NEW_BIN_FILE(TRIM(TEMPDIR),FNO_bounds) !original is 62
#ifdef FLAG_TS_T
#ifdef FLAG_TS_S
  IF(FL_RL_TS.EQ.1) THEN
    TEMPDIR=TRIM(WORKDIR)//'/TEMP_GLOBAL/annualevitus'
    CALL OPEN_BIGENDIAN_OLD_BIN_FILE(TRIM(TEMPDIR),FNO_annu) !original is 70
  END IF
#endif
#endif 
  IF (TRIM(CASE_NAME) == 'ILE2D') THEN
    VTMP=0.
    IF (LOPENS==2) THEN	
      ! specify open "southern" inflow
      DO K=1,K1
        DO I=2,I1
          V(I,1,K)=VINFLOW
          V1(I,1,K)=VTMP
        END DO
      END DO	
    END IF
  
    ! match initial outflow to inflow
    ! outflow evolves in model
    IF (LOPENN==2) THEN
      DO K=1,K1
        DO I=2,I1				
          V(I,J1,K)=VINFLOW
          V1(I,J0,K)=VTMP
        END DO
      END DO	
    END IF

    ! WE MAY USE USE boundaries INSTEAD OF inflow TO SAVE THE INFLOW BOUNDARY
    IF (LOPENS==2 .or. LOPENN==2) THEN		
      WRITE(FNO_bounds) V,V1,V1,V1
    END IF  
  END IF

  IF (TRIM(CASE_NAME) == 'DTRAC') THEN
    ! CALCULATE DENSITIES
    DMIN=1.E10		!ZWM	assigned here but later changes, only for output, not used in calculation
    DMAX=-1.E10
    DO K=1,K1
      DO J=1,J0
        DO I=1,I0
          RHO_TEMP=.0002*(20.-T1(I,J,K))
          DMIN=MIN(DMIN,RHO_TEMP)
          DMAX=MAX(DMAX,RHO_TEMP)
          D1(I,J,K)=RHO_TEMP
        END DO
      END DO
    END DO
    !WRITE(17,327) DMIN,DMAX
    WRITE(*,327) DMIN,DMAX
    327  FORMAT('DMIN,DMAX=',2(1X,F9.4))

    ! DOME problem
    ! NORTH
    SVIN=0.		
    N=6.*K1/22.
    ILO=8.*I2/11
    IHI=9.*I2/11
    DO K=2,N		!ZWM edited below for easy readability
      DO I=ILO-1,IHI+1
        DDDX=.25*ODXV(J1)*(D1(I+1,J0,K)-D1(I-1,J0,K)+D1(I+1,J0,K-1)-D1(I-1,J0,K-1))
        VNOR(I,K)=IN(I,J1,K)*(VNOR(I,K-1)+G*DDDX/(F(J2)*ODZW(K)))  !bc? anything to do with bc? user spec
        SVIN=SVIN-VNOR(I,K)*DXV(J1)/ODZ(K)
      END DO
    END DO

    DO K=1,N+1		!ZWM edited for easy readability
      !WRITE(17,342) K,(INT(-10.*VNOR(I,K)),I=ILO,IHI)
      WRITE(*,342) K,(INT(-10.*VNOR(I,K)),I=ILO,IHI)
    END DO
    342  FORMAT('K=',I2/(21I4))
    SVIN=1.E-12*SVIN
    !WRITE(17,344) SVIN
    WRITE(*,344) SVIN
    344  FORMAT('northern inflow velocity (mm/sec)'/   &
    'total northern inflow=',F6.2, 'Sv')

    ! EAST
    SV=0.
    DO J=2,J1		!ZWM edited for easy readability
      UEST(J,1)=0.
    END DO 
    DO K=2,K1
      DO J=3,J2
        DDDY=.25*ODY(J)*(D1(I0,J+1,K)-D1(I0,J-1,K)+D1(I0,J+1,K-1)-D1(I0,J-1,K-1))
        UEST(J,K)=IN(I1,J,K)*(UEST(J,K-1)-G*DDDY/(F(J-1)*ODZW(K)))
        SV=SV-UEST(J,K)*DY(J)/ODZ(K)
      END DO
      UEST(2,K)=UEST(3,K)
      UEST(J1,K)=UEST(J2,K)
      SV=SV-(UEST(2,K)*DY(2)+UEST(J1,K)*DY(J1))/ODZ(K)
    END DO

    SV=1.E-12*SV
    SVIN=SVIN+SV
    !WRITE(17,355) SV,SVIN
    WRITE(*,355) SV,SVIN
    355  FORMAT(F4.1,' Sv eastern inflow added for a total of ',F4.1,   &
    'Sv inflow.')
    WRITE(FNO_bounds) UEST,UEST,UWST,UWST,VNOR,VNOR,VSUD,VSUD
  END IF

#ifdef FLAG_TS_T
#ifdef FLAG_TS_S
  IF (FL_RL_TS>=1) THEN
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! THERE ARE SEVEN STEPS IN THE FOLLOWING                              !
    ! 1) READ IN MONTHLY CLIMATOLOGICAL S & T                             !
    ! 2) ASSIGN THE ABOVE VALUE TO SSURF & TSURF                          !
    ! 3) CALCULATE ANNUAL MEAN CLIMATOLOGICAL SXYCLI & TXYCLI             !
    ! 4) ASSIGN SPONGE LAYER CLIMATOLOGY                                  !
    ! 5) INTERPOLATE SSURF & TSURF           (BY SETTING FL_SURF_ITP = 1) !
    ! 6) INTERPOLATE SXYCLI & TXYCLI         (BY SETTING FL_ANN_ITP  = 1) !
    ! 7) INTERPOLATE SSSP, SNSP, TSSP & TNSP (BY SETTING FL_SPG_ITP  = 1) !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FL_SURF_ITP= 1
      FL_ANN_ITP = 1
      FL_SPG_ITP = 1
    
238   format('month',i3,' in bounds')
    DO NM=1,12
      write(*,238) NM
    ! 1) READ IN MONTHLY CLIMATOLOGICAL S & T                 
      READ(FNO_annu) S1,T1
      DO K=1,K1
        DO J=2,J1
          DO I=2,I1
            S1(I,J,K)=MIN(MAX(S1(I,J,K),20.),41.)
            T1(I,J,K)=MIN(MAX(-1.8,T1(I,J,K)),40.)
          END DO
        END DO
      END DO    
    ! 2) ASSIGN THE ABOVE VALUE TO SSURF & TSURF              
      DO J=2,J1
        DO I=2,I1
          TSURF(I-1,J-1,NM)=T1(I,J,1)
          SSURF(I-1,J-1,NM)=S1(I,J,1)
        END DO
      END DO
      
    ! 3) CALCULATE ANNUAL MEAN CLIMATOLOGICAL SXYCLI & TXYCLI 
    ! CALCULATIONS OF MAX AND MIN IS NEGLECTED
      DO K=1,K1
        AVG1=0.
        AVG2=0.
        AREAXY=0.
        DO J=2,J1
          TMP=DX(J)*DY(J)
          DO I=2,I1
            AREAXY=AREAXY+IN(I,J,K)*TMP
            AVG1=AVG1+IN(I,J,K)*TMP*T1(I,J,K)
            AVG2=AVG2+IN(I,J,K)*TMP*S1(I,J,K)
          END DO
        END DO
        TXYCLI(K,NM)=AVG1/AREAXY
        SXYCLI(K,NM)=AVG2/AREAXY
      END DO
      
    ! 4) ASSIGN SPONGE LAYER CLIMATOLOGY                   
      DO K=1,K1
        DO J=1,10
          DO I=2,I1
            TSSP(I-1,J,K,NM)=T1(I,J+1,K)
            TNSP(I-1,J,K,NM)=T1(I,J0-J,K)
            SSSP(I-1,J,K,NM)=S1(I,J+1,K)
            SNSP(I-1,J,K,NM)=S1(I,J0-J,K)
          END DO
        END DO
      END DO
    END DO

  ! 5) INTERPOLATE SSURF & TSURF                            
    TMP=1./3.
    TEM=2./3.
    IF (FL_SURF_ITP == 1) THEN 
      DO N=1,10,3
        NM=N+3
        IF (NM.EQ.13) NM=1
        DO J=1,J2
          DO I=1,I2
            TSURF(I,J,N+1)=TSURF(I,J,N)+TMP*(TSURF(I,J,NM)-TSURF(I,J,N))
            TSURF(I,J,N+2)=TSURF(I,J,N)+TEM*(TSURF(I,J,NM)-TSURF(I,J,N))
            SSURF(I,J,N+1)=SSURF(I,J,N)+TMP*(SSURF(I,J,NM)-SSURF(I,J,N))
            SSURF(I,J,N+2)=SSURF(I,J,N)+TEM*(SSURF(I,J,NM)-SSURF(I,J,N))
          END DO
        END DO
      END DO
      !PRINT*,'T0',SUM(TSSP),TMP,TEM
    END IF  

  ! 6) INTERPOLATE SXYCLI & TXYCLI                          
    IF (FL_ANN_ITP == 1) THEN 
      DO N=1,10,3
        NM=N+3
        IF (NM.EQ.13) NM=1
        DO K=1,K1
          TXYCLI(K,N+1)=TXYCLI(K,N)+TMP*(TXYCLI(K,NM)-TXYCLI(K,N))
          TXYCLI(K,N+2)=TXYCLI(K,N)+TEM*(TXYCLI(K,NM)-TXYCLI(K,N))
          SXYCLI(K,N+1)=SXYCLI(K,N)+TMP*(SXYCLI(K,NM)-SXYCLI(K,N))
          SXYCLI(K,N+2)=SXYCLI(K,N)+TEM*(SXYCLI(K,NM)-SXYCLI(K,N))
        END DO
      END DO
    END IF  

  ! 7) INTERPOLATE SSSP, SNSP, TSSP & TNSP
    IF (FL_SPG_ITP == 1) THEN 
      DO N=1,10,3
        NM=N+3
        IF (NM.EQ.13) NM=1
        DO K=1,K1
          DO J=1,10
            DO I=1,I2
              TSSP(I,J,K,N+1)=TSSP(I,J,K,N)+TMP*(TSSP(I,J,K,NM)-TSSP(I,J,K,N))
              TSSP(I,J,K,N+2)=TSSP(I,J,K,N)+TEM*(TSSP(I,J,K,NM)-TSSP(I,J,K,N))
              TNSP(I,J,K,N+1)=TNSP(I,J,K,N)+TMP*(TNSP(I,J,K,NM)-TNSP(I,J,K,N))
              TNSP(I,J,K,N+2)=TNSP(I,J,K,N)+TEM*(TNSP(I,J,K,NM)-TNSP(I,J,K,N))
              SSSP(I,J,K,N+1)=SSSP(I,J,K,N)+TMP*(SSSP(I,J,K,NM)-SSSP(I,J,K,N))
              SSSP(I,J,K,N+2)=SSSP(I,J,K,N)+TEM*(SSSP(I,J,K,NM)-SSSP(I,J,K,N))
              SNSP(I,J,K,N+1)=SNSP(I,J,K,N)+TMP*(SNSP(I,J,K,NM)-SNSP(I,J,K,N))
              SNSP(I,J,K,N+2)=SNSP(I,J,K,N)+TEM*(SNSP(I,J,K,NM)-SNSP(I,J,K,N))
            END DO
          END DO
        END DO
      END DO
    END IF
   
    print *,"write to bounds" 
    WRITE(FNO_bounds) SSURF,TSURF,SSSP,SNSP,TSSP,TNSP
    WRITE(FNO_bounds) SXYCLI,TXYCLI
    CLOSE(FNO_bounds)
    CLOSE(FNO_annu)
  END IF
#endif
#endif
END SUBROUTINE BOUNDS
