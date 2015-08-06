SUBROUTINE FS_GLOBAL
! ----------------------------------------------------------------------
! MAIN DYNAMICS PROGRAM
! FS is a wind stress and buoyancy driven free stream model, with
! rigid top, variable DZ, and bottom topography.
! FS resolves bottom features using a staircase approximation.
! EXACT conservation is satisfied.
! Bottom fluxes are specified as SOURCES when coupled to TS (thin-shell
! submodel). W is positive downward. W(I,J,1) is at free surface.
! W(I,J,KB(I,J)+1), where 0.LE.KB.LE.K1, is at bottom
! ----------------------------------------------------------------------
#include "namelist.in_GLOBAL"
USE OCN_PARA_GLOBAL
USE GRID_VAR_GLOBAL, ONLY: RINV,RINV1,DUM0,DUM1,DUM2,XX,H,X,S,AL,AB,AC,AR,AT,SRC,CL,CB,CC,OCC,CR,CT,IE
USE GRID_VAR_GLOBAL, ONLY: RHO,Mcrho,McTP,McS,McTPTP,McSTP,P0,PX,PY,U1,U2,V1,V2,S1,S2,T1,T2,P,ULF,VLF,SLF,TLF,EV,HV,F,TANPHI,SUMIN
USE GRID_VAR_GLOBAL, ONLY: U,V,W
USE GRID_VAR_GLOBAL, ONLY: PBAR,PVAR,XBAR,UCLI,VCLI,RMSV,SBAR,TBAR
USE GRID_VAR_GLOBAL, ONLY: Z,ODZ,ODZW
USE GRID_VAR_GLOBAL, ONLY: A,XDEG,YV,YVDEG,YDEG,CS,OCS,DX,ODX,DY,ODY,CSV,OCSV,DXV,ODXV,DYV,ODYV
USE GRID_VAR_GLOBAL, ONLY: KB,IU0,IV0,IN,IU,IV,IW
USE CLIMAT_GLOBAL
USE GRID_VAR_GLOBAL, ONLY: TNUDGE,QAVG,WAVG,SCLI,TCLI,SSURF,TSURF,SSSP,SNSP,TSSP,TNSP
USE CONTROL_GLOBAL
USE SCA_GLOBAL
USE GRID_VAR_GLOBAL, ONLY: SCR
USE GRID_VAR_GLOBAL, ONLY: SXY,TXY,SXYCLI,TXYCLI
USE PHYSICS_GLOBAL
USE SURF_GLOBAL, ONLY:QSURF,QVOLUME,SALSURF
USE WINDS_GLOBAL,ONLY:WIND
USE GRID_VAR_GLOBAL, ONLY: DMX,DMY
USE GRID_VAR_GLOBAL, ONLY: XP,YP,DHX,DHY
USE GLOBALSPLT, ONLY: SC,TC,UC,VC,nexp
#ifdef FLAG_TRACER
USE GRID_VAR_GLOBAL, ONLY: C1,C2,CLF
#endif
USE GRID_VAR_GLOBAL, ONLY: AROUT,IFLT
use GRID_VAR_GLOBAL, only: nudge
USE DENSITY_EOS_GLOBAL

use grid_var_GLOBAL, only: qdot,qdot2,rain,evapo
USE GRID_VAR_GLOBAL, ONLY:TAUX,TAUY !REAL,DIMENSION(I2,J2,12)::TAUX,TAUY
USE GRID_VAR_GLOBAL, ONLY:TAUX1,TAUY1
use getm_meteo

REAL(8)::TMPD,SLTD,PISD,p01,p02,ccrho,ccTP,ccS,ccTPTP,ccSTP
REAL,DIMENSION(I1,J2)::UX,VX,SX,TX
REAL,DIMENSION(I2,J1)::UY,VY,SY,TY
REAL,DIMENSION(I2,J2,2)::UZ,VZ,SZ,TZ
REAL,DIMENSION(I1,J2)::CX
REAL,DIMENSION(I2,J1)::CY
REAL,DIMENSION(I2,J2,2)::CZ
INTEGER::UP_BOUND,M
REAL :: PSM,NSM,TMP_U,TMP_V,TMP_T,TMP_S,TMP_C
REAL(8) :: WFACE(I2,J2,K0)
REAL(8) :: UFACE(I1,J2)
REAL(8) :: VFACE(I2,J1)
logical::space
real :: seavol,rivervol,storet,stores

  UFACE=0.;VFACE=0.;WFACE=0.;TMPBC=0.
  TMPD=0.;SLTD=0.;PISD=0.;p01=0.;p02=0.
  UX=0.;VX=0.;SX=0.;TX=0.;UY=0.;VY=0.;SY=0.;TY=0.;UZ=0.;VZ=0.;SZ=0.;TZ=0.
  O12=1./12.; O16=1./16.; O24=1./24.
  

  ! =============================================================================== !
  ! THIS CODE INCLUDES THE FOLLOWING STEPS:                                         !
  !  1. PHYSICS: VERTICAL MIXING COEFF., MODIFIED PACANOWSKI & PHILANDER 82         !
  !  2. WEAK NUDGING TO CLIMATOLOGY                                                 !
  !  3. HYDROSTAIC EQUATION WITH STATIC BALANCE SUBTRACTED OUT                      !
  !     A. CALCULATE RHO                                                            !
  !     B. CALCULATE FACE AVERAGED P @ BOTTOM OF EACH LAYER                         !
  !     C. CALCULATE VOLUME AVERAGED P (4TH-ORDER ACCURATE)                         !
  !  3bis. ADD RIVERS                                                               !
  !  4. CALCULATE ALL INTERNAL FLUXES AND SUBSTITUTE INTO CONSERVATION LAWS         !
  !     A. HORIZONTAL PRESSURE GRADIENT, 4TH-ORDER-ACCURATE                         !
  !     B. VERTICAL FLUXES                                                          !
  !     C. LONGITUDINAL FLUXES                                                      !
  !     D. LATITUDINAL FLUXES                                                       !
  !  5. MOMENTUM EQUATIONS & TRANSPORT EQUATIONS (TMP U2,V2,T2,S2,C2)               !
  !  6. OPEN BOUNDARY CONDITIONS                                                    !
  !  7. SURFACE EFFECTS                                                             !
  !     A. HEAT FLUX     B. WIND STRESS                                             !
  !     C. EVAPORATION & PRECIPITATION (NO SALT SOURCES)                            !
  !     D. SURFACE NUDGING                                                          !
  !  8. CLIMATOLOGCAL SPONGE LAYER RESTORING                                        !
  !  9. BOTTOM STRESS                                                               !
  ! 10. TRAPEZOIDAL CORIOLIS                                                        !
  ! 11. SET PERIODIC BOUNDARY CONDITIONS IF NEEDED                                  !
  ! 12. INTERPOLATE "A" GRID QUANTITIES TO "C" GRID LOCATIONS, 4TH-ORDER-ACCURATE   !
  ! 13. INFLOW & OUTFLOW CONDITIONS, CASE DEPENDENT                                 !
  ! 14. ITERATE EVP SOLVER TO GET EXACTLY NON-DIVERGENT BAROTROPIC MODE (P0, U, V)  !
  ! 15. CALCULATE W BY CONTINUITY EQUATION                                          !
  ! 16. INTERPOLATE INCOMPRESSIBILITY INDUCED CHANGES BACK TO CELL CENTERS (U2&V2)  !
  ! 17. CHECK INCOMPRESSIBILITY                                                     !
  ! 18. ELIMINATE ARBITRARY PRESSURE CONSTANT SO THAT MEAN SURFACE PRESSURE IS ZERO !
  ! 19. SET PERIODIC BOUNDARY CONDITIONS IF NEEDED                                  !
  ! 20. UPDATE USING FLTW METHOD                                                    !
  ! 21. NUDGE TO CLIMATOLOGY DATA                                                   !
  ! 22. RELAX SWAMP LAYER & MASKED WATER TO INTERIOR TO REDUCE SEEPAGE EFFECTS      !    
  ! 23. BIHARMONIC FILTERS                                                          !
  ! 24. LIMITER APPROACH FOR TRACER CONCENTRATION                                   !
  ! =============================================================================== !
  ALPHA=0.1
  ! ================================================================ !
  ! 0. DETERMINE NUDGING COEFFICIENTS (CURRENT MONTH AND FNEW, FOLD) !
  ! ================================================================ !
  if (FL_FNEWFOLD_CALC.eq.1) then
    NLD=DAYS/30. ! NLD IS AN INTEGER
    FNEW=DAYS/30.-NLD
    NLD=MOD(NLD,12)+1
  else
    nld=month
    fnew=day/30.0 ! this is a small approximation, as current month doesn't necessarily have 30 days
  endif
  FOLD=1.-FNEW
  NEW=NLD+1
  IF (NEW.EQ.13) NEW=1

! PUT IN THE VERY BEGINNING
  ! ======================================================================== !
  ! 1. PHYSICS: VERTICAL MIXING COEFF., MODIFIED PACANOWSKI & PHILANDER 82   !
  ! ======================================================================== !
  IF (MOD(ITF,ITF_PHYS).EQ.1) CALL TURBULENCE(TURB_PP82)   !E.G. ITF_PHY = 4, SPECIFIED IN NAMELIST
#ifdef _DEBUG_
  PRINT*,'Pass Step 1'
#endif
  ! ================================== !
  ! 2. WEAK NUDGING TO CLIMATOLOGY     !
  ! ================================== !
  ! NUDGE HORIZONTAL T,S MEAN TOWARD CLIMATOLOGY HORIZONTAL MEANS
  IF (FL_NUDGE_HMEAN_T .GE. 1) THEN !E.G. GLOBAL
      DO K=2,K1
        TXY(K)=0.
        DO J=2,J1
          TMP=DX(J)*DY(J)
          DO I=2,I1
            TXY(K)=TXY(K)+IN(I,J,K)*TMP*T1(I,J,K)
          END DO
        END DO
        IF (A(K).NE.0.) TXY(K)=TXY(K)/A(K)
      END DO
      TMP=1./(NUDGING_HMEAN*DAODT)
      DO  K=2,K1
        EPS=TMP
        SCR(1,2,K)=FNEW*TXYCLI(K,NEW)+FOLD*TXYCLI(K,NLD)
        TMP2=EPS*(SCR(1,2,K)-TXY(K))
        DO J=2,J1
          DO I=2,I1
            if (a(k).ne.0.)   T1(I,J,K)=T1(I,J,K)+IN(I,J,K)*TMP2
          END DO
        END DO
      END DO
  ENDIF
  IF (FL_NUDGE_HMEAN_S .GE. 1) THEN !E.G. GLOBAL
      DO K=2,K1
        SXY(K)=0.
        DO J=2,J1
          TMP=DX(J)*DY(J)
          DO I=2,I1
            SXY(K)=SXY(K)+IN(I,J,K)*TMP*S1(I,J,K)
          END DO
        END DO
        IF (A(K).NE.0.) SXY(K)=SXY(K)/A(K)
      END DO
      TMP=1./(NUDGING_HMEAN*DAODT)
      DO  K=2,K1
        EPS=TMP
        SCR(1,2,K)=FNEW*SXYCLI(K,NEW)+FOLD*SXYCLI(K,NLD)
        TMP2=EPS*(SCR(1,2,K)-SXY(K))
        DO J=2,J1
          DO I=2,I1
            if (a(k).ne.0.)   S1(I,J,K)=S1(I,J,K)+IN(I,J,K)*TMP2
          END DO
        END DO
      END DO
  ENDIF
#ifdef _DEBUG_
  PRINT*,'Pass Step 2'
#endif

  ! =========================================================== !
  ! 3. HYDROSTAIC EQUATION WITH STATIC BALANCE SUBTRACTED OUT   !
  ! =========================================================== !
  ! A. CALCULATE RHO

!Yu-Chiao 20110928
  CALL EOS_COMPUTE(T2,S2,Z,RHO,FL_EOS_OP)
  ! B. CALCULATE FACE AVERAGED P @ BOTTOM OF EACH LAYER
  DO J=2,J1; JJ=J-1;        !(J:CELL, JJ:SIZE OF WFACE, EASY TO READ)
    DO I=2,I1; II=I-1;      !(I:CELL, II:SIZE OF WFACE, EASY TO READ)
      WFACE(II,JJ,1)=P0(I,J)
      DO K=1,K1             !(K: LAYER, TOP FACE: K, BOTTOM FACE: K+1)
        WFACE(II,JJ,K+1)=WFACE(II,JJ,K)+RHO(I,J,K)*G/ODZ(K)
      END DO
    END DO
  END DO
  ! C. CALCULATE VOLUME AVERAGED P (4TH-ORDER ACCURATE)
  DO J=2,J1; JJ=J-1;        !(J:CELL, JJ:SIZE OF WFACE, EASY TO READ)
    DO I=2,I1; II=I-1;      !(I:CELL, II:SIZE OF WFACE, EASY TO READ)
      DO K=2,K2             !(K:INTERIOR LAYERS, TOP FACE: K, BOTTOM FACE: K+1)
        P(I,J,K)=(WFACE(II,JJ,K)+WFACE(II,JJ,K+1))*12.+    &
                 (WFACE(II,JJ,K)+WFACE(II,JJ,K+1)-WFACE(II,JJ,K-1)-WFACE(II,JJ,K+2))
        P(I,J,K)=P(I,J,K)*O24                 
      END DO
      ! 2ND-ORDER-ACCURATE IN THE TOP & BOTTOM LAYERS
      P(I,J,1) =.5*(WFACE(II,JJ,1) +WFACE(II,JJ,2))
      P(I,J,K1)=.5*(WFACE(II,JJ,K1)+WFACE(II,JJ,K0))
    END DO
!!    IF (LOPENE == 1) THEN   ! LOPENW = LOPENE =1 ! PERIODIC B.C.
#ifdef FLAG_PERIODIC_WE
      DO K=1,K1
        P(1,J,K) =P(I1,J,K)
        P(I0,J,K)=P(2,J,K)
      END DO
#endif
!!    END IF  
  END DO
#ifdef _DEBUG_
  PRINT*,'Pass Step 3'
#endif

  ! ================================== !
  ! 3bis. ADD RIVERS: INFLOW, T and S  !
  ! ================================== !
  !  use month as an index for riverF(k,month) riverT(k,month) riverS(k,month)
  do k=1,rivercount
   seavol=1/ocs(riverj(k))/odx(riverj(k))/ody(riverj(k))/odz(1) ! cm3
   rivervol=riverFlow(k,month)*1000000*2*86400/daodt ! idem
   !write(*,*) "putting riverdata : ",riveri(k),riverj(k),riverD(k),month,riverFlow(k,month)
   !write(*,*) "seavol, river volume=",seavol,rivervol
   if (riverD(k).eq.'W') then
      U(riveri(k),riverj(k),1)=-riverFlow(k,month)*1000000*odz(1)*ocs(riverj(k))*ody(riverj(k))
      U1(riveri(k),riverj(k),1)=-riverFlow(k,month)*1000000*odz(1)*ocs(riverj(k))*ody(riverj(k))
      storet=T2(riveri(k),riverj(k),1); stores=S2(riveri(k),riverj(k),1)
      T2(riveri(k),riverj(k),1)=(seavol*T2(riveri(k),riverj(k),1) + (rivervol*riverT(k,month))) / (seavol+rivervol)
      S2(riveri(k),riverj(k),1)=(seavol*S2(riveri(k),riverj(k),1) + (rivervol*riverS(k,month))) / (seavol+rivervol)
      T1(riveri(k),riverj(k),1)=(seavol*T1(riveri(k),riverj(k),1) + (rivervol*riverT(k,month))) / (seavol+rivervol)
      S1(riveri(k),riverj(k),1)=(seavol*S1(riveri(k),riverj(k),1) + (rivervol*riverS(k,month))) / (seavol+rivervol)
      !write(*,*) "river: ",riverlist(k)," T,S: ",storet,stores," --> ",T2(riveri(k),riverj(k),1), S2(riveri(k),riverj(k),1)
   elseif (riverD(k).eq.'E') then
      U(riveri(k)-1,riverj(k),1)=riverFlow(k,month)*1000000*odz(1)*ocs(riverj(k))*ody(riverj(k))
      U1(riveri(k)-1,riverj(k),1)=riverFlow(k,month)*1000000*odz(1)*ocs(riverj(k))*ody(riverj(k))
      storet=T2(riveri(k),riverj(k),1); stores=S2(riveri(k),riverj(k),1)
      T2(riveri(k),riverj(k),1)=(seavol*T2(riveri(k),riverj(k),1) + (rivervol*riverT(k,month))) / (seavol+rivervol)
      S2(riveri(k),riverj(k),1)=(seavol*S2(riveri(k),riverj(k),1) + (rivervol*riverS(k,month))) / (seavol+rivervol)
      T1(riveri(k),riverj(k),1)=(seavol*T1(riveri(k),riverj(k),1) + (rivervol*riverT(k,month))) / (seavol+rivervol)
      S1(riveri(k),riverj(k),1)=(seavol*S1(riveri(k),riverj(k),1) + (rivervol*riverS(k,month))) / (seavol+rivervol)
   elseif (riverD(k).eq.'N') then
      V(riveri(k),riverj(k)-1,1)=riverFlow(k,month)*1000000*odz(1)*odx(riverj(k))
      V1(riveri(k),riverj(k)-1,1)=riverFlow(k,month)*1000000*odz(1)*odx(riverj(k))
      storet=T2(riveri(k),riverj(k),1); stores=S2(riveri(k),riverj(k),1)
      T2(riveri(k),riverj(k),1)=(seavol*T2(riveri(k),riverj(k),1) + (rivervol*riverT(k,month))) / (seavol+rivervol)
      S2(riveri(k),riverj(k),1)=(seavol*S2(riveri(k),riverj(k),1) + (rivervol*riverS(k,month))) / (seavol+rivervol)
      T1(riveri(k),riverj(k),1)=(seavol*T1(riveri(k),riverj(k),1) + (rivervol*riverT(k,month))) / (seavol+rivervol)
      S1(riveri(k),riverj(k),1)=(seavol*S1(riveri(k),riverj(k),1) + (rivervol*riverS(k,month))) / (seavol+rivervol)
   elseif (riverD(k).eq.'S') then
      V(riveri(k),riverj(k),1)=-riverFlow(k,month)*1000000*odz(1)*odx(riverj(k))
      V1(riveri(k),riverj(k),1)=-riverFlow(k,month)*1000000*odz(1)*odx(riverj(k))
      storet=T2(riveri(k),riverj(k),1); stores=S2(riveri(k),riverj(k),1)
      T2(riveri(k),riverj(k),1)=(seavol*T2(riveri(k),riverj(k),1) + (rivervol*riverT(k,month))) / (seavol+rivervol)
      S2(riveri(k),riverj(k),1)=(seavol*S2(riveri(k),riverj(k),1) + (rivervol*riverS(k,month))) / (seavol+rivervol)
      T1(riveri(k),riverj(k),1)=(seavol*T1(riveri(k),riverj(k),1) + (rivervol*riverT(k,month))) / (seavol+rivervol)
      S1(riveri(k),riverj(k),1)=(seavol*S1(riveri(k),riverj(k),1) + (rivervol*riverS(k,month))) / (seavol+rivervol)
   end if
   IF ( (ITF_GENERAL .LT. TOPTS) .OR. (MOD(ITF_GENERAL,TOPTS) .EQ. 0) .OR. (ITFTO .LE. NITFTO) ) &
     write(*,*) "river: ",riverlist(k)," T,S: ",storet,stores," --> ",T2(riveri(k),riverj(k),1), S2(riveri(k),riverj(k),1)
  end do


  ! ======================================================================== !
  ! 4. CALCULATE ALL INTERNAL FLUXES AND SUBSTITUTE INTO CONSERVATION LAWS   !
  ! ======================================================================== !
  ! USE MOVING WINDOW FOR SMALLER STORAGE & EFFICIENCY I=LB; LB=LT; LT=I (SWITCH)
  LB=1; LT=2 
  DO K=1,K1; KK=K+1;        !(K:LAYER, KK:BOTTOM FACE OF CURRENT LAYER)
    ! A. HORIZONTAL PRESSURE GRADIENT, 4TH-ORDER-ACCURATE
    ! A-1. LONGITUDINAL DIRECTION
    DO J=2,J1; JJ=J-1;      !(J:CELL, JJ:SIZE OF UFACE, EASY TO READ)
      DO II=2,I2; I=II;     !(I:CELL, II:INTERIOR FACE, EASY TO READ)
        IF (IN(I,J,K)*IN(I+1,J,K) == 1) THEN
          TMP=IN(I-1,J,K)*IN(I,J,K)*IN(I+1,J,K)*IN(I+2,J,K)
          UFACE(II,JJ)=6.*(P(I,J,K)+P(I+1,J,K))   &
                         +(P(I,J,K)+P(I+1,J,K)-P(I-1,J,K)-P(I+2,J,K))*TMP
        ELSE
          UFACE(II,JJ)=12.*(IN(I,J,K)*P(I,J,K)+IN(I+1,J,K)*P(I+1,J,K))
        END IF             
      END DO
      ! BOUNDARY FACES (II=1:WBC, II=I1:EBC)
!!      IF (LOPENE == 1) THEN ! LOPENW = LOPENE = 1 ! PERIODIC B.C.
#ifdef FLAG_PERIODIC_WE
        TMP=IN(I2,J,K)*IN(I1,J,K)*IN(2,J,K)*IN(3,J,K)
        UFACE(1,JJ) =6.*(P(I1,J,K)+P(2,J,K))   &
                       +(P(I1,J,K)+P(2,J,K)-P(I2,J,K)-P(3,J,K))*TMP
        UFACE(I1,JJ)=UFACE(1,JJ)
!!      ELSE 
#else
        UFACE(1,JJ) =12.*P(2,J,K); UFACE(I1,JJ)=12.*P(I1,J,K)
#endif
!!      END IF
      DO I=2,I1; II=I;      !(I:CELL, II:FACE, EASY TO READ)
        PX(I,J)=O12*(UFACE(II,JJ)-UFACE(II-1,JJ))*IN(I,J,K)
      END DO 
    END DO
    ! A-2. LATITUDINAL DIRECTION
    DO I=2,I1; II=I-1;      !(I:CELL, II:SIZE OF VFACE, EASY TO READ)
      DO JJ=2,J2; J=JJ;     !(J:CELL, JJ:INTERIOR FACE, EASY TO READ)
        IF (IN(I,J,K)*IN(I,J+1,K) == 1) THEN
          TMP=IN(I,J-1,K)*IN(I,J,K)*IN(I,J+1,K)*IN(I,J+2,K)
          VFACE(II,JJ)=6.*(P(I,J,K)+P(I,J+1,K)) &
                         +(P(I,J,K)+P(I,J+1,K)-P(I,J-1,K)-P(I,J+2,K))*TMP
        ELSE
          VFACE(II,JJ)=12.*(IN(I,J,K)*P(I,J,K)+IN(I,J+1,K)*P(I,J+1,K))          
        END IF             
      END DO
      ! BOUNDARY FACES (JJ=1:SBC, JJ=J1:NBC) (NO PERIODIC B.C. ALONG LAT.-DIR.)
      VFACE(II,1)=12.*P(I,2,K); VFACE(II,J1)=12.*P(I,J1,K)
      DO J=2,J1; JJ=J;      !(J:CELL, JJ:FACE, EASY TO READ)
        PY(I,J)=O12*(VFACE(II,JJ)-VFACE(II,JJ-1))*IN(I,J,K)        
      END DO
    END DO
    ! B. VERTICAL FLUXES   (KK=1&K0:0, KK=2&K1:2ND-ORDER, KK=3,K2:4TH-ORDER)
    IF (K .EQ. 1) THEN      ! TOP FLUX OF TOP LAYER = 0

      DO J=2,J1; JJ=J-1;    !(J:CELL, JJ:SIZE OF UZ, ETC., EASY TO READ)
        DO I=2,I1; II=I-1;  !(I:CELL, II:SIZE OF UZ, ETC., EASY TO READ)
          UZ(II,JJ,LB)=0.
#ifdef FLAG_FS
          UZ(II,JJ,LB)=W(I,J,1)*U2(I,J,1)
#endif
          VZ(II,JJ,LB)=0.
#ifdef FLAG_FS
          VZ(II,JJ,LB)=W(I,J,1)*V2(I,J,1)
#endif
#ifdef FLAG_TS_T
          TZ(II,JJ,LB)=0.
#ifdef FLAG_FS
          TZ(II,JJ,LB)=W(I,J,1)*T2(I,J,1)
#endif
#endif
#ifdef FLAG_TS_S
          SZ(II,JJ,LB)=0.
#ifdef FLAG_FS
          SZ(II,JJ,LB)=W(I,J,1)*S2(I,J,1)
#endif
#endif
#ifdef FLAG_TRACER
          CZ(II,JJ,LB)=0. 
#ifdef FLAG_FS
          CZ(II,JJ,LB)=W(I,J,1)*C2(I,J,1)
#endif
#endif
        END DO
      END DO


    END IF
    IF (KK .LT. K0) THEN    ! BOTTOM FLUX (KK=2~K0) OF EACH LAYER (K=1~K1)
      DO J=2,J1             !(J:CELL, I:CELL) 
        DO I=2,I1           !(K:CURRENT LAYER, K+1:LOWER LAYER)
          SCR(I,J,1)=6.*(U2(I,J,K)+U2(I,J,K+1))
          SCR(I,J,2)=6.*(V2(I,J,K)+V2(I,J,K+1))
#ifdef FLAG_TS_T
          SCR(I,J,3)=6.*(T2(I,J,K)+T2(I,J,K+1))
#endif
#ifdef FLAG_TS_S
          SCR(I,J,4)=6.*(S2(I,J,K)+S2(I,J,K+1))
#endif
#ifdef FLAG_TRACER
          SCR(I,J,5)=6.*(C2(I,J,K)+C2(I,J,K+1))
#endif
        END DO
      END DO
      IF ((KK .GE. 3) .AND. (KK .LE. K2)) THEN !(KK=3,K2:4TH-ORDER)
        DO J=2,J1                              !(J:CELL, I:CELL)
          DO I=2,I1
            TMP=IN(I,J,K-1)*IN(I,J,K)*IN(I,J,K+1)*IN(I,J,K+2)
            !U2(I,J,K-1) IS UPDATED, THUS USE ULF(I,J,K-1)
            SCR(I,J,1)=SCR(I,J,1) &
                    +(-ULF(I,J,K-1)+U2(I,J,K)+U2(I,J,K+1)-U2(I,J,K+2))*TMP 
            SCR(I,J,2)=SCR(I,J,2) &
                    +(-VLF(I,J,K-1)+V2(I,J,K)+V2(I,J,K+1)-V2(I,J,K+2))*TMP
#ifdef FLAG_TS_T
             SCR(I,J,3)=SCR(I,J,3) &
                    +(-TLF(I,J,K-1)+T2(I,J,K)+T2(I,J,K+1)-T2(I,J,K+2))*TMP
#endif
#ifdef FLAG_TS_S
            SCR(I,J,4)=SCR(I,J,4) &
                    +(-SLF(I,J,K-1)+S2(I,J,K)+S2(I,J,K+1)-S2(I,J,K+2))*TMP
#endif
!            IF (FL_TRACER_ON ==1) &
#ifdef FLAG_TRACER
            SCR(I,J,5)=SCR(I,J,5) &
                    +(-CLF(I,J,K-1)+C2(I,J,K)+C2(I,J,K+1)-C2(I,J,K+2))*TMP
#endif
          END DO
        END DO
      END IF
      DO J=2,J1
        DO I=2,I1
          DO N=1,2
            SCR(I,J,N)=O12*SCR(I,J,N)
          END DO
#ifdef FLAG_TS_T
            SCR(I,J,3)=O12*SCR(I,J,3)
#endif
#ifdef FLAG_TS_S
            SCR(I,J,4)=O12*SCR(I,J,4)
#endif
#ifdef FLAG_TRACER
            SCR(I,J,5)=O12*SCR(I,J,5)
#endif
        END DO
      END DO
      DO J=2,J1; JJ=J-1;    !(J:CELL, JJ:SIZE OF UZ, I:CELL, II:SIZE OF UZ)
        DO I=2,I1; II=I-1;  !(KK:BOTTOM FACE, K+1:LOWER LAYER, K:CURRENT LAYER)
            UZ(II,JJ,LT)=W(I,J,KK)*SCR(I,J,1) &
                        -EV(II,JJ,K)*(U1(I,J,K+1)-U1(I,J,K))*IW(I,J,KK)
            VZ(II,JJ,LT)=W(I,J,KK)*SCR(I,J,2) &
                        -EV(II,JJ,K)*(V1(I,J,K+1)-V1(I,J,K))*IW(I,J,KK)
#ifdef FLAG_TS_T
	    TZ(II,JJ,LT)=W(I,J,KK)*SCR(I,J,3)&
                        -HV(II,JJ,K)*(T1(I,J,K+1)-T1(I,J,K))*IW(I,J,KK)
#endif
#ifdef FLAG_TS_S
            SZ(II,JJ,LT)=W(I,J,KK)*SCR(I,J,4)&
                  -ALPHA*HV(II,JJ,K)*(S1(I,J,K+1)-S1(I,J,K))*IW(I,J,KK)
#endif
#ifdef FLAG_TRACER
            CZ(II,JJ,LT)=W(I,J,KK)*SCR(I,J,5)&
                        -HV(II,JJ,K)*(C1(I,J,K+1)-C1(I,J,K))*IW(I,J,KK)
#endif
        END DO
      END DO            
    ELSE IF (KK .EQ. K0) THEN ! BOTTOM FLUX OF THE BOTTOM LAYER = 0 (DRAG LAW AT #9)
      DO J=2,J1; JJ=J-1;    !(J:CELL, JJ:SIZE OF UZ, ETC., EASY TO READ)
        DO I=2,I1; II=I-1;  !(I:CELL, II:SIZE OF UZ, ETC., EASY TO READ)
          UZ(II,JJ,LT)=0.
          VZ(II,JJ,LT)=0.
#ifdef FLAG_TS_T
          TZ(II,JJ,LT)=0.
#endif
#ifdef FLAG_TS_S
          SZ(II,JJ,LT)=0.
#endif
!          IF (FL_TRACER_ON ==1) CZ(II,JJ,LT)=0.
#ifdef FLAG_TRACER
          CZ(II,JJ,LT)=0.
#endif
        END DO
      END DO
    END IF
    ! C. LONGITUDINAL FLUXES
    DO J=2,J1; JJ=J-1;      !(J:CELL, JJ:SIZE OF UZ, ETC., EASY TO READ)
      DO I=1,I1             !(I:CELL, II:FACE, II=I HERE)
        SCR(I,J,1)=6.*(U2(I,J,K)+U2(I+1,J,K))
        SCR(I,J,2)=6.*(V2(I,J,K)+V2(I+1,J,K))
#ifdef FLAG_TS_T
        SCR(I,J,3)=6.*(T2(I,J,K)+T2(I+1,J,K))
#endif
#ifdef FLAG_TS_S
        SCR(I,J,4)=6.*(S2(I,J,K)+S2(I+1,J,K))
#endif
#ifdef FLAG_TRACER
        SCR(I,J,5)=6.*(C2(I,J,K)+C2(I+1,J,K))
#endif
      END DO
      DO I=2,I2             !(I:CELL, II:INTERIOR FACE, II=I HERE)
        TMP=IN(I-1,J,K)*IN(I,J,K)*IN(I+1,J,K)*IN(I+2,J,K)
        SCR(I,J,1)=SCR(I,J,1) &
                +(-U2(I-1,J,K)+U2(I,J,K)+U2(I+1,J,K)-U2(I+2,J,K))*TMP
        SCR(I,J,2)=SCR(I,J,2) &
                +(-V2(I-1,J,K)+V2(I,J,K)+V2(I+1,J,K)-V2(I+2,J,K))*TMP
#ifdef FLAG_TS_T
        SCR(I,J,3)=SCR(I,J,3) &
                +(-T2(I-1,J,K)+T2(I,J,K)+T2(I+1,J,K)-T2(I+2,J,K))*TMP
#endif
#ifdef FLAG_TS_S
        SCR(I,J,4)=SCR(I,J,4) &
                +(-S2(I-1,J,K)+S2(I,J,K)+S2(I+1,J,K)-S2(I+2,J,K))*TMP
#endif
#ifdef FLAG_TRACER
        SCR(I,J,5)=SCR(I,J,5) &
                +(-C2(I-1,J,K)+C2(I,J,K)+C2(I+1,J,K)-C2(I+2,J,K))*TMP
#endif
      END DO
      ! BOUNDARY FACES (II=I=1:WBC, II=I=I1:EBC) 
!!      IF (LOPENE == 1) THEN ! LOPENW = LOPENE = 1 ! PERIODIC
#ifdef FLAG_PERIODIC_WE
        ! WEST
        TMP=IN(I2,J,K)*IN(I1,J,K)*IN(2,J,K)*IN(3,J,K)
        SCR(1,J,1)=SCR(1,J,1)+TMP*(-U2(I2,J,K)+U2(I1,J,K)+U2(2,J,K)-U2(3,J,K))
        SCR(1,J,2)=SCR(1,J,2)+TMP*(-V2(I2,J,K)+V2(I1,J,K)+V2(2,J,K)-V2(3,J,K))
#ifdef FLAG_TS_T
        SCR(1,J,3)=SCR(1,J,3)+TMP*(-T2(I2,J,K)+T2(I1,J,K)+T2(2,J,K)-T2(3,J,K))
#endif
#ifdef FLAG_TS_S
        SCR(1,J,4)=SCR(1,J,4)+TMP*(-S2(I2,J,K)+S2(I1,J,K)+S2(2,J,K)-S2(3,J,K))
#endif
#ifdef FLAG_TRACER
        SCR(1,J,5)=SCR(1,J,5)+TMP*(-C2(I2,J,K)+C2(I1,J,K)+C2(2,J,K)-C2(3,J,K))
#endif
        ! EAST
        SCR(I1,J,1)=SCR(1,J,1)         
        SCR(I1,J,2)=SCR(1,J,2)         
#ifdef FLAG_TS_T
        SCR(I1,J,3)=SCR(1,J,3)
#endif
#ifdef FLAG_TS_S
        SCR(I1,J,4)=SCR(1,J,4)
#endif
!        IF (FL_TRACER_ON ==1) SCR(I1,J,5)=SCR(1,J,5)
#ifdef FLAG_TRACER
        SCR(I1,J,5)=SCR(1,J,5)
#endif

#endif
!!      END IF
      DO I=1,I1
        DO N=1,2
          SCR(I,J,N)=O12*SCR(I,J,N)
        END DO
#ifdef FLAG_TS_T
          SCR(I,J,3)=O12*SCR(I,J,3)
#endif
#ifdef FLAG_TS_S
            SCR(I,J,4)=O12*SCR(I,J,4)
#endif
#ifdef FLAG_TRACER
          SCR(I,J,5)=O12*SCR(I,J,5)
#endif
      END DO
      DO I=1,I1             !(I:CELL, II:FACE, II=I HERE)
        AM=DMX(I,J,K)*ODX(J)
        !AH=DMX(I,J,K)*ODX(J)*(1./PRN)
        ah=dmx(i,j,k)*odx(j)*dm0(i,j,k)/de0(i,j,k)
        UX(I,JJ)=(U(I,J,K)*SCR(I,J,1)       &
                -AM*(U1(I+1,J,K)-U1(I,J,K)))*IU(I,J,K)
        VX(I,JJ)=(U(I,J,K)*SCR(I,J,2)       &
                -AM*(V1(I+1,J,K)-V1(I,J,K)))*IU(I,J,K)
#ifdef FLAG_TS_T
        TX(I,JJ)=(U(I,J,K)*SCR(I,J,3)       &
                -AH*(T1(I+1,J,K)-T1(I,J,K)))*IU(I,J,K)
#endif
#ifdef FLAG_TS_S
        SX(I,JJ)=(U(I,J,K)*SCR(I,J,4)       &
                -AH*(S1(I+1,J,K)-S1(I,J,K)))*IU(I,J,K)
#endif
!        IF (FL_TRACER_ON ==1)               &       
#ifdef FLAG_TRACER
        CX(I,JJ)=(U(I,J,K)*SCR(I,J,5)       &
                -AH*(C1(I+1,J,K)-C1(I,J,K)))*IU(I,J,K)
#endif
      END DO
    END DO
    ! D. LATITUDINAL FLUXES
    DO I=2,I1; II=I-1;      !(I:CELL, II:SIZE OF UY, ETC., EASY TO READ)
      DO J=1,J1             !(J:CELL, JJ:FACE, JJ=J HERE)
        SCR(I,J,1)=6.*(U2(I,J,K)+U2(I,J+1,K))
        SCR(I,J,2)=6.*(V2(I,J,K)+V2(I,J+1,K))
#ifdef FLAG_TS_T
        SCR(I,J,3)=6.*(T2(I,J,K)+T2(I,J+1,K))
#endif
#ifdef FLAG_TS_S
        SCR(I,J,4)=6.*(S2(I,J,K)+S2(I,J+1,K)) 
#endif
#ifdef FLAG_TRACER
        SCR(I,J,5)=6.*(C2(I,J,K)+C2(I,J+1,K))
#endif
      END DO
      DO J=2,J2             !(J:CELL, JJ:INTERIOR FACE, JJ=J HERE)
        TMP=IN(I,J-1,K)*IN(I,J,K)*IN(I,J+1,K)*IN(I,J+2,K)        
        SCR(I,J,1)=SCR(I,J,1) &
                +(-U2(I,J-1,K)+U2(I,J,K)+U2(I,J+1,K)-U2(I,J+2,K))*TMP
        SCR(I,J,2)=SCR(I,J,2)+TMP*&
                +(-V2(I,J-1,K)+V2(I,J,K)+V2(I,J+1,K)-V2(I,J+2,K))*TMP
#ifdef FLAG_TS_T
        SCR(I,J,3)=SCR(I,J,3) &
                +(-T2(I,J-1,K)+T2(I,J,K)+T2(I,J+1,K)-T2(I,J+2,K))*TMP
#endif
#ifdef FLAG_TS_S
        SCR(I,J,4)=SCR(I,J,4) &
                +(-S2(I,J-1,K)+S2(I,J,K)+S2(I,J+1,K)-S2(I,J+2,K))*TMP
#endif
#ifdef FLAG_TRACER
        SCR(I,J,5)=SCR(I,J,5) &
                +(-C2(I,J-1,K)+C2(I,J,K)+C2(I,J+1,K)-C2(I,J+2,K))*TMP
#endif
      END DO
      ! NO PERIODIC B.C.
      DO J=1,J1
        DO N=1,2
          SCR(I,J,N)=O12*SCR(I,J,N)
        END DO
#ifdef FLAG_TS_T
          SCR(I,J,3)=O12*SCR(I,J,3)
#endif
#ifdef FLAG_TS_S
          SCR(I,J,4)=O12*SCR(I,J,4)
#endif
#ifdef FLAG_TRACER
          SCR(I,J,5)=O12*SCR(I,J,5)
#endif
      END DO
      DO J=1,J1             !(J:CELL, JJ:FACE, JJ=J HERE)
        AM=DMY(I,J,K)*ODYV(J)
        AH=DMY(I,J,K)*ODYV(J)*dm0(i,j,k)/de0(i,j,k)
        UY(II,J)=CSV(J)*(V(I,J,K)*SCR(I,J,1)&
                -AM*(U1(I,J+1,K)-U1(I,J,K)))*IV(I,J,K)
        VY(II,J)=CSV(J)*(V(I,J,K)*SCR(I,J,2)&
                -AM*(V1(I,J+1,K)-V1(I,J,K)))*IV(I,J,K)
#ifdef FLAG_TS_T
        TY(II,J)=CSV(J)*(V(I,J,K)*SCR(I,J,3)&
                -AH*(T1(I,J+1,K)-T1(I,J,K)))*IV(I,J,K)
#endif
#ifdef FLAG_TS_S
        SY(II,J)=CSV(J)*(V(I,J,K)*SCR(I,J,4)&
                -AH*(S1(I,J+1,K)-S1(I,J,K)))*IV(I,J,K)
#endif
#ifdef FLAG_TRACER
        CY(II,J)=CSV(J)*(V(I,J,K)*SCR(I,J,5)&
                -AH*(C1(I,J+1,K)-C1(I,J,K)))*IV(I,J,K)
#endif
      END DO
    END DO
    ! ================================================================== !
    ! 5. MOMENTUM EQUATIONS & TRANSPORT EQUATIONS (TMP U2,V2,T2,S2,C2)   !
    ! ================================================================== !
    ! PARTIALLY UPDATED, CORIOLIS & BOUNDARY FORCINGS ADDED LATER
    DO J=2,J1; JJ=J-1;      !(J:CELL & NOR. FACE, JJ:SIZE OF UX, ETC., EASY TO READ)
       ODYJ=OCS(J)*ODY(J)
       DO I=2,I1; II=I-1;   !(I:CELL & EST. FACE, II:SIZE OF UY, ETC., EASY TO READ)
         DTIN=DT*IN(I,J,K)  !(K:LAYER, LT:BOTTOM FACE, LB:TOP FACE)
         ! LONGITUDINAL MOMENTUM
         U2(I,J,K)=U1(I,J,K)-DTIN*                          &
                 ((UX(I,JJ)-UX(I-1,JJ)+PX(I,J))*ODX(J)      &
                 +(UY(II,J)-UY(II,J-1))*ODYJ                &
                 +(UZ(II,JJ,LT)-UZ(II,JJ,LB))*ODZ(K))
         ! LATITUDINAL MOMENTUM
         V2(I,J,K)=V1(I,J,K)-DTIN*                          &
                 ((VX(I,JJ)-VX(I-1,JJ))*ODX(J)              &
                 +(VY(II,J)-VY(II,J-1))*ODYJ+PY(I,J)*ODY(J) &
                 +(VZ(II,JJ,LT)-VZ(II,JJ,LB))*ODZ(K))
         ! TEMPERATURE
#ifdef FLAG_TS_T
         T2(I,J,K)=T1(I,J,K)-DTIN*                          &
                 ((TX(I,JJ)-TX(I-1,JJ))*ODX(J)              &
                 +(TY(II,J)-TY(II,J-1))*ODYJ                &
                 +(TZ(II,JJ,LT)-TZ(II,JJ,LB))*ODZ(K))
#endif
#ifdef FLAG_TS_S
         S2(I,J,K)=S1(I,J,K)-DTIN*                          &
                 ((SX(I,JJ)-SX(I-1,JJ))*ODX(J)              &
                 +(SY(II,J)-SY(II,J-1))*ODYJ                &
                 +(SZ(II,JJ,LT)-SZ(II,JJ,LB))*ODZ(K))
#endif
#ifdef FLAG_TRACER
         ! TRACER
         C2(I,J,K)=C1(I,J,K)-DTIN*                          &
                 ((CX(I,JJ)-CX(I-1,JJ))*ODX(J)              &
                 +(CY(II,J)-CY(II,J-1))*ODYJ                &
                 +(CZ(II,JJ,LT)-CZ(II,JJ,LB))*ODZ(K))
#endif
      END DO
    END DO 
    I=LB; LB=LT; LT=I
  END DO
#ifdef _DEBUG_
  PRINT*,'Pass Step 4,5'
#endif
  ! ============================= !
  ! 6. OPEN BOUNDARY CONDITIONS   !
  ! ============================= !
  TMP=.5*DT
  IF (LOPENW .GE. 2) THEN   ! WEST
    DO K=1,K1
      DO J=2,J1
        TMPIN=IN(2,J,K)*TMP*ODX(J)
        TEMP=ABS(U(1,J,K))
        TEMP1=TMPIN*(TEMP+U(1,J,K))
        TEMP2=TMPIN*(TEMP-U(1,J,K))
        U2(2,J,K)=U2(2,J,K)+TEMP1*U2(1,J,K)-TEMP2*ULF(2,J,K) !U2, V2 HAS BEEN UPDATED
        V2(2,J,K)=V2(2,J,K)+TEMP1*V2(1,J,K)-TEMP2*VLF(2,J,K) !THUS USE ULF, VLF
#ifdef FLAG_TS_T
        T2(2,J,K)=T2(2,J,K)+TEMP1*T2(1,J,K)-TEMP2*TLF(2,J,K)
#endif
#ifdef FLAG_TS_S
        S2(2,J,K)=S2(2,J,K)+TEMP1*S2(1,J,K)-TEMP2*SLF(2,J,K)
#endif
#ifdef FLAG_TRACER
        C2(2,J,K)=C2(2,J,K)+TEMP1*C2(1,J,K)-TEMP2*CLF(2,J,K)
#endif
      END DO
    END DO
  END IF
  IF (LOPENE .GE. 2) THEN   ! EAST
    DO K=1,K1
      DO J=2,J1
        TMPIN=IN(I1,J,K)*TMP*ODX(J)
        TEMP=ABS(U(I1,J,K))
        TEMP1=TMPIN*(TEMP+U(I1,J,K))
        TEMP2=TMPIN*(TEMP-U(I1,J,K))
        U2(I1,J,K)=U2(I1,J,K)-TEMP1*ULF(I1,J,K)+TEMP2*U2(I0,J,K)
        V2(I1,J,K)=V2(I1,J,K)-TEMP1*VLF(I1,J,K)+TEMP2*V2(I0,J,K)
#ifdef FLAG_TS_T
        T2(I1,J,K)=T2(I1,J,K)-TEMP1*TLF(I1,J,K)+TEMP2*T2(I0,J,K)
#endif
#ifdef FLAG_TS_S
        S2(I1,J,K)=S2(I1,J,K)-TEMP1*SLF(I1,J,K)+TEMP2*S2(I0,J,K)
#endif
#ifdef FLAG_TRACER
        C2(I1,J,K)=C2(I1,J,K)-TEMP1*CLF(I1,J,K)+TEMP2*C2(I0,J,K)
#endif
      END DO
    END DO    
  END IF
  IF (LOPENS .GE. 2) THEN   ! SOUTH
    DO K=1,K1
      DO I=2,I1
        TMPIN=IN(I,2,K)*TMP*ODY(2)*CSV(1)*OCS(2)
        TEMP=ABS(V(I,1,K))
        TEMP1=TMPIN*(TEMP+V(I,1,K))
        TEMP2=TMPIN*(TEMP-V(I,1,K))
        U2(I,2,K)=U2(I,2,K)+TEMP1*U2(I,1,K)-TEMP2*ULF(I,2,K)
        V2(I,2,K)=V2(I,2,K)+TEMP1*V2(I,1,K)-TEMP2*VLF(I,2,K)
#ifdef FLAG_TS_T
        T2(I,2,K)=T2(I,2,K)+TEMP1*T2(I,1,K)-TEMP2*TLF(I,2,K)
#endif
#ifdef FLAG_TS_S
        S2(I,2,K)=S2(I,2,K)+TEMP1*S2(I,1,K)-TEMP2*SLF(I,2,K)
#endif
#ifdef FLAG_TRACER
        C2(I,2,K)=C2(I,2,K)+TEMP1*C2(I,1,K)-TEMP2*CLF(I,2,K)
#endif
      END DO
    END DO
  END IF
  IF (LOPENN .GE. 2) THEN   ! NORTH
    DO K=1,K1
      DO I=2,I1
        TMPIN=IN(I,J1,K)*TMP*ODY(J1)*CSV(J1)*OCS(J1)
        TEMP=ABS(V(I,J1,K))
        TEMP1=TMPIN*(TEMP+V(I,J1,K))
        TEMP2=TMPIN*(TEMP-V(I,J1,K))
        U2(I,J1,K)=U2(I,J1,K)-TEMP1*ULF(I,J1,K)+TEMP2*U2(I,J0,K)
        V2(I,J1,K)=V2(I,J1,K)-TEMP1*VLF(I,J1,K)+TEMP2*V2(I,J0,K)
#ifdef FLAG_TS_T
        T2(I,J1,K)=T2(I,J1,K)-TEMP1*TLF(I,J1,K)+TEMP2*T2(I,J0,K)
#endif
#ifdef FLAG_TS_S
        S2(I,J1,K)=S2(I,J1,K)-TEMP1*SLF(I,J1,K)+TEMP2*S2(I,J0,K)
#endif
#ifdef FLAG_TRACER
        C2(I,J1,K)=C2(I,J1,K)-TEMP1*CLF(I,J1,K)+TEMP2*C2(I,J0,K)
#endif
      END DO
    END DO
  END IF
#ifdef _DEBUG_
  PRINT*,'Pass Step 6'
#endif

  ! ==================== ! 
  ! 7. SURFACE EFFECTS   !
  ! ==================== !
  ! preliminary: compute bluk formulae if necessary
  if (bulk_formulae) then
     space=.false.
     call interfgetfield(localWindU10m,modelmjd,tempU10M,space)
     call interfgetfield(localWindV10m,modelmjd,tempV10M,space)
     call interfgetfield(localAirTemperature2m,modelmjd,tempT2M,space)
     call interfgetfield(localDewTemperature2m,modelmjd,tempDT2M,space)
     call interfgetfield(localAtmPressure,modelmjd,tempPMSL,space)
     if (cloudtype.eq.1) then
       call interfgetfield(localCloudCoverage,modelmjd,tempCC,space)
     else
       tempCC=0.0
     end if
     call compute_fluxes(modelmjd,xdeg(2:I1),ydeg(2:J1),tempCC,tempU10M,tempV10M,tempPMSL,T2(:,:,1),tempT2M,tempDT2M,qdot,qdot2,taux1,tauy1,evapo)
     !convert taux1 to the right unit and mask
     DO J=1,J2
       DO I=1,I2
         taux1(I,J)=10.*taux1(I,J)*IN(I+1,J+1,1)
	 tauy1(I,J)=10.*tauy1(I,J)*IN(I+1,J+1,1)
       END DO
     END DO
  end if

  ! A. HEAT FLUXES
  if (heattype.eq.2) then
      !overwrite qdot and qdot2 (obtained from bulk formulae) with values read from disk
      !interpolate Qup, Qsolar in time --> QDOT,QDOT2
      space=.false.
      call interfgetfield(localQup,modelmjd,QDOT,space)
      call interfgetfield(localQsolar,modelmjd,QDOT2,space)
  end if
  if (heattype.eq.3) then
      space=.false.
      !shortwave
      call interfgetfield(localQsolar,modelmjd,QDOT2,space)
      !longwave: emission, sensible and latent from qdot (from compute_fluxes), back-radiation from localQup (from disk)
      call interfgetfield(localQup,modelmjd,tempCC,space) !tempCC now temporarily holds longwave back-radiation
      qdot=qdot+tempCC ! qdot contains 3 contributions from bulk formulae, all positive upwards
		       ! tempCC is the back radiation read from disk, follows the same convention
		       !   but is a negative number as it actually goes down
  end if
  if (heattype.ge.0) then
     !actually add the heat, whether it comes from bulk formulae or fluxes read from disk
     !call QSURF(T1(:,:,1),modelmjd,dt,odz)
     call QSURF(T2(:,:,1),modelmjd,dt,odz)
     !call QSURF(TLF(:,:,1),modelmjd,dt,odz)
     !call QVOLUME(T1,modelmjd,dt,odz)
     call QVOLUME(T2,modelmjd,dt,odz)
     !call QVOLUME(TLF,modelmjd,dt,odz)
  end if

  ! B. WIND STRESS
  if (windtype.eq.-1) then
     IF (FL_WD_ON == 1) CALL WIND(U2(:,:,1),V2(:,:,1),DAYS,DT,ODZ)
  elseif (windtype.eq.1) then
      !overwrite bulk formulae tauX and tauY with values read from disk
      space=.false.
      call interfgetfield(localMomentfluxX,modelmjd,taux1,space)
      call interfgetfield(localMomentfluxY,modelmjd,tauy1,space)
      !convert taux1 to the right unit and mask
      DO J=1,J2
        DO I=1,I2
          taux1(I,J)=10.*taux1(I,J)*IN(I+1,J+1,1)
          tauy1(I,J)=10.*tauy1(I,J)*IN(I+1,J+1,1)
        END DO
      END DO
  end if
  if (windtype.ge.0) then
      !actually add the surface stress to ocean surface current
      !CALL WIND(U1(:,:,1),V1(:,:,1),modelmjd,dt,odz)
      CALL WIND(U2(:,:,1),V2(:,:,1),modelmjd,dt,odz)
      !CALL WIND(ULF(:,:,1),VLF(:,:,1),modelmjd,dt,odz)
  end if

  ! C. EVAPORATION & PRECIPITATION (NO SALT SOURCES)
  if (salttype.eq.0) then
      space=.false.
      !overwrite the evaporation with values read from disk, need to be kg/m2/s
      call interfgetfield(localEvap,modelmjd,evapo,space)
  ! elseif (salttype.eq.-1)  evapo=0  ! we will initialize evapo later, if salttype.eq.-1 AND raintype.eq.1, because in that case, rain must be substracted from it  ! 
  ! elseif (salttype.eq.1) then it is already calculated with the bulk formulae and stored in evapo
  end if
  if (raintype.eq.1) then
      space=.false.
      !precipitation can only be read from disk, and should be in kg/m2/s
      call interfgetfield(localPrecipitation,modelmjd,rain,space)
      !remove precipitation from evaporation
      !write(*,*) "sizes=",size(evapo,1),size(evapo,2),size(rain,1),size(rain,2)
      !write(*,*) "mean evapo,rain=",sum(evapo)/size(evapo,1)/size(evapo,2),sum(rain)/size(rain,1)/size(rain,2)
      !write(*,*) "max  evapo,rain=",maxval(evapo),maxval(rain)
      if (salttype.eq.-1) evapo=0.0   ! initialize evapo
      evapo=evapo-rain
  end if
  if (salttype.ge.0.or.raintype.eq.1) then
      !CALL SalSURF(S1(:,:,1),modelmjd,dt,odz)
      CALL SalSURF(S2(:,:,1),modelmjd,dt,odz)
      !CALL SalSURF(SLF(:,:,1),modelmjd,dt,odz)
  end if

  ! D. CLIMATOLOGICAL SURFACE RESTORING
  IF (FL_SURF_EP == 1) THEN
    ! TAKE OUT MOMENTUM WITH POSITIVE E-P; ADD NONE FOR NEGATIVE E-P
    ! CLIMATOLOGICAL SURFACE RESTORING
    DO J=2,J1; JJ=J-1;
      DO I=2,I1; II=I-1;
        ! QAVG=HEAT/TIME IN DEG. C PER TIME STEP
        TEMP=IN(I,J,1)*NUDGE(I,J,1)*QAVG(II,JJ,NLD)
        T2(I,J,1) =T2(I,J,1) +TEMP
        T1(I,J,1) =T1(I,J,1) +TEMP
        TLF(I,J,1)=TLF(I,J,1)+TEMP
        ! NONLINEARLY NUDGE TOWARD SURFACE T CLIMATOLOGY
        TKLI=FOLD*TSURF(II,JJ,NLD)+FNEW*TSURF(II,JJ,NEW)
        TEMP=T2(I,J,1)
        TMP1=TAU2D(i,j)*(1.+(TEMP-TKLI)**2)
        TMP=1.-TMP1
        !T2(I,J,1)=IN(I,J,1)*NUDGE(I,J,1)*(TMP*TEMP+TMP1*TKLI)+(1-IN(I,J,1)*NUDGE(I,J,1))*TKLI
        T2(I,J,1)=IN(I,J,1)*NUDGE(I,J,1)*(TMP*TEMP+TMP1*TKLI)+IN(I,J,1)*(1-NUDGE(I,J,1))*T2(I,J,1)+(1-IN(I,J,1))*TKLI
        T2(I,J,1)=MAX(-1.8,T2(I,J,1))
        ! ACCUMULATE NUDGINGS FOR PRESENT MONTH
        TNUDGE(II,JJ)=TNUDGE(II,JJ)+T2(I,J,1)-TEMP
      END DO
    END DO
    ! SURFACE FLUXES
    I=30*ITFDAY
    IF (MOD(ITF,I) .EQ. I-1) THEN
      NSOMBO(NLD)=NSOMBO(NLD)+1 ! NSOMBO=1 WHEN QAVG & WAVG ARE INITIALIZED
      ENS1=1./NSOMBO(NLD); ENS2=1.-ENS1
      TEMP=1./(30*ITFDAY) ! TEMP=1/YEARS/(30*ITFDAY), FOR LONG-TERM AVERAGE
      TMP1=TEMP*ENS1; TMP2=Z(3)*ODT*ENS1
      DO J=2,J1; JJ=J-1;
        TMP=DX(J)*DY(J)
        DO I=2,I1; II=I-1;
          ! EVAP IS THE NET SURFACE OUTFLOW (POSITIVE FOR NET POSITIVE E-P)
          ! WAVG=-(DZ*DS/S)/DT, SELINC: CHANGES OF SALINITY PER TIME STEP
          SELINC=-TEMP*(S2(I,J,1)-SSURF(II,JJ,NEW))
          WAVG(II,JJ,NLD)=(WAVG(II,JJ,NLD)&
                         -TMP2*SELINC/S2(I,J,1))*IN(I,J,1)
          QAVG(II,JJ,NLD)=QAVG(II,JJ,NLD)&
                         +TMP1*(TNUDGE(II,JJ)+TSURF(II,JJ,NEW)-T2(I,J,1))
          ! SET MONTHLY SUM OF NUDGINGS = 0 TO INITIALIZE NEXT MONTH
          TNUDGE(II,JJ)=0.
        END DO
      END DO
      ! GIVE REASONABLE ADVECTION ESTIMATE FOR NEXT MONTH, FIRST YEAR ONLY
      IF (DAYS.LE.361.) THEN
        DO JJ=1,J2
          DO II=1,I2
            WAVG(II,JJ,NEW)=WAVG(II,JJ,NLD)
            QAVG(II,JJ,NEW)=QAVG(II,JJ,NLD)
          END DO
        END DO
      END IF
    END IF
    ! NUDGE W TOWARD MONTHLY MEAN
    TMP1=1./(NUDGE_W*DAODT)
    DO J=2,J1; JJ=J-1;
      DO I=2,I1; II=I-1;
        W(I,J,1)=IN(I,J,1)*(W(I,J,1)+TMP1*nudge(i,j,1)*(WAVG(II,JJ,NLD)-W(I,J,1)))
      END DO
    END DO
  ENDIF
    
  !CALCULATE NET SURFACE FLUX
  EVAP=0.
  DO J=2,J1
    TEMP=DX(J)*DY(J)
    DO I=2,I1
      EVAP=EVAP-W(I,J,1)*TEMP
    END DO
  END DO 
  TMPBC=EVAP/A(1)
  IF(FL_INFL_W+FL_INFL_E+FL_INFL_S+FL_INFL_N.EQ.0)THEN
   DO J=2,J1
    DO I=2,I1
      W(I,J,1)=IN(I,J,1)*(W(I,J,1)+TMPBC)
    END DO
   END DO
  ELSE
    TMPBC=EVAP
  ENDIF
  ! ========================================= !
  ! 8. CLIMATOLOGCAL SPONGE LAYER RESTORING   !
  ! ========================================= !
  IF (FL_SPL_ON ==1) THEN
    N=10 ! 10 BOUNDARY LATITUDES ONLY
    DO K=2,K1; KK=K-1;      ! SURFACE SOURCES TAKE CARE OF K=1
      !TEMP=1 ! used to be TAUDTN, which is now tau3D(i,j,k)
      DO J=1,N
        !TEMP=.9*TEMP
        DO I=2,I1; II=I-1;
          tmp=1.-tau3d(i,j,k)
          ! NORTH
          TKLI=FOLD*TNSP(II,J,K,NLD)+FNEW*TNSP(II,J,K,NEW)
          TMP2=T2(I,J0-J,K)
          T2(I,J0-J,K)=TMP*TMP2+TAU3D(i,j,k)*TKLI
          SKLI=FOLD*SNSP(II,J,K,NLD)+FNEW*SNSP(II,J,K,NEW)
          TMP1=S2(I,J0-J,K)
          S2(I,J0-J,K)=TMP*TMP1+TAU3D(i,j,k)*SKLI
          !SOUTH    
          TKLI=FOLD*TSSP(II,J,K,NLD)+FNEW*TSSP(II,J,K,NEW)
          TMP2=T2(I,J+1,K)
          T2(I,J+1,K)=TMP*TMP2+TAU3D(i,j,k)*TKLI
          SKLI=FOLD*SSSP(II,J,K,NLD)+FNEW*SSSP(II,J,K,NEW)
          TMP1=S2(I,J+1,K)
          S2(I,J+1,K)=TMP*TMP1+TAU3D(i,j,k)*SKLI
        END DO
      END DO
    END DO
  END IF
  ! ================== !
  ! 9. BOTTOM STRESS   !
  ! ================== !
  TMP=DT*DRAG
  DO J=2,J1
    DO I=2,I1
      K=KB(I,J)
      ! NOT EXCEED EXPLICIT LIMIT
      DRG=MIN(.2,IN(I,J,K)*TMP*ODZ(K)*SQRT(U1(I,J,K)**2+V1(I,J,K)**2))
      U2(I,J,K)=U2(I,J,K)-DRG*U1(I,J,K)
      V2(I,J,K)=V2(I,J,K)-DRG*V1(I,J,K)
    END DO
  END DO
  ! ========================== !
  ! 10. TRAPEZOIDAL CORIOLIS   !
  ! ========================== !
  DO K=1,K1
    DO J=2,J1
      DO I=2,I1
        TMP=.5*(F(J-1)+ULF(I,J,K)*TANPHI(J-1))*DT
        TEMP=1./(1.+TMP**2)
        QU=U2(I,J,K)+TMP*V1(I,J,K)
        QV=V2(I,J,K)-TMP*U1(I,J,K)
        U2(I,J,K)=TEMP*(QU+TMP*QV)
        V2(I,J,K)=TEMP*(QV-TMP*QU)
      END DO
    END DO
  END DO
  ! ================================================ !
  ! 11. SET PERIODIC BOUNDARY CONDITIONS IF NEEDED   !
  ! ================================================ !
  IF (LOPENE == 1) THEN     !LOPENW = LOPENE = 1
    DO K=1,K1
      DO J=2,J1
       U2(1,J,K) =U2(I1,J,K)
       V2(1,J,K) =V2(I1,J,K)
       U2(I0,J,K)=U2(2,J,K)
       V2(I0,J,K)=V2(2,J,K)
      END DO
    END DO
  END IF  
#ifdef _DEBUG_
  PRINT*,'Pass Step 7-11'
#endif
  ! =============================================================================== !
  ! 12. INTERPOLATE "A" GRID QUANTITIES TO "C" GRID LOCATIONS, 4TH-ORDER-ACCURATE   !
  ! =============================================================================== !
  DO K=1,K1
    ! LONGITUDINAL DIRECTION 
    DO J=2,J1               !(J:CELL, K:LAYER) 
      DO I=2,I2             !(I:CELL, II:INTERIOR FACE, II=I HERE) 
        SCR(I,J,1)=6.*(U2(I,J,K)+U2(I+1,J,K))
        TMP=IN(I-1,J,K)*IN(I,J,K)*IN(I+1,J,K)*IN(I+2,J,K)
        SCR(I,J,1)=SCR(I,J,1)+TMP*&
                 (-U2(I-1,J,K)+U2(I,J,K)+U2(I+1,J,K)-U2(I+2,J,K))
        U(I,J,K)=O12*SCR(I,J,1)*IU(I,J,K)
      END DO
      ! BOUNDARY FACES (II=1:WBC, II=I1:EBC) (INITIAL DATA OR PERIODIC B.C.)
      IF (LOPENW == 1) THEN !LOPENW = LOPENE = 1 ! PERIODIC B.C.
        SCR(1,J,1)=6.*(U2(1,J,K)+U2(2,J,K))&
           -U2(I2,J,K)+U2(1,J,K)+U2(2,J,K)-U2(3,J,K)
        U(1,J,K) =O12*SCR(1,J,1)*IU(1,J,K)
        U(I1,J,K)=U(1,J,K)
      END IF
    END DO
    ! LATITUDINAL DIRECTION
    DO I=2,I1               !(I:CELL, K:LAYER)  
      DO J=2,J2             !(J:CELL, JJ:INTERIOR FACE, JJ=J HERE) 
        SCR(I,J,1)=6.*(V2(I,J,K)+V2(I,J+1,K))
        TMP=IN(I,J-1,K)*IN(I,J,K)*IN(I,J+1,K)*IN(I,J+2,K)
        SCR(I,J,1)=SCR(I,J,1)+TMP*&
                 (-V2(I,J-1,K)+V2(I,J,K)+V2(I,J+1,K)-V2(I,J+2,K))
        V(I,J,K)=O12*SCR(I,J,1)*IV(I,J,K)
      END DO
      ! BOUNDARY FACES (JJ=1:SBC, JJ=J1:NBC) (INITIAL DATA)
    END DO
  END DO
  ! ================================================= !
  ! 13. INFLOW & OUTFLOW CONDITIONS, CASE DEPENDENT   !
  ! ================================================= !

  !FL_INFL_W.GE.1  INCLUDE VELOCITY CORRECTION 
  !FL_INFL_W.EQ.2  INCLUDE VELOCITY CORRECTION BUT NOT ADD IT BACK
  !FL_INFL_W.EQ.3  INCLUDE VELOCITY CORRECTION AND ADD IT BACK
  !FL_INFL_W.EQ.4  ENFORCE OUTFLOW (NO INFLOW ALLOWED)
  !FL_INFL_W.EQ.5  ENFORCE INFLOW (NO OUTFLOW ALLOWED)

  IF(FL_INFL_W+FL_INFL_E+FL_INFL_S+FL_INFL_N.NE.0)THEN
  !IMPOSE LATER FOR E,S,N
  IF(FL_INFL_W.EQ.4)THEN ! NO INFLOW ALLOWED
    DO K=1,K1
      DO J=2,J1
        U(1,J,K)=MIN(0.,U(2,J,K));U1(1,J,K)=U(1,J,K);U2(1,J,K)=U(1,J,K)
      END DO
    END DO
  ELSEIF(FL_INFL_W.EQ.5)THEN 
    DO K=1,K1
      DO J=2,J1
        U(1,J,K)=MAX(0.,U(2,J,K));U1(1,J,K)=U(1,J,K);U2(1,J,K)=U(1,J,K)
      END DO
    END DO
  ENDIF

    ! SUBTRACT MEAN BOUNDARY OUTFLOW AS REQUIRED BY INCOMPRESSIBILITY
  IF(FL_INFL_W.GE.1)THEN 
    DO K=1,K1
      DO J=2,J1
        TMPBC=TMPBC-U(1,J,K)*IN(2,J,K)/(ODY(J)*ODZ(K)) 
      END DO
    END DO
  ENDIF
  IF(FL_INFL_E.GE.1)THEN 
    DO K=1,K1
      DO J=2,J1
        TMPBC=TMPBC+U(I1,J,K)*IN(I1,J,K)/(ODY(J)*ODZ(K)) 
      END DO
    END DO
  ENDIF
  IF(FL_INFL_S.GE.1)THEN 
    DO K=1,K1
      DO I=2,I1
        TMPBC=TMPBC-V(I,1,K)*IN(I,2,K)/(ODXV(1)*ODZ(K))                    
      END DO
    END DO
  ENDIF
  IF(FL_INFL_N.GE.1)THEN 
    DO K=1,K1
      DO I=2,I1
        TMPBC=TMPBC+V(I,J1,K)*IN(I,J1,K)/(ODXV(J1)*ODZ(K))                    
      END DO
    END DO
  ENDIF
  ! ADJUST OUTFLOWS
  TMPBC=TMPBC/AROUT
  IF ( (ITF_GENERAL .LT. TOPTS) .OR. (MOD(ITF_GENERAL,TOPTS) .EQ. 0) .OR. (ITFTO .LE. NITFTO) ) THEN
  !IF (MOD(ITF,ITFDAY).LT.3) WRITE(*,648) TMPBC
      WRITE(*,648) TMPBC
  ENDIF
648 FORMAT('outflow vel correction= ',1PE9.2,' cm/sec')
  IF(FL_INFL_W.GE.3)THEN 
    DO K=1,K1
      DO J=2,J1
        U(1,J,K)=(U(1,J,K)+TMPBC)*IN(2,J,K)
      END DO
    END DO
  ENDIF
  IF(FL_INFL_E.GE.3)THEN 
    DO K=1,K1
      DO J=2,J1
        U(I1,J,K)=(U(I1,J,K)-TMPBC)*IN(I1,J,K)
      END DO
    END DO
  ENDIF
  IF(FL_INFL_S.GE.3)THEN 
    DO K=1,K1
      DO I=2,I1
        V(I,1,K)=(V(I,1,K)+TMPBC)*IN(I,2,K)
      END DO
    END DO
  ENDIF
  IF(FL_INFL_N.GE.3)THEN 
    DO K=1,K1
      DO I=2,I1
        V(I,J1,K)=(V(I,J1,K)-TMPBC)*IN(I,J1,K)
      END DO
    END DO
  ENDIF
  ENDIF
#ifdef _DEBUG_
  PRINT*,'Pass Step 12-13'
#endif

  ! ============================================================================== !
  ! 14. ITERATE EVP SOLVER TO GET EXACTLY NON-DIVERGENT BAROTROPIC MODE (P0,U,V)   !
  ! ============================================================================== !
  ! USE SCR TO STORE CHANGES FOR INCOMPRESSIBILTY
  DO K=1,4; DO J=1,J0; DO I=1,I0; SCR(I,J,K)=0.; END DO; END DO; END DO
  !MXMASK=18 ! MXMASK: MAXIMUM ITERATIONS ALLOWED (SPECIFIED IN NML)
  ITMASK=0
  ITMASKLOOP: DO WHILE(.true.)
    ITMASK=ITMASK+1
    ! ZERO OUT ADVECTION VELOCITY OVER LAND
    TMP=0.
    DO J=2,J1               !(I&J:CELL, II;INTERIOR FACE, II=I HERE)
      DO I=2,I2 
        U(I,J,1)=IU0(I,J)*U(I,J,1)
      END DO
      IF (LOPENW == 1) THEN ! LOPENW = LOPENE = 1 ! PERIODIC B.C.
        U(1,J,1) =IU0(1,J) *U(1,J,1)
        U(I1,J,1)=IU0(I1,J)*U(I1,J,1)
      END IF
    END DO
    DO J=2,J2               !(I&J:CELL, JJ;INTERIOR FACE, JJ=J HERE)
      DO I=2,I1
        V(I,J,1)=IV0(I,J)*V(I,J,1)
      END DO
    END DO

!!!!hyedit 
    U = 1.0
    V = 0.0
    W(:,:,K0) = 0.0

#ifdef FLAG_FS
    DO K=K1,1,-1
      TEMP1=1./ODZ(K)
      DO J=2,J1
        TEMP=OCS(J)*ODY(J)
        DO I=2,I1
          W(I,J,K)=W(I,J,K+1)+((U(I,J,K)-U(I-1,J,K))*ODX(J)&
                    +(CSV(J)*V(I,J,K)-CSV(J-1)*V(I,J-1,K))*TEMP)*TEMP1
        END DO
      END DO
    END DO
    DO J=2,J1; JJ=J-1;      !(J:CELL, JJ:SIZE OF S, EASY TO READ)
      DO I=2,I1; II=I-1;    !(I:CELL, II:SIZE OF S, EASY TO READ)
        S(II,JJ)=W(I,J,1)
      END DO
    END DO
#else
    DO K=1,K1
      TEMP1=1./ODZ(K)
      DO J=2,J1
        TEMP=OCS(J)*ODY(J)
        DO I=2,I1
          W(I,J,K+1)=W(I,J,K)-((U(I,J,K)-U(I-1,J,K))*ODX(J)&
                    +(CSV(J)*V(I,J,K)-CSV(J-1)*V(I,J-1,K))*TEMP)*TEMP1
        END DO
      END DO
    END DO
    DO J=2,J1; JJ=J-1;      !(J:CELL, JJ:SIZE OF S, EASY TO READ)
      DO I=2,I1; II=I-1;    !(I:CELL, II:SIZE OF S, EASY TO READ)
        S(II,JJ)=-W(I,J,KB(I,J)+1)
      END DO
    END DO
#endif

    ! SOLVE SYSTEM MATRIX BY EVP METHOD
    !CALL REPBIR(AL,AB,AC,AR,AT,RINV,RINV1,DUM0,DUM1,DUM2,S,H,X,IE,I0,I2,I0,I2,NB0)
    CALL FS_SOLVER 
    ! RIGID-LID PRESSURE ADJUSTMENT
    DO J=2,J1               !(J&I:CELL)
      DO I=2,I1
        P0(I,J)=P0(I,J)+ODT*X(I,J)
      END DO
    END DO
    ! BT MODE ADJUSTMENT (FROM PRESSURE GRADIENT CHANGE)
    DO J=2,J1 ! LON.-DIR. FOR U COMPONENT (J&I:CELL, II:FACE, II=I HERE)
      DO I=2,I2             !(I:INTERIOR FACES, BC FACES(1&I1) = 0 OR PERIODIC) 
        SCR(I,J,1)=-(X(I+1,J)-X(I,J))*ODX(J) ! CURRENT CORRECTION
      END DO                
      SCR(1,J,1)=0; SCR(I1,J,1)=0                ! IF LOPENE /= 1 
!!      IF (LOPENE == 1) THEN ! LOPENW = LOPENE =1 ! PERIODIC B.C.
#ifdef FLAG_PERIODIC_WE
        SCR(1,J,1) =-(X(2,J)-X(I1,J))*ODX(J)  
        SCR(I1,J,1)=SCR(1,J,1)
#endif
!!      END IF  
      DO I=1,I1
        SCR(I,J,3)=SCR(I,J,3)+SCR(I,J,1)      ! OVERALL CORRECTIONS
        U(I,J,1)=U(I,J,1)+SCR(I,J,1)          ! TOP LAYER
        DO K=2,K1                             ! REMAINING LAYERS
        U(I,J,K)=U(I,J,K)+SCR(I,J,1)*IU(I,J,K)
        END DO
      END DO 
    END DO
    DO J=2,J2 ! LAT.-DIR. FOR FOR V COMPONENT (I&J:CELL, JJ:FACE, JJ=J HERE)
      DO I=2,I1             !(J:INTERIOR FACES, BC FACES(1&J1) = 0 NO PERIODIC)
        SCR(I,J,2)=-(X(I,J+1)-X(I,J))*ODYV(J) ! CURRENT CORRECTION
      END DO                
    END DO
    DO I=2,I1
      SCR(I,1,2)=0; SCR(I,J1,2)=0
    END DO
    DO J=1,J1
      DO I=2,I1
        SCR(I,J,4)=SCR(I,J,4)+SCR(I,J,2)      ! OVERALL CORRECTIONS
        V(I,J,1)=V(I,J,1)+SCR(I,J,2)          ! TOP LAYER
      END DO
    END DO
    DO K=2,K1                             ! REMAINING LAYERS
      DO J=1,J1
        DO I=2,I1
          V(I,J,K)=V(I,J,K)+SCR(I,J,2)*IV(I,J,K)
        END DO
      END DO
    END DO
    ! CHECK CONVERGENCE TO ZERO FLOW OVER LAND
    TMP=0.
    DO J=2,J1               ! FOR U COMPONENT (J&I:CELL; II:FACE, II=I HERE)
      DO I=2,I2             !(I:INTERIOR FACES)
        TMP=MAX(TMP,(1-IU0(I,J))*ABS(U(I,J,1)))
      END DO
    END DO
    DO J=2,J2               ! FOR V COMPONENT (I&J:CELL; JJ:FACE, JJ=J HERE)
      DO I=2,I1             !(J:INTERIOR FACES)
        TMP=MAX(TMP,(1-IV0(I,J))*ABS(V(I,J,1)))
      END DO
    END DO

    IF ( (ITF_GENERAL .LT. TOPTS) .OR. (MOD(ITF_GENERAL,TOPTS) .EQ. 0) .OR. (ITFTO .LE. NITFTO) ) THEN
    !IF ( (ITF-IT0 .LT. ITFDAY) .OR. ( MOD(ITF,ITFDAY) .eq. 0) ) THEN
       WRITE(*,684)ITF-IT0,ITMASK,TMP
       WRITE(FNO_LOG,684)ITF-IT0,ITMASK,TMP
    ENDIF
!@    if(itf-it0.lt.itfday.or.mod(itf,itfday).eq.0) write(14,684) itf-it0,itmask,tmp
 684 FORMAT('@itf-it0=',i6,',itmask=',I2,', Vmx on land=',F9.5,' cm/s')

     IF ((TMP .LE. TOLERANCE) .OR. (ITMASK .EQ. MXMASK)) EXIT ! USE TOLERANCE VALUE
   END DO ITMASKLOOP
#ifdef _DEBUG_
  PRINT*,'Pass Step 14'
#endif
   ! ======================================== !
   ! 15. CALCULATE W BY CONTINUITY EQUATION   !
   ! ======================================== !
#ifdef FLAG_FS
   DO J=2,J1              !(J:CELL & NORTH FACE FOR V)
     TEMP=OCS(J)*ODY(J)
     DO I=2,I1            !(I:CELL & EAST  FACE FOR U)
       DO K=KB(I,J),1,-1
         TMP=1./ODZ(K)
         W(I,J,K)=(W(I,J,K+1)+((U(I,J,K)-U(I-1,J,K))*ODX(J)&
                   +(CSV(J)*V(I,J,K)-CSV(J-1)*V(I,J-1,K))*TEMP)*TMP)
       END DO
     END DO
   END DO
#else
   DO K=1,K1                !(K:LAYER& TOP   FACE FOR W, K+1:BOTTOM FACE FOR W)
     TMP=1./ODZ(K)
     DO J=2,J1              !(J:CELL & NORTH FACE FOR V)
       TEMP=OCS(J)*ODY(J)
       DO I=2,I1            !(I:CELL & EAST  FACE FOR U)
         W(I,J,K+1)=IW(I,J,K+1)*(W(I,J,K)-((U(I,J,K)-U(I-1,J,K))*ODX(J)&
                   +(CSV(J)*V(I,J,K)-CSV(J-1)*V(I,J-1,K))*TEMP)*TMP)
       END DO
     END DO
   END DO
#endif
   ! ============================================================================== !
   ! 16. INTERPOLATE INCOMPRESSIBILITY INDUCED CHANGES BACK TO CELL CENTERS (U2&V2) !
   ! ============================================================================== !
   DO J=2,J1 ! LON.-DIR. FOR U2 COMPONENT (J&I:CELL, I=2:WBC, I=I1:EBC)
     DO I=3,I2 ! 4TH-ORDER ENHANCEMENT (SCR(I,J,1):CELL, SCR(I,J,3):EAST FACE)
       SCR(I,J,1) =-SCR(I-2,J,3)+SCR(I-1,J,3)+SCR(I,J,3) -SCR(I+1,J,3) 
     END DO                 
     SCR(2,J,1)=0; SCR(I1,J,1)=0                  ! IF LOPENE /= 1
   END DO
!!   IF (LOPENE == 1) THEN  ! LOPENW = LOPENE = 1 ! PERIODIC B.C.
#ifdef FLAG_PERIODIC_WE
     DO J=2,J1
       SCR(2,J,1) =-SCR(I2,J,3) +SCR(1,J,3)  +SCR(2,J,3) -SCR(3,J,3) ! 0  = I2
       SCR(I1,J,1)=-SCR(I3,J,3) +SCR(I2,J,3) +SCR(I1,J,3)-SCR(2,J,3) ! I0 = 2
     END DO
#endif
!!   END IF 
   DO K=1,K1
     DO J=2,J1
       DO I=3,I2            !(I:INTERIOR CELLS, I=2:WBC, I=I1:EBC)
         TMP1=IU(I-1,J,K) 
         TMP2=IU(I,J,K) 
         TMP =IU(I-2,J,K)*IU(I-1,J,K)*IU(I,J,K)*IU(I+1,J,K)
         U2(I,J,K)=IN(I,J,K)*(U2(I,J,K) &
                  +O24*(12.*(TMP1*SCR(I-1,J,3)+TMP2*SCR(I,J,3))+TMP*SCR(I,J,1)))
       END DO
       TMP1=IU(1,J,K) 
       TMP2=IU(2,J,K)
       TMP=1 
       U2(2,J,K) =IN(2,J,K) *(U2(2,J,K)  &
                 +O24*(12.*(TMP1*SCR(1,J,3) +TMP2*SCR(2,J,3)) +TMP*SCR(2,J,1)))
       TMP1=IU(I2,J,K) 
       TMP2=IU(I1,J,K) 
       U2(I1,J,K)=IN(I1,J,K)*(U2(I1,J,K) &
                 +O24*(12.*(TMP1*SCR(I2,J,3)+TMP2*SCR(I1,J,3))+TMP*SCR(I1,J,1)))
     END DO
   END DO
   IF ((TRIM(CASE_NAME) .EQ. 'NPBTAI') .AND. (LOPENW .GE. 2)) THEN
     DO K=1,K1
       DO J=2,J1  
         U2(2,J,K)=.5*(U(1,J,K)+U(2,J,K))    !SPECIAL TREATMENT FOR NPBTAI COUPLING
       END DO
     END DO
   END IF
   DO J=3,J2 ! 4TH-ORDER ENHANCEMENT (SCR(I,J,2):CELL, SCR(I,J,4):NORTH FACE)
     DO I=2,I1 ! LAT.-DIR. FOR V2 COMPONENT (J&I:CELL; J=2:SBC; J=J1:NBC)
       SCR(I,J,2) =-SCR(I,J-2,4)+SCR(I,J-1,4)+SCR(I,J,4)-SCR(I,J+1,4)
     END DO
   END DO
   DO I=2,I1
     SCR(I,2,2)=0; SCR(I,J1,2)=0                  ! NO PERIODIC B.C. IN LAT.-DIR.
   END DO
   DO K=1,K1
     DO J=3,J2            !(J:INTERIOR CELLS, J=2:SBC, J=J1:NBC)
       DO I=2,I1
         TMP1=IV(I,J-1,K)
         TMP2=IV(I,J,K)
         TMP=IV(I,J-2,K)*IV(I,J-1,K)*IV(I,J,K)*IV(I,J+1,K)
         V2(I,J,K)=IN(I,J,K)*(V2(I,J,K) &
                  +O24*(12.*(TMP1*SCR(I,J-1,4)+TMP2*SCR(I,J,4))+TMP*SCR(I,J,2)))
       END DO
     END DO
     DO I=2,I1 
       TMP1=IV(I,1,K)
       TMP2=IV(I,2,K)
       TMP=1
       V2(I,2,K) =IN(I,2,K) *(V2(I,2,K) &
                 +O24*(12.*(TMP1*SCR(I,1,4) +TMP2*SCR(I,2,4)) +TMP*SCR(I,2,2)))
       TMP1=IV(I,J2,K)
       TMP2=IV(I,J1,K)
       V2(I,J1,K)=IN(I,J1,K)*(V2(I,J1,K) &
                 +O24*(12.*(TMP1*SCR(I,J2,4)+TMP2*SCR(I,J1,4))+TMP*SCR(I,J1,2)))
     END DO
   END DO
   IF (TRIM(CASE_NAME) .EQ. 'DTRAC') THEN
     DO K=1,K1
       DO I=2,I1  
         V2(I,2,K) =.5*(V(I,1,K) +V(I,2,K))
         V2(I,J1,K)=.5*(V(I,J2,K)+V(I,J1,K))
       END DO
     END DO
   END IF
   ! ============================= !
   ! 17. CHECK INCOMPRESSIBILITY   !
   ! ============================= !
   IF ((ITF-IT0 .LE. ITFDAY) .OR. (MOD(ITF,ITFDAY) .EQ. 0) .OR. (ITFTO .LE. NITFTO) ) THEN
     TMP=0.; ERR=0.
     DO K=1,K1              !(K:LAYER& TOP   FACE FOR W, K+1:BOTTOM FACE FOR W)
       DO J=2,J1            !(J:CELL & NORTH FACE FOR V)
         DO I=2,I1          !(I:CELL & EAST  FACE FOR U)
           TMP1=(U(I,J,K)-U(I-1,J,K))*ODX(J)
           TMP2=(CSV(J)*V(I,J,K)-CSV(J-1)*V(I,J-1,K))*OCS(J)*ODY(J)
           TMP3=(W(I,J,K+1)-W(I,J,K))*ODZ(K)
           TMP=TMP+MAX(ABS(TMP1),ABS(TMP2),ABS(TMP3))*IN(I,J,K)
           ERR=ERR+ABS(TMP1+TMP2+TMP3)*IN(I,J,K)
         END DO
       END DO
     END DO
     ERR=ERR/TMP
     WRITE(*,7111) ERR
     7111  FORMAT(' *** NORMALIZED mean incompressibility error = ',1PE9.2)
   END IF
#ifdef _DEBUG_
  PRINT*,'Pass Step 15-17'
#endif

  ! ========================================================================== !
  ! 18. ELIMINATE ARBITRARY PRESSURE CONSTANT, MEAN SURFACE PRESSURE IS ZERO   !
  ! ========================================================================== !
  if (FL_ARBR_P0==1) then
    PSM=0.
    NSM=0.
    DO J=2,J1
      DO I=2,I1
        NSM=NSM+IN(I,J,1)
        PSM=PSM+P0(I,J)*IN(I,J,1)
      END DO
    END DO
    PSM=PSM/NSM
    DO J=2,J1
      DO I=2,I1
        P0(I,J)=(P0(I,J)-PSM)*IN(I,J,1)
      END DO
    END DO
    IF (TRIM(CASE_NAME) == 'NPBTAI') THEN ! NOTICE!!! CASE DEPENDENT TREATMENT 
      TMP=P0(2,127)
      DO J=2,J1
        DO I=2,I1
          P0(I,J)=P0(I,J)-TMP
        END DO
      END DO
    END IF
  end if

  ! ================================================ !
  ! 19. SET PERIODIC BOUNDARY CONDITIONS IF NEEDED   !
  ! ================================================ !
  IF (LOPENE == 1) THEN ! LOPENW = LOPENE = 1 ! PERIODIC B.C.
    DO K=1,K1
      DO J=2,J1
        ! WEST
        U2(1,J,K)=U2(I1,J,K)
        V2(1,J,K)=V2(I1,J,K)
#ifdef FLAG_TS_T
        T2(1,J,K)=T2(I1,J,K)
#endif
#ifdef FLAG_TS_S
        S2(1,J,K)=S2(I1,J,K)
#endif
#ifdef FLAG_TRACER
        C2(1,J,K)=C2(I1,J,K)
#endif
        ! EAST
        U2(I0,J,K)=U2(2,J,K)
        V2(I0,J,K)=V2(2,J,K)
#ifdef FLAG_TS_T
        T2(I0,J,K)=T2(2,J,K)
#endif
#ifdef FLAG_TS_S
        S2(I0,J,K)=S2(2,J,K)
#endif
#ifdef FLAG_TRACER
        C2(I0,J,K)=C2(2,J,K)
#endif
      END DO
    END DO
  END IF
  ! ============================== !
  ! 20. UPDATE USING FLTW METHOD   !
  ! ============================== !
  TMP=FLTW
  OFLTW=1.-2.*TMP

!  PRINT*,'FLTW,WRAF,U,V,T,S',ITF,FLTW,WRAF,SUM(U2),SUM(V2),SUM(T2),SUM(S2)
!  PRINT*,'U,V,T,S',MAXVAL(U2),MAXVAL(V2),MAXVAL(T2),MAXVAL(S2)
!  PAUSE


  DO K=1,K1
    DO J=1,J0
      DO I=1,I0
        TMP_U=FLTW*(.5*U1(I,J,K)-ULF(I,J,K)+.5*U2(I,J,K))
        U1(I,J,K)=ULF(I,J,K)+WRAF*TMP_U
        ULF(I,J,K)=U2(I,J,K)+(WRAF-1.)*TMP_U
        TMP_V=FLTW*(.5*V1(I,J,K)-VLF(I,J,K)+.5*V2(I,J,K))
        V1(I,J,K)=VLF(I,J,K)+WRAF*TMP_V
        VLF(I,J,K)=V2(I,J,K)+(WRAF-1.)*TMP_V
#ifdef FLAG_TS_T
        TMP_T=FLTW*(.5*T1(I,J,K)-TLF(I,J,K)+.5*T2(I,J,K))
        T1(I,J,K)=TLF(I,J,K)+WRAF*TMP_T
        TLF(I,J,K)=T2(I,J,K)+(WRAF-1.)*TMP_T
#endif
#ifdef FLAG_TS_S
        TMP_S=FLTW*(.5*S1(I,J,K)-SLF(I,J,K)+.5*S2(I,J,K))
        S1(I,J,K)=SLF(I,J,K)+WRAF*TMP_S
        SLF(I,J,K)=S2(I,J,K)+(WRAF-1.)*TMP_S
#endif
#ifdef FLAG_TRACER
        TMP_C=FLTW*(.5*C1(I,J,K)-CLF(I,J,K)+.5*C2(I,J,K))
        C1(I,J,K)=CLF(I,J,K)+WRAF*TMP_C
        CLF(I,J,K)=C2(I,J,K)+(WRAF-1.)*TMP_C
#endif
      END DO
    END DO
  END DO
  
#ifdef _DEBUG_
  PRINT*,'Pass Step 18-20'
#endif

  ! =============================== !
  ! 21. NUDGE TO CLIMATOLOGY DATA   !
  ! =============================== !
  IF (FL_NUDGE_OP_T .GE. 2) THEN
    DO K=2,K1 ! FROM K = 2 ! K = 1 (SURFACE FORCING HEAT FLUX) 
      summ=0.
      DO J=2,J1
        TMP1=-MIN(.1*(J-2)**2,10.)    
        TMP2=-MIN(.001*(J1-J)**2,10.) 
        !TMP1=TAUDTN*(EXP(TMP1)+EXP(TMP2))
        TMP3=EXP(TMP1)+EXP(TMP2)
        !TMP2=1.-TMP1
        AREA=DX(J)*DY(J)
        DO I=2,I1
          TMP1=TAU3D(i,j,k)*TMP3
          TMP2=1.-TMP1
          N=1-IN(I,J,K)*NUDGE(I,J,K)
          TEMP=T1(I,J,K)
          TMP=N*TEMP+IN(I,J,K)*NUDGE(I,J,K)*(TMP2*TEMP+TMP1*TCLI(I,J,K))
          EPS=TMP-TEMP
          T1(I,J,K)=TEMP+EPS
          T2(I,J,K)=T2(I,J,K)+EPS
          TLF(I,J,K)=TLF(I,J,K)+EPS
          summ=summ+AREA*EPS
        END DO
      END DO
      IF(FL_NUDGE_OP_T.GE.3)THEN
        if(A(K).EQ.0) then
            SUMM=0. 
        else 
            summ=summ/A(K)
        end if
        ! AVOID MEAN T CHANGES
        DO J=2,J1
          DO I=2,I1
            TEMP=IN(I,J,K)*summ
            T1(I,J,K)=T1(I,J,K)-TEMP
            T2(I,J,K)=T2(I,J,K)-TEMP
            TLF(I,J,K)=TLF(I,J,K)-TEMP
          END DO
        END DO
      ENDIF
    END DO
  END IF
  IF (FL_NUDGE_OP_S .GE. 2) THEN
    DO K=2,K1 ! FROM K = 2 ! K = 1 (SURFACE FORCING HEAT FLUX) 
      summ=0.
      DO J=2,J1
        TMP1=-MIN(.1*(J-2)**2,10.)    
        TMP2=-MIN(.001*(J1-J)**2,10.) 
        TMP3=EXP(TMP1)+EXP(TMP2)
        AREA=DX(J)*DY(J)
        DO I=2,I1
          tmp1=tau3d(i,j,k)*tmp3
          TMP2=1.-TMP1
          N=1-IN(I,J,K)*NUDGE(I,J,K)
          TEMP=S1(I,J,K)
          TMP=N*TEMP+IN(I,J,K)*NUDGE(I,J,K)*(TMP2*TEMP+TMP1*SCLI(I,J,K))
          EPS=TMP-TEMP
          S1(I,J,K)=TEMP+EPS
          S2(I,J,K)=S2(I,J,K)+EPS
          SLF(I,J,K)=SLF(I,J,K)+EPS
          summ=summ+AREA*EPS
        END DO
      END DO
      IF(FL_NUDGE_OP_S.GE.3)THEN
        if(A(K).EQ.0) then
          SUMM=0.
        else
          summ=summ/A(K)
        end if
        ! AVOID MEAN S CHANGES
        DO J=2,J1
          DO I=2,I1
            TEMP=IN(I,J,K)*summ
            S1(I,J,K)=S1(I,J,K)-TEMP
            S2(I,J,K)=S2(I,J,K)-TEMP
            SLF(I,J,K)=SLF(I,J,K)-TEMP
          END DO
        END DO
      ENDIF
    END DO
  END IF

  ! ============================================================================ !
  ! 22. RELAX SWAMP LAYER & MASKED WATER TO INTERIOR TO REDUCE SEEPAGE EFFECTS   !    
  ! ============================================================================ !
!  IF (FL_SWAMP_RE == 1) THEN
!    DO J=2,J1
!      T1(1,J,1) =T1(2,J,1)  ! JUST FOR SMOOTHING T1(1,J,1) = T1(2,J,1) IS O.K. 
!      T1(I0,J,1)=T1(I1,J,1)
!      S1(1,J,1) =S1(2,J,1)
!      S1(I0,J,1)=S1(I1,J,1)
!    END DO
!    DO I=2,I1
!      T1(I,1,1) =T1(I,2,1)
!      T1(I,J0,1)=T1(I,J1,1)
!      S1(I,1,1) =S1(I,2,1)
!      S1(I,J0,1)=S1(I,J1,1)
!    END DO
!    DO J=2,J1
!      DO I=2,I1
!        TMP=.25*(1-IN(I,J,1))
!        T1(I,J,1)=IN(I,J,1)*T1(I,J,1)+TMP*&
!                 (T1(I-1,J,1)+T1(I+1,J,1)+T1(I,J+1,1)+T1(I,J-1,1))
!        S1(I,J,1)=IN(I,J,1)*S1(I,J,1)+TMP*&
!                 (S1(I-1,J,1)+S1(I+1,J,1)+S1(I,J+1,1)+S1(I,J-1,1))
!      END DO
!    END DO
!  END IF

  IF (FL_SWAMP_RE == 1) THEN
    IF (LOPENE == 1) THEN ! LOPENW = LOPENE = 1 ! PERIODIC
      DO J=2,J1
        T1(1,J,1) =T1(I1,J,1)
        T1(I0,J,1)=T1(2,J,1)
        S1(1,J,1) =S1(I1,J,1)
        S1(I0,J,1)=S1(2,J,1)
      END DO
    ELSE
      DO J=2,J1
        T1(1,J,1) =T1(2,J,1)  ! JUST FOR SMOOTHING T1(1,J,1) = T1(2,J,1) IS O.K.
        T1(I0,J,1)=T1(I1,J,1)
        S1(1,J,1) =S1(2,J,1)
        S1(I0,J,1)=S1(I1,J,1)
      END DO
    ENDIF
    DO I=2,I1
      T1(I,1,1) =T1(I,2,1)
      T1(I,J0,1)=T1(I,J1,1)
      S1(I,1,1) =S1(I,2,1)
      S1(I,J0,1)=S1(I,J1,1)
    END DO
    DO J=2,J1
      DO I=2,I1
        if (in(i,j,1).eq.0) then
          TMP=IN(I-1,J,1)+IN(I+1,J,1)+IN(I,J+1,1)+IN(I,J-1,1) ! USE ONLY THE WATER POINTS
          if (TMP.ne.0) then
            TMP=(1-IN(I,J,1))/TMP
            T1(I,J,1)=TMP*(IN(I-1,J,1)*T1(I-1,J,1)+IN(I+1,J,1)*T1(I+1,J,1)+IN(I,J+1,1)*T1(I,J+1,1)+IN(I,J-1,1)*T1(I,J-1,1))
            S1(I,J,1)=TMP*(IN(I-1,J,1)*S1(I-1,J,1)+IN(I+1,J,1)*S1(I+1,J,1)+IN(I,J+1,1)*S1(I,J+1,1)+IN(I,J-1,1)*S1(I,J-1,1))
          end if
        end if
      END DO
    END DO
    IF (LOPENE == 1) THEN ! LOPENW = LOPENE = 1 ! PERIODIC
      DO J=2,J1
        T1(1,J,1) =T1(I1,J,1)
        T1(I0,J,1)=T1(2,J,1)
        S1(1,J,1) =S1(I1,J,1)
        S1(I0,J,1)=S1(2,J,1)
      END DO
    ELSE
      DO J=2,J1
        T1(1,J,1) =T1(2,J,1)  ! JUST FOR SMOOTHING T1(1,J,1) = T1(2,J,1) IS O.K.
        T1(I0,J,1)=T1(I1,J,1)
        S1(1,J,1) =S1(2,J,1)
        S1(I0,J,1)=S1(I1,J,1)
      END DO
    ENDIF
  END IF 

  ! ======================== !
  ! 23. BIHARMONIC FILTERS   !
  ! ======================== !
  IF (FL_BI_FIL == 1) THEN
    K=MOD(ITF,2*K1)+1
    IF (K.GT.K1) K=2*K1-K+1
    BF=0.975**K
#ifdef FLAG_TS_T
    CALL BFLT_GLOBAL(BF,T1(1,1,K),K)
#endif
#ifdef FLAG_TS_S
    CALL BFLT_GLOBAL(BF,S1(1,1,K),K)
#endif
    BF=0.975
#ifdef FLAG_TS_T
    CALL BFLT_GLOBAL(BF,T1,1)
#endif
#ifdef FLAG_TS_S
    CALL BFLT_GLOBAL(BF,S1,1)
#endif
    IF (LOPENE == 1) THEN ! LOPENW = LOPENE = 1 ! PERIODIC
      DO J=2,J1
#ifdef FLAG_TS_T
        T1(1,J,K) =T1(I1,J,K)
        T1(I0,J,K)=T1(2,J,K)
        T1(1,J,1) =T1(I1,J,1)
        T1(I0,J,1)=T1(2,J,1)
#endif
#ifdef FLAG_TS_S
        S1(1,J,K) =S1(I1,J,K)
        S1(I0,J,K)=S1(2,J,K)
        S1(1,J,1) =S1(I1,J,1)
        S1(I0,J,1)=S1(2,J,1)
#endif
      END DO
    END IF
  END IF
#ifdef _DEBUG_
  PRINT*,'Complete Step 23'
#endif

  ! =============================================== !
  ! 24. LIMITER APPROACH FOR TRACER CONCENTRATION   !
  ! =============================================== !
#ifdef FLAG_TRACER
    if (mod(itf,itfday).ne.0) go to 841
    write(*,830) int(days+0.01)
830 format('day',i3,', before fct')
    k=kb(i0/2,j1) ! CHECK SILL LEVEL
    CALL RANGER(C2(1,1,K),IN,I0,2,2,I1,J1,ILO,JLO,IHI,JHI,CLO,CHI)
840 write(*,802) K,CLO,ILO,JLO,CHI,IHI,JHI
802 format('level',i3,' tracer: mn=',ES12.5,'@ (',I3,',',I3,'), mx=',ES12.5,'@ (',I3,',',I3,')')
841 DO 843 K=1,K1
    DO 843 J=2,J1
    DO 843 I=2,I1
      IFLT(I,J,K)=0
      TMP=(C2(I,J,K)-0.)*(100.-C2(I,J,K))
      IF (TMP.GE.0.) GO TO 843
      IFLT(I,J,K)=1
      IF (C2(I,J,K).LT.0.) IFLT(I,J,K)=-1
    843 CONTINUE  
    DO 850 K=1,K1
    DO 850 J=2,J1
    DO 850 I=2,I1
      IF (IFLT(I,J,K).EQ.0) GO TO 850
      IL=MAX(2,I-1); IH=MIN(I1,I+1)
      JL=MAX(2,J-1); JH=MIN(J1,J+1)
      KL=MAX(1,K-1); KH=MIN(K1,K+1)
      IF (IFLT(I,J,K).EQ.1) GO TO 847
      TMP=0. ! FOR C2.LT.0
      DO KK=KL,KH
        DO JJ=JL,JH
          DO II=IL,IH
            TMP=(1-IN(II,JJ,KK))*TMP+IN(II,JJ,KK)*MAX(TMP,C2(II,JJ,KK))
          END DO
        END DO
      END DO
      DO 846 KK=KL,KH
      DO 846 JJ=JL,JH
      DO 846 II=IL,IH
        IF (C2(II,JJ,KK).NE.TMP) GO TO 846
        CLF(II,JJ,KK)=C2(II,JJ,KK)+C2(I,J,K)
        CLF(I,J,K)=0.
        GO TO 850
      846 CONTINUE
      GO TO 850
      847  TMP=100.! FOR C2.GT.100 
      DO KK=KL,KH
        DO JJ=JL,JH
          DO II=IL,IH
            TMP=(1-IN(II,JJ,KK))*TMP+IN(II,JJ,KK)*MIN(TMP,C2(II,JJ,KK))
          END DO
        END DO
      END DO
      DO 849 KK=KL,KH
      DO 849 JJ=JL,JH
      DO 849 II=IL,IH                                        
        IF (C2(II,JJ,KK).NE.TMP) GO TO 849
        CLF(II,JJ,KK)=C2(II,JJ,KK)+C2(I,J,K)-100.
        CLF(I,J,K)=100.
        GO TO 850
      849 CONTINUE
    850 CONTINUE
    DO K=1,K1
      DO J=2,J1
        DO I=2,I1
          C2(I,J,K)=CLF(I,J,K)
        END DO
      END DO
    END DO
    if (mod(itf,itfday).ne.0) return
    write(*,930) days
930 format('after fct',f7.3)
    k=kb(i0/2,j1) ! CHECK SILL LEVEL
    CALL RANGER(C2(1,1,K),IN,I0,2,2,I1,J1,ILO,JLO,IHI,JHI,CLO,CHI)
940 write(*,802) k,clo,ilo,jlo,chi,ihi,jhi 
#endif

  CONTAINS

  SUBROUTINE FS_SOLVER
    ! local variables for compare bicg and EVP 
    real, parameter  :: pi=3.1415926535897932384626433832795029
    real*8 :: err1,err2
    real*8 :: mineig, maxeig
    real*8 :: arbit, arbit1
    real*8, dimension(i0,j0) :: x1,x0

    DO J=1,J2+2
      DO I=1,I2+2
        !X(I,J)= cos((I-2)/real(I2)*2.0*pi ) +1.0 ! TS
        !X(I,J)= cos((I-2)/real(I2)*2.0*pi )
        X(I,J)= 0.0
        !X(I,J)=.5*sin((I-2)/real(I2)*2.0*pi)*cos((J-2)/real(I2))
      END DO
    END DO
    x0 = x
    write(*,*) "==============intial value========="
    write(*,*) X(:,1)
    write(*,*) "==============intial value========="
    write(*,'(3e13.3)') X(1,10),X(2,10),X(3,10)
    write(*,'(3e13.3)') X(I2,10),X(I2+1,10),X(I2+2,10)
    !X = 1.0
    !AL = -0.25
    !AB = -0.25
    !AC =  1.0
    !AR = -0.25
    !AT = -0.25

    S(:,:) = 0.0
    DO J=2,J2+1
      DO I=2,I2+1
        S(I-1,J-1)=AL(I-1,J-1)*X(I-1,J)+AB(I-1,J-1)*X(I,J-1)+AC(I-1,J-1)*X(I,J)+AR(I-1,J-1)*X(I+1,J)+AT(I-1,J-1)*X(I,J+1)
      END DO
    END DO

    DO J=3,J2-1
      DO I=3,I2-1
        X(I,J) = 1.0
      END DO
    END DO
    !!X(3:I2-1,3:J2-1) = 1.0
    X1(:,:) = X(:,:)

    IF (FL_EVP_STP == 0) THEN ! DTRAC, NPBTAI
      !CALL REPBIR_PERIOD(AL,AB,AC,AR,AT,RINV,RINV1,DUM0,DUM1,DUM2,S,H,X,IE,I0,I2,I0,I2,NB0)        
      CALL REP1(AL,AB,AC,AR,AT,RINV,RINV1,DUM0,DUM1,DUM2,S,H,X,IE,I0,I2,NB0)        
      call p_bicgstab_le(al,ab,ac,ar,at,s,x1,al,ab,ac,ar,at,i2,j2)
     
      arbit  = (x (10,10)-x0(10,10))
      arbit1 = (x1(10,10)-x0(10,10))
      x(:,2:j2+1) = x(:,2:j2+1) -arbit
      x1(:,2:j2+1) = x1(:,2:j2+1) -arbit1
      !do j = 2,j2+1
      !  do i = 2,i2+1
      !    err1  = x1(i,j) - x(i,j) 
      !    err2  = min(abs(x1(i,j)),abs(x(i,j)))
      !    if ( abs(err1) .ge. 1.0 .and. err1/err2 .ge. 1.0e-1 ) then
      !    ! write(*,'(2I5,5e13.3)') i,j,x1(i,j),x(i,j),oac(i-1,j-1),err1,err1/err2
      !     write(*,'(2I5,5e13.3)') i,j,x1(i,j),x(i,j)
      !    end if 
      !  end do 
      !end do 
      do j = 1,j2+2
        do i = 1,i2+2
          err1  = x1(i,j) - x0(i,j) 
          if ( abs(err1) .ge. 1.0e-3 ) then
           write(*,'(A5,2I5,5e13.3)') "bicg0 ", i,j,x0(i,j),x1(i,j)
          end if 
        end do 
      end do 
      do j = 1,j2+2
        do i = 1,i2+2
          err1  = x(i,j) - x0(i,j) 
          if ( abs(err1) .ge. 1.0e-3 ) then
           write(*,'(A5,2I5,5e13.3)') "evp0  ", i,j,x0(i,j), x(i,j)
          end if 
        end do 
      end do 
      err1 = maxval(abs(x1-x))
      err2 = norm2(x1-x)
      
      write(*,'(A20,3f20.11)') "EVP-BICG PointErr NormErr NormX ",err1,err2,norm2(x)
      write(*,'(A20,2I)') "Max error at", maxloc(abs(x1-x))
      write(*,'(A20,2I)') "Dimensions", I2,J2

      open(unit=110,file="evp-solution.bin",action="write",form="unformatted")
      write(110) x  -x0
      close(110)
      open(unit=110,file="bicg-solution.bin",action="write",form="unformatted")
      write(110) x1 -x0
      close(110)
      open(unit=110,file="evp-bicg-sol.bin",action="write",form="unformatted")
      write(110) x1-x
      close(110)
      stop
    ELSE IF (FL_EVP_STP == 1) THEN ! GLOBAL

      NITERLOOP: DO NITER=1,5 ! NITER IS THE BIR SWEEP NUMBER
        write(*,*) 'NITER =', NITER
        INC=IBIR/2
        DO N=1,NBIR ! N IS THE BIR LONGITUDINAL BAND #
          IL=(N-1)*INC+1
          IF (N.NE.NBIR) THEN
            IR=IL+IBIR-1
            DO J=1,J2
              DO I=1,IBIR
                SRC(I,J)=S(IL+I-1,J)
              END DO
            END DO
            CALL REPBIR(AL(IL,1),AB(IL,1),AC(IL,1),AR(IL,1),AT(IL,1),&
               RINV(1,1,N),RINV1(1,1,N),DUM0,DUM1,DUM2,SRC,H,X(IL,1),IE,&
               I0,I2,IBIR+2,IBIR,NB0)
          ELSE ! LAST BIR STRIP
            IR=IL+INC
            IR2=IL+IBIR-1
            DO J=1,J2
              XX(1,J)=X(IL,J)         ! SET SCRATCH ARRAY XX EAST & WEST BOUNDARIES
              XX(IBIR+2,J)=X(INC+2,J)
              DO I=IL,IR-1
                SRC(I-IL+1,J)=S (I,J)
                CL (I-IL+1,J)=AL(I,J)
                CR (I-IL+1,J)=AR(I,J)
                CB (I-IL+1,J)=AB(I,J)
                CT (I-IL+1,J)=AT(I,J)
                CC (I-IL+1,J)=AC(I,J)
                !OCC (I-IL+1,J)=OAC(I,J)
              END DO
              DO I=IR,IR2
                SRC(I-IL+1,J)=S (I-IR+1,J)
                CL (I-IL+1,J)=AL(I-IR+1,J)
                CR (I-IL+1,J)=AR(I-IR+1,J)
                CB (I-IL+1,J)=AB(I-IR+1,J)
                CT (I-IL+1,J)=AT(I-IR+1,J)
                CC (I-IL+1,J)=AC(I-IR+1,J)
                !OCC (I-IL+1,J)=OAC(I-IR+1,J)
              END DO
            END DO
            CALL REPBIR(CL,CB,CC,CR,CT,&
               RINV(1,1,N),RINV1(1,1,N),DUM0,DUM1,DUM2,SRC,H,XX,IE,&
               IBIR+2,IBIR,IBIR+2,IBIR,NB0)
            DO J=2,J1
              DO I=IL,IR-1
                X(I+1,J)=XX(I-IL+2,J)
              END DO
              DO I=IR,IR2
                X(I-IR+2,J)=XX(I-IL+2,J)
              END DO
            END DO
!!            IF (LOPENE == 1) THEN ! LOPENW = LOPENE = 1 ! PERIODIC
#ifdef FLAG_PERIODIC_WE
              DO J=2,J1
                X(1,J)=X(I1,J)
                X(I0,J)=X(2,J)
              END DO
              WRITE(*,*)'TEST',NITER,MAXVAL(abs(X(2:i2+1,2:j2+1)-X1(2:i2+1,2:j2+1)))
#endif
!!            END IF
          END IF
        END DO
      END DO NITERLOOP      

      call p_bicgstab_le(al,ab,ac,ar,at,s,x1,al,ab,ac,ar,at,i2,j2)
      arbit  = (x (10,10)-x0(10,10))
      arbit1 = (x1(10,10)-x0(10,10))
      !x(:,2:j2+1) = x(:,2:j2+1) -arbit
      x1(:,2:j2+1) = x1(:,2:j2+1) -arbit1


      !do j = 2,j2+1
      !  do i = 2,i2+1
      !    err1  = x1(i,j) - x(i,j) 
      !    err2  = max(abs(x1(i,j)),abs(x(i,j)))
      !    if ( abs(err1) .ge. 1.0e-3 ) then
      !     !write(*,'(2I5,5e13.3)') i,j,x1(i,j),x(i,j),oac(i-1,j-1),err1,err1/err2
      !     write(*,'(2I5,5e13.3)') i,j,x1(i,j),x(i,j)
      !    end if 
      !  end do 
      !end do 
      !err1 = maxval(abs(x1-x))
      !err2 = norm2(x1-x)
      
      do j = 1,j2+2
        do i = 1,i2+2
          err1  = x1(i,j) - x0(i,j) 
          if ( abs(err1) .ge. 1.0e-3 ) then
           write(*,'(A5,2I5,5e13.3)') "bicg1 ", i,j,x0(i,j),x1(i,j)
          end if 
        end do 
      end do 
      do j = 1,j2+2
        do i = 1,i2+2
          err1  = x(i,j) - x0(i,j) 
          if ( abs(err1) .ge. 1.0e-3 ) then
           write(*,'(A5,2I5,5e13.3)') "evp1  ", i,j,x0(i,j), x(i,j)
          end if 
        end do 
      end do 
      err1 = maxval(abs(x1-x))
      err2 = norm2(x1-x)
      
      write(*,'(A20,3f20.11)') "EVP-BICG PointErr NormErr NormX ",err1,err2,norm2(x)
      write(*,'(A20,2I)') "Max error at", maxloc(abs(x1-x))
      write(*,'(A20,2I)') "Dimensions", I2,J2

      open(unit=110,file="evp-solution.bin",action="write",form="unformatted")
      write(110) x  -x0
      close(110)
      open(unit=110,file="bicg-solution.bin",action="write",form="unformatted")
      write(110) x1 -x0
      close(110)
      open(unit=110,file="evp-bicg-sol.bin",action="write",form="unformatted")
      write(110) x1-x
      close(110)
      stop


    else if (FL_EVP_STP == 2 ) then ! call bicg
      call p_bicgstab_le(al,ab,ac,ar,at,s,x,al,ab,ac,ar,at,i2,j2)
      call lanczos(al,ab,ac,ar,at,mineig,maxeig,i2,j2)
      call p_csi(al,ab,ac,ar,at,mineig,maxeig,s,x1,al,ab,ac,ar,at,i2,j2)
      do j = 2,j2+1
        do i = 2,i2+1
          err1  = x1(i,j) - x(i,j) 
          err2  = max(abs(x1(i,j)),abs(x(i,j)))
          if ( abs(err1) .ge. 1.0e-3 ) then
           !write(*,'(2I5,5e13.3)') i,j,x1(i,j),x(i,j),oac(i-1,j-1),err1,err1/err2
            write(*,'(2I5,5e13.3)') i,j,x1(i,j),x(i,j)
          end if 
        end do 
      end do 
      err1 = maxval(abs(x1-x))
      err2 = norm2(x1-x)

      write(*,'(A20,3e15.8)') "BICG-TURE Max Norm2",err1,err2,norm2(x)
      write(*,'(A20,2I)') "Max error at", maxloc(abs(x1-x))
      write(*,'(A20,2I)') "Dimensions", I2,J2
      open(unit=110,file="pcsi-bicg.bin",action="write",form="unformatted")
      write(110) x1-x
      close(110)

      stop

    END IF    
  END SUBROUTINE FS_SOLVER

INCLUDE './turbulence.f90'

END SUBROUTINE FS_GLOBAL

! ----------------------------------------------------------------------

SUBROUTINE BFLT_GLOBAL(PWT,P,K)
USE OCN_PARA_GLOBAL,ONLY:I0,J0,I1,J1,I2,J2
USE GRID_VAR_GLOBAL, ONLY: INFX,INFY,VAR_ALLOCATE_WINDMX_MAIN
IMPLICIT NONE
REAL,DIMENSION(I0,J0),INTENT(INOUT)::P
REAL,DIMENSION(I0,J0)::S
REAL::PWT,TMP
INTEGER::I,J,K
  TMP=(1.-PWT)/16.
  S=0.0
  DO J=2,J1
    DO I=2,I1
      S(I,J)=P(I-1,J)-2.*P(I,J)+P(I+1,J)
    END DO
    DO I=2,I1
      P(I,J)=P(I,J)-TMP*INFX(I,J,K)*(S(I-1,J)+S(I+1,J)-2.*S(I,J))
    END DO
  END DO
  DO I=2,I1
    DO J=3,J2
      S(I,J)=P(I,J-1)-2.*P(I,J)+P(I,J+1)
    END DO
    DO J=4,J0-3
      P(I,J)=P(I,J)-TMP*INFY(I,J,K)*(S(I,J-1)+S(I,J+1)-2.*S(I,J))
    END DO
  END DO
END SUBROUTINE BFLT_GLOBAL

