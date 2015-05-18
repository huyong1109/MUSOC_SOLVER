MODULE density_EOS_GLOBAL
 IMPLICIT NONE

 CONTAINS

SUBROUTINE EOS_COMPUTE(T,S,Z,RHO,FL_EOS_OP)
 USE OCN_PARA_GLOBAL
 USE GRID_VAR_GLOBAL,ONLY: KB,IU0,IV0,IN,IU,IV,IW
 USE GRID_VAR_GLOBAL, ONLY: Mcrho,McTP,McS,McTPTP,McSTP
 IMPLICIT NONE
 REAL,INTENT(IN) :: T(:,:,:),S(:,:,:)
 REAL,INTENT(IN) :: Z(:)
 REAL(8),INTENT(INOUT) :: RHO(:,:,:)
 INTEGER,INTENT(IN) :: FL_EOS_OP
 REAL(8) :: TMPD,SLTD,PISD,p02,p01,ccrho,ccTP,ccS,ccTPTP,ccSTP,TEMP,TMP1,TMP2,TMP3
 INTEGER :: I,J,K

 TEMP=1.D5 
 IF (FL_EOS_OP.EQ.1) THEN
  ! LINEAR EOS
   DO K=1,K1
     DO J=2,J1
       DO I=2,I1
         TMP1=DBLE(IN(I,J,K))
         TMP2=1.D+0-TMP1
         RHO(I,J,K)=2.D-4*(1.D+1-DBLE(T(I,J,K)))
         TEMP=TMP2*TEMP+TMP1*MIN(TEMP,RHO(I,J,K))
       END DO
     END DO
   END DO
 ELSE IF (FL_EOS_OP.EQ.2) THEN
 ! DAN WRIGHT'S FULL EOS, RHO(I,J,K)=r(TMPD,SLTD,PISD)
 ! Yu-Chiao 20110919
!   DO K=1,K1
!     DO J=2,J1
!       DO I=2,I1
!         TMP1=DBLE(IN(I,J,K))
!         TMP2=1.D+0-TMP1
!         TMPD=DBLE(MAX(-1.8,T(I,J,K)))
!         SLTD=DBLE(MIN(41.,S(I,J,K)))
!         SLTD=MAX(0.D0,SLTD)
!         PISD=1.D-3*DBLE(Z(2*K))
!         RHO(I,J,K)=Wright_density(TMPD,SLTD,PISD)
!         TEMP=TMP2*TEMP+TMP1*MIN(TEMP,RHO(I,J,K))
!       END DO
!     END DO
!   END DO
   DO K=1,K1
     DO J=2,J1
       DO I=2,I1
         TMPD=MAX(-1.8,T(I,J,K))
         SLTD=MIN(41.,S(I,J,K))
         SLTD=MAX(0.D0,SLTD)
         PISD=.001*Z(2*K)
         RHO(I,J,K)=Wright_density(TMPD,SLTD,PISD)
         TEMP=(1-IN(I,J,K))*TEMP+IN(I,J,K)*MIN(TEMP,RHO(I,J,K))
       END DO
     END DO
   END DO
 ELSE IF (FL_EOS_OP.EQ.3) THEN
 ! UNESCO
   DO K=1,K1
     DO J=2,J1
       DO I=2,I1
         TMP1=DBLE(IN(I,J,K))
         TMP2=1.D+0-TMP1
         TMPD=DBLE(MAX(-1.8,T(I,J,K)))
         SLTD=DBLE(MIN(41.,S(I,J,K)))
         SLTD=MAX(0.D0,SLTD)
         PISD=1.D-3*DBLE(Z(2*K))
         TMPD=temperature(TMPD,SLTD,PISD)
         RHO(I,J,K)=densityTSP(TMPD,SLTD,PISD)*1.D-3
         TEMP=TMP2*TEMP+TMP1*MIN(TEMP,RHO(I,J,K))
       END DO
     END DO
   END DO
 ELSE IF (FL_EOS_OP.EQ.4) THEN
 ! local linear scheme
   DO K=1,K1
     DO J=2,J1
       DO I=2,I1
         TMP1=DBLE(IN(I,J,K))
         TMP2=1.D+0-TMP1
         TMPD=DBLE(MAX(-1.8,T(I,J,K)))
         SLTD=DBLE(MIN(41.,S(I,J,K)))
         SLTD=MAX(0.D0,SLTD)
         PISD=1.D-3*DBLE(Z(2*K))
         ccrho=Mcrho(I,J,K)
         ccTP=McTP(I,J,K)
         ccS=McS(I,J,K)
         RHO(I,J,K)=local_linear_density(TMPD,SLTD,ccrho,ccTP,ccS)*1.D-3
         TEMP=TMP2*TEMP+TMP1*MIN(TEMP,RHO(I,J,K))
       END DO
     END DO
   END DO
 ELSE IF (FL_EOS_OP.EQ.5) THEN
 ! local nonlinear scheme 1
   DO K=1,K1
     DO J=2,J1
       DO I=2,I1
         TMP1=DBLE(IN(I,J,K))
         TMP2=1.D+0-TMP1 
         TMPD=DBLE(MAX(-1.8,T(I,J,K)))
         SLTD=DBLE(MIN(41.,S(I,J,K)))
         SLTD=MAX(0.D0,SLTD)
         PISD=1.D-3*DBLE(Z(2*K))
         ccrho=Mcrho(I,J,K)
         ccTP=McTP(I,J,K)
         ccS=McS(I,J,K)
         ccTPTP=McTPTP(I,J,K)
         RHO(I,J,K)=local_density2(TMPD,SLTD,ccrho,ccTP,ccS,ccTPTP)*1.D-3
         TEMP=TMP2*TEMP+TMP1*MIN(TEMP,RHO(I,J,K))
       END DO
     END DO
   END DO
 ELSE IF (FL_EOS_OP.EQ.6) THEN
 ! local nonlinear scheme 2
   DO K=1,K1
     DO J=2,J1
       DO I=2,I1
         TMP1=DBLE(IN(I,J,K))
         TMP2=1.D+0-TMP1
         TMPD=DBLE(MAX(-1.8,T(I,J,K)))
         SLTD=DBLE(MIN(41.,S(I,J,K)))
         SLTD=MAX(0.D0,SLTD)
         PISD=1.D-3*DBLE(Z(2*K))
         ccrho=Mcrho(I,J,K)
         ccTP=McTP(I,J,K)
         ccS=McS(I,J,K)
         ccTPTP=McTPTP(I,J,K)
         ccSTP=McSTP(I,J,K)
         RHO(I,J,K)=local_density3(TMPD,SLTD,ccrho,ccTP,ccS,ccTPTP,ccSTP)*1.D-3
         TEMP=TMP2*TEMP+TMP1*MIN(TEMP,RHO(I,J,K))
       END DO
     END DO
   END DO
 END IF

 DO K=1,K1
   DO J=2,J1
     DO I=2,I1
       RHO(I,J,K)=RHO(I,J,K)-TEMP
     END DO
   END DO
 END DO

END SUBROUTINE EOS_COMPUTE


subroutine EOS_COMPUTATION_STABILITY1(T,S,ZZ,RHO,FL_EOS_OP,FL_LEV_OP,EOS_FILE)
 USE OCN_PARA_GLOBAL !,ONLY: I0,I1,I2,J2
 USE GRID_VAR_GLOBAL, ONLY: KB,IU0,IV0,IN,IU,IV,IW,VAR_ALLOCATE_BATHY
 IMPLICIT NONE
 REAL,INTENT(INOUT),DIMENSION(:,:,:) :: T,S
 REAL,INTENT(IN),DIMENSION(:) :: ZZ
 REAL(8),INTENT(OUT),DIMENSION(:,:,:) :: RHO
 INTEGER,INTENT(IN) :: FL_EOS_OP,FL_LEV_OP,EOS_FILE
 INTEGER :: I,J,K,N,Ia1,Ia2,Ib1,Ib2,Ja1,Ja2,Jb1,Jb2,Ka1,Ka2,Kb1,Kb2,FNO_eos
 REAL(8) :: ccrho,ccTP,ccS,ccTPTP,ccSTP,TMPD,SLTD,PISD
 REAL(8),ALLOCATABLE :: Mcrho(:,:,:),McTP(:,:,:),McS(:,:,:), McTPTP(:,:,:),McSTP(:,:,:)
 REAL,DIMENSION(2*SIZE(ZZ)) :: Z
 CHARACTER(160) :: TEMPDIR
 INTEGER :: ISTAT
 
 ALLOCATE(Mcrho(I0,J0,K1),STAT=ISTAT) ; Mcrho=0.
 ALLOCATE(McTP(I0,J0,K1),STAT=ISTAT) ; McTP=0.
 ALLOCATE(McS(I0,J0,K1),STAT=ISTAT) ; McS=0.
 ALLOCATE(McTPTP(I0,J0,K1),STAT=ISTAT)	; McTPTP=0.
 ALLOCATE(McSTP(I0,J0,K1),STAT=ISTAT) ; McSTP=0.

 IF (EOS_FILE.EQ.0) THEN
 IF (FL_LEV_OP.EQ.1) THEN
   Ka1=1;Ka2=33;Kb1=2;Kb2=33
   Ja1=1;Ja2=180;Jb1=1;Jb2=180
   Ia1=1;Ia2=360;Ib1=1;Ib2=360
   DO K=1,SIZE(ZZ)
     Z(2*K)=ZZ(K)*100.
   END DO
 ELSE
   Ka1=1;Ka2=K1;Kb1=2;Kb2=K1
   Ja1=1;Ja2=J1;Jb1=2;Jb2=J1
   Ia1=2;Ia2=I0;Ib1=2;Ib2=I1
   DO K=1,SIZE(ZZ)
     Z(K)=ZZ(K)
   END DO
 END IF
 IF (FL_EOS_OP.EQ.1) THEN
   DO K=Ka1,Ka2
     DO J=Ja1,Ja2
       DO I=Ia1,Ia2
         TMPD=DBLE(T(I,J,K))
         SLTD=DBLE(S(I,J,K))
         PISD=1.D-3*DBLE(Z(2*K))
         RHO(I,J,K)=2.D-4*(1.D+1-TMPD)
       END DO
     END DO
   END DO
   DO K=Kb1,Kb2
     N=0
     DO J=Jb1,Jb2
       DO I=Ib1,Ib2
         IF (RHO(I,J,K)<RHO(I,J,K-1)) THEN
           N=N+1
           T(I,J,K)=T(I,J,K-1)
           S(I,J,K)=S(I,J,K-1)
           TMPD=DBLE(T(I,J,K))
           SLTD=DBLE(S(I,J,K))
           PISD=1.D-3*DBLE(Z(2*K))
           RHO(I,J,K)=2.D-4*(1.D+1-TMPD)
         END IF
       END DO
     END DO
     WRITE(*,58) N,K
   END DO
 ELSE IF (FL_EOS_OP.EQ.2) THEN
   DO K=Ka1,Ka2
     DO J=Ja1,Ja2
       DO I=Ia1,Ia2
         TMPD=T(I,J,K)
         SLTD=S(I,J,K)
         PISD=1.D-3*DBLE(Z(2*K))
         RHO(I,J,K)=Wright_density(TMPD,SLTD,PISD)
       END DO
     END DO
   END DO
   DO K=Kb1,Kb2
     N=0
     DO J=Jb1,Jb2
       DO I=Ib1,Ib2
         IF (RHO(I,J,K)<RHO(I,J,K-1)) THEN
           N=N+1
           T(I,J,K)=T(I,J,K-1)
           S(I,J,K)=S(I,J,K-1)
           TMPD=DBLE(T(I,J,K))
           SLTD=DBLE(S(I,J,K))
           PISD=1.D-3*DBLE(Z(2*K))
           RHO(I,J,K)=Wright_density(TMPD,SLTD,PISD)
         END IF
       END DO
     END DO
     WRITE(*,58) N,K
   END DO
 ELSE IF (FL_EOS_OP.EQ.3) THEN
   DO K=Ka1,Ka2
     DO J=Ja1,Ja2
       DO I=Ia1,Ia2
         TMPD=DBLE(T(I,J,K))
         SLTD=DBLE(S(I,J,K))
         PISD=1.D-3*DBLE(Z(2*K))
         TMPD=temperature(TMPD,SLTD,PISD)
         RHO(I,J,K)=densityTSP(TMPD,SLTD,PISD)*1.D-3
       END DO
     END DO
   END DO
   DO K=Kb1,Kb2
     N=0
     DO J=Jb1,Jb2
       DO I=Ib1,Ib2
         IF (RHO(I,J,K)<RHO(I,J,K-1)) THEN
           N=N+1
           T(I,J,K)=T(I,J,K-1)
           S(I,J,K)=S(I,J,K-1)
           TMPD=DBLE(T(I,J,K))
           SLTD=DBLE(S(I,J,K))
           PISD=1.D-3*DBLE(Z(2*K))
           TMPD=temperature(TMPD,SLTD,PISD)
           RHO(I,J,K)=densityTSP(TMPD,SLTD,PISD)*1.D-3
         END IF
       END DO
     END DO
     WRITE(*,58) N,K
   END DO
 ELSE IF (FL_EOS_OP.EQ.4) THEN
   DO K=Ka1,Ka2
     DO J=Ja1,Ja2
       DO I=Ia1,Ia2
         TMPD=DBLE(T(I,J,K))
         SLTD=DBLE(S(I,J,K))
         PISD=1.D-3*DBLE(Z(2*K))
         CALL make_local_coefficients(TMPD,SLTD,PISD,ccrho,ccTP,ccS)
         RHO(I,J,K)=local_linear_density(TMPD,SLTD,ccrho,ccTP,ccS)*1.D-3
       END DO
     END DO
   END DO
   DO K=Kb1,Kb2
     N=0
     DO J=Jb1,Jb2
       DO I=Ib1,Ib2
         IF (RHO(I,J,K)<RHO(I,J,K-1)) THEN
           N=N+1
           T(I,J,K)=T(I,J,K-1)
           S(I,J,K)=S(I,J,K-1)
           TMPD=DBLE(T(I,J,K))
           SLTD=DBLE(S(I,J,K))
           PISD=1.D-3*DBLE(Z(2*K))
           CALL make_local_coefficients(TMPD,SLTD,PISD,ccrho,ccTP,ccS)
           RHO(I,J,K)=local_linear_density(TMPD,SLTD,ccrho,ccTP,ccS)*1.D-3
         END IF
       END DO
     END DO
     WRITE(*,58) N,K
   END DO
 ELSE IF (FL_EOS_OP.EQ.5) THEN
   DO K=Ka1,Ka2
     DO J=Ja1,Ja2
       DO I=Ia1,Ia2
         TMPD=DBLE(T(I,J,K))
         SLTD=DBLE(S(I,J,K))
         PISD=1.D-3*DBLE(Z(2*K))
         CALL make_local_coefficients2(TMPD,SLTD,PISD,ccrho,ccTP,ccS,ccTPTP)
         RHO(I,J,K)=local_density2(TMPD,SLTD,ccrho,ccTP,ccS,ccTPTP)*1.D-3
       END DO
     END DO
   END DO
   DO K=Kb1,Kb2
     N=0
     DO J=Jb1,Jb2
       DO I=Ib1,Ib2
         IF (RHO(I,J,K)<RHO(I,J,K-1)) THEN
           N=N+1
           T(I,J,K)=T(I,J,K-1)
           S(I,J,K)=S(I,J,K-1)
           TMPD=DBLE(T(I,J,K))
           SLTD=DBLE(S(I,J,K))
           PISD=1.D-3*DBLE(Z(2*K))
           CALL make_local_coefficients2(TMPD,SLTD,PISD,ccrho,ccTP,ccS,ccTPTP)
           RHO(I,J,K)=local_density2(TMPD,SLTD,ccrho,ccTP,ccS,ccTPTP)*1.D-3
         END IF
       END DO
     END DO
     WRITE(*,58) N,K
   END DO
 ELSE IF (FL_EOS_OP.EQ.6) THEN
   DO K=Ka1,Ka2
     DO J=Ja1,Ja2
       DO I=Ia1,Ia2
         TMPD=DBLE(T(I,J,K))
         SLTD=DBLE(S(I,J,K))
         PISD=1.D-3*DBLE(Z(2*K))
         CALL make_local_coefficients3(TMPD,SLTD,PISD,ccrho,ccTP,ccS,ccTPTP,ccSTP)
         RHO(I,J,K)=local_density3(TMPD,SLTD,ccrho,ccTP,ccS,ccTPTP,ccSTP)*1.D-3
       END DO
     END DO
   END DO
   DO K=Kb1,Kb2
     N=0
     DO J=Jb1,Jb2
       DO I=Ib1,Ib2
         IF (RHO(I,J,K)<RHO(I,J,K-1)) THEN
           N=N+1
           T(I,J,K)=T(I,J,K-1)
           S(I,J,K)=S(I,J,K-1)
           TMPD=DBLE(T(I,J,K))
           SLTD=DBLE(S(I,J,K))
           PISD=1.D-3*DBLE(Z(2*K))
           CALL make_local_coefficients3(TMPD,SLTD,PISD,ccrho,ccTP,ccS,ccTPTP,ccSTP)
           RHO(I,J,K)=local_density3(TMPD,SLTD,ccrho,ccTP,ccS,ccTPTP,ccSTP)*1.D-3
         END IF
       END DO
     END DO
     WRITE(*,58) N,K
   END DO
 END IF

 ELSE
 write(*,*) "WHY DIDNT YOU PUT EOS_FILE=0 IN STABILITY1 ARGUMENTS???"
!! Open file for storing coefficients of EOS scheme
!  DO K=1,SIZE(ZZ)
!     Z( K )=ZZ ( K)
!  END DO
!   TEMPDIR=TRIM(WORKDIR)//'/INPUT_GLOBAL/eos'
!   CALL OPEN_BIGENDIAN_NEW_BIN_FILE(TRIM(TEMPDIR),FNO_eos)
!
!   IF (FL_EOS_OP.EQ.4) THEN
!     DO K=1,K1
!       PISD=1.D-3*DBLE(Z(2*K))
!       DO J=1,J0
!         DO I=1,I0
!           TMPD=DBLE(T(I,J,K))
!           SLTD=DBLE(S(I,J,K))
!           PISD=1.D-3*DBLE(Z(2*K))
!           CALL make_local_coefficients(TMPD,SLTD,PISD,ccrho,ccTP,ccS)
!           Mcrho(I,J,K)=ccrho
!           McTP(I,J,K)=ccTP
!           McS(I,J,K)=ccS
!         END DO
!       END DO
!     END DO
!     WRITE(FNO_eos) Mcrho,McTP,McS
!     CLOSE(FNO_eos)
!   ELSE IF (FL_EOS_OP.EQ.5) THEN
!     DO K=1,K1
!       PISD=1.D-3*DBLE(Z(2*K))
!       DO J=1,J0
!         DO I=1,I0
!           TMPD=DBLE(T(I,J,K))
!           SLTD=DBLE(S(I,J,K))
!!          PISD=1.D-3*DBLE(Z(2*K))
!           CALL make_local_coefficients2(TMPD,SLTD,PISD,ccrho,ccTP,ccS,ccTPTP)
!           Mcrho(I,J,K)=ccrho
!           McTP(I,J,K)=ccTP
!           McS(I,J,K)=ccS
!           McTPTP(I,J,K)=ccTPTP
!         END DO
!       END DO
!     END DO
!     WRITE(FNO_eos) Mcrho,McTP,McS,McTPTP
!     CLOSE(FNO_eos)
!   ELSE IF (FL_EOS_OP.EQ.6) THEN
!     DO K=1,K1
!       PISD=1.D-3*DBLE(Z(2*K))
!       DO J=1,J0
!         DO I=1,I0
!           TMPD=DBLE(T(I,J,K))
!           SLTD=DBLE(S(I,J,K))
!!          PISD=1.D-3*DBLE(Z(2*K))
!           CALL make_local_coefficients3(TMPD,SLTD,PISD,ccrho,ccTP,ccS,ccTPTP,ccSTP)
!           Mcrho(I,J,K)=ccrho
!           McTP(I,J,K)=ccTP
!           McS(I,J,K)=ccS
!           McTPTP(I,J,K)=ccTPTP
!           McSTP(I,J,K)=ccSTP
!         END DO
!       END DO
!     END DO
!     WRITE(FNO_eos) Mcrho,McTP,McS,McTPTP,McSTP
!     CLOSE(FNO_eos)
!   END IF
 END IF

 58    FORMAT(I6,' unstable points at level',I3)

end subroutine EOS_COMPUTATION_STABILITY1
	


SUBROUTINE make_local_coefficients3(TP,S,P,crho,cTP,cS,cTPTP,cSTP)
 ! Calculates local coefficients
 ! TP    potential temperature [degrees, Celcius]
 ! S     salinity              [psu, effectively ppt]
 ! P     pressure              [bars, 10^5 Pascals]
 ! crho  local density when TP=S=0
 ! cTP   gradient of local density wrt TP
 ! cS    gradient of local density wrt S
 ! cTPTP .5*second derivative of local density wrt TP
 ! cSTP  cross-derivative of local density wrt S and TP [kg/m^3/C/ppt]
! USE UNESCO_density, ONLY: density_T_S_TT_TS, temperature, &
!                           potential_temperature_T, potential_temperature_TT, &
!                           potential_temperature_TS
 IMPLICIT NONE
 REAL(8),INTENT(IN) :: TP,S,P 
 REAL(8),INTENT(OUT) :: crho,cTP,cS,cTPTP,cSTP 
 REAL(8) :: T, dum
 ! find coeffs for rho = crho + cTP*(TP-TPr) + cS*(S-Sr) 
 !                            + cTPTP*(TP-TPr)**2 + cSTP*(S-Sr)*(TP-TPr)
 T=temperature(TP,S,P)
 CALL density_T_S_TT_TS(T,S,P,crho,cTP,cS,cTPTP,cSTP)
 dum=1.D0/potential_temperature_T(T,S,P)
 cTP=cTP*dum
 cSTP = dum*(cSTP - cTP*potential_temperature_TS(P))
 cTPTP=.5*dum*dum*(cTPTP-cTP*potential_temperature_TT(T,P))
 ! modify coeffs for rho = crho + cTP*TP + cS*S +cTPTP*TP**2 + cSTP*S*TP
 crho=crho-(cTP-cTPTP*TP)*TP-cS*S+cSTP*S*TP
 cTP=cTP-2*cTPTP*TP-cSTP*S
 cS=cS-cSTP*TP
END SUBROUTINE make_local_coefficients3               


SUBROUTINE make_local_coefficients2(TP,S,P,crho,cTP,cS,cTPTP)
 ! Calculates local coefficients
 ! TP    potential temperature [degrees, Celcius]
 ! S     salinity              [psu, effectively ppt]
 ! P     pressure              [bars, 10^5 Pascals]
 ! crho  local density when TP=S=0
 ! cTP   gradient of local density wrt TP
 ! cS    gradient of local density wrt S
 ! cTPTP .5*second derivative of local density wrt TP
! USE UNESCO_density, ONLY: density_T_S_TT_TS, temperature, &
!                           potential_temperature_T, potential_temperature_TT
 IMPLICIT NONE
 REAL(8),INTENT(IN) :: TP,S,P
 REAL(8),INTENT(OUT) :: crho,cTP,cS,cTPTP
 REAL(8) :: T, dum
 ! find coeffs for rho = crho + cTP*(TP-TPr) + cS*(S-Sr) + cTPTP*(TP-TPr)**2 
 T=temperature(TP,S,P)
 CALL density_T_S_TT_TS(T,S,P,crho,cTP,cS,cTPTP)
 dum=1.D0/potential_temperature_T(T,S,P)
 cTP=cTP*dum
!YC ??
 cTPTP=.5*dum*dum*(cTPTP-cTP*potential_temperature_TT(T,P))
 ! modify coeffs for rho = crho + cTP*TP + cS*S +cTPTP*TP**2 
 crho=crho-(cTP-cTPTP*TP)*TP-cS*S
 cTP=cTP-2*cTPTP*TP
END SUBROUTINE make_local_coefficients2               

SUBROUTINE make_local_coefficients(TP,S,P,crho,cTP,cS)
 ! Calculates local coefficients 
 ! TP    potential temperature [degrees, Celcius]
 ! S     salinity              [psu, effectively ppt]
 ! P     pressure              [bars, 10^5 Pascals]
 ! crho  local density when TP=S=0
 ! cTP   gradient of local density wrt TP
 ! cS    gradient of local density wrt S
! USE UNESCO_density, ONLY: density_T_S_TT_TS,temperature,potential_temperature_T
 IMPLICIT NONE
 REAL(8),INTENT(IN)::TP,S,P
 REAL(8),INTENT(OUT)::crho,cTP,cS
 REAL(8)::T
 ! find coeffs for rho = crho + cTP*(TP-TPr) + cS*(S-Sr)  
 T=temperature(TP,S,P)
!CALL density_T_S(T,S,P,crho,cTP,cS)
 CALL density_T_S_TT_TS(T,S,P,crho,cTP,cS)
 cTP=cTP/potential_temperature_T(T,S,P)
 ! modify coeffs for rho = crho + cTP*TP + cS*S 
 crho=crho-cTP*TP-cS*S
END SUBROUTINE make_local_coefficients

SUBROUTINE density_T_S_TT_TS(T,S,P,rho,drho_dT,drho_dS,drho_dT_dT, &
                                                       drho_dT_dS)
 ! Uses equation (A3.1) of Gill 1982 Atmosphere Ocean Dynamics to
 ! calculate the density of fresh water and then uses (A3.2) to
 ! calculate the density of saline water at one standard atmosphere.
 ! Pressure effects then accounted for using (A3.3).
 !
 ! The derivatives of density with respect to temperature T and
 ! salinity S are also calculated. The second derivative wrt T
 ! and the TS cross-derivative are calculated.
 !
 ! INPUT VARIABLES:
 ! T   temperature, [degrees Celcius]
 ! S   salinity, [psu, which effectively equals ppt]
 ! P   pressure, [bars, 10^5 Pascals]
 !
 ! OUTPUT VARIABLES:
 ! rho         density  at T,S,P            [kg/m^3]
 ! drho_dT     derivative of density wrt T  [kg/m^3/C]
 ! drho_dS     derivative of density wrt S  [kg/m^3/ppt]
 !
 ! OPTIONAL OUTPUT VARIABLES:
 ! drho_dT_dT  2nd derivative of density wrt T  [kg/m^3/C/C]
 ! drho_dT_dS  cross derivative of density wrt T and S  [kg/m^3/C/ppt]
 !
 ! test values are rho(T=5, S=0, P=0)      =  999.96675 - 1000
 !                 rho(T=5, S=35,P=0)      = 1027.67547 - 1000
 !                 rho(T=25, S=35, P=1000) = 1062.53817 - 1000
 !
 ! This algorithm is modified so it works well at single precision.
 ! This is done by extracting an offset density of 1000 kg/m^3
! USE global_parameters, ONLY: ks
 IMPLICIT NONE
 REAL(8),INTENT(IN) :: T,S,P
 REAL(8),INTENT(OUT) :: rho,drho_dT,drho_dS
 REAL(8),INTENT(OUT),OPTIONAL :: drho_dT_dT,drho_dT_dS
 REAL(8) :: K,dum, dKdT, dKdS, r,dKdTdT, rdT, rdS,dKdSdT
 ! calculate the density of fresh water at temperature T
 rho = -0.1574060 + (6.793952E-2 - (9.095290E-3 - (1.001685E-4 &
            - (1.120083E-6 - 6.536332E-9*T)*T)*T)*T)*T;
 
 ! calculate the derivative wrt T of density of fresh water at temperature T
 drho_dT =  6.793952E-2 - (1.819058E-2 - (3.005055E-4  &
            -(4.480332E-6 -3.268166E-8*T)*T)*T)*T;

 ! calculate the 2nd derivative wrt T of density of fresh water at temperature T
 IF (PRESENT(drho_dT_dT)) drho_dT_dT = - (1.819058E-2 - (6.01011E-4         &
                                       - (1.3440996E-5 - 1.3072664E-7*T)*T)*T);
 
 dum=SQRT(S)
 
 ! calculate the density of saline water at 1 atmosphere
 rho = rho + S*(0.824493 +4.8314E-4*S - (4.0899E-3 - (7.6438E-5  &
               - (8.2467E-7 - 5.3875E-9*T)*T)*T)*T)              &
         + S*dum*(-5.72466E-3 + (1.0227E-4 - 1.6546E-6*T)*T)
 
 ! calculate the derivative wrt T of density of saline water at 1 atmosphere
 drho_dT = drho_dT                                                    &
           + S*(-4.0899E-3 + (1.52876E-4 - (2.47401E-6 - 2.155E-8*T)*T)*T)  &
           + S*dum*(1.0227E-4 - 3.3092E-6*T);
 
 ! calculate the 2nd derivative wrt T of density of saline water at 1 atmosphere
 IF (PRESENT(drho_dT_dT)) drho_dT_dT = drho_dT_dT                           &
                            + S*(1.52876E-4 - (4.94802E-6 - 6.465E-8*T)*T)  &
                            - S*dum*3.3092E-6;
 
 ! calculate the derivative wrt S of density of saline water at 1 atmosphere
 drho_dS = 0.824493 - (4.0899E-3 - (7.6438E-5                &
               - (8.2467E-7 - 5.3875E-9*T)*T)*T)*T           &
         + dum*(-8.58699E-3 + (1.53405E-4 - 2.4819E-6*T)*T)  &
         + 9.6628E-4*S;

 ! calculate the density cross-derivative wrt S and T of saline water at 1 atmos
 IF (PRESENT(drho_dT_dS))  drho_dT_dS = - (4.0899E-3 - (1.52876E-4        &
                                        - (2.47401E-6 - 2.155E-8*T)*T)*T) &
                                        + dum*(1.53405E-4 - 4.9638E-6*T)  

 
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 !%%%%%%%%  PRESSURE EFFECTS via secant bulk modulus  %%%%%%%%%%%%%%%%%%%
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 ! calculate the secant bulk modulus of pure water
 K = 19652.21 + T*(148.4206 - T*(2.327105 - T*(1.360477E-2  &
       - T*5.155288E-5)));
 
 ! calculate the derivative wrt T of the secant bulk modulus of pure water
 dKdT = 148.4206 - T*(4.65421 - T*(4.081431E-2 - T*2.0621152E-4));
 
 ! calculate the 2nd derivative wrt T of the secant bulk modulus of pure water
 dKdTdT =  - (4.65421 - T*(8.162862E-2 - T*6.1863456E-4));
 
 ! calculate the bulk modulus at one atmosphere
 K = K +S*(54.6746 - T*(.603459 - T*(1.09987E-2 - T*6.1670E-5))) &
       +S*dum*(7.944E-2 + T*(1.6483E-2 - T*5.3009E-4));
 
 ! calculate derivative wrt T of the bulk modulus at one atmosphere
 dKdT = dKdT +S*(- .603459 + T*(2.19974E-2 - T*1.8501E-4)) &
       +S*dum*(1.6483E-2 - 1.06018E-3*T);

 ! calculate 2nd derivative wrt T of the bulk modulus at one atmosphere
 IF (PRESENT(drho_dT_dT))  &
     dKdTdT = dKdTdT +S*(2.19974E-2 - T*3.7002E-4 - dum*1.06018E-3);
 
 ! calculate the derivative wrt S of the bulk modulus at one atmosphere
 dKdS = 54.6746 - T*(.603459 - T*(1.09987E-2 - T*6.1670E-5)) &
       +dum*(1.1916E-1 + T*(2.47245E-2 - T*7.95135E-4));

 ! calculate the cross derivative wrt S and T 
 ! of the bulk modulus at one atmosphere
 IF (PRESENT(drho_dT_dS)) dKdSdT = - (.603459 - T*(2.19974E-2 - T*1.8501E-4)) &
                                   + dum*(2.47245E-2 - T*1.59027E-3);
 
 ! correct the bulk modulus for pressure p
 K = K + P*(  (3.239908 + T*(1.43713E-3 + T*(1.16092E-4             &
                - T*5.77905E-7)))                                   &
            + S*(2.2838E-3 - T*(1.0981E-5 + T*1.6078E-6))           &
            + 1.91075E-4*S*dum                                      &
            + P*(   (8.50935E-5 - T*(6.12293E-6 - T*5.2787E-8))     &
                    +S*(-9.9348E-7 + T*(2.0816E-8 + T*9.1697E-10)) ));
 
 ! correct the derivative wrt T of the bulk modulus for pressure p
 dKdT = dKdT + P*(  (1.43713E-3 + T*(2.32184E-4 - T*1.733715E-6))   &
                  + S*(- 1.0981E-5 - 3.2156E-6*T)                   &
                  + P*(  (- 6.12293E-6 + 1.05574E-7*T)              &
                        + S*(2.0816E-8 + 1.83394E-9*T)   ));

 ! correct the 2nd derivative wrt T of the bulk modulus for pressure p
 IF (PRESENT(drho_dT_dT)) dKdTdT = dKdTdT + P*(  (2.32184E-4 - T*3.46743E-6)  &
                                                 - S*3.2156E-6                &
                                                 + P*(   1.05574E-7           &
                                                       + S*1.83394E-9   ));
 
 ! correct the derivative wrt S of the bulk modulus for pressure p
 dKdS = dKdS + P*(  (2.2838E-3 - T*(1.0981E-5 + T*1.6078E-6))       &
                   + 2.866125E-4*dum                                &
                   + P*(-9.9348E-7 + T*(2.0816E-8 + T*9.1697E-10))  );
 
 ! correct the cross derivative wrt S and T of the bulk modulus for pressure p
 IF (PRESENT(drho_dT_dS)) dKdSdT = dKdSdT - P*(  1.0981E-5 + T*3.2156E-6  &
                                                - P*(2.0816E-8 + T*1.83394E-9));
 
 ! correct the density for pressure effects
 r=rho;
 rdT=drho_dT
 rdS=drho_dS

 dum=1.D0/(K-P)
 rho = (K*rho+1000*P)*dum;
 drho_dT = (K*drho_dT + (r - rho)*dKdT)*dum ;
 drho_dS = (K*drho_dS + (r - rho)*dKdS)*dum ;

 IF (PRESENT(drho_dT_dT)) drho_dT_dT = dum*(r*dKdTdT+2*rdT*dKdT+K*drho_dT_dT  &
                                            + rho*(2*dum*dKdT*dKdT-dKdTdT)    &
                                            - 2*dum*(r*dKdT+K*rdT)*dKdT);               

 IF (PRESENT(drho_dT_dS)) drho_dT_dS = dum*( dKdS*rdT+K*drho_dT_dS+rdS*dKdT  &
                                           +r*dKdSdT-rho*dKdSdT-drho_dS*dKdT &
                                           -drho_dT*dKdS)
 
END SUBROUTINE density_T_S_TT_TS

FUNCTION local_linear_density(TP,S,crho,cTP,cS)
 ! Calculates the density from local coefficients.
 ! TP    potential temperature [degrees, Celcius]
 ! S     salinity              [psu, effectively ppt]
 ! crho  local density when TP=S=0
 ! cTP   gradient of local density wrt TP [kg/m^3/C]
 ! cS    gradient of local density wrt S  [kg/m^3/ppt]
 IMPLICIT NONE
! REAL(KIND=ks), DIMENSION(:,:,:) :: TP,S,crho,cTP,cS
! REAL(KIND=ks), DIMENSION(SIZE(S,1),SIZE(S,2),SIZE(S,3)) :: local_linear_density
 REAL(8) :: TP,S,crho,cTP,cS,local_linear_density
 local_linear_density=crho+cTP*TP+cS*S
END FUNCTION local_linear_density

FUNCTION local_density2(TP,S,crho,cTP,cS,cTPTP)
 ! Calculates the density from local coefficients.
 ! TP    potential temperature [degrees, Celcius]
 ! S     salinity              [psu, effectively ppt]
 ! crho  local density when TP=S=0
 ! cTP   gradient of local density wrt TP [kg/m^3/C]
 ! cS    gradient of local density wrt S  [kg/m^3/ppt]
 ! cTPTP .5*second derivative of local density wrt TP [kg/m^3/C]
 ! Operation count = 6
 IMPLICIT NONE
 REAL(8) :: TP,S,crho,cTP,cS,cTPTP,local_density2
 local_density2=crho+(cTP+cTPTP*TP)*TP+cS*S
END FUNCTION local_density2                             

FUNCTION local_density3(TP,S,crho,cTP,cS,cTPTP,cSTP)
 ! Calculates the density from local coefficients.
 ! TP    potential temperature [degrees, Celcius]
 ! S     salinity              [psu, effectively ppt]
 ! crho  local density when TP=S=0
 ! cTP   gradient of local density wrt TP [kg/m^3/C]
 ! cS    gradient of local density wrt S  [kg/m^3/ppt]
 ! cTPTP .5*second derivative of local density wrt TP [kg/m^3/C/C]
 ! cSTP  cross-derivative of local density wrt S and TP [kg/m^3/C/ppt]
 ! Operation count = 8
 IMPLICIT NONE
 REAL(8) :: TP,S,crho,cTP,cS,cTPTP,cSTP,local_density3
 local_density3=crho+(cTP+cTPTP*TP+cSTP*S)*TP+cS*S
END FUNCTION local_density3


FUNCTION potential_temperature_T(T,S,P)
 ! Calculates the derivative of potential temperature wrt insitu temperature T
 ! T  temperature in degrees Celcius
 ! S  salinity in parts per thousand
 ! P  pressure in bars (10^5 Pascals)
 !
 ! Note: the derivative of T wrt potential temperature is 
 !       1/potential_temperature_T
 !
 ! References:
 ! Gill, A. E. (1981) Atmosphere-Ocean Dynamics. Academic Press
 ! Bryden, H. L. (1973) New polynomials for thermal expansion,
 !        adiabatic temperature gradient and geopotential
 !        temperature gradient of sea water. Deep-Sea Res., 20, 401-408.
 IMPLICIT NONE
 REAL(8)::T,S,P,potential_temperature_T
 potential_temperature_T = 1. +                                  &
             P*( -8.3198e-5 + T*(1.0813e-6 - T*1.20822e-8)       &
                +(S-35.)*2.9778e-7                               &
                +P*( (3.1628e-8 - T*4.3974e-10) - P*5.0484e-12 ) );
END FUNCTION potential_temperature_T

FUNCTION potential_temperature_TT(T,P)
 ! Calculates the 2nd derivative of potential temperature wrt temperature T
 ! T  temperature in degrees Celcius
 ! P  pressure in bars (10^5 Pascals)
 !
  ! References:
 ! Gill, A. E. (1981) Atmosphere-Ocean Dynamics. Academic Press
 ! Bryden, H. L. (1973) New polynomials for thermal expansion,
 !        adiabatic temperature gradient and geopotential
 !        temperature gradient of sea water. Deep-Sea Res., 20, 401-408.
 IMPLICIT NONE
 REAL(8) :: T,P,potential_temperature_TT
 potential_temperature_TT = P*( 1.0813e-6 - T*2.41644e-8 - P*4.3974e-10 );
END FUNCTION potential_temperature_TT

FUNCTION potential_temperature_TS(P)
 ! Calculates the second derivative of potential temperature wrt 
 ! insitu temperature T and salinity S
 ! P  pressure in bars (10^5 Pascals)
 !
 ! References:
 ! Gill, A. E. (1981) Atmosphere-Ocean Dynamics. Academic Press
 ! Bryden, H. L. (1973) New polynomials for thermal expansion,
 !        adiabatic temperature gradient and geopotential
 !        temperature gradient of sea water. Deep-Sea Res., 20, 401-408.
 IMPLICIT NONE
 REAL(8) :: P,potential_temperature_TS
  potential_temperature_TS = P*2.9778e-7 
END FUNCTION potential_temperature_TS


FUNCTION Wright_density(TP,S,P)
! Dan Wright's full equation of state
! TP potential temperature [degrees Celcius]
! S  salinity [ppt, parts per thousand]
! P  pressure .001/cm depth change [bars, 10^5 Pascals]
! Wright_density [grams/cm^3] multiply by 1000 to get kg/m^3
 REAL(8)::TP,S,P,p01,p02,Wright_density 

 p02=1747.4508988+TP*(11.51588-0.046331033*TP)-S*(3.85429655+0.01353985*TP)
 p01=P+5884.81703666+TP*(39.803732+TP*(-0.3191477 &
                    +TP*0.0004291133))+2.6126277*S
 Wright_density=p01/(p02+0.7028423*p01)
END FUNCTION Wright_density

FUNCTION densityTSP(T,S,p)
  ! Uses equation (A3.1) of Gill 1982 Atmosphere Ocean Dynamics to 
  ! calculate the density of fresh water and then uses (A3.2) to
  ! calculate the density of saline water at one standard atmosphere.
  ! Pressure effects then accounted for using (A3.3).
  ! Test value is   densityTSP(T=25, S=35, p=1000) = 1062.53817-1000
  IMPLICIT NONE
  REAL(8)::T,S,p,K,S3o2,densityTSP

  densityTSP = densityTS(T,S) ! density of saline water at 1 atmosphere
  S3o2=dsqrt(S**3)
  ! calculate the secant bulk modulus of pure water
  K = 19652.21D0 + T*(148.4206D0 - T*(2.327105D0 - T*(1.360477D-2 &
      - T*5.155288D-5)));
  ! calculate the bulk modulus at one atmosphere
  K = K +S*(54.6746D0 - T*(.603459D0 - T*(1.09987D-2 - T*6.1670D-5))) &
        +S3o2*(7.944D-2 + T*(1.6483D-2 - T*5.3009D-4));
  ! correct the bulk modulus for pressure p
  K = K + p*(                                                         &
         3.239908D0 + T*(1.43713D-3 + T*(1.16092D-4 - T*5.77905D-7)) &
        + S*(2.2838D-3 - T*(1.0981D-5 + T*1.6078D-6))                 &
        + 1.91075D-4*S3o2                                             &
        + p*(8.50935D-5 - T*(6.12293D-6 - T*5.2787D-8))               &
        + p*S*(-9.9348D-7 + T*(2.0816D-8 + T*9.1697D-10))             &
              );
  ! correct the density for pressure effects
! densityTSP = densityTSP/(1.D0-p/K) ! original UNESCO, density is not offset
  ! Brian Sanderson recommends the following modification to the UNESCO
  ! equation of state. This modification enables much better accuracy
  ! to be obtained using single precision arithmetic.
! K=p/K
! densityTSP = 1000.D0*K/(1.D0-K) + densityTSP/(1.D0-K) 
  densityTSP = (1000.D0*P+densityTSP*K)/(K-P)
 END FUNCTION densityTSP

FUNCTION densityTS(T,S)
  ! Uses equation (A3.1) of Gill 1982 Atmosphere Ocean Dynamics to 
  ! calculate the density of fresh water and then uses (A3.2) to
  ! calculate the density of saline water at one standard atmosphere.
  ! Test value is    densityTS(T=5, S=35) = 1027.67547-1000
  IMPLICIT NONE
  REAL(8)::T,S,densityTS

  densityTS=densityT(T)    ! density of fresh water at temperature T
  densityTS = densityTS + S*(0.824493D0 - T*(4.0899D-3 - T*(7.6438D-5 &
             - T*(8.2467D-7 - T*5.3875D-9))))                       &
       + S*DSQRT(S)*(-5.72466D-3 + T*(1.0227D-4 - T*1.6546D-6))      &
       + 4.8314D-4*S*S;    ! density of saline water at 1 atmosphere
 END FUNCTION densityTS

FUNCTION densityT(T)
  ! Uses equation (A3.1) of Gill 1982 Atmosphere Ocean Dynamics to 
  ! calculate the density of fresh water at temperature T.
  ! Test value is    densityT(T=5) =  999.96675-1000
  ! ***NOTE*** for accuracy purposes we subtract 1000 from the density
  IMPLICIT NONE
  REAL(8)::T,densityT

! densityT=999.842594 + T*(6.793952E-2 - T*(9.095290E-3     &
! densityT=-20.1574060 + T*(6.793952E-2 - T*(9.095290E-3    &
  densityT=-0.1574060D0 + T*(6.793952D-2 - T*(9.095290D-3    &
          - T*(1.001685D-4 -  T*(1.120083D-6 - T*6.536332D-9))))
 END FUNCTION densityT

FUNCTION temperature(TP,S,P)
  ! Calculates insitu temperature from potential temperature
  ! TP potential temperature [degrees Celcius]
  ! S  salinity [ppt]
  ! P  pressure [bars, 10^5 pascals]
  !
  ! The UNESCO equation for potential_temperature is inverted iteratively.
  !
  ! See also: potential_temperature, density
  IMPLICIT NONE
  REAL(8)::TP,S,P,T,temperature
  INTEGER :: iterate
  T=TP;             
  temperature=T;     ! first `guess'
  ! under relax to obtain temperature from potential temperatuer
  DO iterate=1,100   
   T=T+.999D0*(TP-potential_temperature(T,S,P));
!  T=T+(TP-potential_temperature(T,S,P));
   IF (DABS(temperature-T)<1.D-5) THEN
!     PRINT *,'number of iterations = ',iterate
      EXIT
   END IF
   temperature=T;
  END DO
  IF (iterate>=100) THEN
    PRINT *,'ERROR: temperature did not converge in 100 iterations'
    PRINT *,'Program stopped by FUNCTION temperature in UNESCO_density.f90'
    STOP
  END IF
 END FUNCTION temperature


FUNCTION potential_temperature(T,S,P)
  ! Calculates the potential temperature from
  ! T  temperature  [degrees Celcius]
  ! S  salinity in [parts per thousand (ppt)]
  ! P  pressure in [bars (10^5 Pascals)]
  ! EXAMPLE:
  !        t=potential_temperature(T=10,S=25,P=1000)
  !        gives t=8.46785160000000
  ! References:
  ! Gill, A. E. (1981) Atmosphere-Ocean Dynamics. Academic Press
  ! Bryden, H. L. (1973) New polynomials for thermal expansion,
  !        adiabatic temperature gradient and geopotential
  !        temperature gradient of sea water. Deep-Sea Res., 20, 401-408.
  IMPLICIT NONE
  REAL(8)::T,S,p,potential_temperature
  potential_temperature =                                             &
       T-P*(3.6504D-4 + T*(8.3198D-5 - T*(5.4065D-7 - T*4.0274D-9)))  &
        -P*(S-35.D0)*(1.7439D-5 - T*2.9778D-7)                          & 
        -P*P*(8.9309D-7 - T*(3.1628D-8 - T*2.1987D-10))               &
        +4.1057D-9*(S-35.D0)*P*P                                        &
        -P*P*P*(-1.6056D-10 + T*5.0484D-12);
 END FUNCTION potential_temperature


END MODULE density_EOS_GLOBAL
