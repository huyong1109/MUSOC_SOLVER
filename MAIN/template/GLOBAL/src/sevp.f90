! ----------------------------------------------------------------------
! SEVP elliptic solver arrays
! SEVP subregion boundary array "IE" is a user-defined array
! Double precision SEVP elliptic solver arrays
!
! ******************
! *ELLIPTIC SOLVER *
! ******************
! ----------------------------------------------------------------------
SUBROUTINE REPBIR(AX,AY,BB,CX,CY,RINV,RINV1,DUM0,DUM1,DUM2,F,H,X,IE,I0,I2,M0,M2,NBLK)
! ----------------------------------------------------------------------
  INTEGER, INTENT (IN) :: I0,I2,M0,M2,NBLK
  REAL(8) :: RINV(M2,M2,*),RINV1(M2,M2,*),DUM0(M2,*),DUM1(*),DUM2(*),H(M0,*),X(I0,*)
  REAL :: AX(I2,*),AY(I2,*),BB(I2,*),CX(I2,*),CY(I2,*),F(M2,*)
  INTEGER :: IE(*)
  
  ! diagonalize R
  !write(*,*) "Dimsize F",size(F)
  !write(*,*) "Dimsize OBB",size(OBB)
  !write(*,*) "IE(NBLK) = ", IE(NBLK), " J2 = ",J2
  !write(*,*) AX(4,4),AY(4,4),BB(4,4),OBB(4,4),CX(4,4),CY(4,4)
  !F(1:M2,1:IE(NBLK)-2) = F(1:M2,1:IE(NBLK)-2) *OBB(1:M2,1:IE(NBLK)-2)

  
  JS=1
  DO NB=1,NBLK
    JF=IE(NB)-2
    DO J=JS,JF
    DO I=1,M2
      X(I+1,J+2)=(F(I,J)-AX(I,J)*X(I,J+1)-AY(I,J)*X(I+1,J)-BB(I,J)*X(I+1,J+1)-CX(I,J)*X(I+2,J+1))/CY(I,J)
    ENDDO
    ENDDO
    IF (NB.EQ.NBLK) GO TO 150
    J=IE(NB)-1
    DO I=1,M2
      DUM1(I)=F(I,J)-AX(I,J)*X(I,J+1)-AY(I,J)*X(I+1,J)-BB(I,J)*X(I+1,J+1)-CX(I,J)*X(I+2,J+1)-CY(I,J)*X(I+1,J+2)
    ENDDO
    J=IE(NB)
    DO N=1,M2
      DUM2(N)=0.
      DO M=1,M2
        DUM2(N)=DUM2(N)+DUM1(M)*RINV1(M,N,NB)
      ENDDO
      DUM0(N,NB)=X(N+1,J)
      X(N+1,J)=X(N+1,J)-DUM2(N)
    ENDDO
150 JS=IE(NB)
  ENDDO
  DO NBS=1,NBLK
    NB=NBLK-NBS+1
    JS=1
    IF (NB.NE.1) JS=IE(NB-1)
    JF=IE(NB)-2
    IF (NB.EQ.NBLK) GO TO 201
    J=IE(NB)
    DO N=1,M2
      X(N+1,J)=DUM0(N,NB)
    ENDDO
201 N=IE(NB)
    DO J=JS,N
    DO I=1,M0
      H(I,J)=0.
    ENDDO
    ENDDO
    J=IE(NB)-1
    DO I=1,M2
      DUM1(I)=F(I,J)-AX(I,J)*X(I,J+1)-AY(I,J)*X(I+1,J)-BB(I,J)*X(I+1,J+1)-CX(I,J)*X(I+2,J+1)-CY(I,J)*X(I+1,J+2)
    ENDDO
    DO N=1,M2
      DUM2(N)=0.
      DO M=1,M2
        DUM2(N)=DUM2(N)+DUM1(M)*RINV(M,N,NB)
      ENDDO
      H(N+1,JS+1)=DUM2(N)
      X(N+1,JS+1)=X(N+1,JS+1)+DUM2(N)
    ENDDO
    IF (NB.EQ.1) GO TO 250
    DO M=1,M2
      DUM1(M)=H(M+1,JS+1)*CY(M,JS-1)
    ENDDO
    J=IE(NB-1)
    DO N=1,M2
      DUM2(N)=0.
      DO M=1,M2
        DUM2(N)=DUM2(N)+DUM1(M)*RINV1(M,N,NB-1)
      ENDDO
      H(N+1,J)=DUM2(N)
    ENDDO
250 DO J=JS,JF
    DO I=1,M2
      H(I+1,J+2)=(-AX(I,J)*H(I,J+1)-AY(I,J)*H(I+1,J)-BB(I,J)*H(I+1,J+1)-CX(I,J)*H(I+2,J+1))/CY(I,J)
      X(I+1,J+2)=X(I+1,J+2)+H(I+1,J+2)
    ENDDO
    ENDDO
  ENDDO
END SUBROUTINE REPBIR


! *******************
! * ELLIPTIC SOLVER *
! *******************
! ----------------------------------------------------------------------
SUBROUTINE REP(AX,AY,BB,CX,CY,RINV,RINV1,DUM0,DUM1,DUM2,F,H,X,IE, I0,I2,NBLK)
! ----------------------------------------------------------------------
! NOTES
! REP is the SEVP Poisson solver for rigid pressure against rigid lid
! (see "Elliptic Marching Methods and Domain Decomposition" book by
! Patrick J. Roache).

! The Poisson approach is equivalent to solving implicitly free surface
! gravity wave terms in the limit as time step goes to infinity, and is
! also equivalent to assuming that the divergence of the barotropic mode
! is zero (as does the rigid lid approximation for incompressible flow).

! The big influence coefficient arrays, RINV and RINV1, are calculated
! in the PREPXXX/prepxxx.f codes for the various grids (XXX). These
! depend only on the grid resolutions and bathymetry data calculated in
! PREPXXX/DATA/inmets.f. The SEVP solver efficiency and roundoff error
! in the solution depend on the choice of NB0 and the IE vector. The
! NB0 value and the spacing of the IE values should ideally be large
! enough to allow reasonable efficiency while small enough for tolerable
! roundoff error. It is recommended that the IE spacings be about 5-8 for
! 64-bit word length (REAL*8) of the RINV, RINV1 and other arrays
! associated with their derivation and usage. Longitudinal multi-grid
! (domain decomposition) used herein for all the field operations
! reduces the size of RINV and RINV1, while being nearly seamless and
! accurate. A truly 4th-order-accurate domain decomposition is being
! developed. To further reduce storage, the SEVP solver may itself be
! decomposed within any of the grids, but that is not implemented here.

  REAL*8 RINV,RINV1,DUM0,DUM1,DUM2,X,H
  DIMENSION AX(I2,*),AY(I2,*),BB(I2,*),CX(I2,*),CY(I2,*), RINV(I2,I2,*),RINV1(I2,I2,*),H(I0,*),IE(*),DUM0(I2,*),DUM1(*), DUM2(*),F(I2,*),X(I0,*)

  JS=1
  DO NB=1,NBLK
    JF=IE(NB)-2
    DO J=JS,JF
    DO I=1,I2
      X(I+1,J+2)=(F(I,J)-AX(I,J)*X(I,J+1)-AY(I,J)*X(I+1,J)-BB(I,J)* X(I+1,J+1)-CX(I,J)*X(I+2,J+1))/CY(I,J)
    ENDDO
    ENDDO
    IF (NB.EQ.NBLK) GO TO 150
    J=IE(NB)-1
    DO I=1,I2
      DUM1(I)=F(I,J)-AX(I,J)*X(I,J+1)-AY(I,J)*X(I+1,J)-BB(I,J)* X(I+1,J+1)-CX(I,J)*X(I+2,J+1)-CY(I,J)*X(I+1,J+2)
    ENDDO
    J=IE(NB)
    DO N=1,I2
      DUM2(N)=0.
      DO M=1,I2
        DUM2(N)=DUM2(N)+DUM1(M)*RINV1(M,N,NB)
      ENDDO
      DUM0(N,NB)=X(N+1,J)
      X(N+1,J)=X(N+1,J)-DUM2(N)
    ENDDO
150 JS=IE(NB)
  ENDDO
  
  DO NBS=1,NBLK
    NB=NBLK-NBS+1
    JS=1
    IF (NB.NE.1) JS=IE(NB-1)
    JF=IE(NB)-2
    IF (NB.EQ.NBLK) GO TO 201
    J=IE(NB)
    DO N=1,I2
      X(N+1,J)=DUM0(N,NB)
    ENDDO
201 N=IE(NB)
    DO J=JS,N
    DO I=1,I0
      H(I,J)=0.
    ENDDO
    ENDDO
    J=IE(NB)-1
    DO I=1,I2
      DUM1(I)=F(I,J)-AX(I,J)*X(I,J+1)-AY(I,J)*X(I+1,J)-BB(I,J)* X(I+1,J+1)-CX(I,J)*X(I+2,J+1)-CY(I,J)*X(I+1,J+2)
    ENDDO
    DO N=1,I2
      DUM2(N)=0.
      DO M=1,I2
        DUM2(N)=DUM2(N)+DUM1(M)*RINV(M,N,NB)
      ENDDO
      H(N+1,JS+1)=DUM2(N)
      X(N+1,JS+1)=X(N+1,JS+1)+DUM2(N)
    ENDDO
    IF (NB.EQ.1) GO TO 250
    DO M=1,I2
      DUM1(M)=H(M+1,JS+1)*CY(M,JS-1)
    ENDDO
    J=IE(NB-1)
    DO N=1,I2
      DUM2(N)=0.
      DO M=1,I2
        DUM2(N)=DUM2(N)+DUM1(M)*RINV1(M,N,NB-1)
      ENDDO
      H(N+1,J)=DUM2(N)
    ENDDO
250 DO J=JS,JF
      DO I=1,I2
        H(I+1,J+2)=(-AX(I,J)*H(I,J+1)-AY(I,J)*H(I+1,J)-BB(I,J)* H(I+1,J+1)-CX(I,J)*H(I+2,J+1))/CY(I,J)
        X(I+1,J+2)=X(I+1,J+2)+H(I+1,J+2)
      END DO
    ENDDO
  ENDDO
END

! *****************
! ELLIPTIC SOLVER *
! *****************
! ----------------------------------------------------------------------
      SUBROUTINE P_BICGSTAB_LE(AB,AL,AC,AR,AT,B,X,CB,CL,CC,CR,CT,R,RH,P,V,S,T,PH,SH)

      INCLUDE 'mpif.h'
      INCLUDE 'resglo.h'
      PARAMETER(I0J0=I0*J0)
      REAL*8 X,R,RH,P,V,S,T,PH,SH
C REAL*8 AB,AL,AC,AR,AT,X,R,RH,P,V,S,T,PH,SH,CB,CL,CC,CR,CT,B
      REAL*8 RHO,RHN,ALPHA,BETA,W,RES,TMP,TP,TOL
      DIMENSION AB(I2,J2),AL(I2,J2),AC(I2,J2),AR(I2,J2),AT(I2,0:J2),
     1 B(I2,J2),X(I0,J0),R(I0,J0),RH(I0,J0),P(I0,J0),V(I0,J0),
     2 S(I0,J0),T(I0,J0),PH(I0,J0),SH(I0,J0),
     3 CB(I2,J2),CL(I2,J2),CC(I2,J2),CR(I2,J2),CT(I2,J2)
      COMMON/MPI/M_CART,M_CLON,M_CLAT,MYID,MYLON,MYLAT,M_N,M_E,
     1 M_S,M_W,M_VLON,M_VLAT,M_V8LON,M_V8LAT,JS,JF,IERR,
     1 NPX(14),NPY(14),ISTAT(MPI_STATUS_SIZE)

      REAL*8 TIMER
      COMMON/TIMER/TIMER(10)

      ITER=0
 10   RHO=1.d0
      ALPHA=1.d0
      W=1.d0
      CALL P_SQPROD(AB,AL,AC,AR,AT,X,R)

      DO 50 J=2,J1
      DO 50 I=2,I1
 50   R(I,J)=B(I-1,J-1)-R(I,J)

      CALL DCOPY(I0J0,R,1,RH,1)
      CALL P_INPROD(R,R,TP)

      RES=SQRT(TP)
      IF (RES .LT. 1e-7) THEN
      CALL MPI_SENDRECV(X(2,1),1,M_V8LON,M_W,1,
     1 X(I0,1),1,M_V8LON,M_E,1,M_CART,ISTAT,IERR)
      CALL MPI_SENDRECV(X(I1,1),1,M_V8LON,M_E,1,
     1 X(1,1),1,M_V8LON,M_W,1,M_CART,ISTAT,IERR)
      CALL MPI_SENDRECV(X(1,2),1,M_V8LAT,M_S,1,
     1 X(1,J0),1,M_V8LAT,M_N,1,M_CART,ISTAT,IERR)
      CALL MPI_SENDRECV(X(1,J1),1,M_V8LAT,M_N,1,
     1 X(1,1),1,M_V8LAT,M_S,1,M_CART,ISTAT,IERR)
      RETURN
      ENDIF

 100  ITER=ITER+1
      CALL P_INPROD(R,RH,RHN)
      BETA=(RHN/RHO)*(ALPHA/W)

C PRINT*,ITER,'BETA=',BETA

      CALL DAXPY(I0J0,-1.d0*W,V,1,P,1)
      CALL DSCAL(I0J0,BETA,P,1)
      CALL DAXPY(I0J0,1.d0,R,1,P,1)

c DO 200 J=2,J1
C DO 200 I=2,I1
C 200 PH(I,J)=P(I,J)
C/AC(I-1,J-1)

C CALL P_SQELIM(CB,CL,CC,CR,CT,P,PH)
C CALL P_LUSQELIM(CB,CL,CC,CR,CT,P,PH)
      CALL P_SIPSQELIM(CB,CL,CC,CR,CT,P,PH)
      CALL P_SQPROD(AB,AL,AC,AR,AT,PH,V)
      CALL P_INPROD(RH,V,TMP)
      ALPHA=RHN/TMP

C PRINT*,ITER,'ALPHA=', ALPHA

      CALL DCOPY(I0J0,R,1,S,1)
      CALL DAXPY(I0J0,-1.d0*ALPHA,V,1,S,1)
C CALL DCOPY(I0J0,S,1,SH,1)

C DO 300 J=2,J1
C DO 300 I=2,I1
C 300 SH(I,J)=S(I,J)
C/AC(I-1,J-1)

      CALL P_SIPSQELIM(CB,CL,CC,CR,CT,S,SH)
C CALL P_LUSQELIM(CB,CL,CC,CR,CT,S,SH)
C CALL P_SQELIM(CB,CL,CC,CR,CT,S,SH)
      CALL P_SQPROD(AB,AL,AC,AR,AT,SH,T)
      CALL P_INPROD(T,T,TMP)
      CALL P_INPROD(T,S,W)
      W=W/TMP

C PRINT*,ITER,'WN=',W

      TMP=0.d0

      CALL DAXPY(I0J0,ALPHA,PH,1,X,1)
      CALL DAXPY(I0J0,W,SH,1,X,1)
      CALL DCOPY(I0J0,S,1,R,1)
      CALL DAXPY(I0J0,-1.d0*W,T,1,R,1)

      CALL P_INPROD(R,R,TP)

      RHO=RHN
C IF (MYID .EQ. 0) PRINT*,ITER,SQRT(TP)
C,ALPHA,RHN,W,SQRT(TP)
      IF (SQRT(TP)/RES .GT. 1e-3 .AND. MOD(ITER,100) .NE. 0) GOTO 100
C IF (MYID .EQ. 0) PRINT*, ITER,SQRT(TP)
      IF (SQRT(TP) .GT. 1e-7 .AND. ITER .LT. 10000) GOTO 10

      CALL MPI_BARRIER(M_CART,IERR)
      TIMER(8)=MPI_WTIME()
      CALL MPI_SENDRECV(X(2,1),1,M_V8LON,M_W,1,
     1 X(I0,1),1,M_V8LON,M_E,1,M_CART,ISTAT,IERR)
      CALL MPI_SENDRECV(X(I1,1),1,M_V8LON,M_E,1,
     1 X(1,1),1,M_V8LON,M_W,1,M_CART,ISTAT,IERR)
      CALL MPI_SENDRECV(X(1,2),1,M_V8LAT,M_S,1,
     1 X(1,J0),1,M_V8LAT,M_N,1,M_CART,ISTAT,IERR)
      CALL MPI_SENDRECV(X(1,J1),1,M_V8LAT,M_N,1,
     1 X(1,1),1,M_V8LAT,M_S,1,M_CART,ISTAT,IERR)

      CONTINUE
      END
