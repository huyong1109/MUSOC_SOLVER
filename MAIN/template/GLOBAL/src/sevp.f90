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
  !write(*,*) "IE = ",IE(1:NBLK)
  !write(*,*) AX(4,4),AY(4,4),BB(4,4),OBB(4,4),CX(4,4),CY(4,4)
  !F(1:M2,1:IE(NBLK)-2) = F(1:M2,1:IE(NBLK)-2) *OBB(1:M2,1:IE(NBLK)-2)
  JF = IE(NBLK)-2
  !write(*,*) 'ax'
  !write(*,*) ax(:,1:jf)
  !write(*,*) 'al'
  !write(*,*) ay(:,1:jf)
  !write(*,*) 'ac'
  !write(*,*) bb(:,1:jf)
  !write(*,*) 'ar'
  !write(*,*) cx(:,1:jf)
  !write(*,*) 'at'
  !write(*,*) cy(:,1:jf)

  !JF = IE(NBLK)-2
  !write(*,*) "check symmetry, JF=", JF
  !do j = 2,JF
  !  do i = 2,I2+1
  !    if (AX(i,j) /= CX(i-1,j))  write(*,*) "AX", i,j,AX(i,j),CX(i-1,j)
  !    if (CX(i,j) /= AX(i+1,j))  write(*,*) "CX", i,j,CX(i,j),AX(i+1,j)
  !    if (AY(i,j) /= CY(i,j-1))  write(*,*) "AY", i,j,AY(i,j),CY(i,j-1)
  !    if (CY(i,j) /= AY(i,j+1))  write(*,*) "CY", i,j,CY(i,j),AY(i,j+1)
  !  end do 
  !end do 
  !write(*,*) "END check symmetry"

  
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
subroutine p_bicgstab_le(al,ab,ac,ar,at,b,x,cl,cb,cc,cr,ct,i2,j2)
  ! input & output
  integer, intent(in) :: i2,j2 ! dimensions
  real,intent(in),dimension(i2,j2) ::  &
      ab,al,ac,ar,at,cl,cc,cr,ct
  real,intent(in),dimension(i2,j2) :: b
  real*8,intent(inout),dimension(i2+2,j2+2) :: x
  ! local variables 
  integer :: i0,j0 ! dimensions
  real*8 :: rho,rhn,alpha,beta,w,res,tmp,tp,tol
  real*8, dimension(i2+2,j2+2) ::r,rh,p,v,s,t,ph,sh 
  
  i0 = i2 + 2
  j0 = j2 + 2
  iter=0
  tp = 1.0d0
  !write(*,*) 'ab'
  !write(*,*) ab(:,:)
  !write(*,*) 'al'
  !write(*,*) al(:,:)
  !write(*,*) 'ac'
  !write(*,*) ac(:,:)
  !write(*,*) 'ar'
  !write(*,*) ar(:,:)
  !write(*,*) 'at'
  !write(*,*) at(:,:)
  
  call solver_boundary_update(x,i2,j2)
  call printVarSummaryR8(x,i0,j0,"init x")
  call printVarSummaryR4(b,i2,j2,"init b")
  r = 0.
  rh = 0.
  p = 0. 
  v = 0. 
  s = 0. 
  t = 0.
  ph = 0.
  sh = 0.

  do while(sqrt(tp) .gt. 1e-7 .and. iter .lt. 1000) 
     rho=1.d0
     alpha=1.d0
     w=1.d0
     print*,iter,'matrix-vector multiplication'
     call p_sqprod(al,ab,ac,ar,at,x,r,i2,j2)

     r(2:i2+1,2:j2+1)  = b(:,:)-r(2:i2+1,2:j2+1)
     call solver_boundary_update(r,i2,j2)

     rh(:,:) = r(:,:)
     print*,iter,'vector-vector'
     call p_inprod(r,r,tp,i2,j2) ! norm2 of r

     res=sqrt(tp)
     print*,iter,'res = ', res
     if (res .lt. 1e-7) then ! periodic boundary 
       return
     endif
     do while (sqrt(tp)/res .gt. 1e-3 )
        iter=iter+1
        call p_inprod(r,rh,rhn,i2,j2)
        beta=(rhn/rho)*(alpha/w)

        print*,iter,'beta=',beta
        ! P_i = R_t +beta(P_t - omega_t V_t)
        p(:,:) = r(:,:) + beta*(p(:,:) -w*v(:,:))

        call preconditioning(cl,cb,cc,cr,ct,p,ph,i2,j2)
        call solver_boundary_update(ph,i2,j2)
        call p_sqprod(al,ab,ac,ar,at,ph,v,i2,j2)
        call solver_boundary_update(v,i2,j2)
        call p_inprod(rh,v,tmp,i2,j2)
        alpha=rhn/tmp

        print*,iter,'alpha=', alpha
        s(:,:) = r(:,:) - alpha*v(:,:) 

        call preconditioning(cl,cb,cc,cr,ct,s,sh,i2,j2)
        call solver_boundary_update(sh,i2,j2)
        call p_sqprod(al,ab,ac,ar,at,sh,t,i2,j2)
        call solver_boundary_update(t,i2,j2)
        call p_inprod(t,t,tmp,i2,j2)
        call p_inprod(t,s,w,i2,j2)
        w=w/tmp


        tmp=0.d0
        x(:,:) = x(:,:) + alpha*ph(:,:) + w*sh(:,:)
        r(:,:) = s(:,:) -w*t(:,:)

        call p_inprod(r,r,tp,i2,j2)
        print*,iter,'tp', tp

        rho=rhn
        if( mod(iter,100) .eq. 0 ) then 
          exit
        end if 
      end do !if (sqrt(tp)/res .gt. 1e-3 .and. mod(iter,100) .ne. 0) goto 100
      
  end do !while (sqrt(tp) .gt. 1e-7 .and. iter .lt. 10000) goto 10

      call solver_boundary_update(x,i2,j2)
      continue
end subroutine 
subroutine p_csi(al,ab,ac,ar,at,mineig,maxeig,b,x,cl,cb,cc,cr,ct,i2,j2)
  ! input & output
  implicit none
  integer, intent(in) :: i2,j2 ! dimensions
  real*8, intent(in) :: mineig,maxeig
  real,intent(in),dimension(i2,j2) ::  &
      ab,al,ac,ar,at,cl,cb,cc,cr,ct
  real,intent(in),dimension(i2,j2) :: b
  real*8,intent(inout),dimension(i2+2,j2+2) :: x
  ! local variables 
  integer :: i0,j0,iter,i ! dimensions
  real*8 :: alpha,beta,csy,omega,tmp,tp
  real*8, dimension(i2+2,j2+2) ::r,s,q,work0,work1,a0r
  
  i0 = i2 + 2
  j0 = j2 + 2
  iter=0
  tp = 1.0d0
  r(:,:) = 0.
  s(:,:) = 0. 
  q(:,:) = 0. 
  work1 =0.
  call printVarSummaryR4(b,i0,j0,"pcsi b")
  call printVarSummaryR8(x,i0,j0,"pcsi x")
  !write(*,'(A15,2f10.5)')'Max-min eigs', mineig,maxeig 
  write(*,*)'Min-Max eigs', mineig,maxeig 
  alpha = 2.0/(maxeig-mineig)
  beta = (maxeig+mineig)/(maxeig-mineig)
  csy = beta/alpha
  omega = 2.0/csy

  call solver_boundary_update(x,i2,j2)
  call p_sqprod(al,ab,ac,ar,at,x,s,i2,j2)
  call solver_boundary_update(s,i2,j2)
  r(2:i2+1,2:j2+1)  = b(:,:)-s(2:i2+1,2:j2+1)
  call preconditioning(cl,cb,cc,cr,ct,r,work1,i2,j2)
  call solver_boundary_update(work1,i2,j2)
  q(:,:) = (1.0/csy)*work1(:,:)
  x(:,:) = x(:,:) +q(:,:)
  call p_sqprod(al,ab,ac,ar,at,x,s,i2,j2)
  call solver_boundary_update(s,i2,j2)
  r(2:i2+1,2:j2+1)  = b(:,:)-s(2:i2+1,2:j2+1)
  do while(sqrt(tp) .gt. 1e-7 .and. iter .lt. 1000) 
     omega = 1.0/(csy-omega/(4.0*alpha*alpha))
     call preconditioning(cl,cb,cc,cr,ct,r,work1,i2,j2)
     call solver_boundary_update(work1,i2,j2)
     q(:,:) = omega*work1(:,:)+ (csy*omega-1.0)*q(:,:)
     x(:,:) = x(:,:) + q(:,:)
     call p_sqprod(al,ab,ac,ar,at,x,s,i2,j2)
     r(2:i2+1,2:j2+1)  = b(:,:)-s(2:i2+1,2:j2+1)
     call solver_boundary_update(r,i2,j2)
     call p_inprod(r,r,tp,i2,j2)
     call printVarSummaryR8(r,i0,j0,"pcsi r")
     print*,iter,'tp', tp
     iter = iter +1

  end do !while (sqrt(tp) .gt. 1e-7 .and. iter .lt. 10000) goto 10
  call solver_boundary_update(x,i2,j2)
  continue
end subroutine 

subroutine lanczos(al,ab,ac,ar,at,mineig,maxeig,i2,j2)
  implicit none
  ! input & output
  integer, intent(in) :: i2,j2 ! dimensions
  real,intent(in),dimension(i2,j2) ::  &
      ab,al,ac,ar,at
  real*8,dimension(i2+2,j2+2) :: v,w
  ! local variables 
  integer :: i0,j0,ierr,iter,i ! dimensions
  real*8 :: mineig,maxeig,tp
  integer, parameter :: m = 10
  real*8,dimension(m) :: alpha,beta,sigma,mcsa,mcsb,mcsb2
  real*8,dimension(m,3) :: t
  real*8, dimension(i2+2,j2+2) ::vp,wp,s,q,work
  
  i0 = i2 + 2
  j0 = j2 + 2
  v = 1.
  vp = 0.
  wp = 0.
  alpha = 0.
  beta = 0.
  sigma = 0.
  call p_inprod(v,v,tp,i2,j2)
  write(*,*) 'lanczos tp', tp
  v = 1/sqrt(tp)*v
  w = v

  iter=0
  do while(iter .lt. m)
    iter = iter +1
    call p_sqprod(al,ab,ac,ar,at,v,s,i2,j2)
    call p_inprod(s,w,tp,i2,j2)
    alpha(iter) = tp
    work =beta(iter)*vp
    vp = v
    v = s- alpha(iter)*v - work
    call p_sqprod_T(al,ab,ac,ar,at,w,q,i2,j2)
    work = sigma(iter)*wp
    wp = w
    w = q- alpha(iter)*w - work
    call p_inprod(v,w,tp,i2,j2)
    sigma(iter+1) = sqrt(abs(tp))
    beta(iter+1) = sign(sigma(iter+1),tp)
    v = v/sigma(iter+1)
    w = w/beta(iter+1)
  end do 

  call p_sqprod(al,ab,ac,ar,at,v,s,i2,j2)
  call p_inprod(s,w,tp,i2,j2)
  alpha(m) = tp
  write(*,*) 'alpha '
  write(*,'(10f8.5)') alpha(:)
  write(*,*) 'beta  '
  write(*,'(10f8.5)') beta(:)
  write(*,*) 'sigma '
  write(*,'(10f8.5)') sigma(:)

  t(:,:) = 0.0
  t(2:m,1)   = sigma(2:m)
  t(1:m,2)   = alpha(1:m)
  t(1:m-1,3) = sigma(2:m)
  call figi ( m, t, mcsa, mcsb,mcsb2, ierr)
  write(*,*) 'mcsa'
  write(*,*) mcsa(:)
  write(*,*) 'mcsb'
  write(*,*) mcsb(:)

  maxeig = mcsa(1) + mcsb(1)
  do i = 2,m
    maxeig = max(maxeig, mcsa(i) + mcsb(i) + mcsb(i-1))
  end do 
  call ratqr(m, 1.0d-6, mcsa,mcsb,mineig,ierr)
  write(*,*) 'EIGS ierr: ', ierr
  write(*,*) 'Min-Max eigs: ', mineig,maxeig
      continue
end subroutine 

subroutine p_sqprod(al,ab,ac,ar,at,x,r,i2,j2)
  ! input & output
  integer, intent(in) :: i2,j2 ! dimensions
  real,intent(in),dimension(i2,j2) ::  &
      ab,al,ac,ar,at
  real*8,intent(in),dimension(i2+2,j2+2) :: x
  real*8,intent(inout),dimension(i2+2,j2+2) :: r
  ! local variables 
  integer :: i,j 

  do j = 1,j2
    do i = 1,i2
      r(i+1,j+1) =ab(i,j)*x(i+1,j)+al(i,j)*x(i,j+1)+ac(i,j)*x(i+1,j+1) & 
                 +ar(i,j)*x(i+2,j+1)+at(i,j)*x(i+1,j+2)
    end do
  end do 

end subroutine p_sqprod
subroutine p_sqprod_T(al,ab,ac,ar,at,x,r,i2,j2)
  ! input & output
  integer, intent(in) :: i2,j2 ! dimensions
  real,intent(in),dimension(i2,j2) ::  &
      ab,al,ac,ar,at
  real,dimension(i2+2,j2+2) ::  &
      tab,tal,tar,tat
  real*8,intent(in),dimension(i2+2,j2+2) :: x
  real*8,intent(inout),dimension(i2+2,j2+2) :: r
  ! local variables 
  integer :: i,j 
  
  tab = 0.
  tal = 0.
  tar = 0.
  tat = 0.
  tab(2:i2+1,2:j2+1) = ab
  tal(2:i2+1,2:j2+1) = al
  tar(2:i2+1,2:j2+1) = ar
  tat(2:i2+1,2:j2+1) = at

  do j = 1,j2
    do i = 1,i2
      r(i+1,j+1) =tat(i+1,j)*x(i+1,j)+tar(i,j+1)*x(i,j+1)+ac(i,j)*x(i+1,j+1) & 
                 +tal(i+2,j+1)*x(i+2,j+1)+tab(i+1,j+2)*x(i+1,j+2)
    end do
  end do 

end subroutine p_sqprod_T
subroutine p_inprod(r,p,tp,i2,j2)
  integer, intent(in) :: i2,j2 ! dimensions
  real*8,intent(in),dimension(i2+2,j2+2) :: r,p
  real*8,intent(inout) :: tp

  !local variables
  integer :: i,j
  tp = 0.d0
  do j = 2,j2+1
    do i = 2,i2+1
      tp = tp + r(i,j)*p(i,j) 
    end do
  end do 

end subroutine p_inprod

subroutine preconditioning(cl,cb,cc,cr,ct,p,ph,i2,j2)
  integer, intent(in) :: i2,j2 ! dimensions
  real,intent(in),dimension(i2,j2) :: cb,cl,cc,cr,ct
  real*8,intent(in),dimension(i2+2,j2+2) :: p
  real*8,intent(inout),dimension(i2+2,j2+2) :: ph

  !local variables
  integer :: i,j
  !!! on preconditioning
  ph(:,:) = p(:,:)
  !!! diagonal preconditioning
  !ph(:,:) = p(:,:)/cc(:,:)

end subroutine
subroutine solver_boundary_update(x,i2,j2)
  integer, intent(in) :: i2,j2 ! dimensions
  real*8,intent(in),dimension(i2+2,j2+2) :: x

  !local variables
  integer :: i,j
#ifdef FLAG_PERIODIC_WE
              DO J=2,J2+1
                X(1,J)=X(I2+1,J)
                X(I0,J)=X(2,J)
              END DO
#endif
end subroutine solver_boundary_update
subroutine printVarSummaryR4(x,m,n,varname)
  integer, intent(in) :: m,n
  real,intent(in),dimension(m,n) :: x
  character(len=8), intent(in) :: varname

  !local variables
  integer :: i,j,imin(2),imax(2)
  real :: fmax,fmin,infinity 

  infinity = huge(0.0d0)
  
  imin = minloc(x)
  fmin = x(imin(1),imin(2))
  imax = minloc(-x)
  fmax = x(imax(1),imax(2))

  write(*,'(2A8,2I4,f15.5)') trim(varname), "min :",imin,fmin 
  write(*,'(2A8,2I4,f15.5)') trim(varname), "max :",imax,fmax
  if (isnan(fmax) .or. isnan(fmin) .or. & 
    fmax .ge. infinity .or. fmin .le.  -infinity ) then
    stop
  end if 

end subroutine printVarSummaryR4

subroutine printVarSummaryR8(x,m,n,varname)
  integer, intent(in) :: m,n
  real*8,intent(in),dimension(m,n) :: x
  character(len=8), intent(in) :: varname

  !local variables
  integer :: i,j,imin(2),imax(2)
  real :: fmax,fmin,infinity 

  infinity = huge(0.0d0)
  
  imin = minloc(x)
  fmin = x(imin(1),imin(2))
  imax = minloc(-x)
  fmax = x(imax(1),imax(2))

  write(*,'(2A8,2I4,f15.5)') trim(varname), "min :",imin,fmin 
  write(*,'(2A8,2I4,f15.5)') trim(varname), "max :",imax,fmax
  if (isnan(fmax) .or. isnan(fmin) .or. & 
    fmax .ge. infinity .or. fmin .le.  -infinity ) then
    stop
  end if 

end subroutine printVarSummaryR8

subroutine figi ( n, t, d, e, e2, ierr )
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) T(N,3) contains the input matrix.  Its subdiagonal
!    is stored in the last N-1 positions of the first column, its diagonal in
!    the N positions of the second column, and its superdiagonal in the
!    first N-1 positions of the third column.  T(1,1) and T(N,3) are arbitrary.
!
!    Output, real ( kind = 8 ) D(N), the diagonal elements of the symmetric 
!    matrix.
!
!    Output, real ( kind = 8 ) E(N), contains the subdiagonal elements of
!    the symmetric matrix in E(2:N).  E(1) is not set.
!
!    Output, real ( kind = 8 ) E2(N), the squares of the corresponding elements 
!    of E.  E2 may coincide with E if the squares are not needed.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, for normal return,
!    N+I, if T(I,1) * T(I-1,3) is negative,
!    -(3*N+I), if T(I,1) * T(I-1,3) is zero with one factor non-zero.  In
!      this case, the eigenvectors of the symmetric matrix are not simply
!      related to those of T and should not be sought.
!
implicit none

integer ( kind = 4 ) n

real ( kind = 8 ) d(n)
real ( kind = 8 ) e(n)
real ( kind = 8 ) e2(n)
integer ( kind = 4 ) i
integer ( kind = 4 ) ierr
real ( kind = 8 ) t(n,3)

ierr = 0

do i = 1, n

  if ( 1 < i ) then

    e2(i) = t(i,1) * t(i-1,3)

    if ( e2(i) < 0.0D+00 ) then
      ierr = n + i
      return
    else if ( e2(i) == 0.0D+00 ) then
      if ( t(i,1) /= 0.0D+00 .or. t(i-1,3) /= 0.0D+00 ) then
        ierr = - 3* n - i
        return
      end if
      e(i) = 0.0D+00

    else

      e(i)   = sqrt ( e2(i))
    end if
  end if

  d(i)  = t(i,2)

end do
return
end

 subroutine  ratqr(n, eps1, d, e, mineig,ierr )
! !DESCRIPTION:
!  This subroutine finds the algebraically smallest 
!  eigenvalue of a positive definite symmetric tridiagonal matrix by the
!  rational QR method with Newton corrections.
!  Follow EISPACK lib subroutine RATQR 
!
! ! INPUT VARIABLES
   integer*4, intent(in) ::  n    ! matrix order
   real*8,intent(in) :: d(n)      ! diagonal elements
   real*8,intent(in) :: e(n)      ! off-diagonal elements, e(1) is arbitrary
   real*8,intent(in) :: eps1      ! convergence tolerance

!  ! OUTPUT VARIABLES
   integer*4, intent(inout) ::  ierr  ! status flag 
   real*8, intent(inout) ::  mineig   !smallest eigevalue

!  ! LOCAL  VARIABLES
   real*8, dimension(n) ::  e2, bd, w
   real*8 :: f, ep, delta, err
   integer*4 i, ii ! local counters
   real*8 p, q, qp,r,s,tot
    
   ierr = 0
   w(1:n) = d(1:n)
   err = 0.0 
   s = 0.0 
   tot = w(1)
   q = 0.0 

   do i = 1, n
      p = q
      bd(i) = e(i)*e(i)
      q = 0.0 

      if ( i /= n ) then
        q = abs ( e(i+1) )
      endif

      tot = min ( w(i) - p - q, tot )
   end do
   bd(1) = 0.0 

   if ( tot < 0.0  ) then
     tot = 0.0 
   else 
     w(1:n) = w(1:n) - tot
   endif

!  QR transformation.

   do while (.true.)
     tot = tot + s
     delta = w(n) - s
     i = n
     if ( delta <= eps1 ) then
       exit
     endif 

     f = bd(n) / delta
     qp = delta + f
     p = 1.0

     do ii = 1, n - 1
       i = n - ii
       q = w(i) - s - f
       r = q / qp
       p = p * r + 1.0
       ep = f * r
       w(i+1) = qp + ep
       delta = q - ep

       ! check convergence 
       if ( delta <= eps1 ) exit

       f = bd(i) / q
       qp = delta + f
       bd(i+1) = qp * ep
     end do
     
     ! check convergence 
     if ( delta <= eps1 ) exit

     w(1) = qp
     s = qp / p

     !  Set error: irregular end of iteration.
     if ( tot + s <= tot ) then
       ierr = 1
       return
     endif

   end do 

   w(1) = tot
   err = err + abs ( delta)
   bd(1) = err
   mineig = w(1)

   return
 end subroutine ratqr
