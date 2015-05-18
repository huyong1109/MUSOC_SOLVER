MODULE PMOVE_GLOBAL 
USE TIMCOM_GENERAL,ONLY : WORKDIR
USE OCN_PARA_GLOBAL
! grid_var.f90
!USE B_GLOBAL       ! for U2 V2
USE GRID_VAR_GLOBAL, ONLY: RHO,Mcrho,McTP,McS,McTPTP,McSTP,P0,PX,PY,U1,U2,V1,V2,S1,S2,T1,T2,P,ULF,VLF,SLF,TLF,EV,HV,F,TANPHI,SUMIN
!USE CGRID_GLOBAL   ! for W
USE GRID_VAR_GLOBAL, ONLY: U,V,W
!USE XYFS_GLOBAL    ! for X0DEG,DXDEG
!USE ZFS_GLOBAL     ! for Z 
USE GRID_VAR_GLOBAL, ONLY: Z,ODZ,ODZW
!USE METRIC_GLOBAL  ! for ODX,ODY,ODZ,YVDEG
USE GRID_VAR_GLOBAL, ONLY: A,XDEG,YV,YVDEG,YDEG,CS,OCS,DX,ODX,DY,ODY,CSV,OCSV,DXV,ODXV,DYV,ODYV
! control.f90
USE CONTROL_GLOBAL ! for DT 
USE SCA_GLOBAL     ! for N360,NYR,IT0,ITF,ITFDAY 
!@USE GRID_VAR_GLOBAL, ONLY: UL,VL,WL
USE DIR_GLOBAL, ONLY: IDIR,ODIR

REAL,PRIVATE,ALLOCATABLE :: UL(:,:,:),VL(:,:,:),WL(:,:,:)

!LOGICAL 
INTEGER,PRIVATE :: L0,LOC                              !total number of tracer
INTEGER,PRIVATE :: LGRP                                !number of tracer group
INTEGER,PARAMETER,PRIVATE :: LGRP_MAX=256              ! maximum tracer group number that system can be process
INTEGER,PRIVATE :: LXSTA(LGRP_MAX),LXEND(LGRP_MAX), &  !store the max and min grid no of x-dim of each tracer group
                   LYSTA(LGRP_MAX),LYEND(LGRP_MAX), &
                   LZSTA(LGRP_MAX),LZEND(LGRP_MAX)
REAL,PRIVATE,ALLOCATABLE :: XG(:),YG(:),ZG(:)          !store the tracer location in logical coordinate 

! scratch arrays
REAL,PRIVATE,ALLOCATABLE,SAVE :: &
XLD(:),XC(:),XNU(:),YLD(:),YC(:),YNU(:), &
ZLD(:),ZC(:),ZNU(:),XCW(:),YCS(:), &
DXCU(:),DYCU(:),DXCV(:),DYCV(:), &
DL(:),XD(:),YD(:),ZM(:), &
USW(:),UES(:),UNW(:),UNE(:),UFL(:), &
VSW(:),VSE(:),VNW(:),VNE(:),VFL(:), &
WSW(:),WSE(:),WNW(:),WNE(:),WFL(:)
INTEGER,PRIVATE,ALLOCATABLE,SAVE :: ICW(:),JCS(:),KC(:)
INTEGER,PRIVATE :: FNO_XYZLOGI=177,FNO_XYZREAL=178
CHARACTER(LEN=256),PRIVATE :: FNM_XYZLOGI,FNM_XYZREAL
! private subroutine
PRIVATE :: PMOVE_RESET, PMOVE_SCRATCH_ARRAY, PMOVE_CONVERT

! public subroutine for calling interface
PUBLIC :: PMOVE_INIT, PMOVE_CHECK, PMOVE_RUN, PMOVE_FINALIZE

CONTAINS

!######################################################################

SUBROUTINE PMOVE_INIT
INTEGER :: ISTAT(3)
  FNM_XYZLOGI=TRIM(WORKDIR)//"/"//TRIM(ODIR)//"/ltracer_logi.txt"
  !#FNM_XYZREAL=TRIM(WORKDIR)//"/"//TRIM(ODIR)//"/ltracer_real.txt"
  OPEN(FNO_XYZLOGI,file=FNM_XYZLOGI,form="formatted")
  !#OPEN(FNO_XYZREAL,file=FNM_XYZREAL,form="formatted")
  allocate (UL(I0,J0,K1),STAT=ISTAT(1))
  allocate (VL(I0,J0,K1),STAT=ISTAT(2))
  allocate (WL(I0,J0,K0),STAT=ISTAT(3))
  if (ISTAT(1)/=0 .or. ISTAT(2)/=0 .or. ISTAT(3)/=0) STOP "Cannot allocate memory"
  UL=0.;VL=0.;WL=0.
  CALL PMOVE_RESET
END SUBROUTINE PMOVE_INIT

!######################################################################

SUBROUTINE PMOVE_FINALIZE
  CLOSE(FNO_XYZLOGI)
  !#CLOSE(FNO_XYZREAL)
END SUBROUTINE PMOVE_FINALIZE

!######################################################################

SUBROUTINE PMOVE_RESET
  INTEGER :: ixs,ixe,iys,iye,izs,ize
  INTEGER :: L
  CHARACTER(256) :: FN
  LGRP=0
  ! "tracer in ltracers.txt must be in format as following:"
  ! <start grid pt in x> <end grid pt in x> <start grid pt in y> <end grid pt in y> <start grid pt in z>
  FN=TRIM(WORKDIR)//"/"//TRIM(IDIR)//"/ltracers.txt"
  open(23,file=TRIM(FN),err=23)
  do while( .true. ) 
   read(23,*,iostat=ios)ixs,ixe,iys,iye,izs,ize
   if ( ios < 0 ) exit
   if ( ios > 0 ) cycle
   LGRP = LGRP + 1
   if ( LGRP > LGRP_MAX ) then
      write(*,*)"number of tracer group can't excced 256!"
      stop
   endif
   LXSTA(LGRP) = ixs
   LXEND(LGRP) = ixe
   LYSTA(LGRP) = iys
   LYEND(LGRP) = iye
   LZSTA(LGRP) = izs
   LZEND(LGRP) = ize
  enddo
  close(23)
  L0=0
  do L = 1, LGRP
     L0 = L0 + (LXEND(L)-LXSTA(L)+1)*(LYEND(L)-LYSTA(L)+1)*(LZEND(L)-LZSTA(L)+1)
  enddo
  if (.not. allocated(XG)) then     !test for allocation
    ALLOCATE(XG(L0),YG(L0),ZG(L0))
  end if

  LOC=0
  do L = 1, LGRP
   do K = LZSTA(L),LZEND(L)
   do J = LYSTA(L),LYEND(L)
   do I = LXSTA(L),LXEND(L)
      LOC = LOC + 1
      XG(LOC) = REAL(I)-1.5
      YG(LOC) = REAL(J)-1.5
      ZG(LOC) = REAL(K)-0.5
   enddo
   enddo
   enddo
  enddo
  if ( LOC /= L0 ) then
    write(*,*)"total number of tracer does NOT consistent!"
    write(*,*)"In module PMOVE_GLOBAL : LOC /= L0"
    write(*,*)"I'd better stop it!"
    stop 1
  endif

  WRITE(FNO_XYZLOGI,7) NYR,N360,LOC,(XG(N),YG(N),ZG(N),N=1,LOC)
7 FORMAT(3I5/(3F7.2))
  CALL FLUSH(FNO_XYZLOGI)
!  WRITE(FNO_XYZREAL,8) NYR,N360,LOC,(XD(N),YD(N),ZM(N),N=1,LOC)
!8 FORMAT(3I5/(2F7.2,F7.1))
!  CALL FLUSH(FNO_XYZREAL)

  CALL PMOVE_SCRATCH_ARRAY
  RETURN 
 23 write(*,*)"error occured when reading ltracers.txt"
  STOP
END SUBROUTINE PMOVE_RESET

!######################################################################

SUBROUTINE PMOVE_SCRATCH_ARRAY
integer,save :: L0_CUR = -1

if ( (L0_CUR /= L0) .and. (L0 > 0) ) then
   if (allocated(XLD)) then
      deallocate &
      (XLD,XC,XNU,YLD,YC,YNU, &
       ZLD,ZC,ZNU,XCW,YCS, &
       DXCU,DYCU,DXCV,DYCV,DL,XD,YD,ZM, &
       USW,UES,UNW,UNE,UFL, &
       VSW,VSE,VNW,VNE,VFL, &
       WSW,WSE,WNW,WNE,WFL, &
       ICW,JCS,KC)
   endif
   L0_CUR = L0
   allocate &
      (XLD(L0),XC(L0),XNU(L0),YLD(L0),YC(L0),YNU(L0), &
       ZLD(L0),ZC(L0),ZNU(L0),XCW(L0),YCS(L0), &
       DXCU(L0),DYCU(L0),DXCV(L0),DYCV(L0),DL(L0),XD(L0),YD(L0),ZM(L0), &
       USW(L0),UES(L0),UNW(L0),UNE(L0),UFL(L0), &
       VSW(L0),VSE(L0),VNW(L0),VNE(L0),VFL(L0), &
       WSW(L0),WSE(L0),WNW(L0),WNE(L0),WFL(L0), &
       ICW(L0),JCS(L0),KC(L0))
endif
END SUBROUTINE PMOVE_SCRATCH_ARRAY 

!######################################################################

SUBROUTINE PMOVE_CHECK(NEWDAY)
  INTEGER :: NEWDAY

! NEWDAY = 0 : timestep is still in old day
! NEWDAY = 1 : timestep is update to new day, it might be output the tracer location
! NEWDAY =-1 : new tracers need to be re-put 

! new tracer diagnostics
  NEWDAY=0
  IF (MOD(ITF,ITFDAY).EQ.1) NEWDAY=1
  IF (ITF.EQ.IT0+1.OR.MOD(ITF,30*ITFDAY).EQ.1) NEWDAY=-1

! reset or initialize tracer locations if NEWDAY = -1
! save tracer locations if NEWDAY /= 0 ( = 1 )
! keep walking if NEWDAY = 0
  IF (NEWDAY.EQ.-1) THEN
    CALL PMOVE_RESET
  ENDIF

  !IF (NEWDAY.NE.0) THEN
  IF (NEWDAY.EQ.1) THEN
    ! convert logical to real coordinate
    !#  CALL PMOVE_CONVERT
    ! write initial or present tracer locations at start of new day
    ! write LOGICAL grid locations
      WRITE(FNO_XYZLOGI,7) NYR,N360,LOC,(XG(N),YG(N),ZG(N),N=1,LOC)
    7 FORMAT(3I5/(3F7.2))
      CALL FLUSH(FNO_XYZLOGI)
   !   WRITE(FNO_XYZREAL,8) NYR,N360,LOC,(XD(N),YD(N),ZM(N),N=1,LOC)
   ! 8 FORMAT(3I5/(2F7.2,F7.1))
   !   CALL FLUSH(FNO_XYZREAL)

  ENDIF

END SUBROUTINE PMOVE_CHECK

!######################################################################

SUBROUTINE PMOVE_CONVERT
! convert logical coordinate to real spherical coordinate
    ! write LOGICAL and SPHERICAL coordinate locations daily
    ! XG,YG,ZG is logical coordinate
    ! XD,YD,ZM is real spherical coordinate
      DO N=1,L0
      ! convert to deg lon (X=0 at west face of westernmost INTERIOR zone)
        !XD(N)=WESTDEG+XG(N)*DXDEG
        !XD(N)=X0DEG+XG(N)*DXDEG
        XD(N)=XG(N)*DXDEG
      ! convert to deg lat (Y=0 at south face of southernmost INTERIOR zone)
        TMP=YG(N)+1.
        J=TMP
        YD(N)=YVDEG(J)+(TMP-J)*(YVDEG(J+1)-YVDEG(J))
      ! convert to meters. However, logical vertical coordinate may be better!
      ! NOTES
      ! Z is LOGICAL GRID tracer location (Z=0 at rigid lid)
      ! K is the layer containing the tracer
      ! ZZ(2*K-1) is first z-level layer interface above the tracer
      ! TMP is the fraction (0<TMP<1) toward the next layer interface ZZ(2*K+1)
      ! ZM is distance from the tracer to the rigid lid in meters.
        K=ZG(N)+1.
        TMP=ZG(N)+1.-K
        ZM(N)=.01*(Z(2*K-1)+TMP*(Z(2*K+1)-Z(2*K-1)))
      ENDDO
END SUBROUTINE PMOVE_CONVERT

!######################################################################

! ======================================================================
SUBROUTINE PMOVE_RUN(MXITP)
! PM moves tracers in a 3-d lagrangian manner
! ======================================================================
! MXITP= number of iterations used to get new tracer locations
! L0= size of tracer and scratch arrays
!     used to organize computation into sequential array operations
!     L0.GE.(number of tracers to be followed)
! UL,VL,WL are at cell center locations
! UL,VL,WL are XG,YG,ZG LOGICAL displacements(grid intervals per time step)
! UL,VL,WL are interpolated to logical tracer locations to get tracer
! displacements in fractional grid intervals per time step
! XG,YG,ZG are tracer LOGICAL grid locations, and are saved on disk only at
! Z is z-level grid array
! XG is logical distance from west face of the westernmost INTERIOR cell
! "logical locations" are in grid increments
! end of each model day

!HANGE UNITS to fractional grid interval displacement per time interval
! time interval is DT/2, tracer leapfrog eliminates extra 0.5 factor.
  DO K=1,K1
    DO J=1,J0
      TMP1=DT*ODX(J)
      TMP2=DT*ODY(J)
      TMP3=.5*DT*ODZ(K)
      DO I=1,I0
      ! mid-layer cell center values
      ! assume vertically invarient flow in each layer
      ! to avoid needing vertical interpolations
        UL(I,J,K)=TMP1*U2(I,J,K)
        VL(I,J,K)=TMP2*V2(I,J,K)
        WL(I,J,K)=TMP3*.5*(W(I,J,K)+W(I,J,K+1))
      ENDDO
    ENDDO
  ENDDO

! "initial" tracer location at start of time step
  DO N=1,LOC
    ZLD(N)=ZG(N)
    XLD(N)=XG(N)
    YLD(N)=YG(N)
! initial guess for tracer location at end of time step
    ZNU(N)=ZLD(N)
    XNU(N)=XLD(N)
    YNU(N)=YLD(N)
  ENDDO
! iterate to improve accuracy
  DO ITP=1,MXITP
    DO N=1,LOC
    ! two-level time average tracer location
    ! fluid velocity is at same time level
      ZC(N)=.5*(ZLD(N)+ZNU(N))
      XC(N)=.5*(XLD(N)+XNU(N))
      YC(N)=.5*(YLD(N)+YNU(N))
    ENDDO
  ! no need to re-initialize tracers near dry points or open boundaries
    DO N=1,LOC
    ! first scalar cell center west of tracer #N
    ! ghost zones have value either IC=1, IC=I0, JC=1 or JC=J0
    ! westernmost interior cell center is at XC=0.5
    ! 0 < XC < I2 and 0 < YC < J2 for interior points
      ICW(N)=XC(N)+1.5
      JCS(N)=YC(N)+1.5
    ! no ghost zones in vertical, so it's simplest to take nearest z-level
    ! to avoid interpolating data from outside vertical domain
    ! (ok if vertical variations of velocity are not big)
    ! 0 < ZC < K1
      KC(N)=ZC(N)+1.
    ! nearest cell center line west of tracer
      XCW(N)=ICW(N)-1.5
    ! nearest cell center line south of tracer
      YCS(N)=JCS(N)-1.5
    ENDDO
! ZC not used because no vertical interpolations are used
! velocity at cell corners
    DO N=1,LOC
      I=ICW(N)
      J=JCS(N)
      K=KC(N)
    ! cell corners
      USW(N)=UL(I,J,K)
      UES(N)=UL(I+1,J,K)
      UNW(N)=UL(I,J+1,K)
      UNE(N)=UL(I+1,J+1,K)
      VSW(N)=VL(I,J,K)
      VSE(N)=VL(I+1,J,K)
      VNW(N)=VL(I,J+1,K)
      VNE(N)=VL(I+1,J+1,K)
    ! W is at cell center, and assumed constant vertically within cell
      WSW(N)=WL(I,J,K)
      WSE(N)=WL(I+1,J,K)
      WNW(N)=WL(I,J+1,K)
      WNE(N)=WL(I+1,J+1,K)
    ENDDO
! interpolate nearest corner values to tracer location
    DO N=1,LOC
      DXCU(N)=XC(N)-XCW(N)
      DYCU(N)=YC(N)-YCS(N)
      DXCV(N)=XC(N)-XCW(N)
      DYCV(N)=YC(N)-YCS(N)
      DL(N)=UNW(N)-USW(N)
! add z-interpolations!
      UFL(N)=USW(N)+(UES(N)-USW(N))*DXCU(N)+DYCU(N)* &
            (DL(N)+(UNE(N)-UES(N)-DL(N))*DXCU(N))
      DL(N)=VNW(N)-VSW(N)
      VFL(N)=VSW(N)+(VSE(N)-VSW(N))*DXCV(N)+DYCV(N)* &
            (DL(N)+(VNE(N)-VSE(N)-DL(N))*DXCV(N))
      DL(N)=WNW(N)-WSW(N)
      WFL(N)=WSW(N)+(WSE(N)-WSW(N))*DXCV(N)+DYCU(N)* &
            (DL(N)+(WNE(N)-WSE(N)-DL(N))*DXCV(N))
    ENDDO
    DO N=1,LOC
      ZNU(N)=ZLD(N)+WFL(N)
      XNU(N)=XLD(N)+UFL(N)
      YNU(N)=YLD(N)+VFL(N)
    ENDDO
  ENDDO
  DO N=1,L0
    ZG(N)=.5*(ZNU(N)+ZLD(N))
    XG(N)=.5*(XNU(N)+XLD(N))
    YG(N)=.5*(YNU(N)+YLD(N))
  ENDDO
END SUBROUTINE PMOVE_RUN


END MODULE PMOVE_GLOBAL 

