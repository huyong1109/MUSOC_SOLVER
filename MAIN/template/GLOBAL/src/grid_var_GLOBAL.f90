MODULE OCN_PARA_GLOBAL
#include "namelist.in_GLOBAL"
USE TIMCOM_GENERAL
IMPLICIT NONE
INTEGER        ::I0,I1,I2,I3,I4
INTEGER        ::J0,J1,J2,J3
INTEGER        ::K0,K1,K2,K01
INTEGER        ::NB0,NB1,NBIR,IBIR,IBR2,N0
type(grid2D) ::localgrid2D
type(field2d)::localSSTfield,localHeatflux,localSolarHeatflux, &
     localMomentfluxX,localMomentfluxY,localPrecipitation,localAtmPressure, &
     localWindU10m,localWindV10m,localAirTemperature2m,localDewTemperature2m, &
     localCloudCoverage,localQsolar,localQup, localEvap
logical :: bulk_formulae, hum_type
real,allocatable,dimension(:,:) :: tempU10M,tempV10M,tempT2M,tempDT2M,tempPMSL,tempCC

CONTAINS


SUBROUTINE VAR_SPECIFY_OCN_PARA
  INTEGER :: I,J,K
  INTEGER :: STATUS(4)

  I1=I0-1;I2=I0-2;I3=I0-3;I4=I0-4
  J1=J0-1;J2=J0-2;J3=J0-3
  K1=K0-1;K2=K0-2;K01=K0+K1
  NB0=ANINT(REAL(J2)/REAL(N0))
  IBIR=2*I2/NBIR;IBR2=IBIR*IBIR
  NB1=NB0-1
END SUBROUTINE VAR_SPECIFY_OCN_PARA

END MODULE OCN_PARA_GLOBAL


!######################################################################


MODULE WINDID_GLOBAL
INTEGER :: WINDID
END MODULE WINDID_GLOBAL


!######################################################################


MODULE GRID_VAR_GLOBAL
!Horizontal resolution parameters are I0,I1,I2,I3, J0,J1,J2,J3
!Vertical resolution parameters are K0,K1,K2
!NB0,NB1 are dimension parameters for SEVP elliptic solver

! "A" grid arrays
USE OCN_PARA_GLOBAL

IMPLICIT NONE

SAVE

REAL :: AROUT

REAL(8),ALLOCATABLE :: RHO(:,:,:),Mcrho(:,:,:),McTP(:,:,:),McS(:,:,:), &
                       McTPTP(:,:,:),McSTP(:,:,:)
REAL,ALLOCATABLE :: P0(:,:),PX(:,:),PY(:,:)
REAL,ALLOCATABLE :: U1(:,:,:),U2(:,:,:),V1(:,:,:),V2(:,:,:),S1(:,:,:),S2(:,:,:),T1(:,:,:),T2(:,:,:),P(:,:,:),ULF(:,:,:),VLF(:,:,:),SLF(:,:,:),TLF(:,:,:)
REAL,ALLOCATABLE :: EV(:,:,:),HV(:,:,:)
REAL,ALLOCATABLE :: F(:),TANPHI(:),SUMIN(:)

#ifdef FLAG_TRACER
REAL,ALLOCATABLE::C1(:,:,:),C2(:,:,:),CLF(:,:,:) 
#endif
!use OCN_PARA_GLOBAL,ONLY:I0,J0,K1

!use OCN_PARA_GLOBAL,ONLY:I0,J0,J1,K0,K1
REAL,ALLOCATABLE :: U(:,:,:),V(:,:,:),W(:,:,:)

!use OCN_PARA_GLOBAL,ONLY:K0,K1
! Vertical grid arrays
! "Z" array contains cell center AND interface depths
real :: myZ
REAL,ALLOCATABLE :: Z(:),ODZ(:),ODZW(:)
REAL,ALLOCATABLE :: transmitted(:,:,:)

!use OCN_PARA_GLOBAL,ONLY:I0,J0,J1,K0,K1
! Lower precision logical depth and masking arrays
! Logical depth KB and associated masking arrays
INTEGER(2),ALLOCATABLE :: KB(:,:),IU0(:,:),IV0(:,:)
INTEGER(2),ALLOCATABLE :: IN(:,:,:),IU(:,:,:)
INTEGER(2),ALLOCATABLE :: IV(:,:,:)
INTEGER(2),ALLOCATABLE :: IW(:,:,:)

INTEGER(2),ALLOCATABLE :: IFLT(:,:,:)

!use OCN_PARA_GLOBAL,ONLY:I0,J0,J1,K1
REAL,ALLOCATABLE :: A(:)
REAL,ALLOCATABLE :: XDEG(:)
REAL,ALLOCATABLE :: YV(:),YVDEG(:),YDEG(:),CS(:),OCS(:),DX(:),ODX(:),DY(:),ODY(:)
REAL,ALLOCATABLE :: CSV(:),OCSV(:),DXV(:),ODXV(:),DYV(:),ODYV(:)

REAL,ALLOCATABLE :: Y(:)!(J0)

! Time averaged arrays
REAL,ALLOCATABLE :: PBAR(:,:),PVAR(:,:),XBAR(:,:),UCLI(:,:),VCLI(:,:),RMSV(:,:) !(I0,J0)
REAL,ALLOCATABLE :: SBAR(:,:,:),TBAR(:,:,:) !(I0,J0,K1)

! time avg BCS for duo grid NAB model
! split the following data between GLO2NAB and GLO2IBE
! and add NAB2GLO and IBE2GLO common blocks
REAL,ALLOCATABLE :: VGLO(:,:),UAVG(:,:),VAVG(:,:),SAVG(:,:),TAVG(:,:) !(I0,K1)

! Scratch arrays
REAL,ALLOCATABLE :: SCR(:,:,:) !(I0,J0,K0+1)

!USE OCN_PARA_GLOBAL,ONLY:I0,J0,I1,J1,K1,K2
REAL,ALLOCATABLE :: DMX(:,:,:),DMY(:,:,:)
REAL,ALLOCATABLE :: VBK(:,:,:),HBK(:)
!FILTERGLO
INTEGER(2),ALLOCATABLE :: INFX(:,:,:),INFY(:,:,:)

!use OCN_PARA_GLOBAL,ONLY:K1
REAL,ALLOCATABLE :: SXY(:),TXY(:)
REAL,ALLOCATABLE :: SXYCLI(:,:),TXYCLI(:,:)

!use OCN_PARA_GLOBAL,ONLY:K1
!INTEGER,DIMENSION( 8)::ILOC
!INTEGER,DIMENSION(10)::JLOC
INTEGER,ALLOCATABLE :: KP(:)

!use OCN_PARA_GLOBAL,ONLY:I0,J0,I2,J2,K1,K2
! Derived model climate source arrays
! Annual cycle sources in all sponge layer locations
REAL,ALLOCATABLE :: TNUDGE(:,:)
REAL,ALLOCATABLE :: QAVG(:,:,:),WAVG(:,:,:)
REAL,ALLOCATABLE :: SCLI(:,:,:),TCLI(:,:,:),QDAVG(:,:,:),SDAVG(:,:,:)
!INTEGER,DIMENSION(12)::NSOMBO
! climate data for winter
! northern and southern sponge layer Levitus Climatology
! COMMON BLOCK/KLIMATGLO/
! COMMON BLOCK/BCSGLO/
!REAL,DIMENSION(12)::TSURFM,SSURFM
REAL,ALLOCATABLE :: SSURF(:,:,:),TSURF(:,:,:)
REAL,ALLOCATABLE :: SSSP(:,:,:,:),SNSP(:,:,:,:),TSSP(:,:,:,:),TNSP(:,:,:,:)

! nudging control matrix
integer*2, allocatable :: nudge(:,:,:)

REAL,ALLOCATABLE :: TAUX(:,:,:),TAUY(:,:,:),taux1(:,:),tauy1(:,:),qdot(:,:),qdot2(:,:),rain(:,:),evapo(:,:)

REAL(8),ALLOCATABLE::RINV(:,:,:)
REAL(8),ALLOCATABLE::RINV1(:,:,:)
REAL(8),ALLOCATABLE::DUM0(:,:)
REAL(8),ALLOCATABLE::DUM1(:),DUM2(:)
REAL(8),ALLOCATABLE::XX(:,:),H(:,:)
REAL(8),ALLOCATABLE::X(:,:)
REAL,ALLOCATABLE::S(:,:),AL(:,:),AB(:,:),AC(:,:),OAC(:,:),AR(:,:),AT(:,:)
REAL,ALLOCATABLE::SRC(:,:),CL(:,:),CB(:,:),CC(:,:),OCC(:,:), CR(:,:),CT(:,:)
INTEGER,ALLOCATABLE ::IE(:)

REAL,ALLOCATABLE :: DEPTH(:,:)

!PG_GLOBAL
REAL,ALLOCATABLE :: DHX(:,:),DHY(:,:),XP(:,:,:),YP(:,:)
REAL,ALLOCATABLE :: tlev(:,:,:),slev(:,:,:),ssp(:,:,:),tsp(:,:,:),ADD(:,:,:)

CONTAINS

SUBROUTINE VAR_ALLOCATE_B

INTEGER :: STATUS(4)

allocate (RHO(I0,J0,K1),STAT=status(1)); allocate (P0(I0,J0),STAT=status(2)); allocate (PX(I0,J0),STAT=status(3)); allocate (PY(I0,J0),STAT=status(4))

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot allocate memory"

allocate (Mcrho(I0,J0,K1),STAT=status(1)); allocate (McTP(I0,J0,K1),STAT=status(2)); allocate (McS(I0,J0,K1),STAT=status(3)); allocate (McTPTP(I0,J0,K1),STAT=status(4))

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot allocate memory"

allocate (McSTP(I0,J0,K1),STAT=status(1))
if (status(1)/=0) STOP "Cannot allocate memory"


allocate (U1(I0,J0,K1),STAT=status(1)); allocate (U2(I0,J0,K1),STAT=status(2)); allocate (V1(I0,J0,K1),STAT=status(3)); allocate (V2(I0,J0,K1),STAT=status(4))

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot allocate memory"

allocate (S1(I0,J0,K1),STAT=status(1)); allocate (S2(I0,J0,K1),STAT=status(2)); allocate (T1(I0,J0,K1),STAT=status(3)); allocate (T2(I0,J0,K1),STAT=status(4))

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot allocate memory"

allocate (P(I0,J0,K1),STAT=status(1)); allocate (ULF(I0,J0,K1),STAT=status(2)); allocate (VLF(I0,J0,K1),STAT=status(3)); allocate (SLF(I0,J0,K1),STAT=status(4))

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot allocate memory"

allocate (TLF(I0,J0,K1),STAT=status(1)); allocate (EV(I2,J2,K2),STAT=status(2)); allocate (HV(I2,J2,K2),STAT=status(3))

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0) STOP "Cannot allocate memory"

allocate (F(J2),STAT=status(1)); allocate (TANPHI(J2),STAT=status(2)); allocate (SUMIN(K1),STAT=status(3))

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0) STOP "Cannot allocate memory"

allocate (DEPTH(I0,J0),STAT=status(1))

if (status(1)/=0) STOP "Cannot allocate memory"


RHO=0.;P0=0.;PX=0.;PY=0.;U1=0.;U2=0.;V1=0.;V2=0.;S1=0.;S2=0.;T1=0.;T2=0.;P=0.;ULF=0.;VLF=0.;SLF=0.;TLF=0.;EV=0.;HV=0.;F=0.;TANPHI=0.;SUMIN=0.;DEPTH=0.;Mcrho=0.d0;McTP=0.d0;McS=0.d0;McTPTP=0.d0;McSTP=0.d0

END SUBROUTINE VAR_ALLOCATE_B

!======================================================

SUBROUTINE VAR_ALLOCATE_GLOBAL

INTEGER :: STATUS(4)

#ifdef FLAG_TRACER
allocate (C1(I0,J0,K1),STAT=status(1)); allocate (C2(I0,J0,K1),STAT=status(2)); allocate (CLF(I0,J0,K1),STAT=status(3))
if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0) STOP "Cannot allocate memory"
C1=0.;C2=0.;CLF=0.
#endif
allocate (IFLT(I0,J0,K1),STAT=status(2))
if (status(1)/=0 .or. status(2)/=0) STOP "Cannot allocate memory"
IFLT=0.

END SUBROUTINE VAR_ALLOCATE_GLOBAL

!======================================================

SUBROUTINE VAR_ALLOCATE_CGRID

INTEGER :: STATUS(4)

allocate (U(I0,J0,K1),STAT=status(1)); allocate (V(I0,J1,K1),STAT=status(2)); allocate (W(I0,J0,K0),STAT=status(3))

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0) STOP "Cannot allocate memory"

U=0.;V=0.;W=0.

END SUBROUTINE VAR_ALLOCATE_CGRID

!======================================================

SUBROUTINE VAR_ALLOCATE_ZFS

INTEGER :: STATUS(3)

allocate (Z(K0+K1),STAT=status(1)); allocate (ODZ(K1),STAT=status(2)); allocate (ODZW(K0),STAT=status(3))

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0) STOP "Cannot allocate memory during VAR_ALLOCATE_ZFS"

Z=0.;ODZ=0.;ODZW=0.

END SUBROUTINE VAR_ALLOCATE_ZFS

!============================================================

SUBROUTINE VAR_ALLOCATE_ABSORPTION
INTEGER :: STATUS
allocate(transmitted(i2,j2,k1),stat=status);
if (status/=0) STOP "Cannot allocate memory for transmitted"
END SUBROUTINE VAR_ALLOCATE_ABSORPTION

!=============================================================

SUBROUTINE CALC_ABSORPTION_RATE
integer :: i,j,k,n
real::R=0.58,xi1=0.35,xi2=23
real(kind=8)::total_transmitted

do i=1,i2
do j=1,j2

  total_transmitted=0.0
  do k=1,k1
    myZ=-z(2*K)/100.      ! convert cm to m , make depth negative
    transmitted(i,j,k)=(R*exp(myZ/xi1)+(1-R)*exp(myZ/xi2)) * in(i+1,j+1,k)*odz(k)
    total_transmitted=total_transmitted+transmitted(i,j,k)
  end do

  if (total_transmitted.gt.1e-6) then
    do k=1,k1
      transmitted(i,j,k)=transmitted(i,j,k)/total_transmitted
    end do    
  else
      transmitted(i,j,:)=0.0
  end if
end do
end do

END SUBROUTINE CALC_ABSORPTION_RATE

!=============================================================

SUBROUTINE VAR_ALLOCATE_BATHY

INTEGER :: STATUS(4)

allocate (KB(I0,J0),STAT=status(1)); allocate (IU0(I0,J0),STAT=status(2)); allocate (IV0(I0,J0),STAT=status(3))

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0) STOP "Cannot allocate memory"

allocate (IN(I0,J0,K1),STAT=status(1)); allocate (IU(I0,J0,K1),STAT=status(2)); allocate (IV(I0,J1,K1),STAT=status(3)); allocate (IW(I0,J0,K0),STAT=status(4))

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot allocate memory"

KB=0;IU0=0;IV0=0;IN=0;IU=0;IV=0;IW=0

END SUBROUTINE VAR_ALLOCATE_BATHY

!=======================================================

SUBROUTINE VAR_ALLOCATE_METRIC

INTEGER :: STATUS(4)

allocate (A(K1),STAT=status(1)); allocate (XDEG(I0),STAT=status(2)); allocate (YV(J0),STAT=status(3)); allocate (YVDEG(J0),STAT=status(4))

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot allocate memory"

allocate (YDEG(J0),STAT=status(1)); allocate (CS(J0),STAT=status(2)); allocate (OCS(J0),STAT=status(3)); allocate (DX(J0),STAT=status(4))

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot allocate memory"

allocate (ODX(J0),STAT=status(1)); allocate (DY(J0),STAT=status(2)); allocate (ODY(J0),STAT=status(3))

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0) STOP "Cannot allocate memory"

allocate (CSV(J1),STAT=status(1)); allocate (OCSV(J1),STAT=status(2)); allocate (DXV(J1),STAT=status(3))

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0) STOP "Cannot allocate memory"

allocate (ODXV(J1),STAT=status(1)); allocate (DYV(J1),STAT=status(2)); allocate (ODYV(J1),STAT=status(3))

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0) STOP "Cannot allocate memory"

A=0.;XDEG=0.;YV=0.;YVDEG=0.;YDEG=0.;CS=0.;OCS=0.;DX=0.;ODX=0.;DY=0.;CSV=0.;OCSV=0.;DXV=0.;ODXV=0.;DYV=0.;ODYV=0.;ODY=0.

END SUBROUTINE VAR_ALLOCATE_METRIC

!======================================================

SUBROUTINE VAR_ALLOCATE_METRIC_PRE

INTEGER :: STATUS(4)

allocate (Y(J0),STAT=status(1))

if (status(1)/=0) STOP "Cannot allocate memory"

Y=0.

END SUBROUTINE VAR_ALLOCATE_METRIC_PRE

!======================================================

SUBROUTINE VAR_ALLOCATE_TAVG_MAIN

INTEGER :: STATUS(4)

allocate (PBAR(I0,J0),STAT=status(1)); allocate (PVAR(I0,J0),STAT=status(2)); allocate (XBAR(I0,J0),STAT=status(3)); allocate (UCLI(I0,J0),STAT=status(4))

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot allocate memory"

allocate (VCLI(I0,J0),STAT=status(1)); allocate (RMSV(I0,J0),STAT=status(2)); allocate (SBAR(I0,J0,K1),STAT=status(3)); allocate (TBAR(I0,J0,K1),STAT=status(4))

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot allocate memory"

PBAR=0.;PVAR=0.;XBAR=0.;UCLI=0.;VCLI=0.;RMSV=0.;SBAR=0.;TBAR=0.

END SUBROUTINE VAR_ALLOCATE_TAVG_MAIN

!======================================================

SUBROUTINE VAR_ALLOCATE_GLO2NAB_MAIN

INTEGER :: STATUS(4)

allocate (VGLO(I0,K1),STAT=status(1)); allocate (UAVG(I0,K1),STAT=status(2)); allocate (VAVG(I0,K1),STAT=status(3)); allocate (SAVG(I0,K1),STAT=status(4))

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot allocate memory"

allocate (TAVG(I0,K1),STAT=status(1))

if (status(1)/=0) STOP "Cannot allocate memory"

VGLO=0.;UAVG=0.;VAVG=0.;SAVG=0.;TAVG=0.

END SUBROUTINE VAR_ALLOCATE_GLO2NAB_MAIN

!======================================================

SUBROUTINE VAR_ALLOCATE_SCRATCH_MAIN

INTEGER :: STATUS(4)

allocate (SCR(I0,J0,K0+1),STAT=status(1))

if (status(1)/=0) STOP "Cannot allocate memory"

SCR=0.

END SUBROUTINE VAR_ALLOCATE_SCRATCH_MAIN

!======================================================

SUBROUTINE VAR_ALLOCATE_WINDMX_MAIN

INTEGER :: STATUS(4)

allocate (DMX(I1,J0,K1),STAT=status(1)); allocate (DMY(I0,J1,K1),STAT=status(2))

if (status(1)/=0 .or. status(2)/=0) STOP "Cannot allocate memory"

allocate (HBK(K2),STAT=status(2)); allocate (INFX(I0,J0,K1),STAT=status(3)); allocate (INFY(I0,J0,K1),STAT=status(4))

if (status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot allocate memory"
!vbk(:,:,:) will be allocated when read from namelist, in module input_XXX , subroutine read_input_parameters

DMX=0.;DMY=0.;HBK=0.;INFX=0;INFY=0

END SUBROUTINE VAR_ALLOCATE_WINDMX_MAIN

!======================================================

SUBROUTINE VAR_ALLOCATE_XYMEANS_MAIN

INTEGER :: STATUS(4)

allocate (SXY(K1),STAT=status(1)); allocate (TXY(K1),STAT=status(2)); allocate (SXYCLI(K1,12),STAT=status(3)); allocate (TXYCLI(K1,12),STAT=status(4))

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot allocate memory"

SXY=0.;TXY=0.;SXYCLI=0.;TXYCLI=0.

END SUBROUTINE VAR_ALLOCATE_XYMEANS_MAIN

!======================================================

SUBROUTINE VAR_ALLOCATE_IJKPLOT_MAIN

INTEGER :: STATUS(4)


allocate (KP(K1),STAT=status(1))

if (status(1)/=0) STOP "Cannot allocate memory"


KP=0

KP=(/1,2*0,1,3*0,1,3*0,1,3*0,1,3*0,1,10*0/)

END SUBROUTINE VAR_ALLOCATE_IJKPLOT_MAIN

!======================================================

SUBROUTINE VAR_ALLOCATE_CLIMAT

INTEGER :: STATUS(4)

allocate (TNUDGE(I2,J2),STAT=status(1)); allocate (QAVG(I2,J2,12),STAT=status(2)); allocate (WAVG(I2,J2,12),STAT=status(3))

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0) STOP "Cannot allocate memory"

allocate (SCLI(I0,J0,K1),STAT=status(1)); allocate (TCLI(I0,J0,K1),STAT=status(2)); allocate (QDAVG(I0,J0,K1),STAT=status(3)); allocate (SDAVG(I0,J0,K1),STAT=status(4))

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot allocate memory"

allocate (SSURF(I2,J2,12),STAT=status(3)); allocate (TSURF(I2,J2,12),STAT=status(4))

if (status(3)/=0 .or. status(4)/=0) STOP "Cannot allocate memory"

allocate (SSSP(I2,10,K1,12),STAT=status(1)); allocate (SNSP(I2,10,K1,12),STAT=status(2)); allocate (TSSP(I2,10,K1,12),STAT=status(3)); allocate (TNSP(I2,10,K1,12),STAT=status(4))

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot allocate memory"

TNUDGE=0.;QAVG=0.;WAVG=0.;SCLI=0.;TCLI=0.;QDAVG=0.;SDAVG=0.;SSURF=0.;TSURF=0.;SSSP=0.;SNSP=0.;TSSP=0.;TNSP=0.

allocate(nudge(i0,j0,k1),stat=status(1))
if (status(1)/=0) STOP "cannot allocate memory for variable nudge(i,j,k)"

END SUBROUTINE VAR_ALLOCATE_CLIMAT

!======================================================

SUBROUTINE VAR_ALLOCATE_CLIMAT2

INTEGER :: STATUS(4)

allocate (SSURF(I2,J2,12),STAT=status(3)); allocate (TSURF(I2,J2,12),STAT=status(4))

if (status(3)/=0 .or. status(4)/=0) STOP "Cannot allocate memory"

allocate (SSSP(I2,10,K1,12),STAT=status(1)); allocate (SNSP(I2,10,K1,12),STAT=status(2)); allocate (TSSP(I2,10,K1,12),STAT=status(3)); allocate (TNSP(I2,10,K1,12),STAT=status(4))

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot allocate memory"

SSURF=0.;TSURF=0.;SSSP=0.;SNSP=0.;TSSP=0.;TNSP=0.

END SUBROUTINE VAR_ALLOCATE_CLIMAT2

!======================================================

SUBROUTINE VAR_ALLOCATE_TAUXY
use timcom_general, only: windtype

INTEGER :: STATUS(4)

if (windtype.eq.-1) then
  allocate (TAUX(I2,J2,12),STAT=status(3)); allocate (TAUY(I2,J2,12),STAT=status(4))
  if (status(3)/=0 .or. status(4)/=0) STOP "Cannot allocate memory"
  TAUX=0.;TAUY=0.
elseif (windtype.ge.0) then
  allocate(TAUX1(I2,J2),stat=status(3)); allocate(tauy1(I2,J2),stat=status(4))
  if (status(3)/=0 .or. status(4)/=0) STOP "Cannot allocate memory"
end if

if (heattype.ge.0.or.windtype.ge.0.or.salttype.ge.0) then
  allocate(qdot(I2,J2),stat=status(3))
  allocate(qdot2(I2,J2),stat=status(4))
  if (status(3)/=0 .or. status(4)/=0) STOP "Cannot allocate memory"
end if

if (raintype.eq.1) then
  allocate(rain(I2,J2),stat=status(4))
  if (status(4).ne.0) STOP "Cannot allocate memory"
end if
if (raintype.eq.1.or.salttype.ge.0.or.windtype.ge.0.or.heattype.ge.0) then
  allocate(evapo(I2,J2),stat=status(4))
  if (status(4).ne.0) STOP "Cannot allocate memory"
end if

END SUBROUTINE VAR_ALLOCATE_TAUXY

!======================================================

SUBROUTINE VAR_ALLOCATE_SEVP
USE CONTROL_GLOBAL

INTEGER :: STATUS(4)

IF (FL_EVP_STP == 1) THEN

  allocate (RINV(IBR2,NB0,NBIR),STAT=status(1)); allocate (RINV1(IBR2,NB1,NBIR),STAT=status(2)); allocate (DUM0(IBIR,NB1),STAT=status(3)) !new
  allocate (DUM1(IBIR),STAT=status(1)); allocate (DUM2(IBIR),STAT=status(2)); allocate (XX(IBIR+2,J0),STAT=status(3)); allocate (H(IBIR+2,J0),STAT=status(4))

ELSE if (FL_EVP_STP == 0) then

  allocate (RINV(I2,I2,NB0),STAT=status(1)); allocate (RINV1(I2,I2,NB1),STAT=status(2)); allocate (DUM0(I2,NB1),STAT=status(3)) !old
  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0) STOP "Cannot allocate memory"
  allocate (DUM1(I2),STAT=status(1)); allocate (DUM2(I2),STAT=status(2)); allocate (XX(I0,J0),STAT=status(3)); allocate (H(I0,J0),STAT=status(4))

ELSE if (FL_EVP_STP == 2) then  !!! hyedit, save for bicg
  allocate (RINV(IBR2,NB0,NBIR),STAT=status(1)); allocate (RINV1(IBR2,NB1,NBIR),STAT=status(2)); allocate (DUM0(IBIR,NB1),STAT=status(3)) !new
  allocate (DUM1(IBIR),STAT=status(1)); allocate (DUM2(IBIR),STAT=status(2)); allocate (XX(IBIR+2,J0),STAT=status(3)); allocate (H(IBIR+2,J0),STAT=status(4))

END IF

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot allocate memory"

allocate (X(I0,J0),STAT=status(1))

if (status(1)/=0) STOP "Cannot allocate memory"

allocate (S(I2,J2),STAT=status(1)); allocate (AL(I2,J2),STAT=status(2)); allocate (AB(I2,J2),STAT=status(3)); allocate (AC(I2,J2),STAT=status(4))
if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot allocate memory"

allocate (OAC(I2,J2),STAT=status(1)) 
if (status(1)/=0) STOP "Cannot allocate memory"

allocate (AR(I2,J2),STAT=status(1)); allocate (AT(I2,J2),STAT=status(2))

if (status(1)/=0 .or. status(2)/=0) STOP "Cannot allocate memory"

allocate (SRC(IBIR,J2),STAT=status(1)); allocate (CL(IBIR,J2),STAT=status(2)); allocate (CB(IBIR,J2),STAT=status(3)); allocate (CC(IBIR,J2),STAT=status(4))

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot allocate memory"

allocate (OCC(IBIR,J2),STAT=status(1))
if (status(1)/=0 ) STOP "Cannot allocate memory"

allocate (CR(IBIR,J2),STAT=status(1)); allocate (CT(IBIR,J2),STAT=status(2))

if (status(1)/=0 .or. status(2)/=0) STOP "Cannot allocate memory"

allocate (IE(NB0),STAT=status(1))

if (status(1)/=0) STOP "Cannot allocate memory"

RINV=0.;RINV1=0.;DUM0=0.;DUM1=0.;DUM2=0.;XX=0.;H=0.;X=0.;S=0.;AL=0.;AB=0.;AC=0.;OAC=0.;AR=0.;AT=0.;SRC=0.;CL=0.;CB=0.;CC=0.;OCC=0.;CR=0.;CT=0.;IE=0

END SUBROUTINE VAR_ALLOCATE_SEVP

!======================================================

SUBROUTINE VAR_ALLOCATE_PG
INTEGER :: status(4)
allocate(DHX(I1,J0),STAT=status(1))
allocate(DHY(I0,J1),STAT=status(2))
allocate(XP(I0,J0,K1),STAT=status(3))
allocate(YP(I0,J0),STAT=status(4))

if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot allocate memory"

DHX=0.
DHY=0.
XP=0.
YP=0.

END SUBROUTINE VAR_ALLOCATE_PG

!======================================================

SUBROUTINE VAR_ALLOCATE_ADD
INTEGER :: status(4)
!allocate(ADD(2:I1,2:J1,K2),STAT=status(1))
allocate(ADD(1:I0-2,1:J0-2,K2),STAT=status(1))
if (status(1)/=0) STOP "Cannot allocate memory"
ADD=0.
END SUBROUTINE VAR_ALLOCATE_ADD

SUBROUTINE DE_INI_VAR

INTEGER :: STATUS(4)

if (allocated(RHO)) then

  deallocate (RHO,STAT=status(1)); deallocate (P0,STAT=status(2)); deallocate (PX,STAT=status(3)); deallocate (PY,STAT=status(4))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot deallocate memory"

  deallocate (Mcrho,STAT=status(1)); deallocate (McTP,STAT=status(2)); deallocate (McS,STAT=status(3)); deallocate (McTPTP,STAT=status(4))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot deallocate memory"

  deallocate (McSTP,STAT=status(1))
  if (status(1)/=0) STOP "Cannot deallocate memory"


  deallocate (U1,STAT=status(1)); deallocate (U2,STAT=status(2)); deallocate (V1,STAT=status(3)); deallocate (V2,STAT=status(4))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot deallocate memory"

  deallocate (S1,STAT=status(1)); deallocate (S2,STAT=status(2)); deallocate (T1,STAT=status(3)); deallocate (T2,STAT=status(4))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot deallocate memory"

  deallocate (P,STAT=status(1)); deallocate (ULF,STAT=status(2)); deallocate (VLF,STAT=status(3)); deallocate (SLF,STAT=status(4))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot deallocate memory"

  deallocate (TLF,STAT=status(1)); deallocate (EV,STAT=status(2)); deallocate (HV,STAT=status(3))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0) STOP "Cannot deallocate memory"

  deallocate (F,STAT=status(1)); deallocate (TANPHI,STAT=status(2)); deallocate (SUMIN,STAT=status(3))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0) STOP "Cannot deallocate memory"

end if

#ifdef FLAG_TRACER
if (allocated(C1)) then
  deallocate (C1,STAT=status(1)); deallocate (C2,STAT=status(2)); deallocate (CLF,STAT=status(3))
  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0) STOP "Cannot deallocate memory"
endif
#endif

  deallocate (IFLT,STAT=status(1))
  if (status(1)/=0) STOP "Cannot deallocate memory"

if (allocated(U)) then

  deallocate (U,STAT=status(1)); deallocate (V,STAT=status(2)); deallocate (W,STAT=status(3))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0) STOP "Cannot deallocate memory"

end if

if (allocated(Z)) then

  deallocate (Z,STAT=status(1)); deallocate (ODZ,STAT=status(2)); deallocate (ODZW,STAT=status(3))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0) STOP "Cannot deallocate memory"

end if

if (allocated(KB)) then

  deallocate (KB,STAT=status(1)); deallocate (IU0,STAT=status(2)); deallocate (IV0,STAT=status(3))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0) STOP "Cannot deallocate memory"

  deallocate (IN,STAT=status(1)); deallocate (IU,STAT=status(2)); deallocate (IV,STAT=status(3)); deallocate (IW,STAT=status(4))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot deallocate memory"

end if

if (allocated(A)) then

  deallocate (A,STAT=status(1)); deallocate (XDEG,STAT=status(2)); deallocate (YV,STAT=status(3)); deallocate (YVDEG,STAT=status(4))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot deallocate memory"

  deallocate (YDEG,STAT=status(1)); deallocate (CS,STAT=status(2)); deallocate (OCS,STAT=status(3)); deallocate (DX,STAT=status(4))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot deallocate memory"

  deallocate (ODX,STAT=status(1)); deallocate (DY,STAT=status(2)); deallocate (ODY,STAT=status(3))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0) STOP "Cannot deallocate memory"

  deallocate (CSV,STAT=status(1)); deallocate (OCSV,STAT=status(2)); deallocate (DXV,STAT=status(3))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0) STOP "Cannot deallocate memory"

  deallocate (ODXV,STAT=status(1)); deallocate (DYV,STAT=status(2)); deallocate (ODYV,STAT=status(3))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0) STOP "Cannot deallocate memory"

end if

if (allocated(Y)) then

  deallocate (Y,STAT=status(1))

  if (status(1)/=0) STOP "Cannot deallocate memory"

end if

if (allocated(PBAR)) then

  deallocate (PBAR,STAT=status(1)); deallocate (PVAR,STAT=status(2)); deallocate (XBAR,STAT=status(3)); deallocate (UCLI,STAT=status(4))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot deallocate memory"

  deallocate (VCLI,STAT=status(1)); deallocate (RMSV,STAT=status(2)); deallocate (SBAR,STAT=status(3)); deallocate (TBAR,STAT=status(4))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot deallocate memory"

end if

if (allocated(VGLO)) then

  deallocate (VGLO,STAT=status(1)); deallocate (UAVG,STAT=status(2)); deallocate (VAVG,STAT=status(3)); deallocate (SAVG,STAT=status(4))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot deallocate memory"

  deallocate (TAVG,STAT=status(1))

  if (status(1)/=0) STOP "Cannot deallocate memory"

end if

if (allocated(SCR)) then

  deallocate (SCR,STAT=status(1))

  if (status(1)/=0) STOP "Cannot deallocate memory"

end if

if (allocated(DMX)) then

  deallocate (DMX,STAT=status(1)); deallocate (DMY,STAT=status(2))

  if (status(1)/=0 .or. status(2)/=0) STOP "Cannot deallocate memory"

  deallocate (VBK,STAT=status(1)); deallocate (HBK,STAT=status(2)); deallocate (INFX,STAT=status(3)); deallocate (INFY,STAT=status(4))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot deallocate memory"

end if

if (allocated(SXY)) then

  deallocate (SXY,STAT=status(1)); deallocate (TXY,STAT=status(2)); deallocate (SXYCLI,STAT=status(3)); deallocate (TXYCLI,STAT=status(4))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot deallocate memory"

end if

if (allocated(KP)) then

  deallocate (KP,STAT=status(1))

  if (status(1)/=0) STOP "Cannot deallocate memory"

end if

if (allocated(TNUDGE)) then

  deallocate (TNUDGE,STAT=status(1)); deallocate (QAVG,STAT=status(2)); deallocate (WAVG,STAT=status(3))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0) STOP "Cannot deallocate memory"

  deallocate (SCLI,STAT=status(1)); deallocate (TCLI,STAT=status(2)); deallocate (QDAVG,STAT=status(3)); deallocate (SDAVG,STAT=status(4))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot deallocate memory"

end if

if (allocated(SSURF)) then

   deallocate (SSURF,STAT=status(3)); deallocate (TSURF,STAT=status(4))

  if (status(3)/=0 .or. status(4)/=0) STOP "Cannot deallocate memory"

  deallocate (SSSP,STAT=status(1)); deallocate (SNSP,STAT=status(2)); deallocate (TSSP,STAT=status(3)); deallocate (TNSP,STAT=status(4))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot deallocate memory"

end if

if (allocated(TAUX)) then

   deallocate (TAUX,STAT=status(3)); deallocate (TAUY,STAT=status(4))

  if (status(3)/=0 .or. status(4)/=0) STOP "Cannot deallocate memory"  

end if

if (allocated(RINV)) then

  deallocate (RINV,STAT=status(1)); deallocate (RINV1,STAT=status(2)); deallocate (DUM0,STAT=status(3))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0) STOP "Cannot deallocate memory"

  deallocate (DUM1,STAT=status(1)); deallocate (DUM2,STAT=status(2)); deallocate (XX,STAT=status(3)); deallocate (H,STAT=status(4))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot deallocate memory"

  deallocate (X,STAT=status(1)); deallocate (S,STAT=status(2)); deallocate (AL,STAT=status(3)); deallocate (AB,STAT=status(4))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot deallocate memory"

  deallocate (AC,STAT=status(1)); deallocate (AR,STAT=status(2)); deallocate (AT,STAT=status(3));deallocate (OAC,STAT=status(4))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot deallocate memory"

  deallocate (SRC,STAT=status(3)); deallocate (CL,STAT=status(4))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot deallocate memory"

  deallocate (CB,STAT=status(1)); deallocate (CC,STAT=status(2)); deallocate (CR,STAT=status(3)); deallocate (CT,STAT=status(4))

  if (status(1)/=0 .or. status(2)/=0 .or. status(3)/=0 .or. status(4)/=0) STOP "Cannot deallocate memory"

  deallocate (OCC,STAT=status(1))
  if (status(1)/=0 ) STOP "Cannot deallocate memory"

  deallocate (IE,STAT=status(1))

  if (status(1)/=0) STOP "Cannot deallocate memory"

end if

if (allocated(DEPTH)) then

  deallocate (DEPTH,STAT=status(1))

  if (status(1)/=0) STOP "Cannot deallocate memory"

end if

if (allocated(ADD)) then

  deallocate (ADD,STAT=status(1))

  if (status(1)/=0) STOP "Cannot deallocate memory"

end if

END SUBROUTINE DE_INI_VAR

END MODULE GRID_VAR_GLOBAL


!######################################################################


MODULE XYFS_GLOBAL
CHARACTER(6):: COORD_SYS
REAL::DYDX,DXMNUT
REAL::Y0DEG,X0DEG
REAL::DXDEG,DYDEG
INTEGER::DX0STD
REAL::R0=6.4E8
END MODULE XYFS_GLOBAL


MODULE ZCO_GLOBAL    ! LI
USE OCN_PARA_GLOBAL,ONLY:K0
INTEGER::ZOP
REAL::ZTOP   !  ONLY APPLIED TO LIENAR AND STRATCHED COORDINATE
REAL::C=0.01,D=0.03    !  ONLY APPLIED TO LINEAR EXPONENTIAL Z COORDINATE
!REAL,DIMENSION(K0)::Z_CUSTOM   !  ONLY APPLIED TO USER-DEFINED Z COORDINATE
END MODULE ZCO_GLOBAL


MODULE GLOBALSPLT
USE OCN_PARA_GLOBAL,ONLY:I0,J0,K1
INTEGER::NEXP
REAL::DTS
REAL,DIMENSION(:,:,:),ALLOCATABLE::SC,TC,UC,VC
PRIVATE :: ALLOCATE_STATUS

CONTAINS

SUBROUTINE VAR_ALLOCATE_GLOBALSPLT
integer :: istat
ALLOCATE(SC(I0,J0,K1),STAT=istat)
CALL ALLOCATE_STATUS(istat)
ALLOCATE(TC(I0,J0,K1),STAT=istat)
CALL ALLOCATE_STATUS(istat)
ALLOCATE(UC(I0,J0,K1),STAT=istat)
CALL ALLOCATE_STATUS(istat)
ALLOCATE(VC(I0,J0,K1),STAT=istat)
CALL ALLOCATE_STATUS(istat)
END SUBROUTINE VAR_ALLOCATE_GLOBALSPLT

SUBROUTINE ALLOCATE_STATUS(istat)
INTEGER :: istat
if (istat/=0) STOP "Cannot allocate memory"
END SUBROUTINE ALLOCATE_STATUS

END MODULE GLOBALSPLT

