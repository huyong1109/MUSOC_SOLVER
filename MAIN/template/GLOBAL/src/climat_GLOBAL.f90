MODULE CLIMAT_GLOBAL
use OCN_PARA_GLOBAL,ONLY:I0,J0,I2,J2,K1,K2
IMPLICIT NONE
! Derived model climate source arrays
! Annual cycle sources in all sponge layer locations
!REAL,DIMENSION(I2,J2)::TNUDGE
!REAL,DIMENSION(I2,J2,12)::QAVG,WAVG
!REAL,DIMENSION(I0,J0,K1)::SCLI,TCLI,QDAVG,SDAVG
INTEGER,DIMENSION(12)::NSOMBO
! climate data for winter
! northern and southern sponge layer Levitus Climatology
! COMMON BLOCK/KLIMAT_GLOBAL/
REAL::FNEW,FOLD
INTEGER::NEW,NLD
! COMMON BLOCK/BCS_GLOBAL/
REAL,DIMENSION(12)::TSURFM,SSURFM
!REAL,DIMENSION(I2,J2,12)::SSURF,TSURF
!REAL,DIMENSION(I2,10,K2,12)::SSSP,SNSP,TSSP,TNSP
END MODULE CLIMAT_GLOBAL


!######################################################################


MODULE SURF_GLOBAL
use OCN_PARA_GLOBAL,ONLY:I0,J0,I2,J2,K1, localQsolar, localQup
USE GRID_VAR_GLOBAL,ONLY:QDOT,qdot2,evapo,transmitted,IN
USE TIMCOM_GENERAL,ONLY:heattype
use interfext, interfinit => init, interfgetfield => getfield
IMPLICIT NONE
INTEGER,PRIVATE::I,J
REAL,PRIVATE:: TMP,DTEMP

CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE QSURF(T,days,DT,ODZ)
! ----------------------------------------------------------------------
REAL,DIMENSION(I0,J0):: T
REAL,DIMENSION(*),INTENT(IN):: ODZ
REAL,INTENT(IN) :: DT
REAL, intent(in) :: days  ! current modified julian time
INTEGER:: I,J
logical:: space

	! The temperature array T  has a permiter "ghost zone".
! Thus, we set atmospheric heat exchange only in the interior zones.
!      QDOT=50. W/m2
! TMP is conversion factor from Watts per square meter
! to deg C per time step in top model layer.
!     TMP=DT*ODZ(1)/RHO/CP
! All units are cgs, so we use rho=1.
! CP for water is 2.5E4 ergs/gram/deg C.
TMP=DT*ODZ(1)/2.5E4
	      !add QDOT to T
	      DO J=1,J2
	      DO I=1,I2
	        T(I+1,J+1)=T(I+1,J+1) - (TMP*QDOT(i,j)*IN(i+1,j+1,1)) ! substract because heat is going OUT
	      END DO
	      END DO
	
END SUBROUTINE QSURF

!-----------------------------------------------------------------------------
SUBROUTINE QVOLUME(T,days,DT,ODZ) ! <<--- solar heat flux coming in the volume
!-----------------------------------------------------------------------------
REAL,DIMENSION(I0,J0,K1):: T
REAL,DIMENSION(*),INTENT(IN):: ODZ
REAL,INTENT(IN)::DT
REAL, intent(in) :: days  ! current modified julian time

INTEGER::I,J,K
logical::space

! transmitted*QDOT is added to T
DO K=1,k1
  TMP=DT*ODZ(k)/2.5E4 ! unit conversation
  DO I=1,I2
    do J=1,J2
      T(I+1,J+1,K)=T(I+1,J+1,K) + (TMP*QDOT2(i,j)*transmitted(i,j,k)*IN(i+1,j+1,k)) ! add because heat is going IN
    end do
  END DO
END DO

END SUBROUTINE QVOLUME

! ----------------------------------------------------------------------
SUBROUTINE SalSURF(S,days,DT,ODZ)
! ----------------------------------------------------------------------
REAL,DIMENSION(I0,J0)::	 S
REAL,DIMENSION(*),INTENT(IN)::	 ODZ
REAL,INTENT(IN) ::	 DT
REAL, intent(in) ::	 days  ! current modified julian time
INTEGER::	 I,J
logical::	 space

! Salinity= Salinity + [ tmp * sdot ]     (sdot==evapo)
! tmp converts precipitation and evaporation rate (kg/m2/s) to psu
! hence tmp is DT/DZ /waterdensity * salinity
! tmp =        dt*odz/1000 * psu

TMP=DT*ODZ(1)/1000*100 !the latter *100 because ODZ is in 1/cm rather than 1/m
!write(*,*) "mean salinity before=",sum(S)/size(S,1)/size(S,2)
DO J=1,J2
  DO I=1,I2
    S(I+1,J+1)=S(I+1,J+1)* (1 + (TMP*EVAPO(i,j)*IN(i+1,j+1,1))) ! EVAPO REPRESENTS EVAP-PRECIP. POSITIVE EVAPO --> SALINITY INCREASE
  END DO
END DO
!write(*,*) "mean salinity after =",sum(S)/size(S,1)/size(S,2)
	
END SUBROUTINE SalSURF
	
			
END MODULE SURF_GLOBAL

!######################################################################

 
MODULE WINDS_GLOBAL
use OCN_PARA_GLOBAL,ONLY:I0,J0,I1,J1,I2,J2
USE GRID_VAR_GLOBAL, ONLY:TAUX,TAUY !REAL,DIMENSION(I2,J2,12)::TAUX,TAUY
USE GRID_VAR_GLOBAL, ONLY:TAUX1,TAUY1
IMPLICIT NONE

CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE WIND(U,V,DAYS,DT,ODZ)
! ----------------------------------------------------------------------
use CLIMAT_GLOBAL,ONLY:FNEW,FOLD,NEW,NLD
USE TIMCOM_GENERAL,ONLY:windtype
use interfext, interfinit => init, interfgetfield => getfield
USE GRID_VAR_GLOBAL,ONLY:iu,iv
REAL,DIMENSION(I0,J0)::U,V
REAL,DIMENSION(*),INTENT(IN)::ODZ
REAL::DT,DAYS,TMP
INTEGER::I,J
logical::space	
! The velocity arrays U,V have a perimeter "ghost zone".
! Thus, we set wind forcing only in the interior zones.
! TAUX,TAUY are surface wind stress components.
! TAUX,TAUY units are force per unit area (i.e., energy per unit volume).
!     TMP=ODZ(1)/RHO
! All units are cgs, so we use RHO=1.
TMP=DT*ODZ(1)
if (windtype.eq.-1) then  !traditional timcom climatology
      DO J=2,J1
      DO I=2,I1
        U(I,J)=U(I,J)+TMP*(FOLD*TAUX(I-1,J-1,NLD)+FNEW*TAUX(I-1,J-1,NEW))
        V(I,J)=V(I,J)+TMP*(FOLD*TAUY(I-1,J-1,NLD)+FNEW*TAUY(I-1,J-1,NEW))
      END DO
      END DO
else ! when using windype=0 or 1, taux1 is already at the right time (thus only 1 layer)
      DO J=2,J1
        DO I=2,I1
          U(I,J)=U(I,J)+TMP*TAUX1(I-1,J-1)*iu(i,j,1)
          V(I,J)=V(I,J)+TMP*TAUY1(I-1,J-1)*iv(i,j,1)
        END DO
      END DO
endif
	
END SUBROUTINE WIND



END MODULE WINDS_GLOBAL
