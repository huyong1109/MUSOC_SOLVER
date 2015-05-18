MODULE PHYSICS_GLOBAL
USE OCN_PARA_GLOBAL
!use B_GLOBAL
USE GRID_VAR_GLOBAL, ONLY: RHO,Mcrho,McTP,McS,McTPTP,McSTP,P0,PX,PY,U1,U2,V1,V2,S1,S2,T1,T2,P,ULF,VLF,SLF,TLF,EV,HV,F,TANPHI,SUMIN
!use CGRID_GLOBAL
USE GRID_VAR_GLOBAL, ONLY: U,V,W
!use ZFS_GLOBAL
USE GRID_VAR_GLOBAL, ONLY: Z,ODZ,ODZW
!use BATHY_GLOBAL
USE GRID_VAR_GLOBAL, ONLY: KB,IU0,IV0,IN,IU,IV,IW
use CONTROL_GLOBAL
use SCA_GLOBAL
!use WINDMX_GLOBAL_MAIN
USE GRID_VAR_GLOBAL, ONLY: VBK,HBK,INFX,INFY
USE GRID_VAR_GLOBAL, ONLY: ADD

INTEGER,PARAMETER :: TURB_PP82=0,TURB_P81=1

CONTAINS

SUBROUTINE PP82
! ----------------------------------------------------------------------
!calculate gradient Ri based vertical mixing coefficents
! so called "eddy viscosity and diffusivity"
! as per Pacanowski and Philander(1982)
! APPROACH:
! Set background vertical diffusivities to DMZ0 (O(0.1) cm-cm/s) and
! add Richardson number based vertical diffusivity plus numerical
! contribution according to vertical cell Reynolds number.
! Vertical cell Re is O(100), except during winter cooling conditions
! and in wind-blown surface mixed layer, as internal waves, which
! dominate W below the SML in the model results during summer, do not
! mix T or S in nature.
! NOTE: one cannot depend on diffusive closure by itself if one wants to
!       model possible contra-diffusive effects
     
! In loop 750 below, EV,HV units are cm-cm/s
! but EV,HV normalization by DZ is done 
      DO 750 K=1,K2
      L=K+1
      HBK0=HBK(K)
      TMPW=ORZMX/ODZW(L)
      DO 750 J=2,J1
      DO 750 I=2,I1
! TEMP must have units of cm-cm/s
      TEMP=TMPW*ABS(W(I,J,L))
! rho is double precision
!     RI=MAX(-0.5D0,980.*(RHO(I,J,L)-RHO(I,J,K))*ODZW(L)/(ODZW(L)**2*
      RI=MAX(-0.9D0,980.*(RHO(I,J,L)-RHO(I,J,K))*ODZW(L)/(ODZW(L)**2* &
      (0.001D0+(U2(I,J,L)-U2(I,J,K))**2+(V2(I,J,L)-V2(I,J,K))**2)))
      TMP=1./(1.+RI)
      TEMP=TMP*TEMP
! we add ODZW factor & apply explicit stability limit 
      EVISC=MIN(20.*TMP**2,100.)
      HV(I-1,J-1,K)=(EVISC*TMP+HBK0+TEMP)*ODZW(L)
! ADD includes VBK0
 750  EV(I-1,J-1,K)=(EVISC+TEMP+ADD(I-1,J-1,K))*ODZW(L)
 
END SUBROUTINE PP82

END MODULE PHYSICS_GLOBAL
