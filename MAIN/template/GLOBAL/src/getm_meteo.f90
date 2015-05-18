#define DEBUG

module getm_meteo

  real, public	:: www,L,rho_air,qs,qa,ea,cd_heat,cd_mom
  real,public,parameter :: cpa=1008.,KELVIN=273.15,emiss=0.97,bolz=5.67e-8

  !intel fortran special numbers (NaN, Inf...)
  !INTEGER(4), PARAMETER :: FOR_K_FP_SNAN       = '00000000'X
  !INTEGER(4), PARAMETER :: FOR_K_FP_QNAN       = '00000001'X
  INTEGER(4), PARAMETER :: FOR_K_FP_POS_INF    = '00000002'X
  INTEGER(4), PARAMETER :: FOR_K_FP_NEG_INF    = '00000003'X
  !INTEGER(4), PARAMETER :: FOR_K_FP_POS_NORM   = '00000004'X
  !INTEGER(4), PARAMETER :: FOR_K_FP_NEG_NORM   = '00000005'X
  !INTEGER(4), PARAMETER :: FOR_K_FP_POS_DENORM = '00000006'X
  !INTEGER(4), PARAMETER :: FOR_K_FP_NEG_DENORM = '00000007'X
  !INTEGER(4), PARAMETER :: FOR_K_FP_POS_ZERO   = '00000008'X
  !INTEGER(4), PARAMETER :: FOR_K_FP_NEG_ZERO   = '00000009'X

contains

subroutine compute_fluxes(Tmjd,gridx,gridy,cc,u10,v10,airp,sst,t2,hum,shf,swr,tausx,tausy,evapoo)
  use date
  use TIMCOM_GENERAL,ONLY:heattype
  implicit none

#include "./physic.h"
  real, intent(in) :: Tmjd
  real, dimension(:), intent(in) :: gridx,gridy
  real, dimension(:,:), intent(in) :: u10,v10,airp,sst,t2,hum
  real, dimension(:,:), intent(out) :: shf,swr,tausx,tausy,evapoo
  real, dimension(:,:) :: cc

  integer :: i,j, y,m,d
  integer :: yearday
  real :: s,lat,lon,hh,windmult 
  logical :: isnan,hasproblem


  !coefficient to artificially increase obtained momentum stress
  windmult=1.0

  if (any(cc.gt.1.0)) then
    cc=cc/100.0
  end if

  call gregd(Tmjd,y,m,d,s)
  yearday = Tmjd-mjd(y,1,1,0.)
  hh = s/3600

  do j=1,size(cc,2)   ! J2
      lat = gridy(j)
    do i=1,size(cc,1) ! I2
      lon = gridx(i)
      !we used to substract 360° from model longitude; now this is done with the "rotateforcing" init variable
      !if (lon.gt.180) lon=lon-360

      !short wave
      if (heattype.ne.3.and.heattype.ne.-1) &
        swr(i,j) = short_wave_radiation(yearday,hh,lon,lat,cc(i,j))

      !long wave, evaporation, momentum
      call exchange_coefficients(windmult,u10(i,j),v10(i,j),t2(i,j),airp(i,j),sst(i,j),hum(i,j))
      call fluxes(u10(i,j),v10(i,j),t2(i,j),cc(i,j),sst(i,j),shf(i,j),tausx(i,j),tausy(i,j),evapoo(i,j) )
    end do
  end do

! dimensions at this point
!  tausx,tausy: N/m2
!  t2: deg C
!  evap: kg/m2/s
! convertion into model fluxes

! luc 2012, for TIMCOM I don't need this variable fles
!  fles = -1.8e-9 *  sqrt(1e3 * sqrt(tausx**2 + tausy**2))**3

!In TIMCOM, all these conversions are done later in the code.
!at this point, I just need fluxes in  W/m2 and evap in kg/m2/s
! so i don't do anything anymore here

!  tausx = tausx * UNSRO0
!  tausy = tausy * UNSRO0
!  evap  = -  evap * UNSRO0
! normalement evap est maintenant négatif si il y a evaporation 
! luc 2012, j'enleve le signe moins ci-dessus, je veux evap positif s'il y a evaporation
! en outre je ne convertis plus avec UNSRO0 donc je commente carrement toute la ligne

! 1/(rho_0 * heat_capacity) = 0.25e-6
!  shf = -UNSRO0/HeatCapacity * shf  
!  swr = -UNSRO0/HeatCapacity * swr

  end subroutine 

!$Id: fluxes.F90,v 1.1.1.1 2002/05/02 14:01:39 gotm Exp $
!#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Heat and momentum fluxes.
!
! !INTERFACE:
   subroutine fluxes(u10,v10,airt,cc,sst,hf,taux,tauy,evappoint)
   use TIMCOM_GENERAL,ONLY:heattype

! !DESCRIPTION:
!  The sum of the latent and sensible heat fluxes + longwave
!  back-radiation is calculated and returned in \emph{hf} [$W/m^2$]. Also the
!  sea surface stresses are calculated and returned in \emph{taux} and
!  \emph{tauy} [$N/m^2$], repsectively. The wind velocities are following the
!  meteorological convention (from where) and are in $m/s$. The
!  temperatures \emph{airt} and \emph{sst} can be in Kelvin or Celcius -
!  if they are $>$ 100 - Kelvin is assumed. \emph{cc} - the cloud cover -
!  is specified as fraction between 0 and 1.
!
! !SEE ALSO:
!  meteo.F90, exchange_coefficients.F90
!
! !USES:
!   use meteo, only: cpa,emiss,bolz,KELVIN
!   use meteo, only: w,L,rho_air,qs,qa,ea,cd_heat,cd_mom
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   real, intent(in)		:: u10,v10,airt,cc,sst
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   real, intent(out)	:: hf,taux,tauy,evappoint
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding and Hans Burchard
!
!  $Log: fluxes.F90,v $
!  Revision 1.1.1.1  2002/05/02 14:01:39  gotm
!  recovering after CVS crash
!
!  Revision 1.1  2001/07/26 14:35:18  bbh
!  initial import into CVS
!
!
! !DEFINED PARAMETERS:
   integer, parameter   :: clark=1      ! Clark et. al, 1974
   integer, parameter   :: hastenrath=2 ! Hastenrath and Lamb, 1978
!
! !LOCAL VARIABLES:
   real		:: tmp
   real		:: qe,qh,qb
   real		:: ta,tw,tw_k
   integer      :: back_radiation_method=clark
!
!EOP
!-----------------------------------------------------------------------
!BOC
!   if (sst .lt. 100.) then
      tw  = sst
      tw_k= sst+KELVIN
!   else
!      tw  = sst-KELVIN
!      tw_k= sst
!   end if

!   if (airt .gt. 100.) then
!      ta  = airt - KELVIN
!   else
      ta = airt
!   end if

   evappoint=cd_heat*rho_air*www*(qs-qa)          ! evaporation in kg/m²s      
   qe=L*evappoint			          ! latent
   qh=cd_heat*cpa*rho_air*www*(tw-ta)             ! sensible heat flux, positive going up

   select case(back_radiation_method)		  ! back radiation
      case(clark)
         qb=emiss*bolz*(tw_k**4)*(0.39-0.05*sqrt(ea/100.0))	&
            +4.0*emiss*bolz*(tw_k**3)*(tw-ta)
      case(hastenrath) ! qa in g(water)/kg(wet air)
         qb=emiss*bolz*(tw_k**4)*(0.39-0.056*sqrt(1000*qa))		&
           +4.0*emiss*bolz*(tw_k**3)*(tw-ta)
      case default
   end select
   if (heattype.ne.3.and.heattype.ne.-1) &
     qb=(1.0-0.8*cc*cc)*qb
     !when heattype=3, we don't want to remove the back-scattering from the outgoing longwave
     !it will be read from disk, and removed later (in the fs routine)

   !total outgoing heat:
   hf = qe+qh+qb ! for timcom, outgoing heat is positive going up

   tmp   = cd_mom*rho_air*www
   taux  = tmp*u10
   tauy  = tmp*v10

   return
   end subroutine fluxes
!EOC

!-----------------------------------------------------------------------
!Copyright (C) 2001 - Karsten Bolding & Hans Burchard
!-----------------------------------------------------------------------

!$Id: short_wave_radiation.F90,v 1.1.1.1 2002/05/02 14:01:39 gotm Exp $
!#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Short wave radiation.
!
! !INTERFACE:
   real function short_wave_radiation(yday,hour,lat,lon,cc)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Short wave radiation is calculated based on the following input
!  parameters - year day, hour of day, latitude and longitude and cloud cover.
!  The albedo monthly values are from Payne (1972) as means of the values
!  at 40N and 30N for the Atlantic Ocean ( hence the same latitudinal
!  band of the Mediterranean Sea ) :
!  The radiation is returned as $W/m^2$.
!
! !SEE ALSO:
!  meteo.F90
!
! !USES:
!
! !INPUT PARAMETERS:
   integer, intent(in)	:: yday
   real, intent(inout)	:: hour
   real, intent(in)	:: lat,lon,cc
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding and Hans Burchard
!
!  $Log: short_wave_radiation.F90,v $
!  Revision 1.1.1.1  2002/05/02 14:01:39  gotm
!  recovering after CVS crash
!
!  Revision 1.1  2001/07/26 14:35:18  bbh
!  initial import into CVS
!
!
! !DEFINED PARAMETERS:
   real, parameter  :: pi=3.1415926535897932384626433832795029
   real, parameter  :: deg2rad=pi/180.,rad2deg=180./pi

   real, parameter	:: solar=1350.
   real, parameter	:: eclips=23.439*deg2rad
   real, parameter	:: tau=0.7
   real, parameter	:: aozone=0.09

   real, parameter	:: yrdays(2)= (/365.,366./)

! element alb1(21) is not used, but should exist
! modif Mahdia alb1(21) devient alb1(20)

   real, parameter	:: alb1(20) = 	&
                    (/.719,.656,.603,.480,.385,.300,.250,.193,.164, &
                      .131,.103,.084,.071,.061,.054,.039,.036,.032,.031,.030 /)

   real, parameter	:: za(20) = (/90.,88.,86.,84.,82.,80.,78.,76.,74.,70., &
                                      66.,62.,58.,54.,50.,40.,30.,20.,10.,0.0 /)
   real	:: dza(19)
   data            dza/8*2.0, 6*4.0, 5*10.0/

! !LOCAL VARIABLES:
   real	:: alat,alon,eqnx
   real	:: th0,th02,th03,sundec
   real	:: thsun,coszen,zen,dzen,sunbet
   real	:: qatten,qzer,qdir,qdiff,qtot,qshort
   real	:: albedo
   integer	:: jab
   logical :: isnanf, isinff
!
! !TO DO:
!
! !BUGS:
!  yrdays should be set depending on leap year or not.
!
!EOP
!-----------------------------------------------------------------------
!BOC
   th0 = 2.*pi*yday/yrdays(1)
   th02 = 2.*th0
   th03 = 3.*th0

!  sun declination :
   sundec = 0.006918 - 0.399912*cos(th0) + 0.070257*sin(th0)	&
           - 0.006758*cos(th02) + 0.000907*sin(th02)		&
           - 0.002697*cos(th03) + 0.001480*sin(th03)

   alon = deg2rad*lon
   alat = deg2rad*lat

!  sun hour angle :
   thsun = (hour-12.)*15.*deg2rad + alon

!  cosine of the solar zenith angle :
   coszen =sin(alat)*sin(sundec)+cos(alat)*cos(sundec)*cos(thsun)

   if (coszen .le. 0.0) then
      coszen = 0.0
      qatten = 0.0
   else
      qatten = tau**(1./coszen)
   end if
   qzer  = coszen * solar
   qdir  = qzer * qatten
   qdiff = ((1.-aozone)*qzer - qdir) * 0.5
   qtot  =  qdir + qdiff

   eqnx = (yday-81.)/yrdays(1)*2.*pi

!  sin of the solar noon altitude in radians :
   sunbet=sin(alat)*sin(eclips*sin(eqnx))+cos(alat)*cos(eclips*sin(eqnx))
!  necessary in simple precision 
!  alex
   sunbet = min(max(sunbet,-1.),1.)
   sunbet = asin(sunbet)*rad2deg


!  calculates the albedo as a function of the solar zenith angle :
!  (after Payne jas 1972)
!  solar zenith angle in degrees :

! necessary in simple precision 
! alex
   if (coszen.ge.1.) then
     zen = 0.
   else
     zen=(180./pi)*acos(coszen)
   end if

   if(zen .ge. 74.)then
      jab=.5*(90.-zen)+1.
   else if (zen .ge. 50.) then
      jab=.23*(74.-zen)+9.
   else
      jab=.10*(50.-zen)+15.
   endif

   if (jab.lt.1.or.jab+1.gt.size(alb1)) then
     !write(99,*) 'getm_meteo: Warning jab ',jab,coszen
     jab = min(size(alb1)-1,max(1,jab))
     !write(99,*) 'getm_meteo: changed jab ',jab
   end if

   dzen=(za(jab)-zen)/dza(jab)
   albedo=alb1(jab)+dzen*(alb1(jab+1)-alb1(jab))

!  radiation as from Reed(1977), Simpson and Paulson(1979)
!  calculates SHORT WAVE FLUX ( watt/m*m )
!  Rosati,Miyakoda 1988 ; eq. 3.8

!HB I do not know where this "if(cc .lt. 0.3) then" comes from. It is not
!HB included in Rosati,Miyakoda 1988 ; eq. 3.8 ...
!HB   if(cc .lt. 0.3) then
!HB      short_wave_radiation  = qtot
!HB   else
      short_wave_radiation  = qtot*(1-0.62*cc + .0019*sunbet)*(1.-albedo)
!HB   endif


   return
   end function short_wave_radiation
!EOC


!-----------------------------------------------------------------------
!Copyright (C) 2001 - Karsten Bolding & Hans Burchard
!-----------------------------------------------------------------------

! !ROUTINE: Air/sea exchange coefficients
!
! !INTERFACE:
   subroutine exchange_coefficients(windmult,u10,v10,airt,airp,sst,hum)
!
! !DESCRIPTION:
!  Various variables for calculating meteorological forcing is calculated
!  here. The input variables are given in SI units ($m/s$, $^oC | Kelvin$,
!  and $Pa$). The following formulaes are used:
!  \begin{itemize}
!     \item Latent heat: $L = (2.5-0.00234T)10^6$
!     \item Specific vapor pressure: $e_s = a_0 + a_1 T^1 + a_2 T^2 + a_3 T^3
!                                    + a_4 T^4 + a_5 T^5 + a_6 T^6 + a_7 T^7$
!     \item Specific humidity: $q_s = const06 e_s / (airp-0.377 e_s) $
!     \item Absolute humidity: $q_a = 0.01  rh q_s $
!     \item Absolute vapor pressure: $q_s = q_a airp / (airp + 0.377 q_a) $
!     \item Virtual temperature: $T_v = T_k ( 1 + q_a/const06 )/( 1 + q_a)$
!     \item Air density: $\rho_a = airp/(287.05 T_v)$
!     \item Air stability: $S_0 = 0.25(t_w - t_a)/W and s = S_0 |S_0|
!                                                         / (|S_0| +0.01)$
!  \end{itemize}
!  The variables $L, e_s, q_s, e_a$ and $q_a$ are later used for calculating
!  the latent and sensible heat fluxes.
!  The exchange coefficients $C_{D_{mom}}$ and $C_{D_{heat}}$ are functions of
!  emperical parameters - which a actual wind speed dependend - and the
!  air stability.
!
! !SEE ALSO:
!  meteo.F90, fluxes.F90
!
! !USES:
!   use meteo, only: cpa,KELVIN
!   use meteo, only: L,rho_air,w,qs,qa,cd_heat,cd_mom
  use TIMCOM_GENERAL,ONLY:humtype
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   real, intent(in)	:: windmult,u10,v10,airt,airp,sst,hum
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  $Log: exchange_coefficients.F90,v $
!  Revision 1.1.1.1  2002/05/02 14:01:38  gotm
!  recovering after CVS crash
!
!  Revision 1.1  2001/07/26 14:35:18  bbh
!  initial import into CVS
!
!
! !DEFINED PARAMETERS:
   real, parameter	:: a1=6.107799961
   real, parameter	:: a2=4.436518521e-1
   real, parameter	:: a3=1.428945805e-2
   real, parameter	:: a4=2.650648471e-4
   real, parameter	:: a5=3.031240396e-6
   real, parameter	:: a6=2.034080948e-8
   real, parameter	:: a7=6.136820929e-11
   real, parameter	:: const06=0.62198
   real, parameter	:: eps=1.0e-12
!
! !LOCAL VARIABLES:
   real		:: tvirt
   real(kind=8) :: s,s0
   real		:: ae_h,be_h,ce_h,pe_h
   real		:: ae_m,be_m,ce_m,pe_m
   real		:: x,es,ea
   real		:: cee_heat,cee_mom
   real		:: ta,ta_k,tw,tw_k

   real		:: twet,rh
   real		:: d_tmp
!
!  !TO DO:
!
!  !BUGS:
!   How to calculate absolute humidity based on the input variables.
!
!EOP
!-----------------------------------------------------------------------
!BOC

!   if (sst .lt. 100.) then
      tw  = sst
      tw_k= sst+KELVIN
!   else
!      tw  = sst-KELVIN
!      tw_k= sst
!   end if

!   if (airt .lt. 100.) then
      ta   = airt
      ta_k = airt+KELVIN
!   else
!      ta   = airt-KELVIN
!      ta_k = airt
!   end if

   www = sqrt(u10*u10+v10*v10)
   L = (2.5-0.00234*tw)*1.e6
   es = a1 +tw*(a2+tw*(a3+tw*(a4+tw*(a5+tw*(a6+tw*a7)))))
   es = es * 100.0 ! Conversion millibar --> Pascal
   qs = const06*es/(airp-0.377*es)

   if (humtype.ne.2) then
      d_tmp = hum
!     Piece of code taken from HAMSOM for calculating relative
!     humidity from dew point temperature and dry air temperature.
!     ARTHUR changé ta_k en tw_k pour le calcul de es selon une correction observée dans GETM
      ea  = 611.21*exp((18.729 - (min(d_tmp,300.)-273.15)/227.3)*       &
              (min(d_tmp,300.)-273.15)/(max(d_tmp,200.)-273.15+257.87))
      es  = 611.21*exp((18.729 - (min(tw_k,300.)-273.15)/227.3)*       &
              (min(tw_k,300.)-273.15)/(max(tw_k,200.)-273.15+257.87))
      !correction for seawater following Kraus 1972
      es = 0.98 * es
      rh = ea/es * 100.

      qa = CONST06*ea/(airp-0.377*ea);
   else
      qa=hum !hum actually comes from DewTemperature2m which itself actually containted specific humidity @ 10m
      ea=qa*airp/(const06+0.377*qa) ! I inverted the formula from above linking qa and ea
      es  = 611.21*exp((18.729 - (min(tw_k,300.)-273.15)/227.3)*       &
              (min(tw_k,300.)-273.15)/(max(tw_k,200.)-273.15+257.87))
      rh = ea/es * 100.      
   end if

   tvirt = ta_k*(1+qa/const06)/(1+qa)
   rho_air = airp/(287.05*Tvirt)

!  Stability
   s0=0.25*(tw-ta)/(www+1.0e-10)**2
   s=s0/(abs(s0)+0.01)*abs(s0)

!  Transfer coefficient for heat and momentum
   if (www.lt.2.2) then
      ae_h=0.0;   be_h=1.23;   ce_h=0.0;      pe_h=-0.16;
      ae_m=0.0;   be_m=1.185;  ce_m=0.0;      pe_m=-0.157;
   else if (www.lt.5.0) then
      ae_h=0.969; be_h=0.0521; ce_h=0.0;      pe_h=1.0;
      ae_m=0.927; be_m=0.0546; ce_m=0.0;      pe_m=1.0;
   else if (www.lt.8.0) then
      ae_h=1.18;  be_h=0.01;   ce_h=0.0;      pe_h=1.0;
      ae_m=1.15;  be_m=0.01;   ce_m=0.0;      pe_m=1.0;
   else if (www.lt.25.0) then
      ae_h=1.196; be_h=0.008;  ce_h=-0.0004;  pe_h=1.0
      ae_m=1.17;  be_m=0.0075; ce_m=-0.00044; pe_m=1.0;
   else
      ae_h=1.68;  be_h=-0.016; ce_h=0;        pe_h=1;
      ae_m=1.652; be_m=-0.017; ce_m=0.0;      pe_m=1.0;
   end if

   cee_heat=(ae_h+be_h*exp(pe_h*log(www+eps))+ce_h*(www-8.0)**2)*1.0e-3
   cee_mom =(ae_m+be_m*exp(pe_m*log(www+eps)))*1.0e-3

   if(s .lt. 0.) then
      x = 0.1+0.03*s+0.9*exp(4.8*s)
      cd_heat=x*cee_heat
      cd_mom =x*cee_mom
   else
      cd_heat=cee_heat*(1.0+0.63*sqrt(s))
      cd_mom =cee_mom *(1.0+0.47*sqrt(s))
   end if

   cd_heat=cd_heat*windmult
   return
   end subroutine exchange_coefficients
!EOC

!-----------------------------------------------------------------------
!Copyright (C) 2001 - Karsten Bolding & Hans Burchard
!-----------------------------------------------------------------------



end module
