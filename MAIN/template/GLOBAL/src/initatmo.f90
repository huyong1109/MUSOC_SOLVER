SUBROUTINE INITATMO()
! ----------------------------------------------------------------------
use TIMCOM_GENERAL 
use OCN_PARA
use INPUT

use interfext, interfinit => init, interfgetfield => getfield
use initfile 

!============================================
! Load the data and interpolate on model grid
!============================================
write(*,*) "Loading atmospheric data specified in: ",trim(atmosfile)
bulk=.false.
humtype=0
windtype=-1
cloudtype=0
if (presentInitValue(atmosfile,'wind')) then
  call getInitValue(atmosfile,'wind',windtype)
endif
if (windtype.eq.0) then
  write(*,*) "You asked to specify wind speed; momentum will be calculated"
  bulk=.true.
elseif (windtype.eq.1) then
  write(*,*) "You asked to specify surface momentum stress immediately"
  call interfinit(MomentfluxX,atmosfile,'tauX10m.')
  call interfinit(MomentfluxY,atmosfile,'tauY10m.')
else
  write(*,*) "You asked for (traditional) timcom wind momentum forcings"
end if

heattype=-1
if (presentInitValue(atmosfile,'heat')) &
  call getInitValue(atmosfile,'heat',heattype)
if (presentInitValue(atmosfile,'humtype')) &
  call getInitValue(atmosfile,'humtype',humtype)
if (heattype.eq.0.and.humtype.eq.0) then
  bulk=.true.
  cloudtype=1
  write(*,*) "You asked to specify T2M, DT2M and CC; heat fluxes will be calculated"
  call interfinit(DewTemperature2m,atmosfile,'DewTemp2m.')
  !call interfinit(AtmPressure,atmosfile,'AtmPressure.')
  !call interfinit(PME,atmosfile,'PME.')
elseif (heattype.eq.0.and.humtype.eq.1) then !before, this resumed to condition heattype.eq.1
  bulk=.true.
  cloudtype=1
  write(*,*) "You asked to specify T2M, RelHum and CC; heat fluxes will be calculated"
  call interfinit(RelHum2m,atmosfile,'RelHum2m.')
  !call interfinit(AtmPressure,atmosfile,'AtmPressure.')
  !call interfinit(PME,atmosfile,'PME.')
elseif (heattype.eq.2) then
  write(*,*) "You want to use prescribed heat fluxes, they will be loaded from disk now"
  call interfinit(shortwavedown,atmosfile,'SWRdown.')
  call interfinit(longwaveup,atmosfile,'LWRup.')
  call interfinit(sensibleup,atmosfile,'SSenHF.')
  call interfinit(latentup,atmosfile,'SLatHF.')
  if (size(shortwavedown%day).ne.size(longwaveup%day).or.size(sensibleup%day).ne.size(latentup%day)) &
        write(*,*) "PLEASE USE HEAT FLUX INPUTS FOR shortwave_down, longwave_up, sensible_up, latent_up WITH IDENTICAL TIME STEPS"
  Qup%initfname="sum of 3 other field2D variables"
  Qup%key="Qup"
  Qup%g2D=longwaveup%g2D
  Qup%initialised=longwaveup%initialised
  allocate(Qup%array(size(longwaveup%array,1),size(longwaveup%array,2),size(longwaveup%array,3)))
  Qup%array=(longwaveup%array)+(sensibleup%array)+(latentup%array)
  allocate(Qup%day(size(longwaveup%day)))
  Qup%day=longwaveup%day
  Qup%scale=1
  Qup%index=1
  deallocate(longwaveup%array,latentup%array,sensibleup%array)
  deallocate(longwaveup%day,latentup%day,sensibleup%day)
  !write(*,*) shortwavedown%array(1,1,1),longwaveup%array(1,1,1),sensibleup%array(1,1,1),latentup%array(1,1,1)
elseif (heattype.eq.3) then
  bulk=.true.
  !cloudtype stays at 0
  humtype=2
  write(*,*) "Shortwave- (solar down) and longwave back- (downward) radiation is loaded from disk now; longwave radiation, sensible and latent heat fluxes will be calculated later (bulk)"
  call interfinit(shortwavedown,atmosfile,'SWRdown.')
  call interfinit(longwaveup,atmosfile,'LWRdown.')
  Qup%initfname="sum of 3 upward radiation fields"
  Qup%key="Qup"
  Qup%g2D=longwaveup%g2D
  Qup%initialised=longwaveup%initialised
  allocate(Qup%array(size(longwaveup%array,1),size(longwaveup%array,2),size(longwaveup%array,3)))
  Qup%array=-1*(longwaveup%array)
  allocate(Qup%day(size(longwaveup%day)))
  Qup%day=longwaveup%day
  Qup%scale=1
  Qup%index=1
  deallocate(longwaveup%array)
  deallocate(longwaveup%day)
  call interfinit(DewTemperature2m,atmosfile,'SpecificHum2m.')
  DewTemperature2m%initfname="Specific Humidity"
else
  write(*,*) "You asked for (traditional) timcom heat fluxes"
end if  

salttype=-1
if (presentInitValue(atmosfile,'salt')) then
     call getInitValue(atmosfile,'salt',salttype)
endif
if (salttype.eq.-1) then
  write(*,*) "You did NOT ask to compute evaporation, I will do nothing about evaporation"
elseif (salttype.eq.1) then
  bulk=.true.
  write(*,*) "You asked to compute evaporation interactively, surface salinity will be modified accordingly"
elseif (salttype.eq.0) then
  call interfinit(EvapField,atmosfile,'EVAP.')
  write(*,*) "You asked to load pre-computated evaporation, surface salinity will be modified accordingly"
end if


raintype=-1
if (presentInitValue(atmosfile,'rain')) then
     call getInitValue(atmosfile,'rain',raintype)
endif
if (raintype.eq.1) then
  write(*,*) "You asked to use rain forcing, they will be loaded from disk now"
  call interfinit(Precip,atmosfile,'PRECIP.')
else
  write(*,*) "You did not ask for rain forcings..."
end if

if (bulk) then
  call interfinit(WindU10m,atmosfile,'WindU10m.')
  call interfinit(WindV10m,atmosfile,'WindV10m.')
  call interfinit(AirTemperature2m,atmosfile,'AirTemperature2m.')
  if (heattype.ne.3.and.heattype.ne.-1) call interfinit(CloudCoverage,atmosfile,'CloudCoverage.')
  call interfinit(PresMSL,atmosfile,'AtmPressure.')
  if (heattype.eq.-1.or.heattype.eq.2) then !because we need some kind of humidity to use bulk formulae (even if we don't use those for heat fluxes)
    if (humtype.eq.1) &
      call interfinit(RelHum2m,atmosfile,'RelHum2m.')
    if (humtype.eq.2) then
      call interfinit(DewTemperature2m,atmosfile,'SpecificHum2m.')
      DewTemperature2m%initfname="Specific Humidity"
    end if
  end if
  if (humtype.eq.1) then ! before this condition was: (heattype.eq.1.or.heattype.eq.2)
    write(*,*) "Converting RelHum to DewPointTemperature"
    !!convert from RELHUM to DT2M, a good formula:
    !param_a=17.271
    !param_b=237.7
    !gamma=(param_a*T2M)/(param_b+T2M) + log(RELHUM/100)
    !DT2M=(param_b*gamma)/(param_a-gamma)
    !!convert from RELHUM to DT2M, a lousy formula:
    !DT2M=T2M-(100-RH)/5
    DewTemperature2m%initfname="DT2M obtained from RelHum"
    DewTemperature2m%key="DewTemperature at 2m"
    DewTemperature2m%g2D=RelHum2m%g2D
    DewTemperature2m%initialised=RelHum2m%initialised
    allocate(DewTemperature2m%array(size(RelHum2m%array,1),size(RelHum2m%array,2),size(RelHum2m%array,3)))
    DewTemperature2m%array=AirTemperature2m%array-((100-RelHum2m%array)/5.0)+273.15
    deallocate(RelHum2m%array)
    allocate(DewTemperature2m%day(size(RelHum2m%day)))
    DewTemperature2m%day=RelHum2m%day
    deallocate(RelHum2m%day)
    DewTemperature2m%scale=1
    DewTemperature2m%index=1
  end if
end if
write(*,*) "InitAtmo determined that cloudtype=",cloudtype
END SUBROUTINE INITATMO
