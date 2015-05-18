SUBROUTINE INITFS_GLOBAL(I0_GLOBAL,J0_GLOBAL,K0_GLOBAL)
#include "namelist.in_GLOBAL"
! ****************
! INITIALIZATION *
! ****************
USE TIMCOM_GENERAL
USE OCN_PARA_GLOBAL
USE WINDID_GLOBAL
USE GRID_VAR_GLOBAL, ONLY: RINV,RINV1,DUM0,DUM1,DUM2,XX,H,X,S,AL,AB,AC,OAC,AR,AT,SRC,CL,CB,CC,CR,CT,IE
USE GRID_VAR_GLOBAL, ONLY: RHO,Mcrho,McTP,McS,McTPTP,McSTP,P0,PX,PY,U1,U2,V1,V2,S1,S2,T1,T2,P,ULF,VLF,SLF,TLF,EV,HV,F,TANPHI,SUMIN
USE GRID_VAR_GLOBAL, ONLY: U,V,W
USE GRID_VAR_GLOBAL, ONLY: PBAR,PVAR,XBAR,UCLI,VCLI,RMSV,SBAR,TBAR
USE GRID_VAR_GLOBAL, ONLY: Z,ODZ,ODZW
USE GRID_VAR_GLOBAL, ONLY: VGLO,UAVG,VAVG,SAVG,TAVG
USE GRID_VAR_GLOBAL, ONLY: A,XDEG,YV,YVDEG,YDEG,CS,OCS,DX,ODX,DY,ODY,CSV,OCSV,DXV,ODXV,DYV,ODYV
USE GRID_VAR_GLOBAL, ONLY: KB,IU0,IV0,IN,IU,IV,IW
USE CLIMAT_GLOBAL
USE GRID_VAR_GLOBAL, ONLY: TNUDGE,QAVG,WAVG,SCLI,TCLI,QDAVG,SDAVG,SSURF,TSURF,SSSP,SNSP,TSSP,TNSP
USE CONTROL_GLOBAL
USE SCA_GLOBAL
USE GRID_VAR_GLOBAL, ONLY: SCR
USE GRID_VAR_GLOBAL, ONLY:TAUX,TAUY
USE GRID_VAR_GLOBAL, ONLY: SXY,TXY,SXYCLI,TXYCLI
USE GRID_VAR_GLOBAL, ONLY: DMX,DMY,VBK,HBK,INFX,INFY
USE GRID_VAR_GLOBAL, ONLY: ADD,VAR_ALLOCATE_ADD
USE GRID_VAR_GLOBAL, ONLY: AROUT
USE GRID_VAR_GLOBAL, ONLY: nudge
USE GLOBALSPLT
USE INPUT_GLOBAL
USE INIT_VAR_GLOBAL
USE RESTART_GLOBAL
USE PMOVE_GLOBAL
USE REALWINDS_GLOBAL
#ifdef FLAG_TRACER
USE GRID_VAR_GLOBAL, ONLY: C1,C2,CLF
#endif
USE GRID_VAR_GLOBAL, ONLY:       &
    VAR_SPECIFY_OCN_PARA,       &
    VAR_ALLOCATE_B,             &
    VAR_ALLOCATE_CGRID,         &
    VAR_ALLOCATE_ZFS,           &
    VAR_ALLOCATE_BATHY,         &
    VAR_ALLOCATE_METRIC,        &
    VAR_ALLOCATE_TAVG_MAIN,     &
    VAR_ALLOCATE_GLO2NAB_MAIN,  &
    VAR_ALLOCATE_CLIMAT,        &
    VAR_ALLOCATE_SEVP,          &
    VAR_ALLOCATE_PG,            &
    VAR_ALLOCATE_SCRATCH_MAIN,  &
    VAR_ALLOCATE_WINDMX_MAIN,   &
    VAR_ALLOCATE_XYMEANS_MAIN,  &
    VAR_ALLOCATE_GLOBAL,         &
    VAR_ALLOCATE_ABSORPTION,    &
    calc_absorption_rate

use grids
use initfile
use ufileformat
USE DENSITY_EOS_GLOBAL

  integer :: imax,jmax,kmax,lucstat(6)

  REAL ZC, ZD, TMP3, TMP4 !NELSON, PLEASE PUT THIS DECLARATION IN A PROPER PLACE
  INTEGER :: I0_GLOBAL,J0_GLOBAL,K0_GLOBAL

  character(64) :: key
  real, pointer :: temparray(:)
  character(256) :: tempstring
  character(3) :: mstring
  real, allocatable :: tempfield2D(:,:),tempfield3D(:,:,:)

  integer :: icsmN

  print *,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
  print *, "initialize the FS for GLOBAL domain"
  WRITE(FNO_LOG,*),"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
  WRITE(FNO_LOG,*),"initialize the FS for GLOBAL domain"

  I0=I0_GLOBAL
  J0=J0_GLOBAL
  K0=K0_GLOBAL

  CALL READ_INPUT_PARAMETERS

  CALL VAR_SPECIFY_OCN_PARA
  CALL VAR_ALLOCATE_B
  CALL VAR_ALLOCATE_CGRID
  CALL VAR_ALLOCATE_ZFS
  CALL VAR_ALLOCATE_BATHY
  CALL VAR_ALLOCATE_METRIC
  CALL VAR_ALLOCATE_TAVG_MAIN
  CALL VAR_ALLOCATE_GLO2NAB_MAIN
  CALL VAR_ALLOCATE_CLIMAT
  CALL VAR_ALLOCATE_SEVP
  CALL VAR_ALLOCATE_PG
  CALL VAR_ALLOCATE_SCRATCH_MAIN
  CALL VAR_ALLOCATE_WINDMX_MAIN
  CALL VAR_ALLOCATE_XYMEANS_MAIN
  CALL VAR_ALLOCATE_GLOBALSPLT
  CALL VAR_ALLOCATE_ADD

  CALL VAR_ALLOCATE_GLOBAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1. read input files: namelist files & initial files ! 
! 2. assign S and T                                   !
! 3. assign wind shear stress                         !
! 4. specify surface evaporation and precipitation    !
! 5. specify horizontal and vertical viscosity        !
!    inclusive case-dependant special points          !
! 6. mask to eliminate biharmonic filter near shore   !
! 7. particle tracking, tracer, pmove                 !
! 8bis. compute absorption coefficients               !
! 8tris. rivers                                       !
! 9. check restart                                    !
! 10. check alternative IC                            !
! 11. initialize nudge(:,:,:)                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! step 1: INPUT DATA FILES
  CALL READ_INPUT_FILES

! step 2: assign S and T 
#ifdef FLAG_TS_T
    DO K=1,K1
      DO J=1,J0
        DO I=1,I0
          T1(I,J,K)=TCLI(I,J,K) 
          TLF(I,J,K)=TCLI(I,J,K)
          T2(I,J,K)=TCLI(I,J,K)
        END DO
      END DO
    END DO
#endif
#ifdef FLAG_TS_S
    DO K=1,K1
      DO J=1,J0
        DO I=1,I0
          S1(I,J,K)=SCLI(I,J,K)
          SLF(I,J,K)=SCLI(I,J,K)
          S2(I,J,K)=SCLI(I,J,K)
        END DO
      END DO
    END DO
#endif

! step 2.5: assign tracer
#ifdef FLAG_TRACER
  DO K=KB(I0/2,J1)/2,KB(I0/2,J1)
    DO I=2,I1
      C1(I,J0,K)=100.
      CLF(I,J0,K)=100.
      C2(I,J0,K)=100.
    END DO
  END DO
#endif

! step 3: assign wind shear stress
  bulk_formulae=.false.
  if (windtype.ge.0.or.heattype.ge.0.or.raintype.eq.1.or.salttype.ge.0) then
    !create the "grid2D" variable which corresponds to this domain's grid
    localgrid2D%typex=parall
    localgrid2D%typey=parall
    imax=size(XDEG)-2 ; jmax=size(YDEG)-2
    localgrid2D%imax=imax ; localgrid2D%jmax=jmax
    allocate(localgrid2d%x(imax,jmax),localgrid2d%y(imax,jmax),localgrid2d%mask(imax,jmax))
    do i=1,imax
      do j=1,jmax
         localgrid2D%x(i,j)=xdeg(i+1)
         localgrid2D%y(i,j)=ydeg(j+1)
         localgrid2D%mask(i,j)=1
      end do
    end do
  end if

  if (windtype.eq.-1) then
    IF (FL_WD_ON == 1) THEN
      IF (FL_RL_WD == 1) CALL REALWIND(TIMCOM_INIT)
      ! SURFACE PHYSICS INPUT DATA PROCESS
      ! Screen out land sources and convert from Newtons/m-m to dynes/cm-cm
      DO N=1,12
        DO J=1,J2
          DO I=1,I2
            TAUX(I,J,N)=10.*TAUX(I,J,N)*IN(I+1,J+1,1)
            TAUY(I,J,N)=10.*TAUY(I,J,N)*IN(I+1,J+1,1)
          END DO
        END DO
      END DO
    END IF
  elseif(windtype.eq.0) then
    bulk_formulae=.true.
  elseif (windtype.eq.1) then  ! atmos init files -> taux, tauy
    !build tauX field2D
    localMomentfluxX%initfname="interpolated from another file"
    localMomentfluxX%scale=1
    localMomentfluxX%g2D=localgrid2D
    localMomentfluxX%initialised=.false.
    localMomentfluxX%index=1
    localMomentfluxX%key="Wind stress tauX"
    allocate(localMomentfluxX%day(size(MomentfluxX%day)))
    localMomentfluxX%day=MomentfluxX%day
    !build tauY field2D
    localMomentfluxY%initfname="interpolated from another file"
    localMomentfluxY%scale=1
    localMomentfluxY%g2D=localgrid2D
    localMomentfluxY%initialised=.false.
    localMomentfluxY%index=1
    localMomentfluxY%key="Wind stress tauY"
    allocate(localMomentfluxY%day(size(MomentfluxY%day)))
    localMomentfluxY%day=MomentfluxY%day
    kmax=size(MomentfluxX%day)
    !write(*,*) "local grid minx miny maxx maxy imax jmax kmax = ",localMomentfluxX%g2D%x(1,1),localMomentfluxX%g2D%x(imax,jmax),localMomentfluxX%g2D%y(1,1),localMomentfluxX%g2D%y(imax,jmax),localMomentfluxX%g2D%imax,localMomentfluxX%g2D%jmax,kmax
    allocate(localMomentfluxX%array(imax,jmax,kmax),localMomentfluxY%array(imax,jmax,kmax))
    !interpolate from MomentfluxX to localMomentfluxX
    write(*,*) "Interpolating (in space) wind momentum --> local grid (for all timesteps)"
    do n=1,kmax
      call interp_field2D(MomentfluxX%g2D,MomentfluxX%array(:,:,n),localMomentfluxX%g2D,localMomentfluxX%array(:,:,n))
      call interp_field2D(MomentfluxY%g2D,MomentfluxY%array(:,:,n),localMomentfluxY%g2D,localMomentfluxY%array(:,:,n))
      DO J=1,J2
        DO I=1,I2
          localMomentfluxX%array(I,J,N)=10.*localMomentfluxX%array(I,J,N)*IN(I+1,J+1,1)
          localMomentfluxY%array(I,J,N)=10.*localMomentfluxY%array(I,J,N)*IN(I+1,J+1,1)
        END DO
      END DO
    end do
    !write(*,*) "        FOR THIS DOMAIN, localMomentfluxX%array(1,1,1)   = ",localMomentfluxX%array(1,1,1)
    !write(*,*) "        FOR THIS DOMAIN, localMomentfluxX%array(100,80,1) = ",localMomentfluxX%array(100,80,1)
  else
    write(*,*) "The wind forcing TYPE is not set correctly in the atmospheric forcings init file ... SERIOUS PROBLEM!!!"
  end if

  if (heattype.eq.-1) then     ! traditional TIMCOM heat flux
  elseif (heattype.eq.0) then  ! heatflux computed from T, DEWT,   CLOUDCOVERAGE, U,V and PRESSURE
    bulk_formulae=.true.
    cloudtype=1
!  elseif (heattype.eq.1) then  ! heatflux computed from T, RELHUM, CLOUDCOVERAGE, U,V and PRESSURE
!    bulk_formulae=.true.
  elseif (heattype.eq.2) then  ! heatfluxes read from files
    ! Qsolar
    localQsolar%initfname="interpolated from another file"
    localQsolar%scale=1
    localQsolar%g2D=localgrid2D
    localQsolar%initialised=.false.
    localQsolar%index=1
    localQsolar%key="Solar shortwave heat flux (downward +)"
    allocate(localQsolar%day(size(shortwavedown%day)))
    localQsolar%day=shortwavedown%day
    kmax=size(shortwavedown%day)
    allocate(localQsolar%array(imax,jmax,kmax))
    ! Qup
    localQup%initfname="interpolated from another file"
    localQup%scale=1
    localQup%g2D=localgrid2D
    localQup%initialised=.false.
    localQup%index=1
    localQup%key="Upward heat flux = Qlongwave + Qlatent + Qsensible (upward +)"
    allocate(localQup%day(size(Qup%day)))
    localQup%day=Qup%day
    kmax=size(Qup%day)
    ! interpolate on the local domain grid
    allocate(localQup%array(imax,jmax,kmax))
    write(*,*) "Interpolating (in space) heat fluxes --> local grid (for all timesteps)"
    do n=1,kmax
      call interp_field2D(shortwavedown%g2D,shortwavedown%array(:,:,n),localQsolar%g2D,localQsolar%array(:,:,n))
      call interp_field2D(Qup%g2D,Qup%array(:,:,n),localQup%g2D,localQup%array(:,:,n))
    end do
    !write(*,*) "        FOR THIS DOMAIN, localQsolar%array(100,80,1)   = ",localQsolar%array(100,80,1)
    !write(*,*) "        FOR THIS DOMAIN, localQup%array(100,80,1)   = ",localQup%array(100,80,1)
  elseif (heattype.eq.3) then ! Qshort, Qlong read from file ; Qlatent, Qsensible computed (but no need for TCC)
    bulk_formulae=.true.
    ! Qsolar
    localQsolar%initfname="interpolated from another file"
    localQsolar%scale=1
    localQsolar%g2D=localgrid2D
    localQsolar%initialised=.false.
    localQsolar%index=1
    localQsolar%key="Solar shortwave heat flux (downward +)"
    localQsolar%model%imax=imax ; localQsolar%model%jmax=jmax
    allocate(localQsolar%day(size(shortwavedown%day)))
    localQsolar%day=shortwavedown%day
    kmax=size(shortwavedown%day)
    allocate(localQsolar%array(imax,jmax,kmax))
    ! Qup
    localQup%initfname="longwave radiation interpolated from file, latent and sensible radiation computed"
    localQup%scale=1
    localQup%g2D=localgrid2D
    localQup%initialised=.false.
    localQup%index=1
    localQup%key="Longwave back-radiation flux (positive downward)"
    localQup%model%imax=imax; localQup%model%jmax=jmax
    allocate(localQup%day(size(Qup%day)))
    localQup%day=Qup%day
    kmax=size(Qup%day)
    ! interpolate Qsolar and Qup on the local domain grid
    allocate(localQup%array(imax,jmax,kmax))
    write(*,*) "Interpolating (in space) heat fluxes --> local grid (for all timesteps)"
    do n=1,kmax
      call interp_field2D(shortwavedown%g2D,shortwavedown%array(:,:,n),localQsolar%g2D,localQsolar%array(:,:,n))
      call interp_field2D(Qup%g2D,Qup%array(:,:,n),localQup%g2D,localQup%array(:,:,n))
    end do
  else
    write(*,*) "The heat forcing TYPE is not set correctly in the atmospheric forcings init file ... SERIOUS PROBLEM!!!"
  end if

! step 4: specify surface evaporation and precipitation
  if (raintype.eq.1) then
    localPrecipitation%initfname="interpolated from another file"
    localPrecipitation%scale=1
    localPrecipitation%g2D=localgrid2D
    localPrecipitation%initialised=.false.
    localPrecipitation%index=1
    localPrecipitation%key="Precipitation [kg/m2/s]"
    localPrecipitation%model%imax=imax; localPrecipitation%model%jmax=jmax
    kmax=size(Precip%day)
    allocate(localPrecipitation%day(kmax))
    localPrecipitation%day=Precip%day
    allocate(localPrecipitation%array(imax,jmax,kmax))
    ! interpolate on the local domain grid
    write(*,*) "Interpolating (in space) rain fluxes --> local grid (for all timesteps)"
    do n=1,kmax
      call interp_field2D(Precip%g2D,Precip%array(:,:,n),localPrecipitation%g2D,localPrecipitation%array(:,:,n))
    end do
  end if

  if (salttype.eq.0) then
    localEvap%initfname="interpolated from another file"
    localEvap%scale=1
    localEvap%g2D=localgrid2D
    localEvap%initialised=.false.
    localEvap%index=1
    localEvap%key="Evaporation [kg/m2/s]"
    kmax=size(EvapField%day)
    allocate(localEvap%day(kmax))
    localEvap%day=EvapField%day
    allocate(localEvap%array(imax,jmax,kmax))
    write(*,*) "Interpolating (in space) evaporation fluxes --> local grid (for all timesteps)"
    do n=1,kmax
      call interp_field2D(EvapField%g2D,EvapField%array(:,:,n),localEvap%g2D,localEvap%array(:,:,n))
    end do
  elseif (salttype.eq.1) then
    bulk_formulae=.true.
  end if

  if (bulk_formulae) then
    if (.not.allocated(taux1)) allocate(taux1(i2,j2),tauy1(i2,j2))
    localWindU10m%initfname="interpolated from another file"
    localWindU10m%scale=1
    localWindU10m%g2D=localgrid2D
    localWindU10m%initialised=.false.
    localWindU10m%index=1
    localWindU10m%key="Wind speed U10m"
    localWindU10m%model%imax=imax ; localWindU10m%model%jmax=jmax
    allocate(localWindU10m%day(size(WindU10m%day)))
    localWindU10m%day=WindU10m%day
    kmax=size(WindU10m%day)
    allocate(localWindU10m%array(imax,jmax,kmax))
    localWindV10m%initfname="interpolated from another file"
    localWindV10m%scale=1
    localWindV10m%g2D=localgrid2D
    localWindV10m%initialised=.false.
    localWindV10m%index=1
    localWindV10m%key="Wind speed V10m"
    localWindV10m%model%imax=imax ; localWindV10m%model%jmax=jmax
    allocate(localWindV10m%day(size(WindV10m%day)))
    localWindV10m%day=WindV10m%day
    kmax=size(WindV10m%day)
    allocate(localWindV10m%array(imax,jmax,kmax))
    localAirTemperature2m%initfname="interpolated from another file"
    localAirTemperature2m%scale=1
    localAirTemperature2m%g2D=localgrid2D
    localAirTemperature2m%initialised=.false.
    localAirTemperature2m%index=1
    localAirTemperature2m%key="Air Temperature at 2m"
    localAirTemperature2m%model%imax=imax ; localAirTemperature2m%model%jmax=jmax
    allocate(localAirTemperature2m%day(size(AirTemperature2m%day)))
    localAirTemperature2m%day=AirTemperature2m%day
    kmax=size(AirTemperature2m%day)
    allocate(localAirTemperature2m%array(imax,jmax,kmax))
    localDewTemperature2m%initfname="interpolated from another file"
    localDewTemperature2m%scale=1
    localDewTemperature2m%g2D=localgrid2D
    localDewTemperature2m%initialised=.false.
    localDewTemperature2m%index=1
    localDewTemperature2m%key="Dew point temperature at 2m"
    if (heattype.eq.3) localDewTemperature2m%key="Specific humidity at 2m"
    localDewTemperature2m%model%imax=imax ; localDewTemperature2m%model%jmax=jmax
    allocate(localDewTemperature2m%day(size(DewTemperature2m%day)))
    localDewTemperature2m%day=DewTemperature2m%day
    kmax=size(DewTemperature2m%day)
    allocate(localDewTemperature2m%array(imax,jmax,kmax))
    if (cloudtype.eq.1) then !(heattype.ne.3) then
      localCloudCoverage%initfname="interpolated from another file"
      localCloudCoverage%scale=1
      localCloudCoverage%g2D=localgrid2D
      localCloudCoverage%initialised=.false.
      localCloudCoverage%index=1
      localCloudCoverage%key="Cloud Coverage"
      localCloudCoverage%model%imax=imax ; localCloudCoverage%model%jmax=jmax
      allocate(localCloudCoverage%day(size(CloudCoverage%day)))
      localCloudCoverage%day=CloudCoverage%day
      kmax=size(CloudCoverage%day)
      allocate(localCloudCoverage%array(imax,jmax,kmax))
    end if
    localAtmPressure%initfname="interpolated from another file"
    localAtmPressure%scale=1
    localAtmPressure%g2D=localgrid2D
    localAtmPressure%initialised=.false.
    localAtmPressure%index=1
    localAtmPressure%key="Atmospheric Pressure"
    localAtmPressure%model%imax=imax ; localAtmPressure%model%jmax=jmax
    allocate(localAtmPressure%day(size(PresMSL%day)))
    localAtmPressure%day=PresMSL%day
    kmax=size(PresMSL%day)
    allocate(localAtmPressure%array(imax,jmax,kmax))

    write(*,*) "Interpolating remaining atmospheric fields from original grid --> local grid"
    do n=1,kmax
      call interp_field2D(WindU10m%g2D,WindU10m%array(:,:,n),localWindU10m%g2D,localWindU10m%array(:,:,n))
      call interp_field2D(WindV10m%g2D,WindV10m%array(:,:,n),localWindV10m%g2D,localWindV10m%array(:,:,n))
      call interp_field2D(AirTemperature2m%g2D,AirTemperature2m%array(:,:,n),localAirTemperature2m%g2D,localAirTemperature2m%array(:,:,n))
      call interp_field2D(PresMSL%g2D,PresMSL%array(:,:,n),localAtmPressure%g2D,localAtmPressure%array(:,:,n))
      call interp_field2D(DewTemperature2m%g2D,DewTemperature2m%array(:,:,n),localDewTemperature2m%g2D,localDewTemperature2m%array(:,:,n))
      if (cloudtype.eq.1) then
        call interp_field2D(CloudCoverage%g2D,CloudCoverage%array(:,:,n),localCloudCoverage%g2D,localCloudCoverage%array(:,:,n))
      end if
    end do

    !allocate 6 temporary fields
    allocate(tempU10M(imax,jmax),STAT=lucstat(1))
    allocate(tempV10M(imax,jmax),STAT=lucstat(2))
    allocate(tempT2M(imax,jmax),STAT=lucstat(3))
    allocate(tempDT2M(imax,jmax),STAT=lucstat(4))
    allocate(tempPMSL(imax,jmax),STAT=lucstat(5))
    allocate(tempCC(imax,jmax),STAT=lucstat(6))
    if(any(lucstat.ne.0)) STOP "Cannot allocate memory"
  end if

    DO K=1,K1
      A(K)=0.
      DO J=2,J1
        TMP=DX(J)*DY(J)
        DO I=2,I1
          A(K)=A(K)+IN(I,J,K)*TMP
        END DO
      END DO
    END DO
  IF (FL_SURF_EP == 1) THEN
    ! conversion of QAVG (deg C per model time step) to watts/m-m
    ! use factors for deg C per year and g/kg/yr
    WATTS=360.*24.*3600.*0.00134/(DT*ODZ(1))

    ! Initialize QAVG and WAVG
    N=30*DAODT
    TEMP=Z(3)/(30.*24.*3600.)
    DO M=1,12
      NSOMBO(M)=1
      MP=M+1
      IF (M.EQ.12) MP=1      
      EVAP=0.
      DO J=1,J2
        TMP=DX(J+1)*DY(J+1)
        DO I=1,I2
          ! This initialization is arbitrary and optional
          QAVG(I,J,M)=IN(I+1,J+1,1)*(TSURF(I,J,MP)-TSURF(I,J,M))/N
          WAVG(I,J,M)=-IN(I+1,J+1,1)*TEMP*(SSURF(I,J,MP)-SSURF(I,J,M))/&
                     (.5*(SSURF(I,J,M)+SSURF(I,J,MP)))
          ! Zero guess for e-p during first year (at least until month 2,
          ! when month 1 may give a good estimate)
          EVAP=EVAP-WAVG(I,J,M)*TMP
        END DO
      END DO      
      TEMP=EVAP/A(1)
      DO J=1,J2
        DO I=1,I2
          WAVG(I,J,M)=WAVG(I,J,M)+IN(I+1,J+1,1)*TEMP
        END DO
      END DO
      ! CHECK: THE FUNCTION OF THIS PART
! Yu-Chiao 10192011
!      IF ((M == 1) .AND. (FL_EVAP1 == 1)) THEN
      IF (M==1) THEN
        DO J=1,J2
          TMP=DX(J+1)*DY(J+1)
          DO I=1,I2
            EVAP=EVAP-WAVG(I,J,M)*TMP
          END DO
        END DO
        EVAP1=EVAP
      END IF       
    END DO
    ! CHECK: THE FUNCTION OF THIS PART
! Yu-Chiao 10192011
!    IF (FL_EVAP1 == 1) EVAP=EVAP1      
    DO J=1,J2
      DO I=1,I2
        W(I+1,J+1,1)=WAVG(I,J,1)
      END DO
    END DO    
  END IF

! step 5: specify horizontal and vertical viscosity
  ! step 5-1: horizontal viscosity
  ! ILE2D:  FL_TURB_H=1, 
  ! DTRAC:  FL_TURB_H=1, 
  ! GLOBAL: FL_TURB_H=1, 
  ! NPBTAI: FL_TURB_H=2, 
  
  ! IF (FL_TURB_H >= 1) THEN
  ! DMX=DM0
  ! DMY=DM0
  DMX=DM0(1:i1,1:j0,1:k1)
  DMY=DM0(1:i0,1:j1,1:k1)
  ! END IF 
  
  ! CASE DEPENDENT TREATMENT
  IF (TRIM(CASE_NAME) == 'GLOBAL') THEN
    ! increase mixing near equator and in high latitudes
    ! big northern latitude dissipation fixes initialization data problem
    ! Increasing diffusivity where KB(I,J)=1 is now superfluous, because
    ! KB > 1 except over land to allow baroclinic coastal ventillation
    ! equatorial damping parameterizes strong instabilities and augments
    ! Ri-based vertical mixing
    DO K=1,K1
      DO J=1,J1
        TMP=5.E7*(EXP(-.01*(J-J2)**2))&      ! far north
           +1.E7*(EXP(-.001*(J-J0/2-1)**2))& ! equatorial
           +5.E7*(EXP(-.01*(2-J )**2))       ! far south
        DO I=1,I1
          ! isolated bad points
          TEMP=TMP+5.E7*EXP(-.001*((I-148)**2+(J-59)**2))  &
                  +5.E7*EXP(-.001*((I-97)**2+(J-27)**2))   &
                  +2.E8*EXP(-.001*((I-21)**2+(J-32)**2))   &
                  +2.E8*EXP(-.001*((I-54)**2+(J-33)**2))   &
                  +2.E8*EXP(-.001*((I-148)**2+(J-58)**2))  &
                  +5.E8*EXP(-.0001*((I-180)**2+(J-34)**2)) &
                  +5.E8*EXP(-.0001*((I-2)**2+(J-35)**2))
          DMX(I,J,K)=DMX(I,J,K)+TEMP
          DMY(I,J,K)=DMY(I,J,K)+TEMP
        END DO
      END DO
    END DO
  END IF
  IF (TRIM(CASE_NAME) == 'MEDINA') THEN
    ! increase mixing around "bad" points
    DO J=1,J1
      DO I=1,I1
        ! isolated bad point
        TEMP=2.E8*EXP(-.001*((I-256)**2+(J-53)**2))
        do k=1,k1
          DMX(I,J,K)=DMX(I,J,K)+TEMP
          DMY(I,J,K)=DMY(I,J,K)+TEMP
        end do
      END DO
    END DO
  END IF

  IF (FL_TURB_H == 2) THEN
  ! big coastal lateral mixing
    DO JJ=2,J1; J=MAX(JJ,3)
      TMP1=2.E7; ! IF (J.GT.127) TMP1=2.E7*EXP(-.004*(J-127)**2) ! equatorial boundary where f=0 allows big wind forced flow and one must dissipate deep western undercurrent     
      TMP2=5.E5*EXP(-.004*(J-J1)**2)                           ! stormy high latitudes where natural scales are small and vigorous
      DO II=2,I1; I=MAX(II,3)
        ZA=0.; ZB=0.; ZC=0.; ZD=0.
        IF (KB(I+1,JJ)/=0) ZA=Z(2*KB(I+1,JJ))
        IF (KB(I-1,JJ)/=0) ZB=Z(2*KB(I-1,JJ))
        IF (KB(II,J+1)/=0) ZC=Z(2*KB(II,J+1))
        IF (KB(II,J-1)/=0) ZD=Z(2*KB(II,J-1))
        DZDX=(ZA-ZB)*ODX(J); DZDY=(ZC-ZD)*ODY(J)  
        TMP=SQRT(DZDX**2+DZDY**2)

        ZBOT=Z(2*KB(I,J)+1)
        DO K=1,KB(I,J)
          !SCR(II,JJ,K)=DM0+TMP1+TMP2+5.E6*EXP(-1.E-9*(ZBOT-Z(2*K+1))**2)*TMP
          SCR(II,JJ,K)=TMP1+TMP2+5.E6*EXP(-1.E-9*(ZBOT-Z(2*K+1))**2)*TMP
        END DO
      END DO
    END DO
    
    DO K=1,K1; DO N=1,1; DO J=2,J1; DO I=3,I1
      M=IN(I+1,J,K)*IN(I-1,J,K)*IN(I,J+1,K)*IN(I,J-1,K)
      IF (M.EQ.0.AND.IN(I,J,K).EQ.1) SCR(I,J,K)=SCR(I,J,K)+1.E6
    END DO; END DO; END DO; END DO
    
    ! average, 10 times
    DO K=1,K1; DO N=1,10; DO J=2,J1; DO I=3,I1
      M=IN(I+1,J,K)*IN(I-1,J,K)*IN(I,J+1,K)*IN(I,J-1,K)
      IF (M.EQ.1) SCR(I,J,K) =.25*(SCR(I+1,J,K)+SCR(I-1,J,K)+SCR(I,J+1,K)+SCR(I,J-1,K))
    END DO; END DO; END DO; END DO;
          
    DO K=1,K1
      DO J=2,J1
        DMX(1,J,K) =DMX(1,J,K) +SCR(2,J,K)
        DO I=2,I2
        DMX(I,J,K) =DMX(I,J,K) +.5*(SCR(I,J,K)+SCR(I+1,J,K))
        END DO
        DMX(I1,J,K)=DMX(I1,J,K)+SCR(I2,J,K)
      END DO
      ! note: except in hi lats, DMY slightly violates ang. mom. conservation
      DO I=2,I1
        DMY(I,1,K) =DMY(I,1,K) +SCR(I,2,K)
        DO J=2,J2
        DMY(I,J,K) =DMY(I,J,K) +.5*(SCR(I,J,K)+SCR(I,J+1,K))
        END DO
        DMY(I,J1,K)=DMY(I,J1,K)+SCR(I,J2,K)
      END DO
    END DO     
  END IF    
      
  ! step 5-2: vertical viscosity
  ! CURRENTLY, ONLY FOR NPBTAI
  IF (FL_TURB_V_ADD == 1) THEN
    DO JJ=2,J1; J=MAX(JJ,3)
      TMP3=500.; ! IF (JJ.GT.127) TMP3=500.*EXP(-.004*(JJ-127)**2)
      DO II=2,I1; I=MAX(II,3)
        ZA=0.; ZB=0.; ZC=0.; ZD=0.
        IF (KB(I+1,JJ)/=0) ZA=Z(2*KB(I+1,JJ))
        IF (KB(I-1,JJ)/=0) ZB=Z(2*KB(I-1,JJ))
        IF (KB(II,J+1)/=0) ZC=Z(2*KB(II,J+1))
        IF (KB(II,J-1)/=0) ZD=Z(2*KB(II,J-1))
        DZDX=(ZA-ZB)*ODX(J); DZDY=(ZC-ZD)*ODY(J)  
        TMP=SQRT(DZDX**2+DZDY**2)

        ZBOT=Z(2*KB(I,J)+1)
        DO K=1,K2
          TMP4=99.*EXP(-1.E-9*(ZBOT-Z(2*K+1))**2+3.E-5*(ZBOT-4.E5))
          ADD(II-1,JJ-1,K)=0.2*VBK(II-1,JJ-1,K)+TMP3+TMP4+99.*EXP(-1.E-9*(ZBOT-Z(2*K+1))**2)*TMP
        END DO          
      END DO
    END DO      
  END IF

! step 6: mask to eliminate biharmonic filter
  IF (FL_BI_FIL == 1) THEN
    ! Mask to eliminate biharmonic filter application near shore
    DO K=1,K1
      DO J=2,J1
        IF (LOPENE == 1) THEN ! PERIODIC B.C.
          DO I=3,I0-2
            INFX(I,J,K)=IN(I-2,J,K)*IN(I-1,J,K)*IN(I,J,K)*IN(I+1,J,K)*IN(I+2,J,K)
          END DO
        ELSE ! CLOSED B.C. OR OPEN B.C.
          DO I=4,I0-3
            INFX(I,J,K)=IN(I-2,J,K)*IN(I-1,J,K)*IN(I,J,K)*IN(I+1,J,K)*IN(I+2,J,K)
          END DO
        END IF
      END DO      
      DO J=4,J0-3
        DO I=2,I1
          INFY(I,J,K)=IN(I,J-2,K)*IN(I,J-1,K)*IN(I,J,K)*IN(I,J+1,K)*IN(I,J+2,K)
        END DO
      END DO
    END DO
  END IF   
! step 7: calculate AROUT for to avoid incompressibility error
  AROUT=0.
  IF(FL_INFL_W.GE.3)THEN 
    DO K=1,K1
      DO J=2,J1
        AROUT=AROUT+IN(2,J,K)/(ODY(J)*ODZ(K)) 
      END DO
    END DO
  ENDIF
  IF(FL_INFL_E.GE.3)THEN 
    DO K=1,K1
      DO J=2,J1
        AROUT=AROUT+IN(I1,J,K)/(ODY(J)*ODZ(K)) 
      END DO
    END DO
  ENDIF
  IF(FL_INFL_S.GE.3)THEN 
    DO K=1,K1
      DO I=2,I1
        AROUT=AROUT+IN(I,2,K)/(ODXV(1)*ODZ(K))                    
      END DO
    END DO
  ENDIF
  IF(FL_INFL_N.GE.3)THEN 
    DO K=1,K1
      DO I=2,I1
        AROUT=AROUT+IN(I,J1,K)/(ODXV(J1)*ODZ(K))                    
      END DO
    END DO
  ENDIF
  WRITE(FNO_LOG,*)'AROUT=',AROUT
! step 8: particle tracking
  IF (FL_PMOVE_ON == 1) CALL PMOVE(TIMCOM_INIT) !GLOBAL, NPB, TAI  
! Write formatted run description and parameters
! Yu-Chiao 10192011
!  WRITE(FNO_LOG,1) 'GLOBAL',DSCRIB,MXIT,I0,J0,K0,KTRM,DT,FLTW,G,ORZMX,DRAG
! 1  FORMAT('CASE ',A5,' (',A63,')'/'----------'/&
!       ' MXIT    I0    J0   K0   KTRM'/1I7,4I5//&
!       2X,'DT',8X,'FLTW',6X,'G',6X,'ORZMX',6X,'DRAG'/1P,5(1X,E9.2))
  WRITE(FNO_LOG,1) 'GLOBAL',DSCRIB,MXIT,I0,J0,K0,DT,FLTW,G,ORZMX,DRAG
 1  FORMAT('CASE ',A5,' (',A63,')'/'----------'/&
       ' MXIT    I0    J0   K0   '/1I7,3I5//&
       2X,'DT',8X,'FLTW',6X,'G',6X,'ORZMX',6X,'DRAG'/1P,5(1X,E9.2))

! step 8bis; absorption
  call var_allocate_absorption
  call CALC_ABSORPTION_RATE

! step 8tris: river data
  write(*,*) "Loading river data specified in :",trim(riverfile)
  rivercount=0
  if (presentInitValue(riverfile,'rivers')) then
     call getInitValue(riverfile,'rivers',riverlist)
     rivercount=size(riverlist)
     allocate(riveri(rivercount),riverj(rivercount),riverD(rivercount),riverFlow(rivercount,12),riverT(rivercount,12),riverS(rivercount,12))
  end if
  do i=1,rivercount
     key=trim(adjustl(riverlist(i)))//'.i'
       call getInitValue(riverfile,key,riveri(i))
     key=trim(adjustl(riverlist(i)))//'.j'
       call getInitValue(riverfile,key,riverj(i))
     key=trim(adjustl(riverlist(i)))//'.direction'
       call getInitValue(riverfile,key,riverD(i))
     key=trim(adjustl(riverlist(i)))//'.flow'
       call getInitValue(riverfile,key,temparray)
       riverFlow(i,:)=temparray
       deallocate(temparray)
     key=trim(adjustl(riverlist(i)))//'.temp'
       call getInitValue(riverfile,key,temparray)
       riverT(i,:)=temparray
       deallocate(temparray)
     key=trim(adjustl(riverlist(i)))//'.sal'
       call getInitValue(riverfile,key,temparray)
       riverS(i,:)=temparray
       deallocate(temparray)
  end do

! step 9: Check restart
#ifndef FLAG_TS_S
    ! NEED FURTHER MODIFICATION
    ! THERE WAS A BUG IN TIMCOM, I.E. NOT CONSISTENT W/ DIECAST
    IF (RESTARTSIG .NE. 0) THEN
      CALL RD_RESTART
      WRITE(*,92) ITF,DAYS
      I=DAYS*DAODT+.5
      IF (I /= ITF) THEN
        WRITE(*,93) ITF,I
        ITF=I
      END IF
    END IF
#else
#ifdef FLAG_TS_T
    RESTART_CONDITION:IF (RESTARTSIG.NE.0) THEN 
      CALL RD_RESTART
      WRITE(*,92) ITF,DAYS
      WRITE(FNO_LOG,40)
   40 FORMAT('climatological means')
      DO K=1,K1
        TMP=0.
        TEMP=0.
        N=0
        DO J=2,J1
          DO I=2,I1
            N=N+IN(I,J,K)
            TMP=TMP+IN(I,J,K)*T1(I,J,K)
            TEMP=TEMP+IN(I,J,K)*S1(I,J,K)
          END DO
        END DO
        IF (N.NE.0) TMP=TMP/N
        IF (N.NE.0) TEMP=TEMP/N
        WRITE(FNO_LOG,44) K,TMP,TEMP
   44 FORMAT('level',I3,' mean T,S=',2F6.2)
      END DO    
      N=DAYS*DAODT+.5
      IF (N /= ITF) THEN 
        WRITE(*,93) ITF,N
        DAODT_NEW=DAODT
        DAODT_OLD=ITF/DAYS
        TMP=DAODT_OLD/DAODT_NEW
        DO NM=1,12
          DO J=1,J2
            DO I=1,I2
              QAVG(I,J,NM)=TMP*QAVG(I,J,NM)
            END DO
          END DO
        END DO
        ITF=N
      END IF
      IF (DAYS.LT.3241.) AV=0.
      ! Make sure W and WAVG have zero horizontal mean
      ! (this allows one to change surface geometry "on the fly")
      DO N=1,12
        EVAP=0.
        EVAP2=0.
        DO J=1,J2
          TMP=DX(J+1)*DY(J+1)
          DO I=1,I2
            EVAP2=EVAP2-IN(I+1,J+1,1)*W(I+1,J+1,1)*TMP
            EVAP=EVAP-IN(I+1,J+1,1)*WAVG(I,J,N)*TMP
          END DO
        END DO
        ! Modify W(I,J,1) to assure no net volume inflow (avg EVAP=0)
        TMP=EVAP/A(1)
        TMP2=EVAP2/A(1)
        DO J=1,J2
          DO I=1,I2
            W(I+1,J+1,1)=IN(I+1,J+1,1)*(W(I+1,J+1,1)+TMP2)
            WAVG(I,J,N)=IN(I+1,J+1,1)*(WAVG(I,J,N)+TMP)
          END DO
        END DO
      END DO
      IF (DAYS.LE.3601.) THEN 
        ! WAVG & QAVG include advection effects
        ! WAVG is surface freshwater sources/sinks (e-p) consistent with model
        ! and climatological data
        ! NSOMBO(M)=number of ensemble-averaging years
        !           assumed already included in QAVG and WAVG
        ! However, we weight each new year equally to the previously completed
        ! year's ensemble average, until year 11, thus given more weight to
        ! recent results and thus acclerating convergence to a multi-year ensemble
        ! average.
        DO M=1,12
          NSOMBO(M)=1
        END DO
      END IF
      WRITE(*,90) ITF,DAYS,NSOMBO
      WRITE(FNO_LOG,90) ITF,DAYS,NSOMBO
   90 FORMAT(/"Restart data read for GLOBAL's ITF ",I7,", day ",F8.2,", NSOMBO:",12I3)
    END IF RESTART_CONDITION 
! Yu-Chiao 10152011
!  END IF
#endif
#endif
 92   FORMAT("Restart data read for GLOBAL's ITF ",I6,", day ",F8.2)
 93   FORMAT(/'GLOBAL ITF # changed from',I7,' to',I7,' at restart to reflect DT change')

  IF (ITF.GE.MXIT) THEN
    WRITE(*,*)"INITFS_GLOBAL: ITF ",ITF, ">= MAX TIMESTEPS", MXIT
    WRITE(*,*)"I'd be stop better!"
    STOPSIG=TIMCOM_FINISH
    CALL TIMCOM_STOPSIG_SIGNAL
  END IF

!step 10. check alternative IC
if (icSIG.eq.1) then
  ! load from icfile
  write(mstring,'(i3.3)')member
  write(*,*) "Now loading IC from file: ",trim(adjustl(icfile)),'/GLOBAL'//mstring//'.nc, variables : T S U V P0'
  tempstring=trim(icfile)//'/GLOBAL'//mstring//'.nc#T'
  call uload(tempstring,tempfield3D)
  if (size(tempfield3D,1).eq.I0.and.size(tempfield3D,2).eq.J0.and.size(tempfield3D,3).eq.K1) then
    T1=tempfield3D
    T2=tempfield3D
    TLF=tempfield3D
    deallocate(tempfield3D)
  else
    write(*,*) "IC file, variable T: variable size do not correspond"
    stopsig=TIMCOM_FINISH
    call TIMCOM_STOPSIG_SIGNAL
  end if
  tempstring=trim(icfile)//'/GLOBAL'//mstring//'.nc#S'
  call uload(tempstring,tempfield3D)
  if (size(tempfield3D,1).eq.I0.and.size(tempfield3D,2).eq.J0.and.size(tempfield3D,3).eq.K1) then
    S1=tempfield3D
    S2=tempfield3D
    SLF=tempfield3D
    deallocate(tempfield3D)
  else
    write(*,*) "IC file, variable S: variable size do not correspond"
    stopsig=TIMCOM_FINISH
    call TIMCOM_STOPSIG_SIGNAL
  end if
  tempstring=trim(icfile)//'/GLOBAL'//mstring//'.nc#U'
  call uload(tempstring,tempfield3D)
  if (size(tempfield3D,1).eq.I0.and.size(tempfield3D,2).eq.J0.and.size(tempfield3D,3).eq.K1) then
    U1=tempfield3D
    U2=tempfield3D
    ULF=tempfield3D
    deallocate(tempfield3D)
  else
    write(*,*) "IC file, variable U: variable size do not correspond"
    stopsig=TIMCOM_FINISH
    call TIMCOM_STOPSIG_SIGNAL
  end if
  tempstring=trim(icfile)//'/GLOBAL'//mstring//'.nc#V'
  call uload(tempstring,tempfield3D)
  if (size(tempfield3D,1).eq.I0.and.size(tempfield3D,2).eq.J0.and.size(tempfield3D,3).eq.K1) then
    V1=tempfield3D
    V2=tempfield3D
    VLF=tempfield3D
    deallocate(tempfield3D)
  else
    write(*,*) "IC file, variable V: variable size do not correspond"
    stopsig=TIMCOM_FINISH
    call TIMCOM_STOPSIG_SIGNAL
  end if
  tempstring=trim(icfile)//'/GLOBAL'//mstring//'.nc#P0'
  call uload(tempstring,tempfield2D)
  if (size(tempfield2D,1).eq.I0.and.size(tempfield2D,2).eq.J0) then
    !P0=tempfield2D*100.0*G      !loading SSH into P0
    !P=tempfield3D                !this works but we usually don't have P
    P0=tempfield2D
    deallocate(tempfield2D)
  else
    write(*,*) "IC file, variable SSH: variable size do not correspond"
    stopsig=TIMCOM_FINISH
    call TIMCOM_STOPSIG_SIGNAL
  end if

  if (smoothSIG.eq.1) then
  !smoothing the initial condition
    do k=1,k1
      icsmN=ics_top+(k-1)*(1.0*ics_bottom-ics_top)/(1.0*k1-1.0)
      WRITE(*,20) k,icsmN
20    FORMAT('smoothing T and S, layer ',I3,' ',I3,' iterations')
      DO N=1,icsmN
        SMIN=42.
        SMAX=0.
        TMIN=42.
        TMAX=0.
        do J=2,J1
          do I=2,I1
            M=IN(I-1,J,K)+IN(I+1,J,K)+IN(I,J-1,K)+IN(I,J+1,K)
            if (M*IN(I,J,K).NE.0) then
              TMP=1./M
              T1(I,J,K)=TMP*(IN(I+1,J,K)*T1(I+1,J,K)+IN(I-1,J,K)*T1(I-1,J,K) &
                +IN(I,J+1,K)*T1(I,J+1,K)+IN(I,J-1,K)*T1(I,J-1,K))
              S1(I,J,K)=TMP*(IN(I+1,J,K)*S1(I+1,J,K)+IN(I-1,J,K)*S1(I-1,J,K) &
                +IN(I,J+1,K)*S1(I,J+1,K)+IN(I,J-1,K)*S1(I,J-1,K))
              SMIN=MIN(SMIN,S1(I,J,K))
              SMAX=MAX(SMAX,S1(I,J,K))
              TMIN=MIN(TMIN,T1(I,J,K))
              TMAX=MAX(TMAX,T1(I,J,K))
            end if
          end do
        end do
!       WRITE(17,31) N,K,SMIN,SMAX,TMIN,TMAX
!31      FORMAT('   N,K,SMIN,SMAX,TMIN,TMAX:',I4,I3,4F6.2)
      end do
    end do
    T2=T1
    S2=S1
    TLF=T1
    SLF=S1
  end if

  if (stabilizeSIG.eq.1) then
    call EOS_COMPUTATION_STABILITY1(T1,S1,Z,RHO,FL_EOS_OP,0,0)
    if (icSIG.eq.1) then
      T2=T1 ; TLF=T1 ; S2=S1 ; SLF=S1
    end if
    write(*,*) "Loaded T and S were modified for water columns to be stable"
  end if

  write(*,*) "  ... loading succesfull !"
end if

!step 11: initialize nudging
  nudge=1

  print *,"End of initfs for GLOBAL domain"
  print *,">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  WRITE(FNO_LOG,*),"End of initfs for GLOBAL domain"
  WRITE(FNO_LOG,*),">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

CONTAINS

  INCLUDE './realwind.f90'
  INCLUDE './pmove.f90'

END SUBROUTINE INITFS_GLOBAL
