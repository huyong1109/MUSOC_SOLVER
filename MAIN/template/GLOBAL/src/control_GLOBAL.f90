MODULE CONTROL_GLOBAL
! User-defined (scalar and array) control parameters

REAL::EVAP,QSUM,OLDE,OLDQ,AV,AVN,NAVG ! ,TAUDT,TAUDTN
INTEGER :: MXIT,MXITN


INTEGER :: DIMEN                 !NAMELIST:SPATIAL

REAL :: DAODT,DT                 !NAMELIST:TEMPORAL

INTEGER :: LEVP,MXMASK,FL_EVP_STP !NAMELIST:SOLVER
REAL :: TOLERANCE                 !NAMELIST:SOLVER

REAL :: DMZ0,TL,RZMX,DRAG    !NAMELIST:PHYSICS
REAL*4, allocatable :: DM0(:,:,:),DE0(:,:,:)
	
INTEGER :: FL_TURB_H, FL_TURB_V, FL_TURB_V_ADD, ITF_PHYS  ! PHYSICS_FLAGS

INTEGER :: FL_RL_TS,FL_EOS_OP,FL_TRACER_ON    ! TSC_FLAGS
REAL :: TMAX

INTEGER :: FL_PMOVE_ON,FL_PMOVE_NEWDAY           ! PMOVE_FLAGS

INTEGER :: FL_SURF_EP,FL_EVAP1,FL_WD_ON, FL_RL_WD    ! SURFACE_BOUND_FLAGS

!NAMELIST:	NUDGE_FACTOR
INTEGER :: FL_FNEWFOLD_CALC
INTEGER :: FL_NUDGE_HMEAN_T,FL_NUDGE_HMEAN_S
INTEGER :: FL_NUDGE_OP_T,FL_NUDGE_OP_S,FL_ARBR_P0
REAL :: NUDGING_HMEAN,NUDGE_W
real, allocatable ::	tau2D(:,:), tau3D(:,:,:)
	
INTEGER :: LOPENW,LOPENE,LOPENS,LOPENN             !NAMELIST:LATERAL_BOUND_FLAGS 
INTEGER :: FL_INFL_W,FL_INFL_E,FL_INFL_S,FL_INFL_N !NAMELIST:LATERAL_BOUND_FLAGS

REAL :: UINFLOW,VINFLOW                            !NAMELIST:INFLOW_BOUND_COND


INTEGER :: FL_SPL_ON,FL_SWAMP_RE                   !NAMELIST:SPONGE_AND_SWAMP_LAYERS

REAL :: FLTW                    !NAMELIST:FILTERS
REAL :: WRAF                    !NAMELIST:FILTERS
INTEGER :: FL_BI_FIL            !NAMELIST:FILTERS

INTEGER :: FL_INI_VIS, FL_BK_VIS

REAL ::	 modelmjd
integer :: year,month,day,hour,minut
	
character(256) :: RIVERFILE
integer :: rivercount
character(len=10), pointer :: riverlist(:)
integer, pointer :: riveri(:),riverj(:)
character,pointer:: riverD(:)
real,pointer     :: riverFlow(:,:),riverT(:,:),riverS(:,:)
	
END MODULE CONTROL_GLOBAL


!######################################################################


MODULE SCA_GLOBAL
! Derived scalars
REAL :: ODT   ! NAMELIST:TEMPORAL
!REAL :: PRN   ! NAMELIST:PHYSICS
! prn used to be de0/dm0, but is not stored anymore since dm0 and de0 are 3D arrays
! prn was used only in fs_GLOBAL and is replaced by its value dm0(i,j,)/de0(i,j,k)
REAL::DAYS,ORZMX,OFLTW,WATTS
INTEGER::N360=0,NYR=0,IT0=0,ITF=0,ITFDAY=0
END MODULE SCA_GLOBAL


!######################################################################


MODULE DIR_GLOBAL
CHARACTER(LEN=30) :: IDIR="/INPUT_GLOBAL/", &
                     ODIR="/OUTPUT_GLOBAL/",&
                     TDIR="/TEMP_GLOBAL/"
END MODULE DIR_GLOBAL
