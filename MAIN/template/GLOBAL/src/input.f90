MODULE INPUT
USE TIMCOM_GENERAL
USE OCN_PARA
USE date
USE initfile
CONTAINS
SUBROUTINE READ_CASE_INFO
OPEN(20,FILE=".case_info",CONVERT="BIG_ENDIAN",STATUS="OLD",FORM="UNFORMATTED")
READ(20)CASE_NAME
READ(20)DSCRIB
CLOSE(20)
END SUBROUTINE READ_CASE_INFO

SUBROUTINE READ_CASE_SIZE
OPEN(20,FILE="./INPUT_GLOBAL/GLOBAL",CONVERT="BIG_ENDIAN",STATUS="OLD",FORM="UNFORMATTED")
READ(20)I0_GLOBAL,J0_GLOBAL,K0_GLOBAL
CLOSE(20)
END SUBROUTINE READ_CASE_SIZE

SUBROUTINE READ_NAMELIST_INPUT
CHARACTER(15)::tempname
write(tempname,"(A,i3.3,A)") "namelist",member,".run"
write(*,*) "reading general configuration: ",tempname
call getInitValue(tempname,'WORKDIR',WORKDIR)
call getInitValue(tempname,'RUNS',RUNS)
call getInitValue(tempname,'ISAV',ISAV)
call getInitValue(tempname,'VAR_I0J0K1_OUTPUT',VAR_I0J0K1_OUTPUT)
call getInitValue(tempname,'VAR_I0J0K1_HDF5',VAR_I0J0K1_HDF5)
call getInitValue(tempname,'UVAR_OUTPUT',UVAR_OUTPUT)
call getInitValue(tempname,'VVAR_OUTPUT',VVAR_OUTPUT)
call getInitValue(tempname,'SVAR_OUTPUT',SVAR_OUTPUT)
call getInitValue(tempname,'TVAR_OUTPUT',TVAR_OUTPUT)
call getInitValue(tempname,'PVAR_OUTPUT',PVAR_OUTPUT)
call getInitValue(tempname,'VAR_I0J0K0_OUTPUT',VAR_I0J0K0_OUTPUT)
call getInitValue(tempname,'WVAR_OUTPUT',WVAR_OUTPUT)
if (presentInitValue(tempname,'VAR_FLUXES_OUTPUT')) &
  call getInitValue(tempname,'VAR_FLUXES_OUTPUT',VAR_FLUXES_OUTPUT)
call getInitValue(tempname,'SCRNOUT',SCRNOUT)
call getInitValue(tempname,'CPSAV',CPSAV)
call getInitValue(tempname,'ATMOSFILE',ATMOSFILE)
call getInitValue(tempname,'T0YEAR',T0YEAR)
call getInitValue(tempname,'T0MONTH',T0MONTH)
call getInitValue(tempname,'T0DAY',T0DAY)
call getInitValue(tempname,'T0HOUR',T0HOUR)
call getInitValue(tempname,'T0MINUT',T0MINUT)
t0mjd=mjd2(T0YEAR,T0MONTH,T0DAY,T0HOUR,T0MINUT)
END SUBROUTINE READ_NAMELIST_INPUT
END MODULE INPUT
