SUBROUTINE INITFS()
! ----------------------------------------------------------------------
use TIMCOM_GENERAL 
use OCN_PARA
use INPUT

!===============================
! READ THE CASE INFO 
!===============================
  CALL READ_CASE_INFO
  CALL READ_CASE_SIZE
!===============================
! READ THE GENERAL NAMELIST FILE
!===============================
  CALL READ_NAMELIST_INPUT

END SUBROUTINE INITFS
