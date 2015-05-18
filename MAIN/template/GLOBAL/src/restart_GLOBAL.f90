MODULE RESTART_GLOBAL
use OCN_PARA_GLOBAL
use SCA_GLOBAL
use CLIMAT_GLOBAL
use CONTROL_GLOBAL
!use B_GLOBAL
USE GRID_VAR_GLOBAL, ONLY: RHO,Mcrho,McTP,McS,McTPTP,McSTP,P0,PX,PY,U1,U2,V1,V2,S1,S2,T1,T2,P,ULF,VLF,SLF,TLF,EV,HV,F,TANPHI,SUMIN
!use CGRID_GLOBAL
USE GRID_VAR_GLOBAL, ONLY: U,V,W
!use TAVG_GLOBAL_MAIN
USE GRID_VAR_GLOBAL, ONLY: PBAR,PVAR,XBAR,UCLI,VCLI,RMSV,SBAR,TBAR
USE DIR_GLOBAL,ONLY: TDIR
USE GRID_VAR_GLOBAL, ONLY:X,QAVG,TNUDGE,WAVG
use TIMCOM_GENERAL, only : member
use ufileformat

IMPLICIT NONE

CHARACTER(160),PRIVATE::FN,FN_CHK
CHARACTER(256),PRIVATE::FNnc
CHARACTER(1024),PRIVATE :: CMD
CHARACTER(3) :: member_string
INTEGER,PRIVATE::FNO=19
REAL,PRIVATE::	valex
INTEGER,PRIVATE:: valexi
REAL, ALLOCATABLE:: tempreal(:)		
CONTAINS

! ------------------
! WRITE RESTART DATA
! ------------------
SUBROUTINE WR_RESTART

  print *, "---------------------------------------------------------------------------"
  print *, "Checkpoint to RESTART file:"

  write(member_string,"(i3.3)") member
  FN=TRIM(WORKDIR)//TRIM(TDIR)//'restart'//member_string//'.nc'
  FN_CHK=TRIM(WORKDIR)//TRIM(TDIR)//'restart'//member_string//'.completed'

! mv restart.complete to restart.complete.bak
  CMD='if [ -f "'//TRIM(FN_CHK)//'" ]; then mv '//TRIM(FN_CHK)//' '//TRIM(FN_CHK)//'.bak ; fi'
  CALL SYSTEM(TRIM(CMD))

! mv restart to restart.bak
  CMD='if [ -f "'//TRIM(FN)//'" ]; then mv '//TRIM(FN)//' '//TRIM(FN)//'.bak ; fi'
  CALL SYSTEM(TRIM(CMD))

!! A VERSION WITH BINARY RESTART
!  OPEN(FNO,CONVERT='BIG_ENDIAN',FILE=FN,form='unformatted')
!! dynamic restart data
!  WRITE(FNO) ITF,NSOMBO,DAYS,AV,X,&
!            U1,U2,V1,V2,S1,S2,T1,T2,P0,P,ULF,VLF,SLF,TLF,U,V,W,&
!            QAVG,TNUDGE,WAVG,EVAP,QSUM,OLDE,OLDQ
!! model time mean diagnostic data
!  WRITE(FNO) PBAR,PVAR,XBAR,SBAR,TBAR,UCLI,VCLI,RMSV,SSURFM,TSURFM
!  CLOSE(FNO)
!
!! A VERSION WITH NETCDF RESTART
   valex=9999.
   FNnc=trim(FN)//'#MJD' ;    call usave(FNnc,modelmjd,valex)
   FNnc=trim(FN)//'#ITF' ;    call usave(FNnc,ITF,valex)
   FNnc=trim(FN)//'#NSOMBO' ;    call usave(FNnc,NSOMBO,valex)
   FNnc=trim(FN)//'#DAYS' ;    call usave(FNnc,DAYS,valex)
   FNnc=trim(FN)//'#AV'	;    call usave(FNnc,AV,valex)
   FNnc=trim(FN)//'#X'	;    call usave(FNnc,X,valex)
   FNnc=trim(FN)//'#U1'	;    call usave(FNnc,U1,valex)
   FNnc=trim(FN)//'#U2'	;    call usave(FNnc,U2,valex)
   FNnc=trim(FN)//'#V1'	;    call usave(FNnc,V1,valex)
   FNnc=trim(FN)//'#V2'	;    call usave(FNnc,V2,valex)
   FNnc=trim(FN)//'#S1'	;    call usave(FNnc,S1,valex)
   FNnc=trim(FN)//'#S2'	;    call usave(FNnc,S2,valex)
   FNnc=trim(FN)//'#T1'	;    call usave(FNnc,T1,valex)
   FNnc=trim(FN)//'#T2'	;    call usave(FNnc,T2,valex)
   FNnc=trim(FN)//'#P0'	;    call usave(FNnc,P0,valex)
   FNnc=trim(FN)//'#P'	;    call usave(FNnc,P,valex)
   FNnc=trim(FN)//'#ULF' ;    call usave(FNnc,ULF,valex)
   FNnc=trim(FN)//'#VLF' ;    call usave(FNnc,VLF,valex)
   FNnc=trim(FN)//'#SLF' ;    call usave(FNnc,SLF,valex)
   FNnc=trim(FN)//'#TLF' ;    call usave(FNnc,TLF,valex)
   FNnc=trim(FN)//'#U'	;    call usave(FNnc,U,valex)
   FNnc=trim(FN)//'#V'	;    call usave(FNnc,V,valex)
   FNnc=trim(FN)//'#W'	;    call usave(FNnc,W,valex)
   FNnc=trim(FN)//'#QAVG'   	;    call usave(FNnc,QAVG,valex)
   FNnc=trim(FN)//'#TNUDGE' ;    call usave(FNnc,TNUDGE,valex)
   FNnc=trim(FN)//'#WAVG' ;    call usave(FNnc,WAVG,valex)
   FNnc=trim(FN)//'#EVAP' ;    call usave(FNnc,EVAP,valex)
   FNnc=trim(FN)//'#QSUM' ;    call usave(FNnc,QSUM,valex)
   FNnc=trim(FN)//'#OLDE' ;    call usave(FNnc,OLDE,valex)
   FNnc=trim(FN)//'#OLDQ' ;    call usave(FNnc,OLDQ,valex)
   FNnc=trim(FN)//'#PBAR' ;    call usave(FNnc,PBAR,valex)
   FNnc=trim(FN)//'#PVAR' ;    call usave(FNnc,PVAR,valex)
   FNnc=trim(FN)//'#XBAR' ;    call usave(FNnc,XBAR,valex)
   FNnc=trim(FN)//'#SBAR' ;    call usave(FNnc,SBAR,valex)
   FNnc=trim(FN)//'#TBAR' ;    call usave(FNnc,TBAR,valex)
   FNnc=trim(FN)//'#UCLI' ;    call usave(FNnc,UCLI,valex)
   FNnc=trim(FN)//'#VCLI' ;    call usave(FNnc,VCLI,valex)
   FNnc=trim(FN)//'#RMSV' ;    call usave(FNnc,RMSV,valex)
   FNnc=trim(FN)//'#SSURFM' ;    call usave(FNnc,SSURFM,valex)
   FNnc=trim(FN)//'#TSURFM' ;    call usave(FNnc,TSURFM,valex)

	! make new restart.complete when restart checkpoint complete
  WRITE(CMD,'(A,I10,A,A)')'echo ',ITF_GENERAL,' >& ',TRIM(FN_CHK)
  CALL SYSTEM(TRIM(CMD))

  print *, "End of checkpoint !"
  print *, "---------------------------------------------------------------------------"

END SUBROUTINE WR_RESTART


	
	
! -----------------
! READ RESTART DATA
! -----------------
SUBROUTINE RD_RESTART

integer :: i,j,k
real :: d_r_var1(k0,2),d_r_var2(J0,K1),d_r_var3(i0,j0,12)
! d_r_var = dummy real variable

  print *, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
  print *, "Reading from RESTART:"

  write(member_string,"(i3.3)") member
  FN=TRIM(WORKDIR)//TRIM(TDIR)//'restart'//member_string//'.nc'

!! THE VERSION READS FROM BINARY RESTART FILE
!  OPEN(FNO,CONVERT='BIG_ENDIAN',FILE=FN,form='unformatted')
!! AV is for long term averages. AVG is for nabcs 5-day averages.
!  READ(FNO) ITF,NSOMBO,DAYS,AV,X,&
!      U1,U2,V1,V2,S1,S2,T1,T2,P0,P,ULF,VLF,SLF,TLF,U,V,W,&
!      QAVG,TNUDGE,WAVG,EVAP,QSUM,OLDE,OLDQ
!! model time mean diagnostic data
!  READ(FNO) PBAR,PVAR,XBAR,SBAR,TBAR,UCLI,VCLI,RMSV,SSURFM,TSURFM
!  CLOSE(FNO)
!!
!! THE VERSION READS FROM NETCDF FILE
   FNnc=trim(FN)//'#ITF'    ;     call uload(FNnc,ITF)
write(*,*) "ITF=",ITF
   FNnc=trim(FN)//'#NSOMBO' ;     call uload(FNnc,tempreal)
            NSOMBO=tempreal ;     deallocate(tempreal)
   FNnc=trim(FN)//'#DAYS'   ;     call uload(FNnc,DAYS)
   FNnc=trim(FN)//'#AV'	    ;     call uload(FNnc,AV)
   FNnc=trim(FN)//'#X'	    ;     deallocate(X)    ;call uload(FNnc,X,valex)
write(*,*) "X(180,50)=",X(180,50)
   FNnc=trim(FN)//'#U1'	    ;     deallocate(U1)   ; call uload(FNnc,U1,valex)
   FNnc=trim(FN)//'#U2'	    ;     deallocate(U2)   ; call uload(FNnc,U2,valex)
   FNnc=trim(FN)//'#V1'	    ;     deallocate(V1)   ; call uload(FNnc,V1,valex)
   FNnc=trim(FN)//'#V2'	    ;     deallocate(V2)   ; call uload(FNnc,V2,valex)
   FNnc=trim(FN)//'#S1'	    ;     deallocate(S1)   ; call uload(FNnc,S1,valex)
   FNnc=trim(FN)//'#S2'	    ;     deallocate(S2)   ; call uload(FNnc,S2,valex)
   FNnc=trim(FN)//'#T1'	    ;     deallocate(T1)   ; call uload(FNnc,T1,valex)
   FNnc=trim(FN)//'#T2'	    ;     deallocate(T2)   ; call uload(Fnnc,T2,valex)
write(*,*) "T2(180,50,10)=",T2(180,50,10)
   FNnc=trim(FN)//'#P0'	    ;     deallocate(P0)   ; call uload(FNnc,P0,valex)
   FNnc=trim(FN)//'#P'	    ;     deallocate(P)    ; call uload(FNnc,P,valex)
   FNnc=trim(FN)//'#ULF'    ;     deallocate(ULF)  ; call uload(FNnc,ULF,valex)
   FNnc=trim(FN)//'#VLF'    ;     deallocate(VLF)  ; call uload(FNnc,VLF,valex)
   FNnc=trim(FN)//'#SLF'    ;     deallocate(SLF)  ; call uload(FNnc,SLF,valex)
   FNnc=trim(FN)//'#TLF'    ;     deallocate(TLF)  ; call uload(FNnc,TLF,valex)
   FNnc=trim(FN)//'#U'	    ;     deallocate(U)    ; call uload(FNnc,U,valex)
   FNnc=trim(FN)//'#V'	    ;     deallocate(V)    ; call uload(FNnc,V,valex)
   FNnc=trim(FN)//'#W'	    ;     deallocate(W)    ; call uload(FNnc,W,valex)
   FNnc=trim(FN)//'#QAVG'   ;     deallocate(QAVG) ; call uload(FNnc,QAVG,valex)
   FNnc=trim(FN)//'#TNUDGE' ;     deallocate(TNUDGE);call uload(FNnc,TNUDGE,valex)
   FNnc=trim(FN)//'#WAVG'   ;     deallocate(WAVG) ; call uload(FNnc,WAVG,valex)
   FNnc=trim(FN)//'#EVAP'   ;     call uload(FNnc,EVAP,valex)
   FNnc=trim(FN)//'#QSUM'   ;     call uload(FNnc,QSUM,valex)
   FNnc=trim(FN)//'#OLDE'   ;     call uload(FNnc,OLDE,valex)
   FNnc=trim(FN)//'#OLDQ'   ;     call uload(FNnc,OLDQ,valex)
   FNnc=trim(FN)//'#PBAR'   ;     deallocate(PBAR) ; call uload(FNnc,PBAR,valex)
   FNnc=trim(FN)//'#PVAR'   ;     deallocate(PVAR) ; call uload(FNnc,PVAR,valex)
   FNnc=trim(FN)//'#XBAR'   ;     deallocate(XBAR) ; call uload(FNnc,XBAR,valex)
   FNnc=trim(FN)//'#SBAR'   ;     deallocate(SBAR) ; call uload(FNnc,SBAR,valex)
   FNnc=trim(FN)//'#TBAR'   ;     deallocate(TBAR) ; call uload(FNnc,TBAR,valex)
   FNnc=trim(FN)//'#UCLI'   ;     deallocate(UCLI) ; call uload(FNnc,UCLI,valex)
   FNnc=trim(FN)//'#VCLI'   ;     deallocate(VCLI) ; call uload(FNnc,VCLI,valex)
   FNnc=trim(FN)//'#RMSV'   ;     deallocate(RMSV) ; call uload(FNnc,RMSV,valex)
   FNnc=trim(FN)//'#SSURFM' ;     call uload(FNnc,tempreal,valex)
            ssurfm=tempreal ;     deallocate(tempreal)
   FNnc=trim(FN)//'#TSURFM' ;     call uload(FNnc,tempreal,valex)
            tsurfm=tempreal ;     deallocate(tempreal)

write(*,*) "SSURFM(5)=",SSURFM(5)
	

#ifdef _DEBUG_
  print *, "MED pbar", maxval(pbar),minval(pbar),maxloc(pbar),minloc(pbar),sum(pbar)
  print *, "MED pvar", maxval(pvar),minval(pvar),maxloc(pvar),minloc(pvar),sum(pvar)
  print *, "MED xbar", maxval(xbar),minval(xbar),maxloc(xbar),minloc(xbar),sum(xbar)
  print *, "MED sbar", maxval(sbar),minval(sbar),maxloc(sbar),minloc(sbar),sum(sbar)
  print *, "MED tbar", maxval(tbar),minval(tbar),maxloc(tbar),minloc(tbar),sum(tbar)
  print *, "MED ucli", maxval(ucli),minval(ucli),maxloc(ucli),minloc(ucli),sum(ucli)
  print *, "MED vcli", maxval(vcli),minval(vcli),maxloc(vcli),minloc(vcli),sum(vcli)
  print *, "MED rmsv", maxval(rmsv),minval(rmsv),maxloc(rmsv),minloc(rmsv),sum(rmsv)
  print *, "MED ssurfm", maxval(ssurfm),minval(ssurfm),maxloc(ssurfm),minloc(ssurfm),sum(ssurfm)
  print *, "MED tsurfm", maxval(tsurfm),minval(tsurfm),maxloc(tsurfm),minloc(tsurfm),sum(tsurfm)
  print *, "MED u1", maxval(u1),minval(u1),maxloc(u1),minloc(u1),sum(u1)
  print *, "MED u2", maxval(u2),minval(u2),maxloc(u2),minloc(u2),sum(u2)
  print *, "MED v1", maxval(v1),minval(v1),maxloc(v1),minloc(v1),sum(v1)
  print *, "MED v2", maxval(v2),minval(v2),maxloc(v2),minloc(v2),sum(v2)
  print *, "MED s1", maxval(s1),minval(s1),maxloc(s1),minloc(s1),sum(s1)
  print *, "MED s2", maxval(s2),minval(s2),maxloc(s2),minloc(s2),sum(s2)
  print *, "MED t1", maxval(t1),minval(t1),maxloc(t1),minloc(t1),sum(t1)
  print *, "MED t2", maxval(t2),minval(t2),maxloc(t2),minloc(t2),sum(t2)
  print *, "MED u", maxval(u),minval(u),maxloc(u),minloc(u),sum(u)
  print *, "MED v", maxval(v),minval(v),maxloc(v),minloc(v),sum(v)
  print *, "MED w", maxval(w),minval(w),maxloc(w),minloc(w),sum(w)
  print *, "MED p0", maxval(p0),minval(p0),maxloc(p0),minloc(p0),sum(p0)
  print *, "MED p", maxval(p),minval(p),maxloc(p),minloc(p),sum(p)
  print *, "MED rho", maxval(rho),minval(rho),maxloc(rho),minloc(rho),sum(rho)
  print *, "MED x", maxval(x),minval(x),maxloc(x),minloc(x),sum(x)
#endif

! all grid coupling data
      DO K=1,K1
        DO J=1,J0
          DO I=1,I0
            ULF(I,J,K)=U2(I,J,K)
            VLF(I,J,K)=V2(I,J,K)
            SLF(I,J,K)=S2(I,J,K)
            TLF(I,J,K)=T2(I,J,K)
          END DO
        END DO
      END DO
		
 print *, "End of reading !"
 print *, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

END SUBROUTINE RD_RESTART

END MODULE RESTART_GLOBAL
