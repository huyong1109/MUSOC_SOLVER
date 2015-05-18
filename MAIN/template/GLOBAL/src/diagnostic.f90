FUNCTION DIAG_BLOWUP()
!
! diagnostic with main procedure 
! This file is included in main.f90
!
LOGICAL DIAG_BLOWUP

! ==============
! Monitor solution progress in detail and check for possible instability
! NOTE: DIECAST has been found to ALWAYS BE STABLE unless DT is too large
! or user-specified i.c.s or b.c.s are UNPHYSICAL)
  CALL RANGER(V2,IN,I0,2,2,I1,J1,ILO,JLO,IHI,JHI,VLO,VHI)
  TEMP=DAYS/360.

  IF ( (ITF_GENERAL .LT. TOPTS) .OR. (MOD(ITF_GENERAL,TOPTS) .EQ. 0) .OR. (ITFTO .LE. NITFTO) ) THEN
  !IF ((ITF-IT0) <= ITFDAY .OR. MOD(ITF,5*ITFDAY) == 0) THEN
   IF (ITFTO .LE. NITFTO ) THEN
   ! IF ( (ITF_GENERAL .LT. TOPTS) .OR. (MOD(ITF_GENERAL,TOPTS) .EQ. 0) .OR. (ITFTO .LE. NITFTO) ) THEN
   !IF (MOD(ITF,5*ITFDAY).EQ.0) THEN
     J=J0/8
     WRITE(FNO_LOG,1401) (U2(I,J,1),I=I0-3,I0), (U2(I,J,1),I=3,6)
   ENDIF
    TMP=0.
    DO J=2,J2
    DO I=2,I2
      TMP=MAX(TMP,(1.-IU0(I,J))*ABS(U(I,J,1)),(1.-IV0(I,J))*ABS(V(I,J,1)))
    ENDDO
    ENDDO
    WRITE(6,1402) ITF,TEMP,ILO,JLO,VLO,IHI,JHI,VHI,TMP
    WRITE(FNO_LOG,1402) ITF,TEMP,ILO,JLO,VLO,IHI,JHI,VHI,TMP
    CALL FLUSH(FNO_LOG)
  ENDIF

   IF ((VHI-VLO) >= 2000.) THEN
     WRITE(6,1402) ITF,TEMP,ILO,JLO,VLO,IHI,JHI,VHI,TMP
     WRITE(6,1403) ITF
     WRITE(FNO_LOG,1402) ITF,TEMP,ILO,JLO,VLO,IHI,JHI,VHI,TMP
     WRITE(FNO_LOG,1403) ITF
       DIAG_BLOWUP=.true.
   ELSE
       DIAG_BLOWUP=.false. 
   ENDIF
1401  FORMAT('periodic conditions check '/(8F7.2))
1402  FORMAT('Blowup-diagnostic at itf',I8,', yr',F7.2,/'(ILO,JLO,VLO,IHI,JHI,VHI)=(',&
             2I4,1X,F9.4,2I4,1X,F9.4,')'/'max vel on land=',F8.4/)
1403  FORMAT('STOP. UNIX-pectedly/ large velocity at ITF=',I7)

END FUNCTION DIAG_BLOWUP

SUBROUTINE DIAG_CROSS
!cross-equatorial exchange
      WRITE(6,1404) INT(DAYS)
      WRITE(FNO_LOG,1404) INT(DAYS)
        TMP1=0.
        TMP2=0.
        J=112
        DO K=1,K1
          TMP=1.E-12*DXV(J)/ODZ(K)
        ! this range must be modified
          DO I=301,345
            TMP1=TMP1+TMP*MAX(V(I,J,K),0.)
            TMP2=TMP2-TMP*MIN(V(I,J,K),0.)
          ENDDO
          WRITE(6,1405) TMP1,TMP2,TMP1-TMP2,.01*Z(2*K+1)
          WRITE(FNO_LOG,1405) TMP1,TMP2,TMP1-TMP2,.01*Z(2*K+1)
        ENDDO
1404  FORMAT('at day',I5/'NAB summed southern boundary source (Sv)',6X,'to depth (m)')
1405  FORMAT('in=',F5.1,', out=',F5.1,', net=',F5.1,8X,F7.2)
END SUBROUTINE DIAG_CROSS

