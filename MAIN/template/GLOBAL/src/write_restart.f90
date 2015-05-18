SUBROUTINE WRITE_RESTART(SIG)
INTEGER :: SIG
INTEGER :: BASE

  SELECT CASE (SIG)
  CASE(0)
  !===================================
  ! Re-arrange the checkpoint timestep
  !===================================
     IF ( CPSAV < 0 ) THEN
       BASE=ABS(CPSAV)*INT(DAODT)
     ELSEIF ( CPSAV > 0 ) THEN
       BASE=CPSAV
     ELSE
       RETURN
     ENDIF
     IF ( MOD(ITF,BASE) == 0 .OR. ITF == MXIT ) THEN
       CALL WR_RESTART
       CHKPTSIG=1
     ELSE
       RETURN
     ENDIF
  CASE(1)
     CALL WR_RESTART
     CHKPTSIG=1
  CASE DEFAULT
     WRITE(6,*)"WRONG WRITE_RESTART SIGNAL"
     WRITE(FNO_LOG,*)"WRONG WRITE_RESTART SIGNAL"
     RETURN
  END SELECT

      WRITE(6,4040) ITF
      WRITE(6,1406) TRIM(CASE_NAME),N360,NYR
      WRITE(FNO_LOG,4040) ITF
      WRITE(FNO_LOG,1406) TRIM(CASE_NAME),N360,NYR
4040  FORMAT('Restart file written at execution STEPS = ',I6)
1406  FORMAT(A,' save at day',I4,', year',I3)
END SUBROUTINE WRITE_RESTART

