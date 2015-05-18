SUBROUTINE PMOVE(SIG)
  INTEGER :: SIG,NEWDAY,NEWDAY_SIG=0
  LOGICAL :: alive
  INTEGER :: BASE
  CHARACTER(LEN=256) :: FNAME

  IF ( FL_PMOVE_ON == 0 ) RETURN
  FNAME=TRIM(WORKDIR)//"/"//TRIM(IDIR)//"/newtracer"
  INQUIRE(FILE=TRIM(FNAME),EXIST=alive)
  IF ( alive ) THEN
    OPEN(43,file=TRIM(FNAME))
    READ(43,"(I7)",IOSTAT=istat) NEWDAY
    CLOSE(43,STATUS="DELETE")
    IF ( istat /= 0 ) THEN
       WRITE(6,*)"I got wrong newtracer signal!"
       WRITE(6,*)"The value of new tracer signal must be"
       WRITE(6,*)" <value> > 0 : re-new the FL_PMOVE_NEWDAY value to every FL_PMOVE_NEWDAY timesteps re-put the tracer"
       WRITE(6,*)" <value> < 0 : re-new the FL_PMOVE_NEWDAY value to every abs(FL_PMOVE_NEWDAY) days re-put the tracer"
       WRITE(6,*)" <value> = 0 : new tracers need to be re-put"
       WRITE(6,*)""
       WRITE(6,*)"I'll continue to run the job with the original FL_PMOVE_NEWDAY value"
    ENDIF
    IF ( NEWDAY == 0 ) THEN
       NEWDAY_SIG=-1
    ELSE
        FL_PMOVE_NEWDAY=NEWDAY
    ENDIF 
  ENDIF
  SELECT CASE (SIG)
  CASE(0)
    CALL PMOVE_INIT
  CASE(1)
    ! NEWDAY_SIG = 0 : timestep is still in old day
    ! NEWDAY_SIG = 1 : timestep is update to new day, it might be output the tracer location
    ! NEWDAY_SIG =-1 : new tracers need to be re-put
    IF ( .NOT. (NEWDAY_SIG == -1) ) THEN
      IF ( FL_PMOVE_NEWDAY < 0 ) THEN
        BASE=ABS(FL_PMOVE_NEWDAY)*INT(DAODT)
      ELSEIF ( FL_PMOVE_NEWDAY == 1) THEN
        BASE=FL_PMOVE_NEWDAY
      ELSE
        RETURN
      ENDIF
      IF ( MOD(ITF,BASE) == 0 ) THEN
        NEWDAY_SIG=1
      ELSE
        NEWDAY_SIG=0
      ENDIF
    ENDIF
    CALL PMOVE_CHECK(NEWDAY_SIG)
    CALL PMOVE_RUN(1)
  CASE(2)
    CALL PMOVE_FINALIZE
  CASE DEFAULT
    WRITE(6,*)"WRONG PMOVE SIGNAL"
    WRITE(FNO_LOG,*)"WRONG PMOVE SIGNAL"
  END SELECT
END SUBROUTINE PMOVE

