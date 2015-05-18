SUBROUTINE SNAPSHOT_OUTPUT(SIG)
  INTEGER :: SIG
  INTEGER :: BASE,NN,MM
  CHARACTER(LEN=160) :: FILENAME,DIR
  DIR=TRIM(WORKDIR)//"/"//TRIM(ODIR)//"/"
! Compute output quotient base
  IF ( ISAV < 0 ) THEN
    BASE=ABS(ISAV)*INT(DAODT)
  ELSEIF ( ISAV > 0 ) THEN
    BASE=ISAV
  ELSE
    RETURN
  ENDIF
  SELECT CASE (SIG)
  CASE(0)
    IF (MOD(ITF,BASE)==0) THEN
      NN=ITF/ITFDAY
      MM=MOD(ITF,ITFDAY)
      IF (MM == 0) THEN
        IF ( ISAV < 0 ) THEN
          WRITE(FILENAME,"(A,A,i3.3,A,i5.5,A)")TRIM(DIR),"member",member,"_T1_",NN,".nc"
          CALL TIMCOM_NC_OUTPUT_TYPE1(FILENAME)
          WRITE(FILENAME,"(A,A,i3.3,A,i5.5,A)")TRIM(DIR),"member",member,"_T2_",NN,".nc"
          CALL TIMCOM_NC_OUTPUT_TYPE2(FILENAME)
          WRITE(FILENAME,"(A,A,i3.3,A,i5.5,A)")TRIM(DIR),"member",member,"_fluxes_",NN,".nc"
          CALL TIMCOM_NC_OUTPUT_FLUXES(FILENAME)
        ELSEIF ( ISAV > 0 ) THEN
          WRITE(FILENAME,"(A,A,i3.3,A,i5.5,A)")TRIM(DIR),"member",member,"_T1_",NN,"-00000.nc"
          CALL TIMCOM_NC_OUTPUT_TYPE1(FILENAME)
          WRITE(FILENAME,"(A,A,i3.3,A,i5.5,A)")TRIM(DIR),"member",member,"_T2_",NN,"-00000.nc"
          CALL TIMCOM_NC_OUTPUT_TYPE2(FILENAME)
          WRITE(FILENAME,"(A,A,i3.3,A,i5.5,A)")TRIM(DIR),"member",member,"_fluxes_",NN,"-00000.nc"
          CALL TIMCOM_NC_OUTPUT_FLUXES(FILENAME)
        ELSE
          RETURN
        ENDIF
      ELSE
        WRITE(FILENAME,"(A,A,i3.3,A,i5.5,A,i5.5,A)")TRIM(DIR),"member",member,"_T1_",NN,"-",MM,".nc"
        CALL TIMCOM_NC_OUTPUT_TYPE1(FILENAME)
        WRITE(FILENAME,"(A,A,i3.3,A,i5.5,A,i5.5,A)")TRIM(DIR),"member",member,"_T2_",NN,"-",MM,".nc"
        CALL TIMCOM_NC_OUTPUT_TYPE2(FILENAME)
        WRITE(FILENAME,"(A,A,i3.3,A,i5.5,A,i5.5,A)")TRIM(DIR),"member",member,"_fluxes_",NN,"-",MM,".nc"
        CALL TIMCOM_NC_OUTPUT_FLUXES(FILENAME)
      ENDIF
    ENDIF
  CASE(1)
    NN=ITF/ITFDAY
    MM=MOD(ITF,ITFDAY)
    IF (MM == 0) THEN
      IF ( ISAV < 0 ) THEN
        WRITE(FILENAME,"(A,A,i3.3,A,i5.5,A)")TRIM(DIR),"member",member,"_T1_",NN,".nc"
        CALL TIMCOM_NC_OUTPUT_TYPE1(FILENAME)
        WRITE(FILENAME,"(A,A,i3.3,A,i5.5,A)")TRIM(DIR),"member",member,"_T2_",NN,".nc"
        CALL TIMCOM_NC_OUTPUT_TYPE2(FILENAME)
        WRITE(FILENAME,"(A,A,i3.3,A,i5.5,A)")TRIM(DIR),"member",member,"_fluxes_",NN,".nc"
        CALL TIMCOM_NC_OUTPUT_FLUXES(FILENAME)
      ELSEIF ( ISAV > 0 ) THEN
        WRITE(FILENAME,"(A,A,i3.3,A,i5.5,A)")TRIM(DIR),"member",member,"_T1_",NN,"-00000.nc"
        CALL TIMCOM_NC_OUTPUT_TYPE1(FILENAME)
        WRITE(FILENAME,"(A,A,i3.3,A,i5.5,A)")TRIM(DIR),"member",member,"_T2_",NN,"-00000.nc"
        CALL TIMCOM_NC_OUTPUT_TYPE2(FILENAME)
        WRITE(FILENAME,"(A,A,i3.3,A,i5.5,A)")TRIM(DIR),"member",member,"_fluxes_",NN,"-00000.nc"
        CALL TIMCOM_NC_OUTPUT_FLUXES(FILENAME)
      ELSE
        RETURN
      ENDIF
    ELSE 
      WRITE(FILENAME,"(A,A,i3.3,A,i5.5,A,i5.5,A)")TRIM(DIR),"member",member,"_T1_",NN,"-",MM,".nc"
      CALL TIMCOM_NC_OUTPUT_TYPE1(FILENAME)
      WRITE(FILENAME,"(A,A,i3.3,A,i5.5,A,i5.5,A)")TRIM(DIR),"member",member,"_T2_",NN,"-",MM,".nc"
      CALL TIMCOM_NC_OUTPUT_TYPE2(FILENAME)
      WRITE(FILENAME,"(A,A,i3.3,A,i5.5,A,i5.5,A)")TRIM(DIR),"member",member,"_fluxes_",NN,"-",MM,".nc"
      CALL TIMCOM_NC_OUTPUT_FLUXES(FILENAME)
    ENDIF
  CASE DEFAULT
      RETURN
  END SELECT
END SUBROUTINE SNAPSHOT_OUTPUT

