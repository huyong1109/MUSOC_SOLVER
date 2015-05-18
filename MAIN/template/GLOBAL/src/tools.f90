
! ----------------------------------------------------------------------
!Nelson notice - find max and min value and their locations
!              - IN is land information
SUBROUTINE RANGER(FLD,IN,IR,IL,JL,IH,JH,ILO,JLO,IHI,JHI,FMIN,FMAX)
! ----------------------------------------------------------------------
INTEGER(2) IN(IR,*)
REAL FLD(IR,*)
  FMIN=1.E20
  FMAX=-1.E20
  DO J=JL,JH
  DO I=IL,IH
    FMIN=(1-IN(I,J))*FMIN+IN(I,J)*MIN(FMIN,FLD(I,J))
    FMAX=(1-IN(I,J))*FMAX+IN(I,J)*MAX(FMAX,FLD(I,J))
  ENDDO
  ENDDO
  DO J=JL,JH
  DO I=IL,IH
    IF (FMAX.EQ.FLD(I,J)) THEN
      IHI=I
      JHI=J
    ELSE
      IF (FMIN.EQ.FLD(I,J)) THEN 
        ILO=I
        JLO=J
      ENDIF
    ENDIF
  ENDDO
  ENDDO
END SUBROUTINE RANGER

