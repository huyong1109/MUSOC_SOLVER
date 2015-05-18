SUBROUTINE AVERAGE

! ====================================================
! This is exactly at southern boundary of MEDiNA model
! ====================================================
      DO K=1,K1
      DO I=2,I1
        VGLO(I,K)=V(I,J2,K)
        UAVG(I,K)=U2(I,J2,K)
        VAVG(I,K)=V2(I,J2,K)
        SAVG(I,K)=S2(I,J2,K)
        TAVG(I,K)=T2(I,J2,K)
      ENDDO
      ENDDO

! =============
! time averages
! =============
      IF (MOD(ITF,ITFDAY).EQ.1) THEN
        AV=AV+1.
        WT=1./AV
        OWT=1.-WT
        DO J=2,J1
        DO I=2,I1
      ! Mean surface pressure
          PBAR(I,J)=OWT*PBAR(I,J)+WT*P(I,J,1)
      ! Mean squared deviation of surface pressure from mean
          PVAR(I,J)=OWT*PVAR(I,J)+WT*(P(I,J,1)-PBAR(I,J))**2
        ENDDO
        ENDDO
      ! Integrate velocity to get barotropic streamfunction and
      ! evaluate at natural psi locations (grid line intersects)
        SCR(:,:,1:2)=0.
        DO K=1,K1
          TMP=1.E-12/ODZ(K)
          DO I=2,I1
            SCR(I,1,1)=SCR(I-1,1,1)+V(I,1,K)*TMP*DXV(1)
          ENDDO
          DO J=2,J1
          DO I=1,I1
            SCR(I,J,2)=SCR(I,J,2)+U(I,J,K)*TMP*DY(J)
          ENDDO
          ENDDO
        ENDDO
        DO I=1,I1
        DO J=2,J1
          SCR(I,J,1)=SCR(I,J-1,1)-SCR(I,J,2)
        ENDDO
        ENDDO
      ! Interpolate streamfunction to "P" cell centers
        DO J=1,J2
        DO I=1,I2
          SCR(I+1,J+1,2)=0.25*(SCR(I,J,1)+SCR(I+1,J,1)+SCR(I,J+1,1)+SCR(I+1,J+1,1))
        ENDDO
        ENDDO
      ! vertically averaged vel (cell centered: not for xi integration)
        DO K=1,K1
          TMP=1./ODZ(K)
          DO J=2,J1
          DO I=2,I1
            SCR(I,J,1)=SCR(I,J,1)+U2(I,J,K)*TMP
            SCR(I,J,2)=SCR(I,J,2)+V2(I,J,K)*TMP
          ENDDO
          ENDDO
        ENDDO
        DO J=2,J1
        DO I=2,I1
          TMP=1./Z(2*KB(I,J)+1)
          SCR(I,J,1)=TMP*SCR(I,J,1)
          SCR(I,J,2)=TMP*SCR(I,J,2)
        ENDDO
        ENDDO
      ! time mean streamfunction
        DO I=2,I1
        DO J=2,J1
          XBAR(I,J)=OWT*XBAR(I,J)+WT*SCR(I,J,2)
        ENDDO
        ENDDO
        IF(K1>=3)THEN
      ! 700m depth rms eddy velocity
        K=0
        DO WHILE(.TRUE.)
          K=K+1
          IF (Z(2*K).LT.7.E4) CYCLE
          EXIT
        ENDDO
        KL1=K-1
        KL2=K
        WT2=(7.E4-Z(2*KL1))/(Z(2*KL2)-Z(2*KL1))
        WT1=1.-WT2
        DO J=2,J1
        DO I=2,I1
          TMP= WT1*U2(I,J,KL1)+WT2*U2(I,J,KL2)
          TEMP=WT1*V2(I,J,KL1)+WT2*V2(I,J,KL2)
        ! Climate mean U and V at depth 700 meters
          UCLI(I,J)=OWT*UCLI(I,J)+WT*TMP
          VCLI(I,J)=OWT*VCLI(I,J)+WT*TEMP
          RMSV(I,J)=OWT*RMSV(I,J)+WT*((UCLI(I,J)-TMP)**2+(VCLI(I,J)-TMP)**2)
        ENDDO
        ENDDO
        ENDIF
      ! time mean salinity, temperature and longitudinally integrated
      ! latitudinal volume transport (latitudnal velocity * area)
        DO K=1,K1
          DO J=2,J1
          DO I=2,I1
            SBAR(I,J,K)=OWT*SBAR(I,J,K)+WT*S2(I,J,K)
          ENDDO
          ENDDO
          DO J=2,J1
          DO I=2,I1
            TBAR(I,J,K)=OWT*TBAR(I,J,K)+WT*T2(I,J,K)
          ENDDO
          ENDDO
        ENDDO
      ENDIF


END SUBROUTINE AVERAGE
