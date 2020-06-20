C
C********* ********* ********* ********* ********* ********* *********C
C
C     timedevelop.f
C
C     Program for homogeneous shear turbulent flow by K-e model
C     MADE BY MASAYOSHI OKAMOTO
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON
     & /BAS/TK,TE,TAU
     & /RHS/FTK1,FTE1,FTA1,FTK2,FTE2,FTA2,FTK3,FTE3,FTA3,FTK4,FTE4
      COMMON
     & /ERR/ERRTK,ERRTE
      COMMON
     & /CONMODEL/CE1,CE2,CNU
     & /CON1/DT,SH
     & /CON2/ITYPE,IDEG,ITINC,NOUT
C
      CALL INIT
      OPEN(12,FILE='TIME-DATA.txt',FORM='FORMATTED')
      WRITE(12,*) '       Numerical result     Exact solution'
      WRITE(12,*) 'time    k     e    tau     k      e    tau'
      IF(ITYPE.EQ.0) THEN
C
C     EXPLICIT METHOD
C
        DO 100 IT=1,ITINC
          CALL ADBS(IT)
          IF (MOD(IT,NOUT).EQ.0) THEN
            CALL ENPR(IT)
          ENDIF
  100   CONTINUE
      ELSEIF(ITYPE.EQ.1) THEN
C
C     EXPLICIT METHOD OF RUNGE-KUTTA
C
        DO 110 IT=1,ITINC
          CALL RNKT
          IF (MOD(IT,NOUT).EQ.0) THEN
            CALL ENPR(IT)
          ENDIF
  110   CONTINUE
      ELSE
C
C     IMPLICIT METHOD
C
        DO 120 IT=1,ITINC
          CALL IMCK(IT)
          IF (MOD(IT,NOUT).EQ.0) THEN
            CALL ENPR(IT)
          ENDIF
  120   CONTINUE
      ENDIF
      CLOSE(12)
      WRITE(6,*) 'ERROR OF TURBULENCE ENERGY'
      WRITE(6,*) ERRTK
      WRITE(6,*) 'ERROR OF DISSIPATION RATE'
      WRITE(6,*) ERRTE
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE INIT
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON
     & /BAS/TK,TE,TAU
     & /RHS/FTK1,FTE1,FTA1,FTK2,FTE2,FTA2,FTK3,FTE3,FTA3,FTK4,FTE4
      COMMON
     & /ERR/ERRTK,ERRTE
      COMMON
     & /CONMODEL/CE1,CE2,CNU
     & /CON1/DT,SH
     & /CON2/ITYPE,IDEG,ITINC,NOUT
C
      TK=1.D0
      TE=1.D0
      TAU=1.D0
      FTK1=0.D0
      FTE1=0.D0
      FTA1=0.D0
      FTK2=0.D0
      FTE2=0.D0
      FTA2=0.D0
      FTK3=0.D0
      FTE3=0.D0
      FTA3=0.D0
      FTK4=0.D0
      FTE4=0.D0
      CE1=1.44D0
      CE2=1.94D0
      CNU=0.09D0
      ERRTK=0.D0
      ERRTE=0.D0
      WRITE(6,*) '0: EXPLICIT METHOD OF ADAMS-BASHFORTH'
      WRITE(6,*) '1: EXPLICIT METHOD OF RUNGE-KUTTA'
      WRITE(6,*) '2: IMPLICIT METHOD 0F ADAMS-MOULTON'
      READ(5,*) ITYPE
      IF(ITYPE.EQ.0) THEN
        WRITE(6,*) '0: EXPLICIT EULER METHOD'
        WRITE(6,*) '1: 2ND ADAMS-BASHFORTH METHOD'
        WRITE(6,*) '2: 3RD ADAMS-BASHFORTH METHOD'
        WRITE(6,*) '3: 4TH ADAMS-BASHFORTH METHOD'
      ELSEIF(ITYPE.EQ.1) THEN
        WRITE(6,*) '0: MODIFIED EULER METHOD'
        WRITE(6,*) '1: IMPROVED EULER (HEUN) METHOD'
        WRITE(6,*) '2: 3RD RUNGE-KUTTA METHOD'
        WRITE(6,*) '3: 4TH RUNGE-KUTTA METHOD'
      ELSEIF(ITYPE.EQ.2) THEN
        WRITE(6,*) '0: IMPLICIT EULER METHOD'
        WRITE(6,*) '1: CRANK NICOLSON METHOD'
        WRITE(6,*) '2: 3RD ADAMS-MOULTON METHOD'
        WRITE(6,*) '3: 4TH ADAMS-MOULTON METHOD'
      ENDIF
      READ(5,*) IDEG
      WRITE(6,*) 'SHEAR RATE'
      READ(5,*) SH
      WRITE(6,*) 'TIME INTERVAL'
      READ(5,*) DT
      WRITE(6,*) 'TIME STEP NUMBER'
      READ(5,*) ITINC
      WRITE(6,*) 'TOTAL TIME'
      WRITE(6,*) DT*DFLOAT(ITINC)
      WRITE(6,*) 'OUTPUT INTERVAL NUMBER'
      READ(5,*) NOUT
      WRITE(6,*) 'OUTPUT INTERVAL TIME'
      WRITE(6,*) DT*DFLOAT(NOUT)
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE ADBS(IT)
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON
     & /BAS/TK,TE,TAU
     & /RHS/FTK1,FTE1,FTA1,FTK2,FTE2,FTA2,FTK3,FTE3,FTA3,FTK4,FTE4
      COMMON
     & /CONMODEL/CE1,CE2,CNU
     & /CON1/DT,SH
     & /CON2/ITYPE,IDEG,ITINC,NOUT
C
      IF(IDEG.EQ.0) THEN
        S1=1.D0
        S2=0.D0
        S3=0.D0
        S4=0.D0
      ELSEIF(IDEG.EQ.1) THEN
        S1=1.5D0
        S2=-0.5D0
        S3=0.D0
        S4=0.D0
        IF(IT.EQ.1) THEN
          CALL MKSTAB(-DT,FTK2,FTE2)
        ENDIF
      ELSEIF(IDEG.EQ.2) THEN
        S1=23.D0/12.D0
        S2=-4.D0/3.D0
        S3=5.D0/12.D0
        S4=0.D0
        IF(IT.EQ.1) THEN
          CALL MKSTAB(-DT,FTK2,FTE2)
          CALL MKSTAB(-2.D0*DT,FTK3,FTE3)
        ENDIF
      ELSEIF(IDEG.EQ.3) THEN
        S1=55.D0/24.D0
        S2=-59.D0/24.D0
        S3=37.D0/24.D0
        S4=-3.D0/8.D0
        IF(IT.EQ.1) THEN
          CALL MKSTAB(-DT,FTK2,FTE2)
          CALL MKSTAB(-2.D0*DT,FTK3,FTE3)
          CALL MKSTAB(-3.D0*DT,FTK4,FTE4)
        ENDIF
      ENDIF
      FTK1=CNU*SH**2*TK**2/TE-TE
      FTE1=CE1*CNU*SH**2*TK-CE2*TE**2/TK
      TK=TK+DT*(S1*FTK1+S2*FTK2+S3*FTK3+S4*FTK4)
      TE=TE+DT*(S1*FTE1+S2*FTE2+S3*FTE3+S4*FTE4)
      FTK4=FTK3
      FTE4=FTE3
      FTK3=FTK2
      FTE3=FTE2
      FTK2=FTK1
      FTE2=FTE1
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE RNKT
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON
     & /BAS/TK,TE,TAU
     & /RHS/FTK1,FTE1,FTA1,FTK2,FTE2,FTA2,FTK3,FTE3,FTA3,FTK4,FTE4
      COMMON
     & /CONMODEL/CE1,CE2,CNU
     & /CON1/DT,SH
     & /CON2/ITYPE,IDEG,ITINC,NOUT
C
      IF(IDEG.EQ.0) THEN
        FTK1=CNU*SH**2*TK**2/TE-TE
        FTE1=CE1*CNU*SH**2*TK-CE2*TE**2/TK
        TKA=TK+0.5D0*DT*FTK1
        TEA=TE+0.5D0*DT*FTE1
        FTK2=CNU*SH**2*TKA**2/TEA-TEA
        FTE2=CE1*CNU*SH**2*TKA-CE2*TEA**2/TKA
        TK=TK+DT*FTK2
        TE=TE+DT*FTE2
      ELSEIF(IDEG.EQ.1) THEN
        FTK1=CNU*SH**2*TK**2/TE-TE
        FTE1=CE1*CNU*SH**2*TK-CE2*TE**2/TK
        TKA=TK+DT*FTK1
        TEA=TE+DT*FTE1
        FTK2=CNU*SH**2*TKA**2/TEA-TEA
        FTE2=CE1*CNU*SH**2*TKA-CE2*TEA**2/TKA
        TK=TK+0.5D0*DT*(FTK1+FTK2)
        TE=TE+0.5D0*DT*(FTE1+FTE2)
      ELSEIF(IDEG.EQ.2) THEN
        FTK1=CNU*SH**2*TK**2/TE-TE
        FTE1=CE1*CNU*SH**2*TK-CE2*TE**2/TK
        TKA=TK+0.5D0*DT*FTK1
        TEA=TE+0.5D0*DT*FTE1
        FTK2=CNU*SH**2*TKA**2/TEA-TEA
        FTE2=CE1*CNU*SH**2*TKA-CE2*TEA**2/TKA
        TKA=TK-DT*FTK1+2.D0*DT*FTK2
        TEA=TE-DT*FTE1+2.D0*DT*FTE2
        FTK3=CNU*SH**2*TKA**2/TEA-TEA
        FTE3=CE1*CNU*SH**2*TKA-CE2*TEA**2/TKA
        TK=TK+DT*(FTK1+4.D0*FTK2+FTK3)/6.D0
        TE=TE+DT*(FTE1+4.D0*FTE2+FTE3)/6.D0
      ELSEIF(IDEG.EQ.3) THEN
        FTK1=CNU*SH**2*TK**2/TE-TE
        FTE1=CE1*CNU*SH**2*TK-CE2*TE**2/TK
        TKA=TK+0.5D0*DT*FTK1
        TEA=TE+0.5D0*DT*FTE1
        FTK2=CNU*SH**2*TKA**2/TEA-TEA
        FTE2=CE1*CNU*SH**2*TKA-CE2*TEA**2/TKA
        TKA=TK+0.5D0*DT*FTK2
        TEA=TE+0.5D0*DT*FTE2
        FTK3=CNU*SH**2*TKA**2/TEA-TEA
        FTE3=CE1*CNU*SH**2*TKA-CE2*TEA**2/TKA
        TKA=TK+DT*FTK3
        TEA=TE+DT*FTE3
        FTK4=CNU*SH**2*TKA**2/TEA-TEA
        FTE4=CE1*CNU*SH**2*TKA-CE2*TEA**2/TKA
        TK=TK+DT*(FTK1+2.D0*FTK2+2.D0*FTK3+FTK4)/6.D0
        TE=TE+DT*(FTE1+2.D0*FTE2+2.D0*FTE3+FTE4)/6.D0
      ENDIF
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE IMCK(IT)
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON
     & /BAS/TK,TE,TAU
     & /RHS/FTK1,FTE1,FTA1,FTK2,FTE2,FTA2,FTK3,FTE3,FTA3,FTK4,FTE4
      COMMON
     & /CONMODEL/CE1,CE2,CNU
     & /CON1/DT,SH
     & /CON2/ITYPE,IDEG,ITINC,NOUT
C
      IF(IDEG.EQ.0) THEN
        S0=1.D0
        S1=0.D0
        S2=0.D0
        S3=0.D0
      ELSEIF(IDEG.EQ.1) THEN
        S0=1.D0/2.D0
        S1=1.D0/2.D0
        S2=0.D0
        S3=0.D0
      ELSEIF(IDEG.EQ.2) THEN
        IF(IT.EQ.1) THEN
          S0=1.D0/2.D0
          S1=1.D0/2.D0
          S2=0.D0
          S3=0.D0
        ELSE
          S0=5.D0/12.D0
          S1=8.D0/12.D0
          S2=-1.D0/12.D0
          S3=0.D0
        ENDIF
      ELSEIF(IDEG.EQ.3) THEN
        IF(IT.EQ.1) THEN
          S0=1.D0/2.D0
          S1=1.D0/2.D0
          S2=0.D0
          S3=0.D0
        ELSEIF(IT.EQ.2) THEN
          S0=5.D0/12.D0
          S1=8.D0/12.D0
          S2=-1.D0/12.D0
          S3=0.D0
        ELSE
          S0=9.D0/24.D0
          S1=19.D0/24.D0
          S2=-5.D0/24.D0
          S3=1.D0/24.D0
        ENDIF
      ENDIF
      TAU=TK/TE
      FTA1=-CNU*(CE1-1.D0)*(TAU*SH)**2+(CE2-1.D0)
      FTK1=CNU*TK*TAU*SH**2-TK/TAU
      TAU=(-1.D0+DSQRT(1.D0
     &     +4.D0*S0*CNU*(CE1-1.D0)*SH**2*DT*(S0*(CE2-1.D0)*DT
     &                 +TAU+DT*(S1*FTA1+S2*FTA2+S3*FTA3))))
     &      /(2.D0*S0*CNU*(CE1-1.D0)*DT*SH**2)
      TK=(TK+DT*(S1*FTK1+S2*FTK2+S3*FTK3))
     &  /(1.D0-S0*DT*(CNU*TAU*SH**2-1.D0/TAU))
      TE=TK/TAU
      FTA3=FTA2
      FTK3=FTK2
      FTA2=FTA1
      FTK2=FTK1
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE ENPR(IT)
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON
     & /BAS/TK,TE,TAU
     & /RHS/FTK1,FTE1,FTA1,FTK2,FTE2,FTA2,FTK3,FTE3,FTA3,FTK4,FTE4
      COMMON
     & /ERR/ERRTK,ERRTE
      COMMON
     & /CONMODEL/CE1,CE2,CNU
     & /CON1/DT,SH
     & /CON2/ITYPE,IDEG,ITINC,NOUT
C
      TIME=DFLOAT(IT)*DT
      TAU=TK/TE
      CALL TASV(TIME,TTK,TTE,TTAU)
      WRITE(12,'( 1H ,1P,7E24.14 )') TIME,TK,TE,TAU,TTK,TTE,TTAU
      ETK=DABS(TK-TTK)/TTK
      ETE=DABS(TE-TTE)/TTE
      ERRTK=DMAX1(ERRTK,ETK)
      ERRTE=DMAX1(ERRTE,ETE)
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE TASV(TIME,TTK,TTE,TTAU)
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON
     & /CONMODEL/CE1,CE2,CNU
     & /CON1/DT,SH
     & /CON2/ITYPE,IDEG,ITINC,NOUT
C
      TK0=1.D0
      TE0=1.D0
      TAU0=1.D0
      A=DSQRT((CE2-1.D0)/(CE1-1.D0)/CNU/SH**2)
      B=2.D0*DSQRT(CNU*(CE1-1.D0)*(CE2-1.D0)*SH**2)
      C=DLOG(DABS((TAU0-A)/(TAU0+A)))
      TTAU=A*(1.D0-DEXP(-B*TIME+C))/(1.D0+DEXP(-B*TIME+C))
      TTK=TK0*DEXP(-TIME/A+(CE2-1.D0)*TIME/A/(CE1-1.D0)
     & +1.D0/(CE1-1.D0)*DLOG((1.D0+DEXP(-B*TIME+C))/(1.D0+DEXP(C)))
     & -1.D0/(CE2-1.D0)*DLOG((1.D0-DEXP(-B*TIME+C))/(1.D0-DEXP(C))))
      TTE=TTK/TTAU
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE MKSTAB(TIME,FK,FE)
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON
     & /CONMODEL/CE1,CE2,CNU
     & /CON1/DT,SH
     & /CON2/ITYPE,IDEG,ITINC,NOUT
C
      TK0=1.D0
      TE0=1.D0
      TAU0=1.D0
      A=DSQRT((CE2-1.D0)/(CE1-1.D0)/CNU/SH**2)
      B=2.D0*DSQRT(CNU*(CE1-1.D0)*(CE2-1.D0)*SH**2)
      C=DLOG(DABS((TAU0-A)/(TAU0+A)))
      TTAU=A*(1.D0-DEXP(-B*TIME+C))/(1.D0+DEXP(-B*TIME+C))
      TTK=TK0*DEXP(-TIME/A+(CE2-1.D0)*TIME/A/(CE1-1.D0)
     & +1.D0/(CE1-1.D0)*DLOG((1.D0+DEXP(-B*TIME+C))/(1.D0+DEXP(C)))
     & -1.D0/(CE2-1.D0)*DLOG((1.D0-DEXP(-B*TIME+C))/(1.D0-DEXP(C))))
      TTE=TTK/TTAU
      FK=CNU*SH**2*TTK**2/TTE-TTE
      FE=CE1*CNU*SH**2*TTK-CE2*TTE**2/TTK
      RETURN
      END
