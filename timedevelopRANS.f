C
C********* ********* ********* ********* ********* ********* *********C
C
C     timedevelopRANS.f
C
C     Program for homogeneous shear turbulent flow by RANS
C     MADE BY MASAYOSHI OKAMOTO
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON
     & /BAS/TK,TE,UU,VV,WW,UV
      COMMON
     & /CONMODEL/CE1,CE2
     & /CON1/DT,SH,ANU
     & /CON2/ITYPE,ITINC,NOUT
C
      CALL INIT
      OPEN(12,FILE='TIME-DATA.txt',FORM='FORMATTED')
      WRITE(12,*) 'time    k      e      uu     vv      ww      uv'
      DO 100 IT=1,ITINC
        IF(ITYPE.EQ.0) THEN
          CALL RKLS
        ELSEIF(ITYPE.EQ.1) THEN
          CALL RKLRR
        ELSEIF(ITYPE.EQ.2) THEN
          CALL RKSSG
        ELSEIF(ITYPE.EQ.3) THEN
          CALL RKASM
        ENDIF
        IF (MOD(IT,NOUT).EQ.0) THEN
          CALL ENPR(IT)
        ENDIF
  100 CONTINUE
      CLOSE(12)
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
     & /BAS/TK,TE,UU,VV,WW,UV
      COMMON
     & /CONMODEL/CE1,CE2
     & /CON1/DT,SH,ANU
     & /CON2/ITYPE,ITINC,NOUT
C
      ANU=0.04D0
      TE=0.33597D0
      UU=0.57507D0
      VV=0.36902D0
      WW=0.5108D0
      TK=0.5D0*(UU+VV+WW)
      UV=-0.22575D0
      WRITE(6,*) '0: PREDICTION OF Launder-Sharma model'
      WRITE(6,*) '1: PREDICTION OF LRR model'
      WRITE(6,*) '2: PREDICTION 0F SSG model'
      WRITE(6,*) '3: PREDICTION 0F ASM model'
      READ(5,*) ITYPE
      IF(ITYPE.EQ.0) THEN
        CE1=1.44D0
        CE2=1.94D0
      ELSEIF(ITYPE.EQ.1) THEN
        CE1=1.45D0
        CE2=1.90D0
      ELSEIF(ITYPE.EQ.2) THEN
        CE1=1.44D0
        CE2=1.83D0
      ELSEIF(ITYPE.EQ.3) THEN
        CE1=1.45D0
        CE2=1.90D0
      ENDIF
      SH=2.D0
      WRITE(6,*) 'TIME STEP NUMBER N, DT=5/N'
      READ(5,*) ITINC
      DT=5.D0/DFLOAT(ITINC)
      WRITE(6,*) 'OUTPUT INTERVAL NUMBER'
      READ(5,*) NOUT
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE RKLS
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON
     & /BAS/TK,TE,UU,VV,WW,UV
      COMMON
     & /CONMODEL/CE1,CE2
     & /CON1/DT,SH,ANU
     & /CON2/ITYPE,ITINC,NOUT
C
      RT=TK**2/ANU/TE
      ANT=0.09*TK**2/TE*DEXP(-3.4D0/(1.D0+RT/50.D0)**2)
      FE=1.D0-0.3D0*DEXP(-RT**2)
      FTK1=ANT*SH**2-TE
      FTE1=CE1*ANT*SH**2*TE/TK-CE2*TE**2/TK*FE
      TKA=TK+0.5D0*DT*FTK1
      TEA=TE+0.5D0*DT*FTE1
C
      RT=TKA**2/ANU/TEA
      ANT=0.09*TKA**2/TEA*DEXP(-3.4D0/(1.D0+RT/50.D0)**2)
      FE=1.D0-0.3D0*DEXP(-RT**2)
      FTK2=ANT*SH**2-TEA
      FTE2=CE1*ANT*SH**2*TEA/TKA-CE2*TEA**2/TKA*FE
      TKA=TK+0.5D0*DT*FTK2
      TEA=TE+0.5D0*DT*FTE2
C
      RT=TKA**2/ANU/TEA
      ANT=0.09*TKA**2/TEA*DEXP(-3.4D0/(1.D0+RT/50.D0)**2)
      FE=1.D0-0.3D0*DEXP(-RT**2)
      FTK3=ANT*SH**2-TEA
      FTE3=CE1*ANT*SH**2*TEA/TKA-CE2*TEA**2/TKA*FE
      TKA=TK+DT*FTK3
      TEA=TE+DT*FTE3
C
      RT=TKA**2/ANU/TEA
      ANT=0.09*TKA**2/TEA*DEXP(-3.4D0/(1.D0+RT/50.D0)**2)
      FE=1.D0-0.3D0*DEXP(-RT**2)
      FTK4=ANT*SH**2-TEA
      FTE4=CE1*ANT*SH**2*TEA/TKA-CE2*TEA**2/TKA*FE
      TK=TK+DT*(FTK1+2.D0*FTK2+2.D0*FTK3+FTK4)/6.D0
      TE=TE+DT*(FTE1+2.D0*FTE2+2.D0*FTE3+FTE4)/6.D0
C
      RT=TK**2/ANU/TE
      ANT=0.09*TK**2/TE*DEXP(-3.4D0/(1.D0+RT/50.D0)**2)
      UU=2.D0*TK/3.D0
      VV=2.D0*TK/3.D0
      WW=2.D0*TK/3.D0
      UV=-ANT*SH
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE RKLRR
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON
     & /BAS/TK,TE,UU,VV,WW,UV
      COMMON
     & /CONMODEL/CE1,CE2
     & /CON1/DT,SH,ANU
     & /CON2/ITYPE,ITINC,NOUT
C
      TK=0.5*(UU+VV+WW)
      FUU1=-1.2D0*SH*UV+8.D0/15.D0*TE-1.8D0*TE*UU/TK
      FVV1=-0.4D0*SH*UV+8.D0/15.D0*TE-1.8D0*TE*VV/TK
      FWW1=-0.4D0*SH*UV+8.D0/15.D0*TE-1.8D0*TE*WW/TK
      FUV1=-0.4D0*SH*VV-1.8D0*TE*UV/TK
      FTE1=-CE1*SH*UV*TE/TK-CE2*TE**2/TK
      UUA=UU+0.5D0*DT*FUU1
      VVA=VV+0.5D0*DT*FVV1
      WWA=WW+0.5D0*DT*FWW1
      UVA=UV+0.5D0*DT*FUV1
      TEA=TE+0.5D0*DT*FTE1
C
      TK=0.5*(UUA+VVA+WWA)
      FUU2=-1.2D0*SH*UVA+8.D0/15.D0*TEA-1.8D0*TEA*UUA/TK
      FVV2=-0.4D0*SH*UVA+8.D0/15.D0*TEA-1.8D0*TEA*VVA/TK
      FWW2=-0.4D0*SH*UVA+8.D0/15.D0*TEA-1.8D0*TEA*WWA/TK
      FUV2=-0.4D0*SH*VVA-1.8D0*TEA*UVA/TK
      FTE2=-CE1*SH*UVA*TEA/TK-CE2*TEA**2/TK
      UUA=UU+0.5D0*DT*FUU2
      VVA=VV+0.5D0*DT*FVV2
      WWA=WW+0.5D0*DT*FWW2
      UVA=UV+0.5D0*DT*FUV2
      TEA=TE+0.5D0*DT*FTE2
C
      TK=0.5*(UUA+VVA+WWA)
      FUU3=-1.2D0*SH*UVA+8.D0/15.D0*TEA-1.8D0*TEA*UUA/TK
      FVV3=-0.4D0*SH*UVA+8.D0/15.D0*TEA-1.8D0*TEA*VVA/TK
      FWW3=-0.4D0*SH*UVA+8.D0/15.D0*TEA-1.8D0*TEA*WWA/TK
      FUV3=-0.4D0*SH*VVA-1.8D0*TEA*UVA/TK
      FTE3=-CE1*SH*UVA*TEA/TK-CE2*TEA**2/TK
      UUA=UU+DT*FUU3
      VVA=VV+DT*FVV3
      WWA=WW+DT*FWW3
      UVA=UV+DT*FUV3
      TEA=TE+DT*FTE3
C
      TK=0.5*(UUA+VVA+WWA)
      FUU4=-1.2D0*SH*UVA+8.D0/15.D0*TEA-1.8D0*TEA*UUA/TK
      FVV4=-0.4D0*SH*UVA+8.D0/15.D0*TEA-1.8D0*TEA*VVA/TK
      FWW4=-0.4D0*SH*UVA+8.D0/15.D0*TEA-1.8D0*TEA*WWA/TK
      FUV4=-0.4D0*SH*VVA-1.8D0*TEA*UVA/TK
      FTE4=-CE1*SH*UVA*TEA/TK-CE2*TEA**2/TK
      UU=UU+DT*(FUU1+2.D0*FUU2+2.D0*FUU3+FUU4)/6.D0
      VV=VV+DT*(FVV1+2.D0*FVV2+2.D0*FVV3+FVV4)/6.D0
      WW=WW+DT*(FWW1+2.D0*FWW2+2.D0*FWW3+FWW4)/6.D0
      UV=UV+DT*(FUV1+2.D0*FUV2+2.D0*FUV3+FUV4)/6.D0
      TE=TE+DT*(FTE1+2.D0*FTE2+2.D0*FTE3+FTE4)/6.D0
      TK=0.5*(UU+VV+WW)
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE RKSSG
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON
     & /BAS/TK,TE,UU,VV,WW,UV
      COMMON
     & /CONMODEL/CE1,CE2
     & /CON1/DT,SH,ANU
     & /CON2/ITYPE,ITINC,NOUT
C
      A=4.D0/15.D0
      TK=0.5*(UU+VV+WW)
      B11=0.5D0*UU/TK-1.D0/3.D0
      B22=0.5D0*VV/TK-1.D0/3.D0
      B33=0.5D0*WW/TK-1.D0/3.D0
      B12=0.5D0*UV/TK
      FUU1=-191.D0/60.D0*SH*TK*B12-2.D0/3.D0*TE-3.4D0*TE*B11
     &     +3.6D0*SH*TK*B11*B12
     &     +1.4D0*TE*(2.D0*B11**2-B22**2-B33**2+B12**2)
      FVV1=1.D0/60.D0*SH*TK*B12-2.D0/3.D0*TE-3.4D0*TE*B22
     &     +3.6D0*SH*TK*B22*B12
     &     +1.4D0*TE*(-B11**2+2.D0*B22**2-B33**2+B12**2)
      FWW1=-5.D0/6.D0*SH*TK*B12-2.D0/3.D0*TE-3.4D0*TE*B33
     &     +3.6D0*SH*TK*B33*B12
     &     +1.4D0*TE*(-B11**2-B22**2+2.D0*B33**2-2.D0*B12**2)
      FUV1=-3.4D0*TE*B12+3.6D0*SH*TK*B12**2+4.2D0*TE*B12*(B11+B22)
     &     -(A+0.65D0*DSQRT(B11**2+B22**2+B33**2+2.D0*B12**2))*SH*TK
     &     +0.425D0*SH*TK*B11-1.175D0*SH*TK*B22
      FTE1=-CE1*SH*UV*TE/TK-CE2*TE**2/TK
      UUA=UU+0.5D0*DT*FUU1
      VVA=VV+0.5D0*DT*FVV1
      WWA=WW+0.5D0*DT*FWW1
      UVA=UV+0.5D0*DT*FUV1
      TEA=TE+0.5D0*DT*FTE1
C
      TK=0.5*(UUA+VVA+WWA)
      B11=0.5D0*UUA/TK-1.D0/3.D0
      B22=0.5D0*VVA/TK-1.D0/3.D0
      B33=0.5D0*WWA/TK-1.D0/3.D0
      B12=0.5D0*UVA/TK
      FUU2=-191.D0/60.D0*SH*TK*B12-2.D0/3.D0*TEA-3.4D0*TEA*B11
     &     +3.6D0*SH*TK*B11*B12
     &     +1.4D0*TEA*(2.D0*B11**2-B22**2-B33**2+B12**2)
      FVV2=1.D0/60.D0*SH*TK*B12-2.D0/3.D0*TEA-3.4D0*TEA*B22
     &     +3.6D0*SH*TK*B22*B12
     &     +1.4D0*TEA*(-B11**2+2.D0*B22**2-B33**2+B12**2)
      FWW2=-5.D0/6.D0*SH*TK*B12-2.D0/3.D0*TEA-3.4D0*TEA*B33
     &     +3.6D0*SH*TK*B33*B12
     &     +1.4D0*TEA*(-B11**2-B22**2+2.D0*B33**2-2.D0*B12**2)
      FUV2=-3.4D0*TEA*B12+3.6D0*SH*TK*B12**2+4.2D0*TEA*B12*(B11+B22)
     &     -(A+0.65D0*DSQRT(B11**2+B22**2+B33**2+2.D0*B12**2))*SH*TK
     &     +0.425D0*SH*TK*B11-1.175D0*SH*TK*B22
      FTE2=-CE1*SH*UVA*TEA/TK-CE2*TEA**2/TK
      UUA=UU+0.5D0*DT*FUU2
      VVA=VV+0.5D0*DT*FVV2
      WWA=WW+0.5D0*DT*FWW2
      UVA=UV+0.5D0*DT*FUV2
      TEA=TE+0.5D0*DT*FTE2
C
      TK=0.5*(UUA+VVA+WWA)
      B11=0.5D0*UUA/TK-1.D0/3.D0
      B22=0.5D0*VVA/TK-1.D0/3.D0
      B33=0.5D0*WWA/TK-1.D0/3.D0
      B12=0.5D0*UVA/TK
      FUU3=-191.D0/60.D0*SH*TK*B12-2.D0/3.D0*TEA-3.4D0*TEA*B11
     &     +3.6D0*SH*TK*B11*B12
     &     +1.4D0*TEA*(2.D0*B11**2-B22**2-B33**2+B12**2)
      FVV3=1.D0/60.D0*SH*TK*B12-2.D0/3.D0*TEA-3.4D0*TEA*B22
     &     +3.6D0*SH*TK*B22*B12
     &     +1.4D0*TEA*(-B11**2+2.D0*B22**2-B33**2+B12**2)
      FWW3=-5.D0/6.D0*SH*TK*B12-2.D0/3.D0*TEA-3.4D0*TEA*B33
     &     +3.6D0*SH*TK*B33*B12
     &     +1.4D0*TEA*(-B11**2-B22**2+2.D0*B33**2-2.D0*B12**2)
      FUV3=-3.4D0*TEA*B12+3.6D0*SH*TK*B12**2+4.2D0*TEA*B12*(B11+B22)
     &     -(A+0.65D0*DSQRT(B11**2+B22**2+B33**2+2.D0*B12**2))*SH*TK
     &     +0.425D0*SH*TK*B11-1.175D0*SH*TK*B22
      FTE3=-CE1*SH*UVA*TEA/TK-CE2*TEA**2/TK
      UUA=UU+DT*FUU3
      VVA=VV+DT*FVV3
      WWA=WW+DT*FWW3
      UVA=UV+DT*FUV3
      TEA=TE+DT*FTE3
C
      TK=0.5*(UUA+VVA+WWA)
      B11=0.5D0*UUA/TK-1.D0/3.D0
      B22=0.5D0*VVA/TK-1.D0/3.D0
      B33=0.5D0*WWA/TK-1.D0/3.D0
      B12=0.5D0*UVA/TK
      FUU4=-191.D0/60.D0*SH*TK*B12-2.D0/3.D0*TEA-3.4D0*TEA*B11
     &     +3.6D0*SH*TK*B11*B12
     &     +1.4D0*TEA*(2.D0*B11**2-B22**2-B33**2+B12**2)
      FVV4=1.D0/60.D0*SH*TK*B12-2.D0/3.D0*TEA-3.4D0*TEA*B22
     &     +3.6D0*SH*TK*B22*B12
     &     +1.4D0*TEA*(-B11**2+2.D0*B22**2-B33**2+B12**2)
      FWW4=-5.D0/6.D0*SH*TK*B12-2.D0/3.D0*TEA-3.4D0*TEA*B33
     &     +3.6D0*SH*TK*B33*B12
     &     +1.4D0*TEA*(-B11**2-B22**2+2.D0*B33**2-2.D0*B12**2)
      FUV4=-3.4D0*TEA*B12+3.6D0*SH*TK*B12**2+4.2D0*TEA*B12*(B11+B22)
     &     -(A+0.65D0*DSQRT(B11**2+B22**2+B33**2+2.D0*B12**2))*SH*TK
     &     +0.425D0*SH*TK*B11-1.175D0*SH*TK*B22
      FTE4=-CE1*SH*UVA*TEA/TK-CE2*TEA**2/TK
      UU=UU+DT*(FUU1+2.D0*FUU2+2.D0*FUU3+FUU4)/6.D0
      VV=VV+DT*(FVV1+2.D0*FVV2+2.D0*FVV3+FVV4)/6.D0
      WW=WW+DT*(FWW1+2.D0*FWW2+2.D0*FWW3+FWW4)/6.D0
      UV=UV+DT*(FUV1+2.D0*FUV2+2.D0*FUV3+FUV4)/6.D0
      TE=TE+DT*(FTE1+2.D0*FTE2+2.D0*FTE3+FTE4)/6.D0
      TK=0.5*(UU+VV+WW)
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE RKASM
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON
     & /BAS/TK,TE,UU,VV,WW,UV
      COMMON
     & /CONMODEL/CE1,CE2
     & /CON1/DT,SH,ANU
     & /CON2/ITYPE,ITINC,NOUT
C
      TKI=TK
      TEI=TE
      VVI=VV
      UVI=UV
      CALL SASM(TKI,TEI,VVI,UVI)
      FTK1=-UVI*SH-TE
      FTE1=-CE1*UVI*SH*TE/TK-CE2*TE**2/TK
      TKA=TK+0.5D0*DT*FTK1
      TEA=TE+0.5D0*DT*FTE1
C
      TKI=TKA
      TEI=TEA
      CALL SASM(TKI,TEI,VVI,UVI)
      FTK2=-UVI*SH-TEA
      FTE2=-CE1*UVI*SH*TEA/TKA-CE2*TEA**2/TKA
      TKA=TK+0.5D0*DT*FTK2
      TEA=TE+0.5D0*DT*FTE2
C
      TKI=TKA
      TEI=TEA
      CALL SASM(TKI,TEI,VVI,UVI)
      FTK3=-UVI*SH-TEA
      FTE3=-CE1*UVI*SH*TEA/TKA-CE2*TEA**2/TKA
      TKA=TK+DT*FTK3
      TEA=TE+DT*FTE3
C
      TKI=TKA
      TEI=TEA
      CALL SASM(TKI,TEI,VVI,UVI)
      FTK4=-UVI*SH-TEA
      FTE4=-CE1*UVI*SH*TEA/TKA-CE2*TEA**2/TKA
      TK=TK+DT*(FTK1+2.D0*FTK2+2.D0*FTK3+FTK4)/6.D0
      TE=TE+DT*(FTE1+2.D0*FTE2+2.D0*FTE3+FTE4)/6.D0
C
      VV=VVI
      UV=UVI
      CALL SASM(TK,TE,VV,UV)
      WW=VV
      UU=2.D0*TK-VV-WW
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE SASM(TKI,TEI,VVI,UVI)
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON
     & /CONMODEL/CE1,CE2
     & /CON1/DT,SH,ANU
     & /CON2/ITYPE,ITINC,NOUT
C
      C=32.D0/75.D0
      DO 100 ITER=1,200
        VVO=VVI
        VVI=VVI-(VVI**3-0.8D0*TKI*VVI**2
     &          +(0.16D0*TKI**2+C*TEI**2/SH**2)*VVI
     &          -64.D0/225.D0*TEI**2*TKI/SH**2)
     &   /(3.D0*VVI**2-1.6D0*TKI*VVI+0.16D0*TKI**2+C*TEI**2/SH**2)
        ERR=DABS((VVI-VVO)/VVO)
        IF(ERR.LT.1.D-5) GOTO 1000
 100  CONTINUE
1000  UVI=0.8D0*TEI*(2.D0/3.D0-VVI/TKI)/SH/(0.4D0-VVI/TKI)
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
     & /BAS/TK,TE,UU,VV,WW,UV
      COMMON
     & /CONMODEL/CE1,CE2
     & /CON1/DT,SH,ANU
     & /CON2/ITYPE,ITINC,NOUT
C
      TIME=DFLOAT(IT)*DT
      WRITE(12,'( 1H ,1P,7E14.4 )') TIME,TK,TE,UU,VV,WW,UV
      RETURN
      END
