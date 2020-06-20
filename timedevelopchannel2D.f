C
C********* ********* ********* ********* ********* ********* *********C
C
C     PROGRAM timedevelopchannel2D.f
C
C     Program for 2D channel flow by FDM
C     MADE BY MASAYOSHI OKAMOTO
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NX=40,NY=40)
      COMMON
     & /BVL/U(0:NX+1,0:NY+1),V(0:NX+1,0:NY+1)
      COMMON
     & /TVL/UA(0:NX+1,0:NY+1),VA(0:NX+1,0:NY+1)
      COMMON
     & /PRS/PC(0:NX+1,0:NY+1),PA(0:NX+1,0:NY+1)
      COMMON
     & /XGD/X(0:NX),XX(0:NX)
     & /YGD/Y(0:NY),YY(0:NY)
      COMMON
     & /CON/DT,DX,DY,ANU,OMG
     & /CNM/NIT,IPTER
C
      CALL INIT
      CALL BDCD
      DO 100 IT=1,NIT
        CALL SOLVEVEL
        CALL SOLA(IT)
        IF(MOD(IT,200).EQ.0) THEN
          IFILE=IT/200
          WRITE(6,*) U(20,20)
          CALL ENPR(IFILE)
        ENDIF
 100  CONTINUE
      WRITE(11) U,V,PA
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE INIT
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NX=40,NY=40)
      COMMON
     & /BVL/U(0:NX+1,0:NY+1),V(0:NX+1,0:NY+1)
      COMMON
     & /TVL/UA(0:NX+1,0:NY+1),VA(0:NX+1,0:NY+1)
      COMMON
     & /PRS/PC(0:NX+1,0:NY+1),PA(0:NX+1,0:NY+1)
      COMMON
     & /XGD/X(0:NX),XX(0:NX)
     & /YGD/Y(0:NY),YY(0:NY)
      COMMON
     & /CON/DT,DX,DY,ANU,OMG
     & /CNM/NIT,IPTER
C
      DO 100 J=0,NY+1
      DO 100 I=0,NX+1
        U(I,J)=0.D0
        V(I,J)=0.D0
        UA(I,J)=0.D0
        VA(I,J)=0.D0
        PA(I,J)=0.D0
        PC(I,J)=0.D0
 100  CONTINUE
      DO 110 I=0,NX+1
        X(I)=0.D0
        XX(I)=0.D0
 110  CONTINUE
      DO 120 J=0,NY+1
        Y(J)=0.D0
        YY(J)=0.D0
 120  CONTINUE
      OMG=1.7D0
      RE=100.D0
      ANU=1.D0/RE
      DT=1.D-2
      NIT=5000
      IPTER=50000
      DX=10.D0/DFLOAT(NX)
      DY=1.D0/DFLOAT(NY)
      DO 130 I=0,NX
        XX(I)=DX*DFLOAT(I)
        YY(I)=DY*DFLOAT(I)
 130  CONTINUE
      DO 140 I=1,NY
        X(I)=XX(I)-0.5D0*DX
        Y(I)=YY(I)-0.5D0*DY
 140  CONTINUE
C      READ(11) U,V,PA
C      REWIND 11
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE BDCD
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NX=40,NY=40)
      COMMON
     & /BVL/U(0:NX+1,0:NY+1),V(0:NX+1,0:NY+1)
C
      DO 100 I=0,NX+1
        U(I,NY+1)=-U(I,NY)
        V(I,NY)=0.D0
        U(I,0)=-U(I,1)
        V(I,0)=0.D0
 100  CONTINUE
      DO 110 J=0,NY+1
        U(NX,J)=U(NX-1,J)
        V(NX+1,J)=V(NX,J)
        U(0,J)=1.D0
        V(0,J)=-V(1,J)
 110  CONTINUE
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE SOLVEVEL
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NX=40,NY=40)
      COMMON
     & /BVL/U(0:NX+1,0:NY+1),V(0:NX+1,0:NY+1)
      COMMON
     & /TVL/UA(0:NX+1,0:NY+1),VA(0:NX+1,0:NY+1)
      COMMON
     & /PRS/PC(0:NX+1,0:NY+1),PA(0:NX+1,0:NY+1)
      COMMON
     & /CON/DT,DX,DY,ANU,OMG
     & /CNM/NIT,IPTER
C
C********** STATEMENT FUNCTION **********
C
      UO(I,J)=0.5D0*(U(I,J)+U(I-1,J))
      VO(I,J)=0.5D0*(V(I,J)+V(I,J-1))
      UXY(I,J)=0.5D0*(U(I,J+1)+U(I,J))
      VXY(I,J)=0.5D0*(V(I+1,J)+V(I,J))
C
      DO 100 J=1,NY
      DO 100 I=1,NX-1
        UA(I,J)=U(I,J)
     &    +DT*(-(UO(I+1,J)**2-UO(I,J)**2)/DX
     &         -(VXY(I,J)*UXY(I,J)-VXY(I,J-1)*UXY(I,J-1))/DY
     &         -(PA(I+1,J)-PA(I,J))/DX
     &         +ANU*(U(I+1,J)-2.D0*U(I,J)+U(I-1,J))/DX**2
     &         +ANU*(U(I,J+1)-2.D0*U(I,J)+U(I,J-1))/DY**2)
 100  CONTINUE
      DO 110 J=1,NY-1
      DO 110 I=1,NX
        VA(I,J)=V(I,J)
     &    +DT*(-(UXY(I,J)*VXY(I,J)-UXY(I-1,J)*VXY(I-1,J))/DX
     &         -(VO(I,J+1)**2-VO(I,J)**2)/DY
     &         -(PA(I,J+1)-PA(I,J))/DY
     &         +ANU*(V(I+1,J)-2.D0*V(I,J)+V(I-1,J))/DX**2
     &         +ANU*(V(I,J+1)-2.D0*V(I,J)+V(I,J-1))/DY**2)
 110  CONTINUE
      DO 120 I=0,NX+1
        UA(I,NY+1)=-UA(I,NY)
        VA(I,NY)=0.D0
        UA(I,0)=-UA(I,1)
        VA(I,0)=0.D0
 120  CONTINUE
      DO 130 J=0,NY+1
        UA(NX,J)=UA(NX-1,J)
        VA(NX+1,J)=VA(NX,J)
        UA(0,J)=1.D0
        VA(0,J)=-VA(1,J)
 130  CONTINUE
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE SOLA(IT)
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NX=40,NY=40)
      COMMON
     & /BVL/U(0:NX+1,0:NY+1),V(0:NX+1,0:NY+1)
      COMMON
     & /TVL/UA(0:NX+1,0:NY+1),VA(0:NX+1,0:NY+1)
      COMMON
     & /PRS/PC(0:NX+1,0:NY+1),PA(0:NX+1,0:NY+1)
      COMMON
     & /CON/DT,DX,DY,ANU,OMG
     & /CNM/NIT,IPTER
      DIMENSION DV(0:NX+1,0:NY+1)
C
      DO 100 ITER=1,IPTER
        DVMAX=0.D0
        DO 110 J=1,NY
        DO 110 I=1,NX
          DV(I,J)=(UA(I,J)-UA(I-1,J))/DX+(VA(I,J)-VA(I,J-1))/DY
          DVMAX=DMAX1(DABS(DV(I,J)),DVMAX)
          PC(I,J)=-0.5D0*OMG*DV(I,J)/DT/(1.D0/DX**2+1.D0/DY**2)
          UA(I,J)=UA(I,J)+DT/DX*PC(I,J)
          UA(I-1,J)=UA(I-1,J)-DT/DX*PC(I,J)
          VA(I,J)=VA(I,J)+DT/DY*PC(I,J)
          VA(I,J-1)=VA(I,J-1)-DT/DY*PC(I,J)
          PA(I,J)=PA(I,J)+PC(I,J)
 110    CONTINUE
        IF(DVMAX.LT.1.D-10) GOTO 200
 100  CONTINUE
 200  IF(MOD(IT,1000).EQ.0) THEN
        WRITE(6,*) ITER,DVMAX
      ENDIF
      DO 160 J=1,NY
      DO 160 I=1,NX-1
        U(I,J)=UA(I,J)
 160  CONTINUE
      DO 170 J=1,NY-1
      DO 170 I=1,NX
        V(I,J)=VA(I,J)
 170  CONTINUE
      CALL BDCD
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE ENPR(IFILE)
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NX=40,NY=40)
      COMMON
     & /BVL/U(0:NX+1,0:NY+1),V(0:NX+1,0:NY+1)
      COMMON
     & /PRS/PC(0:NX+1,0:NY+1),PA(0:NX+1,0:NY+1)
      COMMON
     & /XGD/X(0:NX),XX(0:NX)
     & /YGD/Y(0:NY),YY(0:NY)
      COMMON
     & /CON/DT,DX,DY,ANU,OMG
     & /CNM/NIT,IPTER
      CHARACTER FILENAME*128
C
      WRITE(FILENAME,'("DATA",I4.4,".txt")') IFILE
      OPEN(40,FILE=FILENAME,STATUS='REPLACE')
      DO 100 J=1,NY
      DO 100 I=1,NX
        WRITE(40,'( 1H ,1P,5E14.4 )') X(I),Y(J)
     &  ,0.5D0*(U(I,J)+U(I-1,J)),0.5D0*(V(I,J)+V(I,J-1))
     &  ,PA(I,J)
 100  CONTINUE
      CLOSE(40)
      RETURN
      END
