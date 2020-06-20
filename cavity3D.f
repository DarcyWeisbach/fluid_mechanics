C
C********* ********* ********* ********* ********* ********* *********C
C
C     PROGRAM cavity3D.f
C
C     Program for 3D cavity flow by FVM
C     MADE BY MASAYOSHI OKAMOTO
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NX=40,NY=40,NZ=40)
      COMMON
     & /BVL/UA(0:NX+1,0:NY+1,0:NZ+1),VA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,WA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /PRS/PC(0:NX+1,0:NY+1,0:NZ+1),PA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /UCF/AEUA(0:NX+1,0:NY+1,0:NZ+1),AWUA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ANUA(0:NX+1,0:NY+1,0:NZ+1),ASUA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ATUA(0:NX+1,0:NY+1,0:NZ+1),ABUA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,APUA(0:NX+1,0:NY+1,0:NZ+1),BUA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /VCF/AEVA(0:NX+1,0:NY+1,0:NZ+1),AWVA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ANVA(0:NX+1,0:NY+1,0:NZ+1),ASVA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ATVA(0:NX+1,0:NY+1,0:NZ+1),ABVA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,APVA(0:NX+1,0:NY+1,0:NZ+1),BVA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /WCF/AEWA(0:NX+1,0:NY+1,0:NZ+1),AWWA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ANWA(0:NX+1,0:NY+1,0:NZ+1),ASWA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ATWA(0:NX+1,0:NY+1,0:NZ+1),ABWA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,APWA(0:NX+1,0:NY+1,0:NZ+1),BWA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /PCF/AEPC(0:NX+1,0:NY+1,0:NZ+1),AWPC(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ANPC(0:NX+1,0:NY+1,0:NZ+1),ASPC(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ATPC(0:NX+1,0:NY+1,0:NZ+1),ABPC(0:NX+1,0:NY+1,0:NZ+1)
     &     ,APPC(0:NX+1,0:NY+1,0:NZ+1),BPC(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /XGD/X(0:NX),XX(0:NX)
     & /YGD/Y(0:NY),YY(0:NY)
     & /ZGD/Z(0:NZ),ZZ(0:NZ)
      COMMON
     & /CON/DX,DY,DZ,ANU,ALFUA,ALFVA,ALFWA,ALFPA
     & /CNM/NITER,IPTER
C
      CALL INIT
      CALL BDCD
      DO 100 ITER=1,NITER
        CALL SOLVEUA
        CALL SOLVEVA
        CALL SOLVEWA
        CALL SOLVEPC
        CALL UPDATE
        IF(MOD(ITER,1000).EQ.0) THEN
          WRITE(6,*) ITER
          CALL CONVERG(REST)
          IF(REST.LT.1.D-5) GOTO 110
        ENDIF
 100  CONTINUE
 110  CALL ENPR
      WRITE(11) UA,VA,WA,PA
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE INIT
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NX=40,NY=40,NZ=40)
      COMMON
     & /BVL/UA(0:NX+1,0:NY+1,0:NZ+1),VA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,WA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /PRS/PC(0:NX+1,0:NY+1,0:NZ+1),PA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /UCF/AEUA(0:NX+1,0:NY+1,0:NZ+1),AWUA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ANUA(0:NX+1,0:NY+1,0:NZ+1),ASUA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ATUA(0:NX+1,0:NY+1,0:NZ+1),ABUA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,APUA(0:NX+1,0:NY+1,0:NZ+1),BUA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /VCF/AEVA(0:NX+1,0:NY+1,0:NZ+1),AWVA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ANVA(0:NX+1,0:NY+1,0:NZ+1),ASVA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ATVA(0:NX+1,0:NY+1,0:NZ+1),ABVA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,APVA(0:NX+1,0:NY+1,0:NZ+1),BVA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /WCF/AEWA(0:NX+1,0:NY+1,0:NZ+1),AWWA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ANWA(0:NX+1,0:NY+1,0:NZ+1),ASWA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ATWA(0:NX+1,0:NY+1,0:NZ+1),ABWA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,APWA(0:NX+1,0:NY+1,0:NZ+1),BWA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /PCF/AEPC(0:NX+1,0:NY+1,0:NZ+1),AWPC(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ANPC(0:NX+1,0:NY+1,0:NZ+1),ASPC(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ATPC(0:NX+1,0:NY+1,0:NZ+1),ABPC(0:NX+1,0:NY+1,0:NZ+1)
     &     ,APPC(0:NX+1,0:NY+1,0:NZ+1),BPC(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /XGD/X(0:NX),XX(0:NX)
     & /YGD/Y(0:NY),YY(0:NY)
     & /ZGD/Z(0:NZ),ZZ(0:NZ)
      COMMON
     & /CON/DX,DY,DZ,ANU,ALFUA,ALFVA,ALFWA,ALFPA
     & /CNM/NITER,IPTER
C
      DO 100 K=0,NZ+1
      DO 100 J=0,NY+1
      DO 100 I=0,NX+1
        UA(I,J,K)=0.D0
        VA(I,J,K)=0.D0
        WA(I,J,K)=0.D0
        PA(I,J,K)=0.D0
        PC(I,J,K)=0.D0
        APUA(I,J,K)=0.D0
        AEUA(I,J,K)=0.D0
        AWUA(I,J,K)=0.D0
        ANUA(I,J,K)=0.D0
        ASUA(I,J,K)=0.D0
        ATUA(I,J,K)=0.D0
        ABUA(I,J,K)=0.D0
        BUA(I,J,K)=0.D0
        APVA(I,J,K)=0.D0
        AEVA(I,J,K)=0.D0
        AWVA(I,J,K)=0.D0
        ANVA(I,J,K)=0.D0
        ASVA(I,J,K)=0.D0
        ATVA(I,J,K)=0.D0
        ABVA(I,J,K)=0.D0
        BVA(I,J,K)=0.D0
        APWA(I,J,K)=0.D0
        AEWA(I,J,K)=0.D0
        AWWA(I,J,K)=0.D0
        ANWA(I,J,K)=0.D0
        ASWA(I,J,K)=0.D0
        ATWA(I,J,K)=0.D0
        ABWA(I,J,K)=0.D0
        BWA(I,J,K)=0.D0
        APPC(I,J,K)=0.D0
        AEPC(I,J,K)=0.D0
        AWPC(I,J,K)=0.D0
        ANPC(I,J,K)=0.D0
        ASPC(I,J,K)=0.D0
        ATPC(I,J,K)=0.D0
        ABPC(I,J,K)=0.D0
        BPC(I,J,K)=0.D0
 100  CONTINUE
      DO 110 I=0,NX+1
        X(I)=0.D0
        XX(I)=0.D0
 110  CONTINUE
      DO 120 J=0,NY+1
        Y(J)=0.D0
        YY(J)=0.D0
 120  CONTINUE
      DO 130 K=0,NZ+1
        Z(K)=0.D0
        ZZ(K)=0.D0
 130  CONTINUE
      ALFUA=0.1D0
      ALFVA=0.1D0
      ALFWA=0.1D0
      ALFPA=0.1D0
      RE=100.D0
      ANU=1.D0/RE
      NITER=100000
      IPTER=5
      DX=1.D0/FLOAT(NX)
      DY=1.D0/FLOAT(NY)
      DZ=1.D0/FLOAT(NZ)
      DO 140 I=0,NX
        XX(I)=DX*FLOAT(I)
 140  CONTINUE
      DO 150 I=1,NX
        X(I)=XX(I)-0.5D0*DX
 150  CONTINUE
      DO 160 J=0,NY
        YY(J)=DY*FLOAT(J)
 160  CONTINUE
      DO 170 J=1,NY
        Y(J)=YY(J)-0.5D0*DY
 170  CONTINUE
      DO 180 K=0,NZ
        ZZ(K)=DZ*FLOAT(K)
 180  CONTINUE
      DO 190 K=1,NZ
        Z(K)=ZZ(K)-0.5D0*DZ
 190  CONTINUE
C      READ(11) UA,VA,WA,PA
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
      PARAMETER(NX=40,NY=40,NZ=40)
      COMMON
     & /BVL/UA(0:NX+1,0:NY+1,0:NZ+1),VA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,WA(0:NX+1,0:NY+1,0:NZ+1)
C
      DO 100 K=0,NZ+1
      DO 100 I=0,NX+1
        UA(I,NY+1,K)=1.D0
        VA(I,NY,K)=0.D0
        WA(I,NY+1,K)=0.D0
        UA(I,0,K)=0.D0
        VA(I,0,K)=0.D0
        WA(I,0,K)=0.D0
 100  CONTINUE
      DO 110 K=0,NZ+1
      DO 110 J=0,NY+1
        UA(NX,J,K)=0.D0
        VA(NX+1,J,K)=0.D0
        WA(NX+1,J,K)=0.D0
        UA(0,J,K)=0.D0
        VA(0,J,K)=0.D0
        WA(0,J,K)=0.D0
 110  CONTINUE
      DO 120 J=0,NY+1
      DO 120 I=0,NX+1
        UA(I,J,NZ+1)=0.D0
        VA(I,J,NZ+1)=0.D0
        WA(I,J,NZ)=0.D0
        UA(I,J,0)=0.D0
        VA(I,J,0)=0.D0
        WA(I,J,0)=0.D0
 120  CONTINUE
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE SOLVEUA
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NX=40,NY=40,NZ=40)
      COMMON
     & /BVL/UA(0:NX+1,0:NY+1,0:NZ+1),VA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,WA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /PRS/PC(0:NX+1,0:NY+1,0:NZ+1),PA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /UCF/AEUA(0:NX+1,0:NY+1,0:NZ+1),AWUA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ANUA(0:NX+1,0:NY+1,0:NZ+1),ASUA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ATUA(0:NX+1,0:NY+1,0:NZ+1),ABUA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,APUA(0:NX+1,0:NY+1,0:NZ+1),BUA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /CON/DX,DY,DZ,ANU,ALFUA,ALFVA,ALFWA,ALFPA
     & /CNM/NITER,IPTER
C
      DO 100 K=1,NZ
      DO 100 J=1,NY
      DO 100 I=1,NX-1
        AWUA(I,J,K)=ANU*DY*DZ/DX
     &             +0.25D0*DY*DZ*(UA(I,J,K)+UA(I-1,J,K))
        AEUA(I,J,K)=ANU*DY*DZ/DX
     &             -0.25D0*DY*DZ*(UA(I+1,J,K)+UA(I,J,K))
        IF(J.EQ.1) THEN
          ASUA(I,J,K)=2.D0*ANU*DX*DZ/DY
        ELSE
          ASUA(I,J,K)=ANU*DX*DZ/DY
     &               +0.25D0*DX*DZ*(VA(I+1,J-1,K)+VA(I,J-1,K))
        ENDIF
        IF(J.EQ.NY) THEN
          ANUA(I,J,K)=2.D0*ANU*DX*DZ/DY
        ELSE
          ANUA(I,J,K)=ANU*DX*DZ/DY
     &               -0.25D0*DX*DZ*(VA(I+1,J,K)+VA(I,J,K))
        ENDIF
        IF(K.EQ.1) THEN
          ABUA(I,J,K)=2.D0*ANU*DX*DY/DZ
        ELSE
          ABUA(I,J,K)=ANU*DX*DY/DZ
     &               +0.25D0*DX*DY*(WA(I+1,J,K-1)+WA(I,J,K-1))
        ENDIF
        IF(K.EQ.NZ) THEN
          ATUA(I,J,K)=2.D0*ANU*DX*DY/DZ
        ELSE
          ATUA(I,J,K)=ANU*DX*DY/DZ
     &               -0.25D0*DX*DY*(WA(I+1,J,K)+WA(I,J,K))
        ENDIF
        APUA(I,J,K)=AEUA(I,J,K)+AWUA(I,J,K)+ANUA(I,J,K)+ASUA(I,J,K)
     &             +ATUA(I,J,K)+ABUA(I,J,K)
        BUA(I,J,K)=-DY*DZ*(PA(I+1,J,K)-PA(I,J,K))
        UA(I,J,K)=UA(I,J,K)*(1.D0-ALFUA)
     &        +ALFUA*(AEUA(I,J,K)*UA(I+1,J,K)+AWUA(I,J,K)*UA(I-1,J,K)
     &               +ANUA(I,J,K)*UA(I,J+1,K)+ASUA(I,J,K)*UA(I,J-1,K)
     &               +ATUA(I,J,K)*UA(I,J,K+1)+ABUA(I,J,K)*UA(I,J,K-1)
     &               +BUA(I,J,K))/APUA(I,J,K)
 100  CONTINUE
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE SOLVEVA
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NX=40,NY=40,NZ=40)
      COMMON
     & /BVL/UA(0:NX+1,0:NY+1,0:NZ+1),VA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,WA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /PRS/PC(0:NX+1,0:NY+1,0:NZ+1),PA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /VCF/AEVA(0:NX+1,0:NY+1,0:NZ+1),AWVA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ANVA(0:NX+1,0:NY+1,0:NZ+1),ASVA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ATVA(0:NX+1,0:NY+1,0:NZ+1),ABVA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,APVA(0:NX+1,0:NY+1,0:NZ+1),BVA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /CON/DX,DY,DZ,ANU,ALFUA,ALFVA,ALFWA,ALFPA
     & /CNM/NITER,IPTER
C
      DO 100 K=1,NZ
      DO 100 J=1,NY-1
      DO 100 I=1,NX
        IF(I.EQ.1) THEN
          AWVA(I,J,K)=2.0D0*ANU*DY*DZ/DX
        ELSE
          AWVA(I,J,K)=ANU*DY*DZ/DX
     &               +0.25D0*DY*DZ*(UA(I-1,J+1,K)+UA(I-1,J,K))
        ENDIF
        IF(I.EQ.NX) THEN
          AEVA(I,J,K)=2.D0*ANU*DY*DZ/DX
        ELSE
          AEVA(I,J,K)=ANU*DY*DZ/DX
     &               -0.25D0*DY*DZ*(UA(I,J+1,K)+UA(I,J,K))
        ENDIF
        ANVA(I,J,K)=ANU*DX*DZ/DY
     &             -0.25D0*DX*DZ*(VA(I,J+1,K)+VA(I,J,K))
        ASVA(I,J,K)=ANU*DX*DZ/DY
     &             +0.25D0*DX*DZ*(VA(I,J,K)+VA(I,J-1,K))
        IF(K.EQ.1) THEN
          ABVA(I,J,K)=2.D0*ANU*DX*DY/DZ
        ELSE
          ABVA(I,J,K)=ANU*DX*DY/DZ
     &               +0.25D0*DX*DY*(WA(I,J+1,K-1)+WA(I,J,K-1))
        ENDIF
        IF(K.EQ.NZ) THEN
          ATVA(I,J,K)=2.D0*ANU*DX*DY/DZ
        ELSE
          ATVA(I,J,K)=ANU*DX*DY/DZ
     &               -0.25D0*DX*DY*(WA(I,J+1,K)+WA(I,J,K))
        ENDIF
        APVA(I,J,K)=AEVA(I,J,K)+AWVA(I,J,K)+ANVA(I,J,K)+ASVA(I,J,K)
     &             +ATVA(I,J,K)+ABVA(I,J,K)
        BVA(I,J,K)=-DX*DZ*(PA(I,J+1,K)-PA(I,J,K))
        VA(I,J,K)=VA(I,J,K)*(1.D0-ALFVA)
     &        +ALFVA*(AEVA(I,J,K)*VA(I+1,J,K)+AWVA(I,J,K)*VA(I-1,J,K)
     &               +ANVA(I,J,K)*VA(I,J+1,K)+ASVA(I,J,K)*VA(I,J-1,K)
     &               +ATVA(I,J,K)*VA(I,J,K+1)+ABVA(I,J,K)*VA(I,J,K-1)
     &               +BVA(I,J,K))/APVA(I,J,K)
 100  CONTINUE
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE SOLVEWA
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NX=40,NY=40,NZ=40)
      COMMON
     & /BVL/UA(0:NX+1,0:NY+1,0:NZ+1),VA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,WA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /PRS/PC(0:NX+1,0:NY+1,0:NZ+1),PA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /WCF/AEWA(0:NX+1,0:NY+1,0:NZ+1),AWWA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ANWA(0:NX+1,0:NY+1,0:NZ+1),ASWA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ATWA(0:NX+1,0:NY+1,0:NZ+1),ABWA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,APWA(0:NX+1,0:NY+1,0:NZ+1),BWA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /CON/DX,DY,DZ,ANU,ALFUA,ALFVA,ALFWA,ALFPA
     & /CNM/NITER,IPTER
C
      DO 100 K=1,NZ-1
      DO 100 J=1,NY
      DO 100 I=1,NX
        IF(I.EQ.1) THEN
          AWWA(I,J,K)=2.0D0*ANU*DY*DZ/DX
        ELSE
          AWWA(I,J,K)=ANU*DY*DZ/DX
     &               +0.25D0*DY*DZ*(UA(I-1,J,K+1)+UA(I-1,J,K))
        ENDIF
        IF(I.EQ.NX) THEN
          AEWA(I,J,K)=2.D0*ANU*DY*DZ/DX
        ELSE
          AEWA(I,J,K)=ANU*DY*DZ/DX
     &               -0.25D0*DY*DZ*(UA(I,J,K+1)+UA(I,J,K))
        ENDIF
        IF(J.EQ.1) THEN
          ASWA(I,J,K)=2.D0*ANU*DX*DZ/DY
        ELSE
          ASWA(I,J,K)=ANU*DX*DZ/DY
     &               +0.25D0*DX*DZ*(VA(I,J-1,K+1)+VA(I,J-1,K))
        ENDIF
        IF(J.EQ.NY) THEN
          ANWA(I,J,K)=2.D0*ANU*DX*DZ/DY
        ELSE
          ANWA(I,J,K)=ANU*DX*DZ/DY
     &               -0.25D0*DX*DZ*(VA(I,J,K+1)+VA(I,J,K))
        ENDIF
        ABWA(I,J,K)=ANU*DX*DY/DZ
     &             +0.25D0*DX*DY*(WA(I,J,K)+WA(I,J,K-1))
        ATWA(I,J,K)=ANU*DX*DY/DZ
     &             -0.25D0*DX*DY*(WA(I,J,K+1)+WA(I,J,K))
        APWA(I,J,K)=AEWA(I,J,K)+AWWA(I,J,K)+ANWA(I,J,K)+ASWA(I,J,K)
     &             +ATWA(I,J,K)+ABWA(I,J,K)
        BWA(I,J,K)=-DX*DY*(PA(I,J,K+1)-PA(I,J,K))
        WA(I,J,K)=WA(I,J,K)*(1.D0-ALFWA)
     &        +ALFWA*(AEWA(I,J,K)*WA(I+1,J,K)+AWWA(I,J,K)*WA(I-1,J,K)
     &               +ANWA(I,J,K)*WA(I,J+1,K)+ASWA(I,J,K)*WA(I,J-1,K)
     &               +ATWA(I,J,K)*WA(I,J,K+1)+ABWA(I,J,K)*WA(I,J,K-1)
     &               +BWA(I,J,K))/APWA(I,J,K)
 100  CONTINUE
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE SOLVEPC
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NX=40,NY=40,NZ=40)
      COMMON
     & /BVL/UA(0:NX+1,0:NY+1,0:NZ+1),VA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,WA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /PRS/PC(0:NX+1,0:NY+1,0:NZ+1),PA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /UCF/AEUA(0:NX+1,0:NY+1,0:NZ+1),AWUA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ANUA(0:NX+1,0:NY+1,0:NZ+1),ASUA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ATUA(0:NX+1,0:NY+1,0:NZ+1),ABUA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,APUA(0:NX+1,0:NY+1,0:NZ+1),BUA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /VCF/AEVA(0:NX+1,0:NY+1,0:NZ+1),AWVA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ANVA(0:NX+1,0:NY+1,0:NZ+1),ASVA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ATVA(0:NX+1,0:NY+1,0:NZ+1),ABVA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,APVA(0:NX+1,0:NY+1,0:NZ+1),BVA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /WCF/AEWA(0:NX+1,0:NY+1,0:NZ+1),AWWA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ANWA(0:NX+1,0:NY+1,0:NZ+1),ASWA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ATWA(0:NX+1,0:NY+1,0:NZ+1),ABWA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,APWA(0:NX+1,0:NY+1,0:NZ+1),BWA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /PCF/AEPC(0:NX+1,0:NY+1,0:NZ+1),AWPC(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ANPC(0:NX+1,0:NY+1,0:NZ+1),ASPC(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ATPC(0:NX+1,0:NY+1,0:NZ+1),ABPC(0:NX+1,0:NY+1,0:NZ+1)
     &     ,APPC(0:NX+1,0:NY+1,0:NZ+1),BPC(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /CON/DX,DY,DZ,ANU,ALFUA,ALFVA,ALFWA,ALFPA
     & /CNM/NITER,IPTER
C
      DO 100 K=0,NZ+1
      DO 100 J=0,NY+1
      DO 100 I=0,NX+1
        PC(I,J,K)=0.D0
 100  CONTINUE
      DO 110 K=1,NZ
      DO 110 J=1,NY
      DO 110 I=1,NX
        IF(I.EQ.1) THEN
          AWPC(I,J,K)=0.D0
        ELSE
          AWPC(I,J,K)=(DY*DZ)**2/APUA(I-1,J,K)
        ENDIF
        IF(I.EQ.NX) THEN
          AEPC(I,J,K)=0.D0
        ELSE
          AEPC(I,J,K)=(DY*DZ)**2/APUA(I,J,K)
        ENDIF
        IF(J.EQ.1) THEN
          ASPC(I,J,K)=0.D0
        ELSE
          ASPC(I,J,K)=(DX*DZ)**2/APVA(I,J-1,K)
        ENDIF
        IF(J.EQ.NY) THEN
          ANPC(I,J,K)=0.D0
        ELSE
          ANPC(I,J,K)=(DX*DZ)**2/APVA(I,J,K)
        ENDIF
        IF(K.EQ.1) THEN
          ABPC(I,J,K)=0.D0
        ELSE
          ABPC(I,J,K)=(DX*DY)**2/APWA(I,J,K-1)
        ENDIF
        IF(K.EQ.NZ) THEN
          ATPC(I,J,K)=0.D0
        ELSE
          ATPC(I,J,K)=(DX*DY)**2/APWA(I,J,K)
        ENDIF
        APPC(I,J,K)=AEPC(I,J,K)+AWPC(I,J,K)+ANPC(I,J,K)+ASPC(I,J,K)
     &             +ATPC(I,J,K)+ABPC(I,J,K)
        BPC(I,J,K)=-DY*DZ*(UA(I,J,K)-UA(I-1,J,K))
     &             -DX*DZ*(VA(I,J,K)-VA(I,J-1,K))
     &             -DX*DY*(WA(I,J,K)-WA(I,J,K-1))
 110  CONTINUE
      DO 120 IP=1,IPTER
        DO 130 K=1,NZ
        DO 130 J=1,NY
        DO 130 I=1,NX
          PC(I,J,K)=PC(I,J,K)*(1.D0-ALFPA)
     &      +ALFPA*(AEPC(I,J,K)*PC(I+1,J,K)+AWPC(I,J,K)*PC(I-1,J,K)
     &             +ANPC(I,J,K)*PC(I,J+1,K)+ASPC(I,J,K)*PC(I,J-1,K)
     &             +ATPC(I,J,K)*PC(I,J,K+1)+ABPC(I,J,K)*PC(I,J,K-1)
     &             +BPC(I,J,K))/APPC(I,J,K)
 130    CONTINUE
 120  CONTINUE
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE UPDATE
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NX=40,NY=40,NZ=40)
      COMMON
     & /BVL/UA(0:NX+1,0:NY+1,0:NZ+1),VA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,WA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /PRS/PC(0:NX+1,0:NY+1,0:NZ+1),PA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /UCF/AEUA(0:NX+1,0:NY+1,0:NZ+1),AWUA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ANUA(0:NX+1,0:NY+1,0:NZ+1),ASUA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ATUA(0:NX+1,0:NY+1,0:NZ+1),ABUA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,APUA(0:NX+1,0:NY+1,0:NZ+1),BUA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /VCF/AEVA(0:NX+1,0:NY+1,0:NZ+1),AWVA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ANVA(0:NX+1,0:NY+1,0:NZ+1),ASVA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ATVA(0:NX+1,0:NY+1,0:NZ+1),ABVA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,APVA(0:NX+1,0:NY+1,0:NZ+1),BVA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /WCF/AEWA(0:NX+1,0:NY+1,0:NZ+1),AWWA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ANWA(0:NX+1,0:NY+1,0:NZ+1),ASWA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ATWA(0:NX+1,0:NY+1,0:NZ+1),ABWA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,APWA(0:NX+1,0:NY+1,0:NZ+1),BWA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /CON/DX,DY,DZ,ANU,ALFUA,ALFVA,ALFWA,ALFPA
     & /CNM/NITER,IPTER
C
      DO 100 K=1,NZ
      DO 100 J=1,NY
      DO 100 I=1,NX-1
        UA(I,J,K)=UA(I,J,K)
     &           -DY*DZ*(PC(I+1,J,K)-PC(I,J,K))/APUA(I,J,K)
 100  CONTINUE
      DO 110 K=1,NZ
      DO 110 J=1,NY-1
      DO 110 I=1,NX
        VA(I,J,K)=VA(I,J,K)
     &           -DX*DZ*(PC(I,J+1,K)-PC(I,J,K))/APVA(I,J,K)
 110  CONTINUE
      DO 120 K=1,NZ-1
      DO 120 J=1,NY
      DO 120 I=1,NX
        WA(I,J,K)=WA(I,J,K)
     &           -DX*DY*(PC(I,J,K+1)-PC(I,J,K))/APWA(I,J,K)
 120  CONTINUE
      DO 130 K=1,NZ
      DO 130 J=1,NY
      DO 130 I=1,NX
        PA(I,J,K)=PA(I,J,K)+ALFPA*PC(I,J,K)
 130  CONTINUE
      PREF=PA(1,1,1)
      DO 140 K=1,NZ
      DO 140 J=1,NY
      DO 140 I=1,NX
        PA(I,J,K)=PA(I,J,K)-PREF
 140  CONTINUE
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE CONVERG(REST)
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NX=40,NY=40,NZ=40)
      COMMON
     & /BVL/UA(0:NX+1,0:NY+1,0:NZ+1),VA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,WA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /UCF/AEUA(0:NX+1,0:NY+1,0:NZ+1),AWUA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ANUA(0:NX+1,0:NY+1,0:NZ+1),ASUA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ATUA(0:NX+1,0:NY+1,0:NZ+1),ABUA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,APUA(0:NX+1,0:NY+1,0:NZ+1),BUA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /VCF/AEVA(0:NX+1,0:NY+1,0:NZ+1),AWVA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ANVA(0:NX+1,0:NY+1,0:NZ+1),ASVA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ATVA(0:NX+1,0:NY+1,0:NZ+1),ABVA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,APVA(0:NX+1,0:NY+1,0:NZ+1),BVA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /WCF/AEWA(0:NX+1,0:NY+1,0:NZ+1),AWWA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ANWA(0:NX+1,0:NY+1,0:NZ+1),ASWA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,ATWA(0:NX+1,0:NY+1,0:NZ+1),ABWA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,APWA(0:NX+1,0:NY+1,0:NZ+1),BWA(0:NX+1,0:NY+1,0:NZ+1)
C
      RESTUA=0.D0
      RESTVA=0.D0
      RESTWA=0.D0
      DO 100 K=1,NZ
      DO 100 J=1,NY
      DO 100 I=1,NX-1
        RESTUA=RESTUA+DABS(APUA(I,J,K)*UA(I,J,K)-BUA(I,J,K)
     &        -AEUA(I,J,K)*UA(I+1,J,K)-AWUA(I,J,K)*UA(I-1,J,K)
     &        -ANUA(I,J,K)*UA(I,J+1,K)-ASUA(I,J,K)*UA(I,J-1,K)
     &        -ATUA(I,J,K)*UA(I,J,K+1)-ABUA(I,J,K)*UA(I,J,K-1))
 100  CONTINUE
      DO 110 K=1,NZ
      DO 110 J=1,NY-1
      DO 110 I=1,NX
        RESTVA=RESTVA+DABS(APVA(I,J,K)*VA(I,J,K)-BVA(I,J,K)
     &        -AEVA(I,J,K)*VA(I+1,J,K)-AWVA(I,J,K)*VA(I-1,J,K)
     &        -ANVA(I,J,K)*VA(I,J+1,K)-ASVA(I,J,K)*VA(I,J-1,K)
     &        -ATVA(I,J,K)*VA(I,J,K+1)-ABVA(I,J,K)*VA(I,J,K-1))
 110  CONTINUE
      DO 120 K=1,NZ-1
      DO 120 J=1,NY
      DO 120 I=1,NX
        RESTWA=RESTWA+DABS(APWA(I,J,K)*WA(I,J,K)-BWA(I,J,K)
     &        -AEWA(I,J,K)*WA(I+1,J,K)-AWWA(I,J,K)*WA(I-1,J,K)
     &        -ANWA(I,J,K)*WA(I,J+1,K)-ASWA(I,J,K)*WA(I,J-1,K)
     &        -ATWA(I,J,K)*WA(I,J,K+1)-ABWA(I,J,K)*WA(I,J,K-1))
 120  CONTINUE
      WRITE(6,'( 1H ,1P,3E14.4)') RESTUA,RESTVA,RESTWA
      REST=DMAX1(RESTUA,RESTVA)
      REST=DMAX1(REST,RESTWA)
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE ENPR
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NX=40,NY=40,NZ=40)
      COMMON
     & /BVL/UA(0:NX+1,0:NY+1,0:NZ+1),VA(0:NX+1,0:NY+1,0:NZ+1)
     &     ,WA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /PRS/PC(0:NX+1,0:NY+1,0:NZ+1),PA(0:NX+1,0:NY+1,0:NZ+1)
      COMMON
     & /XGD/X(0:NX),XX(0:NX)
     & /YGD/Y(0:NY),YY(0:NY)
     & /ZGD/Z(0:NZ),ZZ(0:NZ)
C
      OPEN(12,FILE='DATA-XY.txt',FORM='FORMATTED')
      OPEN(13,FILE='DATA-YZ.txt',FORM='FORMATTED')
      OPEN(14,FILE='DATA-ZX.txt',FORM='FORMATTED')
      K=NZ/2
      DO 100 J=1,NY
      DO 100 I=1,NX
        WRITE(12,'( 1H ,1P,6E14.4 )') X(I),Y(J)
     &  ,0.25D0*(UA(I,J,K)+UA(I-1,J,K)+UA(I,J,K+1)+UA(I-1,J,K+1))
     &  ,0.25D0*(VA(I,J,K)+VA(I,J-1,K)+VA(I,J,K+1)+VA(I,J-1,K+1))
     &  ,WA(I,J,K),0.5D0*(PA(I,J,K+1)+PA(I,J,K))
 100  CONTINUE
      I=NX/2
      DO 110 K=1,NZ
      DO 110 J=1,NY
        WRITE(13,'( 1H ,1P,6E14.4 )') Z(K),Y(J),UA(I,J,K)
     &  ,0.25D0*(WA(I,J,K)+VA(I,J-1,K)+WA(I+1,J,K)+VA(I+1,J-1,K))
     &  ,0.25D0*(WA(I,J,K)+WA(I,J,K-1)+WA(I+1,J,K)+WA(I+1,J,K-1))
     &  ,0.5D0*(PA(I+1,J,K)+PA(I,J,K))
 110  CONTINUE
      J=NY/2
      DO 120 K=1,NZ
      DO 120 I=1,NX
        WRITE(14,'( 1H ,1P,6E14.4 )') X(I),Z(K)
     &  ,0.25D0*(UA(I,J,K)+UA(I-1,J,K)+UA(I,J+1,K)+UA(I-1,J+1,K))
     &  ,VA(I,J,K)
     &  ,0.25D0*(WA(I,J,K)+WA(I,J,K-1)+WA(I,J+1,K)+WA(I,J+1,K-1))
     &  ,0.5D0*(PA(I,J+1,K)+PA(I,J,K))
 120  CONTINUE
      CLOSE(12)
      CLOSE(13)
      CLOSE(14)
      RETURN
      END
