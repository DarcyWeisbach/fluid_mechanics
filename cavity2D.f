C
C********* ********* ********* ********* ********* ********* *********C
C
C     PROGRAM cavity2D.f
C
C     Program for 2D cavity flow by FVM
C     MADE BY MASAYOSHI OKAMOTO
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NX=40,NY=40)
      COMMON
     & /BVL/UA(0:NX+1,0:NY+1),VA(0:NX+1,0:NY+1)
      COMMON
     & /PRS/PC(0:NX+1,0:NY+1),PA(0:NX+1,0:NY+1)
      COMMON
     & /UCF/APUA(0:NX+1,0:NY+1),AEUA(0:NX+1,0:NY+1),AWUA(0:NX+1,0:NY+1)
     &     ,ANUA(0:NX+1,0:NY+1),ASUA(0:NX+1,0:NY+1),BUA(0:NX+1,0:NY+1)
      COMMON
     & /VCF/APVA(0:NX+1,0:NY+1),AEVA(0:NX+1,0:NY+1),AWVA(0:NX+1,0:NY+1)
     &     ,ANVA(0:NX+1,0:NY+1),ASVA(0:NX+1,0:NY+1),BVA(0:NX+1,0:NY+1)
      COMMON
     & /PCF/APPC(0:NX+1,0:NY+1),AEPC(0:NX+1,0:NY+1),AWPC(0:NX+1,0:NY+1)
     &     ,ANPC(0:NX+1,0:NY+1),ASPC(0:NX+1,0:NY+1),BPC(0:NX+1,0:NY+1)
      COMMON
     & /XGD/X(0:NX),XX(0:NX)
     & /YGD/Y(0:NY),YY(0:NY)
      COMMON
     & /CON/DX,DY,ANU,ALFUA,ALFVA,ALFPA
     & /CNM/NITER,IPTER
C
      CALL INIT
      CALL BDCD
      DO 100 ITER=1,NITER
        CALL SOLVEUA
        CALL SOLVEVA
        CALL SOLVEPC
        CALL UPDATE
        IF(MOD(ITER,1000).EQ.0) THEN
          WRITE(6,*) ITER
          CALL CONVERG(REST)
          IF(REST.LT.1.D-5) GOTO 110
        ENDIF
 100  CONTINUE
 110  CALL ENPR
      WRITE(11) UA,VA,PA
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
     & /BVL/UA(0:NX+1,0:NY+1),VA(0:NX+1,0:NY+1)
      COMMON
     & /PRS/PC(0:NX+1,0:NY+1),PA(0:NX+1,0:NY+1)
      COMMON
     & /UCF/APUA(0:NX+1,0:NY+1),AEUA(0:NX+1,0:NY+1),AWUA(0:NX+1,0:NY+1)
     &     ,ANUA(0:NX+1,0:NY+1),ASUA(0:NX+1,0:NY+1),BUA(0:NX+1,0:NY+1)
      COMMON
     & /VCF/APVA(0:NX+1,0:NY+1),AEVA(0:NX+1,0:NY+1),AWVA(0:NX+1,0:NY+1)
     &     ,ANVA(0:NX+1,0:NY+1),ASVA(0:NX+1,0:NY+1),BVA(0:NX+1,0:NY+1)
      COMMON
     & /PCF/APPC(0:NX+1,0:NY+1),AEPC(0:NX+1,0:NY+1),AWPC(0:NX+1,0:NY+1)
     &     ,ANPC(0:NX+1,0:NY+1),ASPC(0:NX+1,0:NY+1),BPC(0:NX+1,0:NY+1)
      COMMON
     & /XGD/X(0:NX),XX(0:NX)
     & /YGD/Y(0:NY),YY(0:NY)
      COMMON
     & /CON/DX,DY,ANU,ALFUA,ALFVA,ALFPA
     & /CNM/NITER,IPTER
C
      DO 100 J=0,NY+1
      DO 100 I=0,NX+1
        UA(I,J)=0.D0
        APUA(I,J)=0.D0
        AEUA(I,J)=0.D0
        AWUA(I,J)=0.D0
        ANUA(I,J)=0.D0
        ASUA(I,J)=0.D0
        BUA(I,J)=0.D0
        VA(I,J)=0.D0
        APVA(I,J)=0.D0
        AEVA(I,J)=0.D0
        AWVA(I,J)=0.D0
        ANVA(I,J)=0.D0
        ASVA(I,J)=0.D0
        BVA(I,J)=0.D0
        PA(I,J)=0.D0
        PC(I,J)=0.D0
        APPC(I,J)=0.D0
        AEPC(I,J)=0.D0
        AWPC(I,J)=0.D0
        ANPC(I,J)=0.D0
        ASPC(I,J)=0.D0
        BPC(I,J)=0.D0
 100  CONTINUE
      DO 110 I=0,NX+1
        X(I)=0.D0
        XX(I)=0.D0
 110  CONTINUE
      DO 120 J=0,NY+1
        Y(J)=0.D0
        YY(J)=0.D0
 120  CONTINUE
      ALFUA=0.1D0
      ALFVA=0.1D0
      ALFPA=0.1D0
      RE=100.D0
      ANU=1.D0/RE
      NITER=100000
      IPTER=5
      DX=1.D0/DFLOAT(NX)
      DY=1.D0/DFLOAT(NY)
      DO 130 I=0,NX
        XX(I)=DX*DFLOAT(I)
        YY(I)=DY*DFLOAT(I)
 130  CONTINUE
      DO 140 I=1,NY
        X(I)=XX(I)-0.5D0*DX
        Y(I)=YY(I)-0.5D0*DY
 140  CONTINUE
C      READ(11) UA,VA,PA
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
     & /BVL/UA(0:NX+1,0:NY+1),VA(0:NX+1,0:NY+1)
C
      DO 100 I=0,NX+1
        UA(I,NY+1)=1.D0
        VA(I,NY)=0.D0
        UA(I,0)=0.D0
        VA(I,0)=0.D0
 100  CONTINUE
      DO 110 J=0,NY+1
        UA(NX,J)=0.D0
        VA(NX+1,J)=0.D0
        UA(0,J)=0.D0
        VA(0,J)=0.D0
 110  CONTINUE
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
      PARAMETER(NX=40,NY=40)
      COMMON
     & /BVL/UA(0:NX+1,0:NY+1),VA(0:NX+1,0:NY+1)
      COMMON
     & /PRS/PC(0:NX+1,0:NY+1),PA(0:NX+1,0:NY+1)
      COMMON
     & /UCF/APUA(0:NX+1,0:NY+1),AEUA(0:NX+1,0:NY+1),AWUA(0:NX+1,0:NY+1)
     &     ,ANUA(0:NX+1,0:NY+1),ASUA(0:NX+1,0:NY+1),BUA(0:NX+1,0:NY+1)
      COMMON
     & /CON/DX,DY,ANU,ALFUA,ALFVA,ALFPA
     & /CNM/NITER,IPTER
C
      DO 100 J=1,NY
      DO 100 I=1,NX-1
        AEUA(I,J)=ANU*DY/DX-0.25D0*DY*(UA(I+1,J)+UA(I,J))
        AWUA(I,J)=ANU*DY/DX+0.25D0*DY*(UA(I,J)+UA(I-1,J))
        IF(J.EQ.1) THEN
          ASUA(I,J)=2.D0*ANU*DX/DY
        ELSE
          ASUA(I,J)=ANU*DX/DY+0.25D0*DX*(VA(I+1,J-1)+VA(I,J-1))
        ENDIF
        IF(J.EQ.NY) THEN
          ANUA(I,J)=2.D0*ANU*DX/DY
        ELSE
          ANUA(I,J)=ANU*DX/DY-0.25D0*DX*(VA(I+1,J)+VA(I,J))
        ENDIF
        APUA(I,J)=AEUA(I,J)+AWUA(I,J)+ANUA(I,J)+ASUA(I,J)
        BUA(I,J)=-DY*(PA(I+1,J)-PA(I,J))
        UA(I,J)=UA(I,J)*(1.D0-ALFUA)+ALFUA*(AEUA(I,J)*UA(I+1,J)
     &         +AWUA(I,J)*UA(I-1,J)+ANUA(I,J)*UA(I,J+1)
     &         +ASUA(I,J)*UA(I,J-1)+BUA(I,J))/APUA(I,J)
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
      PARAMETER(NX=40,NY=40)
      COMMON
     & /BVL/UA(0:NX+1,0:NY+1),VA(0:NX+1,0:NY+1)
      COMMON
     & /PRS/PC(0:NX+1,0:NY+1),PA(0:NX+1,0:NY+1)
      COMMON
     & /VCF/APVA(0:NX+1,0:NY+1),AEVA(0:NX+1,0:NY+1),AWVA(0:NX+1,0:NY+1)
     &     ,ANVA(0:NX+1,0:NY+1),ASVA(0:NX+1,0:NY+1),BVA(0:NX+1,0:NY+1)
      COMMON
     & /CON/DX,DY,ANU,ALFUA,ALFVA,ALFPA
     & /CNM/NITER,IPTER
C
      DO 100 J=1,NY-1
      DO 100 I=1,NX
        IF(I.EQ.1) THEN
          AWVA(I,J)=2.0D0*ANU*DY/DX
        ELSE
          AWVA(I,J)=ANU*DY/DX+0.25D0*DY*(UA(I-1,J+1)+UA(I-1,J))
        ENDIF
        IF(I.EQ.NX) THEN
          AEVA(I,J)=2.D0*ANU*DY/DX
        ELSE
          AEVA(I,J)=ANU*DY/DX-0.25D0*DY*(UA(I,J+1)+UA(I,J))
        ENDIF
        ANVA(I,J)=ANU*DX/DY-0.25D0*DX*(VA(I,J+1)+VA(I,J))
        ASVA(I,J)=ANU*DX/DY+0.25D0*DX*(VA(I,J)+VA(I,J-1))
        APVA(I,J)=AEVA(I,J)+AWVA(I,J)+ANVA(I,J)+ASVA(I,J)
        BVA(I,J)=-DX*(PA(I,J+1)-PA(I,J))
        VA(I,J)=VA(I,J)*(1.D0-ALFVA)+ALFVA*(AEVA(I,J)*VA(I+1,J)
     &         +AWVA(I,J)*VA(I-1,J)+ANVA(I,J)*VA(I,J+1)
     &         +ASVA(I,J)*VA(I,J-1)+BVA(I,J))/APVA(I,J)
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
      PARAMETER(NX=40,NY=40)
      COMMON
     & /BVL/UA(0:NX+1,0:NY+1),VA(0:NX+1,0:NY+1)
      COMMON
     & /PRS/PC(0:NX+1,0:NY+1),PA(0:NX+1,0:NY+1)
      COMMON
     & /UCF/APUA(0:NX+1,0:NY+1),AEUA(0:NX+1,0:NY+1),AWUA(0:NX+1,0:NY+1)
     &     ,ANUA(0:NX+1,0:NY+1),ASUA(0:NX+1,0:NY+1),BUA(0:NX+1,0:NY+1)
      COMMON
     & /VCF/APVA(0:NX+1,0:NY+1),AEVA(0:NX+1,0:NY+1),AWVA(0:NX+1,0:NY+1)
     &     ,ANVA(0:NX+1,0:NY+1),ASVA(0:NX+1,0:NY+1),BVA(0:NX+1,0:NY+1)
      COMMON
     & /PCF/APPC(0:NX+1,0:NY+1),AEPC(0:NX+1,0:NY+1),AWPC(0:NX+1,0:NY+1)
     &     ,ANPC(0:NX+1,0:NY+1),ASPC(0:NX+1,0:NY+1),BPC(0:NX+1,0:NY+1)
      COMMON
     & /CON/DX,DY,ANU,ALFUA,ALFVA,ALFPA
     & /CNM/NITER,IPTER
C
      DO 100 J=1,NY
      DO 100 I=1,NX
        PC(I,J)=0.D0
 100  CONTINUE
      DO 110 J=1,NY
      DO 110 I=1,NX
        IF(I.EQ.1) THEN
          AWPC(I,J)=0.D0
        ELSE
          AWPC(I,J)=DY*DY/APUA(I-1,J)
        ENDIF
        IF(I.EQ.NX) THEN
          AEPC(I,J)=0.D0
        ELSE
          AEPC(I,J)=DY*DY/APUA(I,J)
        ENDIF
        IF(J.EQ.1) THEN
          ASPC(I,J)=0.D0
        ELSE
          ASPC(I,J)=DX*DX/APVA(I,J-1)
        ENDIF
        IF(J.EQ.NY) THEN
          ANPC(I,J)=0.D0
        ELSE
          ANPC(I,J)=DX*DX/APVA(I,J)
        ENDIF
        APPC(I,J)=AEPC(I,J)+AWPC(I,J)+ANPC(I,J)+ASPC(I,J)
        BPC(I,J)=-DY*(UA(I,J)-UA(I-1,J))-DX*(VA(I,J)-VA(I,J-1))
 110  CONTINUE
      DO 120 IP=1,IPTER
        DO 130 J=1,NY
        DO 130 I=1,NX
          PC(I,J)=PC(I,J)*(1.D0-ALFPA)+ALFPA*(AEPC(I,J)*PC(I+1,J)
     &           +AWPC(I,J)*PC(I-1,J)+ANPC(I,J)*PC(I,J+1)
     &           +ASPC(I,J)*PC(I,J-1)+BPC(I,J))/APPC(I,J)
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
      PARAMETER(NX=40,NY=40)
      COMMON
     & /BVL/UA(0:NX+1,0:NY+1),VA(0:NX+1,0:NY+1)
      COMMON
     & /PRS/PC(0:NX+1,0:NY+1),PA(0:NX+1,0:NY+1)
      COMMON
     & /UCF/APUA(0:NX+1,0:NY+1),AEUA(0:NX+1,0:NY+1),AWUA(0:NX+1,0:NY+1)
     &     ,ANUA(0:NX+1,0:NY+1),ASUA(0:NX+1,0:NY+1),BUA(0:NX+1,0:NY+1)
      COMMON
     & /VCF/APVA(0:NX+1,0:NY+1),AEVA(0:NX+1,0:NY+1),AWVA(0:NX+1,0:NY+1)
     &     ,ANVA(0:NX+1,0:NY+1),ASVA(0:NX+1,0:NY+1),BVA(0:NX+1,0:NY+1)
      COMMON
     & /PCF/APPC(0:NX+1,0:NY+1),AEPC(0:NX+1,0:NY+1),AWPC(0:NX+1,0:NY+1)
     &     ,ANPC(0:NX+1,0:NY+1),ASPC(0:NX+1,0:NY+1),BPC(0:NX+1,0:NY+1)
      COMMON
     & /CON/DX,DY,ANU,ALFUA,ALFVA,ALFPA
     & /CNM/NITER,IPTER
C
      DO 100 J=1,NY
      DO 100 I=1,NX-1
        UA(I,J)=UA(I,J)-DY*(PC(I+1,J)-PC(I,J))/APUA(I,J)
 100  CONTINUE
      DO 110 J=1,NY-1
      DO 110 I=1,NX
        VA(I,J)=VA(I,J)-DX*(PC(I,J+1)-PC(I,J))/APVA(I,J)
 110  CONTINUE
      DO 120 J=1,NY
      DO 120 I=1,NX
        PA(I,J)=PA(I,J)+ALFPA*PC(I,J)
 120  CONTINUE
      PREF=PA(1,1)
      DO 130 J=1,NY
      DO 130 I=1,NX
        PA(I,J)=PA(I,J)-PREF
 130  CONTINUE
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
      PARAMETER(NX=40,NY=40)
      COMMON
     & /BVL/UA(0:NX+1,0:NY+1),VA(0:NX+1,0:NY+1)
      COMMON
     & /PRS/PC(0:NX+1,0:NY+1),PA(0:NX+1,0:NY+1)
      COMMON
     & /UCF/APUA(0:NX+1,0:NY+1),AEUA(0:NX+1,0:NY+1),AWUA(0:NX+1,0:NY+1)
     &     ,ANUA(0:NX+1,0:NY+1),ASUA(0:NX+1,0:NY+1),BUA(0:NX+1,0:NY+1)
      COMMON
     & /VCF/APVA(0:NX+1,0:NY+1),AEVA(0:NX+1,0:NY+1),AWVA(0:NX+1,0:NY+1)
     &     ,ANVA(0:NX+1,0:NY+1),ASVA(0:NX+1,0:NY+1),BVA(0:NX+1,0:NY+1)
C
      RESTUA=0.D0
      RESTVA=0.D0
      DO 100 J=1,NY
      DO 100 I=1,NX-1
        RESTUA=RESTUA+DABS(APUA(I,J)*UA(I,J)-BUA(I,J)
     &        -AEUA(I,J)*UA(I+1,J)-AWUA(I,J)*UA(I-1,J)
     &        -ANUA(I,J)*UA(I,J+1)-ASUA(I,J)*UA(I,J-1))
 100  CONTINUE
      DO 110 J=1,NY-1
      DO 110 I=1,NX
        RESTVA=RESTVA+DABS(APVA(I,J)*VA(I,J)-BVA(I,J)
     &        -AEVA(I,J)*VA(I+1,J)-AWVA(I,J)*VA(I-1,J)
     &        -ANVA(I,J)*VA(I,J+1)-ASVA(I,J)*VA(I,J-1))
 110  CONTINUE
      WRITE(6,'( 1H ,1P,2E14.4)') RESTUA,RESTVA
      REST=DMAX1(RESTUA,RESTVA)
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
      PARAMETER(NX=40,NY=40)
      COMMON
     & /BVL/UA(0:NX+1,0:NY+1),VA(0:NX+1,0:NY+1)
      COMMON
     & /PRS/PC(0:NX+1,0:NY+1),PA(0:NX+1,0:NY+1)
      COMMON
     & /XGD/X(0:NX),XX(0:NX)
     & /YGD/Y(0:NY),YY(0:NY)
      COMMON
     & /CON/DX,DY,ANU,ALFUA,ALFVA,ALFPA
     & /CNM/NITER,IPTER
C
      OPEN(12,FILE='DATA.txt',FORM='FORMATTED')
      DO 100 J=1,NY
      DO 100 I=1,NX
        WRITE(12,'( 1H ,1P,5E14.4 )') X(I),Y(J)
     &  ,0.5D0*(UA(I,J)+UA(I-1,J)),0.5D0*(VA(I,J)+VA(I,J-1))
     &  ,PA(I,J)
 100  CONTINUE
      CLOSE(12)
      RETURN
      END
