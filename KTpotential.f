C
C********* ********* ********* ********* ********* ********* *********C
C
C     KTpotential.f
C
C     WRITTEN BY MASAYOSHI OKAMOTO
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NR=50,NT=72)
      COMPLEX*16 F(0:NT+1,0:NR+1),W(0:NT+1,0:NR+1),Z(0:NT+1,0:NR+1)
     &          ,TETA(0:NT+1,0:NR+1),ZETA(0:NT+1,0:NR+1)
     &          ,ZETA0(0:NT+1,0:NR+1),DZ(0:NT+1,0:NR+1)
      DIMENSION
     &  POR(0:NT+1,0:NR+1),PHR(0:NT+1,0:NR+1),X(0:NT+1,0:NR+1)
     & ,Y(0:NT+1,0:NR+1),UR(0:NT+1,0:NR+1),VR(0:NT+1,0:NR+1)
     & ,PR(0:NT+1,0:NR+1),XZ(0:NT+1,0:NR+1),YZ(0:NT+1,0:NR+1)
      DIMENSION THETA(0:NT+1)
      DIMENSION RAMDA(0:NR+1),RAMDA0(0:NR+1)
      COMPLEX*16 DALM,DALP,TETA0,AE,AP,AM,CAP,CAM,CC,CAPM
C
      PI=4.D0*DATAN(1.D0)
      R1=1.D0
      WRITE(6,*) 'U    : Attack flow velocity'
      READ(5,*) U
      WRITE(6,*) 'ALPHA: Attack angle (degree)'
      READ(5,*) ALPHA
      ALPHA=ALPHA/360.D0*2.D0*PI
      WRITE(6,*) 'K     : KT parameter 2>K>1'
      READ(5,*) TK
      WRITE(6,*) 'ITYPE:0.Biconvex'
      WRITE(6,*) 'ITYPE:1.Symmetric KT wing'
      WRITE(6,*) 'ITYPE:2.Asymmetric KT wing'
      READ(5,*) ITYPE
      IF(ITYPE.EQ.0) THEN
        R=1.D0
      ELSE
        WRITE(6,*) 'R    : Radius rate of zeta circle R>1'
        READ(5,*) R
      ENDIF
      BETA=0.D0
      IF(ITYPE.EQ.2) THEN
        WRITE(6,*) 'BETA : Central angle of Zeta circle  (degree) '
        READ(5,*) BETA
        BETA=BETA/360.D0*2.D0*PI
      ENDIF
      IF(ITYPE.EQ.0) THEN
        TETA0=(0.D0,0.D0)
      ELSEIF(ITYPE.EQ.1) THEN
        TETA0=-(1.D0,0.D0)*(R-R1)
      ELSEIF(ITYPE.EQ.2) THEN
        TETA0=R1-(1.D0,0.D0)*R*DCOS(BETA)
     &          +(0.D0,1.D0)*R*DSIN(BETA)
      ENDIF
      WRITE(6,*) 'IJK  : 0.with JK assumption'
      WRITE(6,*) 'IJK  : 1.without JK assumption'
      READ(5,*) IJK
      IF(IJK.EQ.1) THEN
        WRITE(6,*) 'GAMMA: Circulation'
        READ(5,*) GAMMA
      ENDIF
      IF(IJK.EQ.0) THEN
        GAMMA=4.D0*PI*U*R*DSIN(ALPHA+BETA)
      ENDIF
C
      DR=5.D0/DFLOAT(NR)
      DO 100 I=0,NT
        THETA(I)=2.D0*PI*DFLOAT(I)/DFLOAT(NT)
  100 CONTINUE
      DO 110 J=0,NR
        RAMDA(J)=R1+DR*DFLOAT(J)
  110 CONTINUE
      AE=R*DCOS(ALPHA)+(0.D0,1.D0)*R*DSIN(ALPHA)
      DO 120 J=0,NR
      DO 120 I=0,NT
        ZETA0(I,J)=RAMDA(J)*DCOS(THETA(I))
     & +(0.D0,1.D0)*RAMDA(J)*DSIN(THETA(I))
        F(I,J)=U*(ZETA0(I,J)+R1**2/ZETA0(I,J))
     &        +(0.D0,1.D0)*GAMMA/2.D0/PI*CDLOG(ZETA0(I,J))
        TETA(I,J)=R*RAMDA(J)*DCOS(THETA(I)+ALPHA)
     &    +(0.D0,1.D0)*R*RAMDA(J)*DSIN(THETA(I)+ALPHA)
        ZETA(I,J)=TETA(I,J)+TETA0
        XZ(I,J)=DREAL(ZETA(I,J))
        YZ(I,J)=DIMAG(ZETA(I,J))
        AP=ZETA(I,J)+R1*(1.D0,0.D0)
        AM=ZETA(I,J)-R1*(1.D0,0.D0)
        RAP=CDABS(AP)
        RAM=CDABS(AM)
        TAP=DACOS(DREAL(AP)/RAP)
        TAM=DACOS(DREAL(AM)/RAM)
        IF(DIMAG(AP).LT.0.D0) TAP=2.D0*PI-TAP
        IF(DIMAG(AM).LT.0.D0) TAM=2.D0*PI-TAM
        CAP=RAP**TK*((1.D0,0.D0)*DCOS(TAP*TK)
     &              +(0.D0,1.D0)*DSIN(TAP*TK))
        CAM=RAM**TK*((1.D0,0.D0)*DCOS(TAM*TK)
     &              +(0.D0,1.D0)*DSIN(TAM*TK))
        CAPM=(RAP*RAM)**(TK-1.D0)
     &      *((1.D0,0.D0)*DCOS((TAP+TAM)*(TK-1.D0))
     &       +(0.D0,1.D0)*DSIN((TAP+TAM)*(TK-1.D0)))
        CC=CAP*CAP-2.D0*CAP*CAM+CAM*CAM
        Z(I,J)=TK*R*(CAP+CAM)/(CAP-CAM)
        DZ(I,J)=4.D0*(TK*R1)**2*AE*CAPM/CC
        X(I,J)=DREAL(Z(I,J))
        Y(I,J)=DIMAG(Z(I,J))
        POR(I,J)=DREAL(F(I,J))
        PHR(I,J)=DIMAG(F(I,J))
        W(I,J)=(U*(1.D0-R1**2/ZETA0(I,J)**2)
     &        +(0.D0,1.D0)*GAMMA/2.D0/PI/ZETA(I,J))/DZ(I,J)
        UR(I,J)=DREAL(W(I,J))
        VR(I,J)=-DIMAG(W(I,J))
        PR(I,J)=-0.5D0*(UR(I,J)**2+VR(I,J)**2)
  120 CONTINUE
      J=0
      XMAX=0.D0
      XMIN=0.D0
      DO 130 I=0,NT
        IF(X(I,J).LT.XMIN) THEN
          XMIN=X(I,J)
        ENDIF
        IF(X(I,J).GT.XMAX) THEN
          XMAX=X(I,J)
        ENDIF
  130 CONTINUE
      CHORD=XMAX-XMIN
      OPEN(15,FILE='WING.txt',FORM='FORMATTED')
      OPEN(16,FILE='FORCE.txt',FORM='FORMATTED')
      OPEN(17,FILE='CIRCLE.txt',FORM='FORMATTED')
      OPEN(20,FILE='STREAMLINE.txt',FORM='FORMATTED')
      OPEN(21,FILE='PRESSURE.txt',FORM='FORMATTED')
      OPEN(22,FILE='VELOCITY.txt',FORM='FORMATTED')
      J=0
      DO 140 I=0,NT
        WRITE(15,'( 1H ,1P,2E14.4 )') X(I,J),Y(I,J)
        WRITE(17,'( 1H ,1P,2E14.4 )') XZ(I,J),YZ(I,J)
  140 CONTINUE
      WRITE(16,*) 'Velocity',U
      WRITE(16,*) 'Angle',ALPHA
      WRITE(16,*) 'Circulation',GAMMA
      WRITE(16,*) 'Fy',GAMMA*U
      WRITE(16,*) 'Fx',0.D0
      WRITE(16,*) 'CL',2.D0*GAMMA/U/CHORD
      WRITE(16,*) 'CD',0.D0
      CLOSE(15)
      CLOSE(16)
      CLOSE(17)
      DO 150 J=0,NR
      DO 150 I=0,NT
        WRITE(20,'( 1H ,1P,3E14.4 )') X(I,J),Y(I,J),PHR(I,J)
        IF(PR(I,J).GT.-1000.D0) THEN
          WRITE(21,'( 1H ,1P,3E14.4 )') X(I,J),Y(I,J),PR(I,J)
        ENDIF
  150 CONTINUE
      DO 160 J=2,NR,4
      DO 160 I=2,NT,4
        WRITE(22,'( 1H ,1P,4E14.4 )') X(I,J),Y(I,J),UR(I,J),VR(I,J)
  160 CONTINUE
      CLOSE(20)
      CLOSE(21)
      CLOSE(22)
      END
