C
C********* ********* ********* ********* ********* ********* *********C
C
C     Karman-Trefftz.f
C
C     WRITTEN BY MASAYOSHI OKAMOTO
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NT=720,NNN=200)
      COMPLEX*16 Z(0:NT+1),TETA(0:NT+1),ZETA(0:NT+1)
      DIMENSION X(0:NT+1),Y(0:NT+1)
      DIMENSION THETA(0:NT+1)
      COMPLEX*16 TETA0
      DIMENSION XX(0:NNN),Y1(0:NNN),Y2(0:NNN),CAM(0:NNN),THI(0:NNN)
C
      PI=4.D0*DATAN(1.D0)
      WRITE(6,*) 'ITYPE:0. Biconvex, 1.Symmetry KT wing, 2.KT wing'
      READ(5,*) ITYPE
      IF(ITYPE.EQ.0) THEN
        R=1.D0
      ELSE
        WRITE(6,*) 'R    : Radius rate of zeta circle R>1'
        READ(5,*) R
      ENDIF
      WRITE(6,*) 'K    : KT parameter 2>K>1'
      READ(5,*) TK
      BETA=0.D0
      IF(ITYPE.EQ.2) THEN
        WRITE(6,*) 'BETA : Center angle of zeta circle  (degree) '
        READ(5,*) BETA
        BETA=BETA/360.D0*2.D0*PI
      ENDIF
      A=1.D0
      TETA0=A-(1.D0,0.D0)*R*DCOS(BETA)+(0.D0,1.D0)*R*DSIN(BETA)
C
      DO 100 I=0,NT
        THETA(I)=2.D0*PI*(DFLOAT(I)+0.5D0)/DFLOAT(NT)
  100 CONTINUE
      DO 110 I=0,NT
        TETA(I)=R*DCOS(THETA(I))+(0.0,1.D0)*R*DSIN(THETA(I))
        ZETA(I)=TETA(I)+TETA0
C
C 8888888888    KT Transformation    8888888888888
C
        Z(I)=TK*A*((ZETA(I)+A)**TK+(ZETA(I)-A)**TK)
     &           /((ZETA(I)+A)**TK-(ZETA(I)-A)**TK)
C
C 888888888888888888888888888888888888888888888888
C
        X(I)=DREAL(Z(I))
        Y(I)=DIMAG(Z(I))
  110 CONTINUE
      XMAX=0.D0
      XMIN=0.D0
      DO 120 I=0,NT
        IF(X(I).LT.XMIN) THEN
          XMIN=X(I)
          IMIN=I
        ENDIF
        IF(X(I).GT.XMAX) THEN
          XMAX=X(I)
          IMAX=I
        ENDIF
  120 CONTINUE
      CHORD=XMAX-XMIN
      DX=CHORD/DFLOAT(NNN)
      DO 130 II=0,NNN
        XX(II)=XMIN+DX*DFLOAT(II)
  130 CONTINUE
      DO 140 II=1,NNN-1
        DO 150 I=0,NT
        IF((X(I).LT.XX(II)).AND.(X(I+1).GT.XX(II))) THEN
          I1=I
          GOTO 1000
        ELSEIF((X(I).GT.XX(II)).AND.(X(I+1).LT.XX(II))) THEN
          I1=I
          GOTO 1000
        ENDIF
  150   CONTINUE
 1000   CONTINUE
        DO 180 I=I1+1,NT
        IF((X(I).LT.XX(II)).AND.(X(I+1).GT.XX(II))) THEN
          I2=I
          GOTO 1100
        ELSEIF((X(I).GT.XX(II)).AND.(X(I+1).LT.XX(II))) THEN
          I2=I
          GOTO 1100
        ENDIF
  180   CONTINUE
 1100   CONTINUE
        Y1(II)=(Y(I1)*(X(I1+1)-XX(II))+Y(I1+1)*(XX(II)-X(I1)))
     &        /(X(I1+1)-X(I1))
        Y2(II)=(Y(I2)*(X(I2+1)-XX(II))+Y(I2+1)*(XX(II)-X(I2)))
     &        /(X(I2+1)-X(I2))
        CAM(II)=(Y1(II)+Y2(II))*0.5D0
        THI(II)=DABS(Y1(II)-Y2(II))
  140 CONTINUE
C
      OPEN(11,FILE='KT-CIRCLE.txt',FORM='FORMATTED')
      OPEN(12,FILE='KT-WING.txt',FORM='FORMATTED')
      OPEN(13,FILE='KT-CAMBER.txt',FORM='FORMATTED')
      OPEN(14,FILE='KT-THICK.txt',FORM='FORMATTED')
      OPEN(15,FILE='KT-CHORD.txt',FORM='FORMATTED')
      DO 160 I=0,NT
        WRITE(11,*) DREAL(ZETA(I)),DIMAG(ZETA(I))
        WRITE(12,*) X(I),Y(I)
  160 CONTINUE
      DO 170 II=0,NNN
        WRITE(13,*) XX(II),CAM(II)
        WRITE(14,*) XX(II),THI(II)
  170 CONTINUE
      WRITE(15,*) 'Chord length',CHORD
      CLOSE(11)
      CLOSE(12)
      CLOSE(13)
      CLOSE(14)
      CLOSE(15)
      END
