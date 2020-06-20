C
C********* ********* ********* ********* ********* ********* *********C
C
C     Joukowski.f
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
      WRITE(6,*) 'ITYPE:1.Flat plate wing, 2.Elliptical wing 3.Arc 
     & wing, 4.symmetry JK wing, 5.JK wing'
      READ(5,*) ITYPE
      IF(ITYPE.NE.1) THEN
        WRITE(6,*) 'R    : Radius rate of zeta circle R>1'
        READ(5,*) R
      ENDIF
      IF(ITYPE.EQ.5) THEN
        WRITE(6,*) 'BETA : Center angle of zeta circle  (degree) '
        READ(5,*) BETA
        BETA=BETA/360.D0*2.D0*PI
      ENDIF
      A=1.D0
      IF(ITYPE.EQ.1) THEN
        TETA0=(0.D0,0.D0)
        R=1.D0
        BETA=0.D0
      ELSEIF(ITYPE.EQ.2) THEN
        TETA0=(0.D0,0.D0)
        BETA=0.D0
      ELSEIF(ITYPE.EQ.3) THEN
        TETA0=(0.D0,1.D0)*DSQRT(R**2-A**2)
        BETA=DACOS(A/R)
      ELSEIF(ITYPE.EQ.4) THEN
        TETA0=-(1.D0,0.D0)*(R-A)
        BETA=0.D0
      ELSEIF(ITYPE.EQ.5) THEN
        TETA0=A-(1.D0,0.D0)*R*DCOS(BETA)
     &         +(0.D0,1.D0)*R*DSIN(BETA)
      ENDIF
C
      DO 100 I=0,NT
        THETA(I)=2.D0*PI*(DFLOAT(I)+0.5D0)/DFLOAT(NT)
  100 CONTINUE
      DO 110 I=0,NT
        TETA(I)=R*DCOS(THETA(I))+(0.0,1.D0)*R*DSIN(THETA(I))
        ZETA(I)=TETA(I)+TETA0
C
C 8888888888    JK Transformation    8888888888888
C
        Z(I)=(ZETA(I)+A**2/ZETA(I))
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
      OPEN(11,FILE='JK-CIRCLE.txt',FORM='FORMATTED')
      OPEN(12,FILE='JK-WING.txt',FORM='FORMATTED')
      OPEN(13,FILE='JK-CAMBER.txt',FORM='FORMATTED')
      OPEN(14,FILE='JK-THICK.txt',FORM='FORMATTED')
      OPEN(15,FILE='JK-CHORD.txt',FORM='FORMATTED')
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
