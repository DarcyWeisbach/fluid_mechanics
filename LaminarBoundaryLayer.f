C
C********* ********* ********* ********* ********* ********* *********C
C
C     LaminarBoundaryLayer.f
C
C     Program for Laminar boundary layer
C     MADE BY MASAYOSHI OKAMOTO
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NY=500)
      DIMENSION Y(0:NY+1),YY(0:NY+1),U(0:NY+1),FU(0:NY+1),UO(0:NY+1)
      DIMENSION V(0:NY+1),VO(0:NY+1)
      DIMENSION UE(0:NY)
C
      UY(J)=0.5D0*(U(J+1)+U(J))
      FLO(J)=0.5D0*(FU(J)+FU(J-1))
      VY(J)=0.5D0*(V(J+1)+V(J))
C
      DY=0.02D0
      DO 100 J=0,NY
        YY(J)=DY*DFLOAT(J)
        UE(J)=DERF(0.5D0*YY(J))
 100  CONTINUE
      DO 110 J=1,NY
        Y(J)=(YY(J)+YY(J-1))*0.5D0
        U(J)=DERF(Y(J))
        V(J)=DERF(Y(J))
 110  CONTINUE
      U(0)=-U(1)
      U(NY+1)=1.D0
      V(0)=-V(1)
      V(NY+1)=1.D0
      DO 120 ITER=1,1000000
        FU(0)=0.D0
        DO 130 J=1,NY
          FU(J)=FU(J-1)+U(J)*DY
          UO(J)=U(J)
          VO(J)=V(J)
 130    CONTINUE
        U(1)=(4.D0+DY*FLO(1))/(12.D0-DY*FLO(1))*U(2)
        V(1)=(4.D0+DY*Y(1))/(12.D0-DY*Y(1))*V(2)
        DO 140 J=2,NY
          U(J)=(0.5D0+0.125D0*DY*FLO(J))*U(J+1)
     &        +(0.5D0-0.125D0*DY*FLO(J))*U(J-1)
          V(J)=(0.5D0+0.125D0*DY*Y(J))*V(J+1)
     &        +(0.5D0-0.125D0*DY*Y(J))*V(J-1)
 140    CONTINUE
        ERROR=0.D0
        DO 150 J=1,NY
          ERROR=DMAX1(ERROR,DABS((U(J)-UO(J))/U(J)))
 150    CONTINUE
        IF(ERROR.LT.1.D-10) GOTO 180
 120  CONTINUE
 180  WRITE(6,*) ITER
      U(0)=-U(1)
      U(NY)=1.D0
      V(0)=-V(1)
      V(NY)=1.D0
      OPEN(40,FILE='LaminarBoundary.txt',FORM='FORMATTED')
      WRITE(40,*) 'Y   U    V    U(Linear)'
      WRITE(40,'(1H ,1P,4E24.10)') (YY(J),UY(J),VY(J),UE(J),J=0,NY)
      CLOSE(40)
      END
