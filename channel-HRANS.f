C
C********* ********* ********* ********* ********* ********* *********C
C
C     PROGRAM channel-HRANS.f
C
C     Program for turbulent channel flow 
C             by high-Reynolds-number K-E models
C     MADE BY MASAYOSHI OKAMOTO
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NY=100)
      COMMON
     & /VAL/UA(0:NY+1),TK(0:NY+1),TE(0:NY+1),OM(0:NY+1)
     & /IMP/TVIS(0:NY+1),UV(0:NY+1),SH(0:NY+1)
     & /OLD/UAO(0:NY+1),TKO(0:NY+1),TEO(0:NY+1),OMO(0:NY+1)
      COMMON
     & /LAU/APUA(NY),ANUA(NY),ASUA(NY),BUA(NY)
     & /LAK/APTK(NY),ANTK(NY),ASTK(NY),BTK(NY)
     & /LAE/APTE(NY),ANTE(NY),ASTE(NY),BTE(NY)
     & /LAE/APOM(NY),ANOM(NY),ASOM(NY),BOM(NY)
      COMMON
     & /YPOS/Y(0:NY+1),YW(0:NY+1)
      COMMON
     & /CON1/NITER,IMODEL
     & /CON2/VIS,PG,DY
     & /CON3/CNU,CE1,CE2,CO1,CO2,SGK,SGE,SGO
C
      CALL SETINIT(N)
      DO 100 IT=1,NITER
        CALL SETREY(N)
        CALL SOLVER(N)
        CALL CONVER(COV,N)
        IF(COV.LT.1.D-5) GOTO 1000
 100  CONTINUE
      GOTO 1100
 1000 WRITE(6,*) 'FINISH'
      WRITE(6,*) IT
      OPEN(12,FILE='DATA.txt',FORM='FORMATTED')
      OPEN(13,FILE='DATA-WALLUNIT.txt',FORM='FORMATTED')
      WRITE(12,*) 'Y    UA    UV    TK'
      WRITE(12,'( 1H ,1P,4E14.4 )')
     &     (Y(J),UA(J),UV(J),TK(J),J=1,N)
      WRITE(13,*) 'YW   UA+    UV+    TK+'
      WRITE(13,'( 1H ,1P,4E14.4 )')
     &     (YW(J),UA(J),UV(J),TK(J),J=1,N/2)
      CLOSE(12)
      CLOSE(13)
      GOTO 1200
 1100 WRITE(6,*) 'UNFINISH'
 1200 STOP
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE SETINIT(N)
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NY=100)
      COMMON
     & /VAL/UA(0:NY+1),TK(0:NY+1),TE(0:NY+1),OM(0:NY+1)
     & /IMP/TVIS(0:NY+1),UV(0:NY+1),SH(0:NY+1)
     & /OLD/UAO(0:NY+1),TKO(0:NY+1),TEO(0:NY+1),OMO(0:NY+1)
      COMMON
     & /LAU/APUA(NY),ANUA(NY),ASUA(NY),BUA(NY)
     & /LAK/APTK(NY),ANTK(NY),ASTK(NY),BTK(NY)
     & /LAE/APTE(NY),ANTE(NY),ASTE(NY),BTE(NY)
     & /LAE/APOM(NY),ANOM(NY),ASOM(NY),BOM(NY)
      COMMON
     & /YPOS/Y(0:NY+1),YW(0:NY+1)
      COMMON
     & /CON1/NITER,IMODEL
     & /CON2/VIS,PG,DY
     & /CON3/CNU,CE1,CE2,CO1,CO2,SGK,SGE,SGO
C
C ***** ZERO CLEAR *****
C
      DO 100 J=0,NY+1
        UA(J)=0.D0
        TK(J)=0.D0
        TE(J)=0.D0
        OM(J)=0.D0
        UV(J)=0.D0
        TVIS(J)=0.D0
        SH(J)=0.D0
        UAO(J)=0.D0
        TKO(J)=0.D0
        TEO(J)=0.D0
        OMO(J)=0.D0
        Y(J)=0.D0
        YW(J)=0.D0
 100  CONTINUE
      DO 110 J=1,NY
        APUA(J)=0.D0
        ANUA(J)=0.D0
        ASUA(J)=0.D0
        BUA(J)=0.D0
        APTK(J)=0.D0
        ANTK(J)=0.D0
        ASTK(J)=0.D0
        BTK(J)=0.D0
        APTE(J)=0.D0
        ANTE(J)=0.D0
        ASTE(J)=0.D0
        BOM(J)=0.D0
        APOM(J)=0.D0
        ANOM(J)=0.D0
        ASOM(J)=0.D0
        BOM(J)=0.D0
 110  CONTINUE
C
C     CONSTANTS
C
      NITER=500000
C
C     MODEL CONSTANTS
C
      WRITE(6,*) '0: K-Epsilon model'
      WRITE(6,*) '1: K-Omega model'
      READ(5,*) IMODEL
      WRITE(6,*) 'N: Number of grids'
      READ(5,*) N
C
      IF(IMODEL.EQ.0) THEN
        CNU=0.09D0
        CE1=1.44D0
        CE2=1.92D0
        SGK=1.00D0
        SGE=1.30D0
      ELSEIF(IMODEL.EQ.1) THEN
        CNU=0.09D0
        CO1=5.0D0/9.D0
        CO2=3.0D0/40.D0
        SGK=0.50D0
        SGO=0.50D0
      ENDIF
C
C     FLOW STATE
C
      PG=-1.D0
      WRITE(6,*) 'REYNOLDS NUMBER'
      READ(5,*) RE
      VIS=1.D0/RE
C
C     SET GRID
C
      CALL SETGRID(N)
C
C     INITIAL FIELD
C
      DO 130 J=1,N-1
        UA(J)=2.5D0*DLOG(YW(J))+5.D0
        TK(J)=1.D0
        TE(J)=10.D0
        OM(J)=TE(J)/TK(J)/0.09D0
 130  CONTINUE
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE SETGRID(N)
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NY=100)
      COMMON
     & /YPOS/Y(0:NY+1),YW(0:NY+1)
      COMMON
     & /CON1/NITER,IMODEL
     & /CON2/VIS,PG,DY
     & /CON3/CNU,CE1,CE2,CO1,CO2,SGK,SGE,SGO
C
      YMIN=-1.D0+100.D0*VIS
      YMAX=1.D0-100.D0*VIS
      YH=YMAX-YMIN
      DY=YH/DFLOAT(N-2)
      Y(0)=-1.D0
      Y(1)=YMIN
      YW(1)=(1.D0+Y(1))/VIS
      DO 100 J=2,N-2
        Y(J)=Y(J-1)+DY
        YW(J)=DMIN1(1.D0+Y(J),1.D0-Y(J))/VIS
 100  CONTINUE
      Y(N-1)=YMAX
      YW(N-1)=(1.D0-Y(N-1))/VIS
      Y(N)=1.D0
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE SETREY(N)
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NY=100)
      COMMON
     & /VAL/UA(0:NY+1),TK(0:NY+1),TE(0:NY+1),OM(0:NY+1)
     & /IMP/TVIS(0:NY+1),UV(0:NY+1),SH(0:NY+1)
     & /OLD/UAO(0:NY+1),TKO(0:NY+1),TEO(0:NY+1),OMO(0:NY+1)
      COMMON
     & /YPOS/Y(0:NY+1),YW(0:NY+1)
      COMMON
     & /CON1/NITER,IMODEL
     & /CON2/VIS,PG,DY
     & /CON3/CNU,CE1,CE2,CO1,CO2,SGK,SGE,SGO
C
C********** STATEMENT FUNCTION **********
C
      UAC(J)=0.5D0*(UA(J)+UA(J-1))
C
      DO 100 J=2,N-2
        SH(J)=(UAC(J+1)-UAC(J))/DY
        IF(IMODEL.EQ.0) THEN
          TVIS(J)=CNU*TK(J)**2/TE(J)
        ELSEIF(IMODEL.EQ.1) THEN
          TVIS(J)=TK(J)/OM(J)
        ENDIF
        UV(J)=TVIS(J)*SH(J)
 100  CONTINUE
      IF(IMODEL.EQ.0) THEN
        TVIS(1)=CNU*TK(1)**2/TE(1)
        UV(1)=2.D0*TVIS(1)*(UAC(2)-UA(1))/DY
        TVIS(N-1)=CNU*TK(N-1)**2/TE(N-1)
        UV(N-1)=2.D0*TVIS(N-1)*(UA(N-1)-UAC(N-1))/DY
      ELSEIF(IMODEL.EQ.1) THEN
        TVIS(1)=TK(1)/OM(1)
        UV(1)=2.D0*TVIS(1)*(UAC(2)-UA(1))/DY
        TVIS(N-1)=TK(N-1)/OM(N-1)
        UV(N-1)=2.D0*TVIS(N-1)*(UA(N-1)-UAC(N-1))/DY
      ENDIF
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE BOUND(N)
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NY=100)
      COMMON
     & /VAL/UA(0:NY+1),TK(0:NY+1),TE(0:NY+1),OM(0:NY+1)
     & /IMP/TVIS(0:NY+1),UV(0:NY+1),SH(0:NY+1)
     & /OLD/UAO(0:NY+1),TKO(0:NY+1),TEO(0:NY+1),OMO(0:NY+1)
      COMMON
     & /YPOS/Y(0:NY+1),YW(0:NY+1)
      COMMON
     & /CON1/NITER,IMODEL
     & /CON2/VIS,PG,DY
     & /CON3/CNU,CE1,CE2,CO1,CO2,SGK,SGE,SGO
C
      ERR=0.D0
      UTB=1.D0
      UTT=1.D0
      DO 100 IT=1,100
        UTBO=UTB
        UTTO=UTT
        UTB=DSQRT(0.42D0*UA(1)*UTB
     &     /DLOG(9.D0*UTB*(1.D0+Y(1))/VIS))
        UTT=DSQRT(0.42D0*UA(N-1)*UTT
     &     /DLOG(9.D0*UTT*(1.D0-Y(N-1))/VIS))
        ERRB=DABS(UTB-UTBO)/UTB
        ERR=DMAX1(ERR,ERRB)
        ERRT=DABS(UTT-UTTO)/UTT
        ERR=DMAX1(ERR,ERRT)
        IF(ERR.LT.1.D-5) GOTO 110
 100  CONTINUE
 110  UA(1)=UTB/0.42D0*DLOG(9.D0*UTB*(1.D0+Y(1))/VIS)
      UA(N-1)=UTT/0.42D0*DLOG(9.D0*UTB*(1.D0-Y(N-1))/VIS)
      TK(1)=UTB**2/0.3D0
      TK(N-1)=UTT**2/0.3D0
      IF(IMODEL.EQ.0) THEN
        TE(1)=UTB**3/0.42D0/(1.D0+Y(1))
        TE(N-1)=UTT**3/0.42D0/(1.D0-Y(N-1))
      ELSEIF(IMODEL.EQ.1) THEN
        OM(1)=UTB/0.126D0/(1.D0+Y(1))
        OM(N-1)=UTT/0.126D0/(1.D0-Y(N-1))
      ENDIF
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE SOLVER(N)
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NY=100)
      COMMON
     & /VAL/UA(0:NY+1),TK(0:NY+1),TE(0:NY+1),OM(0:NY+1)
     & /IMP/TVIS(0:NY+1),UV(0:NY+1),SH(0:NY+1)
     & /OLD/UAO(0:NY+1),TKO(0:NY+1),TEO(0:NY+1),OMO(0:NY+1)
      COMMON
     & /LAU/APUA(NY),ANUA(NY),ASUA(NY),BUA(NY)
     & /LAK/APTK(NY),ANTK(NY),ASTK(NY),BTK(NY)
     & /LAE/APTE(NY),ANTE(NY),ASTE(NY),BTE(NY)
     & /LAE/APOM(NY),ANOM(NY),ASOM(NY),BOM(NY)
      COMMON
     & /YPOS/Y(0:NY+1),YW(0:NY+1)
      COMMON
     & /CON1/NITER,IMODEL
     & /CON2/VIS,PG,DY
     & /CON3/CNU,CE1,CE2,CO1,CO2,SGK,SGE,SGO
C
C********** STATEMENT FUNCTION **********
C
      ANGU(J)=VIS+TVIS(J)
      ANGK(J)=VIS+TVIS(J)/SGK
      ANGE(J)=VIS+TVIS(J)/SGE
      ANGO(J)=VIS+TVIS(J)/SGO
      ANGUO(J)=0.5D0*(ANGU(J)+ANGU(J-1))
      ANGKO(J)=0.5D0*(ANGK(J)+ANGK(J-1))
      ANGEO(J)=0.5D0*(ANGE(J)+ANGE(J-1))
      ANGOO(J)=0.5D0*(ANGO(J)+ANGO(J-1))
C
      DO 100 J=1,N
        UAO(J)=UA(J)
        TKO(J)=TK(J)
        TEO(J)=TE(J)
        OMO(J)=OM(J)
 100  CONTINUE
      DO 110 J=2,N-2
        ANUA(J)=ANGUO(J+1)/DY
        ASUA(J)=ANGUO(J)/DY
        APUA(J)=(ANGUO(J+1)+ANGUO(J))/DY
        BUA(J)=-DY*PG
        UA(J)=(BUA(J)+ANUA(J)*UA(J+1)+ASUA(J)*UA(J-1))/APUA(J)
        IF(IMODEL.EQ.0) THEN
          ANTK(J)=ANGKO(J+1)/DY
          ASTK(J)=ANGKO(J)/DY
          APTK(J)=(ANGKO(J+1)+ANGKO(J))/DY+DY*TE(J)/TK(J)
          BTK(J)=DY*UV(J)*SH(J)
          TK(J)=(BTK(J)+ANTK(J)*TK(J+1)+ASTK(J)*TK(J-1))/APTK(J)
          ANTE(J)=ANGEO(J+1)/DY
          ASTE(J)=ANGEO(J)/DY
          APTE(J)=(ANGEO(J+1)+ANGEO(J))/DY
     &           +DY*CE2*TE(J)/TK(J)
          BTE(J)=DY*CE1*TE(J)*UV(J)*SH(J)/TK(J)
          TE(J)=(BTE(J)+ANTE(J)*TE(J+1)+ASTE(J)*TE(J-1))/APTE(J)
        ELSEIF(IMODEL.EQ.1) THEN
          ANTK(J)=ANGKO(J+1)/DY
          ASTK(J)=ANGKO(J)/DY
          APTK(J)=(ANGKO(J+1)+ANGKO(J))/DY+DY*CNU*OM(J)
          BTK(J)=DY*UV(J)*SH(J)
          TK(J)=(BTK(J)+ANTK(J)*TK(J+1)+ASTK(J)*TK(J-1))/APTK(J)
          ANOM(J)=ANGOO(J+1)/DY
          ASOM(J)=ANGOO(J)/DY
          APOM(J)=(ANGOO(J+1)+ANGOO(J))/DY+DY*CO2*OM(J)
          BOM(J)=DY*CO1*OM(J)*UV(J)*SH(J)/TK(J)
          OM(J)=(BOM(J)+ANOM(J)*OM(J+1)+ASOM(J)*OM(J-1))/APOM(J)
        ENDIF
 110  CONTINUE
      CALL BOUND(N)
      RETURN
      END
C
C********* ********* ********* ********* ********* ********* *********C
C
      SUBROUTINE CONVER(COV,N)
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NY=100)
      COMMON
     & /VAL/UA(0:NY+1),TK(0:NY+1),TE(0:NY+1),OM(0:NY+1)
     & /IMP/TVIS(0:NY+1),UV(0:NY+1),SH(0:NY+1)
     & /OLD/UAO(0:NY+1),TKO(0:NY+1),TEO(0:NY+1),OMO(0:NY+1)
C
      TESTUA=0.D0
      TESTTK=0.D0
      TESTTE=0.D0
      TESTOM=0.D0
      COV=1.D-14
      DO 100 J=1,N-1
        TESTUA=DABS(UA(J)-UAO(J))/UAO(J)
        TESTTK=DABS(TK(J)-TKO(J))/TKO(J)
        COV=DMAX1(TESTUA,COV)
        COV=DMAX1(TESTTK,COV)
        IF(IMODEL.EQ.0) THEN
          TESTTE=DABS(TE(J)-TEO(J))/TEO(J)
          COV=DMAX1(TESTTE,COV)
        ELSEIF(IMODEL.EQ.1) THEN
          TESTOM=DABS(OM(J)-OMO(J))/OMO(J)
          COV=DMAX1(TESTOM,COV)
        ENDIF
 100  CONTINUE
      RETURN
      END
