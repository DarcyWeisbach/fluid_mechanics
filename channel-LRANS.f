C
C********* ********* ********* ********* ********* ********* *********C
C
C     PROGRAM channel-LRANS.f
C
C     Program for turbulent channel flow 
C             by Low-Reynolds-number K-E models
C     MADE BY MASAYOSHI OKAMOTO
C
C********* ********* ********* ********* ********* ********* *********C
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NY=1000)
      COMMON
     & /VAL/UA(0:NY+1),TK(0:NY+1),TE(0:NY+1),EP(0:NY+1)
     & /IMP/TVIS(0:NY+1),UV(0:NY+1),SH(0:NY+1)
     & /FNC/DTM(0:NY+1),ETM(0:NY+1),FMU(0:NY+1),FE(0:NY+1)
     & /OLD/UAO(0:NY+1),TKO(0:NY+1),TEO(0:NY+1),EPO(0:NY+1)
      COMMON
     & /LAU/APUA(NY),ANUA(NY),ASUA(NY),BUA(NY)
     & /LAK/APTK(NY),ANTK(NY),ASTK(NY),BTK(NY)
     & /LAE/APTE(NY),ANTE(NY),ASTE(NY),BTE(NY)
      COMMON
     & /YPOS/Y(0:NY+1),YY(0:NY+1),YW(0:NY+1)
     & /YGRD/DY(0:NY+1),DYY(0:NY+1)
     & /YINP/WVN(0:NY+1),WVS(0:NY+1)
      COMMON
     & /CON1/NITER,IMODEL
     & /CON2/VIS,PG
     & /CON3/CNU,CE1,CE2,SGK,SGE
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
      IF(IMODEL.NE.3) THEN
        WRITE(12,*) 'Y    UA    UV    TK    EP'
        WRITE(12,'( 1H ,1P,5E14.4 )')
     &       (Y(J),UA(J),UV(J),TK(J),TE(J)+DTM(J),J=1,N)
        WRITE(13,*) 'YW   UA+    UV+    TK+    EP+'
        WRITE(13,'( 1H ,1P,5E14.4 )')
     &       (YW(J),UA(J),UV(J),TK(J),(TE(J)+DTM(J))*VIS,J=1,N/2)
      ELSE
        WRITE(12,*) 'Y    UA    UV    TK    EP'
        WRITE(12,'( 1H ,1P,5E14.4 )')
     &       (Y(J),UA(J),UV(J),TK(J),EP(J),J=1,N)
        WRITE(13,*) 'YW   UA+    UV+    TK+    EP+'
        WRITE(13,'( 1H ,1P,5E14.4 )')
     &       (YW(J),UA(J),UV(J),TK(J),EP(J)*VIS,J=1,N/2)
      ENDIF
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
      PARAMETER(NY=1000)
      COMMON
     & /VAL/UA(0:NY+1),TK(0:NY+1),TE(0:NY+1),EP(0:NY+1)
     & /IMP/TVIS(0:NY+1),UV(0:NY+1),SH(0:NY+1)
     & /FNC/DTM(0:NY+1),ETM(0:NY+1),FMU(0:NY+1),FE(0:NY+1)
     & /OLD/UAO(0:NY+1),TKO(0:NY+1),TEO(0:NY+1),EPO(0:NY+1)
      COMMON
     & /LAU/APUA(NY),ANUA(NY),ASUA(NY),BUA(NY)
     & /LAK/APTK(NY),ANTK(NY),ASTK(NY),BTK(NY)
     & /LAE/APTE(NY),ANTE(NY),ASTE(NY),BTE(NY)
      COMMON
     & /YPOS/Y(0:NY+1),YY(0:NY+1),YW(0:NY+1)
     & /YGRD/DY(0:NY+1),DYY(0:NY+1)
     & /YINP/WVN(0:NY+1),WVS(0:NY+1)
      COMMON
     & /CON1/NITER,IMODEL
     & /CON2/VIS,PG
     & /CON3/CNU,CE1,CE2,SGK,SGE
C
C ***** ZERO CLEAR *****
C
      DO 100 J=0,NY+1
        UA(J)=0.D0
        TK(J)=0.D0
        TE(J)=0.D0
        EP(J)=0.D0
        UV(J)=0.D0
        TVIS(J)=0.D0
        SH(J)=0.D0
        UAO(J)=0.D0
        TKO(J)=0.D0
        TEO(J)=0.D0
        EPO(J)=0.D0
        DTM(J)=0.D0
        ETM(J)=0.D0
        FMU(J)=0.D0
        FE(J)=0.D0
        Y(J)=0.D0
        YY(J)=0.D0
        YW(J)=0.D0
        DY(J)=0.D0
        DYY(J)=0.D0
        WVN(J)=0.D0
        WVS(J)=0.D0
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
        BTE(J)=0.D0
 110  CONTINUE
C
C     CONSTANTS
C
      NITER=500000
C
C     MODEL CONSTANTS
C
      WRITE(6,*) '0: Launder-Sharma model'
      WRITE(6,*) '1: Cotton-Kirwin model'
      WRITE(6,*) '2: Nagano-Hishida model'
      WRITE(6,*) '3: Abe-Kondoh-Nagano model'
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
        CE1=1.44D0
        CE2=1.92D0
        SGK=1.00D0
        SGE=1.30D0
      ELSEIF(IMODEL.EQ.2) THEN
        CNU=0.09D0
        CE1=1.45D0
        CE2=1.90D0
        SGK=1.00D0
        SGE=1.30D0
      ELSEIF(IMODEL.EQ.3) THEN
        CNU=0.09D0
        CE1=1.50D0
        CE2=1.90D0
        SGK=1.40D0
        SGE=1.40D0
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
      DO 130 J=1,N
        UA(J)=DMIN1(YW(J),2.5D0*DLOG(YW(J))+5.D0)
        TK(J)=3.D0*DMIN1(YW(J)/20.D0,1.D0)
        TE(J)=10.D0*DMIN1(YW(J)/20.D0,1.D0)
        EP(J)=10.D0*DMIN1(YW(J)/20.D0,1.D0)
 130  CONTINUE
C
C     BOUNDARY CONDITION
C
      UA(0)=-UA(1)
      UA(N+1)=-UA(N)
      TK(0)=-TK(1)
      TK(N+1)=-TK(N)
      TE(0)=-TE(1)
      TE(N+1)=-TE(N)
      EP(0)=2.D0*(8.D0*VIS*TK(1)/DY(1)**2)-EP(1)
      EP(N+1)=2.D0*(8.D0*VIS*TK(N)/DY(N)**2)-EP(N)
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
      PARAMETER(NY=1000)
      COMMON
     & /YPOS/Y(0:NY+1),YY(0:NY+1),YW(0:NY+1)
     & /YGRD/DY(0:NY+1),DYY(0:NY+1)
     & /YINP/WVN(0:NY+1),WVS(0:NY+1)
      COMMON
     & /CON1/NITER,IMODEL
     & /CON2/VIS,PG
     & /CON3/CNU,CE1,CE2,SGK,SGE
C
      CCC=0.9D0
      CTC=0.5D0*DLOG((1.D0+CCC)/(1.D0-CCC))
      YY(0)=-1.D0
      DO 100 J=1,N-1
        YY(J)=DTANH(CTC*(-1.D0+2.D0*DFLOAT(J)/DFLOAT(N)))/CCC
 100  CONTINUE
      YY(N)=1.D0
      DO 110 J=1,N
        DY(J)=YY(J)-YY(J-1)
 110  CONTINUE
      DY(0)=DY(1)
      DY(N+1)=DY(N)
      DO 120 J=0,N
        Y(J)=YY(J)-0.5D0*DY(J)
        DYY(J)=0.5D0*DY(J)+0.5D0*DY(J+1)
        YW(J)=DMIN1(1.D0+Y(J),1.D0-Y(J))/VIS
 120  CONTINUE
      DO 130 J=0,N
        WVN(J)=0.5D0*DY(J)/DYY(J)
        WVS(J)=0.5D0*DY(J+1)/DYY(J)
 130  CONTINUE
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
      PARAMETER(NY=1000)
      COMMON
     & /VAL/UA(0:NY+1),TK(0:NY+1),TE(0:NY+1),EP(0:NY+1)
     & /IMP/TVIS(0:NY+1),UV(0:NY+1),SH(0:NY+1)
     & /FNC/DTM(0:NY+1),ETM(0:NY+1),FMU(0:NY+1),FE(0:NY+1)
     & /OLD/UAO(0:NY+1),TKO(0:NY+1),TEO(0:NY+1),EPO(0:NY+1)
      COMMON
     & /LAU/APUA(NY),ANUA(NY),ASUA(NY),BUA(NY)
     & /LAK/APTK(NY),ANTK(NY),ASTK(NY),BTK(NY)
     & /LAE/APTE(NY),ANTE(NY),ASTE(NY),BTE(NY)
      COMMON
     & /YPOS/Y(0:NY+1),YY(0:NY+1),YW(0:NY+1)
     & /YGRD/DY(0:NY+1),DYY(0:NY+1)
     & /YINP/WVN(0:NY+1),WVS(0:NY+1)
      COMMON
     & /CON1/NITER,IMODEL
     & /CON2/VIS,PG
     & /CON3/CNU,CE1,CE2,SGK,SGE
      DIMENSION YK(NY),RT(NY),RTN(NY)
C
C********** STATEMENT FUNCTION **********
C
      STKY(J)=DSQRT(DABS(WVN(J)*TK(J+1)+WVS(J)*TK(J)))
      UAY(J)=WVN(J)*UA(J+1)+WVS(J)*UA(J)
C
      DO 100 J=1,N
        SH(J)=(UAY(J)-UAY(J-1))/DY(J)
 100  CONTINUE
      DO 110 J=1,N
        RT(J)=TK(J)**2/VIS/TE(J)
        RTN(J)=TK(J)**2/VIS/EP(J)
        YK(J)=(1.D0-DABS(Y(J)))*EP(J)**0.25D0/VIS**0.75D0
        IF(IMODEL.EQ.0) THEN
          FMU(J)=DEXP(-3.4D0/(1.D0+RT(J)/50.D0)**2)
          FE(J)=1.D0-0.3D0*DEXP(-RT(J)**2)
          TVIS(J)=CNU*TK(J)**2*FMU(J)/TE(J)
          DTM(J)=2.D0*VIS*(STKY(J)-STKY(J-1))**2/DY(J)**2
          ETM(J)=2.D0*VIS*TVIS(J)*(((UA(J+1)-UA(J))/DYY(J)
     &                 -(UA(J)-UA(J-1))/DYY(J-1))/DY(J))**2
        ELSEIF(IMODEL.EQ.1) THEN
          FMU(J)=1.D0-0.97D0*DEXP(-RT(J)/160.D0)
     &          -0.0045D0*RT(J)*DEXP(-(RT(J)/200.D0)**2)
          FE(J)=1.D0-0.3D0*DEXP(-RT(J)**2)
          TVIS(J)=CNU*TK(J)**2*FMU(J)/TE(J)
          DTM(J)=2.D0*VIS*(STKY(J)-STKY(J-1))**2/DY(J)**2
          ETM(J)=0.95D0*VIS*TVIS(J)*(((UA(J+1)-UA(J))/DYY(J)
     &                  -(UA(J)-UA(J-1))/DYY(J-1))/DY(J))**2
        ELSEIF(IMODEL.EQ.2) THEN
          FMU(J)=(1.D0-DEXP(-YW(J)/26.5D0))**2
          FE(J)=1.D0-0.3D0*DEXP(-RT(J)**2)
          TVIS(J)=CNU*TK(J)**2*FMU(J)/TE(J)
          DTM(J)=2.D0*VIS*(STKY(J)-STKY(J-1))**2/DY(J)**2
          ETM(J)=0.95D0*VIS*TVIS(J)*(((UA(J+1)-UA(J))/DYY(J)
     &                  -(UA(J)-UA(J-1))/DYY(J-1))/DY(J))**2
        ELSEIF(IMODEL.EQ.3) THEN
          FMU(J)=(1.D0-DEXP(-YK(J)/14.D0))**2
     &          *(1+5.D0/RTN(J)**0.75D0*DEXP(-(RTN(J)/200.D0)**2))
          FE(J)=(1.D0-DEXP(-YK(J)/3.1D0))**2
     &         *(1.D0-0.3D0*DEXP(-(RTN(J)/6.5D0)**2))
          TVIS(J)=CNU*TK(J)**2*FMU(J)/EP(J)
          DTM(J)=0.D0
          ETM(J)=0.D0
        ENDIF
        UV(J)=TVIS(J)*SH(J)
 110  CONTINUE
      UV(0)=-UV(1)
      UV(N+1)=-UV(N)
      TVIS(0)=-TVIS(1)
      TVIS(N+1)=-TVIS(N)
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
      PARAMETER(NY=1000)
      COMMON
     & /VAL/UA(0:NY+1),TK(0:NY+1),TE(0:NY+1),EP(0:NY+1)
     & /IMP/TVIS(0:NY+1),UV(0:NY+1),SH(0:NY+1)
     & /FNC/DTM(0:NY+1),ETM(0:NY+1),FMU(0:NY+1),FE(0:NY+1)
     & /OLD/UAO(0:NY+1),TKO(0:NY+1),TEO(0:NY+1),EPO(0:NY+1)
      COMMON
     & /LAU/APUA(NY),ANUA(NY),ASUA(NY),BUA(NY)
     & /LAK/APTK(NY),ANTK(NY),ASTK(NY),BTK(NY)
     & /LAE/APTE(NY),ANTE(NY),ASTE(NY),BTE(NY)
      COMMON
     & /YPOS/Y(0:NY+1),YY(0:NY+1),YW(0:NY+1)
     & /YGRD/DY(0:NY+1),DYY(0:NY+1)
     & /YINP/WVN(0:NY+1),WVS(0:NY+1)
      COMMON
     & /CON1/NITER,IMODEL
     & /CON2/VIS,PG
     & /CON3/CNU,CE1,CE2,SGK,SGE
C
C********** STATEMENT FUNCTION **********
C
      ANGU(J)=VIS+TVIS(J)
      ANGK(J)=VIS+TVIS(J)/SGK
      ANGE(J)=VIS+TVIS(J)/SGE
      ANGUY(J)=WVN(J)*ANGU(J+1)+WVS(J)*ANGU(J)
      ANGKY(J)=WVN(J)*ANGK(J+1)+WVS(J)*ANGK(J)
      ANGEY(J)=WVN(J)*ANGE(J+1)+WVS(J)*ANGE(J)
C
      DO 100 J=1,N
        UAO(J)=UA(J)
        TKO(J)=TK(J)
        TEO(J)=TE(J)
        EPO(J)=EP(J)
 100  CONTINUE
      DO 110 J=1,N
        ANUA(J)=ANGUY(J)/DYY(J)
        ASUA(J)=ANGUY(J-1)/DYY(J-1)
        APUA(J)=ANGUY(J)/DYY(J)+ANGUY(J-1)/DYY(J-1)
        BUA(J)=-DY(J)*PG
        UA(J)=(BUA(J)+ANUA(J)*UA(J+1)+ASUA(J)*UA(J-1))/APUA(J)
        IF(IMODEL.NE.3) THEN
          ANTK(J)=ANGKY(J)/DYY(J)
          ASTK(J)=ANGKY(J-1)/DYY(J-1)
          APTK(J)=ANGKY(J)/DYY(J)+ANGKY(J-1)/DYY(J-1)
     &           +DY(J)*(TE(J)+DTM(J))/TK(J)
          BTK(J)=DY(J)*UV(J)*SH(J)
          TK(J)=(BTK(J)+ANTK(J)*TK(J+1)+ASTK(J)*TK(J-1))/APTK(J)
          ANTE(J)=ANGEY(J)/DYY(J)
          ASTE(J)=ANGEY(J-1)/DYY(J-1)
          APTE(J)=ANGEY(J)/DYY(J)+ANGEY(J-1)/DYY(J-1)
     &           +DY(J)*CE2*TE(J)*FE(J)/TK(J)
          BTE(J)=DY(J)*CE1*TE(J)*UV(J)*SH(J)/TK(J)+DY(J)*ETM(J)
          TE(J)=(BTE(J)+ANTE(J)*TE(J+1)+ASTE(J)*TE(J-1))/APTE(J)
        ELSE
          ANTK(J)=ANGKY(J)/DYY(J)
          ASTK(J)=ANGKY(J-1)/DYY(J-1)
          APTK(J)=ANGKY(J)/DYY(J)+ANGKY(J-1)/DYY(J-1)
     &           +DY(J)*EP(J)/TK(J)
          BTK(J)=DY(J)*UV(J)*SH(J)
          TK(J)=(BTK(J)+ANTK(J)*TK(J+1)+ASTK(J)*TK(J-1))/APTK(J)
          ANTE(J)=ANGEY(J)/DYY(J)
          ASTE(J)=ANGEY(J-1)/DYY(J-1)
          APTE(J)=ANGEY(J)/DYY(J)+ANGEY(J-1)/DYY(J-1)
     &           +DY(J)*CE2*EP(J)*FE(J)/TK(J)
          BTE(J)=DY(J)*CE1*EP(J)*UV(J)*SH(J)/TK(J)
          EP(J)=(BTE(J)+ANTE(J)*EP(J+1)+ASTE(J)*EP(J-1))/APTE(J)
        ENDIF
 110  CONTINUE
      UA(0)=-UA(1)
      UA(N+1)=-UA(N)
      TK(0)=-TK(1)
      TK(N+1)=-TK(N)
      TE(0)=-TE(1)
      TE(N+1)=-TE(N)
      EP(0)=2.D0*(8.D0*VIS*TK(1)/DY(1)**2)-EP(1)
      EP(N+1)=2.D0*(8.D0*VIS*TK(N)/DY(N)**2)-EP(N)
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
      PARAMETER(NY=1000)
      COMMON
     & /VAL/UA(0:NY+1),TK(0:NY+1),TE(0:NY+1),EP(0:NY+1)
     & /IMP/TVIS(0:NY+1),UV(0:NY+1),SH(0:NY+1)
     & /FNC/DTM(0:NY+1),ETM(0:NY+1),FMU(0:NY+1),FE(0:NY+1)
     & /OLD/UAO(0:NY+1),TKO(0:NY+1),TEO(0:NY+1),EPO(0:NY+1)
C
      TESTUA=0.D0
      TESTTK=0.D0
      TESTTE=0.D0
      COV=1.D-14
      DO 100 J=1,N
        TESTUA=DABS(UA(J)-UAO(J))/UAO(J)
        TESTTK=DABS(TK(J)-TKO(J))/TKO(J)
        IF(IMODEL.NE.3) THEN
          TESTTE=DABS(TE(J)-TEO(J))/TEO(J)
        ELSE
          TESTTE=DABS(EP(J)-EPO(J))/TEO(J)
        ENDIF
        COV=DMAX1(TESTUA,COV)
        COV=DMAX1(TESTTK,COV)
        COV=DMAX1(TESTTE,COV)
 100  CONTINUE
      RETURN
      END
