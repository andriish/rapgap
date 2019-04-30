      DOUBLE PRECISION FUNCTION BLUEML(x,g2,Q2)
C*******************************************************************C
C     Gluon distribution by the method of J.Bluemlein               C
C      BLUEML = x*gluon(x,kt2,Q2)               DESY 94-121         C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION(A-G,P-Z)
      EXTERNAL BJ0,BI0
      COMMON/CONVOL/cx,cg2,cQ2
      COMMON/CONST/PI,ALP,ALS
C... Parameters and constants
*! pi constant
      PI  = 3.1415926
*! electromagnetic constant
      ALP = 1./137.
*! strong coupling constant
      ALS = 0.2
      BLUEML=0.D0
      cx = x
      cg2= g2
      cQ2= dabs(Q2)
*!Precision of integration
      acc= 2.D-3
      IF(cg2.GE.cQ2) THEN
         BLUEML=DGAUSS(BI0, x, 1.D0, acc)
      ENDIF
      IF(cg2.LT.cQ2) THEN
         BLUEML=DGAUSS(BJ0, x, 1.D0, acc)
      ENDIF
      BLUEML=BLUEML/cg2
      RETURN
      END

      DOUBLE PRECISION FUNCTION BJ0(eta)
      IMPLICIT DOUBLE PRECISION(A-G,P-Z)
      COMMON/CONST/PI,ALP,ALS
      COMMON/CONVOL/cx,cg2,cQ2
*!Leading order GRV
      LO=0
      Alsb= ALS*3./PI
      argj= 2.*dsqrt(-Alsb*dlog(eta)*dabs(dlog(cg2/cQ2)))
      argg= cx/eta
      BJ0 = Alsb/eta * DBESJ0(argj) * GRVg(argg,cQ2,LO)
      RETURN
      END

      DOUBLE PRECISION FUNCTION BI0(eta)
      IMPLICIT DOUBLE PRECISION(A-G,P-Z)
      COMMON/CONST/PI,ALP,ALS
      COMMON/CONVOL/cx,cg2,cQ2
*!Leading order GRV
      LO=0
      Alsb= ALS*3./PI
      argj= 2.*dsqrt(-Alsb*dlog(eta)*dabs(dlog(cg2/cQ2)))
      argg= cx/eta
      BI0 = Alsb/eta * DBESI0(argj) * GRVg(argg,cQ2,LO)
      RETURN
      END
C###################################################################C
C                     Parton distributions                          C
C     M.Gluck & E.Reya & A.Vogt:  DO-TH 94/24;  DESY 94-206         C
C###################################################################C
C
C
      DOUBLE PRECISION FUNCTION GRVg(X,Q2,L)
C*******************************************************************C
C     Gluon distribution function:      GRVg = x*gluon(x,Q2)        C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION(A-G,P-Z)
C...  L=0 >Leading Order;  L=1 >NLO (MS bar); L=2 >NLO (DIS)
      QQ2=DABS(Q2)
      IF(QQ2.LT.0.4) QQ2=0.4D0
      T0=DLOG(.23D0/(0.232**2))
      T =DLOG(QQ2/(0.232**2))
      S =DLOG(T/T0)
      S2=S**2
      S3=S2*S
      SR=DSQRT(S)
      IF(L-1) 10,20,30
C... Leading Order set
   10 AL = 0.524
      BE = 1.088
      AA = 1.742 -0.930*S
      BB =                -0.399*S2
      A  = 7.486 -2.185*S
      B  = 16.69 -22.74*S +5.779*S2
      C  =-25.59 +29.71*S -7.296*S2
      D  = 2.792 +2.215*S +0.422*S2 -0.104*S3
      E  = 0.807 +2.005*S
      EP = 3.841 +0.316*S
      GOTO 40
C... Next-to-Leading Order set, MS bar scheme
   20 AL = 1.014
      BE = 1.738
      AA = 1.724 +0.157*S
      BB = 0.800 +1.016*S
      A  = 7.517 -2.547*S
      B  = 34.09 +17.47*S -52.21*DSQRT(S)
      C  = 4.039 +1.491*S
      D  = 3.404 +0.830*S
      E  =-1.112 +3.438*S -0.302*S2
      EP = 3.256 +0.436*S
      GOTO 40
C... Next-to-Leading order set, DIS scheme
   30 AL = 1.258
      BE = 1.846
      AA = 2.423
      BB = 2.427 +1.311*S -0.153*S2
      A  = 25.09 -7.935*S
      B  =-14.84 +72.18*S           -124.3*SR
      C  = 590.3 -173.8*S
      D  = 5.196 +1.857*S
      E  =-1.648 +3.988*S -0.432*S2
      EP = 3.232 -0.542*S
   40 GRVg =(X**AA *(A +B*X +C*X**2)*(-DLOG(X))**BB
     +     +S**AL *DEXP(-E +DSQRT(-EP*S**BE*DLOG(X))))*(1.-X)**D
      RETURN
      END
