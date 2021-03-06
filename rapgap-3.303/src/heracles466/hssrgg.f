CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  RENORMIERTE SELBSTENERGIEN
C  IN DEN UNTERPROGRAMMEN
C     HSSRGG, HSSRGZ, HSSRZZ, HSSRWW
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      FUNCTION HSSRGG(Q2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       HSSRGG(Q2): BEITRAEGE ZUR PHOTONSELBSTENERGIE  (RENORMIERT)
C              Q2 = Q**2 VIERERIMPULSQUADRAT
C-----------------------------------------------------------------------
C       09.11.83
C       14.05.91 HS: BURKHARD'S PARAMETRIZATION
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX*16 HSSRGG,HSFONE
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      IF (LPAR(7).EQ.1) THEN
C...HADRONIC CONTRIBUTION FROM EFFECTIVE QUARK LOOPS
        HSSRGG=+ALP1PI*(
     G   -DFLOAT(LPAR(15))*((3D0*Q2+4D0*MW2)*HSFONE(Q2,MW,MW)
     G                                       - Q2/1.5D0)/4D0
     L   +(+(Q2+2D0*ME2  )*HSFONE(Q2,ME  ,ME  )
     L     +(Q2+2D0*MMY2 )*HSFONE(Q2,MMY ,MMY )
     L     +(Q2+2D0*MTAU2)*HSFONE(Q2,MTAU,MTAU) - Q2)/3D0
     2   +(+(Q2+2D0*MU2  )*HSFONE(Q2,MU  ,MU  )
     2     +(Q2+2D0*MC2  )*HSFONE(Q2,MC  ,MC  )
     2     +(Q2+2D0*MT2  )*HSFONE(Q2,MT  ,MT  ) - Q2)/2.25D0
     1   +(+(Q2+2D0*MD2  )*HSFONE(Q2,MD  ,MD  )
     1     +(Q2+2D0*MS2  )*HSFONE(Q2,MS  ,MS  )
     1     +(Q2+2D0*MB2  )*HSFONE(Q2,MB  ,MB  ) - Q2)/9D0        )
      ELSEIF(LPAR(7).EQ.2) THEN
C...HADRONIC CONTRIBUTION FROM BURKHARDT'S PARAMETRIZATION
        HSSRGG=+ALP1PI*(
     G   -DFLOAT(LPAR(15))*((3D0*Q2+4D0*MW2)*HSFONE(Q2,MW,MW)
     G                                       - Q2/1.5D0)/4D0
     L   +(+(Q2+2D0*ME2  )*HSFONE(Q2,ME  ,ME  )
     L     +(Q2+2D0*MMY2 )*HSFONE(Q2,MMY ,MMY )
     L     +(Q2+2D0*MTAU2)*HSFONE(Q2,MTAU,MTAU) - Q2)/3D0
     2   +(+(Q2+2D0*MT2  )*HSFONE(Q2,MT  ,MT  ) - Q2/3D0)/2.25D0     )
     H   - DCMPLX(HSHADQ(Q2)*Q2,0D0)
      ELSEIF(LPAR(7).EQ.0) THEN
       HSSRGG=DCMPLX(0D0,0D0)
      ELSE
      ENDIF
      RETURN
      END
