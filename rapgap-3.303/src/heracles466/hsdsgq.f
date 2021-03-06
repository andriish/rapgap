C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      FUNCTION HSDSGQ(Q2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       DIFFERENCE OF THE RENORMALIZED PHOTON SELF ENERGY DIVIDED BY Q2
C       OF THE PERTURBATIVE RESULT WITH EFFECTIVE QUARK MASSES AND
C       THE PARAMETRIZATION OF BURKHARD
C       14.05.91 HS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX*16 SGEFFQ,HSFONE
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSPARL/ LPAR(20),LPARIN(12)

      IF (LPAR(7).LT.2) THEN
C...NO CORRECTION
        HSDSGQ=0D0
      ELSE
        SGEFFQ=+ALP1PI*(
     2   +(+(Q2+2D0*MU2  )*HSFONE(Q2,MU  ,MU  )
     2     +(Q2+2D0*MC2  )*HSFONE(Q2,MC  ,MC  ) - 2D0*Q2/3D0)/2.25D0
     1   +(+(Q2+2D0*MD2  )*HSFONE(Q2,MD  ,MD  )
     1     +(Q2+2D0*MS2  )*HSFONE(Q2,MS  ,MS  )
     1     +(Q2+2D0*MB2  )*HSFONE(Q2,MB  ,MB  ) - Q2        )/9D0   )
        HSDSGQ=-HSHADQ(Q2)-DREAL(SGEFFQ)/Q2
      ENDIF
      RETURN
      END
