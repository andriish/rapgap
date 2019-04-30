      FUNCTION PQPM(YY1,Q2)
      IMPLICIT None
#include "rgpart.inc"
#include "rgpara.inc"
#include "rgluco.inc"
#include "rglujets.inc"
#include "rgrapgki.inc"
#include "rgdisdif.inc"
      REAL XPQ(-6:6)
      REAL SNGL,XS,Q2S
      REAL SCALE
      Integer KPART
      COMMON/PINT/ SCALE,KPART
      Integer ID
      double precision xpom
      COMMON /PQPMID/xpom,ID
      Double Precision YY1,Q2,GEV2NB,X,F2,ALPH_EM,SIGMA,WTGLU,GF,WMAT
      Double Precision PQPM,beta
      Integer I
      DATA GEV2NB/.3893857D+6/

      X = Q2/YY1/SSS
      XS = REAL(X)
      Q2S = REAL(Q2)
c      CALL rgstfu(K(2,2),SNGL(X),SNGL(Q2),XPQ)
      do i=-6,6
         XPQ(i)=0.0
      enddo
      IF(ID.EQ.1) THEN
         CALL rgstfu(2212,SNGL(X),SNGL(Q2),XPQ)
      ELSE
         BETA = X/XPOM
c     write(6,*) xpom,x,beta,id
         IF(BETA.LT.1.and.beta.gt.0.) THEN
c            write(6,*) xpom,x,beta,id
            CALL RASTFU(KINT(2,2),SNGL(beta),SNGL(Q2),XPQ)
c        write(6,*) KINT(2,2),SNGL(beta),SNGL(Q2),XPQ
         ENDIF
      ENDIF
c check for exclusive heavy flavor production
      if(IHF.GT.0) then
         do i=-6,6
            if(iabs(i).ne.IHFLA) xpq(i)=0
         enddo
      endif
c      write(6,*) ' PQPM x,q2,xpq     ',ihf,ihfla,x,q2,(xpq(i),i=0,5)

ccc      write(6,*) ' PQPM x,q2,xpq MXcut',x,q2,(xpq(i),i=0,4)
      F2 = 0.D0
      ALPH_EM = ALPH
      IF(IRUNAEM.EQ.1) ALPH_EM = PYALEM(Q2)
      IF(INTER.LT.2) THEN
         DO 10 I=-NFLAV,NFLAV
   10    F2 = F2 + DBLE(XPQ(I))*DFLOAT(PYCHGE(I))**2/9.D0
cccc         F2 = 2*DBLE(XPQ(4))*DFLOAT(PYCHGE(4))**2/9.D0
         SIGMA = 2.D0*PI*ALPH_EM**2 /YY1 /Q2/Q2 *(1.d0+(1.d0-yy1)**2)*
     +   F2
      ELSEIF(INTER.EQ.2) THEN
         WTGLU = 0.D0
         DO 20  I=1,NFLAV-1,2
            IF(ISIGN(1,K(1,2)).EQ.1) THEN
               WTGLU = WTGLU + (1.D0 - YY1 + YY1**2/2.D0)*2.D0*
     +         DBLE(XPQ( -I)+XPQ(I+1)) + (YY1 - YY1**2/2.D0)*2.D0*
     +         DBLE(-XPQ(-I)+ XPQ(I+1))
            ELSEIF(ISIGN(1,K(1,2)).EQ.-1) THEN
               WTGLU = WTGLU + (1.D0 - YY1 + YY1**2/2.D0)*2.D0*
     +         DBLE(XPQ(I )+XPQ(-I-1)) - (YY1 - YY1**2/2.D0)*2.D0*
     +         DBLE(XPQ(I)-XPQ(- I-1))
            ENDIF
   20    CONTINUE
         GF = PI*ALPH_EM/(SIN2W*XMW2*DSQRT(2.D0))
         WMAT = GF**2/2.D0/PI /(1.D0 + Q2/XMW2)**2 /YY1
capply additional factor 0.5 because only right/or left handed electrons contr.
         WMAT = WMAT * 0.5D0
         SIGMA = WMAT * WTGLU
c         write(6,*) ' here in CC pqpm',sigma*gev2nb
      ENDIF
c      write(6,*) pi,alph,q2,f2,yy1,q2
      PQPM = SIGMA * GEV2NB
      RETURN
      END
