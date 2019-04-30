C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C   INTERFACE FOR CALLS OF PARTON DISTRIBUTIONS FROM HERACLES
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSPVER(X,Q2)
      DOUBLE PRECISION QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
      DOUBLE PRECISION X,Q2
      COMMON /HSPDFQ/ QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      DIMENSION XPQ(-6:6)
C---
      RX=X
      RQ2=Q2
C---
C...USE PYSTFU ROUTINES via LYSTFU from LEPTO 6.5
      CALL LYSTFU(2212,RX,RQ2,XPQ)
      QU=XPQ(2)/X
      QD=XPQ(1)/X
      QS=XPQ(3)/X
      QC=XPQ(4)/X
      QB=XPQ(5)/X
      QT=XPQ(6)/X
      QBU=XPQ(-2)/X
      QBD=XPQ(-1)/X
      QBS=XPQ(-3)/X
      QBC=XPQ(-4)/X
      QBB=XPQ(-5)/X
      QBT=XPQ(-6)/X

      RETURN
      END
