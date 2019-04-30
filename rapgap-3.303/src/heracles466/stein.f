C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE STEIN(Q2,W2,F1,F2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSRADC/ AMP , AMP2, PI , PI2 , ALFA, AML , AML2
      COMMON /HSRADF/ AP,AP2,AP4,AL2,AL4,ALPS,CALPI,W2PIT,W2TR,TTR
      DIMENSION AN(5)
      DATA AN /1.0621D0, -2.2594D0, 10.54D0, -15.8277D0, 6.7931D0/
C
      X=Q2/(Q2+W2-AMP2)
      X1 = 1D0-X
C
C     STEIN'S FIT OF W2 IN THE LOW X REGION
C
      OS=1D0+W2/Q2
      GE2=1D0/((1D0+.61D0*Q2)*(1D0+2.31D0*Q2)*(1D0+.04D0*Q2))**2
      GM2 = 7.7841d0*GE2
      TAU=Q2/AP4
      W2EL=(GE2+TAU*GM2)/(1D0+TAU)
      ANU=(W2+Q2-AMP2)/AP
      SR=0D0
      DO 2 I=1,5
2     SR=SR+AN(I)*(1D0-1D0/OS)**(I+2)
      F2=(1D0-W2EL)*SR
C
C     F1 = 2.*AMP*W1
C
C   SLAC VALUES OF R
C
      R=.18D0
C
      ANU2=ANU**2
      F1=AP*(1D0+ANU2/Q2)/ANU/(1D0+R)*F2
      END
