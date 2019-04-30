C=======================================================================
      FUNCTION HSCLM1(Q2,MF2)
C       PHOTONIC VERTEX CORRECTION
C       IR-FINITE PART
C                                 OHNE LOG**2
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX*16 HSCLM1
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      IF (Q2) 10,20,30
10    AMB1 = DLOG(-Q2/MF2)
     1         +4.*(PI*PI/12. - 1D0)
      HSCLM1=DCMPLX(AMB1)
      GOTO 40
20    HSCLM1 = (0D0,0D0)
      GOTO 40
30    AMB1 = DLOG(Q2/MF2)
     1         +4.*(PI*PI/3. - 1D0)
      HSCLM1=DCMPLX(AMB1)
40    RETURN
      END
