CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        FUNCTION HSSPEN(Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SPENCE-FUNKTION KOMPLEX, FREI NACH HOLLIK
C-----------------------------------------------------------------------
C       20.07.83    UE 20.08.85
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        COMPLEX*16 HSSPEN,HSCLN,W,SUM,Z,U
        DOUBLE PRECISION RZ,AZ,A1
        DOUBLE PRECISION B(9)/
     1   0.1666666666666666666666666667D0,
     2  -0.0333333333333333333333333333D0,
     3   0.0238095238095238095238095238D0,
     4  -0.0333333333333333333333333333D0,
     5   0.0757575757575757575757575758D0,
     6  -0.2531135531135531135531135531D0,
     7   1.1666666666666666666666666667D0,
     8  -7.09215686274509804D0         ,
     9  54.97117794486215539D0         /
C     BEACHTE:                 B(N)=B2N
C     B(1)=1./6.
C     B(2)=-1./30.
C     B(3)=1./42.
C     B(4)=-1./30.
C     B(5)=5./66.
C     B(6)=-691./2730.
C     B(7)=7./6.
C     B(8)=-3617./510.
C     B(9)=43867./798.
C     B(10)=-174611./330.
C     B(11)=854513./138.
C     PI=3.1415926535897932384
C     PI*PI/6.=1.6449..., PI*PI/3=3.28986...
C
      RZ=DREAL(Z)
      AZ=ABS(Z)
      A1=ABS(1D0-Z)
CH
CH FOR VS FORTRAN:
      IF (AZ.LT.1D-15) THEN
         HSSPEN=DCMPLX(0D0,0D0)
         RETURN
      ENDIF
CH
      IF((RZ .EQ. 1D0) .AND. (DIMAG(Z) .EQ. 0D0)) GOTO 40
      IF(RZ.GT.5D-1) GOTO 20
      IF(AZ.GT.1D0) GOTO 10
      W=-HSCLN(1D0-Z)
      SUM=W-0.25D0*W*W
      U=W
      DO 1 K=1,9
      U=U*W*W/FLOAT(2*K*(2*K+1))
      SUM=SUM+U*B(K)
 1    CONTINUE
      HSSPEN=SUM
      RETURN
10    W=-HSCLN(1D0-1D0/Z)
      SUM=W-0.25*W*W
      U=W
      DO 11 K=1,9
      U=U*W*W/DFLOAT(2*K*(2*K+1))
      SUM=SUM+U*B(K)
11    CONTINUE
      HSSPEN=-SUM-1.64493406684822643D0-.5*HSCLN(-Z)**2
      RETURN
20    IF(A1.GT.1D0) GOTO 30
      W=-HSCLN(Z)
      SUM=W-0.25*W*W
      U=W
      DO 21 K=1,9
      U=U*W*W/FLOAT(2*K*(2*K+1))
      SUM=SUM+U*B(K)
21    CONTINUE
      HSSPEN=-SUM+1.64493406684822643D0-HSCLN(Z)*HSCLN(1.-Z)
      RETURN
30    W=HSCLN(1D0-1./Z)
      SUM=W-0.25*W*W
      U=W
      DO 31 K=1,9
      U=U*W*W/FLOAT(2*K*(2*K+1))
      SUM=SUM+U*B(K)
31    CONTINUE
      HSSPEN=SUM+3.28986813369645287D0+.5*HSCLN(Z-1.)**2
     *      -HSCLN(Z)*HSCLN(1.-Z)
50    CONTINUE
      RETURN
40    HSSPEN=1.64493406684822643D0
      RETURN
      END