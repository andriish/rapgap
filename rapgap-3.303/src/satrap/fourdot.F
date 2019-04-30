C----------------------------------------------------------------------
      FUNCTION FourDot(P,Q)
      DOUBLE PRECISION FourDot
C---LORENTZ 4-VECTOR DOT PRODUCT
      DOUBLE PRECISION P(5),Q(5)
      FourDot=P(4)*Q(4)-(P(1)*Q(1)+P(2)*Q(2)+P(3)*Q(3))
      END
*
