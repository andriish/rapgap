C----------------------------------------------------------------------
      FUNCTION THREEDOT(P,Q)
      DOUBLE PRECISION THREEDOT
C---Galilei 3-vector   DOT PRODUCT
      DOUBLE PRECISION P(3),Q(3)
      THREEDOT=P(1)*Q(1)+P(2)*Q(2)+P(3)*Q(3)
      END
