
      SUBROUTINE GENSPH(VEC)
C--------------------------------------------------------C
C GENERATES AN UNIT-LENGTH VECTOR UNIFORMLY IN 4PI.      C
C USES 1 SQRT AND NO SIN,COS'S. HITOSHI OCT-86           C
C--------------------------------------------------------C
      REAL VEC(3)
C
      Double Precision RNDMV(2)
   10 CALL draprnV(RNDMV,2)
      U1=2.*RNDMV(1)-1.
      U2=2.*RNDMV(2)-1.
      S=U1*U1+U2*U2
      IF(S.GE.1.) GOTO 10
C
      CO2=2.*SQRT(1.-S)
      VEC(1)=U1*CO2
      VEC(2)=U2*CO2
      VEC(3)=1.-2.*S
C
      RETURN
      END
