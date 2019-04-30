
      FUNCTION rggamm(X)

C...Gives ordinary Gamma function Gamma(x) for positive, real arguments;
C...see M. Abramowitz, I. A. Stegun: Handbook of Mathematical Functions
C...(Dover, 1965) 6.1.36.
      Double Precision B
      DIMENSION B(8)
      DATA B/-0.577191652D0,0.988205891D0,-0.897056937D0,0.918206857D0,
     +-0.756704078D0,0.482199394D0,-0.193527818D0,0.035868343D0/

      NX=INT(X)
      DX=X-NX

      rggamm=1.
      DXP=1.
      DO 10  I=1,8
         DXP=DXP*DX
   10 rggamm=rggamm+SNGL(B(I))*DXP
      IF(X.LT.1.) THEN
         rggamm=rggamm/X
      ELSE
         DO 20 IX=1,NX-1
   20    rggamm=(X-IX)*rggamm
      ENDIF
      RETURN
      END
