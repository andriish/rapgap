      SUBROUTINE FLUXH1(X_POM,T2,AP,B0,C,P,FLUX)
*#########################################################
* Generic flux template
*#########################################################
      IMPLICIT REAL*8 (A-G,O-Z)

      B = B0 + 2.D0*AP*LOG(1.D0/X_POM)

      BT=B*DABS(T2)
      IF(BT.GT.170D0) THEN
         FLUX=0.D0
      ELSE
         FLUX=C*(X_POM**(-P))*DEXP(-BT)
      ENDIF

      RETURN
      END
