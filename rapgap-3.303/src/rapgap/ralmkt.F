      SUBROUTINE RALMKT(S,KT,PHISPL)
      Implicit None
      REAL KT,KTMIN,KTMAX,PHISPL,S
      Double Precision draprn
      EXTERNAL draprn
      KTMAX = SQRT(S/4.)
      KTMIN = 0.0
      PHISPL = 0.0
      KT = -LOG(draprn()*(EXP(-5.5*KTMAX) - EXP(-5.5*KTMIN))
     &  + EXP(-5.5*KTMIN))
      KT = KT/5.5
c      write(6,*) ' RALMKT ',KT,S,Ktmin,ktmax
      RETURN
      END
