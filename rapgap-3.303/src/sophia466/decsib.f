

      SUBROUTINE DECSIB
C***********************************************************************
C
C   Decay all unstable particle in SIBYLL
C   decayed particle have the code increased by 10000
C
C   (taken from SIBYLL 1.7, R.E. 04/98)
C
C***********************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

       COMMON /S_CSYDEC/ CBR(102), IDB(49), KDEC(612), LBARP(49)
       COMMON /S_PLIST/ P(2000,5), LLIST(2000), NP, Ideb
      COMMON /S_PLIST1/ LLIST1(2000)
      SAVE

      DIMENSION P0(5), LL(10), PD(10,5)

      NN = 1
      DO J=1,NP
         LLIST1(J) = 0
      ENDDO
      DO WHILE (NN .LE. NP)
         L= LLIST(NN)
         IF (IDB(IABS(L)) .GT. 0)  THEN
            DO K=1,5
              P0(K) = P(NN,K)
            ENDDO
            ND = 0
            CALL DECPAR (L,P0,ND,LL,PD)
            LLIST(NN) = LLIST(NN)+ISIGN(10000,LLIST(NN))
            DO J=1,ND
               DO K=1,5
                  P(NP+J,K) = PD(J,K)
               ENDDO
               LLIST(NP+J)=LL(J)
               LLIST1(NP+J)=NN
            ENDDO
            NP=NP+ND
         ENDIF
         NN = NN+1
      ENDDO

      END
