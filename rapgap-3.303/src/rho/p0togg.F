      SUBROUTINE P0TOGG( JPI0, PPI0 )
C--------------------------------------------------------------------
C       LET PI0 DECAY INTO GAMMA GAMMA
C
C  INPUT
      INTEGER JPI0
      Double Precision EG0,MP0
      COMMON/DECA2/PG1,PG2
C                                       ! INDEX OF PI0
c       JPI0 IS 1 IN OUR CASE
      REAL*8 PPI0(4)
C                                       ! FOUR MOMENTUM OF PI0 IN LAB
C  OUTPUT
C       /MCLIST/
C  VAR
      integer k 
      REAL VEC(3)
      real*8 PG1(4), PG2(4)
      DATA EG0,MP0/0.0674867,0.1349734/
C  BEGIN
C
C     GENERATE GAMMAS IN THE PI-SYSTEMS..........................
      CALL GENSPH(VEC)
C
C...EG0 IS THE MASS OF PIZERO /2.
      PG1(4)=EG0
      PG2(4)=EG0

      DO 10 K=1,3
         PG1(K) = EG0 * DBLE(VEC(K))
         PG2(K) = -PG1(K)
   10 CONTINUE
C
C     BOOST THE GAMMAS TO THE LAB SYSTEM..................
      CALL DBOOST(PPI0,MP0,PG1,PG1)
      CALL DBOOST(PPI0,MP0,PG2,PG2)


      RETURN
      END
