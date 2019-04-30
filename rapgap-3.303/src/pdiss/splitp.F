      SUBROUTINE SPLITP (P, Q, PQ, PREMN, IERR)

      IMPLICIT NONE
************************************************************************
*
*     Purpose:  A quark in a hadron with 5-momentum P is struck
*     --------  by a boson with 5-momentum Q;
*               Q is to be negative for spacelike particles.
*               The hadron splits into
*               a remnant with 5-momentum (1-x)*P(5) and
*               a quark with mass P (5) and 3-momentum x*P(1:3)+Q(1:3).
*               x is determined so that 4-vector conservation holds.
*
*     Input:    P:         5-momentum of the hadron
*     ------    Q:         5-momentum of the boson
*
*     Output:   PQ:        5-vector of the quark
*     -------   PREMN:     5-vector of the remnant
*               IERR:      0 if everything is O.K.,
*                         -1 if x<0 or x>1.
*                          1 if no solution was found.
*
*     Author:   Benno List, 20.1.94
*     -------
*
************************************************************************
      Double Precision P (5), Q (5), PQ (5), PREMN (5)
      INTEGER IERR,I
      Double Precision  X, D, P2, Q2
      
      P2= P (5)*ABS (P (5))
      Q2= Q (5)*ABS (Q (5))
      D = (P(4)+Q(4))**2-(P(1)+Q(1))**2-(P(2)+Q(2))**2-(P(3)+Q(3))**2
     +     - P2 - Q2
      IF (D .EQ. 0D0) THEN
         IERR = 1
         PRINT *, '### SPLITP error: D = 0!'
         PRINT *, 'P = (', P(1),',',P(2),',',P(3),',',P(4),'), ',P(5)
         PRINT *, 'Q = (', Q(1),',',Q(2),',',Q(3),',',Q(4),'), ',Q(5)
      ELSE
         X = -Q2 / D
         IF (X .LT. 0D0 .OR. X .GT. 1D0) THEN
            IERR = -1
            write(6,*) 'SPLITP: x<0 or x>1 X,D',X,D
         ELSE
            IERR = 0
         END IF
         DO 10, I = 1, 4
            PREMN (I) = (1D0 - X)*P (I)
 10      CONTINUE
         PREMN (5) = ABS (1D0 - X)*P (5)
         PQ (5) = ABS (X)*P (5)
         PQ (4) = X**2 * P2
         DO 20, I = 1, 3
            PQ (I) = X*P (I) + Q (I)
            PQ (4) = PQ (4) + PQ (I)**2
 20      CONTINUE
         IF (PQ (4) .LT. 0D0) THEN
            IERR = 1
            write(6,*)'SPLITP: No solution found'
         END IF
         PQ (4) = SQRT (ABS (PQ (4)))

c         PRINT *, '   SPLITP: X =             ', X
c         PRINT *, '   SPLITP: (P+Q)**2 =      ', (P (4) + Q (4))**2
c     +        -(P (1)+Q (1))**2-(P (2)+Q (2))**2-(P (3)+Q (3))**2
c         PRINT *, '   SPLITP: (PQ+PREMN)**2 = ',(PQ (4) + PREMN (4))**2
c     +        -(PQ(1)+PREMN(1))**2-(PQ(2)+PREMN(2))**2
c     +        -(PQ(3)+PREMN(3))**2
c         PRINT *, '   SPLITP: PQ (5) =        ', PQ (5)
c         PRINT *, '   SPLITP: PREMN (5) =     ', PREMN (5)
      END IF

      RETURN
      END




