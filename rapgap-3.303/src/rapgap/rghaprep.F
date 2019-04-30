      subroutine rghaprep
      implicit none
#include "rgludat1.inc"
#include "rgparas.inc"
#include "rgrapgki.inc"
#include "rgluco.inc"
#include "rgdisdif.inc"
      Real     remPARJ32
      Real     remPARJ(5)
      LOGICAL FIRST
      DATA FIRST/.TRUE./
c henry's commons
      Double Precision  Stot, Q2_henry, W2, Mx2, AJAC
      COMMON /CQ2W2MX/  Stot, Q2_henry, W2, Mx2, AJAC
      real rgmx2
      IF(FIRST) THEN
         first = .FALSE.
         remPARJ32 = PARJ(32)
      ELSE
         PARJ(32)=remPARJ32
      ENDIF
      if(ipro.eq.30) then
         rgmx2 = real(mx2)
      else
         rgmx2 = shh
      endif
      If (RGMX2.LT.4.0.AND.IDIR.EQ.0) Then
         If(ABS(KPA).LE.2) Then
            PARJ(32) = 0.35
            GOTO 10
         ElseIf(ABS(KPA).EQ.3) Then
            PARJ(32) = 0.1
            GOTO 10
         EndIf
   10    Continue
*
         remPARJ(1) = PARJ(11)
         remPARJ(2) = PARJ(12)
         remPARJ(3) = PARJ(13)
         remPARJ(4) = PARJ(27)
         remPARJ(5) = PARJ(28)
*
         PARJ(11) = 1.0
         PARJ(12) = 1.0
         PARJ(13) = 1.0
         PARJ(27) = 0.1
         PARJ(28) = 0.9
*
         CALL PYPREP(0)
         PARJ(11) = remPARJ(1)
         PARJ(12) = remPARJ(2)
         PARJ(13) = remPARJ(3)
         PARJ(27) = remPARJ(4)
         PARJ(28) = remPARJ(5)
c      ElseIf (RGMX2.GT.9.0.AND.RGMX2.LT.14.AND. ABS(KPA) .EQ.4.AND.IDIR
c     +.EQ.0) Then
      ElseIf (RGMX2.GT.9.0.AND. ABS(KPA) .EQ.4.AND.IDIR
     +.EQ.0) Then
         remPARJ(2) = PARJ(12)
         remPARJ(3) = PARJ(13)
         PARJ(13) = 1.0
         CALL PYPREP(0)
         PARJ(13) = remPARJ(3)
         PARJ(12) = remPARJ(2)

      Else

         CALL PYPREP(0)


      ENDIF
      return
      end
