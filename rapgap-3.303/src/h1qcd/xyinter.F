C========================================================================
C 2D Interpolation to point within square (path dependent but simple)
C========================================================================
      FUNCTION XYINTER(X1,X2,Y1,Y2,XT,YT,F11,F12,F21,F22)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      XEXP = (XT-X1)/(X2-X1)
      YEXP = (YT-Y1)/(Y2-Y1)

      XYINTER = F11 + XEXP*(F21-F11)
     &          + YEXP*(F12 + XEXP*(F22-F12) - (F11 + XEXP*(F21-F11)))
      RETURN
      END
