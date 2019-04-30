
C->
       DOUBLE PRECISION FUNCTION GAUSS (FUN, A,B)
c*********************************************************
C	Returns the  8 points Gauss-Legendre integral
C	of function FUN from A to B
c       this routine was provided by T.Stanev
c*********************************************************
c** Date: 20/01/98   **
c** A.Muecke         **
c**********************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      SAVE

      EXTERNAL FUN

C...........................................................
	DIMENSION X(8), W(8)
	DATA X /.0950125098D0,.2816035507D0,.4580167776D0,.6178762444D0
     +         ,.7554044083D0,.8656312023D0,.9445750230D0,.9894009349D0/
	DATA W /.1894506104D0,.1826034150D0,.1691565193D0,.1495959888D0
     +        ,.1246289712D0,.0951585116D0,.0622535239D0, .0271524594D0/
        if (A-1.D-10.le.B.and.A+1.D-10.ge.B) then
         GAUSS = 0.D0
c         print*,'A=B in integration routine'
         RETURN
        endif 
	XM = 0.5D0*(B+A)
	XR = 0.5D0*(B-A)
	SS = 0.D0
	DO NJ=1,8
	  DX = XR*X(NJ)
	  SS = SS + W(NJ) * (FUN(XM+DX) + FUN(XM-DX))
	ENDDO
	GAUSS = XR*SS
	RETURN
	END
