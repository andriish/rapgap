      SUBROUTINE BLEPIN

*#**********************************************************************
*#
*# to initialise block data from phythia and lepto neede for
*# structure function and parton shower
*#
*# CALLING: UCOPY
*#
*# AUTHOR:                          CREATED AT:
*#
*# CHANGED BY:                              AT:
*# REASON:
*#
*#**********************************************************************
      IMPLICIT NONE
      CHARACTER*23  FLAG,CHECK

C...Give sensible default values to switches and parameters.
      INTEGER I
      DOUBLE PRECISION XARL(30)
	REAL XUT(14)
	Integer KST(40)
      DOUBLE PRECISION PARL
	REAL CUT,X,Y,W2,Q2,U
	Integer LST
      COMMON /RAPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      Double Precision PYPAR,PYVAR
      Double Precision XYPAR(80),XYVAR(80)
      Integer IPY
      Integer KPY(80)
      COMMON /RGPYPARA/ IPY(80),PYPAR(80),PYVAR(80)
 
C...RAPTOU: Cuts, basic switches and parameters.
      DATA XUT/1.E-04,1.,0.,1.,4.,1.E+08,5.,1.E+08,1.,1.E+08,1.,1.E+08,
     +0.,3.1416/
      DATA KST/0,1,5,1,3,1,1,12,5,1,0,4,5,1,1,1,0,2,3,21*0/
**      DATA PARL/1.,1.,0.44,0.75,0.226,0.,0.,0.015,2.,0.,0.01,4.,
**     +0.1,0.44,0.01,7.29735E-03,1.16637E-05,0.044,0.03,1.,10*0./
      DATA XARL/1.,0.7,0.44,0.75,0.226,0.,0.,0.015,2.,0.,0.01,4.,
     +0.1,0.44,0.01,7.29735E-03,1.16637E-05,0.044,0.03,1.,10*0./
      DATA KPY/
     + 0,     0,     2,     2,     6,     1,     1,     4,     3,     1,
     + 3,     1,     1,     2,     1,     1,     4,     1,     1,     1,
     + 0,     1,     1,     1,     1,     1,     1,     0,     0,     0,
     + 1,     2,     1,     1,    30,    33,     1,     1,     7,     0,
     + 0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     + 0,     0,     0,     1,   100,     0,     0,     0,     0,     0,
     + 0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     + 0,     0,     0,     0,     0,     0,     0,     0,     0,     0/
      DATA (XYPAR(I),I=1,40)/
     +   7.299E-03,   2.290E-01,   2.000E-01,   2.500E-01,   4.000E+00,
     +   1.000E+00,   4.400E-01,   4.400E-01,   7.500E-02,   0.000E+00,
     +   2.000E+00,   2.000E+00,   1.000E+00,   0.000E+00,   3.000E+00,
     +   1.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   1.000E+00,
     +   2.500E-01,   1.000E+00,   2.000E+00,   1.000E-03,   1.000E+00,
     +   1.000E+00,   1.000E+00,  -2.000E-02,  -1.000E-02,   0.000E+00,
     +   0.000E+00,   1.600E+00,   0.500E+00,   0.200E+00,   3.894E-01,
     +   1.000E+00,   3.300E-01,   6.600E-01,   0.000E+00,   1.000E+00/
      DATA (XYPAR(I),I=41,80)/
     +   2.260E+00,   1.000E+04,   1.000E-04,   0.000E+00,   0.000E+00,
     +   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     +   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     +   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     +   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     +   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     +   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     +   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00/
      DATA XYVAR/80*0./
      DATA CHECK /'RAPGAP 2.0 INITIALISED'/
      DATA FLAG  /'                      '/

      IF (FLAG .EQ. CHECK) THEN
         WRITE(6,*) ' BLOCK DATA for RAPGAP 2.0 REINITIALISED.'
      ENDIF

      FLAG = CHECK
      WRITE(6,*) ' NEW BLOCK DATA for RAPGAP  initialised'
      
      do i=1,14
      CUT(I) = XUT(I)
      end do
      do i=1,40
      LST(I) = KST(I)
      end do
      do i=1,30
      PARL(I) = XARL(I)
      end do
      
      do i=1,80
      IPY(I) = KPY(I)
      PYPAR(I) = XYPAR(I)
      PYVAR(I) = XYVAR(I)
      end do

      RETURN
      END
