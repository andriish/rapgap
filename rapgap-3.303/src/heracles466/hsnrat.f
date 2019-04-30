C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   NUCLEAR SHADOWING
C   ASSUME Q2-INDEPENDENT SHADOWING (ANTI-SHADOWING) FOR HEAVY NUCLEI
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSNRAT(X)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSNUCL/ HNA,HNZ
      LOGICAL LFIRST
      DATA LFIRST/.TRUE./
      DATA AN1/0.130D0/
      DATA AN2/0.456D0/
      DATA AN3/0.773D0/

      IF (LFIRST) THEN
        LFIRST=.FALSE.
        HNA3=HNA**(1D0/3D0)
        HMI=1D0-1D0/HNA3-1.145D0/HNA3/HNA3+0.93D0/HNA+0.88D0/HNA/HNA3
     *         -0.59D0/HNA/HNA3/HNA3
        HM1=HMI*AN1
        H1M2=1D0+HMI*AN2
        HM3=HMI*AN3
      ENDIF
      IF (HNA.EQ.2D0.AND.HNZ.EQ.1D0) THEN
        HSNRAT=1D0
       ELSE
        HSNRAT=X**HM1*(1D0-HM3*X)*H1M2
      ENDIF
      RETURN
      END
