C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C   Interface to the parametrization of BK for F2, F1
C   B.Badelek and J.Kwiecinski, Nucl.Phys. B295 (1992) 263
C   See comments in HSF2BK below
C
      SUBROUTINE HSSTBK(DX,DQ2,ZF1,ZF2)
      DOUBLE PRECISION DX,DQ2,ZF1,ZF2
      LOGICAL LFIRST
      DATA LFIRST /.TRUE./
c     Example of usage of JKBB F2 routine 'f2pd';
c     for the conditions of usage read the instruction in f2pd.
c     When implementing f2pd in your program use "DATA" statement
c     instead of "READ" to store matrix h(2,20,20)
C     F2PD RENAMED TO HSF2BK (H.S.)
c
      common//h(2,20,20)
c
c     choose "mode, x, Q2" and read in the data from the "mode" file
c
C     PRINT*,'MODE, X, Q2='
C     READ*,MODE,X,Q2
      IF (LFIRST) THEN
      LFIRST=.FALSE.
      lode=21
      mode=3
      CALL HSBKIN
C     DO 200 I=1,2
C     DO 200 IQ=1,20
C     READ(MODE,*)H(I,IQ, 1),H(I,IQ, 2),H(I,IQ, 3),H(I,IQ, 4),H(I,IQ, 5)
C     READ(MODE,*)H(I,IQ, 6),H(I,IQ, 7),H(I,IQ, 8),H(I,IQ, 9),H(I,IQ,10)
C     READ(MODE,*)H(I,IQ,11),H(I,IQ,12),H(I,IQ,13),H(I,IQ,14),H(I,IQ,15)
C     READ(MODE,*)H(I,IQ,16),H(I,IQ,17),H(I,IQ,18),H(I,IQ,19),H(I,IQ,20)
C200   CONTINUE
      ENDIF
c
c     and calculate proton & deuteron F2 (f2p, f2d respectively)
c     for a defined x,Q2 by calling the subroutine f2pd
c
      X=DX
      Q2=DQ2
      CALL HSF2BK(MODE,X,Q2,F2P,F2D)
      ZF2=DBLE(F2P)
      ZF1=ZF2/2D0/DX
      RETURN
      end
