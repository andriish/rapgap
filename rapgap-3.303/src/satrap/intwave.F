      Subroutine  INTWAVE
c
c author H.Kowalski
c Integration of the wave function
c
      Implicit NONE
      Double Precision  Stot, Q2, W2, Mx2, AJAC
      COMMON /CQ2W2MX/  Stot, Q2, W2, Mx2, AJAC
      Double Precision  YBj, XBj, Xpom, Beta, Tdf
      COMMON /GDVARB1/  YBj, XBj, Xpom, Beta, Tdf

      Double Precision  R02, LAM, SIGM0, X0, ALPHAS
      COMMON /SATURM/   R02, LAM, SIGM0, X0, ALPHAS

      Double Precision    Kt2, WT2
      COMMON /C2BODVAR/   Kt2, WT2

      Double Precision  Lt2min, Lt2max

      Double Precision  IntWUT, INTWUL, INTWUG
      COMMON /CINTWAVE/ IntWUT, INTWUL, INTWUG

      Double Precision  SWAVET, SWAVEL, SWAVEG
      External          SWAVET, SWAVEL, SWAVEG

      Double Precision  DGQUAD
      External          DGQUAD

      Integer         NQUAD
      COMMON /CNQUAD/ NQUAD

      Double Precision  Mx

      Integer NQUA
      Mx = dsqrt(Mx2)

      Lt2min = 0.0
      Lt2max = 10.0/R02

      NQUA   = NQUAD
      IntWUT = DGQUAD(SWAVET,Lt2min,Lt2max,NQUA)
      IntWUL = DGQUAD(SWAVEL,Lt2min,Lt2max,NQUA)
      IntWUG = 0.

      IF(Mx.lt.1.) Return
      IntWUG = DGQUAD(SWAVEG,Lt2min,Lt2max,NQUA)

c      Write(*,1000) IntWUT, IntWUL, IntWUG
10000 Format('INTWAVE',3E10.2)

      Return
      END

