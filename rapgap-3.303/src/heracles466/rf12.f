
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C

      SUBROUTINE RF12(TB,HM,RF1,RF2)
      COMMON /HSSINC/ S,AP,AMP2,W2PIT,W2MIN,W2TR
C
      CW=3471.16/389383.4
C
      IF(HM-1.6)1,3,3
1     IF(TB-2.0)2,3,3
2     SRF=0D0
      GO TO 4
3     SRF=.18
4     R=SRF
      XB=S-TB+AMP2-HM
      SN=TB*AMP2/S/XB
      CN=1.-SN
      TN=SN/CN
      ANB=(S-XB)/AP
      ANB2=ANB**2
      VEPS=1./(1.+2.*TN*(1.+ANB2/TB))
      IF(HM-W2MIN)5,5,6
5     STOT=((HM-W2PIT)/(W2MIN-W2PIT))*
     &     SVTOT(TB,W2MIN,VEPS)
      GO TO 7
   6  STOT=SVTOT(TB,HM,VEPS)
7     RF1 =CW*(HM-AMP2)/(1.+VEPS*R)*STOT
      RF2=ANB/AP*TB/(TB+ANB2)*(1.+R)*RF1
      END
