
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C

      FUNCTION SVTOT(TB,HM,VEPS)
      COMMON /HSWABC/ A1,B1,C1,A2,B2,C2,A3,B3,C3
      COMMON /HSSINC/ S,AP,AMP2,W2PIT,W2MIN,W2TR
C
      IF(VEPS-.9)2,1,1
1     A=A1
      B=B1
      C=C1
      GO TO 5
2     IF(VEPS-.6)4,4,3
3     A=A2
      B=B2
      C=C2
      GO TO 5
4     A=A3
      B=B3
      C=C3
5     GD2=1./(1.+TB/.71)**4
      ANB=(HM+TB-AMP2)/AP
      XMQ=SQRT(ANB**2+TB)
      XMQ0=(HM-AMP2)/AP
      AL=LOG(XMQ/XMQ0)
      SVTOT=GD2*EXP(A+B*AL+C*ABS(AL)**3)
      END
