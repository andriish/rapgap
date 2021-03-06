      Subroutine PARTKINE
      IMPLICIT NONE
c
C modified by H. Kowalski 11. Oct 1999
c   author    H.Kowalski
c   Reconstruction of parton four vectors from dynamic (generated) variables
c   Q2,Mx2,W2 and Kt2 (for 2 Body FiPaSt) or QtG2,SM2,KtQ2 (for 3 Body FiPaSt)
c   Results are stored in Common LUJETS
c
      Double Precision  Stot, Q2, W2, Mx2, AJAC
      COMMON /CQ2W2MX/  Stot, Q2, W2, Mx2, AJAC
      Double Precision  YBj, XBj, Xpom, Beta, Tdf
      COMMON /GDVARB1/  YBj, XBj, Xpom, Beta, Tdf
      Double Precision    PGAMCM   , PPOMCM,   PGPOMCM
      Common /CINPKINE/   PGAMCM(5), PPOMCM(5),PGPOMCM(5)
      Double Precision   Qpar2BL   ,QBpar2BL,   Qpar2BC   ,QBpar2BC
      Common /C2BOKINE/  Qpar2BL(4),QBpar2BL(4),Qpar2BC(4),QBpar2BC(4)

      Double Precision    Kpar3BL   , Qpar3BL    ,QBpar3BL
      Double Precision    Kpar3BC   , Qpar3BC    ,QBpar3BC
      Common /C3BOKINE/   Kpar3BL(4), Qpar3BL(4) ,QBpar3BL(4),
     +                    Kpar3BC(4), Qpar3BC(4) ,QBpar3BC(4)

      Double Precision    Kt2, WT2
      COMMON /C2BODVAR/   Kt2, WT2
      Double Precision    KtG2, QtG2, SM2, ZW, KW2, WT3, YCUT
      COMMON /C3BODVAR/   KtG2, QtG2, SM2, ZW, KW2, WT3, YCUT


      Integer QFKF
      Double Precision    MQUARK, MCHARM
      COMMON /CQMASS/     MQUARK, MCHARM, QFKF

      Double Precision    LIMKIND
      Integer             IKIND,      Nevent, Ievent
      Common /CKIND/      LIMKIND(20),IKIND , Nevent, Ievent

#include "rglujets.inc"
#include "rgluco.inc"

      Double Precision  RANUNI, MQU2, MQU
      DATA MQU /0.01/, MQU2 /0.0001/
      Double Precision  TWOPI, Mx
      Double Precision QFRAN, SIGQFRAN

      Integer I, J,  IL, QF, Ip, SIGQF
      Real ECM
      Integer NJOIN, IJOIN(3)

      LOGICAL DEBUG

      DEBUG = .false.
      TWOPI = 2.*3.1415927

c
c store gamma* and pomeron vector
c
      IL = NIA1+2
      Mx = dsqrt(Mx2)

c         Write (*,10003) IL, NIA1,IKIND, MQUARK
10000    Format('PARTKINE IL',3I6,f9.3)

c      Call PYLIST(1)
***************************************************************
c Choice of the quark flavor and its mass
      If(Mx.le.1.) THEN
         QFRAN = RANUNI(0.D0,1.D0)
         QF = 1
         if(QFRAN.le.0.8) QF=2
         MQU = 0.
      EndIf

      If(Mx.gt.1.) THEN
         QFRAN = RANUNI(0.D0,1.D0)
         QF = 3
         if(QFRAN.le.0.66666) QF=2
         if(QFRAN.gt.0.66666.and.QFRAN.le.0.833333) QF=1
         MQU = 0.
      EndIf

      If(IKIND.gt.10) Then
         QF = 4
         MQU = MCHARM
      EndIf
c End of Choice of the quark flavor and its mass
*************************************************************

      MQU2 = MQU**2

      SIGQFRAN = RANUNI(0.D0,1.D0)
      SIGQF = 1
      if(SIGQFRAN.le.0.5) SIGQF = -1

      If(IKIND.eq. 3.or.IKIND.eq. 4) Goto 10
      If(IKIND.eq.13.or.IKIND.eq.14) Goto 10

      ip = il

      If(IKIND.lt.10) Then
         Do i = 1,4
            P(Ip,  i) =   QBpar2BL(i)
            P(Ip+1,i) =   Qpar2BL(i)
         EndDo
      EndIf
      If(IKIND.gt.10) Then
         Do i = 1,4
            P(Ip,  i) =   QBpar2BC(i)
            P(Ip+1,i) =   Qpar2BC(i)
         EndDo
      EndIf
      P(Ip,  5) =  MQU
      P(Ip+1,5) =  MQU

      K(Ip,1)=   2
      K(Ip,2)=  -SIGQF*QF
      K(Ip+1,1)= 1
      K(Ip+1,2)= SIGQF*QF

c added hannes
      KPA = K(Ip,2)
c end added hannes

C prepare string connections for PYJOIN
      DO J=1,5
         P(N+1,J) = P(Ip,J)
         K(N+1,J) = K(Ip,J)
         P(N+2,J) = P(Ip+1,J)
         K(N+2,J) = K(Ip+1,J)
      Enddo
c set these two entrys as documentation lines
      K(Ip,1)=21
      K(Ip+1,1)=21

      if(debug) Call PYLIST(1)
      NJOIN = 2
      IJOIN(1) = N+1
      IJOIN(2) = N+2
      N=N+2
C End of preparation of string connections for PYJOIN

      Call PYJOIN(NJOIN,IJOIN)

      If(DEBUG) THEN
         Write (*,10100) N, IKIND, Q2, sqrt(W2), sqrt(Mx2), sqrt(KT2)
10100 Format('PARTKINE 2BODY- Q2, W, Mx Kt',2I5,4F8.1/)
c        Write (*,1002) (PGAMCM(i),  i=1,5)
c        Write (*,1002) (PPOMCM(i),  i=1,5)
c        Write (*,1002) (PGPOMCM(i),  i=1,5)
10200 Format('   ',5F10.3)

         Call PYLIST(1)
c      Write(*,*) ' '
      EndIf

      Return

   10 Continue
c     3 parton process
c
      ip = il

      If(IKIND.lt.10) Then
         Do i = 1,4
            P(Ip,  i) =   QBpar3BL(i)
            P(Ip+1,i) =   Kpar3BL(i)
            P(Ip+2,i) =   Qpar3BL(i)
         EndDo
      EndIf
      If(IKIND.gt.10) Then
         Do i = 1,4
            P(Ip,  i) =   QBpar3BC(i)
            P(Ip+1,i) =   Kpar3BC(i)
            P(Ip+2,i) =   Qpar3BC(i)
         EndDo
      EndIf

      P(Ip,  5) =  MQU
      P(Ip+1,5) =  0.
      P(Ip+2,5) =  MQU

      P(Ip,3)   = -P(Ip,3)
      P(Ip+1,3) = -P(Ip+1,3)
      P(Ip+2,3) = -P(Ip+2,3)

c
      K(Ip,1)= 2
      K(Ip,2)=-SIGQF*QF
      K(Ip+1,1)= 2
      K(Ip+1,2)= 21

      K(Ip+2,1)= 1
      K(Ip+2,2)= SIGQF*QF

c added hannes
      KPA = K(Ip,2)
c end added hannes

      ECM = sqrt(SM2)
      DO J=1,5
         P(N+1,J) = P(Ip,J)
         K(N+1,J) = K(Ip,J)
         P(N+2,J) = P(Ip+1,J)
         K(N+2,J) = K(Ip+1,J)
         P(N+3,J) = P(Ip+2,J)
         K(N+3,J) = K(Ip+2,J)
      Enddo
      K(Ip,1)=21
      K(Ip+1,1)=21
      K(Ip+2,1)=21

      If(DEBUG) THEN
         Write (*,10300) N, IKIND, Q2, sqrt(W2), sqrt(Mx2), sqrt(KT2)
10300    Format('PARTKINE 3BODY- Q2, W, Mx Kt',2I5,4F8.1/)
         write(6,*) ' QB ',(P(ip  ,i),i=1,4)
         write(6,*) ' K  ',(P(ip+1,i),i=1,4)
         write(6,*) ' Q  ',(P(ip+2,i),i=1,4)
         Call PYLIST(1)
c     Write(*,*) ' '
      EndIf


      NJOIN = 3
      IJOIN(1) = N+1
      IJOIN(2) = N+2
      IJOIN(3) = N+3
      N=N+3
      Call PYJOIN(NJOIN,IJOIN)

      if(debug) Call PYLIST(1)
      Return
      End
