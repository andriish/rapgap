      SUBROUTINE COLRES

      IMPLICIT NONE
#include "rglujets.inc"
#include "rgrapgki.inc"
#include "rgpara1.inc"
#include "rgparas.inc"
#include "rgpart.inc"
#include "rgluco.inc"
#include "rgdisdif.inc"
#include "rgludat1.inc"
#include "rgcolcon.inc"
      INTEGER LST,IRES
      COMMON/EPPARA/LST(30),IRES(2)
c..hju

      INTEGER KS,NS
      COMMON /COLR/ KS(40,5),NS


      K(NS+7,2)=21
ccc      write(6,*) ' colres IRESPRO ',irespro
      if(Ipro.eq.500) then
            K(NS+3,4) = 21
            K(NS+3,5) = 21
            K(NS+5,4)=0
            K(NS+5,5)=0
            K(NS+6,1)=0
            K(NS+6,2)=0
            K(NS+6,3)=0
            K(NS+6,4)=0
            K(NS+6,5)=0
            K(NS+1,4) = 23
            K(NS+1,5) = 23
            K(NS+7,4)=0
            K(NS+7,5)=0
      
      elseif(ipro.eq.501) then
            K(NS+5,1)=1
            K(NS+5,4)=0
            K(NS+5,5)=0
c      write(6,*) ' colorflow Zq ',ipro
c      call pylist(2)
c  q from proton 1
         IF(K(23,2).GT.0.AND.K(21,2).LT.0.and.K(23,2).le.10) THEN
            K(NS+1,4) = 27
            K(NS+1,5) = 27
            K(NS+3,4)=27
            K(NS+3,5)=27
ctest            K(NS+6,4)=(NS+3)*MSTU(5)
            K(NS+6,4)=(NS+7)*MSTU(5)
            K(NS+6,5)=(NS+7)*MSTU(5)
            K(NS+7,4)=(NS+3)*MSTU(5)+26
            K(NS+7,5)=(NS+1)*MSTU(5)+26
c		write(6,*) ' in  scenario 1a'
         ELSEIF(K(23,2).LT.0.AND.K(21,2).GT.0.and.K(21,2).le.10) THEN
            K(NS+1,4) = 27
            K(NS+1,5) = 27
            K(NS+3,4)=27
            K(NS+3,5)=27
ctest            K(NS+6,4)=(NS+3)*MSTU(5)
            K(NS+6,4)=(NS+7)*MSTU(5)
            K(NS+6,5)=(NS+7)*MSTU(5)
            K(NS+7,4)=(NS+1)*MSTU(5)+26
            K(NS+7,5)=(NS+3)*MSTU(5)+26
c		write(6,*) ' in  scenario 1b'
         ELSEIF(K(23,2).GT.0.AND.K(21,2).EQ.21.and.K(23,2).le.10) THEN
            K(NS+1,4) = 27
            K(NS+1,5) = 23
            K(NS+3,4)=21
            K(NS+3,5)=21
ctest            K(NS+6,4)=(NS+3)*MSTU(5)
            K(NS+6,4)=(NS+7)*MSTU(5)
            K(NS+6,5)=(NS+7)*MSTU(5)
            K(NS+7,4)=(NS+1)*MSTU(5)+26
            K(NS+7,5)=(NS+1)*MSTU(5)+26
            K(NS+7,2) = K(NS+6,2)
c		write(6,*) ' in scenario 1 '
c            call pylist(2)
c  anti-q from proton 1
         ELSEIF(K(23,2).LT.0.AND.K(21,2).EQ.21) THEN
            K(NS+1,4) = 23
            K(NS+1,5) = 27
            K(NS+3,4)=21
            K(NS+3,5)=21
ctest            K(NS+6,4)=(NS+3)*MSTU(5)
            K(NS+6,4)=(NS+7)*MSTU(5)
            K(NS+6,5)=(NS+7)*MSTU(5)
            K(NS+7,4)=(NS+1)*MSTU(5)+26
            K(NS+7,5)=(NS+1)*MSTU(5)+26
            K(NS+7,2) = K(NS+6,2)
c		write(6,*) ' in scenario 2 '
c            call pylist(2)
c  q from proton 2
         ELSEIF(K(21,2).GT.0.AND.K(23,2).EQ.21.and.K(21,2).le.10) THEN
            K(NS+3,4) = 27
            K(NS+3,5) = 21
            K(NS+1,4)=23
            K(NS+1,5)=23
c            K(NS+1,4)=0
            K(NS+6,4)=(NS+7)*MSTU(5)
            K(NS+6,5)=(NS+7)*MSTU(5)
c            K(NS+6,4)=0
            K(NS+7,4)=(NS+3)*MSTU(5)+26
            K(NS+7,5)=(NS+3)*MSTU(5)+26
c            K(NS+7,4)=0
            K(NS+7,2) = K(NS+6,2)
c		write(6,*) ' in  scenario 3'
c            call pylist(2)
c  anti-q from proton 2
         ELSEIF(K(21,2).LT.0.AND.K(23,2).EQ.21) THEN
            K(NS+3,4) = 21
            K(NS+3,5) = 27
            K(NS+1,4)=0
            K(NS+1,5)=23
            K(NS+6,4)=(NS+7)*MSTU(5)
c            K(NS+6,4)=0
            K(NS+6,5)=(NS+7)*MSTU(5)
            K(NS+6,5)=0
            K(NS+7,4)=(NS+3)*MSTU(5)+26
c            K(NS+7,4)=0
            K(NS+7,5)=(NS+3)*MSTU(5)+26
            K(NS+7,5)=0
            K(NS+7,2) = K(NS+6,2)
c		write(6,*) ' in scenario 4 '
c            call pylist(2)
         endif
c         write(6,*) ' after colorflow '
c         call pylist(2)
      elseif(IRESPRO.EQ.1) THEN
         K(NS+3,4)=21
         K(NS+3,5)=27
c qqbar event
         IF(K(NS+5,2).GT.0) THEN
            K(NS+5,4)=(NS+7)*MSTU(5)
            K(NS+5,5)=0
            K(NS+6,4)=0
            K(NS+6,5)=(NS+7)*MSTU(5)
            K(NS+1,4) = 27
            K(NS+1,5) = 23
            K(NS+7,4)=(NS+1)*MSTU(5)+25
            K(NS+7,5)=(NS+3)*MSTU(5)+26
         ELSE
            K(NS+6,4)=(NS+7)*MSTU(5)
            K(NS+6,5)=0
            K(NS+5,4)=0
            K(NS+5,5)=(NS+7)*MSTU(5)
            K(NS+1,4) = 27
            K(NS+1,5) = 23
            K(NS+7,4)=(NS+1)*MSTU(5)+26
            K(NS+7,5)=(NS+3)*MSTU(5)+25
         ENDIF
      ELSEIF(IRESPRO.EQ.2) THEN
c glu glu event

         IF(ICOLORA.EQ.1) THEN
c color flow A
            K(NS+1,4) = 27
            K(NS+1,5) = 23
            K(NS+3,4)=21
            K(NS+3,5)=27
            K(NS+5,4)=(NS+7)*MSTU(5)
            K(NS+5,5)=(NS+6)*MSTU(5)
            K(NS+6,4)=(NS+5)*MSTU(5)
            K(NS+6,5)=(NS+7)*MSTU(5)
            K(NS+7,4)=(NS+1)*MSTU(5)+25
            K(NS+7,5)=(NS+3)*MSTU(5)+26
         ELSEIF(ICOLORA.EQ.2) THEN
c color flow B
            K(NS+1,4) = 27
            K(NS+1,5) = 23
            K(NS+3,4)=21
            K(NS+3,5)=27
            K(NS+5,4)=(NS+6)*MSTU(5)
            K(NS+5,5)=(NS+7)*MSTU(5)
            K(NS+6,4)=(NS+7)*MSTU(5)
            K(NS+6,5)=(NS+5)*MSTU(5)
            K(NS+7,4)=(NS+1)*MSTU(5)+26
            K(NS+7,5)=(NS+3)*MSTU(5)+25
         ELSEIF(ICOLORA.EQ.3) THEN
c color flow C
            K(NS+1,4) = 26
            K(NS+1,5) = 25
            K(NS+3,4)=25
            K(NS+3,5)=26
            K(NS+5,4)=(NS+3)*MSTU(5)
            K(NS+5,5)=(NS+1)*MSTU(5)
            K(NS+6,4)=(NS+1)*MSTU(5)
            K(NS+6,5)=(NS+3)*MSTU(5)
            K(NS+7,4)=0
            K(NS+7,5)=0
            K(NS+7,1)=0
         ELSE
            write(6,*) ' wrong color configuartion ',ICOLORA
         ENDIF

      ELSEIF(IRESPRO.EQ.3) THEN
c qg event
C color configuartion A
         IF(ICOLORA.EQ.1) THEN
C antq from photon
            IF(K(21,2).LT.0.AND.K(23,2).EQ.21) THEN
               IF(K(25,2).LT.0.AND.K(26,2).EQ.21) THEN
                  K(NS+1,4) = 0
                  K(NS+1,5) = 23
                  K(NS+3,4)=21
                  K(NS+3,5)=27
                  K(NS+5,4)=0
                  K(NS+5,5)=(NS+6)*MSTU(5)
                  K(NS+6,4)=(NS+5)*MSTU(5)
                  K(NS+6,5)=(NS+7)*MSTU(5)
                  K(NS+7,4)=(NS+3)*MSTU(5)+25
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+5,2)
c               write(6,*) ' 1st '
               ELSEIF(K(25,2).EQ.21.AND.K(26,2).LT.0) THEN
                  K(NS+1,4) = 0
                  K(NS+1,5) = 23
                  K(NS+3,4)=21
                  K(NS+3,5)=27
                  K(NS+6,4)=0
                  K(NS+6,5)=(NS+5)*MSTU(5)
                  K(NS+5,4)=(NS+6)*MSTU(5)
                  K(NS+5,5)=(NS+7)*MSTU(5)
                  K(NS+7,4)=(NS+1)*MSTU(5)+26
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+6,2)
c               write(6,*) ' 3rd '
               ENDIF
c  q from photon
            ELSEIF(K(21,2).GT.0.AND.K(23,2).EQ.21) THEN
               IF(K(25,2).GT.0.AND.K(26,2).EQ.21) THEN
                  K(NS+1,4) = 23
                  K(NS+1,5) = 0
                  K(NS+3,4)=27
                  K(NS+3,5)=21
                  K(NS+5,4)=(NS+6)*MSTU(5)
                  K(NS+5,5)=0
                  K(NS+6,4)=(NS+7)*MSTU(5)
                  K(NS+6,5)=(NS+5)*MSTU(5)
                  K(NS+7,4)=(NS+3)*MSTU(5)+26
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+5,2)
c               write(6,*) ' 2nd '
               ELSEIF(K(25,2).EQ.21.AND.K(26,2).GT.0) THEN
                  K(NS+1,4) = 23
                  K(NS+1,5) = 0
                  K(NS+3,4)=27
                  K(NS+3,5)=21
                  K(NS+6,4)=(NS+5)*MSTU(5)
                  K(NS+6,5)=0
                  K(NS+5,4)=(NS+7)*MSTU(5)
                  K(NS+5,5)=(NS+6)*MSTU(5)
                  K(NS+7,4)=(NS+3)*MSTU(5)+26
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+6,2)
c               write(6,*) ' 4th '
               ENDIF
C antq from proton
            ELSEIF(K(23,2).LT.0.AND.K(21,2).EQ.21) THEN
               IF(K(25,2).LT.0.AND.K(26,2).EQ.21) THEN
                  K(NS+3,4) = 0
                  K(NS+3,5) = 21
                  K(NS+1,4)=23
                  K(NS+1,5)=27
                  K(NS+5,4)=0
                  K(NS+5,5)=(NS+6)*MSTU(5)
                  K(NS+6,4)=(NS+5)*MSTU(5)
                  K(NS+6,5)=(NS+7)*MSTU(5)
                  K(NS+7,4)=(NS+1)*MSTU(5)+25
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+5,2)
c               write(6,*) ' 1st '
               ELSEIF(K(25,2).EQ.21.AND.K(26,2).LT.0) THEN
                  K(NS+3,4) = 0
                  K(NS+3,5) = 21
                  K(NS+1,4)=23
                  K(NS+1,5)=27
                  K(NS+6,4)=0
                  K(NS+6,5)=(NS+5)*MSTU(5)
                  K(NS+5,4)=(NS+6)*MSTU(5)
                  K(NS+5,5)=(NS+7)*MSTU(5)
                  K(NS+7,4)=(NS+1)*MSTU(5)+26
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+6,2)
c               write(6,*) ' 3rd '
               ENDIF
c  q from proton
            ELSEIF(K(23,2).GT.0.AND.K(21,2).EQ.21) THEN
               IF(K(25,2).GT.0.AND.K(26,2).EQ.21) THEN
                  K(NS+3,4) = 21
                  K(NS+3,5) = 0
                  K(NS+1,4)=27
                  K(NS+1,5)=23
                  K(NS+5,4)=(NS+6)*MSTU(5)
                  K(NS+5,5)=0
                  K(NS+6,4)=(NS+7)*MSTU(5)
                  K(NS+6,5)=(NS+5)*MSTU(5)
                  K(NS+7,4)=(NS+1)*MSTU(5)+26
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+5,2)
c               write(6,*) ' 2nd '
               ELSEIF(K(25,2).EQ.21.AND.K(26,2).GT.0) THEN
                  K(NS+3,4) = 21
                  K(NS+3,5) = 0
                  K(NS+1,4)=27
                  K(NS+1,5)=23
                  K(NS+6,4)=(NS+5)*MSTU(5)
                  K(NS+6,5)=0
                  K(NS+5,4)=(NS+7)*MSTU(5)
                  K(NS+5,5)=(NS+6)*MSTU(5)
                  K(NS+7,4)=(NS+1)*MSTU(5)+25
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+6,2)
c               write(6,*) ' 4th '
               ENDIF
            ENDIF
C color configuartion B
         ELSE
C antq from photon
            IF(K(21,2).LT.0.AND.K(23,2).EQ.21) THEN
               IF(K(25,2).LT.0.AND.K(26,2).EQ.21) THEN
                  K(NS+1,4) = 0
                  K(NS+1,5) = 26
                  K(NS+3,4)=26
                  K(NS+3,5)=25
                  K(NS+5,4)=0
                  K(NS+5,5)=(NS+3)*MSTU(5)
                  K(NS+6,4)=(NS+3)*MSTU(5)
                  K(NS+6,5)=(NS+1)*MSTU(5)
                  K(NS+7,4)=0
                  K(NS+7,5)=0
                  K(NS+7,1)=0
                  K(NS+7,2) = K(NS+5,2)
c               write(6,*) ' 1st '
               ELSEIF(K(25,2).EQ.21.AND.K(26,2).LT.0) THEN
                  K(NS+1,4) = 0
                  K(NS+1,5) = 25
                  K(NS+3,4)=25
                  K(NS+3,5)=26
                  K(NS+6,4)=0
                  K(NS+6,5)=(NS+3)*MSTU(5)
                  K(NS+5,4)=(NS+3)*MSTU(5)
                  K(NS+5,5)=(NS+1)*MSTU(5)
                  K(NS+7,4)=0
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+6,2)
                  K(NS+7,1)=0
c               write(6,*) ' 3rd '
               ENDIF
c  q from photon
            ELSEIF(K(21,2).GT.0.AND.K(23,2).EQ.21) THEN
               IF(K(25,2).GT.0.AND.K(26,2).EQ.21) THEN
                  K(NS+1,4) = 26
                  K(NS+1,5) = 0
                  K(NS+3,4)=25
                  K(NS+3,5)=26
                  K(NS+5,4)=(NS+3)*MSTU(5)
                  K(NS+5,5)=0
                  K(NS+6,4)=(NS+1)*MSTU(5)
                  K(NS+6,5)=(NS+3)*MSTU(5)
                  K(NS+7,4)=0
                  K(NS+7,5)=0
                  K(NS+7,1)=0
                  K(NS+7,2) = K(NS+5,2)
c               write(6,*) ' 2nd '
               ELSEIF(K(25,2).EQ.21.AND.K(26,2).GT.0) THEN
                  K(NS+1,4) = 25
                  K(NS+1,5) = 0
                  K(NS+3,4)=26
                  K(NS+3,5)=25
                  K(NS+6,4)=(NS+3)*MSTU(5)
                  K(NS+6,5)=0
                  K(NS+5,4)=(NS+1)*MSTU(5)
                  K(NS+5,5)=(NS+3)*MSTU(5)
                  K(NS+7,4)=0
                  K(NS+7,5)=0
                  K(NS+7,1)=0
                  K(NS+7,2) = K(NS+6,2)
c               write(6,*) ' 4th '
               ENDIF
C antq from proton
            ELSEIF(K(23,2).LT.0.AND.K(21,2).EQ.21) THEN
               IF(K(25,2).LT.0.AND.K(26,2).EQ.21) THEN
                  K(NS+3,4) = 0
                  K(NS+3,5) = 26
                  K(NS+1,4)=26
                  K(NS+1,5)=25
                  K(NS+5,4)=0
                  K(NS+5,5)=(NS+1)*MSTU(5)
                  K(NS+6,4)=(NS+1)*MSTU(5)
                  K(NS+6,5)=(NS+3)*MSTU(5)
                  K(NS+7,4)=0
                  K(NS+7,5)=0
                  K(NS+7,1)=0
                  K(NS+7,2) = K(NS+5,2)
c               write(6,*) ' 1st '
               ELSEIF(K(25,2).EQ.21.AND.K(26,2).LT.0) THEN
                  K(NS+3,4) = 0
                  K(NS+3,5) = 25
                  K(NS+1,4)=25
                  K(NS+1,5)=26
                  K(NS+6,4)=0
                  K(NS+6,5)=(NS+1)*MSTU(5)
                  K(NS+5,4)=(NS+1)*MSTU(5)
                  K(NS+5,5)=(NS+3)*MSTU(5)
                  K(NS+7,4)=0
                  K(NS+7,5)=0
                  K(NS+7,1)=0
                  K(NS+7,2) = K(NS+6,2)
c               write(6,*) ' 3rd '
               ENDIF
c  q from proton
            ELSEIF(K(23,2).GT.0.AND.K(21,2).EQ.21) THEN
               IF(K(25,2).GT.0.AND.K(26,2).EQ.21) THEN
                  K(NS+3,4) = 26
                  K(NS+3,5) = 0
                  K(NS+1,4)=25
                  K(NS+1,5)=26
                  K(NS+5,4)=(NS+1)*MSTU(5)
                  K(NS+5,5)=0
                  K(NS+6,4)=(NS+3)*MSTU(5)
                  K(NS+6,5)=(NS+1)*MSTU(5)
                  K(NS+7,4)=0
                  K(NS+7,5)=0
                  K(NS+7,1)=0
                  K(NS+7,2) = K(NS+5,2)
c               write(6,*) ' 2nd '
               ELSEIF(K(25,2).EQ.21.AND.K(26,2).GT.0) THEN
                  K(NS+3,4) = 25
                  K(NS+3,5) = 0
                  K(NS+1,4)=26
                  K(NS+1,5)=25
                  K(NS+6,4)=(NS+1)*MSTU(5)
                  K(NS+6,5)=0
                  K(NS+5,4)=(NS+3)*MSTU(5)
                  K(NS+5,5)=(NS+1)*MSTU(5)
                  K(NS+7,4)=0
                  K(NS+7,5)=0
                  K(NS+7,1)=0
                  K(NS+7,2) = K(NS+6,2)
c               write(6,*) ' 4th '
               ENDIF
            ENDIF

         ENDIF
ccc         CALL DULIST(2)
      ELSEIF(IRESPRO.EQ.4) THEN
c qqbar --> g g
         IF(ICOLORA.EQ.1) THEN
c color flow A
            IF(K(23,2).GT.0) THEN
               K(NS+1,4) = 0
               K(NS+1,5) = 27
               K(NS+3,4)=27
               K(NS+3,5)=0
               K(NS+5,4)=(NS+6)*MSTU(5)
               K(NS+5,5)=(NS+7)*MSTU(5)
               K(NS+6,4)=(NS+7)*MSTU(5)
               K(NS+6,5)=(NS+5)*MSTU(5)
               K(NS+7,4)=(NS+3)*MSTU(5)+26
               K(NS+7,5)=(NS+1)*MSTU(5)+25
c            write(6,*) ' flow A, quark'
            ELSEIF(K(23,2).LT.0) THEN
               K(NS+1,4) = 27
               K(NS+1,5) = 0
               K(NS+3,4)=0
               K(NS+3,5)=27
               K(NS+5,4)=(NS+7)*MSTU(5)
               K(NS+5,5)=(NS+6)*MSTU(5)
               K(NS+6,4)=(NS+5)*MSTU(5)
               K(NS+6,5)=(NS+7)*MSTU(5)
               K(NS+7,4)=(NS+1)*MSTU(5)+25
               K(NS+7,5)=(NS+3)*MSTU(5)+26
c            write(6,*) ' flow A, antiquark'
            ELSE
               write(6,*) ' wrong flavor',K(23,2),' for IRESPRO=',
     +         IRESPRO
               write(6,*) ' color conf: ',ICOLORA
            ENDIF
         ELSEIF(ICOLORA.EQ.2) THEN
c color flow B
            IF(K(23,2).GT.0) THEN
               K(NS+1,4) = 0
               K(NS+1,5) = 27
               K(NS+3,4)=27
               K(NS+3,5)=0
               K(NS+5,4)=(NS+7)*MSTU(5)
               K(NS+5,5)=(NS+6)*MSTU(5)
               K(NS+6,4)=(NS+5)*MSTU(5)
               K(NS+6,5)=(NS+7)*MSTU(5)
               K(NS+7,4)=(NS+3)*MSTU(5)+25
               K(NS+7,5)=(NS+1)*MSTU(5)+26
c            write(6,*) ' flow B, quark'
            ELSEIF(K(23,2).LT.0) THEN
               K(NS+1,4) = 27
               K(NS+1,5) = 0
               K(NS+3,4)=0
               K(NS+3,5)=27
               K(NS+5,4)=(NS+6)*MSTU(5)
               K(NS+5,5)=(NS+7)*MSTU(5)
               K(NS+6,4)=(NS+7)*MSTU(5)
               K(NS+6,5)=(NS+5)*MSTU(5)
               K(NS+7,4)=(NS+1)*MSTU(5)+26
               K(NS+7,5)=(NS+3)*MSTU(5)+25
c            write(6,*) ' flow B, antiquark'
            ELSE
               write(6,*) ' wrong flavor',K(23,2),' for IRESPRO=',
     +         IRESPRO
               write(6,*) ' color conf: ',ICOLORA
            ENDIF
         ELSE
            write(6,*) ' wrong color configuartion ',ICOLORA
         ENDIF
      ELSEIF(IRESPRO.EQ.5) THEN
c qqbar --> q qbar
         IF(ICOLORA.EQ.1) THEN
c color flow A
            IF(K(23,2).GT.0) THEN
               K(NS+1,4) = 0
               K(NS+1,5) = 23
               K(NS+3,4)=21
               K(NS+3,5)=0
               K(NS+5,4)=0
               K(NS+5,5)=(NS+6)*MSTU(5)
               K(NS+6,4)=(NS+5)*MSTU(5)
               K(NS+6,5)=0
               K(NS+7,4)=0
               K(NS+7,5)=0
               K(NS+7,1)=0
c            write(6,*) ' flow A, quark'
            ELSEIF(K(23,2).LT.0) THEN
               K(NS+1,4) = 23
               K(NS+1,5) = 0
               K(NS+3,4)=0
               K(NS+3,5)=21
               K(NS+5,4)=(NS+6)*MSTU(5)
               K(NS+5,5)=0
               K(NS+6,4)=0
               K(NS+6,5)=(NS+5)*MSTU(5)
               K(NS+7,4)=0
               K(NS+7,5)=0
               K(NS+7,1)=0
c            write(6,*) ' flow A, antiquark'
            ELSE
               write(6,*) ' wrong flavor',K(23,2),' for IRESPRO=',
     +         IRESPRO
               write(6,*) ' color conf: ',ICOLORA
            ENDIF
         ELSEIF(ICOLORA.EQ.2) THEN
c color flow B
            IF(K(23,2).GT.0) THEN
               K(NS+1,4) = 0
               K(NS+1,5) = 27
               K(NS+3,4)=27
               K(NS+3,5)=0
               IF(K(25,2).GT.0) THEN
                  K(NS+5,4)=(NS+7)*MSTU(5)
                  K(NS+5,5)=0
                  K(NS+6,4)=0
                  K(NS+6,5)=(NS+7)*MSTU(5)
                  K(NS+7,4)=(NS+3)*MSTU(5)+25
                  K(NS+7,5)=(NS+1)*MSTU(5)+26
               ELSE
                  K(NS+5,4)=0
                  K(NS+5,5)=(NS+7)*MSTU(5)
                  K(NS+6,4)=(NS+7)*MSTU(5)
                  K(NS+6,5)=0
                  K(NS+7,4)=(NS+3)*MSTU(5)+26
                  K(NS+7,5)=(NS+1)*MSTU(5)+25
               ENDIF
c            write(6,*) ' flow B, quark'
            ELSEIF(K(23,2).LT.0) THEN
               K(NS+1,4) = 27
               K(NS+1,5) = 0
               K(NS+3,4)=0
               K(NS+3,5)=27
               IF(K(25,2).GT.0) THEN
                  K(NS+5,4)=(NS+7)*MSTU(5)
                  K(NS+5,5)=0
                  K(NS+6,4)=0
                  K(NS+6,5)=(NS+7)*MSTU(5)
                  K(NS+7,4)=(NS+1)*MSTU(5)+25
                  K(NS+7,5)=(NS+3)*MSTU(5)+26
               ELSE
                  K(NS+5,4)=0
                  K(NS+5,5)=(NS+7)*MSTU(5)
                  K(NS+6,4)=(NS+7)*MSTU(5)
                  K(NS+6,5)=0
                  K(NS+7,4)=(NS+1)*MSTU(5)+26
                  K(NS+7,5)=(NS+3)*MSTU(5)+25
               ENDIF
c            write(6,*) ' flow B, antiquark'
            ELSE
               write(6,*) ' wrong flavor',K(23,2),' for IRESPRO=',
     +         IRESPRO
               write(6,*) ' color conf: ',ICOLORA
            ENDIF
         ELSE
            write(6,*) ' wrong color configuartion ',ICOLORA
         ENDIF
      ELSEIF(IRESPRO.EQ.6) THEN
c q q --> q q
         IF(ICOLORA.EQ.1) THEN
c color flow A
            IF(K(23,2).GT.0) THEN
               K(NS+1,4) = 26
               K(NS+1,5) = 0
               K(NS+3,4)=25
               K(NS+3,5)=0
               K(NS+5,4)=(NS+3)*MSTU(5)
               K(NS+5,5)=0
               K(NS+6,4)=(NS+1)*MSTU(5)
               K(NS+6,5)=0
               K(NS+7,4)=0
               K(NS+7,5)=0
               K(NS+7,1)=0
c            write(6,*) ' flow A, quark'
            ELSEIF(K(23,2).LT.0) THEN
               K(NS+1,4) = 0
               K(NS+1,5) = 26

               K(NS+3,4)=0
               K(NS+3,5)=25
               K(NS+5,4)=0
               K(NS+5,5)=(NS+3)*MSTU(5)
               K(NS+6,4)=0
               K(NS+6,5)=(NS+1)*MSTU(5)
               K(NS+7,4)=0
               K(NS+7,5)=0
               K(NS+7,1)=0
c            write(6,*) ' flow A, antiquark'
            ELSE
               write(6,*) ' wrong flavor',K(23,2),' for IRESPRO=',
     +         IRESPRO
               write(6,*) ' color conf: ',ICOLORA
            ENDIF
         ELSEIF(ICOLORA.EQ.2) THEN
c color flow B
            IF(K(23,2).GT.0) THEN
               K(NS+1,4) = 25
               K(NS+1,5) = 0

               K(NS+3,4)=26
               K(NS+3,5)=0
               K(NS+5,4)=(NS+1)*MSTU(5)
               K(NS+5,5)=0
               K(NS+6,4)=(NS+3)*MSTU(5)
               K(NS+6,5)=0
               K(NS+7,4)=0
               K(NS+7,5)=0
               K(NS+7,1)=0
c            write(6,*) ' flow B, quark'
            ELSEIF(K(23,2).LT.0) THEN
               K(NS+1,4) = 0
               K(NS+1,5) = 25

               K(NS+3,4)=0
               K(NS+3,5)=26
               K(NS+5,4)=0
               K(NS+5,5)=(NS+1)*MSTU(5)
               K(NS+6,4)=0
               K(NS+6,5)=(NS+3)*MSTU(5)
               K(NS+7,4)=0
               K(NS+7,5)=0
               K(NS+7,1)=0
c            write(6,*) ' flow B, antiquark'
            ELSE
               write(6,*) ' wrong flavor',K(23,2),' for IRESPRO=',
     +         IRESPRO
               write(6,*) ' color conf: ',ICOLORA
            ENDIF
         ELSE
            write(6,*) ' wrong color configuartion ',ICOLORA
         ENDIF
      ELSEIF(IRESPRO.EQ.7) THEN
c q q --> q q color singlet exchange
         K(NS+1,4) = 25
         K(NS+1,5) = 25
         K(NS+3,4)=26
         K(NS+3,5)=26
         K(NS+5,4)=(NS+1)*MSTU(5)
         K(NS+5,5)=(NS+1)*MSTU(5)
         K(NS+6,4)=(NS+3)*MSTU(5)
         K(NS+6,5)=(NS+3)*MSTU(5)
         K(NS+7,4)=0
         K(NS+7,5)=0
         K(NS+7,1)=0
      ELSEIF(IRESPRO.EQ.8) THEN
c qg event
C color configuartion A
c         IF(ICOLORA.EQ.1) THEN
C antq from photon
            IF(K(21,2).LT.0.AND.K(23,2).EQ.21) THEN
                  K(NS+1,4) = 0
                  K(NS+1,5) = 23
                  K(NS+3,4)=21
                  K(NS+3,5)=27
                  K(NS+5,4)=0
                  K(NS+5,5)=(NS+7)*MSTU(5)
                  K(NS+6,4)=0
                  K(NS+6,5)=0
                  K(NS+7,4)=0
                  K(NS+7,5)=(NS+3)*MSTU(5)+25
                  K(NS+7,2) = K(NS+5,2)
c               write(6,*) ' 1st '
c  q from photon
            ELSEIF(K(21,2).GT.0.AND.K(23,2).EQ.21) THEN
                  K(NS+1,4) = 23
                  K(NS+1,5) = 0
                  K(NS+3,4)=27
                  K(NS+3,5)=21
                  K(NS+5,4)=(NS+7)*MSTU(5)
                  K(NS+5,5)=0
                  K(NS+6,4)=0
                  K(NS+6,5)=0
                  K(NS+7,4)=(NS+3)*MSTU(5)+25
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+5,2)
c               write(6,*) ' 2nd '
C antq from proton
            ELSEIF(K(23,2).LT.0.AND.K(21,2).EQ.21) THEN
                  K(NS+3,4) = 0
                  K(NS+3,5) = 21
                  K(NS+1,4)=23
                  K(NS+1,5)=27
                  K(NS+5,4)=0
                  K(NS+5,5)=0
                  K(NS+6,4)=0
                  K(NS+6,5)=(NS+7)*MSTU(5)
                  K(NS+7,4)=(NS+1)*MSTU(5)+26
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+6,2)
c               write(6,*) ' 3rd '
c  q from proton
            ELSEIF(K(23,2).GT.0.AND.K(21,2).EQ.21) THEN
                  K(NS+3,4) = 21
                  K(NS+3,5) = 0
                  K(NS+1,4)=27
                  K(NS+1,5)=23
                  K(NS+5,4)=0
                  K(NS+5,5)=0
                  K(NS+6,4)=(NS+7)*MSTU(5)
                  K(NS+6,5)=0
                  K(NS+7,4)=(NS+1)*MSTU(5)+26
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+6,2)
c               write(6,*) ' 4th '
            ENDIF

c         ENDIF      
      ELSEIF(IRESPRO.EQ.9) THEN
         IF(K(23,2).GT.0) THEN
               K(NS+1,4) = 0
               K(NS+1,5) = 27
               K(NS+3,4)=27
               K(NS+3,5)=0
               K(NS+6,4)=0
               K(NS+6,5)=0
               K(NS+5,4)=(NS+7)*MSTU(5)
               K(NS+5,5)=(NS+7)*MSTU(5)
               K(NS+7,4)=(NS+3)*MSTU(5)+25
               K(NS+7,5)=(NS+1)*MSTU(5)+25
c            write(6,*) ' flow A, quark'
         ELSEIF(K(23,2).LT.0) THEN
               K(NS+1,4) = 27
               K(NS+1,5) = 0
               K(NS+3,4)=0
               K(NS+3,5)=27
               K(NS+6,4)=0
               K(NS+6,5)=0
               K(NS+5,4)=(NS+7)*MSTU(5)
               K(NS+5,5)=(NS+7)*MSTU(5)
               K(NS+7,4)=(NS+1)*MSTU(5)+25
               K(NS+7,5)=(NS+3)*MSTU(5)+25
c            write(6,*) ' flow A, antiquark'
         ELSE
               write(6,*) ' wrong flavor',K(23,2),' for IRESPRO=',
     +         IRESPRO
               write(6,*) ' color conf: ',ICOLORA
         ENDIF
            
      ELSE
         write(6,*) ' colres : irespro not implemented ',irespro
      ENDIF
      RETURN
      END
