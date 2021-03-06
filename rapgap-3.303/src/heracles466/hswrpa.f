C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSWRPA
C---
C   WRITE PARAMETER DEFINITIONS OF THE CURRENT RUN
C   FOR ALL CHANNELS TO UNIT LUNDAT
C---
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     &                IOPLOT,IPRINT,ICUT
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSPARM/ POLARI,LLEPT,LQUA
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSSTRP/ ICODE,ILIB,ILQMOD
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSISGM/ TCUTQ,TCUTQS
C----------------
      WRITE(LUNDAT) XMIN,XMAX,Q2MIN,YMIN,YMAX,WMIN,ICUT
      WRITE(LUNDAT) TCUTQ,TCUTQS
      WRITE(LUNDAT) (LPARIN(I),I=1,12)
      WRITE(LUNDAT) ICODE,ILIB,ILQMOD
      WRITE(LUNDAT) EELE,EPRO,POLARI,LLEPT
      WRITE(LUNDAT) EGMIN
C
      IF(IPRINT.GE.2) THEN
        WRITE(LUNTES,'(//A)') ' *** TEST PRINT HSWRPA'
        WRITE(LUNTES,'(A/A,I3)')
     &      ' *** PARAMETERS OF THE ACTUAL RUN',
     &      ' *** WRITTEN ONTO UNIT', LUNDAT
        WRITE(LUNTES,'(/A)') ' XMIN,XMAX,Q2MIN,YMIN,YMAX,WMIN,ICUT'
        WRITE(LUNTES,'(6(1PE12.3),I4)')
     &               XMIN,XMAX,Q2MIN,YMIN,YMAX,WMIN,ICUT
        WRITE(LUNTES,'(/A)') ' TCUTQ,TCUTQS'
        WRITE(LUNTES,'(2(1PE12.3))') TCUTQ,TCUTQS
        WRITE(LUNTES,'(/A)') ' (LPARIN(I),I=1,12)'
        WRITE(LUNTES,'(12I2)') (LPARIN(I),I=1,12)
        WRITE(LUNTES,'(/A,I3,A,I8,A,I3))')
     &    ' ICODE = ',ICODE,'ILIB = ',ILIB,' ILQMOD = ',ILQMOD
        WRITE(LUNTES,'(/A)') ' EELE,EPRO,POLARI,LLEPT'
        WRITE(LUNTES,'(3(1PE12.3),I5)') EELE,EPRO,POLARI,LLEPT
        WRITE(LUNTES,'(A,1PE12.3)') ' EGMIN=',EGMIN
      ENDIF
      RETURN
      END
