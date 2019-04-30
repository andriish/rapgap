C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSWRSA(ISET,NREGX,NDIMX,
     &                  SIG,SIGE,TGMAX,NDOX,TLMAX,NM,XX,
     &                  FFGOX,DNCGX,FFLOX,DNCLX,GOLDX,
     &                  NTOTX,NCALX,NCA1X,NCA2X,IBIMX,JCORX,
     &                  LGLOX,LLOCX,
     &                  IT,SI,SI2,SWGT,SCHI)
C---
C   WRITE INFORMATION FOR SAMPLING FOR ONE CONTRIBUTION TO UNIT LUNDAT
C   ISET.GT.0 : RADIATIVE CHANNELS (INTEGRATION INFORMATION FROM VEGAS)
C---
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      real*8 ntotx
      LOGICAL LGLOX,LLOCX
C
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     &                IOPLOT,IPRINT,ICUT
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      DIMENSION TLMAX(NREGX),NM(NREGX),XX(NDOX,NDIMX)
C---------------------------------
      WRITE(LUNDAT) SIG,SIGE,TGMAX,NDOX
      IF(IPRINT.GE.6) THEN
        WRITE(LUNTES,'(/A,3I5)') ' HSWRSA:  ISET, NREGX, NDIMX',
     &                                      ISET, NREGX, NDIMX
        WRITE(LUNTES,'(/A)') ' HSWRSA:  SIG,SIGE,TGMAX,NDOX'
        WRITE(LUNTES,*) SIG,SIGE,TGMAX,NDOX
      ENDIF
      WRITE(LUNDAT) (TLMAX(I),I=1,NREGX)
      WRITE(LUNDAT) (NM(I),I=1,NREGX)
C
      IF(ISET.GT.0) WRITE(LUNDAT) IT,SI,SI2,SWGT,SCHI
      WRITE(LUNDAT) ((XX(I,II),I=1,NDOX),II=1,NDIMX)
      WRITE(LUNDAT) FFGOX,DNCGX,FFLOX,DNCLX,GOLDX,
     &              NTOTX,NCALX,NCA1X,NCA2X,IBIMX,JCORX,
     &              LGLOX,LLOCX
      IF(IPRINT.GE.6) THEN
        WRITE(LUNTES,'(/A/5X,5(1PD13.5))')
     &       ' HSWRSA :    FFGOX,DNCGX,FFLOX,DNCLX,GOLDX',
     &                    FFGOX,DNCGX,FFLOX,DNCLX,GOLDX
        WRITE(LUNTES,'(/A/5X,2I8,4I5)')
     &       ' HSWRSA :    NTOTX,NCALX,NCA1X,NCA2X,IBIMX,JCORX',
     &                    NTOTX,NCALX,NCA1X,NCA2X,IBIMX,JCORX
        WRITE(LUNTES,'(/A,2L4)')
     &       ' HSWRSA:     LGLOX,LLOCX',  LGLOX,LLOCX
      ENDIF
      IF(IPRINT.GE.6) THEN
         WRITE(LUNTES,'(/A)') ' HSWRSA: XX(NDOX,NDIMX)'
         WRITE(LUNTES,2) ((XX(I,J),I=1,NDOX),J=1,NDIMX)
         WRITE(LUNTES,'(/A)') ' HSWRSA: TLMAX(NREGX)'
         WRITE(LUNTES,2) (TLMAX(I),I=1,NREGX)
         WRITE(LUNTES,'(/,A)') ' HSWRSA: NM2(NDOX,NDIMX)'
         WRITE(LUNTES,3) (NM(I),I=1,NREGX)
      ENDIF
      RETURN
 2    FORMAT(5(1PD15.5))
 3    FORMAT(10I8)
      END
