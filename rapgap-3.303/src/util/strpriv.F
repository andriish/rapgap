      SUBROUTINE STRPRIV(X,MUF,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL)
C   This routine provides the interface to the pdf-source code 
C   as provided by the various groups. It returns the pdf-values DXPDF
C   for the different quark flavors at a given  phase space point
C   X, MUF. The quark flavors are ordered according to the PDG convention
C   i.e DXPDF(0) = gluon 
C       DXPDF(1) = down valence quark + \bar{down}_{sea} 
C       DXPDF(2) = up valence quark + \bar{up}_{sea} 
C       DXPDF(3,4,5,6) = strange, charm, bottom, top quark
C       DXPDF(-1,..,-6) are the corresponding antiquarks
C  
C       INPUT: IGROUP = Group Number of which provides pdf 
C                       according to PDFLIB convention
C                       e.g. CTEQ-Group = 4, etc       
C              ISET   = Set porvided by group IGROUP 
C                       Attention: here we do not follow
C                        (for the time being) the 'official'
C                       convention, but ISET must be composed
C                       by a package number and a value
C                       indication the chosen parameterisation (mode) 
C                       within that package
C                       i.e. IGROUP = 4 ISET = 51 would result
C                       to the first parameterization within
C                       the CTEQ5 set. 
C                       For MRST we take the year of the
C                       publication together with the
C                       order of calculation (0=LO 1=NLO) and
C                       chosen parameterization (mode) within the 
C                       set.
C                       i.e. IGROUP = 3 ISET = 2001101
C                       would result in the first NLO
C                       parameterization within the
C                       2001 set.
C                       To obtain information 
C                       on the different modes please consult
C                       http://cpt19.dur.ac.uk/hepdata/pdf.html 
C                       and references therein.
C                   
C                  X  = x-value of parton   
C                MUF  = QCD scale in GeV
C   
C
C      OUTPUT: DXPDF (see above)       
C Author R.P.  
C ******************************************************************   
      IMPLICIT NONE
      INTEGER I
      INTEGER IGROUP, ISET,idif
	Common/strnew/IGROUP,ISET,idif
	INTEGER ISERI, IMODE, ICOMP
      DOUBLE PRECISION MUF, X 
      DOUBLE PRECISION UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL 
      DOUBLE PRECISION UPV1,DNV1,USEA1,DSEA1,STR1,CHM1,BOT1,TOP1,GLU1
      DOUBLE PRECISION SEAU1,SEAD1
      DOUBLE PRECISION CHCKSM
      DOUBLE PRECISION CTQ4PDFm, CTQ5PDFm, CTQ6PDFm
      REAL XIN, Q2IN
      REAL GLUON, DVAL, UVAL, DBAR, UBAR, SBAR, CBAR
      INTEGER     IFIRST
      DATA        IFIRST/0/
      SAVE        IFIRST
	Logical first
	data first/.true./
	Double Precision QCDL4,QCDL5
      COMMON/W50512/ QCDL4,QCDL5
      Double Precision qcd4ctq4(13)
	Data qcd4ctq4/.298,.298,.236,.215,.254,.346,.401,.303,.261,
     &	        .298,.298,.385,.292/
      Double Precision qcd5ctq4(13)
	Data qcd5ctq4/.202,.202,.181,.140,.169,.239,.282,.206,.174,
     &	        .202,.202,.385,.292/
      Double Precision qcd4ctq5(9)
	Data qcd4ctq5/.326,.326,.192,.326,.326,.395,.309,.326,.326/
      Double Precision qcd5ctq5(9)
	Data qcd5ctq5/.226,.226,.146,.226,.226,.395,.309,.226,.226/
      Double Precision qcd4ctq6(3)
	Data qcd4ctq6/.326,.326,.326/
      Double Precision qcd5ctq6(3)
	Data qcd5ctq6/.226,.226,.226/
      Character*30 ehkqsnam(3)
	data ehkqsnam/'p-set1table.dat','p-set2atable.dat',
     &	'p-set2btable.dat'/
      Double precision q2
c     igroup = 3:
c            = 4: CTEQ
c            = 5: ehkqs, set 1,2,3
      idif = 0
c      IGROUP = 4
c      ISET   =  61
c         ISET   =  51
c      IGROUP = 5
c      ISET   =  2

      If(first) then 
	  if(IGROUP.EQ.4) then
	    write(6,*) ' CTEQ pdfs. Set = ',Iset
	    ELSEIF(IGROUP.EQ.5) then
	    write(6,*) ' ehkqs pdfs. Set = ',Iset
	    ELSE
	    write(6,*) 'set not available for Igroup,Iset ',Igroup,Iset
	    Stop
	    Endif
	    first=.false.
	Endif
CRP***Return if required PDF-Set is not available      
      GL = 0
      DNV = 0 
      UPV  = 0         
	DSEA = 0
      USEA = 0
      STR = 0
      CHM = 0
      BOT = 0
      TOP = 0
	 
	if(idif.eq.1) then 
	   QCDL4 = qcd4ctq6(1)
C. this is the lambda (NFL=5) used in the calculation
         QCDL5 = qcd5ctq6(1)
c         call diffpdf(X,MUF,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL)
      else  

	   MUF=max(muf,1d0)
         IF(IGROUP.NE.3.AND.IGROUP.NE.4.and.IGROUP.NE.5) THEN
            PRINT*,'Required set ',IGROUP,' not available'
            RETURN
         ENDIF
         IF(IGROUP.EQ.4) THEN 
CRP***ISER is number of published series of group IGROUP 
CRP***i.e. ISER = 5 for CTEQ5                
            ISERI = INT(FLOAT(ISET)/10.)    
CRP***Return if this series is not available
            IF(ISERI.LT.4.OR.ISERI.GT.6) THEN
               PRINT*,'Required series ',ISERI,' for group ',IGROUP,'' 
               PRINT*,'not available'
               RETURN
            ELSE
               IMODE = ISET - ISERI*10
            ENDIF
CRP***Obtain values        
            IF(IFIRST.EQ.0) THEN
               PRINT*,'IGROUP, ISERI, MODE=', IGROUP, ISERI,IMODE 
               IFIRST = 1
            ENDIF
            IF(ISERI.EQ.4) then
		   CALL SetCtq4(IMODE)
		   QCDL4 = qcd4ctq4(imode)
C. this is the lambda (NFL=5) used in the calculation
               QCDL5 = qcd5ctq4(imode)
            Elseif(ISERI.EQ.5) then
		   CALL SetCtq5(IMODE)
		   QCDL4 = qcd4ctq5(imode)
C. this is the lambda (NFL=5) used in the calculation
               QCDL5 = qcd5ctq5(imode)
            Elseif(ISERI.EQ.6) then
		   CALL SetCtq6(IMODE)  
		   QCDL4 = qcd4ctq6(imode)
C. this is the lambda (NFL=5) used in the calculation
               QCDL5 = qcd5ctq6(imode)
	      Endif    

C
C     generate the individual structure fcn calls
C
            IF(ISERI.EQ.4) THEN 
	         UPV1 = Ctq4Pdfm( 1,X,MUF)
               DNV1 = Ctq4Pdfm( 2,X,MUF)
               SEAU1 = Ctq4Pdfm(-1,X,MUF)
               USEA = SEAU1*X
               UPV = (UPV1-SEAU1)*X
               SEAD1 = Ctq4Pdfm(-2,X,MUF)
               DSEA = SEAD1*X
               DNV = (DNV1-SEAD1)*X
               STR1 = Ctq4Pdfm(-3,X,MUF)
               STR = STR1*X
               CHM1 = Ctq4Pdfm(-4,X,MUF)
               CHM = CHM1*X
               BOT1 = Ctq4Pdfm(-5,X,MUF)
               BOT = BOT1*X
               TOP1 = 0.
C              IF (DSCAL.GT.TMAS) TOP = Ctq5Pdf(ISET,-6,X,MUF)
               TOP = TOP1*X
               GLU1 = Ctq4Pdfm( 0,X,MUF)
               GL = GLU1*X
            ELSEIF(ISERI.EQ.5) THEN 
	         UPV1 = Ctq5Pdfm( 1,X,MUF)
               DNV1 = Ctq5Pdfm( 2,X,MUF)
               SEAU1 = Ctq5Pdfm(-1,X,MUF)
               USEA = SEAU1*X
               UPV = (UPV1-SEAU1)*X
               SEAD1 = Ctq5Pdfm(-2,X,MUF)
               DSEA = SEAD1*X
               DNV = (DNV1-SEAD1)*X
               STR1 = Ctq5Pdfm(-3,X,MUF)
               STR = STR1*X
               CHM1 = Ctq5Pdfm(-4,X,MUF)
               CHM = CHM1*X
               BOT1 = Ctq5Pdfm(-5,X,MUF)
               BOT = BOT1*X
               TOP1 = 0.
C              IF (DSCAL.GT.TMAS) TOP = Ctq5Pdfm(ISET,-6,X,MUF)
               TOP = TOP1*X
               GLU1 = Ctq5Pdfm( 0,X,MUF)
               GL = GLU1*X
            ELSEIF(ISERI.EQ.6) THEN 
	         UPV1 = Ctq6Pdfm( 1,X,MUF)
               DNV1 = Ctq6Pdfm( 2,X,MUF)
               SEAU1 = Ctq6Pdfm(-1,X,MUF)
               USEA = SEAU1*X
               UPV = (UPV1-SEAU1)*X
               SEAD1 = Ctq6Pdfm(-2,X,MUF)
               DSEA = SEAD1*X
               DNV = (DNV1-SEAD1)*X
               STR1 = Ctq6Pdfm(-3,X,MUF)
               STR = STR1*X
               CHM1 = Ctq6Pdfm(-4,X,MUF)
               CHM = CHM1*X
               BOT1 = Ctq6Pdfm(-5,X,MUF)
               BOT = BOT1*X
               TOP1 = 0.
C              IF (DSCAL.GT.TMAS) TOP = Ctq5Pdfm(ISET,-6,X,MUF)
               TOP = TOP1*X
               GLU1 = Ctq6Pdfm( 0,X,MUF)
               GL = GLU1*X
            ENDIF 
         ELSEIF(IGROUP.EQ.5) then
	      
	      if(ifirst.eq.0) then
            call setintPDF(ehkqsnam(iset))
		ifirst=1
		endif
		q2=muf**2
            call getintPDF(x,q2,upv,dnv,usea,dsea,str,chm,bot,top,gl)
		QCDL4 = qcd4ctq6(imode)
C. this is the lambda (NFL=5) used in the calculation
            QCDL5 = qcd5ctq6(1)
	   
	   ENDIF 
      ENDIF 
C
c      write(6,*) ' GL,GLU1,X,MUF ',GL,GLU1,X,MUF
CRP***
      RETURN
      END
