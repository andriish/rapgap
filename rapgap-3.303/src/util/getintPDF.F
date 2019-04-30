C*************************************************************************
C
C  getintPDF  subroutine library.
C
C  By:    Vesa Kolhinen <vesa.kolhinen@phys.jyu.fi>
C  Date:  09.10.2003
C
C
C    Subroutines for reading and interpolating PDF's from given 
C    file. Gives interpolated parton distributions at requested 
C    x and Q^2 points. Uses linear/logarithmic interpolation for 
C    x and logarithmic one for Q^2. Log/lin. treshold for x is
C    currently hardcoded to 0.6
C    
C    Using x and Q^2 values outside the limits returns distributions 
C    at the corresponding limit: e.g. xg(x<x_0,Q^2) returns the value
C    xg(x=x_0,Q^2), etc.
C
C
C  Usage:  
C
C    Read and set up the interpolation tables. Interpolate.
C
C    call setintPDF(filename)
C    call getintPDF(xval,Q2val,uv,dv,us,ds,st,ch,bot,top,gl)
C
C    where the variables are:
C
C    character*30 filename
C    double precision xval,Q2val,uv,dv,us,ds,st,ch,bot,top,gl
C
C
C    setintPDF initializes the interpolation routine and tables.
C    "filename" is the data file containing the PDF's used.
C    This needs to be used only once, but can be repeated if
C    the used set is changed during the run.
C
C    getintPDF returns xf_i(x,Q^2) for given x, Q^2 and parton
C    flavour i. Format is the same as that of PDFLIB's structm.
C
C    These routines use common areas named as "getPDFi", "getPDFd" 
C    and "getPDFt". Do not mess these.
C
C
C  Comments about the code:
C    Data is read to 3-dimensional table "dist".
C    First index defines column of the file, i.e. the type
C    of the distribution (x,xgl,xu_val,...)
C    Second index means x and third Q2 values.
C    
C    About interpolation:
C    First, interpolated x values are sought in two x-loops.
C    When right x_min and x_max are found, two similar Q^2 loops are 
C    executed. Interpolation is then done for all requested Q^2
C    values at this x point.
C
C    Q^2 interpolation is fully logarithmic, x interpolation is 
C    logarithmic up to x = xtres and linear above that.
C    xtres is currently hardcoded to be 0.6
C
C
C    About variables: 
C    Index variables are named as ix, nox, maxx  where ix is 
C    running index, nox is actual size of the table and maxx is 
C    the maximum size of the table, and similarly for all the 
C    variables. Indices which end with "letter i" are used for 
C    interpolated values.
C
C
C
C*************************************************************************



C*************************************************************************
C  setintPDF   
C    Subroutine for making necessary setups for getintPDF routine.
C    setintPDF reads given PDF data file and sets up the common area
C    table for interpolation.
C 
C  Usage: setintPDF(filename)
C    where filename is the PDF data file.
C
C*************************************************************************
      subroutine setintPDF(filename)
      implicit none

      integer ix,nox,maxx,iQ2,noQ2,maxQ2,icol,ncol,maxcol
      parameter (maxx=185,maxQ2=750,maxcol=10)
      integer nx,nQ2, ixstart, ixend, iQ2start, iQ2end

      double precision xmin,xmax,Q2min,Q2max
      double precision dist(maxcol,maxx,maxQ2) 
      double precision xtres
      double precision gettagvalue

      character*30 filename

      common/getPDFd/xtres,xmin,xmax,Q2min,Q2max
      common/getPDFt/dist
      common/getPDFi/ncol,nox,noQ2,nx,nQ2,
     +     ixstart,ixend,iQ2start,iQ2end

 
C Threshold for logarithmic/linear x: currently hardcoded to 0.6
      xtres = 0.6


      call readfile(filename,dist,maxcol,maxx,maxQ2,ncol,nox,noQ2)

C Check in which directions tabulated values increase
      ixstart = 1
      ixend   = nox-1
      nx      = 1
      xmin = dist(1,1,1)
      xmax = dist(1,nox,1)

      iQ2start= 1
      iQ2end  = noQ2-1
      nQ2     = 1
      Q2min = dist(2,1,1)
      Q2max = dist(2,1,noQ2)
      if (dist(1,1+1,1) .le. dist(1,1,1) ) then
         ixstart= 2
         ixend  = nox
         nx     = -1
         xmin = dist(1,nox,1)
         xmax = dist(1,1,1)
      endif
      if ( dist(2,1,1+1) .le. dist(2,1,1) ) then
         iQ2start= 2
         iQ2end  = noQ2
         nQ2     = -1
         Q2min = dist(2,1,noQ2)
         Q2max = dist(2,1,1)
      endif


      return 

      end

C*************************************************************************
C  getintPDF
C    Subroutine returns PDF's xf_i(x,Q^2) at given x,Q^2 point. 
C 
C    Before calling this subroutine, the interpolation routine
C    must be initialized by calling setintPDF with the appropriate
C    data file name.
C
C
C*************************************************************************
      subroutine getintPDF(x,Q2,uv,dv,us,ds,st,ch,bot,top,gl)
      implicit none

      double precision x,Q2,uv,dv,us,ds,st,ch,bot,top,gl

      integer ix,nox,maxx,iQ2,noQ2,maxQ2,icol,ncol,maxcol
      parameter (maxx=185,maxQ2=750,maxcol=10)
      integer nx,nQ2, ixstart, ixend, iQ2start, iQ2end

      double precision xval,Q2val
      double precision dist(maxcol,maxx,maxQ2) 
      double precision intdst(maxcol)
      double precision xtres
      double precision xmin,xmax,Q2min,Q2max
      double precision gettagvalue

      double precision interpolate

      external interpolate

      common/getPDFd/xtres,xmin,xmax,Q2min,Q2max
      common/getPDFt/dist
      common/getPDFi/ncol,nox,noQ2,nx,nQ2,
     +     ixstart,ixend,iQ2start,iQ2end


C Primitive check for initialization
      if (ixstart .eq. 0) then
         write(*,*) 'ERROR: getintPDF is not initialized.'
         stop
      endif

C Check bounds
      xval = x
      Q2val = Q2
      if (x .lt. xmin) then
         xval = xmin
      else if (x .gt. xmax) then
         xval = xmax
      endif
      if (Q2 .lt. Q2min) then
         Q2val = Q2min
      else if (Q2 .gt. Q2max) then
         Q2val = Q2max
      endif


C Search loop and actual interpolation
      do ix = ixstart,ixend
         if ( (dist(1,ix,1).le.xval) .and.  
     +        (xval .le. dist(1,ix+nx,1))) then
            do iQ2 = iQ2start,iQ2end
               if ( (dist(2,ix,iQ2).le.Q2val) .and.  
     +              (Q2val .le. dist(2,ix,iQ2+nQ2))) then

                  do icol = 3,10
                     intdst(icol) = 
     +                    interpolate(xval,Q2val,
     +                    dist(1,ix,iQ2),dist(1,ix+nx,iQ2),
     +                    dist(2,ix,iQ2),dist(2,ix,iQ2+nQ2),
     +                    dist(icol,ix,iQ2),
     +                    dist(icol,ix+nx,iQ2),
     +                    dist(icol,ix,iQ2+nQ2),
     +                    dist(icol,ix+nx,iQ2+nQ2),
     +                    xtres)
                  enddo

                  gl = intdst(3)
                  uv = intdst(4)
                  dv = intdst(5)
                  st = intdst(6) 
                  us = intdst(7) 
                  ds = intdst(8) 
                  ch = intdst(9)
                  bot= intdst(10)
                                    
                  return

               endif
            enddo
         endif
      enddo


      end



C*************************************************************************
C
C  readfile  Reads a table of values from a given file
C
C   fname:     file to be read  (input)
C   dist:      table to be filled  (output)
C   maxcol, maxx, maxQ2:   size of table  (input)
C   ncol, nox, noQ2: actual size of read table (output)
C
C*************************************************************************
      subroutine readfile(fname,dist,maxcol,maxx,maxQ2,ncol,nox,noQ2)
      implicit none

      integer ix,nox,maxx,iQ2,noQ2,maxQ2,icol,ncol,maxcol
      double precision dist(maxcol,maxx,maxQ2) 
      double precision Q2
      double precision gettagvalue

      integer nin
      character*30 fname
      character*200 line
      logical iscomment,istagline,isdata


      nin = 91
      open(unit=nin, file=fname,status="old")

      ncol = 10
      ix = 1
      iQ2 = 0
 100  continue

      read(nin,'(a)',end=200,err=199) line

      if (istagline('Q2=',line).or.istagline('! Q2=',line)) then
         Q2 = gettagvalue('Q2=',line)
         if (Q2.eq.0.0d0) then
            Q2 = gettagvalue('! Q2=',line)
         endif
         iQ2 = iQ2+1
         nox = ix
         ix = 1
         if (iQ2 .gt. maxQ2) then
            write(*,*) 'intpoldst: ERROR: too many Q2 values'
            stop
         endif
      else if ( iscomment(line) ) then
C Print comment lines from the data file if necessary
C         write(*,'(a75)') line
      else 
         if (iQ2 .eq. 0) iQ2 = 1
         if (isdata(10,line)) then
            read(unit=line,fmt='(10(e16.6))')
     +           dist(1,ix,iQ2),dist(2,ix,iQ2),
     +           dist(3,ix,iQ2),dist(4,ix,iQ2),
     +           dist(5,ix,iQ2),dist(6,ix,iQ2),
     +           dist(7,ix,iQ2),dist(8,ix,iQ2),
     +           dist(9,ix,iQ2),dist(10,ix,iQ2)
         else if (isdata(9,line)) then
            read(unit=line,fmt='(9(e16.6))')
     +           dist(1,ix,iQ2),dist(3,ix,iQ2),dist(4,ix,iQ2),
     +           dist(5,ix,iQ2),dist(6,ix,iQ2),dist(7,ix,iQ2),
     +           dist(8,ix,iQ2),dist(9,ix,iQ2),dist(10,ix,iQ2)
         else
            write(*,*) 'intpoldst: invalid data file format'
         endif
         dist(2,ix,iQ2) = Q2
         ix = ix+1
         if (ix .gt. maxx) then
            write(*,*) 'intpoldst: ERROR: too many x values'
            stop
         endif
      endif
      

      goto 100
      
 199  continue
      write(*,*) 'intpoldst: read error'

 200  continue
      close(nin)

      noQ2 = iQ2
      nox = ix-1


      return
      end


C*************************************************************************
C
C Function ISCOMMENT tests whether given line is a comment, i.e. if it 
C begins with comment character !
C
C Returns: .true., .false.
C
C*************************************************************************
      logical function iscomment(line)
      implicit none

      character*(*) line
      integer i, maxlen

      i=1         
      maxlen=len(line)

 100  continue
      if (line(i:i).eq.' ') then
         i=i+1
         goto 100
      endif

      if (line(i:i).eq.'!' .or. line(i:i).eq.'#' ) then
         iscomment = .true.
      else 
         iscomment = .false.
      endif

      return
      end

C*************************************************************************
C
C  istagline  Tests whether given line is a tag line.
C             Tag line is of type " <tag>=<value>, e.q. Q2= 2.25
C
C*************************************************************************
      logical function istagline(tag,line)
      implicit none

      character*(*) tag
      character*(*) line
      integer i, maxlen, taglen

      i=1         
      maxlen=len(line)
      taglen=len(tag)

 100  continue
      if (line(i:i).eq.' ') then
         i=i+1
         goto 100
      endif

      if (line(i:i+taglen).eq.tag) then
         istagline = .true.
      else 
         istagline = .false.
      endif

      return
      end

C*************************************************************************
C
C  gettagvalue   Reads and returns the value of the tag.
C                Tag line is of type " <tag>=<value>, e.q. Q2= 2.25
C
C*************************************************************************
      double precision function gettagvalue(tag,line)
      implicit none

      character*(*) tag
      character*(*) line
      integer i, maxlen, taglen

      i=1         
      maxlen=len(line)
      taglen=len(tag)

 100  continue
      if (line(i:i).eq.' ') then
         i=i+1
         goto 100
      endif

      if (line(i:i+taglen).eq.tag) then
         i = i+taglen
 110     continue
         if (line(i:i).eq.' ') then
            i=i+1
            goto 110
         endif
         read(line(i:maxlen),*) gettagvalue
      else 
         gettagvalue = 0.0d0
      endif

      return
      end

C*************************************************************************
C
C Function ISDATA tests whether given line contains no elements of data
C Maximum no is 15. Change the dimension of datain-table if necessary.          
C
C Returns: .true., .false.
C
C*************************************************************************
      logical function isdata(no,line)
      integer no,i
      character*1000 line
      double precision datain(15)

      read(line,*,err=900) (datain(i),i=1,no)
      isdata = .true.
      return
 900  isdata = .false.
 
      return
      end



C*************************************************************************
C
C  interpolate
C
C
C   Function interpolates a two variable function.
C   x-scale is logarithmic up to xtres and linear from that on,
C   Q2 scale is completely logarithmic
C
C   x1 and x2 are the lower and upper x values
C   Q21 and Q22 are the Q2 values.
C   fmn are the values of the function at the given x and Q2 points:
C
C           f12       f22
C             +-------+
C       Q2    |       |
C             +-------+
C           f11       f21
C                 x
C*************************************************************************

      double precision function interpolate 
     +     (x,Q2,x1,x2,Q21,Q22,f11,f21,f12,f22,xtres)
      implicit none

      double precision x,Q2
      double precision xtres
      double precision x1,x2,Q21,Q22
      double precision f11,f21,f12,f22,fQ21,fQ22

      logical equal
      external equal

C 1) interpolate between x1 < x < x2 at Q21:
      if ( equal(x,x1) ) then
         fQ21 = f11 
      else if ( equal(x,x2) ) then
         fQ21 = f21
      else
        
         if (x.lt.xtres) then
            fQ21 = f11 + (f21-f11)*log(x/x1)/log(x2/x1)
         else
            fQ21 = f11 + (f21-f11)*(x-x1)/(x2-x1)
         end if
         
         if (((fQ21.lt.f11).and.(fQ21.lt.f21)).or.
     1        ((fQ21.gt.f11).and.(fQ21.gt.f21))) then
            write (*,*) 'intpoldst: error in x-log interpolation'
            write (*,*) f11,fQ21,f21
            stop
         end if

      endif

c 2) interpolate between x1 < x < x2 at Q22:
      if ( equal(x,x1) ) then
         fQ22 = f12 
      else if ( equal(x,x2) ) then
         fQ22 = f22
      else
         
         if (x.lt.xtres) then
            fQ22 = f12 + (f22-f12)*log(x/x1)/log(x2/x1)
         else
            fQ22 = f12 + (f22-f12)*(x-x1)/(x2-x1)
         end if
         
         if (((fQ22.lt.f12).and.(fQ22.lt.f22)).or.
     1        ((fQ22.gt.f12).and.(fQ22.gt.f22))) then
            write (*,*) 'intpoldst: error in x-lin interpolation'
            write (*,*) f12,fQ22,f22
            stop
         end if

      endif

c 3) interpolate in logQ2 between Q21 < Q2 < Q22
      if ( equal(Q2,Q21) ) then
         interpolate = fQ21
      else if ( equal(Q2,Q22) ) then
         interpolate = fQ22
      else
         
         interpolate = fQ21 + (fQ22-fQ21)*log(Q2/Q21)/log(Q22/Q21)
         
         if (((interpolate.gt.fQ21).and.(interpolate.gt.fQ22)).or.
     1        ((interpolate.lt.fQ21).and.(interpolate.lt.fQ22))) then
            write (*,*) 'intpoldst: error in Q2-log interpolation'
            write (*,*) fQ21,interpolate,fQ22
            stop
         end if
         
      endif
      
      return
      end
      
C*************************************************************************
C
C  equal   Function checks if two values are equal within an
C          accuracy of 10^{-8}.
C
C
C*************************************************************************
      logical function equal(a,b)
      implicit none
      double precision a,b
      double precision abs

      abs = 1.0E-8

      equal = .false.
      if ( ((a-abs) .le. b) .and. (b .le. (a+abs)) ) 
     +     equal = .true.

      return
      end
      
