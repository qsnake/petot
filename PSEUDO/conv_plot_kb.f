      program convertKB                                           
C
C     *********************************************************************** 
C     *                                                                     *  
C     * convertkb.f transforms the outputfile 'kbplot.dat' of the program  *  
C     * 'kbconv.f' into a form that can be used to create the corresponding *
C     * plots with 'gnuplot'                                                *
C     *                                                                     *
C     * written by Peter Schuster, April 1993                               *
C     *                                                                     *  
C     *********************************************************************** 
C
      parameter (nmax=5000,numw=20)
C      
      implicit double precision (a-h,o-z)
C
      dimension rx(nmax,numw), ry(nmax,numw),
     1 qx(nmax,numw), qy(nmax,numw),tx(nmax,numw), ty(nmax,numw), 
     2 numbr(numw), numbq(numw), numbt(numw)
C
      character* 3 howdat 
      character*12 labelr(numw),labelq(numw),labelt(numw)
C
      call input (nmax,numw,
     1 rx,ry,qx,qy,tx,ty,numbr,numbq,numbt,
     2 ir,iq,it,labelr,labelq,labelt)
C
C
      irsave=ir
      iqsave=iq
      itsave=it
C
      call outdat (nmax,numw,
     1 rx,ry,qx,qy,tx,ty,numbr,numbq,numbt,
     2 ir,iq,it,labelr,labelq,labelt,rxmax,txmax)
C
      ir=irsave
      iq=iqsave
      it=itsave
C
C     ask if plots or postscript-files are wanted
C
   5  write(6,1000)
1000  format(/,' Do you want <Gnuplot> to display the plots on the scree
     1n',/,' or to create postscript- or pbm-files out of them?',//,' Ty
     2pe <dis> for displaying or <ps> for creating postscript-files',/,
     3' or <pbm> for creating pbm-files.')
      read(5,1001)howdat
1001  format(a3)
C
C
      if(howdat.eq.'dis'.or.howdat.eq.'Dis'.or.howdat.eq.'DIS')then               
      call outdis(numw,
     1 ir,iq,it,labelr,labelq,labelt,rxmax,txmax)
C      
      elseif(howdat.eq.'ps'.or.howdat.eq.' ps'.or.howdat.eq.'PS '.or.ho
     1wdat.eq.' PS')then               
      call outps(numw,
     1 ir,iq,it,labelr,labelq,labelt,rxmax,txmax)
C
      elseif(howdat.eq.'pbm'.or.howdat.eq.'Pbm'.or.howdat.eq.'PBM')then
      call outpbm(numw,
     1 ir,iq,it,labelr,labelq,labelt,rxmax,txmax)
C
      else
      write(6,1002)
1002  format('You made a mistake. Please chose again.')
      goto 5
      endif
C
C            
      stop
      end
C
C***   ***   ***   ***   ***   ***   ***   ***   ***   ***   ***   *** *   
C
      subroutine input (nmax,numw,
     1 rx,ry,qx,qy,tx,ty,numbr,numbq,numbt,
     2 ir,iq,it,labelr,labelq,labelt)
C
C     reads the input from <kbplot.dat>
C 
      parameter (mmax=5000)
C
      implicit double precision (a-h,o-z)
C
      dimension absc(mmax), ord(mmax), rx(nmax,numw), ry(nmax,numw),
     1 qx(nmax,numw), qy(nmax,numw),tx(nmax,numw), ty(nmax,numw), 
     2 numbr(numw), numbq(numw), numbt(numw)
C
      character*3 marker
      character*10 fname
      character*12 labelr(numw),labelq(numw),labelt(numw)
C
C
C
      ir=0
      iq=0
      it=0
C
C
C     asking for outputfile to be used
C            
      write(6,1001)
1001  format(/,' Which file shall be read from?',/,' Please give the nam 
     1e. (format = a10)',/,' Press <ENTER> to choose the default <kbplot
     2.dat>') 
      read(5,1002)fname
1002  format(a10)
      if(fname.eq.' ')fname='kbplot.dat'
C
C
C     opening file to read from
C     
      open(unit=3,file=fname,status='old',form='formatted')
C
C       
C     Read data columns
C 
   1  number=1
C   
C  5  read(3,1003,err=200,end=300)absc(number),ord(number)
C1003 format(1x,f7.4,3x,f10.6)
C
C     (it is not yet possible to read data formatted, because there
C      is not a unique format of the data in the file, unit3)
C
C       
   5  read(3,*,err=200,end=300)absc(number),ord(number)
         number=number+1          
         if(number.ge.nmax)then
            write(4,1004)
            write(6,1004)
1004  format(/,' NOT ALL NUMBERS READ! SET NEW DIMENSION FOR NMAX!',/)
            stop
         endif
      goto 5
C
C   
C     Backup one record and get marker
C
 200  backspace (unit=3)
      read(3,1005)marker
1005  format(8x,a3)
C
C
C
C     Potential markers 'p*r' for RQ = 'r'
C
      if(marker.eq.'psr')then
         ir=ir+1
         labelr(ir)='s - Kl & By '
         numbr(ir)=number-1
         call field1(nmax,numw,numbr,absc,ord,rx,ry,ir)
C
      elseif(marker.eq.'ppr')then
         ir=ir+1
         labelr(ir)='p - Kl & By '
         numbr(ir)=number-1
         call field1(nmax,numw,numbr,absc,ord,rx,ry,ir)
C
      elseif(marker.eq.'pdr')then
         ir=ir+1
         labelr(ir)='d - Kl & By '
         numbr(ir)=number-1
         call field1(nmax,numw,numbr,absc,ord,rx,ry,ir)
C
      elseif(marker.eq.'pfr')then
         ir=ir+1
         labelr(ir)='f - Kl & By '
         numbr(ir)=number-1
         call field1(nmax,numw,numbr,absc,ord,rx,ry,ir)
C
      elseif(marker.eq.'pgr')then
         ir=ir+1
         labelr(ir)='g - Kl & By '
         numbr(ir)=number-1
         call field1(nmax,numw,numbr,absc,ord,rx,ry,ir)
C
C         
C     Fourier markers 'p*q' for RQ = 'q'
C
      elseif(marker.eq.'psq')then
         iq=iq+1
         labelq(iq)='trans s K&B '
         numbq(iq)=number-1
         call field1(nmax,numw,numbq,absc,ord,qx,qy,iq)
C
      elseif(marker.eq.'ppq')then
         iq=iq+1
         labelq(iq)='trans p K&B '
         numbq(iq)=number-1
         call field1(nmax,numw,numbq,absc,ord,qx,qy,iq)
C
      elseif(marker.eq.'pdq')then
         iq=iq+1
         labelq(iq)='trans d K&B '
         numbq(iq)=number-1
         call field1(nmax,numw,numbq,absc,ord,qx,qy,iq)
C
      elseif(marker.eq.'pfq')then
         iq=iq+1
         labelq(iq)='trans f K&B '
         numbq(iq)=number-1
         call field1(nmax,numw,numbq,absc,ord,qx,qy,iq)
C
      elseif(marker.eq.'pgq')then
         iq=iq+1
         labelq(iq)='trans g K&B '
         numbq(iq)=number-1
         call field1(nmax,numw,numbq,absc,ord,qx,qy,iq)
C
C         
C     Area markers 'p*t' for RQ = 't'
C
      elseif(marker.eq.'pst')then
         it=it+1
         labelt(it)='pst'
         numbt(it)=number-1
         call field1(nmax,numw,numbt,absc,ord,tx,ty,it)
C
      elseif(marker.eq.'ppt')then
         it=it+1
         labelt(it)='ppt'
         numbt(it)=number-1
         call field1(nmax,numw,numbt,absc,ord,tx,ty,it)
C
      elseif(marker.eq.'pdt')then 
         it=it+1
         labelt(it)='pdt'
         numbt(it)=number-1
         call field1(nmax,numw,numbt,absc,ord,tx,ty,it)
C
      elseif(marker.eq.'pft')then
         it=it+1
         labelt(it)='pft'
         numbt(it)=number-1
         call field1(nmax,numw,numbt,absc,ord,tx,ty,it)
C
      elseif(marker.eq.'pgt')then
         it=it+1
         labelt(it)='pgt'
         numbt(it)=number-1
         call field1(nmax,numw,numbt,absc,ord,tx,ty,it)
C
C
      else 
         write(4,1007)
         write(6,1007)fname
1007     format(//,' ERROR WITHIN MARKERS - LOOK AT ',a10,' !',/)
         stop
      endif
      goto 1
 300  close (unit=3)
C
C     information for program-user
C
      write(6,9991)(numbr(i),i=1,ir)
9991  format(/,' numbers of rows in columns of p*r:   ',10i6)
      write(6,9992)(numbq(k),k=1,iq)
9992  format(' numbers of rows in columns of p*q:   ',10i6)
      write(6,9993)(numbt(j),j=1,it)
9993  format(' numbers of rows in columns of p*t:   ',10i6)    
C
      return
      end
C
C      
C***   ***   ***   ***   ***   ***   ***   ***   ***   ***   ***   ***
C
      subroutine field1(nmax,numw,num,x,y,fieldx,fieldy,ifield)
C
C
      implicit double precision (a-h,o-z)
C
      dimension num(numw),x(nmax),y(nmax),fieldx(nmax,numw),
     1 fieldy(nmax,numw)       
C
C
      do 10,i=1,num(ifield)
            fieldx(i,ifield)=x(i)
            fieldy(i,ifield)=y(i)
  10  continue 
      return
      end 
C
C***   ***   ***   ***   ***   ***   ***   ***   ***   ***   ***   ***
C
C ----------------- end of reading file kbplot.dat ---------------------
C
C***   ***   ***   ***   ***   ***   ***   ***   ***   ***   ***   ***
C
      subroutine outdat (nmax,numw,
     1 rx,ry,qx,qy,tx,ty,numbr,numbq,numbt,
     2 ir,iq,it,labelr,labelq,labelt,rxmax,txmax)
C
      parameter (nscale=20)
C
      implicit double precision (a-h,o-z)
C
      dimension  rx(nmax,numw), ry(nmax,numw),
     1 qx(nmax,numw), qy(nmax,numw),tx(nmax,numw), ty(nmax,numw), 
     3 numbr(numw), numbq(numw), numbt(numw),
     4 xmaxr(nscale), xmaxt(nscale)
C
      character*12 labelr(numw),labelq(numw),labelt(numw)
C
C
C     opening files to write to
C
C  1. 'kbop.gp'.....
C  2. 'kbtr.gp'.....
C  3. 'kbar.gp'.....
C       
      open(unit=4,file='kbop.gp',status='unknown',form='formatted')
      open(unit=9,file='kbtr.gp',status='unknown',form='formatted')
      open(unit=7,file='kbar.gp',status='unknown',form='formatted')
C
C 
C     write Kleinman & Bylander operators to         'kbop.gp'
C     output format:
C                                               x;psr;...;x;pgr 
C
      write(4,1108)(labelr(i),i=1,ir)
1108  format('#     r (a.u.)',2x,a12,19(14x,a12))
      write(4,1208)
1208  format('#')
      do 75,i=1,numbr(ir) 
         write(4,1008)(rx(i,kc),ry(i,kc),kc=1,ir)
1008     format(20(3x,f10.6))
  75  continue
C
C     scaling:
C
      do 76 kc=1,ir
         i=numbr(ir)
 760     if(abs(ry(i,kc)).lt.0.0001)then
            xmaxr(kc)=rx(i,kc)
            i=i-1
            goto 760
         else
            goto 76
      endif
  76  continue
      rxmax=xmaxr(1)
      if(ir.gt.1)then
         do 77 kc=2,ir
            if(xmaxr(kc).ge.xmaxr(kc-1))rxmax=xmaxr(kc)
  77     continue
      endif                                                                            
C
C
C     write Kleinman & Bylander Transforms to       'kbtr.gp'
C     output format:
C                                              x;psq;...x;pgq   
C
      write(9,1109)(labelq(i),i=1,iq)
1109  format('#   q (1/a.u.)',3x,a12,19(14x,a12))
      write(9,1209)
1209  format('#')       
C
      do 80,i=1,numbq(iq)
            write(9,1009)(qx(i,kc),qy(i,kc),kc=1,iq)
1009        format(20(3x,f10.6))
  80  continue
C
C  
C     write Kleinman & Bylander Areas to           'kbar.gp'
C     output format:
C                                             x;pst;...x;pgt
C
      write(7,1110)(labelt(i),i=1,it)
1110  format('#      r (a.u.)',8x,a12,19(16x,a12))  
      write(7,1210)
1210  format('#')       
C
         do 85,i=1,numbt(it)
            write(7,1010)(tx(i,kc),ty(i,kc),kc=1,it)
1010        format(11(4x,f10.6))
  85     continue
C
C   
C     scaling:
C
      do 86 kc=1,it
         i=numbt(it)
 860     if(abs(ty(i,kc)).lt.0.0001)then
            xmaxt(kc)=tx(i,kc)
            i=i-1
            goto 860
         else
            goto 86
      endif
  86  continue
C
      txmax=xmaxt(1)
      if(it.gt.1)then
         do 87 kc=2,it
            if(xmaxt(kc).ge.xmaxt(kc-1))txmax=xmaxt(kc)
  87     continue
      endif
C
      return
      end                        
C
C***   ***   ***   ***   ***   ***   ***   ***   ***   ***   ***   *** 
C
      subroutine outdis(numw,
     1 ir,iq,it,labelr,labelq,labelt,rxmax,txmax)
C
C     writes the appropriate command-file to be accepted by gnuplot;
C     the file is called 'command_kb.gp'
C     the plots will be shown on the screen by gnuplot.
C
      implicit double precision (a-h,o-z)
C
      character*12 labelr(numw), labelq(numw), labelt(numw) 
C
C
      open (unit=1,file='command_kb.gp',status='unknown',
     &  form='formatted')
C
C    
C
C
C     create first sheet: 'Kleinman & Bylander Operators'
C
C     scaling:
C
      icount=0
      do 9 i=1,20
         icount=icount+1
         rxmhlp=10.0-0.5*icount
         if(rxmax.lt.rxmhlp)rmaxx=rxmhlp
   9  continue      
C
      write(1,98)
  98  format('set term x11')
      write(1,99)
  99  format('set nolabel',/,41Hset title 'Kleinman & Bylander Operators
     1')   
      write(1,100)rmaxx
 100  format('set nokey',/,'set notime',/,'set noparametric',/,'set size
     1',/,'set autoscale y',/,'set xrange [0:',f4.1,']',/,20Hset format 
     2xy '%.1f',/,'set nogrid',/,21Hset xlabel 'r [a.u.]',/,'set xtics',
     3/,'set key',/,34Hset ylabel 'V(r) [(Ry/a^3)^0.5]' 7)
      write(1,104)
 104  format('plot \')
      do 5,i=1,ir
         write(1,105)2*i-1,2*i,labelr(i)
 105  format(9H'kbop.gp',' using ',i2,':',i2,5H ti ',a12,1H',' w li,\')
   5  continue
C
C
      write(1,111)   
 111  format(15H0 ti ' ' w li 1)     
      write(1,115)
 115  format('pause -1',/,'#')
C
C
C     create second sheet: 'Kleinman & Bylander Transforms'
C
      write(1,119)
 119  format('set nolabel',/,42Hset title 'Kleinman & Bylander Transform
     1s')  
      write(1,120)
 120  format(19Hset format y '%.2f',/,'set xtics',/,'set key',/,
     1'set autoscale y',/,'set autoscale x',/,
     223Hset xlabel 'q [1/a.u.]',/,35Hset ylabel 'Bessel V(q) [Ry^0.5]' 
     37)
      write(1,121)
 121  format('set nolabel')
      write(1,124)
 124  format('plot \')
      do 25,i=1,iq
         write(1,125)2*i-1,2*i,labelq(i)
 125  format(9H'kbtr.gp',' using ',i2,':',i2,5H ti ',a12,1H',' w li,\')
  25  continue
      write(1,111)
      write(1,115)
C
C
C     create third sheet: 'Kleinman & Bylander Area'
C
C     scaling:
C
      icount=0
      do 37 i=1,20
         icount=icount+1
         txmhlp=10.0-0.5*icount
         if(txmax.lt.txmhlp)tmaxx=txmhlp
  37  continue      
C
      write(1,139)
 139  format(36Hset title 'Kleinman & Bylander Area')     
      write(1,140)tmaxx
 140  format('set nolabel',/,'set key',/,'set autoscale y',
     1/,'set xrange [0:',f4.1,']',/,21Hset xlabel 'r [a.u.]',/,
     224Hset ylabel 'V(r) [Ry]' 5,/,20Hset format xy '%.1f',/,
     3'set xtics')
      write(1,144)
 144  format('plot \')     
      do 45,i=1,it
         write(1,145)2*i-1,2*i,labelt(i)
 145     format(9H'kbar.gp',' using ',i2,':', i2,5H ti ',a5,
     1          1H',' w li,\')
  45  continue
      write(1,111)
C
C
C
      return
      end
C
C***   ***   ***   ***   ***   ***   ***   ***   ***   ***   ***   ***
C
      subroutine outps (numw,
     1 ir,iq,it,labelr,labelq,labelt,rxmax,txmax)
C
C     writes the appropriate command-file to be accepted by gnuplot;
C     the file is called 'command_kb.gp'
C     gnuplot will create postscript-files.
C
      implicit double precision (a-h,o-z)
C
      character*12 labelr(numw), labelq(numw), labelt(numw) 
C
C
      open (unit=1,file='command_kb.gp',status='unknown',
     &  form='formatted')
C
C
C     create first sheet: 'Kleinman & Bylander Operators'
C
C     scaling:
C
      icount=0
      do 9 i=1,20
         icount=icount+1
         rxmhlp=10.0-0.5*icount
         if(rxmax.lt.rxmhlp)rmaxx=rxmhlp
   9  continue      
C
      write(1,97)
  97  format(20Hset output 'kbop.ps')    
      write(1,98)
  98  format(42Hset term postscript monochrome 'Helvetica')
      write(1,99)
  99  format('set nolabel',/,41Hset title 'Kleinman & Bylander Operators
     1')   
      write(1,100)rmaxx
 100  format('set nokey',/,'set notime',/,'set noparametric',/,'set size
     1',/,'set autoscale y',/,'set xrange [0:',f4.1,']',/,20Hset format 
     2xy '%.1f',/,'set nogrid',/,21Hset xlabel 'r [a.u.]',/,'set xtics',
     3/,'set key',/,34Hset ylabel 'V(r) [(Ry/a^3)^0.5]' 1)
      write(1,104)
 104  format('plot \')
      do 5,i=1,ir
         write(1,105)2*i-1,2*i,labelr(i)
 105  format(9H'kbop.gp',' using ',i2,':',i2,5H ti ',a12,1H',' w li,\')
   5  continue
C
      write(1,111)   
 111  format(15H0 ti ' ' w li 1)     
      write(1,115)
 115  format('#')
C
C
C     create second sheet: 'Kleinman & Bylander Transforms'
C
      write(1,117)
 117  format(20Hset output 'kbtr.ps')    
      write(1,119)
 119  format('set nolabel',/,42Hset title 'Kleinman & Bylander Transform
     1s')  
      write(1,120)
 120  format(19Hset format y '%.2f',/,'set xtics',/,'set key',/,
     1'set autoscale y',/,'set autoscale x',/,
     223Hset xlabel 'q [1/a.u.]',/,35Hset ylabel 'Bessel V(q) [Ry^0.5]' 
     31)
      write(1,121)
 121  format('set nolabel')
      write(1,124)
 124  format('plot \')
      do 25,i=1,iq
         write(1,125)2*i-1,2*i,labelq(i)
 125  format(9H'kbtr.gp',' using ',i2,':',i2,5H ti ',a12,1H',' w li,\')
  25  continue
      write(1,111)
      write(1,115)
C
C
C     create third sheet: 'Kleinman & Bylander Area'
C
C     scaling:
C
      icount=0
      do 37 i=1,20
         icount=icount+1
         txmhlp=10.0-0.5*icount
         if(txmax.lt.txmhlp)tmaxx=txmhlp
  37  continue      
C
      write(1,137)
 137  format(20Hset output 'kbar.ps')    
      write(1,139)
 139  format(36Hset title 'Kleinman & Bylander Area')     
      write(1,140)tmaxx
 140  format('set nolabel',/,'set key',/,'set autoscale y',
     1/,'set xrange [0:',f4.1,']',/,21Hset xlabel 'r [a.u.]',/,
     224Hset ylabel 'V(r) [Ry]' 5,/,20Hset format xy '%.1f',/,
     3'set xtics')
      write(1,144)
 144  format('plot \')     
      do 45,i=1,it
         write(1,145)2*i-1,2*i,labelt(i)
 145     format(9H'kbar.gp',' using ',i2,':',i2,5H ti ',a5,
     1          1H',' w li,\')
  45  continue
      write(1,111)
C
C
C
      return
      end
C
C***   ***   ***   ***   ***   ***   ***   ***   ***   ***   ***   ***
C
      subroutine outpbm(numw,
     1 ir,iq,it,labelr,labelq,labelt,rxmax,txmax)
C
C     writes the appropriate command-file to be accepted by gnuplot;
C     the file is called 'command_kb.gp'
C     gnuplot will create portable bitmap-files.
C
      implicit double precision (a-h,o-z)
C
      character*12 labelr(numw), labelq(numw), labelt(numw) 
C
C
      open (unit=1,file='command_kb.gp',status='unknown',
     &   form='formatted')
C
C
C     create first sheet: 'Kleinman & Bylander Operators'
C
C     scaling:
C
      icount=0
      do 9 i=1,20
         icount=icount+1
         rxmhlp=10.0-0.5*icount
         if(rxmax.lt.rxmhlp)rmaxx=rxmhlp
   9  continue      
C
      write(1,97)
  97  format(21Hset output 'kbop.pbm')    
      write(1,98)
  98  format('set term pbm')
      write(1,99)
      write(1,99)
  99  format('set nolabel',/,41Hset title 'Kleinman & Bylander Operators
     1')   
      write(1,100)rmaxx
 100  format('set nokey',/,'set notime',/,'set noparametric',/,'set size
     1',/,'set autoscale y',/,'set xrange [0:',f4.1,']',/,20Hset format 
     2xy '%.1f',/,'set nogrid',/,21Hset xlabel 'r [a.u.]',/,'set xtics',
     3/,'set key',/,34Hset ylabel 'V(r) [(Ry/a^3)^0.5]' 1)
      write(1,104)
 104  format('plot \')
      do 5,i=1,ir
         write(1,105)2*i-1,2*i,labelr(i)
 105  format(9H'kbop.gp',' using ',i2,':',i2,5H ti ',a12,1H',' w li,\')
   5  continue
C
      write(1,111)   
 111  format(15H0 ti ' ' w li 1)     
      write(1,115)
 115  format('#')
C
C
C     create second sheet: 'Kleinman & Bylander Transforms'
C
      write(1,117)
 117  format(21Hset output 'kbtr.pbm')    
      write(1,119)
 119  format('set nolabel',/,42Hset title 'Kleinman & Bylander Transform
     1s')  
      write(1,120)
 120  format(19Hset format y '%.2f',/,'set xtics',/,'set key',/,
     1'set autoscale y',/,'set autoscale x',/,
     223Hset xlabel 'q [1/a.u.]',/,35Hset ylabel 'Bessel V(q) [Ry^0.5]' 
     31)
      write(1,121)
 121  format('set nolabel')
      write(1,124)
 124  format('plot \')
      do 25,i=1,iq
         write(1,125)2*i-1,2*i,labelq(i)
 125  format(9H'kbtr.gp',' using ',i2,':',i2,5H ti ',a12,1H',' w li,\')
  25  continue
      write(1,111)
      write(1,115)
C
C
C     create third sheet: 'Kleinman & Bylander Area'
C
C     scaling:
C
      icount=0
      do 37 i=1,20
         icount=icount+1
         txmhlp=10.0-0.5*icount
         if(txmax.lt.txmhlp)tmaxx=txmhlp
  37  continue      
C
      write(1,137)
 137  format(21Hset output 'kbar.pbm')    
      write(1,139)
 139  format(36Hset title 'Kleinman & Bylander Area')     
      write(1,140)tmaxx
 140  format('set nolabel',/,'set key',/,'set autoscale y',
     1/,'set xrange [0:',f4.1,']',/,21Hset xlabel 'r [a.u.]',/,
     224Hset ylabel 'V(r) [Ry]' 5,/,20Hset format xy '%.1f',/,
     3'set xtics')
      write(1,144)
 144  format('plot \')     
      do 45,i=1,it
         write(1,145)2*i-1,2*i,labelt(i)
 145     format(9H'kbar.gp',' using ',i2,':',i2,5H ti ',a5,
     1          1H',' w li,\')
  45  continue
      write(1,111)
C
C
C
      return
      end
C
C***   ***   ***   ***   ***   ***   ***   ***   ***   ***   ***   ***
