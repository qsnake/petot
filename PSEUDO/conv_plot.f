      program convert                                            
C
C     *********************************************************************** 
C     *                                                                     *
C     * convert transforms the outputfile (e.g.'plot.dat**') of the program *  
C     * 'atm.f' into a form that can be used to create the corresponding    * 
C     * plots with 'gnuplot'                                                *
C     *                                                                     *
C     * written by Peter Schuster, April 1993                               *
C     *                                                                     * 
C     ***********************************************************************
C
      parameter (nmax=5000,numw=20,nump=10,nftp=50)
C      
      implicit double precision (a-h,o-z)
C
      dimension  wx(nmax,numw), wy(nmax,numw),
     1 px(nmax,numw), py(nmax,numw),vx(nmax,nump), vy(nmax,nump), 
     2 fx(nftp,nump), fy(nftp,nump), fwx(nmax,nump), fwy(nmax,nump),
     3 numbw(numw), numbp(numw), numbv(nump), numbf(nump),numbfw(nump)
C
      character* 3 howdat 
      character*12 labelw(numw),labelp(numw),labelv(nump),labelf(nump),
     1             labeli(nump)  
C
      call input (nmax,numw,nump,nftp,
     1 wx,wy,px,py,vx,vy,fx,fy,fwx,fwy,numbw,numbp,
     2 numbv,numbf,numbfw,iw,ip,iv,ifp,ifw,labelw,labelp,labelv,
     3 labelf,labeli,zion)      
C
C
      iwsave=iw
      ipsave=ip
C
      call outdat (nmax,numw,nump,nftp,
     1 wx,wy,px,py,vx,vy,fx,fy,fwx,fwy,numbw,numbp,
     2 numbv,numbf,numbfw,iw,ip,iv,ifp,ifw,labelw,labelp,labelv,
     3 labelf,labeli,yminwp,ymaxwp,yminv)      
C
      iw=iwsave
      ip=ipsave
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
      call outdis(nump,
     1 iw,ip,ifw,iv,labelv,ifp,labelf,yminwp,ymaxwp,yminv,zion)
C      
      elseif(howdat.eq.'ps'.or.howdat.eq.' ps'.or.howdat.eq.'PS '.or.ho
     1wdat.eq.' PS')then               
      call outps(nump,
     1 iw,ip,ifw,iv,labelv,ifp,labelf,yminwp,ymaxwp,yminv,zion)
C
      elseif(howdat.eq.'pbm'.or.howdat.eq.'Pbm'.or.howdat.eq.'PBM')then
      call outpbm(nump,
     1 iw,ip,ifw,iv,labelv,ifp,labelf,yminwp,ymaxwp,yminv,zion)
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
      subroutine input (nmax,numw,nump,nftp,
     1 wx,wy,px,py,vx,vy,fx,fy,fwx,fwy,numbw,numbp,
     2 numbv,numbf,numbfw,iw,ip,iv,ifp,ifw,labelw,labelp,labelv,
     3 labelf,labeli,zion)      
C
C     reads the input from <plot.dat**>
C 
      implicit double precision (a-h,o-z)
C
      parameter (mmax=5000)
C
      dimension absc(mmax), ord(mmax), wx(nmax,numw), wy(nmax,numw),
     1 px(nmax,numw), py(nmax,numw),vx(nmax,nump), vy(nmax,nump), 
     2 fx(nftp,nump), fy(nftp,nump), fwx(nmax,nump), fwy(nmax,nump),
     3 numbw(numw), numbp(numw), numbv(nump), numbf(nump),numbfw(nump)
C
      character*3 marker
      character*10 fname
      character*12 labelw(numw),labelp(numw),labelv(nump),labelf(nump),
     1             labeli(nump)  
      
C
C
C
      iw=0
      ip=0
      iv=0
      if=0
      ifw=0
C
C
C     asking for outputfile to be used
C            
      write(*,1001)
1001  format(/,' Which file <plot.dat__> shall be read from?',/,' Please 
     1 give the name. (format = a10)',/,' Press <ENTER> to choose the de
     2fault <plot.dat01>') 
      read(*,1002)fname
1002  format(a10)
      if(fname.eq.' ')fname='plot.dat01'
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
C1003 format(1x,f7.4,3x,f18.14)
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
C     look if 'ion-change marker'
C
      if(marker.eq.'zio')then
         read(3,1006)zion
1006     format(2x,f5.2)
C
C
C     true wave-function markers
C
      elseif(marker.eq.'wst')then
         iw=iw+1
         labelw(iw)='true wf s'
         numbw(iw)=number-1
         call field1(nmax,numw,numbw,absc,ord,wx,wy,iw)
C
      elseif(marker.eq.'wpt')then
         iw=iw+1
         labelw(iw)='true wf p'
         numbw(iw)=number-1
         call field1(nmax,numw,numbw,absc,ord,wx,wy,iw)
C
      elseif(marker.eq.'wdt')then
         iw=iw+1
         labelw(iw)='true wf d'
         numbw(iw)=number-1
         call field1(nmax,numw,numbw,absc,ord,wx,wy,iw)
C
      elseif(marker.eq.'wft')then
         iw=iw+1
         labelw(iw)='true wf p'
         numbw(iw)=number-1
         call field1(nmax,numw,numbw,absc,ord,wx,wy,iw)
C
      elseif(marker.eq.'wgt')then
         iw=iw+1
         labelw(iw)='true wf p'
         numbw(iw)=number-1
         call field1(nmax,numw,numbw,absc,ord,wx,wy,iw)
C
C         
C     pseudo wave-function markers
C
      elseif(marker.eq.'wsp')then
         ip=ip+1
         labelp(ip)='pseudo wf s'
         numbp(ip)=number-1
         call field1(nmax,numw,numbp,absc,ord,px,py,ip)
C
      elseif(marker.eq.'wpp')then
         ip=ip+1
         labelp(ip)='pseudo wf p'
         numbp(ip)=number-1
         call field1(nmax,numw,numbp,absc,ord,px,py,ip)
C
      elseif(marker.eq.'wdp')then
         ip=ip+1
         labelp(ip)='pseudo wf d'
         numbp(ip)=number-1
         call field1(nmax,numw,numbp,absc,ord,px,py,ip)
C
      elseif(marker.eq.'wfp')then
         ip=ip+1
         labelp(ip)='pseudo wf f'
         numbp(ip)=number-1
         call field1(nmax,numw,numbp,absc,ord,px,py,ip)
C
      elseif(marker.eq.'wgp')then
         ip=ip+1
         labelp(ip)='pseudo wf g'
         numbp(ip)=number-1
         call field1(nmax,numw,numbp,absc,ord,px,py,ip)
C
C         
C     potential markers
C
      elseif(marker.eq.'vns')then
         iv=iv+1
         labelv(iv)='s - nonlocal$'
         numbv(iv)=number-1
         call field1(nmax,nump,numbv,absc,ord,vx,vy,iv)
C
      elseif(marker.eq.'vnp')then
         iv=iv+1
         labelv(iv)='p - nonlocal$'
         numbv(iv)=number-1
         call field1(nmax,nump,numbv,absc,ord,vx,vy,iv)
C
      elseif(marker.eq.'vnd')then 
         iv=iv+1
         labelv(iv)='d - nonlocal$'
         numbv(iv)=number-1
         call field1(nmax,nump,numbv,absc,ord,vx,vy,iv)
C
      elseif(marker.eq.'vnf')then
         iv=iv+1
         labelv(iv)='f - nonlocal$'
         numbv(iv)=number-1
         call field1(nmax,nump,numbv,absc,ord,vx,vy,iv)
C
      elseif(marker.eq.'vng')then
         iv=iv+1
         labelv(iv)='g - nonlocal$'
         numbv(iv)=number-1
         call field1(nmax,nump,numbv,absc,ord,vx,vy,iv)
C
C
C     fourier markers - pseudopotentials
C
      elseif(marker.eq.'fn1')then
         ifp=ifp+1
         labelf(ifp)='s - nonlocal$'
         numbf(ifp)=number-1
         call field1(nftp,nump,numbf,absc,ord,fx,fy,ifp)
C
      elseif(marker.eq.'fn2')then
         ifp=ifp+1
         labelf(ifp)='p - nonlocal$'
         numbf(ifp)=number-1
         call field1(nftp,nump,numbf,absc,ord,fx,fy,ifp)
C
      elseif(marker.eq.'fn3')then
         ifp=ifp+1
         labelf(ifp)='d - nonlocal$'
         numbf(ifp)=number-1
         call field1(nftp,nump,numbf,absc,ord,fx,fy,ifp)
C
      elseif(marker.eq.'fn4')then
         ifp=ifp+1
         labelf(ifp)='f - nonlocal$'
         numbf(ifp)=number-1
         call field1(nftp,nump,numbf,absc,ord,fx,fy,ifp)
C
      elseif(marker.eq.'fn5')then
         ifp=ifp+1
         labelf(ifp)='g - nonlocal$'
         numbf(ifp)=number-1
         call field1(nftp,nump,numbf,absc,ord,fx,fy,ifp)
C
C
C     fourier markers - pseudo wave-functions
C
      elseif(marker.eq.'fw0')then
         ifw=ifw+1
         labeli(ifw)='FT-wf s'
         numbfw(ifw)=number-1
         call field1(nmax,nump,numbfw,absc,ord,fwx,fwy,ifw)
C
      elseif(marker.eq.'fw1')then
         ifw=ifw+1
         labeli(ifw)='FT-wf p'
         numbfw(ifw)=number-1
         call field1(nmax,nump,numbfw,absc,ord,fwx,fwy,ifw)
C
      elseif(marker.eq.'fw2')then
         ifw=ifw+1
         labeli(ifw)='FT-wf d'
         numbfw(ifw)=number-1
         call field1(nmax,nump,numbfw,absc,ord,fwx,fwy,ifw)
C
      elseif(marker.eq.'fw3')then
         ifw=ifw+1
         labeli(ifw)='FT-wf f'
         numbfw(ifw)=number-1
         call field1(nmax,nump,numbfw,absc,ord,fwx,fwy,ifw)
C
      elseif(marker.eq.'fw4')then
         ifw=ifw+1
         labeli(ifw)='FT-wf g'
         numbfw(ifw)=number-1
         call field1(nmax,nump,numbfw,absc,ord,fwx,fwy,ifw)
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
      write(6,9991)(numbw(i),i=1,iw)
9991  format(/,' numbers of rows in columns of true wave-functions:   '
     1,10i6) 
      write(6,9992)(numbp(k),k=1,ip)
9992  format(' numbers of rows in columns of pseudo wave-functions: '
     1,10i6)
      write(6,9993)(numbv(j),j=1,iv)
9993  format(' numbers of rows in columns of pseudopotentials:      '     
     1,10i6)
      write(6,9994)(numbf(l),l=1,ifp)
9994  format(' numbers of rows in columns of FT - ppot:             '
     1,10i6)
      write(6,9995)(numbfw(m),m=1,ifw)
9995  format(' numbers of rows in columns of FT - wafe-functions:   '
     1,10i6)
C
      return
      end
C
C      
C***   ***   ***   ***   ***   ***   ***   ***   ***   ***   ***   ***
C
      subroutine field1(nmax,ncol,num,x,y,fieldx,fieldy,ifield)
C
C
      implicit double precision (a-h,o-z)
C
      dimension num(ncol),x(nmax),y(nmax),fieldx(nmax,ncol),
     1 fieldy(nmax,ncol)       
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
C ----------------- end of reading file plot.dat** ---------------------
C
C***   ***   ***   ***   ***   ***   ***   ***   ***   ***   ***   ***
C
      subroutine outdat (nmax,numw,nump,nftp,
     1 wx,wy,px,py,vx,vy,fx,fy,fwx,fwy,numbw,numbp,
     2 numbv,numbf,numbfw,iw,ip,iv,ifp,ifw,labelw,labelp,labelv,
     3 labelf,labeli,yminwp,ymaxwp,yminv)      
C
C
      implicit double precision (a-h,o-z)
C
      dimension  wx(nmax,numw), wy(nmax,numw),
     1 px(nmax,numw), py(nmax,numw),vx(nmax,nump), vy(nmax,nump), 
     2 fx(nftp,nump), fy(nftp,nump), fwx(nmax,nump), fwy(nmax,nump),
     3 numbw(numw), numbp(numw), numbv(nump), numbf(nump),numbfw(nump)
C
      character*12 labelw(numw),labelp(numw),labelv(nump),labelf(nump),
     1             labeli(nump)  
C
C
C     opening files to write to
C
C  1. 'wfct.gp'.....ae and pseudo wave-functions
C  2. 'fwfct.gp'....ae and pseudo wave-functions + FT's(pwf's)
C  3. 'ppot.gp'.....pseudo potentials
C  4. 'fppot.gp'....pseudo potentials + FT's(pp's)
C       
      open(unit=4,file='wfct.gp',status='unknown',form='formatted')
      open(unit=9,file='fwfct.gp',status='unknown',form='formatted')
      open(unit=7,file='ppot.gp',status='unknown',form='formatted')
      open(unit=8,file='fppot.gp',status='unknown',form='formatted')
C
C 
C     write true wave-functions and pseudo wafe-functions to 'wfct.gp'
C     output format:
C                    x; wst;wsp;(wst;wsp);...;wgp;wgt;(wgt;wgp)
C
      if(iw.eq.0.and.ip.eq.0)goto 91
      if(iw.lt.ip)iw=ip
      write(4,1108)(labelw(i),labelp(i),i=1,iw)
1108  format('#     r (a.u.)',3x,20(1x,a12))
      write(4,1208)
1208  format('#')
      yminwp=0.0
      ymaxwp=0.0    
      do 75,i=1,numbw(iw) 
         write(4,1008)wx(i,1),(wy(i,kc),py(i,kc),kc=1,iw)
1008     format(20(3x,f10.6))
         do 750 k=1,iw
            if(wy(i,k).lt.yminwp) yminwp=wy(i,k)
            if(wy(i,k).gt.ymaxwp) ymaxwp=wy(i,k)
            if(py(i,k).lt.yminwp) yminwp=py(i,k)
            if(py(i,k).gt.ymaxwp) ymaxwp=py(i,k)
 750     continue            
  75  continue                                                              
C
C
C     write FT(pseudo wave-function)                     to 'fwfct.gp'
C     output format:
C                    x; fw0;(fw0);...;fw4;(fw4)   
C
      write(9,1109)(labeli(i),i=1,ifw)
1109  format('#   q (1/a.u.)',3x,10(1x,a12))
      write(9,1209)
1209  format('#')       
C
  91  if(ifw.eq.0)then
         write(9,*)' NO NUMBERS'
         goto 92
      else   
         do 80,i=1,numbfw(ifw)
            write(9,1009)fwx(i,1),(fwy(i,kc),kc=1,ifw)
1009        format(20(3x,f10.6))
  80  continue
      endif                  
C
C
C     write ionic pseudopotentials                        to 'ppot.gp'
C     output format:
C                    x; vn0;(vn0);...;vn4;(vn4)
C
      write(7,1110)(labelv(i),i=1,iv)
1110  format('#      r (a.u.)',1x,10(3x,a12))   
      write(7,1210)
1210  format('#')       
C
  92  if(iv.eq.0)then
         write(7,*)' NO NUMBERS'
         goto 93
      else   
         yminv=0.0   
         do 85,i=1,numbv(iv)
            write(7,1010)vx(i,1),(vy(i,kc),kc=1,iv)
1010        format(11(4x,f11.6))
            do 850 k=1,iv
               if(vy(i,k).lt.yminv) yminv=vy(i,k)
 850        continue               
  85     continue 
      endif
C
C
C     write FT(ionic pseudopotentials)                   to 'fppot.gp'
C     output format:
C                    x; fn1;(fn1);...;fn5;(fn5)
C
      write(8,1111)(labelf(i),i=1,ifp)
1111  format('#    q (1/a.u.)',1x,10(2x,a12))  
      write(8,1211)
1211  format('#')            
C
  93  if(ifp.eq.0)then
         write(8,*)' NO NUMBERS'
         goto 94
      else
         do 90,i=1,numbf(ifp)
            write(8,1011)fx(i,1),(fy(i,kc),kc=1,ifp)
1011        format(20(4x,f10.6))
  90     continue          
      endif         
C
C
  94  return
      end                        
C
C***   ***   ***   ***   ***   ***   ***   ***   ***   ***   ***   *** 
C
      subroutine outdis(nump,
     1 iw,ip,ifw,iv,labelv,ifp,labelf,yminwp,ymaxwp,yminv,zion)
C
C     writes the appropriate command-file to be accepted by gnuplot;
C     the file is called 'command.gp'
C     the plots will be shown on the screen by gnuplot.
C
      implicit double precision (a-h,o-z)
C
      character*12 labelv(nump), labelf(nump) 
C
C
      open (unit=1,file='command.gp',status='unknown',form='formatted')
C
C
C     create first sheet: true and pseudo wafe-functions
C
C     scaling:
C
      icount=1
      yhelp=0.5+yminwp
 410  if(yhelp.gt.0)then
         goto 415
      else
         icount=icount+1
         yhelp=0.5*icount+yminwp
         goto 410
      endif   
 415  yminwp=-0.5*icount              
C
      do 420 i=1,20
         yhelp=10.0-0.5*i
         if(ymaxwp.lt.yhelp)yhmwp=yhelp
 420  continue
      ymaxwp=yhmwp         
C
      write(1,98)
  98  format('set term x11')
      write(1,99)
  99  format('set nolabel',/,26Hset title 'Wave Functions')   
      write(1,100)yminwp,ymaxwp
 100  format('set nokey',/,'set notime',/,'set noparametric',/,'set size
     1',/,'set yrange [',f5.1,':',f5.1,']',/,'set xrange [0:5.5]',/,20Hs
     2et format xy '%.1f',/,'set nogrid',/,21Hset xlabel 'r [a.u.]',/,'s
     3et xtics 0,1',/,20Hset ylabel 'rR(r)' 7)
      write(1,101)
 101  format(44Hset label ' true wave functions ' at 3.7,5.6,/,
     144Hset label '_____________________' at 3.7,5.4,/,
     244Hset label 'pseudo wave functions' at 3.7,4.6,/,
     344Hset label '.....................' at 3.7,4.3)      
      write(1,104)
 104  format('plot',9H'wfct.gp',' using 1:2 with lines  ,',1H\)
      if(iw.gt.1)then
      do 5,i=4,2*iw,2
         write(1,105)i
 105     format(9H'wfct.gp',' using 1:',i2,' with lines 1,',1H\)
   5  continue
      else
      goto 1000
      endif
1000  do 10,i=3,2*ip+1,2
         write(1,110)i
 110     format(9H'wfct.gp',' using 1:',i2,' with lines 5,',1H\)
  10  continue
      write(1,111)
 111  format(15H0 ti ' ' w li 5)     
      write(1,115)
 115  format('pause -1',/,'#')
C
C
C     create second sheet: fourier transforms of wave-functions
C
      write(1,119)
 119  format(36Hset title 'Wave Function Transforms')     
      write(1,120)
 120  format(19Hset format y '%.2f',/,'set xtics 0,5',/,
     1'set yrange [-1.5:1.5]',/,'set xrange [0:10]',/,
     223Hset xlabel 'q [1/a.u.]',/,19Hset ylabel 'R(q)' 7)
      write(1,121)
 121  format('set nolabel')
      write(1,124)
 124  format('plot',10H'fwfct.gp',' using 1:2 with lines 4,',1H\)          
      do 15,i=3,ifw+1
         write(1,125)i
 125     format(10H'fwfct.gp',' using 1:',i2,' with lines 4,'1H\)
  15  continue
      write(1,112)
 112  format(15H0 ti ' ' w li 1)     
      write(1,115)
C
C
C     create third sheet: pseudopotentials
C
C     scaling:
C
      iyminv=int(yminv-1.)
      ihelp=-iyminv
      do 500 i=1,5
         if(mod(ihelp,5).ne.0)ihelp=ihelp+1
 500  continue     
      iyminv=-ihelp         
C
      write(1,139)
 139  format(28Hset title 'Pseudopotentials')     
      write(1,140)iyminv/2.,iyminv
 140  format('set nolabel',/,'set key 3.5,',f7.1,/,'set yrange [',i4,':
     15]',/,'set xrange [0:4]',/,21Hset xlabel 'r [a.u.]',/,
     224Hset ylabel 'V(r) [Ry]' 5,/,20Hset format xy '%.1f',/,
     3'set xtics 0,1.0')
      write(1,144)
 144  format('plot ',1H\)     
      do 25,i=2,iv+1
         write(1,145)i,labelv(i-1),i-1
 145     format(9H'ppot.gp',' using 1:',i2,5H ti ',a12,1H',
     1          ' w li',i2,',',1H\)
  25  continue
      write(1,147)zion*2.
 147  format('-',f5.2,'/x',37H ti 'Z - ion' w li 5, 0 ti ' ' w li 1)      
      write(1,115)
C
C
C     create fourth sheet: fourier transforms of pseudopotentials
C
      write(1,160)
 160  format('set nolabel',/,'set key 10,-0.4',/,
     1'set size 2.8/5.,3/3.',/,30Hset title 'Fourier Transforms',
     2/,'set xrange [0:12]',/,'set yrange [-1.5:1.5]',/,23Hset xlabel 'q
     3 [1/a.u.]',/,36Hset ylabel 'q^2/4piZion V(q) [Ry]' 4,/,'set xtics 
     40,5.0')                             
      write(1,144)
      do 35,i=2,ifp+1
        write(1,165)i,labelf(i-1),i-1
 165    format(10H'fppot.gp',' using 1:',i2,5H ti ',a12,1H',
     1         ' w li',i2,',',1H\)
  35  continue
      write(1,112)
C
C
      return
      end
C
C***   ***   ***   ***   ***   ***   ***   ***   ***   ***   ***   ***
C
      subroutine outps(nump,
     1 iw,ip,ifw,iv,labelv,ifp,labelf,yminwp,ymaxwp,yminv,zion)
C
C     writes the appropriate command-file to be accepted by gnuplot;
C     the file is called 'command.gp'
C     gnuplot will create postscript-files <*.ps> out of the plots
C
      implicit double precision (a-h,o-z)
C
      character*12 labelv(nump), labelf(nump) 
C
      open (unit=1,file='command.gp',status='unknown',form='formatted')
C
C
C     create first sheet: true and pseudo wafe-functions
C
C     scaling:
C
      icount=1
      yhelp=0.5+yminwp
 410  if(yhelp.gt.0)then
         goto 415
      else
         icount=icount+1
         yhelp=0.5*icount+yminwp
         goto 410
      endif   
 415  yminwp=-0.5*icount              
C
      do 420 i=1,20
         yhelp=10.0-0.5*i
         if(ymaxwp.lt.yhelp)yhmwp=yhelp
 420  continue
      ymaxwp=yhmwp         
C
      write(1,98)
  98  format('set term postscript monochrome ',11H'Helvetica',/,20Hset o
     1utput 'wfct.ps')    
      write(1,99)
  99  format('set nolabel',/,26Hset title 'Wave Functions')   
      write(1,100)yminwp,ymaxwp
 100  format('set nokey',/,'set notime',/,'set noparametric',/,'set size
     1',/,'set yrange [',f5.1,':',f5.1,']',/,'set xrange [0:5.5]',/,20Hs
     2et format xy '%.1f',/,'set nogrid',/,21Hset xlabel 'r [a.u.]',/,'s
     3et xtics 0,1',/,20Hset ylabel 'rR(r)' 1)
      write(1,101)
 101  format(44Hset label ' true wave functions ' at 3.7,5.6,/,
     144Hset label '_____________________' at 3.7,5.4,/,
     244Hset label 'pseudo wave functions' at 3.7,4.6,/,
     344Hset label '.....................' at 3.7,4.3)      
      write(1,104)
 104  format('plot',9H'wfct.gp',' using 1:2 with lines 1,',1H\)
      if(iw.gt.1)then
      do 5,i=4,2*iw,2
         write(1,105)i
 105     format(9H'wfct.gp',' using 1:',i2,' with lines 1,',1H\)
   5  continue
      else
      goto 1000
      endif
1000  do 10,i=3,2*ip+1,2
         write(1,110)i
 110     format(9H'wfct.gp',' using 1:',i2,' with lines 5,',1H\)
  10  continue
      write(1,111)
 111  format(15H0 ti ' ' w li 5,/,'#')     
C
C
C     create second sheet: fourier transforms of wave-functions
C
      write(1,118)
 118  format(21Hset output 'fwfct.ps')
      write(1,119)
 119  format(36Hset title 'Wave Function Transforms')     
      write(1,120)
 120  format(19Hset format y '%.2f',/,'set xtics 0,5',/,
     1'set yrange [-1.5:1.5]',/,'set xrange [0:10]',/,
     223Hset xlabel 'q [1/a.u.]',/,19Hset ylabel 'R(q)' 1)
      write(1,121)
 121  format('set nolabel')
      write(1,124)
 124  format('plot',10H'fwfct.gp',' using 1:2 with lines 4,',1H\)          
      do 15,i=3,ifw+1
         write(1,125)i
 125     format(10H'fwfct.gp',' using 1:',i2,' with lines 4,',1H\)
  15  continue
      write(1,112)
 112  format(15H0 ti ' ' w li 1,/,'#')     
C
C
C     create third sheet: pseudopotentials
C
C     scaling:
C
      iyminv=int(yminv-1.)
      ihelp=-iyminv
      do 500 i=1,5
         if(mod(ihelp,5).ne.0)ihelp=ihelp+1
 500  continue     
      iyminv=-ihelp         
C
      write(1,138)
 138  format(20Hset output 'ppot.ps')
      write(1,139)
 139  format(28Hset title 'Pseudopotentials')     
      write(1,140)iyminv/2.,iyminv
 140  format('set nolabel',/,'set key 3.5,',f7.1,/,'set yrange [',i4,':
     15]',/,'set xrange [0:4]',/,21Hset xlabel 'r [a.u.]',/,
     224Hset ylabel 'V(r) [Ry]' 1,/,20Hset format xy '%.1f',/,
     3'set xtics 0,1.0')
      write(1,144)
 144  format('plot ',1H\)     
      do 25,i=2,iv+1
         write(1,145)i,labelv(i-1),i-1
 145     format(9H'ppot.gp',' using 1:',i2,5H ti ',
     1         a12,1H',' w li',i2,',',1H\)
  25  continue
      write(1,147)zion*2.
 147  format('-',f5.2,'/x',37H ti 'Z - ion' w li 5, 0 ti ' ' w li 1)      
      write(1,148)
 148  format('#')
C
C
C     create fourth sheet: fourier transforms of pseudopotentials
C
      write(1,158)
 158  format(21Hset output 'fppot.ps')
      write(1,160)
 160  format('set nolabel',/,'set key 10,-0.4',/,
     1'set size 2.8/5.,3/3.',/,30Hset title 'Fourier Transforms',
     2/,'set xrange [0:12]',/,'set yrange [-1.5:1.5]',/,23Hset xlabel 'q
     3 [1/a.u.]',/,36Hset ylabel 'q^2/4piZion V(q) [Ry]' 1,/,'set xtics 
     40,5.0')                             
      write(1,144)
      do 35,i=2,ifp+1
         write(1,165)i,labelf(i-1),i-1
 165     format(10H'fppot.gp',' using 1:',i2,5H ti ',a12,1H',
     1          ' w li',i2,',',1H\)
  35  continue
      write(1,112)
C
C
      return
      end
C***   ***   ***   ***   ***   ***   ***   ***   ***   ***   ***   ***
C
      subroutine outpbm(nump,
     1 iw,ip,ifw,iv,labelv,ifp,labelf,yminwp,ymaxwp,yminv,zion)
C
C     writes the appropriate command-file to be accepted by gnuplot;
C     the file is called 'command.gp'
C     gnuplot will create pbm-files <*.pbm> out of the plots
C
      implicit double precision (a-h,o-z)
C
      character*12 labelv(nump), labelf(nump) 
C
      open (unit=1,file='command.gp',status='unknown',form='formatted')
C
C
C     create first sheet: true and pseudo wafe-functions
C
C     scaling:
C
      icount=1
      yhelp=0.5+yminwp
 410  if(yhelp.gt.0)then
         goto 415
      else
         icount=icount+1
         yhelp=0.5*icount+yminwp
         goto 410
      endif   
 415  yminwp=-0.5*icount              
C
      do 420 i=1,20
         yhelp=10.0-0.5*i
         if(ymaxwp.lt.yhelp)yhmwp=yhelp
 420  continue
      ymaxwp=yhmwp         
C
      write(1,98)
  98  format('set term pbm',/,21Hset output 'wfct.pbm')    
      write(1,99)
  99  format('set nolabel',/,26Hset title 'Wave Functions')   
      write(1,100)yminwp,ymaxwp
 100  format('set nokey',/,'set notime',/,'set noparametric',/,'set size
     1',/,'set yrange [',f5.1,':',f5.1,']',/,'set xrange [0:5.5]',/,20Hs
     2et format xy '%.1f',/,'set nogrid',/,21Hset xlabel 'r [a.u.]',/,'s
     3et xtics 0,1',/,20Hset ylabel 'rR(r)' 1)
      write(1,101)
 101  format(44Hset label ' true wave functions ' at 3.7,5.6,/,
     144Hset label '_____________________' at 3.7,5.4,/,
     244Hset label 'pseudo wave functions' at 3.7,4.6,/,
     344Hset label '.....................' at 3.7,4.3)      
      write(1,104)
 104  format('plot',9H'wfct.gp',' using 1:2 with lines 1,',1H\)
      if(iw.gt.1)then
      do 5,i=4,2*iw,2
         write(1,105)i
 105     format(9H'wfct.gp',' using 1:',i2,' with lines 1,',1H\)
   5  continue
      else
      goto 1000
      endif
1000  do 10,i=3,2*ip+1,2
         write(1,110)i
 110     format(9H'wfct.gp',' using 1:',i2,' with lines 5,',1H\)
  10  continue
      write(1,111)
 111  format(15H0 ti ' ' w li 5,/,'#')     
C
C
C     create second sheet: fourier transforms of wave-functions
C
      write(1,118)
 118  format(22Hset output 'fwfct.pbm')
      write(1,119)
 119  format(36Hset title 'Wave Function Transforms')     
      write(1,120)
 120  format(19Hset format y '%.2f',/,'set xtics 0,5',/,
     1'set yrange [-1.5:1.5]',/,'set xrange [0:10]',/,
     223Hset xlabel 'q [1/a.u.]',/,19Hset ylabel 'R(q)' 1)
      write(1,121)
 121  format('set nolabel')
      write(1,124)
 124  format('plot',10H'fwfct.gp',' using 1:2 with lines 4,',1H\)          
      do 15,i=3,ifw+1
         write(1,125)i
 125     format(10H'fwfct.gp',' using 1:',i2,' with lines 4,',1H\)
  15  continue
      write(1,112)
 112  format(15H0 ti ' ' w li 1,/,'#')     
C
C
C     create third sheet: pseudopotentials
C
C     scaling:
C
      iyminv=int(yminv-1.)
      ihelp=-iyminv
      do 500 i=1,5
         if(mod(ihelp,5).ne.0)ihelp=ihelp+1
 500  continue     
      iyminv=-ihelp         
C
      write(1,138)
 138  format(21Hset output 'ppot.pbm')
      write(1,139)
 139  format(28Hset title 'Pseudopotentials')     
      write(1,140)iyminv/2.,iyminv
 140  format('set nolabel',/,'set key 3.5,',f7.1,/,'set yrange [',i4,':
     15]',/,'set xrange [0:4]',/,21Hset xlabel 'r [a.u.]',/,
     224Hset ylabel 'V(r) [Ry]' 1,/,20Hset format xy '%.1f',/,
     3'set xtics 0,1.0')
      write(1,144)
 144  format('plot ',1H\)     
      do 25,i=2,iv+1
         write(1,145)i,labelv(i-1),i-1
 145     format(9H'ppot.gp',' using 1:',i2,5H ti ',a12,1H',
     1          ' w li',i2,',',1H\)
  25  continue
      write(1,147)zion*2.
 147  format('-',f5.2,'/x',37H ti 'Z - ion' w li 5, 0 ti ' ' w li 1)      
      write(1,148)
 148  format('#')
C
C
C     create fourth sheet: fourier transforms of pseudopotentials
C
      write(1,158)
 158  format(22Hset output 'fppot.pbm')
      write(1,160)
 160  format('set nolabel',/,'set key 10,-0.4',/,
     1'set size 2.8/5.,3/3.',/,30Hset title 'Fourier Transforms',
     2/,'set xrange [0:12]',/,'set yrange [-1.5:1.5]',/,23Hset xlabel 'q
     3 [1/a.u.]',/,36Hset ylabel 'q^2/4piZion V(q) [Ry]' 1,/,'set xtics 
     40,5.0')                             
      write(1,144)
      do 35,i=2,ifp+1
         write(1,165)i,labelf(i-1),i-1
 165     format(10H'fppot.gp',' using 1:',i2,5H ti ',a12,1H',
     1          ' w li',i2,',',1H\)
  35  continue
      write(1,112)
C
C
      return
      end
CCCCC     CCCCC     CCCCC     CCCCC     CCCCC     CCCCC     CCCCC     CC
