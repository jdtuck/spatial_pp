      subroutine smplxBf(x, y, np, ty1, amu1, amu2, anu1, scls11,
     &    scls22, eps, itmax, itmax1, ipmax, fn, mples, xinit, eps1,
     &    f, iter, nip, ipri, ipflag)
c
      include 'NScluster_f.h'
c
c simplx:  simplex minimization subroutine.
c minmax:  called by subroutine simplx.
c first:             "
c centor:            "
c newsim:            "
c update:            "
c reduce:            "
c epslon:            "
c
c define parameter.
c maxh:    max dimension of vector x.
c
c procedure to use this program.
c 1. save this mail to a file. (save file-name.f)
c 2. delete the comments. (vi file name.f)
c 3. change the define parameter and function.
c
c program starts. ------------------------------------------------------
c
      implicit real * 8 (a-h,o-z)
cc      parameter   (maxh=6, maxh5=maxh+5)
      parameter   (n=5)
c
      common/paramscl/sclmu,sclnu,scla1,scls1,scls2
      common /fnmin/ fmin
cc      common / sizes / tx,ty
cc      common /fname/filea
cc      character*50 filea
      real*8 sclmu,sclnu,scla1,scls1,scls2
cc      integer    n
cc      real*8     xinit(maxh), dist, eps, f
cc      external   funct
      real*8     xinit(n,itmax1), dist, eps, f(itmax1)
      external   bfunct
c
      common/skip/iskip
      dimension  x(np), y(np), rr(np**2)
      dimension  eps1(itmax1)
      dimension  ipri(ipmax), fn(ipmax), mples(ipmax,n)
c
      fmin=1.d10
***************************
cc      open(13, file='TypeB.param')
cc      read(13,*) 
      tx=1.0d0
cc      read(13,*) ty
cc      read(13,*) amu1, amu2, anu, scls1, scls2
      ty=ty1
      scls1=scls11
      scls2=scls22
      sclmu=amu1+amu2
      sclnu=anu1
      scla1=amu1/(amu1+amu2)
cc      read(13,1) filea
cc    1 format(a)
***************************
cc      call input
      iskip=1
      call input(x,y,np,tx,ty,rr,nn)
c
cc      n = 5
c
cc      open(8, FILE='TypeB.simplex.print')
cc      write(8,2) '            -log L            mu          nu ',
cc     &'         a           s1          s2'
    2 format(2a)
cc      close(8)
c
cc      xinit(1) = 1.0d0
cc      xinit(2) = 1.0d0
cc      xinit(3) = 1.0d0
cc      xinit(4) = 1.0d0
cc      xinit(5) = 1.0d0 
      xinit(1,1) = 1.0d0
      xinit(2,1) = 1.0d0
      xinit(3,1) = 1.0d0
      xinit(4,1) = 1.0d0
      xinit(5,1) = 1.0d0 
c
      dist = 0.1d0
cc      eps = 1.0d-3
c
cc      call simplx(xinit, n, funct, dist, eps, f)
      nip = 1
      call simplx(xinit, n, rr, nn, bfunct, dist, eps, f,
     & itmax, itmax1, iter, eps1, ipmax, nip, ipri, fn, mples, ipflag)
      if( (ipflag.eq.1) .or. (ipflag.eq.3) ) nip = nip-1
c
cc      stop
      return
      end
c -------------------------------------------------------------------- c
c
cc      subroutine funct(n,b,fn)
      subroutine bfunct(n,b,fn,rr,nn,nip,jpri,ffn,mples,
     &                       ipmax,ipflag)
c-----------------------------------------------------------------------
c     likelihood function of the inverse power poisson process
c-----------------------------------------------------------------------
      implicit real * 8 (a-h,o-z)
cc      common/datpar/ nn
cc      common/xyod/rr(9234567),th(9234567)
      dimension rr(nn)
      common/ddd/ff,aic
      common/range/rmin,rmax
      common/paramscl/sclmu,sclnu,scla1,scls1,scls2
      common/events/np
      common /fnmin/ fmin
c
      data pi/3.14159265358979d0/
cc      dimension b(5),g(5),h(5)
      dimension b(n),g(n),h(n)
c
      real*8   ffn(ipmax), mples(ipmax,n)
      integer  jpri(ipmax)
c
      pi=3.14159265358979d0
      eps=0.1d-5
      amu=b(1)**2 * sclmu
      anu=b(2)**2 * sclnu
      a1=b(3)**2 * scla1
      s1=b(4)**2 * scls1
      s2=b(5)**2 * scls2
      se=s2-eps
      f1=0.0
c-----
      amuanu=amu*anu
      s124=4*(s1**2)
      s224=4*(s2**2)
c-----
      do 30 i=1,nn
cc      if(rr(i).le.rmin.or.rr(i).ge.rmax) go to 30
cc      ramdai = (amu*anu)
cc     &        +(
cc     &          a1*(exp(-((rr(i))**2)/(4*((s1)**2))))/((s1)**2)
cc     &         +(1-a1)*(exp(-((rr(i))**2)/(4*((s2)**2))))/((s2)**2)
cc     &          )
cc     &           *(anu/(4*pi))
      ramdai = amuanu
     &        +(
     &          a1*(exp(-((rr(i))**2)/s124))/((s1)**2)
     &         +(1-a1)*(exp(-((rr(i))**2)/s224))/((s2)**2)
     &          )
     &           *(anu/(4*pi))
      if(ramdai.le.0.0) go to 190
      f1=f1+log(ramdai)
   30 continue
c
      fs1=pi*(rmax**2)*(amu*anu)
c
cc      ainteg1 = a1*anu*(1-(exp(-(rmax**2)/(4*((s1)**2)))))
cc      ainteg2 = (1-a1)*anu*(1-(exp(-(rmax**2)/(4*((s2)**2)))))       
      ainteg1 = a1*anu*(1-(exp(-(rmax**2)/s124)))
      ainteg2 = (1-a1)*anu*(1-(exp(-(rmax**2)/s224)))       
      fn=f1-(fs1+ainteg1+ainteg2)*np
      fn=-fn
      ff=fn
      if(fmin.gt.fn) then
      fmin=fn
      ipri=1
      else
      ipri=0
      end if
cc      open(8, FILE='TypeB.simplex.print',position='APPEND')
cc      if(ipri.eq.0) write(8,2) 'testfn =',fn,amu,anu,a1,s1,s2
cc      if(ipri.eq.1) write(8,2) 'update =',fn,amu,anu,a1,s1,s2
cc      close(8)
cc      write(6,*) 'h(3)=', h(3)
      ffn(nip) = fn
      mples(nip,1) = amu
      mples(nip,2) = anu
      mples(nip,3) = a1
      mples(nip,4) = s1
      mples(nip,5) = s2
      if((ipflag.eq.0) .or. (ipflag.eq.2)) return
      if(ipri.eq.0) jpri(nip) = -1
      if(ipri.eq.1) jpri(nip) = 1
      nip = nip+1
      return
    3 format(1h ,110x,d18.10)
    2 format(1h , a, d18.10,5d12.5)
    1 format(1h ,7d18.10)
c
  190 continue
      fn=1.0d30
      return
      end
