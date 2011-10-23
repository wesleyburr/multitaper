C$$$    The multitaper R package
C$$$    Multitaper and spectral analysis package for R
C$$$    Copyright (C) 2011 Karim J. Rahim David J. Thomson 

C$$$    This file is part of the multitaper package for R.

C$$$    The multitaper package is free software: you can redistribute it and
C$$$    or modify
C$$$    it under the terms of the GNU General Public License as published by
C$$$    the Free Software Foundation, either version 2 of the License, or
C$$$    any later version.

C$$$    The multitaper package is distributed in the hope that it will be 
C$$$    useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
C$$$    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C$$$    GNU General Public License for more details.

C$$$    You should have received a copy of the GNU General Public License
C$$$    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

C$$$    If you wish to report bugs please contact the author. 
C$$$    karim.rahim@gmail.com
C$$$    112 Jeffery Hall, Queen's University, Kingston Ontario
C$$$    Canada, K7L 3N6


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc This files contains modified djt multitaper files originally
cc from libraries written at Bell Labs by David Thomson.


c**************************************************************************
c     mw2wta multitaper using weights

      subroutine mw2wta(sa,wt,nfreq,nord,s,ev,evp
     1     ,dofs,dofav,var,dt,tol, maxadit
     1     , mxiter, aviter)

      implicit none
      integer nfreq, nord, maxadit, mxiter,n, niter, k
      double precision  sa(nfreq,nord),wt(nfreq,nord),dofs(nfreq)
     1     ,s(nfreq),ev(nord),evp(nord),var,dt,tol, aviter, avewt 
     1     ,wtmin, dofmin, sewn, sbar, dk2, sum, cwt, dofav, wmin
     1     ,dk2l


c     Generate Weights
       mxiter = 0
       aviter = 0.d0
       avewt = 0.d0
       wtmin = 1.d0
       dofmin = dble(2*nord)
       cwt = 0.d0
       wmin = 0.d0
       
c     Equivalent white noise level for bias calculation
       sewn = var*dt
      do 265 n=1,nfreq
c        start at estimate based on two best eigenvalues
         sbar = ( sa(n,1) + sa(n,2) )/2.d0
         dk2 = 1.d0
c     iterate
         do 262 niter=1, maxadit
            sum = 0.d0
            cwt = 0.d0
            wmin = 1.d0
            dk2l = dk2
            do 250 k = 1,nord
               dk2 = 
     1              ( ev(k)*sbar/( ev(k)*sbar + evp(k)*sewn ) )**2
               wt(n,k) = dk2
               sum = sum + sa(n,k)*dk2
               wmin = dmin1(wmin,dk2)
               cwt = cwt + dk2
 250        continue
            sbar = sum/cwt
            if(dabs((dk2-dk2l)/(dk2+dk2l)).le.tol) exit
 262  continue
         mxiter = max0(mxiter,niter)
         aviter = aviter + niter
         avewt = avewt + cwt
         wtmin = dmin1(wtmin,wmin)
         dofs(n) = 2.d0*cwt
         dofmin = dmin1(dofmin,dofs(n))
         s(n) = sbar
         aviter = aviter/dble(nfreq)
 265  continue
      dofav = 2.d0*avewt/dble(nfreq)
      
      end subroutine
c*****end mw2wta

cc**********************************************************************
c     multiwindow jacknifed.

c     Multi-Window Weighting, Jackknifed
      subroutine mw2jkw(sa,wt,nfreq,nord,s,ev,evp
     1     ,dofs,dofav,var,dt,tol, sjk,varjk,bcjk,wjk,cwjk,vwj
     1     ,maxadit, mxiter)
      
      implicit none
      integer nfreq, nord, mxiter, n1, n2,j, n, niter, ks, maxadit
     1     ,k
      double precision sa(nfreq,nord),wt(nfreq,nord),dofs(nfreq)
     1     , s(nfreq)
     1     ,ev(nord),evp(nord),sjk(nord+2),varjk(nfreq),bcjk(nfreq)
     2     ,wjk(nord,nord+2),cwjk(nord+2),vwj(nord)
     1     ,total, avewt, wtmin, dofmin, bcor
     1     ,fnord, vnrm, sewn,var,dt, dofav, varcwt, sbar,sum
     1     ,wmin, slast, tol

c     Generate Weights

      mxiter = 0
      total = 0.d0
      avewt = 0.d0
      wtmin = 1.d0
      niter = 0
      sbar = 0.d0
      wmin =0.d0
      dofmin = dble(2*nord)
      bcor = dble(nord-1)
      fnord = nord
      vnrm = dble(nord-1)/fnord
      n1 = nord + 1
      n2 = nord + 2
c     Equivalent white noise level for bias calculation
      sewn = var*dt
c
      do 365 n=1,nfreq
c     iterate
         do 433 ks = 1, nord+1
c     start at estimate based on two best eigenvalues
            sbar = ( sa(n,1) + sa(n,2) )/2.
            do 362 niter=1, maxadit
               sum = 0.
               cwjk(ks) = 0.d0
               wmin = 1.d0
               slast = sbar
               do 350 k = 1,nord
                  if(k.eq.ks) go to 350
                  wjk(k,ks) = ( ev(k)*sbar/
     1                 ( ev(k)*sbar + evp(k)*sewn ) )**2
                  sum = sum + sa(n,k)*wjk(k,ks)
                  wmin = dmin1(wmin,wjk(k,ks))
                  cwjk(ks) = cwjk(ks) + wjk(k,ks)
 350           continue
               sbar = sum/cwjk(ks)
               sjk(ks) = dlog(sbar)
               if(dabs((sbar-slast)/(sbar+slast)).le.tol) exit
 362        continue
 433     continue

c     Jackknife mean, variance of Log S
         sjk(n2) = 0.d0
         cwjk(n2) = 0.d0
         do 490 k = 1, nord
            wjk(k,n2) = 0.d0
 490     continue
         do 500 k = 1, nord
            cwjk(n2) = cwjk(n2) + cwjk(k)
            sjk(n2) = sjk(n2) + sjk(k)
            do 510 j = 1, nord
               wjk(j,n2) = wjk(j,n2) + wjk(j,k)
 510        continue
 500     continue
         sjk(n2) = sjk(n2)/fnord
         cwjk(n2) = cwjk(n2)/fnord
         do 610 j = 1, nord
            vwj(j) = 0.d0
            wjk(j,n2) = wjk(j,n2)/fnord
            wt(n,j) = wjk(j,n2)
 610     continue
         
c     Jackknife Bias Estimate (Log S )
         bcjk(n) = bcor*(sjk(n2) - sjk(n1))
         
c     Variance Estimate
         varjk(n) = 0.d0
         varcwt = 0.d0
         do 550 k = 1, nord
            varjk(n) = varjk(n) + (sjk(k)-sjk(n2))**2
            varcwt = varcwt + (cwjk(k)-cwjk(n2))**2
            do 560 j = 1, nord
               vwj(j) = vwj(j) + (wjk(j,k)-wjk(j,n2))**2
 560        continue
 550     continue
         
         varjk(n) = varjk(n)*vnrm
         mxiter = max0(mxiter,niter)
         total = total + niter
         avewt = avewt + cwjk(n1)
         wtmin = dmin1(wtmin,wmin)
         dofs(n) = 2.d0*cwjk(n1)
         dofmin = dmin1(dofmin,dofs(n))
         s(n) = sbar
 365  continue
      dofav = 2.d0*avewt/float(nfreq)
      
      end subroutine
c ****** end mw2jkw


c     Multi-Window Average Estimation
      subroutine mweave(x,dw,swz,ndata,nord,ssqswz,cntr,dt
     1     ,spz, varc)
      implicit none
      integer ndata, nord, n, k, nnx
      double precision x(ndata),dw(ndata,nord),swz(nord)
     1     ,sm(nord),sum,spz, zero8, dt
     1     ,ssqswz, cntr, varc
      data  zero8/0.d+00/

c     no need for a max of 9 
      nnx = nord

      call setdp(nnx,zero8,sm)
      do 100 k = 1, nnx
         do 110 n = 1, ndata
            sm(k) = sm(k) + dw(n,k)*x(n)
 110     continue
 100  continue
      
      sum = zero8
      spz = zero8
      do 300 k = 1, nnx, 2
         sum = sum + swz(k)*sm(k)
 300  continue
      sum = sum/ssqswz
      do 500 k = 1, nnx
         spz = spz + (sm(k) - sum*swz(k))**2
 500  continue
      spz = spz/dble(nnx)
      varc = spz/(dt*dble(ndata))
      cntr = sum
      end subroutine   
c ******** end    mweave   

c     Set Real*8 array
      subroutine setdp(npts,val,x)
      implicit none
      integer npts, n
      double precision val, x(npts)
      do 100 n = 1, npts
         x(n) = val
 100  continue

      end subroutine
c ******* end setdp

c ************************** helper functions used in coherence calculation
c djt/ts/ adstoa.f  Add Scalar to Array
      subroutine adstoa(x,y,ndata,xinc)
      implicit none
      integer n, ndata
      double precision  x(ndata), y(ndata), xinc
c     djt/ts/tsu1 -2- add scalar to array
      do 3400 n = 1,ndata
         y(n) = x(n)+xinc
 3400 continue
      
      end subroutine
c ********** end adstoa

c djt/ts/sphsed.f   Basic Phase Unwrapping Routine, Degrees
      subroutine sphsed(ph,nfreq)
      implicit none
      integer nfreq, n
      double precision ph(nfreq), q, pinc,d, t
      
      q=0.d0
      pinc=0.d0
      do 2100 n=1,nfreq
         t=ph(n)
         d=q-t
         q=t
         if(dabs(d).gt.180.d0) pinc=pinc+dsign(360.d0,d)
         ph(n)=t+pinc
 2100 continue
      
      end subroutine
c ****** end       sphsed


c*********************************************************************
cc calculated coherence estimates
      
      subroutine  jkcoh1(cft1, cft2, nord, blklof, blkhif
     1     ,fr, tau, phcorr, NTmsc, NTvar
     1     ,msc, ph, phvar, s1, s2, jkmsc, TRmsc, bias
     1     ,cx)
      
      implicit none
      integer n1, n2, ks, nav, phcorr, blklof, blkhif
     1     ,k, kc, n, nord, nfreqs
      double precision  fr(blklof:blkhif), tau
     1     ,ph(blklof:blkhif), NTmsc(blklof:blkhif),s1(nord+2)
     1     ,s2(nord+2)
     1     ,jkmsc(nord+2),TRmsc(nord+2),bias(blklof:blkhif)
     4     ,phvar(blklof:blkhif),NTvar(blklof:blkhif), cdabs2, phsed
     3     ,trnrm, fnavm, varc, RtoD, RtoD2, msc(blklof:blkhif) 
     1     ,C2toF,  xx, FtoMSC, fnav, xsm2, ff, zpref, d1mach
     1     , dphse
      double complex cft1(blklof:blkhif, nord)
     1     ,cft2(blklof:blkhif, nord)
     1     ,cx(nord+2), zz
      logical phzref
c     


      cdabs2(zz) = dreal(zz)**2 + dimag(zz)**2
      phsed(zz) = RtoD*datan2(dimag(zz),dreal(zz))
c              Transforms from MSC to f, inverse
      C2toF(xx)  = trnrm*dlog((1.+dsqrt(xx))/(1.-dsqrt(xx)))/2.
      FtoMSC(ff) = dtanh(ff/trnrm)**2
c      
      zpref = 0.d0
      nfreqs = blkhif + 1 - blklof
      nav = nord
      n1 = nav + 1
      n2 = nav + 2
      trnrm = dsqrt(dble(2*nav-2))
      fnavm = dble(nav-1)
      fnav = dble(nav)
      varc = fnavm/fnav
      RtoD = 45.d0/datan(1.d0)
      RtoD2 = RtoD**2
 
 
      do 6000 n = blklof, blkhif
         do 1400 ks = 1, nav+1
            kc = 0
            cx(ks) = (0.d0,0.d0)
            s1(ks) = 0.d0
            s2(ks) = 0.d0
            do 1300 k = 1, nav
c     do 300 nb = ns1,ns1+nsav-1
               kc = kc + 1
               if(kc.eq.ks) cycle
               cx(ks) = cx(ks) + cft1(n,k)*dconjg(cft2(n,k))
               s1(ks) = s1(ks) + cdabs2(cft1(n,k))
               s2(ks) = s2(ks) + cdabs2(cft2(n,k))
 1300       continue
            xsm2 = cdabs2(cx(ks))
c     Keep phase in (cos,sin) form
            cx(ks) = cx(ks)/dsqrt(xsm2)
c               MSC
            jkmsc(ks) = xsm2/(s1(ks)*s2(ks))
c     Transform MSC
            TRmsc(ks) = C2toF( jkmsc(ks) )
 1400    continue
c             Bias
         TRmsc(n2) = 0.d0
         cx(n2) = (0.d0,0.d0)
         do 1500 k = 1, nav
            cx(n2) = cx(n2) + cx(k)
            TRmsc(n2) = TRmsc(n2) + TRmsc(k)
 1500    continue
c     Phase and Phase Variance
         cx(n2) = cx(n2)/fnav
         if(cdabs(cx(n2)).le.10.*d1mach(1)) then
            if(n.gt.blklof) then
               ph(n) = ph(n-1)
            else
               ph(n) = 0.d0
            endif
         else
            ph(n) = phsed(cx(n2)) + 360.d0*fr(n)*tau
         endif
         phvar(n) = dble(2*(nav-1))*(1.-cdabs(cx(n2)))*RtoD2
c     Jackknife average of transformed delete-one estimates
         TRmsc(n2) = TRmsc(n2)/fnav
         NTmsc(n) = TRmsc(n1)
         bias(n) = fnavm*( TRmsc(n2) - TRmsc(n1) )
c     J.K. Unbiased Normal Transform to msc
         msc(n) = FtoMSC( NTmsc(n) )
c     Variance
         NTvar(n) = 0.d0
         do 1600 k = 1, nav
            NTvar(n) = NTvar(n) + ( TRmsc(k) - TRmsc(n2) ) **2
 1600    continue
         NTvar(n) = NTvar(n)*varc
 6000 continue
c      cx1(0) = 360.d0
      
c     Keep zero-frequency reference
      phzref = (blklof.le.0).and.(blkhif.ge.0)
      if(phcorr .eq. 1) then 
         if(phzref) zpref = ph(0)
         call sphsed(ph,nfreqs)
         if(phzref) then
            dphse = ph(0) - zpref
            call adstoa(ph,ph,nfreqs,-dphse)
         endif
      endif
          
      end subroutine 
c **** end jkcoh

c*********************************************************************
cc 
cc  quickSineF
cc
cc  Simple non-adaptive (possibly weighted) sine taper multitaper
cc  computation program. Explicitly for calling from within 
cc  the adaptive loop of spec.mtm.sine. The R-specific version 
cc  of this (quickSine) runs quickly on its own; it is the adaptive
cc  loops that need speeding up.
cc
cc  Adapted from Robert Parker's 'psd.f'.
cc
c*********************************************************************
      subroutine quickSineF(nFreqs,nFFT,k,cft,useAdapt,kadapt,spec)

      implicit none

      integer nFreqs,nFFT,k, ks, i, j, i2, j1, j2
      logical useAdapt
      complex*16 cft(nFFT), zz
      real*8 spec(1:nFreqs), ck, wt, kadapt(1:nFreqs)

      do 5 j=1,nFreqs
 5    spec = 0.0d+00

      do 6 i=1,nFreqs
        i2 = 2*(i-1)
        if(useAdapt) then
          ks = int(kadapt(i))
        else
          ks = k
        endif
      
        ck = 1/(real(ks)**2)

        do 7 j=1,ks
          j1 = mod(i2+nFFT-j,nFFT)
          j2 = mod(i2+j,nFFT)
          zz = cft(j1+1) - cft(j2+1)
          wt = 1.0d+00 - ck*(j-1)**2
          spec(i) = spec(i) + (real(zz)**2 + aimag(zz)**2)*wt
 7    continue
      
      spec(i) = spec(i)*(6.0d+00 *real(ks))/(4*(real(ks)**2) + 
     * (3*real(ks)) -1)
 6    continue
      return
      end subroutine

c*********************************************************************
cc 
cc  curbF
cc
cc  Rewrites input vector so all points lie below a piecewise
cc  linear function v(k) + abs(j-k); clips strong peaks, keeps
cc  slopes below 1 in magnitude. Taken from 
cc  Robert Parker's 'psd.f'.
cc
c*********************************************************************
 
      subroutine curbF(n, v)
      implicit none
      integer n, j, k
      real*8 v(n), vloc

      do 1500 j=2, n-1
        if (v(j).lt.v(j+1) .and. v(j).lt.v(j-1)) then
           vloc=v(j)
           do 1200 k=1, n
             v(k)=min(v(k), vloc+abs(j-k))
 1200      continue
        endif
 1500 continue
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc northF
cc
cc Performs quadratically-weighted LS fit to some function 's' by
cc a degree-two polynomial in an orthogonal basis; returns
cc d1 and d2, estimates of 1st and 2nd derivatives at center of record
cc
cc  Taken directly from Robert Parker's 'psd.f'.
cc 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine northF(n, i1, i2, s, ds, dds)
      integer i1, i2, n, el, L, kk, i 
      real*8 gamma, s(n), ds, dds, amid, u1sq, u2sq, dot0, dot1, dot2
     *,   ssq
      L = i2 - i1 + 1
      el=L
      gamma = (el**2 - 1.0)/12.0
      u0sq = el
      u1sq = el*(el**2 - 1.0)/12.0
      u2sq = (el*(el**2 - 1.0)*(el**2- 4.0))/180.0
      amid= 0.5*(el + 1.0)
      dot0=0.0
      dot1=0.0
      dot2=0.0
      ssq=0.0
      do 1100 kk=1, L
        i=kk + i1 - 1
c  Negative or excessive index uses even function assumption
        if (i.le. 0) i=2 - i
        if (i.gt. n) i=2*n - i
        dot0 = dot0 + s(i)
        dot1 = dot1 + (kk - amid) * s(i)
        dot2 = dot2 + ((kk - amid)**2 - gamma)*s(i)
*       ssq = ssq + s(i)**2
 1100 continue
*     c0 = dot0/u0sq
*     c1 = dot1/u1sq
*     c2 = dot2/u2sq
*     resq = abs(ssq - u0sq*c0**2- u1sq*c1**2- u2sq*c2**2)
      ds = dot1/u1sq
      dds = 2.0*dot2/u2sq
c      write(*,*) ds, dds
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc adapt
cc
cc Performs adaptive spectral estimation via sine taper approach.
cc From pilot estimate of spectrum, computes estimates of S'' to be used
cc in Eq. (13) of Riedel & Sidorenko (1995).
cc
cc  Adapted somewhat from Robert Parker's 'psd.f'.
cc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine adapt(ntimes, k, nFreqs, sx, nFFT, cft, df,kopt,fact)
      implicit none
      integer k, ntimes, nFreqs, ispan, iter, nFFT, j
      real*8 sx(nFreqs), kopt(nFreqs), y(nFreqs), dy, ddy, R, ak, phi
     *, sigR, opt(nFreqs), fact, df, c1, c2
      complex*16 cft(nFFT)
      data c1/1.2000/, c2/3.437/

      do 5 j=1,nFreqs
 5      kopt(j) = k

c  Adaptive iteration for MSE spectrum
      do 1600 iter=1, ntimes
c
        do 1100 j=1, nFreqs
           y(j)= log(sx(j))
 1100   continue

c  Estimate K, number of tapers at each freq for MSE spectrum
c  R = S"/S -- use R = Y" + (Y')**2 , Y=ln S.
        do 1200 j=1, nFreqs
c
          ispan = int(kopt(j)*1.4)
          call northF(nFreqs, j-ispan, j+ispan, y, dy, ddy)
          R = (ddy  + dy**2)/df**2
          ak=kopt(j)/real(2*ispan)
          phi=720.0*ak**5*(1.0 - 1.286*ak + 0.476*ak**3 - 0.0909*ak**5)
          sigR= sqrt(phi/real(kopt(j))**5) / df**2
          opt(j)=c2/(df**4 *( R**2 + 1.4*sigR**2) /fact**2)** 0.2
 1200   continue

        call curbF(nFreqs, opt)
        do 1550 j=1, nFreqs
          kopt(j)=max(3.0, opt(j))
 1550   continue
c  Recompute spectrum with optimal variable taper numbers
        call quickSineF(nFreqs,nFFT,1,cft,.true.,kopt,sx)
 1600 continue
      return
      end

