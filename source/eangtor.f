c
c
c     ##########################################################
c     ##  COPYRIGHT (C) 2014 by Chao Lv & Jay William Ponder  ##
c     ##                  All Rights Reserved                 ##
c     ##########################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine eangtor  --  angle-torsion cross term energy  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "eangtor" calculates the angle-torsion potential energy
c
c
      subroutine eangtor
      implicit none
      include 'sizes.i'
      include 'angtor.i'
      include 'atoms.i'
      include 'angle.i'
      include 'angpot.i'
      include 'bound.i'
      include 'energi.i'
      include 'group.i'
      include 'math.i'
      include 'torpot.i'
      include 'tors.i'
      include 'usage.i'
      integer i,k,iangtor
      integer ia,ib,ic,id
      real*8 e,rcb,dr,fgrp
      real*8 rt2,ru2,rtru
      real*8 rba2,rcb2,rdc2
      real*8 xt,yt,zt
      real*8 xu,yu,zu
      real*8 xtu,ytu,ztu
      real*8 v1,v2,v3
      real*8 c1,c2,c3
      real*8 s1,s2,s3
      real*8 sine,cosine
      real*8 sine2,cosine2
      real*8 sine3,cosine3
      real*8 phi1,phi2,phi3
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xba,yba,zba
      real*8 xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 angle,dot
      real*8 cosang,force
      logical proceed      
c
c
c     zero out the energy due to extra potential terms
c
      eat = 0.0d0
c
c     calculate the stretch-torsion interaction energy term
c
      do iangtor = 1, nangtor
         i = iat(1,iangtor)
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ib,ic,id,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     compute the value of the torsional angle
c
         if (proceed) then
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xib = x(ib)
            yib = y(ib)
            zib = z(ib)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
            xid = x(id)
            yid = y(id)
            zid = z(id)
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            rba2 = xba*xba + yba*yba + zba*zba
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
            rdc2 = xdc*xdc + ydc*ydc + zdc*zdc
            if (use_polymer) then
               call image (xba,yba,zba)
               call image (xcb,ycb,zcb)
               call image (xdc,ydc,zdc)
            end if
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(rcb2)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
c
c     compute the multiple angle trigonometry and the phase terms
c
               c1 = tors1(3,i)
               s1 = tors1(4,i)
               c2 = tors2(3,i)
               s2 = tors2(4,i)
               c3 = tors3(3,i)
               s3 = tors3(4,i)
               cosine2 = cosine*cosine - sine*sine
               sine2 = 2.0d0 * cosine * sine
               cosine3 = cosine*cosine2 - sine*sine2
               sine3 = cosine*sine2 + sine*cosine2
               phi1 = 1.0d0 + (cosine*c1 + sine*s1)
               phi2 = 1.0d0 + (cosine2*c2 + sine2*s2)
               phi3 = 1.0d0 + (cosine3*c3 + sine3*s3)
c
c     set the angle-torsion parameters for the first angle
c
               v1 = kant(1,iangtor)
               v2 = kant(2,iangtor)
               v3 = kant(3,iangtor)
c
c     calculate the angle-torsion energy for the first angle 
c
               k = iat(2,iangtor)
               dot = xba*xcb + yba*ycb + zba*zcb
               cosang = -dot / sqrt(rba2*rcb2)
               angle = radian * acos(cosang)
               dr = angle - anat(k)
               force = ak(k)
               e = -2.0d0 * angunit * force * dr
     &                * (v1*phi1 + v2*phi2 + v3*phi3)
c
c     set the angle-torsion parameters for the second angle
c
               v1 = kant(4,iangtor)
               v2 = kant(5,iangtor)
               v3 = kant(6,iangtor)
c
c     calculate the angle-torsion energy for the second angle
c
               k = iat(3,iangtor)
               dot = xcb*xdc + ycb*ydc + zcb*zdc
               cosang = -dot / sqrt(rcb2*rdc2)
               angle = radian * acos(cosang)
               dr = angle - anat(k)
               force = ak(k)
               e = e - 2.0d0 * angunit * force * dr
     &                * (v1*phi1 + v2*phi2 + v3*phi3)
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the total bend-torsion energy
c
               eat = eat + e
            end if
         end if
      end do
      return
      end