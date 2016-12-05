      subroutine krige(ix,iy,iz,icut,MAXCTX,MAXCTY,MAXCTZ,MAXKR1,
     + MAXROT,MAXDAT,MXCUT,MAXKR2,MAXNOD,MAXCUT,MAXXYZ,MAXNST,idbg,
     + imbsim,ivtype,ktype,ldbg,mik,nclose,ncnode,nctx,ncty,nctz, 
     + xx,yy,zz,gmean,cmean,atnode,softdat,vr,vra,x,y,z,icnode,
     + iznode,iynode,ixnode,it,nst,aclose,cnodex,cnodey,cnodez,
     + cnodev,cnodet,covtab,thres,aa,c0,cc,cmax,close,beez,r,rr,s,a,
     + rotmat)

      implicit none

      integer ix,iy,iz,icut,MAXCTX,MAXCTY,MAXCTZ,MAXKR1,MAXROT,MAXDAT,
     +        MXCUT,MAXKR2,MAXNOD,MAXCUT,MAXXYZ,MAXNST,idbg,imbsim,
     +        ivtype,ktype,ldbg,mik,nclose,ncnode,nctx,ncty,nctz 
      real    xx,yy,zz,gmean,cmean
      logical atnode(MAXDAT),softdat(MAXKR1)
      real    vr(MAXDAT,MXCUT),vra(MAXKR1),x(MAXDAT),y(MAXDAT),z(MAXDAT)
      integer icnode(MAXNOD),iznode(MAXXYZ),iynode(MAXXYZ),
     +        ixnode(MAXXYZ),it(MAXCUT*MAXNST),nst(MAXCUT),
     +        aclose(MAXKR1)
      real    cnodex(MAXNOD),cnodey(MAXNOD),cnodez(MAXNOD),
     +        cnodev(MAXNOD),cnodet(MAXNOD),
     +        covtab(MAXCTX,MAXCTY,MAXCTZ,MAXCUT),thres(MAXCUT),
     +        aa(MAXCUT*MAXNST),c0(MAXCUT),cc(MAXCUT*MAXNST),
     +        cmax(MAXCUT),close(MAXDAT),beez(MAXCUT)  
      real*8  r(MAXKR1),rr(MAXKR1),s(MAXKR1),a(MAXKR2),
     +        rotmat(MAXROT,3,3) 


      logical krig,somesoft,bothsoft
      integer i,j,index,irot,in,j1,iii,ind,ix1,iy1,iz1,ie,ii,jj,is,ix2,
     +        iy2,iz2,kk,ising,mclose,na,neq,nhd,jnew
      real    x1,y1,z1,x2,y2,z2,sumwt,cov
c
c Size of the kriging system:  Some of the data values may be missing
c which would correspond to a constraint interval.  Note that there
c should not be any missing values if the median approximation is being
c considered.  The variable ``krig'' is used
c to flag whether kriging is to be done or if the previous weights are
c to be used.
c

      somesoft = .false.
      krig     = .true.
      if(mik.eq.1.and.icut.gt.1) krig = .false.
      if(krig) then
            mclose = 0
            do i=1,nclose
                  index     =  int(close(i))
c            if(ix.eq.27.and.iy.eq.1.and.iz.eq.1) then 
c                  print *,'TAG index=',index
c            end if
                  if(.not.atnode(index).and.vr(index,icut).ge.0.0) then
                        mclose = mclose + 1
                        aclose(mclose) = index
c            if(ix.eq.27.and.iy.eq.1.and.iz.eq.1) then 
c                        print *,'TAG aclose=',aclose(mclose)
c            end if
                  endif
            end do
            na  = mclose + ncnode
            neq = na + ktype
c            print *,'na=',na,'neq=',neq
      endif
c
c There are no data yet:
c
      irot   = 1 + (icut-1)*MAXNST
c
c Set up kriging matrices:
c
      in = 0
      j1 = 0
      do 1 j=1,na
            softdat(j) = .false.
c
c Sort out the actual location of point "j"
c
            if(j.le.mclose) then
                  index  = aclose(j)
                  vra(j) = vr(index,icut)
                  x1     = x(index)
                  y1     = y(index)
                  z1     = z(index)
                  if(index.gt.nhd) softdat(j) = .true.
            else
c
c It is a previously simulated node (keep index for table look-up):
c
                  index  = j-mclose
                  x1     = cnodex(index)
                  y1     = cnodey(index)
                  z1     = cnodez(index)
c
c Is this node informed by a hard datum or a soft datum?
c
c            if(ix.eq.27.and.iy.eq.1.and.iz.eq.1) then 
c       print *,'TAG index=',index,
c     + 'cnodet=',cnodet(index),'cnodev=',
c     + cnodev(index),'thres=',thres(icut)
c            end if
                  if(cnodet(index).le.0.5) then
                        if(ivtype.eq.0) then
                           vra(j) = 0.0
                           if(int(cnodev(index)+0.5).eq.
     +                        int(thres(icut)+0.5)) vra(j) = 1.0
                        else
                           vra(j) = 1.0
                           if(cnodev(index).gt.thres(icut)) vra(j) = 0.0
                        end if
                  else
                        iii    = int(cnodet(index)+0.5)
                        vra(j) = vr(iii,icut)
                        softdat(j) = .true.
                  end if
                  ind    = icnode(index)
                  ix1    = ix + (int(ixnode(ind))-nctx-1)
                  iy1    = iy + (int(iynode(ind))-ncty-1)
                  iz1    = iz + (int(iznode(ind))-nctz-1)
            endif
c
c Only set up the matrix and the RHS if kriging:
c
            if(krig) then
               do 2 i=1,j
c
c Sort out the actual location of point "i"
c
                  if(i.le.mclose) then
                        index  = aclose(i)
                        x2     = x(index)
                        y2     = y(index)
                        z2     = z(index)
                  else
c
c It is a previously simulated node (keep index for table look-up):
c
                        index  = i-mclose
                        x2     = cnodex(index)
                        y2     = cnodey(index)
                        z2     = cnodez(index)
                        ind    = icnode(index)
                        ix2    = ix + (int(ixnode(ind))-nctx-1)
                        iy2    = iy + (int(iynode(ind))-ncty-1)
                        iz2    = iz + (int(iznode(ind))-nctz-1)
                  endif
c
c Now, get the covariance value:
c
                  in = in + 1
c
c Decide whether or not to use the covariance look-up table:
c
                  if(j.le.mclose.or.i.le.mclose) then
                        call cova3(x1,y1,z1,x2,y2,z2,icut,nst,MAXNST,
     +                       c0,it,cc,aa,irot,MAXROT,rotmat,cmax,cov)
                        a(in) = dble(cov)
                  else
c
c Try to use the covariance look-up (if the distance is in range):
c
                        ii = nctx + 1 + (ix1 - ix2)
                        jj = ncty + 1 + (iy1 - iy2)
                        kk = nctz + 1 + (iz1 - iz2)
c                        print *,'ii=',ii,'jj=',jj,'kk=',kk
c                        print *,'x1=',x1,'y1=',y1,'z1=',z1
c                        print *,'x2=',x2,'y2=',y2,'z2=',z2
c                        print *,'ix1=',ix1,'iy1=',iy1,'iz1=',iz1
c                        print *,'ix2=',ix2,'iy2=',iy2,'iz2=',iz2
                        if(ii.lt.1.or.ii.gt.MAXCTX.or.
     +                     jj.lt.1.or.jj.gt.MAXCTY.or.
     +                     kk.lt.1.or.kk.gt.MAXCTZ) then
                              call cova3(x1,y1,z1,x2,y2,z2,icut,nst,
     +                                   MAXNST,c0,it,cc,aa,irot,MAXROT,
     +                                   rotmat,cmax,cov)
                              a(in) = dble(cov)
                        else
                              a(in) = dble(covtab(ii,jj,kk,icut))
                        endif
                  endif
 2          continue
c
c Get the RHS value (possibly with covariance look-up table):
c
            if(j.le.mclose) then
                  call cova3(xx,yy,zz,x1,y1,z1,icut,nst,MAXNST,
     +                       c0,it,cc,aa,irot,MAXROT,rotmat,cmax,cov)
                  r(j) = dble(cov)
            else
c
c Try to use the covariance look-up (if the distance is in range):
c
                  ii = nctx + 1 + (ix - ix1)
                  jj = ncty + 1 + (iy - iy1)
                  kk = nctz + 1 + (iz - iz1)
                  if(ii.lt.1.or.ii.gt.MAXCTX.or.
     +               jj.lt.1.or.jj.gt.MAXCTY.or.
     +               kk.lt.1.or.kk.gt.MAXCTZ) then
                        call cova3(xx,yy,zz,x1,y1,z1,icut,nst,MAXNST,
     +                       c0,it,cc,aa,irot,MAXROT,rotmat,cmax,cov)
                        r(j) = dble(cov)
                  else
                        r(j) = dble(covtab(ii,jj,kk,icut))
                  endif
            endif
            rr(j) = r(j)
c
c End ``if'' block (true if kriging)
c
         endif
c
c End loop over all of the nearby data
c
      if(softdat(j)) somesoft = .true.
 1    continue

c
c If we are doing Markov-Bayes are there are soft data we need to
c correct some of the covariance values in the kriging matrix:
c
      if(imbsim.eq.1.and.somesoft) then
            in = 0
            do j=1,na
                  do i=1,j
                        in = in + 1
                        bothsoft = .false.
                        if(softdat(j).and.softdat(i)) bothsoft = .true.
c
c Correct for soft-soft covariance or soft-hard covariance:
c
                        if(bothsoft) then
                              a(in) = a(in)*dble(beez(icut))
                              if(i.ne.j) a(in) = a(in)*dble(beez(icut))
                        else
                              if(softdat(j).or.softdat(i))
     +                        a(in) = a(in)*dble(beez(icut))
                        end if
                  end do
c
c Correct the right hand side for soft-hard covariance:
c
                  if(softdat(j)) then
                        r(j)  = r(j)*dble(beez(icut))
                        rr(j) = r(j)
                  end if
            end do
      end if
c
c Addition of OK constraint:
c
      if(krig.and.ktype.eq.1) then
            do i=1,na
                  in    = in + 1
                  a(in) = 1.0
            end do
            in      = in + 1
            a(in)   = 0.0
            r(neq)  = 1.0
            rr(neq) = 1.0
      endif
c
c Write out the kriging Matrix if Seriously Debugging:
c
c      if(krig.and.idbg.ge.4) then


c      if(krig.and..false..and.ix.eq.329.and.iy.eq.1.and.iz.eq.1) then
cc            write(ldbg,101) ix,iy,iz
c            write(*,*) ix,iy,iz
c            is = 1
c            do i=1,neq
c                  ie = is + i - 1
cc                  write(ldbg,102) i,r(i),(a(j),j=is,ie)
cc                  write(*,*) i,r(i),(a(j),j=is,ie)
c                  write(*,*) 'r(',i,')=',r(i)
cc                  is = is + i
c            end do
c            do i = 1,(neq*(neq+1)/2)
c                  write(*,*) 'a=',a(i)
c            end do
c 101        format(/,'Kriging Matrices for Node: ',3i4,
c     +               ' RHS first')
c 102        format('    r(',i2,') =',f7.4,'  a= ',9(10f7.4))
c      endif


c
c Solve the Kriging System:
c
      if(krig) then
            if(neq.eq.1.and.ktype.eq.0) then
                  s(1)  = r(1) / a(1)
                  ising = 0
            else
c                  call ksol(1,neq,1,a,r,s,ising)
c                  if(neq.ne.24)then
                  call ksolopt(1,neq,1,a,r,s,ising)
c                  else
c                  call ksolunroll24(1,neq,1,a,r,s,ising)
c                  end if
            endif
      endif
c
c Write a warning if the matrix is singular:
c
      if(ising.ne.0) then
            if(idbg.ge.1) then
                  write(ldbg,*) 'WARNING SISIM: singular matrix'
                  write(ldbg,*) '              for block',ix,iy,iz
            endif
            cmean  = 0.0
            return
      endif
c
c Write out the kriging Matrix if Seriously Debugging:
c
c      if(krig.and.idbg.ge.4) then
c            do i=1,na
c                  write(ldbg,140) i,s(i)
c            end do
c 140        format(' Kriging weight for data: ',i4,' = ',f8.4)
c      endif
c
c Compute the estimate, the sum of weights, correct for SK, and return:
c
      cmean = 0.0
      sumwt = 0.0
      do i=1,na
c            if(ix.eq.27.and.iy.eq.1.and.iz.eq.1) then 
c            print *,'TAG X = ',s(i),' VRA = ',vra(i)
c            end if
            cmean = cmean + real(s(i)) * vra(i)
            sumwt = sumwt + real(s(i))
      end do
      if(ktype.eq.0) cmean = cmean + (1.0-sumwt)*gmean
      return
      end


