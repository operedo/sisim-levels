c-----------------------------------------------------------------------
c
c               Search for nearby Simulated Grid nodes
c               **************************************
c
c The idea is to spiral away from the node being simulated and note all
c the nearby nodes that have been simulated.
c
c
c
c INPUT VARIABLES:
c
c   ix,iy,iz        index of the point currently being simulated
c   sim(nx,ny,nz)   the simulation so far
c   nodmax          the maximum number of nodes that we want
c   nlooku          the number of nodes in the look up table
c   i[x,y,z]node    the relative indices of those nodes.
c   [x,y,z]mn       the origin of the global grid netwrok
c   [x,y,z]siz      the spacing of the grid nodes.
c
c
c
c OUTPUT VARIABLES:
c
c   ncnode          the number of close nodes
c   icnode()        the number in the look up table
c   cnode[x,y,z]()  the location of the nodes
c   cnodev()        the values at the nodes
c
c
c
c-----------------------------------------------------------------------

c      use geostat
c      include  'sisim.inc'

      subroutine srchnd(ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,ncnode,
     +        maxsec,nctx,ncty,nctz,nlooku,nodmax,nx,ny,nz,UNEST,
     +        xmn,ymn,zmn,xsiz,ysiz,zsiz,icnode,ixnode,iynode,iznode,
     +        cnodex,cnodey,cnodez,cnodev,cnodet,tmp,sim,cnodeid)


      implicit none

      integer ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,ncnode,maxsec,
     +        nctx,ncty,nctz,nlooku,nodmax,nx,ny,nz 
      real    UNEST,xmn,ymn,zmn,xsiz,ysiz,zsiz
      integer icnode(MAXNOD),ixnode(MAXXYZ),iynode(MAXXYZ),
     +        iznode(MAXXYZ)
      real    cnodex(MAXNOD),cnodey(MAXNOD),cnodez(MAXNOD),
     + cnodev(MAXNOD),cnodet(MAXNOD),tmp(MAXORD),sim(MXYZ)
      integer cnodeid(MAXNOD)

      integer i,j,k,il,indexloc,ncsec
c
c Consider all the nearby nodes until enough have been found:
c
      ncnode = 0
      ncsec  = 0


      do il=1,nlooku

            if(ncnode.eq.nodmax) return
            i = ix + (int(ixnode(il))-nctx-1)
            if(i.lt.1.or.i.gt.nx) go to 1
            j = iy + (int(iynode(il))-ncty-1)
            if(j.lt.1.or.j.gt.ny) go to 1
            k = iz + (int(iznode(il))-nctz-1)
            if(k.lt.1.or.k.gt.nz) go to 1


c
c Check this potentially informed grid node:
c
            indexloc = (k-1)*nx*ny + (j-1)*nx + i

            if(sim(indexloc).gt.UNEST) then
c            if(sim(indexloc).gt.UNEST.or.
c     + tmp(indexloc).gt.0.5) then
c                  if(sim(indexloc).le.UNEST.and.
c     + tmp(indexloc).gt.0.5) then
c                        ncsec  = ncsec + 1
c                        if(ncsec.gt.maxsec) go to 1
c                  end if

                  ncnode         = ncnode + 1
                  cnodeid(ncnode) = indexloc
                  icnode(ncnode) = il
                  cnodex(ncnode) = xmn + real(i-1)*xsiz
                  cnodey(ncnode) = ymn + real(j-1)*ysiz
                  cnodez(ncnode) = zmn + real(k-1)*zsiz
                  cnodev(ncnode) = sim(indexloc)
                  cnodet(ncnode) = tmp(indexloc)

            endif
 1          continue
      end do
c
c Return to calling program:
c
      return
      end

      subroutine srchndPush(ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,ncnode,
     +        maxsec,nctx,ncty,nctz,nlooku,nodmax,nx,ny,nz,UNEST,
     +        xmn,ymn,zmn,xsiz,ysiz,zsiz,icnode,ixnode,iynode,iznode,
     +        cnodex,cnodey,cnodez,cnodev,cnodet,tmp,sim,cnodeid)


      implicit none

      integer ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,ncnode,maxsec,
     +        nctx,ncty,nctz,nlooku,nodmax,nx,ny,nz 
      real    UNEST,xmn,ymn,zmn,xsiz,ysiz,zsiz
      integer icnode(MAXNOD),ixnode(MAXXYZ),iynode(MAXXYZ),
     +        iznode(MAXXYZ)
      real    cnodex(MAXNOD),cnodey(MAXNOD),cnodez(MAXNOD),
     + cnodev(MAXNOD),cnodet(MAXNOD),tmp(MAXORD),sim(MXYZ)
      integer cnodeid(MAXNOD)

      integer i,j,k,il,indexloc,ncsec
c
c Consider all the nearby nodes until enough have been found:
c
      ncnode = 0
      ncsec  = 0


      do 2 il=1,nlooku

            if(ncnode.eq.nodmax) return
            i = ix + (int(ixnode(il))-nctx-1)
            if(i.lt.1.or.i.gt.nx) go to 2
            j = iy + (int(iynode(il))-ncty-1)
            if(j.lt.1.or.j.gt.ny) go to 2
            k = iz + (int(iznode(il))-nctz-1)
            if(k.lt.1.or.k.gt.nz) go to 2

            indexloc = (k-1)*nx*ny + (j-1)*nx + i

            if(sim(indexloc).gt.UNEST) then
                  ncnode         = ncnode + 1
                  cnodeid(ncnode) = indexloc
                  icnode(ncnode) = il
            endif
 2    continue

c      do il=1,nlooku,2
c
c            if(ncnode.eq.nodmax) return
c            i = ix + (int(ixnode(il))-nctx-1)
c            if(i.lt.1.or.i.gt.nx) go to 1
c            j = iy + (int(iynode(il))-ncty-1)
c            if(j.lt.1.or.j.gt.ny) go to 1
c            k = iz + (int(iznode(il))-nctz-1)
c            if(k.lt.1.or.k.gt.nz) go to 1
c
c            indexloc = (k-1)*nx*ny + (j-1)*nx + i
c
c            if(sim(indexloc).gt.UNEST) then
c                  ncnode         = ncnode + 1
c                  cnodeid(ncnode) = indexloc
c                  icnode(ncnode) = il
c            endif
c 1          continue
c
c            if(ncnode.eq.nodmax) return
c            i = ix + (int(ixnode(il+1))-nctx-1)
c            if(i.lt.1.or.i.gt.nx) go to 2
c            j = iy + (int(iynode(il+1))-ncty-1)
c            if(j.lt.1.or.j.gt.ny) go to 2
c            k = iz + (int(iznode(il+1))-nctz-1)
c            if(k.lt.1.or.k.gt.nz) go to 2
c
c            indexloc = (k-1)*nx*ny + (j-1)*nx + i
c
c            if(sim(indexloc).gt.UNEST) then
c                  ncnode         = ncnode + 1
c                  cnodeid(ncnode) = indexloc
c                  icnode(ncnode) = il+1
c            endif
c 2          continue
c
c      end do
c
c Return to calling program:
c
      return
      end

      subroutine srchndPushOpt(ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,
     +        ncnode,
     +        maxsec,nctx,ncty,nctz,nlooku,nodmax,nx,ny,nz,UNEST,
     +        xmn,ymn,zmn,xsiz,ysiz,zsiz,icnode,ixnode,iynode,iznode,
     +        cnodex,cnodey,cnodez,cnodev,cnodet,tmp,sim,cnodeid)


      implicit none

      integer ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,ncnode,maxsec,
     +        nctx,ncty,nctz,nlooku,nodmax,nx,ny,nz 
      real    UNEST,xmn,ymn,zmn,xsiz,ysiz,zsiz
      integer icnode(MAXNOD),ixnode(MAXXYZ),iynode(MAXXYZ),
     +        iznode(MAXXYZ)
      real    cnodex(MAXNOD),cnodey(MAXNOD),cnodez(MAXNOD),
     + cnodev(MAXNOD),cnodet(MAXNOD),tmp(MAXORD),sim(MXYZ)
      integer cnodeid(MAXNOD)

      integer i,j,k,il,indexloc,ncsec,ind
c
c Consider all the nearby nodes until enough have been found:
c

      ncnode = 0
      ncsec  = 0
      do 2 il=1,nodmax-1
            i = ix + (int(ixnode(il))-nctx-1)
            j = iy + (int(iynode(il))-ncty-1)
            k = iz + (int(iznode(il))-nctz-1)
            if(i.lt. 1.or.j.lt. 1.or.k.lt. 1) go to 2
            if(i.gt.nx.or.j.gt.ny.or.k.gt.nz) go to 2
            ind = i + (j-1)*nx + (k-1)*nx*ny
            if(sim(ind).gt.UNEST) then
                  ncnode = ncnode + 1
                  cnodeid(ncnode) = ind
                  icnode(ncnode) = il
            endif
 2    continue
      do 3 il=nodmax,nlooku
            i = ix + (int(ixnode(il))-nctx-1)
            j = iy + (int(iynode(il))-ncty-1)
            k = iz + (int(iznode(il))-nctz-1)
            if(i.lt. 1.or.j.lt. 1.or.k.lt. 1) go to 3
            if(i.gt.nx.or.j.gt.ny.or.k.gt.nz) go to 3
            ind = i + (j-1)*nx + (k-1)*nx*ny
            if(sim(ind).gt.UNEST) then
                  ncnode = ncnode + 1
                  cnodeid(ncnode) = ind
                  icnode(ncnode) = il
                  if(ncnode.eq.nodmax) return
            endif
 3    continue
      return
c
c Return to calling program:
c
      end


      subroutine srchndPop(ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,ncnode,
     +        maxsec,nctx,ncty,nctz,nlooku,nodmax,nx,ny,nz,UNEST,
     +        xmn,ymn,zmn,xsiz,ysiz,zsiz,icnode,ixnode,iynode,iznode,
     +        cnodex,cnodey,cnodez,cnodev,cnodet,tmp,sim,cnodeid)


      implicit none

      integer ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,ncnode,maxsec,
     +        nctx,ncty,nctz,nlooku,nodmax,nx,ny,nz 
      real    UNEST,xmn,ymn,zmn,xsiz,ysiz,zsiz
      integer icnode(MAXNOD),ixnode(MAXXYZ),iynode(MAXXYZ),
     +        iznode(MAXXYZ)
      real    cnodex(MAXNOD),cnodey(MAXNOD),cnodez(MAXNOD),
     + cnodev(MAXNOD),cnodet(MAXNOD),tmp(MAXORD),sim(MXYZ)
      integer cnodeid(MAXNOD)

      integer i,j,k,il,indexloc,ncsec,nxy
c
c Consider all the nearby nodes until enough have been found:
c
      nxy=nx*ny
      do il=1,ncnode

            indexloc = cnodeid(il)
            k = int((indexloc-1)/nxy) + 1
            j = int((indexloc-(k-1)*nxy-1)/nx) + 1
            i = indexloc - (k-1)*nxy - (j-1)*nx

c            ncnode         = ncnode + 1
c            cnodeid(il) = indexloc
c            icnode(il) = il
            cnodex(il) = xmn + real(i-1)*xsiz
            cnodey(il) = ymn + real(j-1)*ysiz
            cnodez(il) = zmn + real(k-1)*zsiz
c            cnodev(il) = sim(indexloc)
            cnodet(il) = tmp(indexloc)

      end do
c
c Return to calling program:
c
      return
      end



