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


      do il=1,nlooku,2

            if(ncnode.eq.nodmax) return
            i = ix + (int(ixnode(il))-nctx-1)
            if(i.lt.1.or.i.gt.nx) go to 1
            j = iy + (int(iynode(il))-ncty-1)
            if(j.lt.1.or.j.gt.ny) go to 1
            k = iz + (int(iznode(il))-nctz-1)
            if(k.lt.1.or.k.gt.nz) go to 1

            indexloc = (k-1)*nx*ny + (j-1)*nx + i

            if(sim(indexloc).gt.UNEST) then
                  ncnode         = ncnode + 1
                  cnodeid(ncnode) = indexloc
                  icnode(ncnode) = il
            endif
 1          continue

            if(ncnode.eq.nodmax) return
            i = ix + (int(ixnode(il+1))-nctx-1)
            if(i.lt.1.or.i.gt.nx) go to 2
            j = iy + (int(iynode(il+1))-ncty-1)
            if(j.lt.1.or.j.gt.ny) go to 2
            k = iz + (int(iznode(il+1))-nctz-1)
            if(k.lt.1.or.k.gt.nz) go to 2

            indexloc = (k-1)*nx*ny + (j-1)*nx + i

            if(sim(indexloc).gt.UNEST) then
                  ncnode         = ncnode + 1
                  cnodeid(ncnode) = indexloc
                  icnode(ncnode) = il+1
            endif
 2          continue

      end do
c
c Return to calling program:
c
      return
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



      subroutine srchndUnroll(ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,ncnode,
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

      integer i,j,k,il,ilini,indexloc,ncsec
c
c Consider all the nearby nodes until enough have been found:
c
      ncnode = 0
      ncsec  = 0


      do il=1,nlooku,8

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

            if(ncnode.eq.nodmax) return
            i = ix + (int(ixnode(il+1))-nctx-1)
            if(i.lt.1.or.i.gt.nx) go to 2
            j = iy + (int(iynode(il+1))-ncty-1)
            if(j.lt.1.or.j.gt.ny) go to 2
            k = iz + (int(iznode(il+1))-nctz-1)
            if(k.lt.1.or.k.gt.nz) go to 2

c
c Check this potentially informed grid node:
c
            indexloc = (k-1)*nx*ny + (j-1)*nx + i

            if(sim(indexloc).gt.UNEST) then
                  ncnode         = ncnode + 1
                  cnodeid(ncnode) = indexloc
                  icnode(ncnode) = il+1
                  cnodex(ncnode) = xmn + real(i-1)*xsiz
                  cnodey(ncnode) = ymn + real(j-1)*ysiz
                  cnodez(ncnode) = zmn + real(k-1)*zsiz
                  cnodev(ncnode) = sim(indexloc)
                  cnodet(ncnode) = tmp(indexloc)

            endif
 2          continue

            if(ncnode.eq.nodmax) return
            i = ix + (int(ixnode(il+2))-nctx-1)
            if(i.lt.1.or.i.gt.nx) go to 3
            j = iy + (int(iynode(il+2))-ncty-1)
            if(j.lt.1.or.j.gt.ny) go to 3
            k = iz + (int(iznode(il+2))-nctz-1)
            if(k.lt.1.or.k.gt.nz) go to 3

c
c Check this potentially informed grid node:
c
            indexloc = (k-1)*nx*ny + (j-1)*nx + i

            if(sim(indexloc).gt.UNEST) then
                  ncnode         = ncnode + 1
                  cnodeid(ncnode) = indexloc
                  icnode(ncnode) = il+2
                  cnodex(ncnode) = xmn + real(i-1)*xsiz
                  cnodey(ncnode) = ymn + real(j-1)*ysiz
                  cnodez(ncnode) = zmn + real(k-1)*zsiz
                  cnodev(ncnode) = sim(indexloc)
                  cnodet(ncnode) = tmp(indexloc)

            endif
 3          continue

            if(ncnode.eq.nodmax) return
            i = ix + (int(ixnode(il+3))-nctx-1)
            if(i.lt.1.or.i.gt.nx) go to 4
            j = iy + (int(iynode(il+3))-ncty-1)
            if(j.lt.1.or.j.gt.ny) go to 4
            k = iz + (int(iznode(il+3))-nctz-1)
            if(k.lt.1.or.k.gt.nz) go to 4

c
c Check this potentially informed grid node:
c
            indexloc = (k-1)*nx*ny + (j-1)*nx + i

            if(sim(indexloc).gt.UNEST) then
                  ncnode         = ncnode + 1
                  cnodeid(ncnode) = indexloc
                  icnode(ncnode) = il+3
                  cnodex(ncnode) = xmn + real(i-1)*xsiz
                  cnodey(ncnode) = ymn + real(j-1)*ysiz
                  cnodez(ncnode) = zmn + real(k-1)*zsiz
                  cnodev(ncnode) = sim(indexloc)
                  cnodet(ncnode) = tmp(indexloc)

            endif
 4          continue

            if(ncnode.eq.nodmax) return
            i = ix + (int(ixnode(il+4))-nctx-1)
            if(i.lt.1.or.i.gt.nx) go to 5
            j = iy + (int(iynode(il+4))-ncty-1)
            if(j.lt.1.or.j.gt.ny) go to 5
            k = iz + (int(iznode(il+4))-nctz-1)
            if(k.lt.1.or.k.gt.nz) go to 5

c
c Check this potentially informed grid node:
c
            indexloc = (k-1)*nx*ny + (j-1)*nx + i

            if(sim(indexloc).gt.UNEST) then
                  ncnode         = ncnode + 1
                  cnodeid(ncnode) = indexloc
                  icnode(ncnode) = il+4
                  cnodex(ncnode) = xmn + real(i-1)*xsiz
                  cnodey(ncnode) = ymn + real(j-1)*ysiz
                  cnodez(ncnode) = zmn + real(k-1)*zsiz
                  cnodev(ncnode) = sim(indexloc)
                  cnodet(ncnode) = tmp(indexloc)

            endif
 5          continue

            if(ncnode.eq.nodmax) return
            i = ix + (int(ixnode(il+5))-nctx-1)
            if(i.lt.1.or.i.gt.nx) go to 6
            j = iy + (int(iynode(il+5))-ncty-1)
            if(j.lt.1.or.j.gt.ny) go to 6
            k = iz + (int(iznode(il+5))-nctz-1)
            if(k.lt.1.or.k.gt.nz) go to 6

c
c Check this potentially informed grid node:
c
            indexloc = (k-1)*nx*ny + (j-1)*nx + i

            if(sim(indexloc).gt.UNEST) then
                  ncnode         = ncnode + 1
                  cnodeid(ncnode) = indexloc
                  icnode(ncnode) = il+5
                  cnodex(ncnode) = xmn + real(i-1)*xsiz
                  cnodey(ncnode) = ymn + real(j-1)*ysiz
                  cnodez(ncnode) = zmn + real(k-1)*zsiz
                  cnodev(ncnode) = sim(indexloc)
                  cnodet(ncnode) = tmp(indexloc)

            endif
 6          continue

            if(ncnode.eq.nodmax) return
            i = ix + (int(ixnode(il+6))-nctx-1)
            if(i.lt.1.or.i.gt.nx) go to 7
            j = iy + (int(iynode(il+6))-ncty-1)
            if(j.lt.1.or.j.gt.ny) go to 7
            k = iz + (int(iznode(il+6))-nctz-1)
            if(k.lt.1.or.k.gt.nz) go to 7

c
c Check this potentially informed grid node:
c
            indexloc = (k-1)*nx*ny + (j-1)*nx + i

            if(sim(indexloc).gt.UNEST) then
                  ncnode         = ncnode + 1
                  cnodeid(ncnode) = indexloc
                  icnode(ncnode) = il+6
                  cnodex(ncnode) = xmn + real(i-1)*xsiz
                  cnodey(ncnode) = ymn + real(j-1)*ysiz
                  cnodez(ncnode) = zmn + real(k-1)*zsiz
                  cnodev(ncnode) = sim(indexloc)
                  cnodet(ncnode) = tmp(indexloc)

            endif
 7          continue

            if(ncnode.eq.nodmax) return
            i = ix + (int(ixnode(il+7))-nctx-1)
            if(i.lt.1.or.i.gt.nx) go to 8
            j = iy + (int(iynode(il+7))-ncty-1)
            if(j.lt.1.or.j.gt.ny) go to 8
            k = iz + (int(iznode(il+7))-nctz-1)
            if(k.lt.1.or.k.gt.nz) go to 8

c
c Check this potentially informed grid node:
c
            indexloc = (k-1)*nx*ny + (j-1)*nx + i

            if(sim(indexloc).gt.UNEST) then
                  ncnode         = ncnode + 1
                  cnodeid(ncnode) = indexloc
                  icnode(ncnode) = il+7
                  cnodex(ncnode) = xmn + real(i-1)*xsiz
                  cnodey(ncnode) = ymn + real(j-1)*ysiz
                  cnodez(ncnode) = zmn + real(k-1)*zsiz
                  cnodev(ncnode) = sim(indexloc)
                  cnodet(ncnode) = tmp(indexloc)

            endif
 8          continue



      end do

      ilini=il+1

      do il=ilini,nlooku

            if(ncnode.eq.nodmax) return
            i = ix + (int(ixnode(il))-nctx-1)
            if(i.lt.1.or.i.gt.nx) go to 1000
            j = iy + (int(iynode(il))-ncty-1)
            if(j.lt.1.or.j.gt.ny) go to 1000
            k = iz + (int(iznode(il))-nctz-1)
            if(k.lt.1.or.k.gt.nz) go to 1000

c
c Check this potentially informed grid node:
c
            indexloc = (k-1)*nx*ny + (j-1)*nx + i

            if(sim(indexloc).gt.UNEST) then
                  ncnode         = ncnode + 1
                  cnodeid(ncnode) = indexloc
                  icnode(ncnode) = il
                  cnodex(ncnode) = xmn + real(i-1)*xsiz
                  cnodey(ncnode) = ymn + real(j-1)*ysiz
                  cnodez(ncnode) = zmn + real(k-1)*zsiz
                  cnodev(ncnode) = sim(indexloc)
                  cnodet(ncnode) = tmp(indexloc)

            endif
 1000       continue
      end do

c
c Return to calling program:
c
      return
      end

      subroutine srchndomp(ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,ncnode,
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

      integer i,j,k,il,indexloc,ncsec,ncnodetmp
c
c Consider all the nearby nodes until enough have been found:
c
      ncnode = 0
      ncnodetmp = 0
      ncsec  = 0

c$omp parallel default(firstprivate)
c$omp& shared(ixnode,iynode,iznode,sim,tmp,ncsec,ncnode)

c$omp do schedule(static)
      do il=1,nlooku

c            if(ncnode.eq.nodmax) return
            if(ncnode.eq.nodmax) then
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

            if(sim(indexloc).gt.UNEST.or.
     + tmp(indexloc).gt.0.5) then
                  if(sim(indexloc).le.UNEST.and.
     + tmp(indexloc).gt.0.5) then

c$omp atomic
                        ncsec  = ncsec + 1

                        if(ncsec.gt.maxsec) go to 1
                  end if

c$omp critical
                  ncnode         = ncnode + 1
                  ncnodetmp = ncnode
c$omp end critical
                  cnodeid(ncnodetmp) = indexloc
                  icnode(ncnodetmp) = il
                  cnodex(ncnodetmp) = xmn + real(i-1)*xsiz
                  cnodey(ncnodetmp) = ymn + real(j-1)*ysiz
                  cnodez(ncnodetmp) = zmn + real(k-1)*zsiz
                  cnodev(ncnodetmp) = sim(indexloc)
                  cnodet(ncnodetmp) = tmp(indexloc)

            endif
            endif
 1          continue
      end do
c$omp end do nowait
c$omp end parallel



c
c Return to calling program:
c
      return
      end


