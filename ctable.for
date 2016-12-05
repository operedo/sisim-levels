      subroutine ctable(MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ,MAXROT,
     + MAXCUT,MAXORD,MAXXYZ,covtab,tmp,order,iznode,iynode,ixnode,
     + MAXNST,aa,c0,cc,cmax,idbg,isrot,it,ldbg,ncut,nlooku,nodmax,
     + nst,nx,ny,nz,xsiz,ysiz,zsiz,radsqd,rotmat,nctx,ncty,nctz
     + )
c-----------------------------------------------------------------------
c
c               Establish the Covariance Look up Table
c               **************************************
c
c The idea is to establish a 3-D network that contains the covariance
c value for a range of grid node offsets that should be at as large
c as twice the search radius in each direction.  The reason it has to
c be twice as large as the search radius is because we want to use it
c to compute the data covariance matrix as well as the data-block
c covariance matrix.
c
c Secondly, we want to establish a search for nearby nodes that 
c in order of closeness as defined by the variogram.
c
c
c
c INPUT VARIABLES:
c
c   xsiz,ysiz,zsiz  Definition of the grid being considered
c   MAXCTX,Y,Z      Number of blocks in covariance table
c
c   covariance table parameters
c
c
c
c OUTPUT VARIABLES:  covtab()         Covariance table
c
c EXTERNAL REFERENCES:
c
c   sqdist          Computes 3-D anisotropic squared distance
c   sortem          Sorts multiple arrays in ascending order
c   cova3           Computes the covariance according to a 3-D model
c
c
c
c-----------------------------------------------------------------------
c      use geostat

      implicit none

      integer MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ,MAXROT,
     + MAXCUT,MAXORD,MAXXYZ
      real    covtab(MAXCTX,MAXCTY,MAXCTZ,MAXCUT),tmp(MAXORD),
     + order(MAXORD)
      integer iznode(MAXXYZ),iynode(MAXXYZ),ixnode(MAXXYZ)
      integer MAXNST
      real    aa(MAXCUT*MAXNST),c0(MAXCUT),cc(MAXCUT*MAXNST)  
      real    cmax(MAXCUT) 
      integer idbg,isrot
      integer it(MAXCUT*MAXNST) 
      integer ldbg,ncut,nlooku,nodmax
      integer nst(MAXCUT)
      integer nx,ny,nz
      real    xsiz,ysiz,zsiz
      real    radsqd
      real*8  rotmat(MAXROT,3,3) 
      integer nctx,ncty,nctz

      real TINY
      parameter(TINY=1.0e-10)
c      include  'sisim.inc'
      real*8    sqdist,hsqd
      integer i,j,k,ilooku,icut,irot,ic,jc,kc
      integer il,loc,ix,iy,iz
      real    xx,yy,zz 
      real    c(1),d(1),e(1),f(1),g(1),h(1)
      real    cbb 
      real    cov
c
c Size of the look-up table:
c
      nctx = min(((MAXCTX-1)/2),(nx-1))
      ncty = min(((MAXCTY-1)/2),(ny-1))
      nctz = min(((MAXCTZ-1)/2),(nz-1))
c
c Initialize the covariance subroutine and cbb at the same time:
c
      call cova3(0.,0.,0.,0.,0.,0.,1,nst,MAXNST,c0,it,cc,aa,1,
     +           MAXROT,rotmat,cmax,cbb)
c
c Now, set up the table and keep track of the node offsets that are
c within the search radius:
c
      ilooku = max((ncut/2),1)
      nlooku = 0
      do icut=1,ncut
      irot = 1 + (icut-1)*MAXNST
      do i=-nctx,nctx
      xx = i * xsiz
      ic = nctx + 1 + i
      do j=-ncty,ncty
      yy = j * ysiz
      jc = ncty + 1 + j
      do k=-nctz,nctz
      zz = k * zsiz
      kc = nctz + 1 + k
            call cova3(0.,0.,0.,xx,yy,zz,icut,nst,MAXNST,c0,it,cc,aa,
     +                 irot,MAXROT,rotmat,cmax,cov)
            covtab(ic,jc,kc,icut) = cov
c            print *,'TAG cov = ',covtab(ic,jc,kc,icut) 
            if(icut.eq.ilooku) then
c               print *,'xx=',xx,'yy=',yy,'zz=',zz,'isrot=',isrot,
c     + 'MAXROT=',MAXROT,'rotmat=',rotmat(:,:,:)
               hsqd = sqdist(0.0,0.0,0.0,xx,yy,zz,isrot,MAXROT,rotmat)
c               print *,'hsqd=',hsqd,'radsqd=',radsqd,'cov=',cov
               if(real(hsqd).le.radsqd) then
                  nlooku = nlooku + 1
c
c We subtract the covariance from a large value so that the ascending
c sort subroutine will accomplish the sort we want.  Furthermore, a
c fraction of the distance is also taken off so that we search by
c anisotropic distance once we are beyond the range:
c
                  tmp(nlooku)   =-(covtab(ic,jc,kc,icut)-TINY*hsqd)
                  order(nlooku) =real((kc-1)*MAXCXY+(jc-1)*MAXCTX+ic)
               endif 
            endif
      end do
      end do
      end do
      end do

c      print *,'order-bef',order(:)


c
c Finished setting up the look-up table, now order the nodes such
c that the closest ones, according to variogram distance, are searched
c first. Note: the "loc" array is used because I didn't want to make 
c special allowance for 2 byte integers in the sorting subroutine:
c
      call sortem(1,nlooku,tmp,1,order,c,d,e,f,g,h)
c      print *,'order-aft',order(:)
      do il=1,nlooku
            loc = int(order(il))
            iz  = int((loc-1)/MAXCXY) + 1
            iy  = int((loc-(iz-1)*MAXCXY-1)/MAXCTX) + 1
            ix  = loc-(iz-1)*MAXCXY - (iy-1)*MAXCTX
            iznode(il) = iz
            iynode(il) = iy
            ixnode(il) = ix
      end do
      if(nodmax.gt.MAXNOD) then
            write(ldbg,*)
            write(ldbg,*) 'The maximum number of close nodes = ',nodmax
            write(ldbg,*) 'this was reset from your specification due '
            write(ldbg,*) 'to storage limitations.'
            nodmax = MAXNOD
      endif
c
c Debugging output if requested:
c
      if(idbg.le.2) return
      write(ldbg,*)
      write(ldbg,*) 'There are ',nlooku,' nearby nodes that will be '
      write(ldbg,*) 'checked until enough close data are found.'
      write(ldbg,*)
      if(idbg.lt.4) return
      do i=1,nlooku
            xx = (ixnode(i) - nctx - 1) * xsiz
            yy = (iynode(i) - ncty - 1) * ysiz
            zz = (iznode(i) - nctz - 1) * zsiz
            write(ldbg,100) i,xx,yy,zz
      end do
 100  format('Point ',i6,' at ',3f18.6)
c
c All finished:
c
      return
      end


