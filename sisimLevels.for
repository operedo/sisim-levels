C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C                                                                      %
C Copyright (C) 2003, Statios Software and Services Incorporated.  All %
C rights reserved.                                                     %
C                                                                      %
C This program has been modified from the one distributed in 1996 (see %
C below).  This version is also distributed in the hope that it will   %
C be useful, but WITHOUT ANY WARRANTY. Compiled programs based on this %
C code may be redistributed without restriction; however, this code is %
C for one developer only. Each developer or user of this source code   %
C must purchase a separate copy from Statios.                          %
C                                                                      %
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C                                                                      %
C Copyright (C) 1996, The Board of Trustees of the Leland Stanford     %
C Junior University.  All rights reserved.                             %
C                                                                      %
C The programs in GSLIB are distributed in the hope that they will be  %
C useful, but WITHOUT ANY WARRANTY.  No author or distributor accepts  %
C responsibility to anyone for the consequences of using them or for   %
C whether they serve any particular purpose or work at all, unless he  %
C says so in writing.  Everyone is granted permission to copy, modify  %
C and redistribute the programs in GSLIB, but only under the condition %
C that this notice and the above copyright notice remain intact.       %
C                                                                      %
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Module to declare dynamic arrays in multiple subroutines:
c
c      module      geostat
c
cc      implicit none
c      
c      real,allocatable   :: x(:),y(:),z(:),vr(:,:),close(:),actloc(:),
c     +                      tmpdat(:),order(:),gcut(:),sim(:),tmp(:),   
c     +                      gcdf(:),covtab(:,:,:,:),cnodex(:),cnodey(:),
c     +                      cnodez(:),cnodev(:),cnodet(:),vra(:),
c     +                      thres(:),cdf(:),ccdf(:),ccdfo(:),beez(:),
c     +                      c0(:),cmax(:),cc(:),aa(:),ang1(:),ang2(:),
c     +                      ang3(:),anis1(:),anis2(:),aviol(:),xviol(:)
c      real*8,allocatable  :: r(:),rr(:),s(:),a(:)
c      integer,allocatable :: ixnode(:),iynode(:),iznode(:),nisb(:),
c     +                       icnode(:),ixsbtosr(:),iysbtosr(:),
c     +                       izsbtosr(:),it(:),nst(:),nviol(:)
c      real*8,allocatable  :: rotmat(:,:,:)
c      logical,allocatable :: atnode(:),softdat(:)
c
c      end module
c
c
c
      program main
c-----------------------------------------------------------------------
c
c           Conditional Simulation of a 3-D Rectangular Grid
c           ************************************************
c
c The output file will be a GEOEAS file containing the simulated values
c The file is ordered by x,y,z, and then simulation (i.e., x cycles
c fastest, then y, then z, then simulation number).
c
c
c
c-----------------------------------------------------------------------

c      use geostat
#ifdef _OPENMP
      use omp_lib
#endif
#ifdef TRACE
      use extrae_module
#endif
c      implicit none


      include  'sisimLevels.inc'



      real,allocatable   :: x(:),y(:),z(:),vr(:,:),close(:),actloc(:),
c     +                      tmpdat(:),order(:),gcut(:),sim(:),tmp(:),   
     +                      tmpdat(:),gcut(:),sim(:),tmp(:),   
     +                      gcdf(:),covtab(:,:,:,:),cnodex(:),cnodey(:),
     +                      cnodez(:),cnodev(:),cnodet(:),vra(:),
     +                      thres(:),cdf(:),ccdf(:),ccdfo(:),beez(:),
     +                      c0(:),cmax(:),cc(:),aa(:),ang1(:),ang2(:),
     +                      ang3(:),anis1(:),anis2(:),aviol(:),xviol(:)

      real,allocatable   :: simThreads(:,:)
      real*8,allocatable  :: r(:),rr(:),s(:),a(:)
c      integer,allocatable :: ixnode(:),iynode(:),iznode(:),nisb(:),
      integer,allocatable :: ixnode(:),iynode(:),iznode(:),nisb(:),
     +                       order(:),icnode(:),ixsbtosr(:),iysbtosr(:),
     +                       izsbtosr(:),it(:),nst(:),nviol(:)

      integer, allocatable:: aclose(:),lock(:)

      real*8,allocatable  :: rotmat(:,:,:)
      logical,allocatable :: atnode(:),softdat(:)

      integer,allocatable :: ncnodeIndex(:),icnodeIndex(:,:)
      integer,allocatable :: mapIndexCount(:)
      real,allocatable    ::cnodexIndex(:,:),cnodeyIndex(:,:),
     +               cnodezIndex(:,:),cnodevIndex(:,:),cnodetIndex(:,:)
c      integer,allocatable :: ncnodeIndex(:,:),icnodeIndex(:,:,:)
c      real,allocatable    ::cnodexIndex(:,:,:),cnodeyIndex(:,:,:),
c     +               cnodezIndex(:,:,:),cnodevIndex(:,:,:),
c     +               cnodetIndex(:,:,:)

      integer nodmax
      integer nctx,ncty,nctz,nlooku,ncnode,maxsec
      integer mapIndexCountCache

      parameter(MV=20)
      real      var(MV)
      real*8    p,acorni
      integer,allocatable :: ivrs(:)
      integer   test
      character datafl*512,tabfl*512,softfl*512,outfl*512,dbgfl*512,
     +          str*512,title*80
      logical   testfl
      integer MAXTHREADS
#ifdef _OPENMP
      parameter (MAXTHREADS=40)
      integer numIterations(MAXTHREADS)
#endif
      real      ntviol,atviol
c      integer ncnodeaux,ncsecaux

      integer BUFFSIZE,MAXOP1AUX
      parameter(BUFFSIZE=2)
      parameter(MAXOP1AUX=13)
      real  simbuffer(BUFFSIZE)

      integer finLoop
      integer extractValue
      integer :: kernelsisimcwrapper 
      integer :: kernelsisimcudawrapper
     
c levels variables
      integer isim,ind,idum,imult
      integer nnx,nny,nnz,jx,jy,jz,ix,iy,iz,index
      integer c,d,e,f,g,h
      integer i,j,id
      real TINY
      logical testind
      real xx,yy,zz
      integer id2
      real test1,test2,test3
      integer nclose, irepo,icut
      integer in
      real zval,cdfval
      integer ic,ierr
      integer nxysim
      integer,allocatable :: cnodeid(:)
      integer,allocatable :: cnodeidIndex(:,:)
c      integer,allocatable :: cnodeidIndex(:,:,:)

      integer lev,maxLevel
      integer levIni, levFin
      integer levIniLocal, levFinLocal
      integer blocknumber
      integer lastEnd
      integer numberOfLevels
      integer, allocatable ::level(:),
     +  indexSort(:)
c,arrayIndex(:) 
c, levelSort(:)
      real,allocatable :: cdfvalIndex(:)
      integer count
      integer lastCount
      integer,allocatable :: levelCount(:)
      integer,allocatable :: levelStart(:)
      integer levelThreshold

      integer clock_max
      integer clock_rate
      integer clock_start
      integer clock_stop

      integer clock_max_loop
      integer clock_rate_loop
      integer clock_start_loop
      integer clock_stop_loop

      external srchndvectorized
      external srchndvectorized2

 
      integer threadId,numThreads,ithread,ilock
      real    invNumThreads
#ifdef _OPENMP
      integer loutThreads(MAXTHREADS)
      character outflThreads(43,MAXTHREADS) , outfltmp*43
#endif

#ifdef GRIDFS
      character var_name*120
      character gridfs_uri*120
      character current_sim*8
      integer*4 status
      character hostname*16
#ifdef _OPENMP
      character threadId_name*8
#endif
      gridfs_uri = GRIDFS // CHAR(0)
#endif

#ifdef TRACE
      call extrae_init()
#endif

      numThreads=1
#ifdef _OPENMP
c$omp parallel
      numThreads = OMP_get_num_threads()
c      print *,numThreads
c$omp end parallel
#endif


c
c Fortran unit numbers needed:
c
      lin  = 1
      lout = 2
#ifdef _OPENMP
      do i=1,MAXTHREADS
            loutThreads(i)=lout*100+i
      end do
#endif
      ldbg = 3
c
c Note VERSION number:
c
      write(*,9999) VERSION
 9999 format(/' SISIM Version: ',f5.3/)
c
c Get the name of the parameter file - try the default name if no input:
c
      do i=1,512
            str(i:i) = ' '
      end do
      call getarg(1,str)
      if(str(1:1).eq.' ')then
            write(*,*) 'Which parameter file do you want to use?'
            read (*,'(a)') str
      end if
      if(str(1:1).eq.' ') str(1:20) = 'sisim.par           '
      inquire(file=str,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the parameter file does not exist,'
            write(*,*) '        check for the file and try again  '
            write(*,*)
            if(str(1:20).eq.'sisim.par           ') then
                  write(*,*) '        creating a blank parameter file'
                  call makepar
                  write(*,*)
            end if
            stop
      endif
      open(lin,file=str,status='OLD')
c
c Find Start of Parameters:
c
 1    read(lin,'(a4)',end=98) str(1:4)
      if(str(1:4).ne.'STAR') go to 1
c
c Read Input Parameters:
c

      read(lin,*,err=98) ivtype
      write(*,*) ' variable type (1=continuous, 0=categorical)= ',ivtype

      read(lin,*,err=98) ncut
      write(*,*) ' number of thresholds / categories = ',ncut
c
c Find the needed parameters:
c
      MAXCUT = ncut
      MAXROT = MAXCUT * MAXNST + 1
      MXCUT = MAXCUT + 1
c
c Allocate the needed memory:
c34
      allocate(thres(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 34: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c35
      allocate(cdf(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 35: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c36
      allocate(ccdf(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 36: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c37
      allocate(ccdfo(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 37: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c38
      allocate(beez(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 38: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c39
      allocate(c0(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 39: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c40
      allocate(cmax(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 40: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c41
      allocate(cc(MAXCUT*MAXNST),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 41: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c42
      allocate(aa(MAXCUT*MAXNST),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 42: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c43
      allocate(ang1(MAXCUT*MAXNST),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 43: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c44
      allocate(ang2(MAXCUT*MAXNST),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 44: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c45
      allocate(ang3(MAXCUT*MAXNST),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 45: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c46
      allocate(anis1(MAXCUT*MAXNST),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 46: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c47
      allocate(anis2(MAXCUT*MAXNST),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 47: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c48
      allocate(aviol(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 48: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c49
      allocate(xviol(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 49: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c50
      allocate(it(MAXCUT*MAXNST),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 50: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c51
      allocate(nst(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 51: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c52
      allocate(nviol(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 52: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c53
      allocate(rotmat(MAXROT,3,3),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 53: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c54
      allocate(ivrs(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 54: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c
      read(lin,*,err=98) (thres(i),i=1,ncut)
      write(*,*) ' thresholds / categories = ',(thres(i),i=1,ncut)

      read(lin,*,err=98) (cdf(i),i=1,ncut)
      write(*,*) ' global cdf / pdf        = ',(cdf(i),i=1,ncut)

      read(lin,'(a512)',err=98) datafl
      call chknam(datafl,512)
      write(*,*) ' data file = ',datafl(1:40)

      read(lin,*,err=98) ixl,iyl,izl,ivrl
      write(*,*) ' input columns = ',ixl,iyl,izl,ivrl

      read(lin,'(a512)',err=98) softfl
      call chknam(softfl,512)
      write(*,*) ' soft data file = ',softfl(1:40)
      inquire(file=softfl,exist=testfl)

      if(testfl) then
            read(lin,*,err=98) ixs,iys,izs,(ivrs(i),i=1,ncut)
            write(*,*) ' columns = ',ixs,iys,izs,(ivrs(i),i=1,ncut)
            read(lin,*,err=98) imbsim
            write(*,*) ' Markov-Bayes simulation = ',imbsim
            if(imbsim.eq.1) then
                  read(lin,*,err=98) (beez(i),i=1,ncut)
            else
                  read(lin,*,err=98)
            end if
      else
            read(lin,*,err=98)
            read(lin,*,err=98)
            read(lin,*,err=98)
      end if

      read(lin,*,err=98) tmin,tmax
      write(*,*) ' trimming limits      ',tmin,tmax

      read(lin,*,err=98) zmin,zmax
      write(*,*) ' data limits (tails)  ',zmin,zmax

      read(lin,*,err=98) ltail,ltpar
      write(*,*) ' lower tail = ',ltail,ltpar

      read(lin,*,err=98) middle,mpar
      write(*,*) ' middle = ',middle,mpar

      read(lin,*,err=98) utail,utpar
      write(*,*) ' upper tail = ',utail,utpar

      read(lin,'(a512)',err=98) tabfl
      call chknam(tabfl,512)
      write(*,*) ' file for tab. quant. ',tabfl(1:40)

      read(lin,*,err=98) itabvr,itabwt
      write(*,*) ' columns for vr wt = ',itabvr,itabwt

      read(lin,*,err=98) idbg
      write(*,*) ' debugging level = ',idbg

      read(lin,'(a512)',err=98) dbgfl
      call chknam(dbgfl,512)
      write(*,*) ' debugging file = ',dbgfl(1:40)

      read(lin,'(a512)',err=98) outfl
      call chknam(outfl,512)
      write(*,*) ' output file = ',outfl(1:40)

#ifdef _OPENMP
c avoid creation of files when GRIDFS is defined
#ifndef GRIDFS
      do i=1,numThreads
            write(outfltmp,"(A40,I3)") trim(outfl(1:40)),
     + loutThreads(i)
            outfltmp = adjustl(trim(outfltmp))
            write(*,*) ' output file (threads) = ',
     + outfltmp
      end do
#endif
#endif

      read(lin,*,err=98) nsim
      write(*,*) ' number of simulations = ',nsim

      read(lin,*,err=98) nx,xmn,xsiz
      write(*,*) ' X grid specification = ',nx,xmn,xsiz

      read(lin,*,err=98) ny,ymn,ysiz
      write(*,*) ' Y grid specification = ',ny,ymn,ysiz

      read(lin,*,err=98) nz,zmn,zsiz
      write(*,*) ' Z grid specification = ',nz,zmn,zsiz
      nxy  = nx*ny
      nxyz = nx*ny*nz

      read(lin,*,err=98) ixv(1)
      write(*,*) ' random number seed = ',ixv(1)
      do i=1,1000
             p = acorni(idum)
      end do

      read(lin,*,err=98) ndmax
      write(*,*) ' ndmax = ',ndmax

      read(lin,*,err=98) nodmax
      write(*,*) ' max prev sim nodes = ',nodmax

      read(lin,*,err=98) maxsec
      write(*,*) ' max soft indicator data = ',maxsec

      read(lin,*,err=98) sstrat
      write(*,*) ' search strategy = ',sstrat

      read(lin,*,err=98) mults,nmult
      write(*,*) ' multiple grid search flag = ',mults,nmult

      read(lin,*,err=98) noct
      write(*,*) ' max per octant = ',noct

      read(lin,*,err=98) radius,radius1,radius2
      write(*,*) ' search radii = ',radius,radius1,radius2
      if(radius.lt.EPSLON) stop 'radius must be greater than zero'
      radsqd = radius  * radius
      sanis1 = radius1 / radius
      sanis2 = radius2 / radius

      read(lin,*,err=98) sang1,sang2,sang3
      write(*,*) ' search anisotropy angles = ',sang1,sang2,sang3

      read(lin,*,err=98) mxctx,mxcty,mxctz
      write(*,*) ' size of covariance lookup = ',mxctx,mxcty,mxctz
      
      read(lin,*,err=98) mik,cutmik
      write(*,*) ' median IK switch = ',mik,cutmik

      read(lin,*,err=98) ktype
      write(*,*) ' kriging type switch = ',ktype

c
c Output now goes to debugging file:
c
      open(ldbg,file=dbgfl,status='UNKNOWN')
      do i=1,ncut
            read(lin,*,err=98) nst(i),c0(i)
            if(ivtype.eq.0)
     +      write(ldbg,100)  i,thres(i),cdf(i),nst(i),c0(i)
            if(ivtype.eq.1)
     +      write(ldbg,101)  i,thres(i),cdf(i),nst(i),c0(i)
            if(nst(i).gt.MAXNST) stop 'nst is too big'
            istart = 1 + (i-1)*MAXNST
            do j=1,nst(i)
                  index = istart + j - 1
                  read(lin,*,err=98) it(index),cc(index),ang1(index),
     +                               ang2(index),ang3(index)
                  if(it(index).eq.3) STOP 'Gaussian Model Not Allowed!'
                  read(lin,*,err=98) aa(index),aa1,aa2
                  write(ldbg,102)  j,it(index),aa(index),cc(index)
                  anis1(index) = aa1 / max(EPSLON,aa(index))
                  anis2(index) = aa2 / max(EPSLON,aa(index))
                  write(ldbg,103) ang1(index),ang2(index),ang3(index),
     +                            anis1(index),anis2(index)
            end do
      end do
      close(lin)
c
c Find the needed parameters:
c
      MAXX = nx
      MAXY = ny
      MAXZ = nz
      MXYZ = MAXX * MAXY * MAXZ
      MAXCTX = mxctx
      MAXCTY = mxcty
      MAXCTZ = mxctz
      MAXCXY = MAXCTX * MAXCTY
      MAXXYZ = MAXCTX * MAXCTY * MAXCTZ
      MAXSAM = ndmax 
      MAXNOD = nodmax
      MAXKR1 = 2 * MAXNOD + 2 * MAXSAM + 1
      MAXKR2 = MAXKR1 * MAXKR1
      MAXSBX = 1
      if(nx.gt.1)then
            MAXSBX = int(nx/2)
            if(MAXSBX.gt.50)MAXSBX=50
      end if
c
      MAXSBY = 1
      if(ny.gt.1)then
            MAXSBY = int(ny/2)
            if(MAXSBY.gt.50)MAXSBY=50
      end if
c
      MAXSBZ = 1
      if(nz.gt.1)then
            MAXSBZ = int(nz/2)
            if(MAXSBZ.gt.50)MAXSBZ=50
      end if
c
      MAXSB = MAXSBX * MAXSBY * MAXSBZ
      MXSXY = 4 * MAXSBX * MAXSBY
      MXSX = 2 * MAXSBX
      av = 0.0
      ss = 0.0
c
c Find the paramater MAXDAT:
c
      MAXDAT = 1
      inquire(file=datafl,exist=testfl)
      if(testfl) then
            open(lin,file=datafl,status='OLD')
            read(lin,*,err=99)
            read(lin,*,err=99)       nvari
            do i=1,nvari
            read(lin,'()',err=99)
            end do
            MAXDAT = 0
 55         read(lin,*,end=66,err=99) (var(j),j=1,nvari)
            if(var(ivrl).lt.tmin.or.var(ivrl).ge.tmax)go to 55
            MAXDAT = MAXDAT + 1
            go to 55
 66         continue
            rewind(lin)
            close(lin)
      end if
      inquire(file=softfl,exist=testfl)
      if(testfl) then
            open(lin,file=softfl,status='OLD')
            read(lin,*,err=97)
            read(lin,*,err=97) nvari
            do i=1,nvari
            read(lin,'()',err=99)
            end do
 56         read(lin,*,end=67,err=99) (var(j),j=1,nvari)
            MAXDAT = MAXDAT + 1
            go to 56
 67         continue
            close(lin)
      endif
c
c Find the paramater MAXTAB:
c
      MAXTAB = ncut
      inquire(file=tabfl,exist=testfl)
      if(testfl) then
            open(lin,file=tabfl,status='OLD')
            read(lin,*,err=99)
            read(lin,*,err=99)       nvari
            do i=1,nvari
                  read(lin,'()',err=99)
            end do
            MAXTAB = 0
 77         read(lin,*,end=88,err=99) (var(j),j=1,nvari)
            if(var(ivrl).lt.tmin.or.var(ivrl).ge.tmax)go to 77
            MAXTAB = MAXTAB + 1
            go to 77
 88         continue
            rewind(lin)
            close(lin)
      end if
c
c Allocate the needed memory:
c1
      allocate(level(nxyz),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 1: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if
c1
      allocate(indexSort(nxyz),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 1: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if
c1
      allocate(cdfvalIndex(nxyz),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 1: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if

c1
      allocate(x(MAXDAT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 1: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if
c2
      allocate(y(MAXDAT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 2: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if
c3
      allocate(z(MAXDAT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 3: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if
c4
      allocate(vr(MAXDAT,MXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 4: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if
c5

c      integer aclose(MAXKR1)
      allocate(aclose(MAXKR1),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 5: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if

      allocate(close(MAXDAT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 5: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if
c6
      allocate(actloc(MAXDAT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 6: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if
c7
      allocate(tmpdat(MAXDAT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 7: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if
c8
      allocate(sim(MXYZ),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 8: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if

      allocate(lock(MXYZ),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 8: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if


#ifdef _OPENMP
c      allocate(simThreads(MXYZ,MAXTHREADS),stat=test)
c            if(test.ne.0)then
c                  write(*,*)'ERROR 8: Allocation failed due to ',
c     +                  'insufficient memory.'
c                  stop
c            end if
#endif

c9
      MAXORD = MXYZ
      if(MXYZ.lt.MAXCXY) MAXORD=MAXCXY
      allocate(order(MAXORD),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 9: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if
c10
      allocate(tmp(MAXORD),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 10: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if
c11
      allocate(gcut(MAXTAB),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 11: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if
c12
      allocate(gcdf(MAXTAB),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 12: Allocation failed due to ',
     +                  'insufficient memory.'
                  stop
            end if
c13
      allocate(covtab(MAXCTX,MAXCTY,MAXCTZ,MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 13: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c14
      allocate(cnodeid(MAXNOD),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 14: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
      allocate(cnodeidIndex(MAXNOD,MXYZ),stat=test)
c      allocate(cnodeidIndex(MAXNOD,MXYZ,16),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 14: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if

c14
      allocate(cnodex(MAXNOD),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 14: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
      allocate(cnodexIndex(MAXNOD,MXYZ),stat=test)
c      allocate(cnodexIndex(MAXNOD,MXYZ,16),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 14: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if

c15
      allocate(cnodey(MAXNOD),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 15: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
      allocate(cnodeyIndex(MAXNOD,MXYZ),stat=test)
c      allocate(cnodeyIndex(MAXNOD,MXYZ,16),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 14: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if


c16
      allocate(cnodez(MAXNOD),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 16: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
      allocate(cnodezIndex(MAXNOD,MXYZ),stat=test)
c      allocate(cnodezIndex(MAXNOD,MXYZ,16),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 14: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if



c17
      allocate(cnodev(MAXNOD),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 17: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
      allocate(cnodevIndex(MAXNOD,MXYZ),stat=test)
c      allocate(cnodevIndex(MAXNOD,MXYZ,16),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 14: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if



c18
      allocate(cnodet(MAXNOD),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 18: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
      allocate(cnodetIndex(MAXNOD,MXYZ),stat=test)
c      allocate(cnodetIndex(MAXNOD,MXYZ,16),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 14: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if



c19
      allocate(vra(MAXKR1),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 19: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c20
      allocate(r(MAXKR1),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 20: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c21
      allocate(rr(MAXKR1),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 21: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c22
      allocate(s(MAXKR1),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 22: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c23
      allocate(a(MAXKR2),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 23: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c24
      allocate(ixnode(MAXXYZ),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 24: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c25
      allocate(iynode(MAXXYZ),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 25: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c26
      allocate(iznode(MAXXYZ),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 26: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if

c27
      allocate(nisb(MAXSB),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 27: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c28
      allocate(icnode(MAXNOD),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 28: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
      allocate(icnodeIndex(MAXNOD,MXYZ),stat=test)
c      allocate(icnodeIndex(MAXNOD,MXYZ,16),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 28: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
      allocate(ncnodeIndex(MXYZ),stat=test)
c      allocate(ncnodeIndex(MXYZ,16),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 28: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if

      allocate(mapIndexCount(MXYZ),stat=test)

c29
      allocate(ixsbtosr(8*MAXSB),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 29: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c30
      allocate(iysbtosr(8*MAXSB),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 30: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c31
      allocate(izsbtosr(8*MAXSB),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 31: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c32
      allocate(atnode(MAXDAT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 32: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c33
      allocate(softdat(MAXKR1),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 33: Allocation failed due to ',
     +                  'insuffiecent memory.'
                  stop
            end if
c
 100  format(/,' Category  number ',i2,' = ',f12.3,/,
     +         '           global prob value = ',f8.4,/,
     +         '           number of structures = ',i3,/,
     +         '           nugget effect        = ',f8.4)
 101  format(/,' Threshold number ',i2,' = ',f12.3,/,
     +         '           global prob value = ',f8.4,/,
     +         '           number of structures = ',i3,/,
     +         '           nugget effect        = ',f8.4)
 102  format(  '           type of structure ',i3,' = ',i3,/,
     +         '           aa parameter         = ',f12.4,/,
     +         '           cc parameter         = ',f12.4)
 103  format(  '           ang1, ang2, ang3     = ',3f6.2,/,
     +         '           anis1, anis2         = ',2f12.4)
c
c Perform some quick error checking:
c
      if(nx.gt.MAXX) stop 'nx is too big - modify .inc file'
      if(ny.gt.MAXY) stop 'ny is too big - modify .inc file'
      if(nz.gt.MAXZ) stop 'nz is too big - modify .inc file'
c
c Check to make sure the data file exists, then either read in the
c data or write a warning:
c
      title = 'SISIM SIMULATIONS:                      '//
     +        '                                        '
      nd = 0
      inquire(file=datafl,exist=testfl)
      if(.not.testfl) then
            write(*,113) datafl
 113        format('WARNING data file ',a40,' does not exist!',/,
     +             ' Hope your intention was to create an',
     +             ' unconditional simulation.')
      else
c
c The data file exists so open the file and read in the header
c information.
c
            write(*,*) 'Reading input data'
            av = 0.0
            ss = 0.0
            open(lin,file=datafl,status='OLD')
            read(lin,'(a60)',err=99) title(21:80)
            read(lin,*,err=99)       nvari
            do i=1,nvari
                  read(lin,'()',err=99)
            end do
c
c Read all the data until the end of the file:
c
 5          read(lin,*,end=6,err=99) (var(j),j=1,nvari)
            vrt = var(ivrl)
            if(vrt.lt.tmin.or.vrt.ge.tmax) go to 5
            nd  = nd + 1
            x(nd) = xmn
            y(nd) = ymn
            z(nd) = zmn
            if(ixl.gt.0) x(nd) = var(ixl)
            if(iyl.gt.0) y(nd) = var(iyl)
            if(izl.gt.0) z(nd) = var(izl)
            av = av + vrt
            ss = ss + vrt*vrt
c
c The indicator data are constructed knowing the thresholds and the
c data value.
c
            if(ivtype.eq.0) then
                  do ic=1,ncut
                        vr(nd,ic) = 0.0
                        if(int(vrt+0.5).eq.int(thres(ic)+0.5)) 
     +                  vr(nd,ic) = 1.0
                  end do
            else
                  do ic=1,ncut
                        vr(nd,ic) = 1.0
                        if(vrt.gt.thres(ic)) vr(nd,ic) = 0.0
                  end do
            end if
            vr(nd,MXCUT) = vrt
            go to 5
 6          close(lin)
c     
c Compute the averages and variances as an error check for the user:
c
            xd = max(real(nd),1.0)
            av = av / xd
            ss =(ss / xd ) - av * av
            write(*,120)    ivrl,nd,av,ss
            write(ldbg,120) ivrl,nd,av,ss
 120        format(/,'  Data for SISIM: Variable number ',i2,
     +             /,'  Number of acceptable data  = ',i8,
     +             /,'  Equal Weighted Average     = ',f12.4,
     +             /,'  Equal Weighted Variance    = ',f12.4,/)
c
c Check to make sure that the grid is compatible with the data:
c
            if(ixl.le.0.and.nx.gt.1) then
               write(*,*) 'ERROR there is no X coordinate in data file'
               write(*,*) '      nx must be set to 1'
               stop
            end if
            if(iyl.le.0.and.ny.gt.1) then
               write(*,*) 'ERROR there is no Y coordinate in data file'
               write(*,*) '      ny must be set to 1'
               stop
            end if
            if(izl.le.0.and.nz.gt.1) then
               write(*,*) 'ERROR there is no Z coordinate in data file'
               write(*,*) '      nz must be set to 1'
               stop
            end if
      endif
c
c Now, if required, read in the tabulated values for details of the dist
c
      if(ltail.eq.3.or.middle.eq.3.or.utail.eq.3) then
            ng = 0
            inquire(file=tabfl,exist=testfl)
            if(.not.testfl) stop 'ERROR tabfl does not exist'
            open(lin,file=tabfl,status='OLD')
            read(lin,*,err=97)
            read(lin,*,err=97) nvari
            do i=1,nvari
                  read(lin,*,err=97)
            end do
            tcdf = 0.0
            ng   = 0
 21         read(lin,*,end=22,err=97) (var(j),j=1,nvari)
            if(var(itabvr).lt.tmin.or.var(itabvr).ge.tmax) go to 21
            ng = ng + 1
            gcut(ng) = var(itabvr)
            gcdf(ng) = 1.0
            if(itabwt.gt.0) gcdf(ng) = var(itabwt)
            tcdf = tcdf + gcdf(ng)
            go to 21
 22         close(lin)
c
c Sort in ascending order and keep track of where the tabulated values
c switch classes:
c
            if(tcdf.le.0.0) then
                  write(*,*) 'ERROR: either the weights are zero or'
                  write(*,*) '       there are no tabulated data.'
                  stop
            endif
            call sortem(1,ng,gcut,1,gcdf,c,d,e,f,g,h)
c
c Set up gcdf for tabulated quantiles:
c
            oldcp = 0.0
            cp    = 0.0
            tcdf  = 1.0 / tcdf
            do i=1,ng
                  cp      = cp + gcdf(i) * tcdf
                  gcdf(i) =(cp + oldcp) * 0.5
                  oldcp   = cp
            end do
      end if
c
c Direct input of indicator data:
c
      nhd = nd
      inquire(file=softfl,exist=testfl)
      if(testfl) then
            write(*,*)
            write(*,*) 'Reading soft indicator data'
            open(lin,file=softfl,status='OLD')
            read(lin,*,err=97)
            read(lin,*,err=97) nvari
            if(ivrs(ncut).gt.nvari) then
                  write(*,*) ' ERROR: too few variables in ',softfl
                  write(*,*) '        inconsistent with parameters'
                  stop
            end if
            do i=1,nvari
                  read(lin,*,err=97)
            end do
 12         read(lin,*,end=13,err=96) (var(j),j=1,nvari)
c
c Don't keep soft data co-located with hard data:
c
            xx = xmn
            yy = ymn
            zz = zmn
            if(ixs.gt.0) xx = var(ixs)
            if(iys.gt.0) yy = var(iys)
            if(izs.gt.0) zz = var(izs)
            do i=1,nhd
                  test = abs(xx-x(i)) + abs(yy-y(i)) + abs(zz-z(i))
                  if(test.le.EPSLON) go to 12
            end do
c
c Accept this data:
c
            nd = nd + 1
            x(nd) = xx
            y(nd) = yy
            z(nd) = zz
            do j=1,ncut
                  i = ivrs(j)
                  vr(nd,j) = var(i)
                  ccdf(j)  = var(i)
            end do
c
c Draw a value for this soft distribution (in case the distribution is
c co-located with a grid node and Markov-Bayes is not used):
c

            cdfval = real(acorni(idum))
c            cdfval = 0.5
            call ordrel(ivtype,ncut,ccdf,ccdfo,nviol,aviol,xviol)
            zval = UNEST
            call beyond(ivtype,ncut,thres,ccdfo,ng,gcut,gcdf,zmin,zmax,
     +                  ltail,ltpar,middle,mpar,utail,utpar,zval,
     +                  cdfval,ierr)
            vr(nd,MXCUT) = zval
c
c If performing median IK then check for missing values:
c
            if(mik.eq.1) then
                  do ic=1,ncut
                        if(vr(nd,ic).lt.0.0) then
                              write(*,150) softfl
                              stop
                        endif
                  end do
 150              format(' Since the median IK approach is being',
     +                   ' considered no missing values are',
     +                   ' allowed',/,' Check file ',a40)
            endif
            go to 12
 13         close(lin)
            write(*,*) ' finished reading soft data'
      endif
c
c Load the right variogram as the first one if performing median IK:
c
      if(mik.eq.1) then
            icut = 1
            clos = abs(cutmik-thres(1))
            do ic=2,ncut
                  test = abs(cutmik-thres(ic))
                  if(test.lt.clos) then
                        icut = ic
                        clos = test
                  end if
            end do
            c0(1)   = c0(icut)
            nst(1)  = nst(icut)
            istart1 = 1
            istarti = 1 + (icut-1)*MAXNST
            do ist=1,nst(1)
                  index1        = istart1 + ist - 1
                  indexi        = istarti + ist - 1
                  it(index1)    = it(indexi)
                  aa(index1)    = aa(indexi)
                  cc(index1)    = cc(indexi)
                  ang1(index1)  = ang1(indexi)
                  ang2(index1)  = ang2(indexi)
                  ang3(index1)  = ang3(indexi)
                  anis1(index1) = anis1(indexi)
                  anis2(index1) = anis2(indexi)
            end do
      end if
c
c Open the output file and return:
c
      open(lout,file=outfl,status='UNKNOWN')

#ifdef _OPENMP
c each thread will write in the files file.out401, file.out402, so on.
c up to MAXTHREADS except when GRIDFS is defined
#ifndef GRIDFS
      do i=1,numThreads
            write(outfltmp,"(A40,I3)") trim(outfl(1:40)),
     + loutThreads(i)
            outfltmp = adjustl(trim(outfltmp))
            open(100*lout+i,file=outfltmp,status='UNKNOWN')
       end do
#endif
#endif

      write(lout,104) title
 104  format(a80)
      write(lout,105) 1,nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz,nsim
 105  format(i2,3(1x,i4),3(1x,g14.8),3(1x,g12.6),i4) 
      write(lout,106)
 106  format('Simulated Value')
      go to 1111
c
c Error in an Input File Somewhere:
c
 96   stop 'ERROR in soft data file!'
 97   stop 'ERROR in table look up file!'
 98   stop 'ERROR in parameter file!'
 99   stop 'ERROR in data file!'





 1111 print *,'END READING PARAMETERS'
      print *,'START SIMULATION'
c
c Set up the rotation/anisotropy matrices that are needed for the
c variogram and search:
c
      write(*,*) 'Setting up rotation matrices for variogram and search'
      do ic=1,ncut
      do is=1,nst(ic)
            ind = is + (ic-1)*MAXNST
            call setrot(ang1(ind),ang2(ind),ang3(ind),anis1(ind),
     +                  anis2(ind),ind,MAXROT,rotmat)
      end do
      end do
      isrot = MAXNST*MAXCUT + 1
      call setrot(sang1,sang2,sang3,sanis1,sanis2,isrot,MAXROT,rotmat)
c
c Set up for super block searching:
c
      if(sstrat.eq.0) then
            write(*,*) 'Setting up super block search strategy'
            do i=1,nd
                  actloc(i) = real(i)
            end do
            nsec = 0
            call setsupr(nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,nd,x,y,z,
     +             actloc,tmp,nsec,sec1,sec2,sec3,MAXSBX,MAXSBY,MAXSBZ,
     +             nisb,nxsup,xmnsup,xsizsup,nysup,ymnsup,ysizsup,nzsup,
     +             zmnsup,zsizsup)
            call picksup(nxsup,xsizsup,nysup,ysizsup,nzsup,zsizsup,
     +             isrot,MAXROT,rotmat,radsqd,nsbtosr,ixsbtosr,
     +             iysbtosr,izsbtosr)
      end if
c
c Set up the covariance table and the spiral search:
c
      call ctable(MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ,MAXROT,
     + MAXCUT,MAXORD,MAXXYZ,covtab,tmp,order,iznode,iynode,ixnode,
     + MAXNST,aa,c0,cc,cmax,idbg,isrot,it,ldbg,ncut,nlooku,nodmax,
     + nst,nx,ny,nz,xsiz,ysiz,zsiz,radsqd,rotmat,nctx,ncty,nctz
     + )
      
      print *,'calling to kernelSisim...'

#ifdef FORTRAN

      lastEnd=0
      numberOfLevels=0

#ifdef _OPENMP
c$omp parallel default(firstprivate)
      threadId=int(OMP_get_thread_num()) + 1
      print *,'threadId=',threadId
c$omp end parallel
#else
      threadId=1
#endif



      finLoop=nsim

      do isim=1,finLoop

            call system_clock(clock_start,clock_rate,clock_max)


c
c Work out a random path for this realization:
c
            print *,'realization=',isim,nxyz

            do ind=1,nxyz
                  sim(ind)   = real(acorni(idum))
                  order(ind) = ind
            end do

c
c The multiple grid search works with multiples of 4 (yes, that is
c somewhat arbitrary):
c
            if(mults.eq.1) then
                  do imult=1,nmult
                        nnz = max(1,nz/(imult*4))
                        nny = max(1,ny/(imult*4))
                        nnx = max(1,nx/(imult*4))
                        jz  = 1
                        jy  = 1
                        jx  = 1
                        do iz=1,nnz
                           if(nnz.gt.1) jz = iz*imult*4
                           do iy=1,nny
                              if(nny.gt.1) jy = iy*imult*4
                              do ix=1,nnx
                                 if(nnx.gt.1) jx = ix*imult*4
                                 index = jx + (jy-1)*nx + (jz-1)*nxy
                                 sim(index) = sim(index) - imult
                              end do
                           end do
                        end do
                  end do
            end if

            call sortem(1,nxyz,sim,1,order,c,d,e,f,g,h)

            call system_clock(clock_stop,clock_rate,clock_max)
            print *,'TIME (random path)=',
     +    real(clock_stop-clock_start,kind=8)/real(clock_rate,kind=8) 

c
c Initialize the simulation:
c
            call system_clock(clock_start,clock_rate,clock_max)
            do i=1,nxyz
                  sim(i) = UNEST
                  tmp(i) = 0.0
            end do
            write(*,*)
            write(*,*) ' Working on realization number: ',isim
c
c Assign the data to the closest grid node:
c
            TINY = 0.0001
            do id=1,nd
                  call getindx(nx,xmn,xsiz,x(id),ix,testind)
                  call getindx(ny,ymn,ysiz,y(id),iy,testind)
                  call getindx(nz,zmn,zsiz,z(id),iz,testind)
                  ind  = ix + (iy-1)*nx + (iz-1)*nxy
                  xx   = xmn + real(ix-1)*xsiz
                  yy   = ymn + real(iy-1)*ysiz
                  zz   = zmn + real(iz-1)*zsiz
                  test1 = abs(xx-x(id)) + abs(yy-y(id)) + abs(zz-z(id))
c
c Assign this data to the node (unless there is a closer data):
c
                  atnode(id) = .false.
                  if(sstrat.eq.1)                  atnode(id) = .true.
                  if(sstrat.eq.0.and.test1.le.TINY) atnode(id) = .true.
                  if(atnode(id)) then
                        if(sim(ind).ge.0.0) then
                              id2 = int(sim(ind)+0.5)
                              test2 = abs(xx-x(id2)) + abs(yy-y(id2))
     +                                               + abs(zz-z(id2))
                              if(test1.le.test2)then
                                   sim(ind) = real(id)
c                              if(idbg.ge.2) write(ldbg,1002) id,id2
                              end if
                        else
                              sim(ind) = real(id)
                        end if
                  end if
            end do
 1002       format(' WARNING data values ',2i5,' are both assigned to ',
     +           /,'         the same node - taking the closest')

c
c Now, enter the hard data values into the "sim" array and keep the
c data number in the "tmp" array (to be reset when a hard value
c is assigned to that node):
            do i=1,nxyz
                  id = int(sim(i)+0.5)
                  if(id.gt.0) then
                        if(id.le.nhd) then
                              sim(i) = vr(id,MXCUT)
                        else
                              tmp(i) = sim(i)
                              sim(i) = UNEST
                        end if
                  end if
            end do
c
c Accumulate the number and magnitude of order relations violations:
c
            nclose = 0
            irepo  = max(1,min((nxyz/10),10000))
            ntviol = 0.0
            atviol = 0.0
            do icut=1,ncut
                  nviol(icut) =  0
                  aviol(icut) =  0.0
                  xviol(icut) = -1.0
            end do

            call system_clock(clock_stop,clock_rate,clock_max)
            print *,'TIME (initialization)=',
     +    real(clock_stop-clock_start,kind=8)/real(clock_rate,kind=8) 

c
c MAIN LOOP OVER ALL THE NODES:
c

            call system_clock(clock_start,clock_rate,clock_max)
            do in=1,nxyz
               level(in) = 0
            end do

            do in=1,nxyz
                  if((int(in/irepo)*irepo).eq.in) write(*,1004) in
 1004             format('   calculating level of node ',i9)
c                  index = int(order(in)+0.5)
                  index = order(in)
c
c Do we really need to simulate this grid node location?
c


                  if(sim(index).ne.UNEST)then
                        level(index) = 0
                        go to 20
                  end if

                  if(imbsim.eq.0.and.tmp(index).ne.0.0) then
                        id = int(tmp(index)+0.5)
                        sim(index) = vr(id,MXCUT)
                        level(index) = 0
                        go to 20
                  end if
c
c Location of the node we are currently working on:
c
                  iz = int((index-1)/nxy) + 1
                  iy = int((index-(iz-1)*nxy-1)/nx) + 1
                  ix = index - (iz-1)*nxy - (iy-1)*nx
                  xx = xmn + real(ix-1)*xsiz
                  yy = ymn + real(iy-1)*ysiz
                  zz = zmn + real(iz-1)*zsiz
c                  if(idbg.ge.3)
c     +            write(ldbg,*) 'Working on grid index:',ix,iy,iz

c
c Now, we'll simulate the point ix,iy,iz.  First, get the close data
c and make sure that there are enough to actually simulate a value,
c we'll only keep the closest "ndmax" data, and look for previously
c simulated grid nodes:
c
                  if(sstrat.eq.0) then
                        call srchsupr(xx,yy,zz,radsqd,isrot,MAXROT,
     +                       rotmat,nsbtosr,ixsbtosr,iysbtosr,
     +                       izsbtosr,noct,nd,x,y,z,tmpdat,nisb,nxsup,
     +                       xmnsup,xsizsup,nysup,ymnsup,ysizsup,
     +                       nzsup,zmnsup,zsizsup,nclose,close,
     +                       infoct)
                        if(nclose.gt.ndmax) nclose = ndmax
c                       do i=1,nclose
c                             iii = int(close(i)+0.5)
c                             close(i) = real(actloc(iii))
c                       end do
                  endif

                  icnode(:)=0
                  cnodeid(:)=-1
                  ncnode=0

c            call srchndUnroll(ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,ncnode,
c              call srchndList(ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,ncnode,
c                  call srchnd(ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,ncnode,
c       call srchndvectorized(ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,ncnode,
c       call srchndvectorized2(ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,ncnode,
              call srchndPush(ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,ncnode,
     +        maxsec,nctx,ncty,nctz,nlooku,nodmax,nx,ny,nz,UNEST,
     +        xmn,ymn,zmn,xsiz,ysiz,zsiz,icnode,ixnode,iynode,iznode,
     +        cnodex,cnodey,cnodez,cnodev,cnodet,tmp,sim,cnodeid)

                  maxLevel = -1
                  do i=1,ncnode
                     if (level(cnodeid(i)).gt.maxLevel) then
                        maxLevel = level(cnodeid(i))
                     end if
                  end do

                  if(ncnode.eq.0)then
                     level(index) = 0
                     ncnodeIndex(index)=0
                     icnodeIndex(:,index)=-1
                     cnodeidIndex(:,index)=-1
                  else
                     level(index) = maxLevel + 1 
                     ncnodeIndex(index)=ncnode
                     icnodeIndex(:,index)=icnode(:)
                     cnodeidIndex(:,index)=cnodeid(:)
                  end if

                  sim(index) = 999999.0

                  cdfval = real(acorni(idum))
                  cdfvalIndex(index)=cdfval

                  if (maxLevel .gt. numberOfLevels) then
                     numberOfLevels = maxLevel
                  end if

c
c END MAIN LOOP OVER NODES:
c
 20         continue

            tmp(index) = 0.0 

            end do
            
            call system_clock(clock_stop,clock_rate,clock_max)
            print *,'TIME (make levels1)=',
     +    real(clock_stop-clock_start,kind=8)/real(clock_rate,kind=8) 

            call system_clock(clock_start,clock_rate,clock_max)

            levelThreshold = 0
            numberOfLevels = numberOfLevels + 1
            count = 0
            lastCount = 0
            allocate(levelCount(numberOfLevels+1))
            allocate(levelStart(numberOfLevels+1))

            do lev = 0,numberOfLevels
               levelStart(lev+1) = count + 1
               levelCount(lev+1) = 0 
               do in = 1,nxyz
                  if(level(in).eq.lev)then
                     count = count +1
                     indexSort(count) = in
                     levelCount(lev+1) = levelCount(lev+1) + 1
                  end if
               end do
               levelThreshold = levelThreshold + levelCount(lev+1)

            end do

            levelThreshold = int(0.1*ceiling(real(levelThreshold)/
     +                           real(numberOfLevels-1)))


            call system_clock(clock_stop,clock_rate,clock_max)
            print *,'TIME (make levels2)=',
     +    real(clock_stop-clock_start,kind=8)/real(clock_rate,kind=8) 


            call system_clock(clock_start,clock_rate,clock_max)
c
c Initialize (again) the simulation:
c
            do i=1,nxyz
                  sim(i) = UNEST
                  tmp(i) = 0.0
                  lock(i) = 0
            end do
c            write(*,*)
c            write(*,*) ' Working on realization number: ',isim
c
c Assign the data to the closest grid node:
c
            TINY = 0.0001
            do id=1,nd
                  call getindx(nx,xmn,xsiz,x(id),ix,testind)
                  call getindx(ny,ymn,ysiz,y(id),iy,testind)
                  call getindx(nz,zmn,zsiz,z(id),iz,testind)
                  ind  = ix + (iy-1)*nx + (iz-1)*nxy
                  xx   = xmn + real(ix-1)*xsiz
                  yy   = ymn + real(iy-1)*ysiz
                  zz   = zmn + real(iz-1)*zsiz
                  test1 = abs(xx-x(id)) + abs(yy-y(id)) + abs(zz-z(id))
c
c Assign this data to the node (unless there is a closer data):
c
                  atnode(id) = .false.
                  if(sstrat.eq.1)                  atnode(id) = .true.
                  if(sstrat.eq.0.and.test.le.TINY) atnode(id) = .true.
                  if(atnode(id)) then
                        if(sim(ind).ge.0.0) then
                              id2 = int(sim(ind)+0.5)
                              test2 = abs(xx-x(id2)) + abs(yy-y(id2))
     +                                               + abs(zz-z(id2))
                              if(test1.le.test2) sim(ind) = real(id)
c                              if(idbg.ge.2) write(ldbg,1002) id,id2
                        else
                              sim(ind) = real(id)
                        end if
                  end if
            end do


c
c Now, enter the hard data values into the "sim" array and keep the
c data number in the "tmp" array (to be reset when a hard value
c is assigned to that node):
            do i=1,nxyz
                  id = int(sim(i)+0.5)
                  if(id.gt.0) then
                        if(id.le.nhd) then
                              sim(i) = vr(id,MXCUT)
                        else
                              tmp(i) = sim(i)
                              sim(i) = UNEST
                        end if
                  end if
            end do
c
c Accumulate the number and magnitude of order relations violations:
c
            nclose = 0
            irepo  = max(1,min((nxyz/10),10000))
            ntviol = 0.0
            atviol = 0.0
            do icut=1,ncut
                  nviol(icut) =  0
                  aviol(icut) =  0.0
                  xviol(icut) = -1.0
            end do
c
c MAIN LOOP OVER ALL THE NODES:
c


            lev=0

            levIni=levelStart(lev+1) 
            levFin=(levelStart(lev+1)+levelCount(lev+1)-1) 
            print *,'level=',lev,'/',numberOfLevels,
     + 'ini=',levIni,'fin=',levFin

c$omp parallel default(firstprivate)
c$omp& shared(x,y,z,vr,level,indexSort,
c$omp&        cdfvalIndex,ncnodeIndex,icnodeIndex,
c$omp&        cnodeidIndex,cnodexIndex,cnodeyIndex,
c$omp&        cnodezIndex,cnodevIndex,cnodetIndex,
c$omp&        tmp,sim,order,covtab,beez,lock)

#ifdef _OPENMP
            threadId=OMP_get_thread_num()+1
#else
            threadId=1
#endif

c$omp do schedule(runtime)
            do count=levIni,levFin

                  index = indexSort(count)
c
c Do we really need to simulate this grid node location?
c
                  if(sim(index).ne.UNEST)then 
                        lock(index) = 1
                        go to 200
                  end if
                  if(imbsim.eq.0.and.tmp(index).ne.0.0) then
                        id = int(tmp(index)+0.5)
                        sim(index) = vr(id,MXCUT)
                        lock(index) = 1
                        go to 200
                  end if
c
c Location of the node we are currently working on:
c
                  iz = int((index-1)/nxy) + 1
                  iy = int((index-(iz-1)*nxy-1)/nx) + 1
                  ix = index - (iz-1)*nxy - (iy-1)*nx
                  xx = xmn + real(ix-1)*xsiz
                  yy = ymn + real(iy-1)*ysiz
                  zz = zmn + real(iz-1)*zsiz
c
c Now, we'll simulate the point ix,iy,iz.  First, get the close data
c and make sure that there are enough to actually simulate a value,
c we'll only keep the closest "ndmax" data, and look for previously
c simulated grid nodes:
c
                  if(sstrat.eq.0) then
                        call srchsupr(xx,yy,zz,radsqd,isrot,MAXROT,
     +                       rotmat,nsbtosr,ixsbtosr,iysbtosr,
     +                       izsbtosr,noct,nd,x,y,z,tmpdat,nisb,nxsup,
     +                       xmnsup,xsizsup,nysup,ymnsup,ysizsup,
     +                       nzsup,zmnsup,zsizsup,nclose,close,
     +                       infoct)
                        if(nclose.gt.ndmax) nclose = ndmax
                  endif
                  icnode(:)=0
                  cnodex(:)=0.0
                  cnodey(:)=0.0
                  cnodez(:)=0.0
                  cnodev(:)=0.0
                  cnodet(:)=0.0
                  cnodeid(:)=0
                  ncnode=0

                 call srchnd(ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,ncnode,
     +        maxsec,nctx,ncty,nctz,nlooku,nodmax,nx,ny,nz,UNEST,
     +        xmn,ymn,zmn,xsiz,ysiz,zsiz,icnode,ixnode,iynode,iznode,
     +        cnodex,cnodey,cnodez,cnodev,cnodet,tmp,sim,cnodeid)
c
c What cdf value are we looking for?
c
                  zval   = UNEST
                  cdfval = cdfvalIndex(index)
c
c Use the global distribution?
c
                  if((nclose+ncnode).le.0) then
                        call beyond(ivtype,ncut,thres,cdf,ng,gcut,gcdf,
     +                              zmin,zmax,ltail,ltpar,middle,mpar,
     +                              utail,utpar,zval,cdfval,ierr)

                  end if

                  sim(index) = zval
                  lock(index) = 1

 200         continue

                  tmp(index) = 0.0 

            end do
c$omp end do 
c$omp end parallel


            call system_clock(clock_stop,clock_rate,clock_max)
            print *,'TIME (lev',lev,')=',
     +    real(clock_stop-clock_start,kind=8)/real(clock_rate,kind=8) 

            call system_clock(clock_start_loop,clock_rate_loop,
     +                        clock_max_loop)

            invNumThreads = 1.0/real(numThreads)

c$omp parallel default(firstprivate)
c$omp& shared(x,y,z,vr,level,indexSort,
c$omp&        cdfvalIndex,ncnodeIndex,icnodeIndex,
c$omp&        cnodeidIndex,cnodexIndex,cnodeyIndex,
c$omp&        cnodezIndex,cnodevIndex,cnodetIndex,mapIndexCount,
c$omp&        tmp,sim,order,covtab,beez,lock,numberOfLevels, 
c$omp&        levelStart,levelCount,numThreads,invNumThreads)

#ifdef _OPENMP
            threadId=OMP_get_thread_num()+1
#else
            threadId=1
#endif

            do lev=1,(numberOfLevels)


            levIni=levelStart(lev+1) 
            levFin=(levelStart(lev+1)+levelCount(lev+1)-1) 

            if(levelCount(lev+1).ge.-1) then

            if(numThreads.gt.1)then
            blocknumber =ceiling(real(levelCount(lev+1)-numThreads+1)*
     +                         invNumThreads)
               levIniLocal = levIni+
     +                    blocknumber*(threadId-1)
               levFinLocal = levIni+
     +                    blocknumber*threadId-1
               if(threadId.eq.numThreads) levFinLocal=levFin

            else
               levIniLocal=levIni
               levFinLocal=levFin
            end if

            do count=levIniLocal,levFinLocal
            
                  index = indexSort(count)
c
c Location of the node we are currently working on:
c
                  iz = int((index-1)/nxy) + 1
                  iy = int((index-(iz-1)*nxy-1)/nx) + 1
                  ix = index - (iz-1)*nxy - (iy-1)*nx
                  xx = xmn + real(ix-1)*xsiz
                  yy = ymn + real(iy-1)*ysiz
                  zz = zmn + real(iz-1)*zsiz
c                  if(idbg.ge.3)
c     +            write(ldbg,*) 'Working on grid index:',ix,iy,iz

c
c Now, we'll simulate the point ix,iy,iz.  First, get the close data
c and make sure that there are enough to actually simulate a value,
c we'll only keep the closest "ndmax" data, and look for previously
c simulated grid nodes:
c
                  if(sstrat.eq.0) then
                        call srchsupr(xx,yy,zz,radsqd,isrot,MAXROT,
     +                       rotmat,nsbtosr,ixsbtosr,iysbtosr,
     +                       izsbtosr,noct,nd,x,y,z,tmpdat,nisb,nxsup,
     +                       xmnsup,xsizsup,nysup,ymnsup,ysizsup,
     +                       nzsup,zmnsup,zsizsup,nclose,close,
     +                       infoct)
                        if(nclose.gt.ndmax) nclose = ndmax
c                       do i=1,nclose
c                             iii = int(close(i)+0.5)
c                             close(i) = real(actloc(iii))
c                       end do
                  endif

c Instead of calling srchnd, we extract neighbour values from 
c arrays ***Index(:), obtained before...

                 ncnode = ncnodeIndex(index)
                 icnode(:) = icnodeIndex(:,index)
                 cnodeid(:) = cnodeidIndex(:,index)

              call srchndPop(ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,ncnode,
     +        maxsec,nctx,ncty,nctz,nlooku,nodmax,nx,ny,nz,UNEST,
     +        xmn,ymn,zmn,xsiz,ysiz,zsiz,icnode,ixnode,iynode,iznode,
     +        cnodex,cnodey,cnodez,cnodev,cnodet,tmp,sim,cnodeid)


c The values of cnodev are extracted from sim array using the
c indexes obtained before


                  ilock=0
                  do while(ilock.eq.0)
                     ilock=1
                     do i=1,ncnode
                        ilock=ilock*lock(cnodeid(i))
                     end do
                  end do


                  do i=1,ncnode
                     cnodev(i)=sim(cnodeid(i))
                  end do

c
c What cdf value are we looking for?
c
                  zval   = UNEST
                  cdfval = cdfvalIndex(index)
c
c Use the global distribution?
c

                  if((nclose+ncnode).le.0) then
                        call beyond(ivtype,ncut,thres,cdf,ng,gcut,gcdf,
     +                              zmin,zmax,ltail,ltpar,middle,mpar,
     +                              utail,utpar,zval,cdfval,ierr)

                  else
c
c Estimate the local distribution by indicator kriging:
c

                        do ic=1,ncut
                              aclose(:)=0

      call krige(ix,iy,iz,ic,MAXCTX,MAXCTY,MAXCTZ,MAXKR1,
     + MAXROT,MAXDAT,MXCUT,MAXKR2,MAXNOD,MAXCUT,MAXXYZ,MAXNST,idbg,
     + imbsim,ivtype,ktype,ldbg,mik,nclose,ncnode,nctx,ncty,nctz, 
     + xx,yy,zz,cdf(ic),ccdf(ic),atnode,softdat,vr,vra,x,y,z,icnode,
     + iznode,iynode,ixnode,it,nst,aclose,cnodex,cnodey,cnodez,
     + cnodev,cnodet,covtab,thres,aa,c0,cc,cmax,close,beez,r,rr,s,a,
     + rotmat)


                        end do
c
c Correct order relations:
c

                        call ordrel(ivtype,ncut,ccdf,ccdfo,nviol,aviol,
     +                              xviol)

c
c Draw from the local distribution:
c

                        call beyond(ivtype,ncut,thres,ccdfo,ng,gcut,
     +                              gcdf,zmin,zmax,ltail,ltpar,middle,
     +                              mpar,utail,utpar,zval,cdfval,ierr)

c
c Write some debugging information:
c

c                        if(idbg.ge.3) then
c                              do ic=1,ncut
c                              write(ldbg,202) ccdf(ic),ccdfo(ic)
c 202                          format('  CDF (original and fixed)',2f7.4)
c                              end do
c                        endif
                  endif

                  sim(index) = zval
                  lock(index) = 1

c
c END MAIN LOOP OVER NODES:
c
            tmp(index) = 0.0 

           end do

           end if

           end do

c$omp end parallel


            call system_clock(clock_stop_loop,clock_rate_loop,
     +                        clock_max_loop)
            print *,'TIME (loop',lev,')=',
     +    real(clock_stop_loop-clock_start_loop,kind=8)/
     +    real(clock_rate_loop,kind=8) 


c
c Write this simulation to the output file:
c
            nxysim = 0
            do ic=1,ncut
                  ccdf(ic) = 0.0
            end do

            call system_clock(clock_start,clock_rate,
     +                        clock_max)

            do ind=1,nxyz
c
c Calculate the cdf of the simulated values (for error checking):
c

c                  if(sim(ind).gt.UNEST) then
c                        nxysim = nxysim + 1
c                        do ic=1,ncut
c                              if(ivtype.eq.0) then
c                                    if(sim(ind).eq.thres(ic))
c     +                                ccdf(ic)=ccdf(ic)+1.0
c                              else
c                                    if(sim(ind).le.thres(ic))
c     +                                ccdf(ic)=ccdf(ic)+1.0
c                              end if
c                        end do
c                  endif
c

                  write(lout,'(g14.8)') sim(ind)

            end do

            call system_clock(clock_stop,clock_rate,clock_max)
            print *,'TIME (output)=',
     +    real(clock_stop-clock_start,kind=8)/real(clock_rate,kind=8) 

c
c END MAIN LOOP OVER SIMULATIONS:
c
      end do


#endif

      print *,'returning from kernelSisim...'

c
c Finished:
c
      close(lout)
#ifdef _OPENMP
      do i=1,numThreads
            close(lout*100+i)
      end do
#endif
      close(ldbg)
      write(*,9998) VERSION
 9998 format(/' SISIM Version: ',f5.3, ' Finished'/)

#ifdef TRACE
      call extrae_fini()
#endif

      stop
      end
 

      subroutine makepar
c-----------------------------------------------------------------------
c
c                      Write a Parameter File
c                      **********************
c
c
c
c-----------------------------------------------------------------------
      lun = 99
      open(lun,file='sisim.par',status='UNKNOWN')
      write(lun,10)
 10   format('                  Parameters for SISIM',/,
     +       '                  ********************',/,/,
     +       'START OF PARAMETERS:')

      write(lun,11)
 11   format('1                             ',
     +       '-1=continuous(cdf), 0=categorical(pdf)')
      write(lun,12)
 12   format('5                             ',
     +       '-number thresholds/categories')
      write(lun,13)
 13   format('0.5   1.0   2.5   5.0   10.0  ',
     +       '-   thresholds / categories')
      write(lun,14)
 14   format('0.12  0.29  0.50  0.74  0.88  ',
     +       '-   global cdf / pdf')
      write(lun,15)
 15   format('../data/cluster.dat           ',
     +       '-file with data')
      write(lun,16)
 16   format('1   2   0   3                 ',
     +       '-   columns for X,Y,Z, and variable')
      write(lun,17)
 17   format('direct.ik                     ',
     +       '-file with soft indicator input')
      write(lun,18)
 18   format('1   2   0   3 4 5 6 7         ',
     +       '-   columns for X,Y,Z, and indicators')
      write(lun,19)
 19   format('0                             ',
     +       '-   Markov-Bayes simulation (0=no,1=yes)')
      write(lun,20)
 20   format('0.61  0.54  0.56  0.53  0.29  ',
     +       '-      calibration B(z) values')
      write(lun,21)
 21   format('-1.0e21    1.0e21             ',
     +       '-trimming limits')
      write(lun,22)
 22   format('0.0   30.0                    ',
     +       '-minimum and maximum data value')
      write(lun,23)
 23   format('1      0.0                    ',
     +       '-   lower tail option and parameter')
      write(lun,24)
 24   format('1      1.0                    ',
     +       '-   middle     option and parameter')
      write(lun,25)
 25   format('1     30.0                    ',
     +       '-   upper tail option and parameter')
      write(lun,26)
 26   format('cluster.dat                   ',
     +       '-   file with tabulated values')
      write(lun,27)
 27   format('3   0                         ',
     +       '-      columns for variable, weight')
      write(lun,28)
 28   format('0                             ',
     +       '-debugging level: 0,1,2,3')
      write(lun,29)
 29   format('sisim.dbg                     ',
     +       '-file for debugging output')
      write(lun,30)
 30   format('sisim.out                     ',
     +       '-file for simulation output')
      write(lun,31)
 31   format('1                             ',
     +       '-number of realizations')
      write(lun,32)
 32   format('50   0.5    1.0               ',
     +       '-nx,xmn,xsiz')
      write(lun,33)
 33   format('50   0.5    1.0               ',
     +       '-ny,ymn,ysiz')
      write(lun,34)
 34   format('1    1.0   10.0               ',
     +       '-nz,zmn,zsiz')
      write(lun,35)
 35   format('69069                         ',
     +       '-random number seed')
      write(lun,36)
 36   format('12                            ',
     +       '-maximum original data  for each kriging')
      write(lun,37)
 37   format('12                            ',
     +       '-maximum previous nodes for each kriging')
      write(lun,38)
 38   format('1                             ',
     +       '-maximum soft indicator nodes for kriging')
      write(lun,39)
 39   format('1                             ',
     +       '-assign data to nodes? (0=no,1=yes)')
      write(lun,40)
 40   format('0     3                       ',
     +       '-multiple grid search? (0=no,1=yes),num')
      write(lun,41)
 41   format('0                             ',
     +       '-maximum per octant    (0=not used)')
      write(lun,42)
 42   format('20.0  20.0  20.0              ',
     +       '-maximum search radii')
      write(lun,43)
 43   format(' 0.0   0.0   0.0              ',
     +       '-angles for search ellipsoid')
      write(lun,44)
 44   format('51    51    11                ',
     +       '-size of covariance lookup table')
      write(lun,47)
 47   format('0    2.5                      ',
     +       '-0=full IK, 1=median approx. (cutoff)')
      write(lun,48)
 48   format('0                             ',
     +       '-0=SK, 1=OK')
      write(lun,49)
 49   format('1    0.15                     ',
     +       '-One   nst, nugget effect')
      write(lun,50)
 50   format('1    0.85 0.0   0.0   0.0     ',
     +       '-      it,cc,ang1,ang2,ang3')
      write(lun,51)
 51   format('         10.0  10.0  10.0     ',
     +       '-      a_hmax, a_hmin, a_vert')
      write(lun,52)
 52   format('1    0.10                     ',
     +       '-Two   nst, nugget effect')
      write(lun,53)
 53   format('1    0.90 0.0   0.0   0.0     ',
     +       '-      it,cc,ang1,ang2,ang3')
      write(lun,54)
 54   format('         10.0  10.0  10.0     ',
     +       '-      a_hmax, a_hmin, a_vert')
      write(lun,55)
 55   format('1    0.10                     ',
     +       '-Three nst, nugget effect')
      write(lun,56)
 56   format('1    0.90 0.0   0.0   0.0     ',
     +       '-      it,cc,ang1,ang2,ang3')
      write(lun,57)
 57   format('         10.0  10.0  10.0     ',
     +       '-      a_hmax, a_hmin, a_vert')
      write(lun,58)
 58   format('1    0.10                     ',
     +       '-Four  nst, nugget effect')
      write(lun,59)
 59   format('1    0.90 0.0   0.0   0.0     ',
     +       '-      it,cc,ang1,ang2,ang3')
      write(lun,60)
 60   format('         10.0  10.0  10.0     ',
     +       '-      a_hmax, a_hmin, a_vert')
      write(lun,61)
 61   format('1    0.15                     ',
     +       '-Five  nst, nugget effect')
      write(lun,62)
 62   format('1    0.85 0.0   0.0   0.0     ',
     +       '-      it,cc,ang1,ang2,ang3')
      write(lun,63)
 63   format('         10.0  10.0  10.0     ',
     +       '-      a_hmax, a_hmin, a_vert')

      close(lun)
      return
      end
