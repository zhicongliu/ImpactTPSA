!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! Distributionclass: Initial distribution of charged beam bunch class in 
!                    Beam module of APPLICATION layer.
! Version: 1.0
! Author: Ji Qiang, LBNL
! Description: This class defines initial distributions for the charged 
!              particle beam bunch information in the accelerator.
! Comments: we have added three attributes to each particle:
!           x,px,y,py,t,pt,charge/mass,charge weight,id
!----------------------------------------------------------------
      module Distributionclass
        use Pgrid2dclass
        use CompDomclass
        use BeamBunchclass      
        use Timerclass
        use NumConstclass
        use PhysConstclass
      contains
        ! sample the particles with intial distribution.
        subroutine sample_Dist(this,distparam,nparam,flagdist,geom,grid,Flagbc,&
                               nchrg,nptlist,qmcclist,currlist)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nparam,Flagbc,nchrg
        double precision, dimension(nparam) :: distparam
        double precision, dimension(nchrg) :: qmcclist,currlist
        integer, dimension(nchrg) :: nptlist
        type (BeamBunch), intent(inout) :: this
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer, intent(in) :: flagdist
        integer :: myid, myidx, myidy,seedsize,i,isize
        real rancheck


        !the following is added new for 2nd random group ....
        do i = 1, 300
          call random_number(rancheck)
        enddo
!        write(6,*)'myid,rancheck=',seedarray,myid,rancheck

        if(flagdist.eq.2) then
          call Gauss3_Dist(this,nparam,distparam,grid,0)
        else if(flagdist.eq.3) then
          call Waterbag_Dist(this,nparam,distparam,grid,0)
        else if(flagdist.eq.19) then
          call readin_Dist(this)
        else
          print*,"Initial distribution not available!!"
          stop
        endif

        end subroutine sample_Dist
       
        subroutine Gauss3_Dist(this,nparam,distparam,grid,flagalloc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,flagalloc
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: sq12,sq34,sq56
        double precision, allocatable, dimension(:,:) :: x1,x2,x3 
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i,j,k,intvsamp,pid
        double precision :: t0,x11

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call random_number(x11)
        !print*,myid,x11

        avgpts = this%Npt/(npx*npy)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        if(flagalloc.eq.1) then
          this%Pts1 = 0.0
        else
          allocate(this%Pts1(9,avgpts))
          this%Pts1 = 0.0
        endif

!        allocate(x1(2,avgpts))
!        allocate(x2(2,avgpts))
!        allocate(x3(2,avgpts))
!        call normVec(x1,avgpts)
!        call normVec(x2,avgpts)
!        call normVec(x3,avgpts)
        
        intvsamp = 1
        allocate(x1(2,intvsamp))
        allocate(x2(2,intvsamp))
        allocate(x3(2,intvsamp))

        do j = 1, avgpts/intvsamp
          call normVec(x1,intvsamp)
          call normVec(x2,intvsamp)
          call normVec(x3,intvsamp)
          do k = 1, intvsamp
            !x-px:
!            call normdv(x1)
!           Correct Gaussian distribution.
            i = (j-1)*intvsamp + k
            this%Pts1(1,i) = xmu1 + sig1*x1(1,k)/sq12
            this%Pts1(2,i) = xmu2 + sig2*(-muxpx*x1(1,k)/sq12+x1(2,k))
            !y-py
!            call normdv(x1)
!           Correct Gaussian distribution.
            this%Pts1(3,i) = xmu3 + sig3*x2(1,k)/sq34
            this%Pts1(4,i) = xmu4 + sig4*(-muypy*x2(1,k)/sq34+x2(2,k))
            !z-pz
!            call normdv(x1)
!           Correct Gaussian distribution.
            this%Pts1(5,i) = xmu5 + sig5*x3(1,k)/sq56
            this%Pts1(6,i) = xmu6 + sig6*(-muzpz*x3(1,k)/sq56+x3(2,k))
          enddo
        enddo
          
        deallocate(x1)
        deallocate(x2)
        deallocate(x3)

        this%Nptlocal = avgpts

        do j = 1, avgpts
          pid = j + myid*avgpts
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo
       
        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Gauss3_Dist

        subroutine normdv(y)
        implicit none
        include 'mpif.h'
        double precision, dimension(2), intent(out) :: y
        double precision :: twopi,x1,x2,epsilon

        epsilon = 1.0e-18

        twopi = 4.0*asin(1.0)
        call random_number(x2)
10      call random_number(x1)
!        x1 = 0.5
!10      x2 = 0.6
        if(x1.eq.0.0) goto 10
!        if(x1.eq.0.0) x1 = epsilon
        y(1) = sqrt(-2.0*log(x1))*cos(twopi*x2)
        y(2) = sqrt(-2.0*log(x1))*sin(twopi*x2)

        end subroutine normdv

        subroutine normVec(y,num)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: num
        double precision, dimension(2,num), intent(out) :: y
        double precision :: twopi,epsilon
        double precision, dimension(num) :: x1,x2
        integer :: i

        epsilon = 1.0e-18

        twopi = 4.0*asin(1.0)
        call random_number(x2)
        call random_number(x1)
        do i = 1, num
          if(x1(i).eq.0.0) x1(i) = epsilon
          y(1,i) = sqrt(-2.0*log(x1(i)))*cos(twopi*x2(i))
          y(2,i) = sqrt(-2.0*log(x1(i)))*sin(twopi*x2(i))
        enddo

        end subroutine normVec

        ! sample the particles with intial distribution 
        ! using rejection method. 
        subroutine Waterbag_Dist(this,nparam,distparam,grid,flagalloc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,flagalloc
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision, dimension(6) :: lcrange 
        double precision, dimension(2) :: gs
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: rootx,rooty,rootz,r1,r2,x1,x2
        double precision :: r3,r4,r5,r6,x3,x4,x5,x6
        integer :: totnp,npy,npx
        integer :: avgpts,numpts,isamz,isamy
        integer :: myid,myidx,myidy,iran,intvsamp,pid,j
        double precision :: t0,x11
        double precision, allocatable, dimension(:) :: ranum6

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call random_number(x11)
        !print*,"x11: ",myid,x11

        avgpts = this%Npt/(npx*npy)
        !if(mod(avgpts,10).ne.0) then
        !  print*,"The number of particles has to be an integer multiple of 10Nprocs" 
        !  stop
        !endif
 
        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        rootx=sqrt(1.-muxpx*muxpx)
        rooty=sqrt(1.-muypy*muypy)
        rootz=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        if(flagalloc.eq.1) then
          this%Pts1 = 0.0
        else
          allocate(this%Pts1(9,avgpts))
          this%Pts1 = 0.0
        endif
        numpts = 0
        isamz = 0
        isamy = 0
        intvsamp = avgpts
        !intvsamp = 10
        allocate(ranum6(6*intvsamp))

        do 
          ! rejection sample.
10        continue 
          isamz = isamz + 1
          if(mod(isamz-1,intvsamp).eq.0) then
            call random_number(ranum6)
          endif
          iran = 6*mod(isamz-1,intvsamp)
          r1 = 2.0*ranum6(iran+1)-1.0
          r2 = 2.0*ranum6(iran+2)-1.0
          r3 = 2.0*ranum6(iran+3)-1.0
          r4 = 2.0*ranum6(iran+4)-1.0
          r5 = 2.0*ranum6(iran+5)-1.0
          r6 = 2.0*ranum6(iran+6)-1.0
          if(r1**2+r2**2+r3**2+r4**2+r5**2+r6**2.gt.1.0) goto 10
          isamy = isamy + 1
          numpts = numpts + 1
          if(numpts.gt.avgpts) exit
!x-px:
          x1 = r1*sqrt(8.0)
          x2 = r2*sqrt(8.0)
          !Correct transformation.
          this%Pts1(1,numpts) = xmu1 + sig1*x1/rootx
          this%Pts1(2,numpts) = xmu2 + sig2*(-muxpx*x1/rootx+x2)
          !Rob's transformation
          !this%Pts1(1,numpts) = (xmu1 + sig1*x1)*xscale
          !this%Pts1(2,numpts) = (xmu2 + sig2*(muxpx*x1+rootx*x2))/xscale
!y-py:
          x3 = r3*sqrt(8.0)
          x4 = r4*sqrt(8.0)
          !correct transformation
          this%Pts1(3,numpts) = xmu3 + sig3*x3/rooty
          this%Pts1(4,numpts) = xmu4 + sig4*(-muypy*x3/rooty+x4)
          !Rob's transformation
          !this%Pts1(3,numpts) = (xmu3 + sig3*x3)*yscale
          !this%Pts1(4,numpts) = (xmu4 + sig4*(muypy*x3+rooty*x4))/yscale
!t-pt:
          x5 = r5*sqrt(8.0)
          x6 = r6*sqrt(8.0)
          !correct transformation
          this%Pts1(5,numpts) = xmu5 + sig5*x5/rootz
          this%Pts1(6,numpts) = xmu6 + sig6*(-muzpz*x5/rootz+x6)
          !Rob's transformation
          !this%Pts1(5,numpts) = (xmu5 + sig5*x5)*zscale
          !this%Pts1(6,numpts) = (xmu6 + sig6*(muzpz*x5+rootz*x6))/zscale
        enddo

        deallocate(ranum6)
          
        this%Nptlocal = avgpts
        !print*,"avgpts: ",avgpts
       
!        print*,avgpts,isamz,isamy
        do j = 1, avgpts
          pid = j + myid*avgpts
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Waterbag_Dist

        subroutine readin_Dist(this)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer :: i,j,jlow,jhigh,avgpts,myid,nproc,ierr,nptot,nleft
        double precision, dimension(9) :: tmptcl
        double precision :: sum1,sum2,x1,x2
 
        call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
 
        open(unit=12,file='partcl.data',status='old')
 
        sum1 = 0.0
        sum2 = 0.0
 
          read(12,*)nptot
          avgpts = nptot/nproc
          nleft = nptot - avgpts*nproc
          if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
          else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
          endif
          allocate(this%Pts1(9,avgpts))
          this%Pts1 = 0.0
          !jlow = myid*avgpts + 1
          !jhigh = (myid+1)*avgpts
          !print*,"avgpts, jlow, and jhigh: ",avgpts,jlow,jhigh
          do j = 1, nptot
            read(12,*)tmptcl(1:9)
            sum1 = sum1 + tmptcl(1)
            sum2 = sum2 + tmptcl(3)
            if( (j.ge.jlow).and.(j.le.jhigh) ) then
              i = j - jlow + 1
              this%Pts1(1:9,i) = tmptcl(1:9)
            endif
!            if(myid.eq.0) print*,i,sum1,sum2
          enddo
          !print*,"sumx1,sumy1: ",sum1/nptot,sum2/nptot
 
          close(12)
 
          this%Nptlocal = avgpts
 
        end subroutine readin_Dist

      end module Distributionclass
