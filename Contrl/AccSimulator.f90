!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! AccSimulatorclass: Linear accelerator simulator class in CONTROL layer.
! Version: symplectic multi-particle tracking beta verion
! Author: Ji Qiang, LBNL
! Description: This class defines functions to set up the initial beam 
!              particle distribution, field information, computational
!              domain, beam line element lattice and run the dynamics
!              simulation through the system.
! Comments:
!----------------------------------------------------------------
      module AccSimulatorclass
        use Pgrid2dclass
        use CompDomclass
        use BeamLineElemclass
        use BeamBunchclass
        use Timerclass
        use Inputclass
        use Outputclass
        use Dataclass
        use PhysConstclass
        use NumConstclass
        use Distributionclass
        use TPSAMod

!        implicit none
        !# of phase dim., num. total and local particles, int. dist. 
        !and restart switch, error study switch, substep for space-charge
        !switch
        integer, private :: Dim, Nplocal,Flagdist,Rstartflg,Flagerr,&
                            Flagsubstep 
        integer, private :: Np

        !# of num. total x, total and local y mesh pts., type of BC, 
        !# of beam elems, type of integrator.
        integer, private :: Nx,Ny,Nz,Nxlocal,Nylocal,Nzlocal,Flagbc,&
                            Nblem,Flagmap,Flagdiag

        !# of processors in column and row direction.
        !integer, private :: npcol, nprow
        integer :: npcol, nprow

        !beam current, kin. energy, part. mass, and charge.
        double precision, private :: BcurrImp,Bkenergy,Bmass,Bcharge,BfreqImp,&
                                     Perdlen,xrad,yrad

        !conts. in init. dist.
        double precision, dimension(21) :: distparamZ

        !1d logical processor array.
        type (Pgrid2d) :: grid2d

        !beam particle object and array.
        type (BeamBunch)  :: Bpts

        !beam charge density and field potential arrays.

        !geometry object.
        type (CompDom), private :: Ageom

        !the following variables are used for restart purpose
        integer, private :: iend,jend,ibalend,nstepend
        double precision, private :: zend,zblengend

        !beam line element array.
        type (DriftTube),private,target,dimension(Ndriftmax) :: beamln1
        type (ConstFoc),private,target,dimension(Ncfmax) :: beamln7
        type (BeamLineElem),private,dimension(Nblemtmax)::Blnelem
        !beam line element period.
        interface construct_AccSimulator
          module procedure init_AccSimulator
        end interface

        !//total # of charge state
        integer, private :: nchrg
        !//current list of charge state.
        !//charge/mass list of charge state.
        double precision, private, dimension(100) :: currlist,qmcclist
        !//number of particles of charge state.
        integer, private, dimension(100) :: nptlist
        integer, private, allocatable, dimension(:) :: nptlist0
        double precision, private, allocatable, dimension(:) :: currlist0,qmcclist0
        !integer :: maxptsZ, nslicefel !for FEL use
        !real*8 :: zminfel,lamdas
      contains
        !set up objects and parameters.
        subroutine init_AccSimulator(time)
        implicit none
        include 'mpif.h'
        integer :: i,test1,test2,j
        integer :: myid,myidx,myidy,ierr,inb,jstp
        integer, allocatable, dimension(:) :: bnseg,bmpstp,bitype
        double precision, allocatable, dimension(:) :: blength,val1,&
        val2, val3,val4,val5,val6,val7,val8,val9,val10,val11,val12,val13,val14,&
        val15,val16,val17,val18,val19,val20,val21,val22,val23,val24
        double precision :: time
        double precision :: t0
        double precision :: z,phsini
        double precision, dimension(3) :: tmpdr 
        double precision, dimension(5) :: tmpcf 
        integer :: iqr,idr,ibpm,iccl,iccdtl,idtl,isc,icf,islrf,isl,idipole,&
                   iemfld,myrank,imultpole,itws,nfileout
        integer, allocatable, dimension(:) :: seedarray
        real rancheck
        integer :: meanpts20,seedsize
        !integer :: iseedglb = -10

        !start up MPI.
        call init_Input(time)

        ! initialize Timer.
        call construct_Timer(0.0d0)

        call starttime_Timer(t0)

        !Flagmap = 0

        nptlist = 0
        currlist = 0.0
        qmcclist = 0.0
!-------------------------------------------------------------------
! get all global input parameters.
        call in_Input(Dim,Np,Nx,Ny,Nz,Flagbc,Flagdist,Rstartflg,&
              Flagmap,distparamZ,21,BcurrImp,Bkenergy,Bmass,Bcharge,&
        BfreqImp,xrad,yrad,Perdlen,Nblem,npcol,nprow,Flagerr,Flagdiag,&
        Flagsubstep,phsini,nchrg,nptlist,currlist,qmcclist)


        allocate(nptlist0(nchrg))
        allocate(currlist0(nchrg))
        allocate(qmcclist0(nchrg))
        do i = 1, nchrg
          nptlist0(i) =  nptlist(i)
          currlist0(i) =  currlist(i)
          qmcclist0(i) =  qmcclist(i)
        enddo
!        print*,"npt0: ",nptlist0
!        print*,"qmcc0: ",qmcclist0
!-------------------------------------------------------------------
! construct 2D logical processor Cartesian coordinate
        call construct_Pgrid2d(grid2d,MPI_COMM_WORLD,nprow,npcol)
        call getpost_Pgrid2d(grid2d,myid,myidy,myidx)


        if(myid.eq.0) then
          print*,"Start ImpactZ simulation: symplectic beta version"
        endif

        !construct Constants class.
        call construct_PhysConst(BfreqImp)


!-------------------------------------------------------------------
! construct computational domain CompDom class and get local geometry 
! information on each processor.
        !if(Rstartflg.eq.1) then
        !  call ingeom_Output(1500,z,inb,jstp,nprow,npcol,Ageom,Nx,Ny,Nz,&
        !                    myidx,myidy)
        !  if(myid.eq.0) print*,"rstart at: ",z,inb,jstp
        !  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !else
          !xrad = 0.1363243029*0.2
          call construct_CompDom(Ageom,distparamZ,21,Flagdist,&
               Nx,Ny,Nz,grid2d,nprow,npcol,Flagbc,xrad,yrad,Perdlen)
        !endif

!-------------------------------------------------------------------
! initialize Data class.
        call init_Data()

!-------------------------------------------------------------------
! construct BeamBunch class.
        call construct_BeamBunch(Bpts,BcurrImp,Bkenergy,Bmass,Bcharge,&
                            Np,phsini)

        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

!-------------------------------------------------------------------
! construct FieldQuant class objects.

!-------------------------------------------------------------------
! construct beam line elements.
        allocate(blength(Nblem),bnseg(Nblem),bmpstp(Nblem))
        allocate(bitype(Nblem))
        allocate(val1(Nblem),val2(Nblem),val3(Nblem),val4(Nblem))
        allocate(val5(Nblem),val6(Nblem),val7(Nblem),val8(Nblem))
        allocate(val9(Nblem),val10(Nblem),val11(Nblem),val12(Nblem))
        allocate(val13(Nblem),val14(Nblem),val15(Nblem),val16(Nblem))
        allocate(val17(Nblem),val18(Nblem),val19(Nblem),val20(Nblem))
        allocate(val21(Nblem),val22(Nblem),val23(Nblem),val24(Nblem))

        call in_Input(Nblem,blength,bnseg,bmpstp,bitype,val1,val2,val3,&
        val4,val5,val6,val7,val8,val9,val10,val11,val12,val13,val14,&
        val15,val16,val17,val18,val19,val20,val21,val22,val23,val24)

        iccl = 0
        iccdtl = 0
        idtl = 0
        isc = 0
        idr = 0
        iqr = 0
        ibpm = 0
        icf = 0
        islrf = 0
        isl = 0
        idipole = 0
        iemfld = 0
        imultpole = 0
        itws = 0
        do i = 1, Nblem
          if(bitype(i).eq.0) then
            idr = idr + 1
            call construct_DriftTube(beamln1(idr),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmpdr(1) = 0.0
            tmpdr(2) = val1(i)
            tmpdr(3) = val2(i) !artificial vacuum acceleration gradient
            call setparam_DriftTube(beamln1(idr),tmpdr)
            Blnelem(i) = assign_BeamLineElem(beamln1(idr))
          else if(bitype(i).eq.2) then
            icf = icf + 1
            call construct_ConstFoc(beamln7(icf),bnseg(i),bmpstp(i),&
            bitype(i),blength(i))
            tmpcf(1) = 0.0
            tmpcf(2) = val1(i)
            tmpcf(3) = val2(i)
            tmpcf(4) = val3(i)
            tmpcf(5) = val4(i)
            call setparam_ConstFoc(beamln7(icf),tmpcf)
            Blnelem(i) = assign_BeamLineElem(beamln7(icf))
          else
          endif 
        enddo
!-------------------------------------------------------------------
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(myid.eq.0) print*,"pass setting up lattice..."
!-------------------------------------------------------------------
! sample initial particle distribution.
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        iend = 0
        jend = 0
        ibalend = 0
        nstepend = 0
        zend = 0.0d0
        zblengend = 0.0d0

        !meanpts20 = (Np/totnp)*20
        call random_seed(SIZE=seedsize)
        allocate(seedarray(seedsize))
        do i = 1, seedsize
            seedarray(i) = 10.0d0 + myid*6*20+i*1.0d0*myid
        enddo
        call random_seed(PUT=seedarray)
        call random_number(rancheck)
        !the following is added new for 2nd random group ....
        !do i = 1, 3000
        do i = 1, 3000*(myid+1)
           call random_number(rancheck)
        enddo
        deallocate(seedarray)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        print*,"check random: ",myid,rancheck
        
        Rstartflg = 0
        if(Rstartflg.eq.1) then
        else
            call sample_Dist(Bpts,distparamZ,21,Flagdist,Ageom,grid2d,Flagbc,&
                         nchrg,nptlist0,qmcclist0,currlist0)
        endif
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(myid.eq.0) print*,"pass generating initial distribution..."

        !get local particle number and mesh number on each processor.
        call getnpt_BeamBunch(Bpts,Nplocal)


        deallocate(blength,bnseg,bmpstp,bitype)
        deallocate(val1,val2,val3,val4,val5,val6,val7,val8,val9)
        deallocate(val10,val11,val12,val13,val14,val15,val16)
        deallocate(val17,val18,val19,val20,val21,val22,val23)
        deallocate(val24)

        t_init = t_init + elapsedtime_Timer(t0)
        if(myid.eq.0) print*,"pass init ImpactZ..."

        end subroutine init_AccSimulator

        !Run beam dynamics simulation through accelerator.
        subroutine run_AccSimulator()
        !implicit none
        include 'mpif.h'
        !....
        integer :: i,j,bnseg,bmpstp,bitype
        integer :: myid,myidx,myidy,comm2d,commcol,commrow,npx,npy,&
                   totnp,ierr,nylcr,nzlcr
        integer :: nbal,ibal,nstep,ifile
        integer :: tmpfile,nfile
        integer, dimension(3) :: lcgrid
        double precision :: z0,z,tau1,tau2,blength,t0
        double precision, dimension(6) :: lcrange, range, ptrange,ptref
        double precision, dimension(3) :: msize
        double precision :: hy,hz,ymin,zmin,piperad,zedge
        double precision :: tmp1,tmp2,tmp3,tmp4,rfile
        integer :: nmod,k,ii,jj
        !double precision :: sumtest, sumtest2, sumtest3
        double precision, dimension(12) :: drange
        double precision, dimension(3) :: al0,ga0,epson0
        double precision :: realSamplePeriod,tg,tv,gam,piperad2
        integer :: nsubstep,integerSamplePeriod,Flagbctmp
        double precision :: zz,vref
        !parameters for stripper modeling
        integer :: itmp,jtmp
        double precision :: xx,yy,t3dstart,rr,tmger,tmpwk
        double precision  :: aawk,ggwk,lengwk,hzwake,ab,tbtwstart
        integer :: flagwake,flagwakeread,iizz,iizz1,kz,kadd,ipt,flagcsr,flagbtw
        !for bending magnet Transport transfer matrix implementation
        double precision, dimension(6) :: ptarry
        double precision, dimension(10) :: dparam
        real*8 :: tmpgamlc,tmpgamgl
        real*8, dimension(2) :: tmp56lc,tmp56gl
        real*8 :: tmppx,tmppy,tmppt,tmph
        
        type (dctps) ,dimension(6) :: tpsaPtc

!-------------------------------------------------------------------
! prepare initial parameters, allocate temporary array.
        call getpost_Pgrid2d(grid2d,myid,myidy,myidx)
        idproc = myid
        call getcomm_Pgrid2d(grid2d,comm2d,commcol,commrow)
        call getsize_Pgrid2d(grid2d,totnp,npy,npx)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        !get local particle number and mesh number on each processor.
        call getnpt_BeamBunch(Bpts,Nplocal)

        nstep = nstepend
        z = zend

          if(Flagdiag.eq.1) then
            call diagnostic1_Output(z,Bpts,nchrg,nptlist0)
          endif

          flagdecomp = 0

          do i = 0, nprow - 1
            do j = 0, npcol - 1
              Ageom%LcTabnm(1,i,j) = Nz/nprow
              Ageom%LcTabnm(2,i,j) = Ny/npcol
            enddo
          enddo
          Ageom%Mshlocal(3) = Ageom%LcTabnm(1,myidx,myidy)
          Ageom%Mshlocal(2) = Ageom%LcTabnm(2,myidx,myidy)
          Ageom%Mshlocal(1) = Ageom%Meshnum(1)

!-------------------------------------------------------------------
! prepare for round pipe, open longitudinal
        if(Flagbc.eq.3) then
        endif

!--------------------------------------------------------------------
!output initial line density and energy spread
        qchg = Bpts%Current/Scfreq
        pmass = Bpts%Mass
        nfile =44
        gamma0 = -Bpts%refptcl(6)
!--------------------------------------------------------------------
!for matching purpose to avoid dispersion
!        Bpts%Pts1(5,1:Nplocal) = 0.0d0
!        Bpts%Pts1(6,1:Nplocal) = 0.0d0

        !set initial flagcoll
        flagcoll = 1
        if(BcurrImp.lt.1.0d-15) then
          flagcoll = 0
        endif

!-------------------------------------------------------------------
! start looping through 'Nblem' beam line elements.
        tmpfile = 0
        bitypeold = 0
        blengthold = 0.0d0
        zbleng = zblengend
          Flagmap = 2
          if(Flagmap == 2) then
            call dctps_Initialize(6,5)
            do j=1,6            
              call assign(tpsaPtc(j),0.d0,j)
            enddo
          endif
        !do i=1,6
        !  write(*,*) tpsaPtc(i)%map(2),tpsaPtc(i)%map(3),tpsaPtc(i)%map(4),tpsaPtc(i)%map(5),tpsaPtc(i)%map(6),tpsaPtc(i)%map(7)
        !enddo        
        do i = iend+1, Nblem
          call getparam_BeamLineElem(Blnelem(i),blength,bnseg,bmpstp,&
                                     bitype)
          call getradius_BeamLineElem(Blnelem(i),piperad,piperad2)
          nsubstep = bmpstp
          if(myid.eq.0) then
            print*,"enter elment: ",i,bitype,flagwake,flagwakeread
          endif

          tau1 = 0.0d0
          if(bitype.ge.0) tau1 = 0.5d0*blength/bnseg
          tau2 = 2.0d0*tau1

!-------------------------------------------------------------------
! print out beam information using BPM 

          if(bitype.eq.-99) then
            exit
          endif

!-------------------------------------------------------------------
! loop through 'bnseg' numerical segments in each beam element
! using 2 step symplectic integeration (ie. leap frog).
          zedge = z
          call setparam_BeamLineElem(Blnelem(i),1,zedge)
          if(myid.eq.0) print*,"zedge: ",zedge,Flagbc

          !//no bend or bend using Transport transfer map
!-------------------------------------------------------------------
! use linear map or nonlinear Lorentz integrator to advance particles.
          do j = 1, bnseg
            if(bitype.ne.4) then
              ! spatial drift.
              !linear map integrator
              if(Flagmap.eq.1) then
                call map1_BeamBunch(Bpts,Blnelem(i),z,tau1,bitype,&
                                    bnseg,j,ihlf)
              else if(Flagmap.eq.2) then
                call map1_BeamBunch(Bpts,Blnelem(i),z,tau1,bitype,&
                                    bnseg,j,ihlf,tpsaPtc)
              else
              endif
            else
            endif

!-------------------------------------------------------------------
! escape the space charge calculation for 0 current case
!            if(BcurrImp.lt.1.0d-15)  then !no space-charge
!            else if(flagcoll.eq.1) then !calculate space charge forces
!!-------------------------------------------------------------------
!              !for gridless symp.
!              call conv0thB_BeamBunch(Bpts,tau2,Nplocal,Np,ptrange,&
!                                   Flagbc,Perdlen,piperad,piperad2)
!
!              ! assign new 'Nplocal' local particles on each processor.
!              call setnpt_BeamBunch(Bpts,Nplocal)
!
!              gam = -Bpts%refptcl(6)
!              !print*, myid
!              call kicklmnsin_Field(Bpts%Nptlocal,Bpts%Pts1,Nx,Ny,Nz,&
!                         xrad,yrad,Bpts%Npt,gam,tau2,Bpts%Current,Bpts%Mass,&
!                         Bpts%Charge,ptrange(5))
!            endif
!-------------------------------------------------------------------

            if(bitype.ne.4) then
              if(Flagmap.eq.1) then
                call map1_BeamBunch(Bpts,Blnelem(i),z,tau1,bitype,&
                                    bnseg,j,ihlf)
              else if(Flagmap.eq.2) then
                call map1_BeamBunch(Bpts,Blnelem(i),z,tau1,bitype,&
                                    bnseg,j,ihlf,tpsaPtc)
              else
              endif
            else
            endif

            nstep = nstep + 1
            if(myid.eq.0) then 
              print*,"j, nstep, z",j,nstep,z
            endif
          end do
          call diagnostic1_Output(z,Bpts,nchrg,nptlist0)
          zbleng = zbleng + blength
        enddo

        do i=1,6
          write(*,*) tpsaPtc(i)%map(2),tpsaPtc(i)%map(3),tpsaPtc(i)%map(4),tpsaPtc(i)%map(5),tpsaPtc(i)%map(6),tpsaPtc(i)%map(7)
        enddo        
        call tpsaPtc(1)%output()
        !call tpsaPtc(2)%output()
        !call tpsaPtc(5)%output()
        !call tpsaPtc(6)%output()

! final output.
        call MPI_BARRIER(comm2d,ierr)
        !output all particles in 6d phase space.
        !call phase_Output(100,Bpts,1)
        t_integ = t_integ + elapsedtime_Timer(t0)
        call showtime_Timer()

        !final output line density and uncorrelated energy spread
 
444             format(4(1x,e16.8))

        end subroutine run_AccSimulator

        subroutine destruct_AccSimulator(time)
        implicit none
        include 'mpif.h'
        double precision :: time
 
        call destruct_Data()
        call destruct_BeamBunch(Bpts)
        call destruct_CompDom(Ageom)

        deallocate(nptlist0)
        deallocate(currlist0)
        deallocate(qmcclist0)
 
        call end_Output(time)

        end subroutine destruct_AccSimulator

      end module AccSimulatorclass
