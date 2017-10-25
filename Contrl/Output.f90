!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! Outputclass: Output class in I/O module of CONTROL layer. 
! Version: 1.0
! Author: Ji Qiang, LBNL
! Description: This class defines functions to print out the charged
!              particle beam information in the accelerator.
! Comments: We have added 3 more attributes to the particle:
!           x,px,y,py,t,pt,charge/mass,charge weight,id
!----------------------------------------------------------------
      module Outputclass
!        use Timer_class
        use BeamBunchclass
        use PhysConstclass
      contains
        ! calculate <x^2>,<xp>,<px^2>,x emittance, <y^2>,<ypy>,
        ! <py^2> and y emittance, <z^2>,<zp>,<pz^2>,z emittance.
        subroutine diagnostic1_Output(z,this,nchrg,nptlist)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: z
        type (BeamBunch), intent(in) :: this
        integer, intent(in) :: nchrg
        integer, dimension(nchrg), intent(in) :: nptlist
        integer :: innp,nptot
        double precision:: den1,den2,sqsum1,sqsum2,sqsum3,sqsum4,&
                          epsx2,epsy2
        double precision:: xpx,ypy,epx,epy,xrms,pxrms,yrms,pyrms,&
                         xpxfac,ypyfac
        double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
        double precision:: xpxlocal,ypylocal,zpzlocal
        double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms,pzrms,zpzfac
        double precision:: sqsum5local,sqsum6local
        double precision:: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
        pz0lc,pz0,x0lc3,x0lc4,px0lc3,px0lc4,y0lc3,y0lc4,py0lc3,py0lc4,&
        z0lc3,z0lc4,pz0lc3,pz0lc4,x03,x04,px03,px04,y03,y04,py03,py04,&
        z03,z04,pz03,pz04
        double precision ::sqx,cubx,fthx,sqpx,cubpx,fthpx,sqy,cuby,fthy,&
        sqpy,cubpy,fthpy,sqz,cubz,fthz,sqpz,cubpz,fthpz
        double precision:: gam,energy,bet,gambet
        integer :: i,my_rank,ierr,j
        double precision:: qmc,xl,xt
        double precision, dimension(6) :: localmax, glmax
        double precision, dimension(27) :: tmplc,tmpgl
        double precision :: t0,lcrmax,glrmax
!        double precision :: alphax,alphay,alphaz
!       for the calculation of bunched beam current.
!        double precision, dimension(3) :: localaa, aa, localbb, bb 
!        double precision :: pi,tmp5,tmp55
!        integer :: ntmp5,nf
        integer :: npctmin,npctmax

        call starttime_Timer(t0)

        qmc = this%Mass/1.0e6
        xl = Scxl
        xt = Rad2deg

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

        innp = this%Nptlocal
        nptot = this%Npt
        !print*,"pts1: ",this%Pts1(1,1),this%Pts1(2,1),my_rank,innp

        den1 = 1.0/float(nptot)
        den2 = den1*den1
        x0lc = 0.0
        px0lc = 0.0
        y0lc = 0.0
        py0lc = 0.0
        z0lc = 0.0
        pz0lc = 0.0
        sqsum1local = 0.0
        sqsum2local = 0.0
        sqsum3local = 0.0
        sqsum4local = 0.0
        sqsum5local = 0.0
        sqsum6local = 0.0
        xpxlocal = 0.0
        ypylocal = 0.0
        zpzlocal = 0.0
        x0lc3 = 0.0
        x0lc4 = 0.0
        px0lc3 = 0.0
        px0lc4 = 0.0
        y0lc3 = 0.0
        y0lc4 = 0.0
        py0lc3 = 0.0
        py0lc4 = 0.0
        z0lc3 = 0.0
        z0lc4 = 0.0
        pz0lc3 = 0.0
        pz0lc4 = 0.0

!        localaa = 0.0
!        localbb = 0.0
!        aa = 0.0
!        bb = 0.0
!        pi = 2*asin(1.0)

        ! for cache optimization.
        if(innp.ne.0) then
          do i = 1, 6
            localmax(i) = abs(this%Pts1(i,1))
          enddo
          lcrmax = this%Pts1(1,1)**2+this%Pts1(3,1)**2
        else
          do i = 1, 6
            localmax(i) = 0.0
          enddo
          lcrmax = 0.0
        endif
        do i = 1, innp
          x0lc = x0lc + this%Pts1(1,i)
          sqsum1local = sqsum1local + this%Pts1(1,i)*this%Pts1(1,i)
          x0lc3 = x0lc3 + this%Pts1(1,i)*this%Pts1(1,i)*this%Pts1(1,i)
          x0lc4 = x0lc4 + this%Pts1(1,i)*this%Pts1(1,i)*this%Pts1(1,i)*&
                  this%Pts1(1,i)
          xpxlocal = xpxlocal + this%Pts1(1,i)*this%Pts1(2,i)
          px0lc = px0lc + this%Pts1(2,i)
          sqsum2local = sqsum2local + this%Pts1(2,i)*this%Pts1(2,i)
          px0lc3 = px0lc3 + this%Pts1(2,i)*this%Pts1(2,i)*this%Pts1(2,i)
          px0lc4 = px0lc4 + this%Pts1(2,i)*this%Pts1(2,i)*this%Pts1(2,i)*&
                   this%Pts1(2,i)
          y0lc = y0lc + this%Pts1(3,i)
          sqsum3local = sqsum3local + this%Pts1(3,i)*this%Pts1(3,i)
          y0lc3 = y0lc3 + this%Pts1(3,i)*this%Pts1(3,i)*this%Pts1(3,i)
          y0lc4 = y0lc4 + this%Pts1(3,i)*this%Pts1(3,i)*this%Pts1(3,i)*&
                  this%Pts1(3,i)
          ypylocal = ypylocal + this%Pts1(3,i)*this%Pts1(4,i)
          py0lc = py0lc + this%Pts1(4,i)
          sqsum4local = sqsum4local + this%Pts1(4,i)*this%Pts1(4,i)
          py0lc3 = py0lc3 + this%Pts1(4,i)*this%Pts1(4,i)*this%Pts1(4,i)
          py0lc4 = py0lc4 + this%Pts1(4,i)*this%Pts1(4,i)*this%Pts1(4,i)*&
                   this%Pts1(4,i)
          z0lc = z0lc + this%Pts1(5,i)
          sqsum5local = sqsum5local + this%Pts1(5,i)*this%Pts1(5,i)
          z0lc3 = z0lc3 + abs(this%Pts1(5,i)*this%Pts1(5,i)*this%Pts1(5,i))
          z0lc4 = z0lc4 + this%Pts1(5,i)*this%Pts1(5,i)*this%Pts1(5,i)*&
                          this%Pts1(5,i)

          !for the calculation of bunched beam current
!          ntmp5 = this%Pts1(5,i)/pi
!          tmp5 = this%Pts1(5,i) - ntmp5*pi
!          tmp55 = tmp5 - mod(ntmp5,2)*pi
!          do nf = 1, 3
!            localaa(nf) = localaa(nf) + cos(nf*tmp55)
!            localbb(nf) = localbb(nf) + sin(nf*tmp55)
!          enddo

          zpzlocal = zpzlocal + this%Pts1(5,i)*this%Pts1(6,i)
          pz0lc = pz0lc + this%Pts1(6,i)
          sqsum6local = sqsum6local + this%Pts1(6,i)*this%Pts1(6,i)
          pz0lc3 = pz0lc3 + this%Pts1(6,i)*this%Pts1(6,i)*this%Pts1(6,i)
          pz0lc4 = pz0lc4 + this%Pts1(6,i)*this%Pts1(6,i)*this%Pts1(6,i)*&
                            this%Pts1(6,i)
          do j = 1, 6
            if(localmax(j).lt.abs(this%Pts1(j,i))) then
               localmax(j) = abs(this%Pts1(j,i))
            endif
          enddo
          if(lcrmax.lt.(this%Pts1(1,i)**2+this%Pts1(3,i)**2)) then
            lcrmax = this%Pts1(1,i)**2 + this%Pts1(3,i)**2
          endif
        enddo

!        do i = 1, 6
!          localmax(i) = maxval(abs(this%Pts1(i,1:innp)))
!        enddo

        tmplc(1) = x0lc
        tmplc(2) = px0lc
        tmplc(3) = y0lc
        tmplc(4) = py0lc
        tmplc(5) = z0lc
        tmplc(6) = pz0lc
        tmplc(7) = sqsum1local
        tmplc(8) = sqsum2local
        tmplc(9) = sqsum3local
        tmplc(10) = sqsum4local
        tmplc(11) = sqsum5local
        tmplc(12) = sqsum6local
        tmplc(13) = xpxlocal
        tmplc(14) = ypylocal
        tmplc(15) = zpzlocal
        tmplc(16) = x0lc3
        tmplc(17) = x0lc4
        tmplc(18) = px0lc3
        tmplc(19) = px0lc4
        tmplc(20) = y0lc3
        tmplc(21) = y0lc4
        tmplc(22) = py0lc3
        tmplc(23) = py0lc4
        tmplc(24) = z0lc3
        tmplc(25) = z0lc4
        tmplc(26) = pz0lc3
        tmplc(27) = pz0lc4
        !print*,"xlc: ",tmplc(1),tmplc(2),sqsum1local,sqsum2local
        
        !for the calculation of bunched beam current.
!        call MPI_REDUCE(localaa,aa,3,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
!        call MPI_REDUCE(localbb,bb,3,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,0,MPI_COMM_WORLD,ierr)

        call MPI_REDUCE(tmplc,tmpgl,27,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(localmax,glmax,6,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(lcrmax,glrmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innp,npctmin,1,MPI_INTEGER,MPI_MIN,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innp,npctmax,1,MPI_INTEGER,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)

        gam = -this%refptcl(6)
        gambet = sqrt(gam**2-1.0)
        bet = sqrt(gam**2-1.0)/gam
        energy = (gam-1.)*qmc

        if(my_rank.eq.0) then
          x0 = tmpgl(1)*den1
          px0 = tmpgl(2)*den1
          y0 = tmpgl(3)*den1
          py0 = tmpgl(4)*den1
          z0 = tmpgl(5)*den1
          pz0 = tmpgl(6)*den1
          sqx = tmpgl(7)*den1
          sqsum1 = sqx - x0*x0
          sqpx = tmpgl(8)*den1
          sqsum2 = sqpx - px0*px0
          sqy = tmpgl(9)*den1
          sqsum3 = sqy - y0*y0
          sqpy = tmpgl(10)*den1
          sqsum4 = sqpy - py0*py0
          sqz = tmpgl(11)*den1
          sqsum5 = sqz - z0*z0
          sqpz = tmpgl(12)*den1
          sqsum6 = sqpz - pz0*pz0
          xpx = tmpgl(13)*den1 - x0*px0
          ypy = tmpgl(14)*den1 - y0*py0
          zpz = tmpgl(15)*den1 - z0*pz0
          cubx = tmpgl(16)*den1
          fthx = tmpgl(17)*den1
          x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0/3.0)
          x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
          cubpx = tmpgl(18)*den1
          fthpx = tmpgl(19)*den1
          px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0/3.0)
          px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
                 3*px0*px0*px0*px0)))
          cuby = tmpgl(20)*den1
          fthy = tmpgl(21)*den1
          y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0/3.0)
          y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
          cubpy = tmpgl(22)*den1
          fthpy = tmpgl(23)*den1
          py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0/3.0)
          py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
                 3*py0*py0*py0*py0)))
          cubz = tmpgl(24)*den1
          fthz = tmpgl(25)*den1
          z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0/3.0)
          z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
          cubpz = tmpgl(26)*den1
          fthpz = tmpgl(27)*den1
          pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0/3.0)
          pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
                 3*pz0*pz0*pz0*pz0)))
          epsx2 = (sqsum1*sqsum2-xpx*xpx)
          epsy2 = (sqsum3*sqsum4-ypy*ypy)
          epsz2 = (sqsum5*sqsum6-zpz*zpz)
          epx = sqrt(max(epsx2,0.0d0))
          epy = sqrt(max(epsy2,0.0d0))
          epz = sqrt(max(epsz2,0.0d0))
          xrms = sqrt(sqsum1)
          pxrms = sqrt(sqsum2)
          yrms = sqrt(sqsum3)
          pyrms = sqrt(sqsum4)
          zrms = sqrt(sqsum5)
          pzrms = sqrt(sqsum6)
          xpxfac = 0.0
          ypyfac = 0.0
          zpzfac = 0.0
          if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
          if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
          if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)
          write(18,99)z,this%refptcl(5),gam,energy,bet,sqrt(glrmax)*xl
          write(24,100)z,x0*xl,xrms*xl,px0/gambet,pxrms/gambet,-xpx/epx,epx*xl
          write(25,100)z,y0*xl,yrms*xl,py0/gambet,pyrms/gambet,-ypy/epy,epy*xl
          write(26,100)z,z0*xt,zrms*xt,pz0*qmc,pzrms*qmc,-zpz/epz,epz*qmc*xt
          write(27,100)z,glmax(1)*xl,glmax(2)/gambet,glmax(3)*xl,&
                       glmax(4)/gambet,glmax(5)*xt,glmax(6)*qmc
          write(28,101)z,npctmin,npctmax,nptot
          write(29,100)z,x03*xl,px03/gambet,y03*xl,py03/gambet,z03*xt,&
                       pz03*qmc
          write(30,100)z,x04*xl,px04/gambet,y04*xl,py04/gambet,z04*xt,&
                       pz04*qmc
! output bunched beam current.
!          write(31,100)z,aa(1)*den1,aa(2)*den1,aa(3)*den1,bb(1)*den1,&
!                        bb(2)*den1,bb(3)*den1
          write(32,*)z,nptlist(1:nchrg)

!          write(25,100)z,y0*xl,yrms*xl,py0,pyrms,ypy*ypyfac,epy*xl
!          write(24,100)z,x0,xrms,px0,pxrms,-xpx/epx,epx
!          write(24,100)z,x0*xl,xrms*xl,px0,pxrms,-xpx/epx,epx*xl/gambet
!          write(25,100)z,y0,yrms,py0,pyrms,-ypy/epy,epy
!          write(25,100)z,y0*xl,yrms*xl,py0,pyrms,-ypy/epy,epy*xl/gambet
!          write(26,100)z,z0*bet*xl,zrms*bet*xl,pz0,pzrms,-zpz/epz,&
!                       epz*xl/gam**3/bet
!          write(26,100)z,z0*bet*xl,zrms*bet*xl,pz0,pzrms,-zpz/epz,epz*qmc*xt
!          write(26,100)z,z0,zrms,pz0,pzrms,-zpz/epz,epz
!          write(26,100)z,z0*bet*xl,zrms*bet*xl,pz0,pzrms,-zpz/epz,epz*qmc*xt
!          write(26,100)z,z0*xt,zrms*xt,pz0,pzrms,-zpz/epz,epz*xl/(gam**3*bet)
!        write(24,99)z,xrms,pxrms,xpx*xpxfac,epx
!        write(25,99)z,yrms,pyrms,ypy*ypyfac,epy
!        write(26,99)z,zrms,pzrms,zpz*zpzfac,epz
!          write(27,100)z,glmax(1)*xl,glmax(2),glmax(3)*xl,glmax(4), &
!                     glmax(5)*xt,glmax(6)
!        write(27,100)z,glmax(1),glmax(2),glmax(3),glmax(4), &
!                     glmax(5),glmax(6)

          call flush(18)
          call flush(24)
          call flush(25)
          call flush(26)
          call flush(27)
          call flush(28)
          call flush(29)
          call flush(30)
          call flush(32)
        endif

99      format(6(1x,e13.6))
100      format(7(1x,e13.6))
101     format(1x,e13.6,3I10)

        t_diag = t_diag + elapsedtime_Timer(t0)

        end subroutine diagnostic1_Output

        subroutine diagnostic2_Output(this,z,nchrg,nptlist)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(in) :: this
        double precision, intent(in) :: z
        integer, intent(in) :: nchrg
        integer, dimension(nchrg), intent(in) :: nptlist
        integer, parameter :: nbin=100
        integer, dimension(nbin) :: tmpbin,glbin
        integer :: innp,nptot,i,my_rank,ierr,j,npctmin,npctmax,nfreq
        integer :: nii,iitmp
        double precision :: qmc,qmass,xl,xt,tg,pi,t0,gb,bgend,blg
        double precision :: ez1,ezp1,ezpp1,gam,energy,bet,gambet
        double precision:: den1,den2,sqsum1,sqsum2,sqsum3,sqsum4,&
                           epsx2,epsy2
        double precision:: xpx,ypy,epx,epy,xrms,pxrms,yrms,pyrms,&
                           xpxfac,ypyfac
        double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
        double precision:: xpxlocal,ypylocal,zpzlocal
        double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms,pzrms,zpzfac
        double precision:: sqsum5local,sqsum6local,x0lc,x0,px0lc,px0
        double precision:: y0lc,y0,py0lc,py0,z0lc,z0,pz0lc,pz0
        double precision, dimension(6) :: localmax, glmax
        double precision, dimension(15) :: tmplc,tmpgl
        double precision :: lcrmax,glrmax,tmp1,tmp2
        double precision :: f90,f95,f99,ex90,ex95,ex99,ey90,ey95,ey99,&
        ez90,ez95,ez99,r90,r95,r99,rrms,ravg,rrmax,rtmp2,rrmaxlc,ravg2,&
        ravglc,ravg2lc
        double precision, dimension(3) :: Ealpha,Ebeta,Egamma
        double precision :: epsmxlc,epsmylc,epsmzlc,xtmp,pxtmp,ytmp,pytmp,&
                   ztmp,pztmp,epsmx,epsmy,epsmz,eps,hxeps,hyeps,hzeps,hreps
        double precision,allocatable,dimension(:) :: epsiontmp

        call starttime_Timer(t0)

        qmc = this%Mass/1.0e6
        qmass =this%Mass
        xl = Scxl
        xt = Rad2deg
        tg = this%refptcl(5)
        pi = 2*asin(1.0)

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

        innp = this%Nptlocal
        nptot = this%Npt

        den1 = 1.0/float(nptot)
        den2 = den1*den1
        x0lc = 0.0
        px0lc = 0.0
        y0lc = 0.0
        py0lc = 0.0
        z0lc = 0.0
        pz0lc = 0.0
        sqsum1local = 0.0
        sqsum2local = 0.0
        sqsum3local = 0.0
        sqsum4local = 0.0
        sqsum5local = 0.0
        sqsum6local = 0.0
        xpxlocal = 0.0
        ypylocal = 0.0
        zpzlocal = 0.0
        ! for cache optimization.
        if(innp.ne.0) then
          do i = 1, 6
            localmax(i) = abs(this%Pts1(i,1))
          enddo
          lcrmax = this%Pts1(1,1)**2+this%Pts1(3,1)**2
        else
          do i = 1, 6
            localmax(i) = 0.0
          enddo
          lcrmax = 0.0
        endif
        do i = 1, innp
          x0lc = x0lc + this%Pts1(1,i)
          sqsum1local = sqsum1local + this%Pts1(1,i)*this%Pts1(1,i)
!          tmp1 = this%Pts1(2,i)-(xl*xl*ezp1*this%Pts1(1,i)*&
!            sin(nfreq*(this%Pts1(5,i)+tg)+theta*pi/180)/(2*qmass*nfreq))
          tmp1 = this%Pts1(2,i)
          xpxlocal = xpxlocal + this%Pts1(1,i)*tmp1
          px0lc = px0lc + tmp1
          sqsum2local = sqsum2local + tmp1*tmp1
          y0lc = y0lc + this%Pts1(3,i)
          sqsum3local = sqsum3local + this%Pts1(3,i)*this%Pts1(3,i)
!          tmp2 = this%Pts1(4,i)-(xl*xl*ezp1*this%Pts1(3,i)*&
!            sin(nfreq*(this%Pts1(5,i)+tg)+theta*pi/180)/(2*qmass*nfreq))
          tmp2 = this%Pts1(4,i)
          ypylocal = ypylocal + this%Pts1(3,i)*tmp2
          py0lc = py0lc + tmp2
          sqsum4local = sqsum4local + tmp2*tmp2
          z0lc = z0lc + this%Pts1(5,i)
          sqsum5local = sqsum5local + this%Pts1(5,i)*this%Pts1(5,i)
          zpzlocal = zpzlocal + this%Pts1(5,i)*this%Pts1(6,i)
          pz0lc = pz0lc + this%Pts1(6,i)
          sqsum6local = sqsum6local + this%Pts1(6,i)*this%Pts1(6,i)
          do j = 1, 6
            if(localmax(j).lt.abs(this%Pts1(j,i))) then
               localmax(j) = abs(this%Pts1(j,i))
            endif
          enddo
          if(lcrmax.lt.(this%Pts1(1,i)**2+this%Pts1(3,i)**2)) then
            lcrmax = this%Pts1(1,i)**2 + this%Pts1(3,i)**2
          endif
        enddo

        tmplc(1) = x0lc
        tmplc(2) = px0lc
        tmplc(3) = y0lc
        tmplc(4) = py0lc
        tmplc(5) = z0lc
        tmplc(6) = pz0lc
        tmplc(7) = sqsum1local
        tmplc(8) = sqsum2local
        tmplc(9) = sqsum3local
        tmplc(10) = sqsum4local
        tmplc(11) = sqsum5local
        tmplc(12) = sqsum6local
        tmplc(13) = xpxlocal
        tmplc(14) = ypylocal
        tmplc(15) = zpzlocal
        
        call MPI_ALLREDUCE(tmplc,tmpgl,15,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(lcrmax,glrmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(localmax,glmax,6,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innp,npctmin,1,MPI_INTEGER,MPI_MIN,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innp,npctmax,1,MPI_INTEGER,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)

        gam = -this%refptcl(6)
        bet = sqrt(gam**2-1.0)/gam
        gambet = sqrt(gam**2-1.0)
        energy = (gam-1.)*qmc

        x0 = tmpgl(1)
        px0 = tmpgl(2)
        y0 = tmpgl(3)
        py0 = tmpgl(4)
        z0 = tmpgl(5)
        pz0 = tmpgl(6)
        sqsum1 = tmpgl(7) - x0*x0*den1
        sqsum2 = tmpgl(8) - px0*px0*den1
        sqsum3 = tmpgl(9) - y0*y0*den1
        sqsum4 = tmpgl(10) - py0*py0*den1
        sqsum5 = tmpgl(11) - z0*z0*den1
        sqsum6 = tmpgl(12) - pz0*pz0*den1
        xpx = tmpgl(13) - x0*px0*den1
        ypy = tmpgl(14) - y0*py0*den1
        zpz = tmpgl(15) - z0*pz0*den1
        epsx2 = (sqsum1*sqsum2-xpx*xpx)*den2
        epsy2 = (sqsum3*sqsum4-ypy*ypy)*den2
        epsz2 = (sqsum5*sqsum6-zpz*zpz)*den2
        epx = sqrt(max(epsx2,0.0d0))
        epy = sqrt(max(epsy2,0.0d0))
        epz = sqrt(max(epsz2,0.0d0))
        xrms = sqrt(sqsum1*den1)
        pxrms = sqrt(sqsum2*den1)
        yrms = sqrt(sqsum3*den1)
        pyrms = sqrt(sqsum4*den1)
        zrms = sqrt(sqsum5*den1)
        pzrms = sqrt(sqsum6*den1)
        xpx = xpx*den1
        ypy = ypy*den1
        zpz = zpz*den1
        xpxfac = 0.0
        ypyfac = 0.0
        zpzfac = 0.0
        if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
        if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
        if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)
        x0 = x0*den1
        px0 = px0*den1
        y0 = y0*den1
        py0 = py0*den1
        z0 = z0*den1
        pz0 = pz0*den1

        Ealpha(1) = -xpx/epx
        Ealpha(2) = -ypy/epy
        Ealpha(3) = -zpz/epz
        Ebeta(1) = xrms*xrms*gambet/epx
        Ebeta(2) = yrms*yrms*gambet/epy
        Ebeta(3) = zrms*zrms*bet*bet/(epz/gam**3/bet)
        Egamma(:) = (1.0+Ealpha(:)*Ealpha(:))/Ebeta(:)

        allocate(epsiontmp(innp))
        epsmxlc = -1.0e10
        do i = 1, innp
          xtmp = this%Pts1(1,i) - x0
!          tmp1 = this%Pts1(2,i)-(xl*xl*ezp1*this%Pts1(1,i)*&
!            sin(nfreq*(this%Pts1(5,i)+tg)+theta*pi/180)/(2*qmass*nfreq))
          tmp1 = this%Pts1(2,i)
          pxtmp = (tmp1 - px0)/gambet
          epsiontmp(i)=Egamma(1)*xtmp*xtmp+2*Ealpha(1)*xtmp*pxtmp+&
                     Ebeta(1)*pxtmp*pxtmp
          if(epsmxlc.le.epsiontmp(i)) epsmxlc = epsiontmp(i)
        enddo

        call MPI_ALLREDUCE(epsmxlc,epsmx,1,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,MPI_COMM_WORLD,ierr)

        eps = 1.0e-8
        hxeps = (epsmx+eps)/nbin

        do i = 1, nbin
          tmpbin(i) = 0
        enddo

        do i = 1, innp
          iitmp = epsiontmp(i)/hxeps
          nii = iitmp+1
          tmpbin(nii) = tmpbin(nii) + 1 
        enddo
        call MPI_REDUCE(tmpbin,glbin,nbin,MPI_INTEGER,&
                           MPI_SUM,0,MPI_COMM_WORLD,ierr)

        f90 = 0.9*nptot
        f95 = 0.95*nptot
        f99 = 0.99*nptot
        if(my_rank.eq.0) then
          do i = 2, nbin
            glbin(i) = glbin(i) + glbin(i-1)
          enddo
          do i = 1, nbin
            if(glbin(i).gt.f90) then
              ex90 = ((f90 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hxeps+&
                    hxeps*(i-1) 
              exit
            endif
          enddo
          !print*,"i1: ",i,nbin,glbin(i-1),glbin(i),f90,f95,f99
          do i = 1, nbin
            if(glbin(i).gt.f95) then
              ex95 = ((f95 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hxeps+&
                    hxeps*(i-1) 
              exit
            endif
          enddo
          !print*,"i2: ",i,nbin,glbin(i-1),glbin(i)
          do i =1, nbin
            if(glbin(i).gt.f99) then
              ex99 = ((f99 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hxeps+&
                    hxeps*(i-1) 
              exit
            endif
          enddo
          !print*,"i3: ",i,nbin,glbin(i-1),glbin(i)
          !print*,"pass hxeps: ",glbin(nbin),nptot,epsmx,ex90,ex95,ex99
          ex90 = ex90*gambet
          ex95 = ex95*gambet
          ex99 = ex99*gambet
        endif

        epsmylc = -1.0e10
        do i = 1, innp
          ytmp = this%Pts1(3,i) - y0
!          tmp2 = this%Pts1(4,i)-(xl*xl*ezp1*this%Pts1(3,i)*&
!            sin(nfreq*(this%Pts1(5,i)+tg)+theta*pi/180)/(2*qmass*nfreq))
          tmp2 = this%Pts1(4,i)
          pytmp = (tmp2 - py0)/gambet
          epsiontmp(i)=Egamma(2)*ytmp*ytmp+2*Ealpha(2)*ytmp*pytmp+&
                     Ebeta(2)*pytmp*pytmp
          if(epsmylc.le.epsiontmp(i)) epsmylc = epsiontmp(i)
        enddo
        call MPI_ALLREDUCE(epsmylc,epsmy,1,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,MPI_COMM_WORLD,ierr)

        hyeps = (epsmy+eps)/nbin
        do i = 1, nbin
          tmpbin(i) = 0
        enddo

        do i = 1, innp
          iitmp = epsiontmp(i)/hyeps
          nii = iitmp+1
          tmpbin(nii) = tmpbin(nii) + 1
        enddo
        call MPI_REDUCE(tmpbin,glbin,nbin,MPI_INTEGER,&
                           MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 2, nbin
            glbin(i) = glbin(i) + glbin(i-1)
          enddo
          do i = 1, nbin
            if(glbin(i).gt.f90) then
              ey90 = ((f90 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hyeps+&
                    hyeps*(i-1)
              exit
            endif
          enddo
          !print*,"i4: ",i,nbin,glbin(i-1),glbin(i)
          do i = 1, nbin
            if(glbin(i).gt.f95) then
              ey95 = ((f95 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hyeps+&
                    hyeps*(i-1)
              exit
            endif
          enddo
          !print*,"i5: ",i,nbin,glbin(i-1),glbin(i)
          do i = 1, nbin
            if(glbin(i).gt.f99) then
              ey99 = ((f99 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hyeps+&
                    hyeps*(i-1)
              exit
            endif
          enddo
          !print*,"i6: ",i,nbin,glbin(i-1),glbin(i)
          !print*,"pass hyeps: ",glbin(nbin),nptot,epsmy,ey90,ey95,ey99
          ey90 = ey90*gambet
          ey95 = ey95*gambet
          ey99 = ey99*gambet
        endif

        epsmzlc = -1.0e10
        do i = 1, innp
          ztmp = (this%Pts1(5,i) - z0)*bet
          pztmp = (this%Pts1(6,i) - pz0)/gam**3/bet**2
          epsiontmp(i)=Egamma(3)*ztmp*ztmp+2*Ealpha(3)*ztmp*pztmp+&
                     Ebeta(3)*pztmp*pztmp
          if(epsmzlc.le.epsiontmp(i)) epsmzlc = epsiontmp(i)
        enddo
        call MPI_ALLREDUCE(epsmzlc,epsmz,1,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,MPI_COMM_WORLD,ierr)

        hzeps = (epsmz+eps)/nbin
        do i = 1, nbin
          tmpbin(i) = 0
        enddo

        do i = 1, innp
          iitmp = epsiontmp(i)/hzeps
          nii = iitmp+1
          tmpbin(nii) = tmpbin(nii) + 1
        enddo
        call MPI_REDUCE(tmpbin,glbin,nbin,MPI_INTEGER,&
                           MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 2, nbin
            glbin(i) = glbin(i) + glbin(i-1)
          enddo
          do i = 1, nbin
            if(glbin(i).gt.f90) then
              ez90 = ((f90 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hzeps+&
                    hzeps*(i-1)
              exit
            endif
          enddo
          !print*,"i7: ",i,nbin,glbin(i-1),glbin(i)
          do i = 1, nbin
            if(glbin(i).gt.f95) then
              ez95 = ((f95 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hzeps+&
                    hzeps*(i-1)
              exit
            endif
          enddo
          !print*,"i8: ",i,nbin,glbin(i-1),glbin(i)
          do i = 1, nbin
            if(glbin(i).gt.f99) then
              ez99 = ((f99 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hzeps+&
                    hzeps*(i-1)
              exit
            endif
          enddo
          !print*,"i9: ",i,nbin,glbin(i-1),glbin(i)
          !print*,"pass hzeps: ",glbin(nbin),nptot,epsmz,ez90,ez95,ez99
          ez90 = ez90*gam**3*bet
          ez95 = ez95*gam**3*bet
          ez99 = ez99*gam**3*bet
        endif

        deallocate(epsiontmp)

        ravglc = 0.0
        ravg2lc = 0.0
        rrmaxlc = -1.0e10
        do i = 1, innp
          xtmp = this%Pts1(1,i) - x0
          ytmp = this%Pts1(3,i) - y0
          rtmp2 = xtmp*xtmp + ytmp*ytmp
          ravglc = ravglc + sqrt(rtmp2)
          ravg2lc = ravg2lc + rtmp2
          if(rrmaxlc.le.rtmp2) rrmaxlc = rtmp2
        enddo
        call MPI_REDUCE(ravglc,ravg,1,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(ravg2lc,ravg2,1,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(rrmaxlc,rrmax,1,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,MPI_COMM_WORLD,ierr)

        hreps = (sqrt(rrmax)+eps)/nbin
        do i = 1, nbin
          tmpbin(i) = 0
        enddo

        do i = 1, innp
          iitmp = sqrt((this%Pts1(1,i)-x0)**2+(this%Pts1(3,i)-y0)**2)/hreps
          nii = iitmp+1
          tmpbin(nii) = tmpbin(nii) + 1
        enddo
        call MPI_REDUCE(tmpbin,glbin,nbin,MPI_INTEGER,&
                           MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 2, nbin
            glbin(i) = glbin(i) + glbin(i-1)
          enddo
          do i = 1, nbin
            if(glbin(i).gt.f90) then
              r90 = ((f90 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hreps+&
                    hreps*(i-1)
              exit
            endif
          enddo
          !print*,"i10: ",i,nbin,glbin(i-1),glbin(i)
          do i = 1, nbin
            if(glbin(i).gt.f95) then
              r95 = ((f95 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hreps+&
                    hreps*(i-1)
              exit
            endif
          enddo
          !print*,"i11: ",i,nbin,glbin(i-1),glbin(i)
          do i = 1, nbin
            if(glbin(i).gt.f99) then
              r99 = ((f99 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hreps+&
                    hreps*(i-1)
              exit
            endif
          enddo
          !print*,"i12: ",i,nbin,glbin(i-1),glbin(i)
          !print*,"pass hreps: ",glbin(nbin),nptot,sqrt(rrmax),r90,r95,r99
        endif

        if(my_rank.eq.0) then
          write(18,99)z,this%refptcl(5),gam,energy,bet,sqrt(glrmax)*xl
          write(24,100)z,x0*xl,xrms*xl,px0/gam/bet,pxrms/gam/bet,-xpx/epx,epx*xl,ex90*xl,&
                       ex95*xl,ex99*xl
          write(25,100)z,y0*xl,yrms*xl,py0/gam/bet,pyrms/gam/bet,-ypy/epy,epy*xl,ey90*xl,&
                       ey95*xl,ey99*xl
          write(26,100)z,z0*xt,zrms*xt,pz0*qmc,pzrms*qmc,-zpz/epz,epz*qmc*xt,&
                       ez90*qmc*xt,ez95*qmc*xt,ez99*qmc*xt
          !write(26,100)z,z0*bet*xl,zrms*bet*xl,pz0,pzrms,-zpz/epz,epz*qmc*xt,&
          !             ez90*qmc*xt,ez95*qmc*xt,ez99*qmc*xt
          write(27,102)z,glmax(1)*xl,glmax(2),glmax(3)*xl,glmax(4), &
                     glmax(5)*xt,glmax(6)
          write(28,101)z,npctmin,npctmax,nptot
          ravg = ravg/nptot
          rrms = sqrt(ravg2/nptot - ravg*ravg)
          write(29,102)z,ravg*xl,rrms*xl,r90*xl,r95*xl,r99*xl,sqrt(rrmax)*xl
          write(32,*)z,nptlist(1:nchrg)

          call flush(18)
          call flush(24)
          call flush(25)
          call flush(26)
          call flush(27)
          call flush(28)
          call flush(29)
          call flush(32)
!          write(24,100)z,x0*xl,xrms*xl,px0,pxrms,xpx*xpxfac,epx*xl
!          write(25,100)z,y0*xl,yrms*xl,py0,pyrms,ypy*ypyfac,epy*xl
!          write(26,100)z,z0*xt,zrms*xt,pz0,pzrms,zpz*zpzfac,epz*qmc*xt
!          write(24,100)z,x0*xl,xrms*xl,px0,pxrms,-xpx/epx,epx*xl/gambet
!          write(25,100)z,y0*xl,yrms*xl,py0,pyrms,-ypy/epy,epy*xl/gambet
!          write(26,100)z,z0*bet*xl,zrms*bet*xl,pz0,pzrms,-zpz/epz,&
!                       epz*xl/gam**3/bet
!          write(26,100)z,z0*xt,zrms*xt,pz0,pzrms,-zpz/epz,epz*qmc*xt
!          write(26,100)z,z0*xt,zrms*xt,pz0,pzrms,-zpz/epz,epz*xl/(gam**3*bet)
!        write(24,99)z,xrms,pxrms,xpx*xpxfac,epx
!        write(25,99)z,yrms,pyrms,ypy*ypyfac,epy
!        write(26,99)z,zrms,pzrms,zpz*zpzfac,epz
!          write(27,100)z,glmax(1)*xl,glmax(2),glmax(3)*xl,glmax(4), &
!                     glmax(5)*bet*xl,glmax(6)
!        write(27,100)z,glmax(1),glmax(2),glmax(3),glmax(4), &
!                     glmax(5),glmax(6)
        endif

99      format(6(1x,e13.6))
100      format(10(1x,e13.6))
101     format(1x,e13.6,3I10)
102      format(7(1x,e13.6))

        t_diag = t_diag + elapsedtime_Timer(t0)

        end subroutine diagnostic2_Output

        subroutine phase_Output(nfile,this,samplePeriod)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile
        type (BeamBunch), intent(in) :: this
        integer :: samplePeriod
        integer :: np,my_rank,ierr
        integer status(MPI_STATUS_SIZE)
        integer :: i,j,sixnpt,mnpt
        integer, allocatable, dimension(:) :: nptlist
        double precision, allocatable,dimension(:,:) :: recvbuf

        if (samplePeriod .eq. 0) then
           samplePeriod = 1
        endif

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)
        call MPI_ALLREDUCE(this%Nptlocal,mnpt,1,MPI_INTEGER,MPI_MAX,&
                        MPI_COMM_WORLD,ierr)

        allocate(nptlist(0:np-1))
        nptlist = 0
        allocate(recvbuf(9,mnpt))
        sixnpt = 9*this%Nptlocal

        call MPI_GATHER(this%Nptlocal,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        nptlist = 9*nptlist

        if(my_rank.eq.0) then
          open(nfile,status='unknown')
          do i = 1, this%Nptlocal,samplePeriod
            write(nfile,100)this%Pts1(1,i),this%Pts1(2,i),this%Pts1(3,i),&
                            this%Pts1(4,i),this%Pts1(5,i),this%Pts1(6,i),&
                            this%Pts1(7,i),this%Pts1(8,i),this%Pts1(9,i)
          enddo
          do i = 1, np-1
            call MPI_RECV(recvbuf(1,1),nptlist(i),MPI_DOUBLE_PRECISION,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
        
            do j = 1, nptlist(i)/9,samplePeriod
              write(nfile,100)recvbuf(1,j),recvbuf(2,j),recvbuf(3,j),&
                              recvbuf(4,j),recvbuf(5,j),recvbuf(6,j),&
                              recvbuf(7,j),recvbuf(8,j),recvbuf(9,j)
            enddo
          enddo
          close(nfile)
        else
          call MPI_SEND(this%Pts1(1,1),sixnpt,MPI_DOUBLE_PRECISION,0,1,&
                        MPI_COMM_WORLD,ierr)
        endif

100     format(9(1x,e14.7))

        deallocate(nptlist)
        deallocate(recvbuf)

        end subroutine phase_Output

        ! Terminate MPI
        subroutine end_Output(time)
        implicit none
        include 'mpif.h'
        double precision, intent(inout) :: time
        double precision :: endtime, mtime
        integer :: my_rank,ierr

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        endtime = MPI_WTIME()
        time = endtime - time
        call MPI_REDUCE(time,mtime,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          close(24)
          close(25)
          print*,"time: ",mtime
        endif

        !for measurement of memory
        !call system_stats()

        call MPI_Finalize(ierr)

        end subroutine end_Output

      end module Outputclass
