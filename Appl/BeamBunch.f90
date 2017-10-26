!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! BeamBunchclass: Charged beam bunch class in Beam module of APPLICATION 
!                 layer.
! Version: 1.0
! Author: Ji Qiang 
! Description: This class defines the charged particle beam bunch 
!              information in the accelerator.
! Comments: 1) We have added the 3 attributes to the particle array:
!           x,px,y,py,t,pt,charge/mass,charge weight,id. We have moved
!           the charge*curr/freq/Ntot into the charge density calculation,
!           which is represented by the "charge weight" of a particle.
!           2) The map integrator does not work for multiple species, only
!              the Lorenze integrator works for the multiple species.
!----------------------------------------------------------------
      module BeamBunchclass
        use CompDomclass
        use Pgrid2dclass
        use BeamLineElemclass
        use Timerclass
        use PhysConstclass
        use NumConstclass
        use TPSAMod
        type BeamBunch
!          private
          !beam freq, current, part. mass and charge.
          double precision :: Current,Mass,Charge
          !# of total global macroparticles and local particles
          integer :: Nptlocal
          integer :: Npt
          !particles type one.
          double precision, pointer, dimension(:,:) :: Pts1
          !reference particle
          double precision, dimension(6) :: refptcl
        end type BeamBunch
        interface map1_BeamBunch
          module procedure drift1_BeamBunch,tpsaMap_BeamBunch
        end interface
        interface sextupole_BeamBunch
          module procedure sextupole_Ptc_BeamBunch,sextupole_tpsa_BeamBunch
        end interface
      contains
        subroutine construct_BeamBunch(this,incurr,inkin,inmass,incharge,innp,&
                                       phasini)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision, intent(in) :: incurr,inkin,inmass,&
                                        incharge,phasini
        integer, intent(in) :: innp
        integer :: myid, myidx, myidy,comm2d,commrow,commcol,ierr
        integer :: nptot,nprocrow,nproccol
   
        this%Current = incurr
        this%Mass = inmass
        this%Charge = incharge
        this%Npt = innp

        this%refptcl = phasini
        this%refptcl(6) = -(inkin/this%Mass + 1.0)

        end subroutine construct_BeamBunch


        ! Drift beam half step using linear map for external field.
        subroutine drift1_BeamBunch(this,beamln,z,tau,bitype,nseg,nst,ihlf)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        type (BeamLineElem), intent(inout) :: beamln
        double precision, intent(inout) :: z
        double precision, intent (in) :: tau
        integer, intent(in) :: bitype,nseg,nst
        integer, intent(inout) :: ihlf
        double precision, dimension(6) :: temp,temp1,temp2
        double precision, dimension(6,6) :: xm
        double precision :: t0
        integer :: i,j,k
        real*8 :: beta0,gam0,gambet0,tmppx,tmppy,tmppt,tmph
!for multple reference for rf
        real*8 :: phmin,phmax,phmintmp,phmaxtmp
        integer, parameter :: nslice = 20
        real*8, dimension(6,6,nslice+1) :: xmnsl
        real*8, dimension(nslice+1) :: rf5,rf6
        real*8 :: ht,tt,tt0,tmin,tmp15,tmp16
        real*8, dimension(6) :: tmp2
        real*8 :: delta
        integer :: it,isl,ierr,it1,it2
        real*8 :: qmass,hh,vacAcc,del1,del2,del3,tt1,tt2,rr
!fore ideal RF cavity model
        integer :: idealrf
        real*8 :: vtmp,vtmpgrad,phi0lc,phi,gam,gambet,gambeto,&
                  gambetz,invfocus,dgam,fact,harm,dgambet,cosphi
        real*8, dimension(12) :: drange

        call starttime_Timer(t0)

        qmass = this%Charge/this%Mass 

!        call random_number(rr)

        if(bitype.eq.99999) then
          gam0 = -this%refptcl(6)
          gambet0 = sqrt(gam0**2-1.0d0)
          beta0 = gambet0/gam0
          do i = 1, this%Nptlocal
            !tmppx = this%Pts1(2,i)/gambet0
            !tmppy = this%Pts1(4,i)/gambet0
            !tmppt = this%Pts1(6,i)/gambet0
            !tmph = sqrt((tmppt-1.d0/beta0)**2-(1./gambet0)**2-tmppx**2-tmppy**2)
            !this%Pts1(1,i) = this%Pts1(1,i)+tmppx*tau/tmph/Scxl 
            !this%Pts1(3,i) = this%Pts1(3,i)+tmppy*tau/tmph/Scxl 
            !this%Pts1(5,i) = this%Pts1(5,i)+(1./beta0+(tmppt-1./beta0)/tmph)*&
            !                 tau/Scxl
            tmppx = this%Pts1(2,i)
            tmppy = this%Pts1(4,i)
            tmppt = this%Pts1(6,i)
            tmph = sqrt((tmppt-gam0)**2-1-tmppx**2-tmppy**2)
            this%Pts1(1,i) = this%Pts1(1,i)+tmppx*tau/tmph/Scxl 
            this%Pts1(3,i) = this%Pts1(3,i)+tmppy*tau/tmph/Scxl 
            this%Pts1(5,i) = this%Pts1(5,i)-(1./beta0+(tmppt-gam0)/tmph)*&
                             tau/Scxl
          enddo
          this%refptcl(5) = this%refptcl(5) + tau/(Scxl*beta0)
          !artifically linear ramping energy for modulation study
          !if(z.le.15.0d0) then
          !  this%refptcl(6) = this%refptcl(6) - tau*5.0880626
          !endif
          !for artificial vacuum acceleration
          vacAcc = beamln%pdrift%Param(3)
          if(abs(vacAcc) > 1.0d-6) then 
            this%refptcl(6) = this%refptcl(6) - tau*vacAcc*qmass
          endif
        else if(bitype.eq.1) then
          print*,"not available yet"
        else
           !standard rf linear map
           call maplinear_BeamLineElem(beamln,z,tau,xm,this%refptcl,this%Charge,&
                       this%Mass)

           do i = 1, this%Nptlocal
            do j = 1, 6
              temp(j) = 0.0
              do k = 1, 6
                temp(j) = temp(j) + this%Pts1(k,i)*xm(j,k)
              enddo
            enddo
            do j = 1,6
              this%Pts1(j,i) = temp(j)
            enddo
           enddo

          !artifically linear ramping energy for modulation study
          !if(z.lt.15.0d0) then
          !  this%refptcl(6) = this%refptcl(6) - tau*5.0880626
          !endif

        endif

        z=z+tau

        !introduce an instant chirp at the end of the drift.
        !if((abs(z-15.0d0).lt.1.0e-2).and.(bitype.eq.0)) then
        !  gam0 = -this%refptcl(6)
        !  hh = -9.0
        !  do j = 1, this%Nptlocal
        !    this%Pts1(6,j) = this%Pts1(6,j) + gam0*hh*Scxl*this%Pts1(5,j)
        !  enddo
        !endif

        ihlf = ihlf + 1

        t_map1 = t_map1 + elapsedtime_Timer(t0)

        end subroutine drift1_BeamBunch

        subroutine tpsaMap_BeamBunch(this,beamln,z,tau,bitype,nseg,nst,ihlf, tpsaPtc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        type (BeamLineElem), intent(inout) :: beamln
        type (dctps), dimension(6) , intent(inout) ::tpsaPtc
        type (dctps), dimension(6) ::tpsaPtcTemp
        double precision, intent(inout) :: z
        double precision, intent (in) :: tau
        integer, intent(in) :: bitype,nseg,nst
        integer, intent(inout) :: ihlf
        double precision, dimension(6) :: temp,temp1,temp2
        double precision, dimension(6,6) :: xm
        double precision :: t0
        integer :: i,j,k
        real*8 :: beta0,gam0,gambet0,tmppx,tmppy,tmppt,tmph
!for multple reference for rf
        real*8 :: phmin,phmax,phmintmp,phmaxtmp
        integer, parameter :: nslice = 20
        real*8, dimension(6,6,nslice+1) :: xmnsl
        real*8, dimension(nslice+1) :: rf5,rf6
        real*8 :: ht,tt,tt0,tmin,tmp15,tmp16
        real*8, dimension(6) :: tmp2
        real*8 :: delta
        integer :: it,isl,ierr,it1,it2
        real*8 :: qmass,hh,vacAcc,del1,del2,del3,tt1,tt2,rr
!fore ideal RF cavity model
        integer :: idealrf
        real*8 :: vtmp,vtmpgrad,phi0lc,phi,gam,gambet,gambeto,&
                  gambetz,invfocus,dgam,fact,harm,dgambet,cosphi
        real*8, dimension(12) :: drange

        call starttime_Timer(t0)

        qmass = this%Charge/this%Mass 

        !standard rf linear map
        call maplinear_BeamLineElem(beamln,z,tau,xm,this%refptcl,this%Charge,&
          this%Mass)

          do j = 1, 6
            tpsaPtcTemp(j) = 0.d0
            do k = 1, 6
                tpsaPtcTemp(j) = tpsaPtcTemp(j) + tpsaPtc(k)*xm(j,k)
            enddo
          enddo
          do j = 1,6
            tpsaPtc(j) = tpsaPtcTemp(j)
          enddo

        z=z+tau
        ihlf = ihlf + 1
        t_map1 = t_map1 + elapsedtime_Timer(t0)

        end subroutine tpsaMap_BeamBunch


        subroutine sextupole_Ptc_BeamBunch(ptc,ksext,gambet,scxl)
        implicit none        
        double precision, dimension(:,:), intent(inout) :: ptc
        double precision, intent(in)   :: ksext, gambet,scxl
        double precision :: scxl2,factor
        integer :: i,nPtc

        nPtc = size(ptc,2)

        scxl2=scxl*scxl
        factor = 1.0*gambet * ksext/2.0 * scxl2

        do i = 1, nPtc
          ptc(2,i) = ptc(2,i) - (ptc(1,i)**2 - ptc(3,i)**2) * factor
          ptc(4,i) = ptc(4,i) - (ptc(1,i)    * ptc(3,i)   ) * factor
        enddo

        end subroutine sextupole_Ptc_BeamBunch


        subroutine sextupole_tpsa_BeamBunch(ptc,ksext,gambet,scxl)
        implicit none        
        type (dctps), dimension(6) , intent(inout) ::ptc
        double precision, intent(in)   :: ksext, gambet,scxl
        double precision :: scxl2,factor
        integer :: i,nPtc


        scxl2=scxl*scxl
        factor = 1.0*gambet * ksext/2.0 * scxl2

          ptc(2) = ptc(2) - (ptc(1)*ptc(1) - ptc(3)*ptc(3)) * factor
          ptc(4) = ptc(4) - (ptc(1)    * ptc(3)   ) * factor
        end subroutine sextupole_tpsa_BeamBunch


        !from z to t beam frame 0th order transformation.
        subroutine conv0thB_BeamBunch(this,tau,nplc,nptot,ptrange,&
                                     Flagbc,perd,xrad,yrad)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision, intent(in) :: tau,xrad,yrad,perd
        integer, intent(in) :: Flagbc
        double precision, dimension(6), intent(out) :: ptrange
        integer, intent(inout) :: nplc
        integer, intent(inout) :: nptot
        integer :: i
        double precision :: xk,xl
        double precision :: pi,gam,bet,bbyk,rcpgammai,betai,rad
        double precision :: twopi,radtest,tmp0,tmpx,tmpy,tmp5,halfperd
        integer :: ilost,i0,ierr,ntmp5
        real*8 :: fnplc,fnptot,gambet

        pi = 2.0d0*asin(1.0d0)
        twopi = 2.0d0*pi
        xl = Scxl
        xk = 1/xl

        gam = -this%refptcl(6)
        bet = sqrt(gam**2-1.0d0)/gam
        gambet = gam*bet
        bbyk = bet/xk
        rad = (xrad+yrad)/2 
      
        do i = 1, 3
          ptrange(2*i-1) = 1.0e20
          ptrange(2*i) = -1.0e20
        enddo

        halfperd = 0.5d0*perd

        if(Flagbc.eq.1) then ! open 3D
          ilost = 0
          do i0 = 1, this%Nptlocal
            i = i0 - ilost
            !if(abs(this%Pts1(5,i0)).ge.pi) then
            !  ilost = ilost + 1
            !  goto 100
            !endif
            ! The following steps go from z to t frame.
            ! 2) zeroth algorithm to transfer from z to t frame.
!            rcpgammai = 1.0d0/(-this%Pts1(6,i0)+gam)
!            betai = sqrt(1.0d0-rcpgammai*rcpgammai*(1+this%Pts1(2,i0)**2+ &
!                        this%Pts1(4,i0)**2) )
            tmp0 = this%Pts1(1,i0)*xl
            this%Pts1(1,i) = this%Pts1(1,i0)*xl

            if(ptrange(1)>tmp0) then
              ptrange(1) = tmp0
            endif
            if(ptrange(2)<tmp0) then
              ptrange(2) = tmp0
            endif

            this%Pts1(2,i) = this%Pts1(2,i0)
            tmp0 = this%Pts1(3,i0)*xl
            this%Pts1(3,i) = this%Pts1(3,i0)*xl

            if(ptrange(3)>tmp0) then
              ptrange(3) = tmp0
            endif
            !tmp0 = max(tmpy,this%Pts1(3,i))
            if(ptrange(4)<tmp0) then
              ptrange(4) = tmp0
            endif
 
            this%Pts1(4,i) = this%Pts1(4,i0)
            this%Pts1(5,i) = -gambet*this%Pts1(5,i0)*xl

            if(ptrange(5)>this%Pts1(5,i)) then
              ptrange(5) = this%Pts1(5,i)
            endif
            if(ptrange(6)<this%Pts1(5,i)) then
              ptrange(6) = this%Pts1(5,i)
            endif

            this%Pts1(6,i) = this%Pts1(6,i0)
            radtest = sqrt(this%Pts1(1,i)*this%Pts1(1,i)+ &
                      this%Pts1(3,i)*this%Pts1(3,i))
            if(radtest.ge.rad) then
              ilost = ilost + 1
            endif
            this%Pts1(7,i) = this%Pts1(7,i0)
            this%Pts1(8,i) = this%Pts1(8,i0)
            this%Pts1(9,i) = this%Pts1(9,i0)
100         continue
          enddo
        else if(Flagbc.eq.2) then ! open 2D, 1D z periodic
          ilost = 0
          do i0 = 1, this%Nptlocal
            i = i0 - ilost
            ! The following steps go from z to t frame.
            ! 2) Linear algorithm to transfer from z to t frame.
!            rcpgammai = 1.0d0/(-this%Pts1(6,i0)+gam)
!            betai = sqrt(1.0d0-rcpgammai*rcpgammai*(1+this%Pts1(2,i0)**2+ &
!                       this%Pts1(4,i0)**2) )

!            ntmp5 = this%Pts1(5,i0)/pi
!            tmp5 = this%Pts1(5,i0) - ntmp5*pi 
!            if(mod(ntmp5,2).eq.0) then
!              this%Pts1(5,i0) = tmp5
!            else
!              if(tmp5.gt.0.0d0) then
!                this%Pts1(5,i0) = tmp5 - pi
!              else
!                this%Pts1(5,i0) = tmp5 + pi
!              endif 
!            endif
 
            this%Pts1(4,i) = this%Pts1(4,i0)
            this%Pts1(5,i) = -gambet*this%Pts1(5,i0)*xl

            ntmp5 = this%Pts1(5,i)/halfperd
            tmp5 = this%Pts1(5,i) - ntmp5*halfperd 
            this%Pts1(5,i) = tmp5 - mod(ntmp5,2)*halfperd
            tmp5 = -this%Pts1(5,i)/(gambet*xl)

            this%Pts1(6,i) = this%Pts1(6,i0)
            tmp0 = this%Pts1(1,i0)*xl
! for perd bunch
            this%Pts1(1,i) = this%Pts1(1,i0)*xl

            if(ptrange(1)>tmp0) then
               ptrange(1) = tmp0
            endif
            if(ptrange(2)<tmp0) then
               ptrange(2) = tmp0
            endif

            this%Pts1(2,i) = this%Pts1(2,i0)
            tmp0 = this%Pts1(3,i0)*xl
! for perd bunch
            this%Pts1(3,i) = this%Pts1(3,i0)*xl

            if(ptrange(3)>tmp0) then
              ptrange(3) = tmp0
            endif
            if(ptrange(4)<tmp0) then
               ptrange(4) = tmp0
            endif

! for perd bunch
!            this%Pts1(4,i) = this%Pts1(4,i0)
!            this%Pts1(5,i) = -gam*betai*this%Pts1(5,i0)*xl
!            this%Pts1(6,i) = this%Pts1(6,i0)

            radtest = sqrt(this%Pts1(1,i)*this%Pts1(1,i)+ &
                      this%Pts1(3,i)*this%Pts1(3,i))
            if(radtest.ge.rad) then
              ilost = ilost + 1
            endif
            this%Pts1(7,i) = this%Pts1(7,i0)
            this%Pts1(8,i) = this%Pts1(8,i0)
            this%Pts1(9,i) = this%Pts1(9,i0)
          enddo
          ptrange(5) = -halfperd
          ptrange(6) = halfperd
        else if(Flagbc.eq.3) then !round pipe, 1D z open
          ptrange(1) = 0.0
          ptrange(2) = xrad
          ptrange(3) = 0.0
          ptrange(4) = 4*asin(1.0)
          ilost = 0
          do i0 = 1, this%Nptlocal
            i = i0 - ilost
            ! The following steps go from z to t frame.
            ! 2) Linear algorithm to transfer from z to t frame.
!            rcpgammai = 1.0d0/(-this%Pts1(6,i0)+gam)
!            betai = sqrt(1.0d0-rcpgammai*rcpgammai*(1+this%Pts1(2,i0)**2+ &
!                        this%Pts1(4,i0)**2) )
            tmpx = this%Pts1(1,i0)*xl
            this%Pts1(1,i) = this%Pts1(1,i0)*xl

            this%Pts1(2,i) = this%Pts1(2,i0)
            tmpy = this%Pts1(3,i0)*xl
            this%Pts1(3,i) = this%Pts1(3,i0)*xl

            this%Pts1(4,i) = this%Pts1(4,i0)
            this%Pts1(5,i) = -gambet*this%Pts1(5,i0)*xl

            if(ptrange(5)>this%Pts1(5,i)) then
              ptrange(5) = this%Pts1(5,i)
            endif
            if(ptrange(6)<this%Pts1(5,i)) then
              ptrange(6) = this%Pts1(5,i)
            endif

            this%Pts1(6,i) = this%Pts1(6,i0)
            radtest = sqrt(this%Pts1(1,i)*this%Pts1(1,i)+ &
                      this%Pts1(3,i)*this%Pts1(3,i))
            tmp5 = sqrt(tmpx*tmpx+tmpy*tmpy)
            if((radtest.ge.rad).or.(tmp5.ge.rad)) then
              ilost = ilost + 1
            endif
            this%Pts1(7,i) = this%Pts1(7,i0)
            this%Pts1(8,i) = this%Pts1(8,i0)
            this%Pts1(9,i) = this%Pts1(9,i0)
          enddo
        else if(Flagbc.eq.4) then
          ptrange(1) = 0.0
          ptrange(2) = xrad
          ptrange(3) = 0.0
          ptrange(4) = 4*asin(1.0)
          ptrange(5) = -halfperd
          ptrange(6) = halfperd
          ilost = 0
          do i0 = 1, this%Nptlocal
            i = i0 - ilost
            ! The following steps go from z to t frame.
            ! 2) Linear algorithm to transfer from z to t frame.
!            rcpgammai = 1.0d0/(-this%Pts1(6,i0)+gam)
!            betai = sqrt(1.0d0-rcpgammai*rcpgammai*(1+this%Pts1(2,i0)**2+ &
!                        this%Pts1(4,i0)**2) )

            this%Pts1(4,i) = this%Pts1(4,i0)
            this%Pts1(5,i) = -gambet*this%Pts1(5,i0)*xl
            ntmp5 = this%Pts1(5,i)/halfperd
            !if(ntmp5.eq.1) then
            !  print*,"ntmp5: ",ntmp5,this%Pts1(5,i),halfperd,xl,gam*betai
            !endif
            tmp5 = this%Pts1(5,i) - ntmp5*halfperd
            this%Pts1(5,i) = tmp5 - mod(ntmp5,2)*halfperd
            this%Pts1(6,i) = this%Pts1(6,i0)

            tmp5 = -this%Pts1(5,i)/(gam*bet*xl)

            tmpx = this%Pts1(1,i0)*xl
            this%Pts1(1,i) = this%Pts1(1,i0)*xl
            this%Pts1(2,i) = this%Pts1(2,i0)
            tmpy = this%Pts1(3,i0)*xl
            this%Pts1(3,i) = this%Pts1(3,i0)*xl

            radtest = sqrt(this%Pts1(1,i)*this%Pts1(1,i)+ &
                      this%Pts1(3,i)*this%Pts1(3,i))
            tmp5 = sqrt(tmpx*tmpx+tmpy*tmpy)
            if((radtest.ge.rad).or.(tmp5.ge.rad)) then
              ilost = ilost + 1
            endif
            this%Pts1(7,i) = this%Pts1(7,i0)
            this%Pts1(8,i) = this%Pts1(8,i0)
            this%Pts1(9,i) = this%Pts1(9,i0)
          enddo
        else if(Flagbc.eq.5) then !2D rectangular finite, 1D z open
          ptrange(1) = -xrad
          ptrange(2) = xrad
          ptrange(3) = -yrad
          ptrange(4) = yrad
          ilost = 0
          do i0 = 1, this%Nptlocal
            i = i0 - ilost
            ! The following steps go from z to t frame.
            ! 2) Linear algorithm to transfer from z to t frame.
!            rcpgammai = 1.0d0/(-this%Pts1(6,i0)+gam)
!            betai = sqrt(1.0d0-rcpgammai*rcpgammai*(1+this%Pts1(2,i0)**2+ &
!                        this%Pts1(4,i0)**2) )
            tmpx = this%Pts1(1,i0)*xl
            this%Pts1(1,i) = this%Pts1(1,i0)*xl
            this%Pts1(2,i) = this%Pts1(2,i0)
            tmpy = this%Pts1(3,i0)*xl
            this%Pts1(3,i) = this%Pts1(3,i0)*xl
            this%Pts1(4,i) = this%Pts1(4,i0)
            this%Pts1(5,i) = -gambet*this%Pts1(5,i0)*xl

            if(ptrange(5)>this%Pts1(5,i)) then
              ptrange(5) = this%Pts1(5,i)
            endif
            if(ptrange(6)<this%Pts1(5,i)) then
              ptrange(6) = this%Pts1(5,i)
            endif

            this%Pts1(6,i) = this%Pts1(6,i0)

            if(tmpx.le.(-xrad)) then
              ilost = ilost + 1
            else if(tmpx.ge.xrad) then
              ilost = ilost + 1
            else if(tmpy.le.(-yrad)) then
              ilost = ilost + 1
            else if(tmpy.gt.yrad) then
              ilost = ilost + 1
            else if(this%Pts1(1,i).le.(-xrad)) then
              ilost = ilost + 1
            else if(this%Pts1(1,i).ge.xrad) then
              ilost = ilost + 1
            else if(this%Pts1(3,i).le.(-yrad)) then
              ilost = ilost + 1
            else if(this%Pts1(3,i).ge.yrad) then
              ilost = ilost + 1
            else
            endif
            this%Pts1(7,i) = this%Pts1(7,i0)
            this%Pts1(8,i) = this%Pts1(8,i0)
            this%Pts1(9,i) = this%Pts1(9,i0)
          enddo
        else if(Flagbc.eq.6) then
          ptrange(1) = -xrad
          ptrange(2) = xrad
          ptrange(3) = -yrad
          ptrange(4) = yrad
          ptrange(5) = -halfperd
          ptrange(6) = halfperd
          ilost = 0
          do i0 = 1, this%Nptlocal
            i = i0 - ilost
            ! The following steps go from z to t frame.
            ! 2) Linear algorithm to transfer from z to t frame.
!            rcpgammai = 1.0d0/(-this%Pts1(6,i0)+gam)
!            betai = sqrt(1.0d0-rcpgammai*rcpgammai*(1+this%Pts1(2,i0)**2+ &
!                        this%Pts1(4,i0)**2) )

            this%Pts1(4,i) = this%Pts1(4,i0)
            this%Pts1(5,i) = -gambet*this%Pts1(5,i0)*xl
            ntmp5 = this%Pts1(5,i)/halfperd
            tmp5 = this%Pts1(5,i) - ntmp5*halfperd 
            this%Pts1(5,i) = tmp5 - mod(ntmp5,2)*halfperd
            this%Pts1(6,i) = this%Pts1(6,i0)

            tmp5 = -this%Pts1(5,i)/(gambet*xl)

            tmpx = this%Pts1(1,i0)*xl
            this%Pts1(1,i) = this%Pts1(1,i0)*xl
            this%Pts1(2,i) = this%Pts1(2,i0)
            tmpy = this%Pts1(3,i0)*xl
            this%Pts1(3,i) = this%Pts1(3,i0)*xl

            if(tmpx.le.(-xrad)) then
              ilost = ilost + 1
            else if(tmpx.ge.xrad) then
              ilost = ilost + 1
            else if(tmpy.le.(-yrad)) then
              ilost = ilost + 1
            else if(tmpy.gt.yrad) then
              ilost = ilost + 1
            else if(this%Pts1(1,i).le.(-xrad)) then
              ilost = ilost + 1
            else if(this%Pts1(1,i).ge.xrad) then
              ilost = ilost + 1
            else if(this%Pts1(3,i).le.(-yrad)) then
              ilost = ilost + 1
            else if(this%Pts1(3,i).ge.yrad) then
              ilost = ilost + 1
            else
            endif
            this%Pts1(7,i) = this%Pts1(7,i0)
            this%Pts1(8,i) = this%Pts1(8,i0)
            this%Pts1(9,i) = this%Pts1(9,i0)
          enddo
        else
          print*,"no such boundary condition!!!",Flagbc
          stop
        endif

!        if(ilost.gt.0) then

        this%Nptlocal = this%Nptlocal - ilost
        nplc = this%Nptlocal
!        call MPI_ALLREDUCE(nplc,nptot,1,MPI_INTEGER,&
!                           MPI_SUM,MPI_COMM_WORLD,ierr)
        fnplc = nplc*1.0d0 
        call MPI_ALLREDUCE(fnplc,fnptot,1,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        nptot = fnptot + 0.1 

        this%Current = this%Current*nptot/this%Npt
        this%Npt = nptot

!        endif

        end subroutine conv0thB_BeamBunch

        ! find the field using the spectral method (Sin+Hermit) with
        ! particle distribution.
        subroutine kicklmnsin_Field(Nptlocal,rays,Nl,Nm,Nn,xaper,yaper,Npt,gam,tau,&
                   curr,mass,chrg,zmin)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: Nptlocal,Nl,Nm,Nn
        integer, intent(in) :: Npt
        real*8, intent (inout), dimension (9, Nptlocal) :: rays
        real*8, intent (in) :: gam,curr,tau,mass,chrg
        real*8, dimension (Nl, Nptlocal) :: tmpl,tmpcl
        real*8, dimension (Nm, Nptlocal) :: tmpm,tmpcm
        real*8, dimension (Nn, Nptlocal) :: tmpn,tmpcn
        real*8, dimension (Nl,Nm,Nn) :: philmn
        real*8, intent(in) :: xaper,yaper
        real*8 :: sig
        real*8, dimension (Nl,Nm,Nn) :: xden
        integer :: i, j,ip,ierr,k
        real*8 :: t0,t1
        real*8 :: aa,bb,cc,pilc,alphal,betam,gamman,xx,yy,zz
        real*8, dimension(4) :: range
        real*8 :: xmin,ymin,zmin,hx,hy,hz
        real*8  :: epx,epy,epz
        real*8, dimension(2) :: tmpsz,tmpszgl
        real*8 :: twopi,tmpscale,fpei,xk,bet,bbyk,gambet,brho,vz0,perv0
        real*8 :: xycon,tcon,qpert
        integer :: n
        real*8 :: zzmin
! for test
        real*8, dimension (Nn) :: tmph
        real*8 :: egxtest
        integer :: Ntest,ii,innx,inny


!find sigz
!        tmpsz = 0.0d0
!        do i = 1, Nptlocal
!          tmpsz(1) = tmpsz(1) + rays(5,i)
!          tmpsz(2) = tmpsz(2) + rays(5,i)*rays(5,i)
!        enddo
!        call MPI_ALLREDUCE(tmpsz,tmpszgl,2,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        tmpsz = tmpszgl/Npt 
!        sig = sqrt(tmpsz(2) - tmpsz(1)**2)
!        print*,"sigz: ",sig

!        zmin = -8*sig
        zzmin = 1.5*zmin
        cc = 2*abs(zzmin)

! for test
        Ntest = 129
        hx = 2*xaper/(Ntest-1)
        hy = 2*yaper/(Ntest-1)
        xmin = -xaper
        ymin = -yaper
        innx = Nl
        inny = Nm
        !hz = 2*abs(zmin)/(Ntest-1)
        hz = cc/(Ntest-1)
!

        tmpscale = curr/Scfreq
        twopi = 2.0d0*Pi
        fpei = Clight*Clight*1.0d-7
        xk = 1/Scxl

        bet = sqrt(gam**2-1.0d0)/gam
        bbyk = bet/xk
        gambet = gam*bet
        brho = gambet/(Clight*chrg)*mass
        vz0 = bet*Clight
        perv0 = 2.0d0*curr*fpei/(brho*vz0*vz0*gam*gam)
! the curr*charge/freq has been included in the charge density
        xycon = 0.5d0*perv0*gambet*bbyk*twopi/tmpscale
        tcon = bet*xycon*gam**2

!        print*,"xycon,tcon:",xycon,tcon,gam,tau,xk,gam,gambet


        pilc = 2*asin(1.0d0)
        range(1) = -xaper
        range(2) = xaper
        range(3) = -yaper
        range(4) = yaper

        aa = range(2)-range(1)
        bb = range(4)-range(3)

        call starttime_Timer( t0 )
        call starttime_Timer( t1 )

        do ip = 1, Nptlocal
          xx = rays(1,ip)-range(1)
          do i = 1, Nl
            alphal = pilc*i/aa
            tmpl(i,ip) = sin(alphal*xx)
            tmpcl(i,ip) = cos(alphal*xx)
          enddo
        enddo
        !print *, tmpcl(13,366)
        do ip = 1, Nptlocal
          yy = rays(3,ip)-range(3)
          do j =1, Nm
            betam = pilc*j/bb
            tmpm(j,ip) = sin(betam*yy)
            tmpcm(j,ip) = cos(betam*yy)
          enddo
        enddo

!sine transform in z
        do ip = 1, Nptlocal
          zz = rays(5,ip)-zzmin
          do k =1, Nn
            gamman = pilc*k/cc
            tmpn(k,ip) = sin(gamman*zz)
            tmpcn(k,ip) = cos(gamman*zz)
          enddo
        enddo
        t_ntrslo = t_ntrslo + elapsedtime_Timer( t1 )
        call starttime_Timer( t1 )

        do k = 1, Nn
          gamman = pilc*k/cc
          do j =1, Nm
            betam = pilc*j/bb
            do i = 1, Nl
              alphal = pilc*i/aa

              philmn(i,j,k) = 0.0d0
              do ip = 1, Nptlocal
                philmn(i,j,k)=philmn(i,j,k)+tmpl(i,ip)*tmpm(j,ip)*tmpn(k,ip)
              enddo
           
              philmn(i,j,k) = philmn(i,j,k)/(alphal**2+betam**2+gamman**2)
            enddo
          enddo
        enddo
        !xden = philmn

        !call MPI_ALLREDUCE(philmn,xden,Nl*Nm*Nn,MPI_DOUBLE_PRECISION,MPI_SUM,&
        !               MPI_COMM_WORLD,ierr)
        philmn = philmn*8/aa/bb/cc/Npt

        qpert = curr/Scfreq*chrg/abs(chrg)

        philmn = philmn*4*pilc*qpert
        !print*,philmn(9,9,9)
        xden = 0.0d0
        call MPI_ALLREDUCE(philmn,xden,Nl*Nm*Nn,MPI_DOUBLE_PRECISION,MPI_SUM,&
                       MPI_COMM_WORLD,ierr)
        philmn = xden


        !print *, philmn(15,11,10)
        !print*,"qpert: ",qpert

        t_fft2dhpf = t_fft2dhpf + elapsedtime_Timer( t1 )
        call starttime_Timer( t1 )

!-----------------------------------------------------------
! for test
!        yy = 0.0d0 - ymin
!        !zz = 0.0d0
!        !zz = -3.126413042203735E-007 - zzmin
!        zz = 7.242835023471561E-007 - zzmin
!
!        print*,xmin,ymin,innx,inny
!        do ii = 1, Ntest
!          egxtest = 0.0d0
!          xx = xmin + (ii-1)*hx - xmin
!
!          do k = 1, Nn
!            gamman = pilc*k/cc
!            do j = 1, inny
!              betam = pilc*j/bb
!              do i = 1, innx
!                alphal = pilc*i/aa
!
!                egxtest = egxtest - philmn(i,j,k)*alphal*cos(alphal*xx)*&
!                          sin(betam*yy)*sin(gamman*zz)
!              enddo
!            enddo
!          enddo
!          write(101,11)xx+xmin,egxtest
!        enddo
!        call flush(101)
!
!        yy = 0.0d0 - ymin
!        xx = 0.0d0 - xmin
!        print*,"zmin2:",zzmin,hz,Ntest
!
!        do ii = 1, Ntest
!          egxtest = 0.0d0
!          zz = zzmin + (ii-1)*hz - zzmin
!
!          do k = 1, Nn
!            gamman = pilc*k/cc
!            do j = 1, inny
!              betam = pilc*j/bb
!              do i = 1, innx
!                alphal = pilc*i/aa
!                egxtest = egxtest - philmn(i,j,k)*sin(alphal*xx)*sin(betam*yy)*&
!                                    cos(gamman*zz)*gamman
!              enddo
!            enddo
!          enddo
!
!          write(102,11)zz+zzmin,egxtest
!        enddo
!        call flush(102)
!
!11      format(2(1x,e17.8))
!        stop
!-----------------------------------------------------------

        do ip = 1, Nptlocal
          epx = 0.0d0
          epy = 0.0d0
          epz = 0.0d0
          do k = 1, Nn
            gamman = pilc*k/cc
            do j = 1, Nm
              betam = pilc*j/bb
              do i = 1, Nl
                alphal = pilc*i/aa

                epx = epx - philmn(i,j,k)*tmpcl(i,ip)*tmpm(j,ip)*&
                                    tmpn(k,ip)*alphal
                epy = epy - philmn(i,j,k)*tmpl(i,ip)*tmpcm(j,ip)*&
                                    tmpn(k,ip)*betam
                epz = epz - philmn(i,j,k)*tmpl(i,ip)*tmpm(j,ip)*&
                                    tmpcn(k,ip)*gamman
              enddo
            enddo
          enddo

!          print*,"epx,epy,epz,",ip,epx,epy,epz
!          stop

          !0th order algorithm to transfer back from t beam frame to z.
          rays(1,ip) = rays(1,ip)*xk
          rays(2,ip) = rays(2,ip)+tau*xycon*epx*gam
          rays(3,ip) = rays(3,ip)*xk
          rays(4,ip) = rays(4,ip)+tau*xycon*epy*gam
          rays(5,ip) = rays(5,ip)*xk/(-gambet)
          rays(6,ip) = rays(6,ip)-tau*tcon*epz
        enddo

        t_mfft_local1 = t_mfft_local1 + elapsedtime_Timer( t1 )
        t_rhofas = t_rhofas + elapsedtime_Timer( t0 )

        end subroutine kicklmnsin_Field

        subroutine setpts_BeamBunch(this,inparticles,nptsin)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision, intent(in), dimension(:,:) :: inparticles
        integer, intent(in) :: nptsin

        deallocate(this%Pts1)
        allocate(this%Pts1(9,nptsin))
        this%Pts1(:,1:nptsin) = inparticles(:,1:nptsin)

        end subroutine setpts_BeamBunch

        ! set local # of particles.
        subroutine setnpt_BeamBunch(this,innpt)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innpt
        type (BeamBunch), intent(inout) :: this

        this%Nptlocal = innpt

        end subroutine setnpt_BeamBunch
        ! get local # of particles.
        subroutine getnpt_BeamBunch(this,outnpt)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(in) :: this
        integer, intent(out) :: outnpt

        outnpt = this%Nptlocal

        end subroutine getnpt_BeamBunch

        subroutine getpts_BeamBunch(this,outparticles)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(in) :: this
        double precision, dimension(:,:), intent(out) :: outparticles

        outparticles(:,1:this%Nptlocal) = this%Pts1(:,1:this%Nptlocal)

        end subroutine getpts_BeamBunch


        subroutine destruct_BeamBunch(this)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(out) :: this

        deallocate(this%Pts1)

        end subroutine destruct_BeamBunch
      end module BeamBunchclass
