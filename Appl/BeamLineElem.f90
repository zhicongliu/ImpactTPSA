!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! BeamLineElemclass: Beam line element base class in Lattice module of 
!                    APPLICATION layer.
! Version: 1.0
! Author: Ji Qiang, LBNL
! Description: This class defines the base beam line element class for
!              different lattice element class.
! Comments:
!----------------------------------------------------------------
      module BeamLineElemclass
        use DriftTubeclass
        use ConstFocclass
        type BeamLineElem
!          private
          type (DriftTube), pointer :: pdrift
          type (ConstFoc), pointer :: pcf
        end type BeamLineElem
        interface assign_BeamLineElem
          module procedure assign_drift,assign_cf
        end interface
        interface getparam_BeamLineElem
          module procedure getparam1_BeamLineElem, &
                           getparam2_BeamLineElem, &
                           getparam3_BeamLineElem
        end interface
        interface setparam_BeamLineElem
          module procedure setparam1_BeamLineElem, &
                           setparam2_BeamLineElem, &
                           setparam3_BeamLineElem
        end interface
      contains

        subroutine getparam1_BeamLineElem(this,i,blparam)
        implicit none 
        type (BeamLineElem), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        if(associated(this%pdrift)) then
          call getparam_DriftTube(this%pdrift,i,blparam)
        elseif(associated(this%pcf)) then
          call getparam_ConstFoc(this%pcf,i,blparam)
        endif

        end subroutine getparam1_BeamLineElem
  
        subroutine getparam2_BeamLineElem(this,blparams)
        implicit none
        type (BeamLineElem), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        if(associated(this%pdrift)) then
          call getparam_DriftTube(this%pdrift,blparams)
        elseif(associated(this%pcf)) then
          call getparam_ConstFoc(this%pcf,blparams)
        endif

        end subroutine getparam2_BeamLineElem

        subroutine getparam3_BeamLineElem(this,blength,bnseg,bmapstp,&
                                          btype)
        implicit none
        type (BeamLineElem), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        if(associated(this%pdrift)) then
          call getparam_DriftTube(this%pdrift,blength,bnseg,bmapstp,&
                                  btype)
        elseif(associated(this%pcf)) then
          call getparam_ConstFoc(this%pcf,blength,bnseg,bmapstp,btype)
        endif

        end subroutine getparam3_BeamLineElem
       
        subroutine getradius_BeamLineElem(this,piperadius,piperadius2)
        implicit none 
        type (BeamLineElem), intent(in) :: this
        double precision, intent(out) :: piperadius,piperadius2

        if(associated(this%pdrift)) then
          call getparam_DriftTube(this%pdrift,2,piperadius)
          piperadius2 = piperadius
        elseif(associated(this%pcf)) then
          call getparam_ConstFoc(this%pcf,5,piperadius)
          piperadius2 = piperadius
        endif

        end subroutine getradius_BeamLineElem
  
        subroutine setparam1_BeamLineElem(this,i,blparam)
        implicit none 
        type (BeamLineElem), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: blparam

        if(associated(this%pdrift)) then
          call setparam_DriftTube(this%pdrift,i,blparam)
        elseif(associated(this%pcf)) then
          call setparam_ConstFoc(this%pcf,i,blparam)
        endif

        end subroutine setparam1_BeamLineElem
  
        subroutine setparam2_BeamLineElem(this,blparams)
        implicit none
        type (BeamLineElem), intent(inout) :: this
        double precision, dimension(:), intent(in) :: blparams

        if(associated(this%pdrift)) then
          call setparam_DriftTube(this%pdrift,blparams)
        elseif(associated(this%pcf)) then
          call setparam_ConstFoc(this%pcf,blparams)
        endif

        end subroutine setparam2_BeamLineElem

        subroutine setparam3_BeamLineElem(this,blength,bnseg,bmapstp,&
                                          btype)
        implicit none
        type (BeamLineElem), intent(inout) :: this
        double precision, intent(in) :: blength
        integer, intent(in) :: bnseg,bmapstp,btype

        if(associated(this%pdrift)) then
          call setparam_DriftTube(this%pdrift,bnseg,bmapstp,&
                                  btype,blength)
        elseif(associated(this%pcf)) then
          call setparam_ConstFoc(this%pcf,bnseg,bmapstp,btype,blength)
        endif

        end subroutine setparam3_BeamLineElem

        function assign_drift(tdrift) result(ppdrift)
        type (BeamLineElem) :: ppdrift
        type (DriftTube), target, intent(in) :: tdrift

        ppdrift%pdrift => tdrift
        nullify(ppdrift%pcf)

        end function assign_drift

        function assign_cf(tcf) result(ppcf)
        type (BeamLineElem) :: ppcf
        type (ConstFoc), target, intent(in) :: tcf

        ppcf%pcf => tcf
        nullify(ppcf%pdrift)

        end function assign_cf

        subroutine maplinear_BeamLineElem(this,t,tau,xm,refpt,chg,mss)
        implicit none
        type (BeamLineElem), intent(inout) :: this
        double precision, intent(in) :: t,tau,chg,mss
        double precision, dimension(6,6), intent(out) ::xm
        double precision, dimension(6), intent(inout) :: refpt
        integer :: i

        if(associated(this%pdrift)) then
          call maplinear_DriftTube(t,tau,xm,this%pdrift,refpt,chg,mss)
        elseif(associated(this%pcf)) then
          call maplinear_ConstFoc(t,tau,xm,this%pcf,refpt,chg,mss)
        endif
        !do i=1,6
        !print*,xm(i,:)
        !enddo


        end subroutine maplinear_BeamLineElem

      end module BeamLineElemclass
