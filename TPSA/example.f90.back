program main

  use PolyMapMod
  use TPSAMod
  Implicit None
!  type (polymap) :: map1
  type (dctps)   :: x,y,z,m
!  integer ::i
  
!  map1 = polymap(5,7)
!  !line = dot2+dot1
!  !write(*,*) line
!  !line = dot2-dot1
!  !write(*,*) line
!  
!  open(unit=10, file='hello.txt')
!  
!  do i=1,map1%totallength
!    write(10,*)  map1%map(:,i)
!  enddo

!  call polymap_destructor(map1)
  !write(*,*) "Polymap Test Done!"
  
  call dctps_Initialize(4,6)
  call assign(x,3.1d0,1)
  call assign(y,4.5d0,3)
  call assign(z,2.0d0,2)

  m = sinh(cosh(sin(cos(tan(pow(sqrt(log(exp(x)+y-x/y*z)),2.3d0))))))
  call m%output()
  
end program main
