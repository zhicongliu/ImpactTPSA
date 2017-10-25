program main

  use PolyMapMod
  use TPSAMod
  Implicit None
  double precision, parameter :: CLIGHT = 3.0d8
  
  type (dctps)        :: x,xp,y,yp
  type (dctps)        :: bx,by,Fx,Fy
  double precision    :: z,vz
  
  double precision    :: gradient,length,Bro,k
  double precision    :: charge, mass, gam, beta
  
  double precision    :: R(2,2)
  integer       :: stepNumber,i
  
  call dctps_Initialize(4,6)
  
  call assign(x,  0.0d0, 1)
  call assign(xp, 0.0d0, 2)
  call assign(y,  0.0d0, 3)
  call assign(yp, 0.0d0, 4)
  
  z           = 0.0
  vz          = 2.d8        !m/s
  gradient    = 2.d0        !T/m
  length      = 0.1         !m
  charge      = 1.602e-19   !C,  electron
  mass        = 9.109e-31   !kg, electron
  
  stepNumber  = INT(1e4)
  
  

  beta  = vz/CLIGHT
  gam   = 1/sqrt(1-beta**2)
  
  print*, gam
  
  do i=1,stepNumber
    bx = gradient * y
    by = gradient * x
    
    Fx = -charge / (gam*mass) * by
    Fy =  charge / (gam*mass) * bx
    
    xp = xp + Fx  * (length / DBLE(stepNumber)  /2 /vz)
    yp = yp + Fy  * (length / DBLE(stepNumber)  /2 /vz)
    
    x  = x  + xp*vz * (length / DBLE(stepNumber) / vz)
    y  = y  + yp*vz * (length / DBLE(stepNumber) / vz)
    
    xp = xp + Fx  * (length / DBLE(stepNumber)  /2 /vz)
    yp = yp + Fy  * (length / DBLE(stepNumber)  /2 /vz)
    
    z  = z  +       length / DBLE(stepNumber)
  enddo
  
  Bro=mass*CLIGHT/charge*gam*beta      !T m
  k=sqrt(abs(gradient/Bro))
  

  R(1,1)=cos(k*length);
  R(2,1)=sin(k*length)/k;
  R(1,2)=-k*sin(k*length);
  R(2,2)=cos(k*length);
  
  !x  = R(1,1) * x + R(2,1) * xp
  !xp = R(1,2) * x + R(2,2) * xp
  call  x%output()
  call xp%output()
  
  print*, R(1,1),   x%map(2)!,   (x%map(2)-R(1,1))/R(1,1)
  print*, R(2,1),   x%map(3)!,   (x%map(3)-R(2,1))/R(2,1)
  print*, R(1,2),  xp%map(2)!,  (xp%map(2)-R(1,2))/R(1,2)
  print*, R(2,2),  xp%map(3)!,  (xp%map(3)-R(2,2))/R(2,2)
  
  print*,    x%map(2)!,   (x%map(2)-R(1,1))/R(1,1)
  print*,    x%map(3)!,   (x%map(3)-R(2,1))/R(2,1)
  print*,   xp%map(2)!,  (xp%map(2)-R(1,2))/R(1,2)
  print*,   xp%map(3)!,  (xp%map(3)-R(2,2))/R(2,2)
  
end program main
