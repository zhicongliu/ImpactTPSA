Module PolyMapMod
  use MathFuncMod
  
  implicit none
  type  polymap
    integer :: dim
    integer, private :: max_order
    integer, allocatable, dimension(:,:) :: map
    integer(kind=4) :: totallength
  contains
    PROCEDURE, PASS :: set => setindexmap
    PROCEDURE, PASS :: dec => decomposite
  end type polymap
  
  interface polymap
    module procedure constructor
  end interface
  
  
contains

    function constructor(dim,order) result(this)
        implicit none
        type(polymap) :: this
        integer, intent(in) :: dim,order
        !write(*,*) 'Constructor works'
        this%dim   = dim
        this%max_order = order
        
        call setindexMap(this)
    end function
    
    subroutine polymap_destructor(this)
        type(polymap) :: this
        !write(*,*) 'Destructor works'
        if (ALLOCATED(this % map)) deallocate(this % map)
    end subroutine
    
    ! implement instance constructor with arguments
    subroutine setindexmap(this)
      class (polymap) :: this
      integer(kind=4) :: totallength
      integer :: i
      totallength = binomial(this%max_order+this%dim, this%dim)
      this%totallength = totallength
      if (.not. ALLOCATED(this % map)) allocate(this % map(this%dim+1,totallength))
      this %map = 0
      
      do i=1,totallength
        call decomposite(this,i-1,this % map(:,i))
      enddo
      
    end subroutine setindexmap
    
    subroutine decomposite(this,n,results)
      class (polymap) :: this
      integer n,itemp,i,k
      integer :: results(this%dim+1)
      
      itemp = n+1
      do i = this%dim,1,-1
        k=i-1
        do while (binomial(k, i)<itemp)
          k=k+1
        enddo
        itemp=itemp-binomial(k-1, i)
        results(this%dim-i+1)=k-i
      enddo
      do i = this%dim,1,-1
        results(i+1)=results(i)-results(i+1)
      enddo
      !if (n==0) write(*,*) results
    end subroutine decomposite

End Module PolyMapMod
