Module MathFuncMod

contains
recursive Function factorial(n) result(f)
    implicit none
    integer :: n
    integer  (kind = 4) :: f
    if (n<0) then
      f=0
    else
      select case (n) 
          case (0)
              f =  1
          case (1)
              f =  1
          case (2)
              f =  2
          case ( 3)
              f =  6
          case ( 4)
              f =  24
          case ( 5)
              f =  120
          case ( 6)
              f =  720
          case ( 7)
              f =  5040
          case ( 8)
              f =  40320
          case ( 9)
              f =  362880
          case ( 10)
              f =  3628800
          case ( 11)
              f =  39916800
          case ( 12)
              f =  479001600
          case default
              f =  n*factorial(n-1)
      end select
    endif
end function factorial

recursive Function doublefactorial(n)  result(f)
    integer, intent(in):: n
    integer  (kind = 4) :: f
    if (n<0) then
      f=0
      return
    else
      select case (n)
        case (0)
            f =  1
        case (1)
            f =  1
        case (2)
            f =  2
        case (3)
            f =  3
        case (4)
            f =  8
        case (5)
            f =  15
        case (6)
            f =  48
        case (7)
            f =  105
        case (8)
            f =  384
        case (9)
            f =  945
        case (10)
            f =  3840
        case (11)
            f =  10395
        case (12)
            f =  46080
        case (13)
            f =  135135
        case (14)
            f =  645120
        case (15)
            f =  2027025
        case (16)
            f =  10321920
        case (17)
            f =  34459425
        case (18)
            f =  185794560
        case default
            f =  n*doublefactorial(n-2)
      end select
    endif
end function doublefactorial


recursive Function binomial(n, m) result(f)
    implicit none
    integer :: n,m,ml
    integer  (kind = 4) :: f
    if (n<=0 .or. m>n .or. m<0) then
        f = 0;
        return
    endif
    
    if (m>n/2) then
      ml=n-m;
    else 
      ml=m;
    endif
    
    select case (ml)
        case (0)
            f = 1
            return
        case (1)
            f = n
            return
        case (2)
            f = n*(n-1)/2
            return
        case (3)
            f = n*(n-1)*(n-2)/6
            return
    end select

    if (n<=12) then
        f = factorial(n)/factorial(m)/factorial(n-m)
    else 
        f = (n*binomial(n-1,ml-1)/ml)
    endif

end function binomial

End Module MathFuncMod
