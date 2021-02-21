    module polynomial_mod
    implicit none

    contains

    real function fle(k,x)
    implicit none
    integer, intent(in) :: k
    real,intent(in) :: x
    ! 关于Legendre正交多项式的函数
    ! v_0=1.0, v_1=(x-x_j)/dx_j, v_2=( (x-x_j)/dx_j )^2-1.0/12.0,
    ! v_3=( (x-x_j)/dx_j )^3- 0.15*(x-x_j)/dx_j
    ! v_4=( ((x-x_j)/dx_j)^2-3.0/14.0 )*((x-x_j)/dx_j)^2 + 3.0/560.0

    if(k.eq.0) then
        fle=1.0
    elseif (k.eq.1) then
        fle= x
    elseif(k.eq.2) then
        fle= x*x-1./12.0
    elseif (k.eq.3) then
        fle=x*x*x-0.15*x
    elseif(k.eq.4) then
        fle=(x**2-3./14.0)*x*x+3.0/560.0
    endif

    return
    end function fle
    !***************************************************************************
    real function polynomial(a00,x16,xc16,dx16,k)
    implicit none
    real,intent(in) :: x16,xc16,dx16
    integer,intent(in) :: k
    real,intent(in) :: a00(k+1)


    if(k==0)then
        polynomial = a00(1)
    elseif(k==1)then
        polynomial = a00(1) + a00(2)*(x16 - xc16 )/dx16
    elseif(k==2)then
        polynomial = a00(1) + a00(2)*(x16 - xc16 )/dx16 +  a00(3)*(  ((x16 - xc16) /dx16)**2 - 1./12. )
    elseif(k==3)then
        polynomial = a00(1) + a00(2)*(x16 - xc16 )/dx16 +  a00(3)*(  ((x16 - xc16) /dx16)**2 - 1./12. ) &
            +a00(4) *( ((x16 - xc16) /dx16)**3 - 3./20.*((x16 - xc16) /dx16) )
    endif

    end function polynomial


    real function polynomial2d(a00,x16,xc16,dx16,y16,yc16,dy16,k)
    implicit none
    real,intent(in) :: x16,xc16,dx16
    real,intent(in) :: y16,yc16,dy16
    integer,intent(in) :: k
    real,intent(in) :: a00(k)


    if(k==1)then
        polynomial2d = a00(1)
    elseif(k==3)then
        polynomial2d = a00(1) + a00(2)*(x16 - xc16 )/dx16 +  a00(3)*(y16 - yc16 )/dy16
    elseif(k==6)then
        polynomial2d = a00(1) + a00(2)*(x16 - xc16 )/dx16 +  a00(3)*(y16 - yc16 )/dy16 &
            + a00(4)*(  ((x16-xc16)/dx16)**2 -1./12. ) + a00(5)*(x16 - xc16 )/dx16*(y16 - yc16 )/dy16 &
            + a00(6)*(  ((y16-yc16)/dy16)**2 -1./12. )
    elseif(k==10)then
        polynomial2d = a00(1) + a00(2)*(x16 - xc16 )/dx16 +  a00(3)*(y16 - yc16 )/dy16 &
            + a00(4)*(  ((x16-xc16)/dx16)**2 -1./12. ) + a00(5)*(x16 - xc16 )/dx16*(y16 - yc16 )/dy16 &
            + a00(6)*(  ((y16-yc16)/dy16)**2 -1./12. ) &
            + a00(7)*(  ((x16-xc16)/dx16)**3 - 0.15*((x16-xc16)/dx16)   ) &
            + a00(8)*(  (  ((x16-xc16)/dx16)**2 -1./12. )*((y16-yc16)/dy16)  ) &
            + a00(9)*(  (  ((y16-yc16)/dy16)**2 -1./12. )*((x16-xc16)/dx16) )  &
            + a00(10)*(  ((y16-yc16)/dy16)**3 - 0.15*((y16-yc16)/dy16)   )
    endif


    end function polynomial2d

    !********************************************************polynomials for LDG of 2D Poisson equation
    real function flex(k,x,y) !  对x的导数
    implicit none
    integer, intent(in) :: k
    real,intent(in) :: x,y

    select case(k)
    case(0,2,5,9)
        flex=0.0
    case(1)
        flex=1.0
    case(3)
        flex=2.0*x
    case(4)
        flex=y
    case(6)
        flex=3.0*x*x-0.15
    case(7)
        flex=2.0*x*y
    case(8)
        flex=y*y-1.0/12.0
        case default
        write(*,*)'Error in function fled!'
        stop
    end select

    end function flex


    real function fley(k,x,y) !  对y的导数
    implicit none
    integer, intent(in) :: k
    real,intent(in) :: x,y

    select case(k)
    case(0,1,3,6)
        fley=0.0
    case(2)
        fley=1.0
    case(4)
        fley=x
    case(5)
        fley=2.0*y
    case(7)
        fley=x*x-1.0/12.0
    case(8)
        fley=2.0*x*y
    case(9)
        fley=3.0*y*y-0.15
        case default
        write(*,*)'Error in function fled!'
        stop
    end select

    end function fley

    real function fle_comp(k,x,y)
    implicit none
    integer, intent(in) :: k
    real,intent(in) :: x,y

    ! Basis function: Legendre polynomial. 正交，但没标准化。x,y都属于[-0.5,0.5]
    select case(k)
    case(0)
        fle_comp=1.0
    case(1)
        fle_comp=x
    case(2)
        fle_comp=y
    case(3)
        fle_comp=x*x-1.0/12.0
    case(4)
        fle_comp=x*y
    case(5)
        fle_comp=y*y-1.0/12.0
    case(6)
        fle_comp=x*x*x-0.15*x
    case(7)
        fle_comp=(x*x-1.0/12.0)*y
    case(8)
        fle_comp=x*(y*y-1.0/12.0)
    case(9)
        fle_comp=y*y*y-0.15*y
    end select

    end function fle_comp

    real function fphi(k,x,y)
    implicit none
    integer,intent(in) :: k
    real,intent(in) :: x,y

    if(k .eq. 1)then
        fphi = 1
    elseif(k .eq. 2) then
        fphi = x
    elseif(k .eq. 3)then
        fphi = y
    elseif(k .eq. 4)then
        fphi = x**2 - 1./12.
    elseif(k .eq. 5)then
        fphi = x*y
    elseif(k .eq. 6)then
        fphi = y**2 - 1./12.
    elseif( k .eq. 7 )then
        fphi = x*x*x-0.15*x
    elseif( k .eq. 8 )then
        fphi = ( x*x-1./12. )*y
    elseif( k .eq. 9 )then
        fphi = ( y*y-1./12. )*x
    elseif( k .eq. 10 )then
        fphi = y*y*y - 0.15*y
    endif


    end function fphi

    end module polynomial_mod