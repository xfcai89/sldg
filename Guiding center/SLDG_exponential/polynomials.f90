    !real function fphi(k,x,y)
    !implicit none
    !integer,intent(in) :: k
    !real,intent(in) :: x,y
    !
    !if(k .eq. 1)then
    !    fphi = 1
    !elseif(k .eq. 2) then
    !    fphi = x
    !elseif(k .eq. 3)then
    !    fphi = y
    !elseif(k .eq. 4)then
    !    fphi = x**2 - 1./12.
    !elseif(k .eq. 5)then
    !    fphi = x*y
    !elseif(k .eq. 6)then
    !    fphi = y**2 - 1./12.
    !endif
    !
    !
    !end function fphi
    !!****************************************
    !real function polynomial(a00,x16,xc16,dx16,y16,yc16,dy16,k)
    !implicit none
    !real,intent(in) :: x16,xc16,dx16
    !real,intent(in) :: y16,yc16,dy16
    !integer,intent(in) :: k
    !real,intent(in) :: a00(k)
    !
    !
    !if(k==1)then
    !    polynomial = a00(1)
    !elseif(k==3)then
    !    polynomial = a00(1) + a00(2)*(x16 - xc16 )/dx16 +  a00(3)*(y16 - yc16 )/dy16
    !elseif(k==6)then
    !    polynomial = a00(1) + a00(2)*(x16 - xc16 )/dx16 +  a00(3)*(y16 - yc16 )/dy16 &
    !        + a00(4)*(  ((x16-xc16)/dx16)**2 -1./12. ) + a00(5)*(x16 - xc16 )/dx16*(y16 - yc16 )/dy16 &
    !        + a00(6)*(  ((y16-yc16)/dy16)**2 -1./12. )
    !endif
    !
    !
    !end function polynomial
    !!*************************************************
    !real function polynomial_1d(a00,x16,xc16,dx16,k)
    !implicit none
    !real,intent(in) :: x16,xc16,dx16
    !integer,intent(in) :: k
    !real,intent(in) :: a00(k+1)
    !
    !
    !if(k==0)then
    !    polynomial_1d = a00(1)
    !elseif(k==1)then
    !    polynomial_1d = a00(1) + a00(2)*(x16 - xc16 )/dx16
    !elseif(k==2)then
    !    polynomial_1d = a00(1) + a00(2)*(x16 - xc16 )/dx16 +  a00(3)*(  ((x16 - xc16) /dx16)**2 - 1./12. )
    !elseif(k==3)then
    !    polynomial_1d = a00(1) + a00(2)*(x16 - xc16 )/dx16 +  a00(3)*(  ((x16 - xc16) /dx16)**2 - 1./12. ) &
    !    +a00(4) *( ((x16 - xc16) /dx16)**3 - 3./20.*((x16 - xc16) /dx16) )
    !endif
    !
    !end function polynomial_1d
    !