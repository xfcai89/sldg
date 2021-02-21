    subroutine green2( upstream_gl ,upstream_gauss,sum16,nsub, dxy,ixy )
    implicit none
    real,intent(in) :: upstream_gl(1:n_gl)
    real,intent(in) :: upstream_gauss(1:n_g)
    
    real,intent(in) :: dxy
    integer,intent(in) :: ixy,nsub
    real,intent(out) :: sum16(2)
    real :: sum,sum1
    integer :: ns
    real :: xrg

    real :: xctemp
    real :: temp_modal

    real :: a11,a22,bb(2),ainv(3,3)

    real :: x_r,x_l
    integer :: nm,ii,idx

    real :: xi1,xi2,xi3,xi_r,xi_l

    ! Psi^*(x) is reconstructed by (xgl(1),-0.5), (xgl(2),0.5)
    ! compute reconstruction polynomial \Psi^*(x) = a + b * x,
    !                          a + b * x1 = -0.5 ,
    !                          a + b * x2 = 0.5 ,
    !
    !     [  1   x1  ]  [ a ]       [ -0.5 ]
    !     [  1   x2  ]  [ b ]  =    [  0.5  ] , i.e. Ax=b
    !
    !
    !     inverse matrix of A =
    !    {   -x2/(x1-x2)         x1/(x1-x2)      }
    !    {     1/(x1-x2)         -1/(x1-x2)      }
    !
    !*********************************************************************************
    ! We want to avoid the case of divide a small number and we use the formula below
    !                / x_{i+1/2}
    !       1/dx int |                 \hat{u} \Psi(x) * \Psi^*(x) dx
    !                \ x_{i-1/2}
    ! =
    !       int_{-1/2}^{1/2} \hat{u}(xi) * Psi(xi) * Psi^*(xi) dxi
    ! Psi^*(xi) is reconstructed by ( (xgl(1)-xc)/dx ,-0.5), (   (xgl(2)-xc)/dx , 0.5)
    !              x1 = (xgl(1)-xc)/dx, x2 = (xgl(2)-xc)/dx
    ! compute reconstruction polynomial \Psi^*(xi) = a + b * xi,
    !                          a + b * x1 = -0.5 ,
    !                          a + b * x2 = 0.5 ,
    !
    !     [  1   x1  ]  [ a ]       [ -0.5 ]
    !     [  1   x2  ]  [ b ]  =    [  0.5  ] , i.e. Ax=b
    !
    !
    !     inverse matrix of A =
    !    {   -x2/(x1-x2)         x1/(x1-x2)      }
    !    {     1/(x1-x2)         -1/(x1-x2)      }

    do nm = 1 , n_moment
        sum = 0.
        do ns = 1 , nsub
            idx = Dn_star_sub(ns)%id

            if( ixy == 1 )then
                xctemp = x( idx )
            elseif( ixy == 2 )then
                xctemp = y( idx )
            endif

            !xi1 = ( upstream_gl(1) - xctemp )/dxy  ; xi2 =  ( upstream_gl(n_gl) - xctemp)/dxy;
            !do ii = 1 , n_moment
            !    bb(ii) = fle(nm-1, gauss(ii,1) )
            !enddo
            
            !.or.
            
            xi1 = ( upstream_gauss(1) - xctemp )/dxy  ; xi2 =  ( upstream_gauss(n_g) - xctemp )/dxy;
            do ii = 1 , n_moment
                bb(ii) = fle(nm-1, x_g(ii) )
            enddo
            
            ainv(1,1) = -xi2/(xi1-xi2);  ainv(1,2) =  xi1/(xi1-xi2);
            ainv(2,1) =   1./(xi1-xi2);  ainv(2,2) =  -1./(xi1-xi2);

            a11 = ainv(1,1)*bb(1) + ainv(1,2)*bb(2)
            a22 = ainv(2,1)*bb(1) + ainv(2,2)*bb(2)

            xi_r = ( Dn_star_sub(ns)%vpoint(2) - xctemp )/dxy
            xi_l = ( Dn_star_sub(ns)%vpoint(1) - xctemp )/dxy

            do  ii = 1 , n_moment
                if(  ixy ==1)then
                    temp_modal = elem(idx,je)%split_modal(ii)
                elseif( ixy == 2 )then
                    temp_modal = elem(ie,idx)%split_modal(ii)
                endif

                sum = sum +  temp_modal*( f2_i(xi_r,a11,a22,ii) - f2_i(xi_l,a11,a22,ii) )

                !print *,sum ,  temp_modal,( f2_i(xi_r,a11,a22,ii) - f2_i(xi_l,a11,a22,ii) )
                !pause
            enddo
          
        enddo
        sum16(nm) = sum*ai(nm)
        
    enddo

    end subroutine Green2
    !
    real function f2_i(xi16,a16,a26,k)
    implicit none
    real,intent(in) :: xi16,a16,a26
    integer,intent(in) :: k

    if(k==1)then
        f2_i = a16*xi16 + 0.5*a26*xi16**2
    elseif(k==2)then
        f2_i = 0.5*a16*xi16**2 + 1./3.*a26*xi16**3
    endif

    end function f2_i