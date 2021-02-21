    subroutine grid_maker_GL
    implicit none



    ! The physical domain
    xleft = 0.
    xright = 2.*pi

    yleft = 0.
    yright = 2.*pi

    ! width of an element
    dx = (xright - xleft)/nel_x
    dy = (yright - yleft)/nel_y

    ! x(i) ----> center
    ! y(j) ----> center
    DO I = 1 - ighost, nel_x + ighost
        X(I) = XLEFT + (I-0.5) * DX
    ENDDO

    DO J =  1 -ighost , Nel_y + ighost
        Y(J) = YLEFT + (J-0.5) * DY
    ENDDO

    DO I =   1 -ighost , nel_x + 1 + ighost 
        Xgrid(I) = XLEFT + (I-1.) * DX
    ENDDO

    DO J =   1 - ighost , Nel_y+ 1 + ighost
        ygrid(J) = YLEFT + (J-1.) * DY
    ENDDO
    
    xg(1)=-0.466234757101576013906150777246997304567d0
    xg(2)=-0.330604693233132256830699797509952673503d0
    xg(3)=-0.119309593041598454315250860840355967709d0
    xg(4)=-xg(3)
    xg(5)=-xg(2)
    xg(6)=-xg(1)
    wg(1)=1.71324492379170345040296142172733d-1/2d0
    wg(2)=3.60761573048138607569833513837716d-1/2d0
    wg(3)=4.67913934572691047389870343989551d-1/2d0
    wg(4)=wg(3)
    wg(5)=wg(2)
    wg(6)=wg(1)
    
    
    if( n_g .eq. 2 )then
        x_GL( 1 ) = - 0.5
        x_GL( 2 ) = 0.5

        w_GL( 1 ) = 0.5
        w_GL( 2 ) = 0.5
        
        
        x_g(1) = -sqrt(3.)/6.
        x_g(2) = -x_g(1)
        
        w_g(1) = 0.5
        w_g(2) = 0.5
       
        ai(1) = 1.
        ai(2) = 12.
    elseif( n_g .eq. 3 )then
        x_GL( 1 ) = -0.5
        x_GL( 2 ) = 0.5
        
        x_g(1) = -sqrt(15.)/10.
        x_g(2) = 0.
        x_g(3) = sqrt(15.)/10.
        
        w_g(1) = 5./18.
        w_g(2) = 4./9.
        w_g(3) = 5./18.

        
        ai( 1 ) = 1.
        ai( 2 ) = 12.
        ai( 3 ) = 180.
    elseif(n_g .eq. 4)then
        x_GL( 1 ) = -0.5
        x_GL( 2 ) = -0.5/sqrt(5.)
        x_GL( 3 ) = 0.5/sqrt(5.)
        x_GL( 4 ) = 0.5

        w_GL( 1 ) = 1./12.
        w_GL( 2 ) = 5./12.
        w_GL( 3 ) = 5./12.
        w_GL( 4 ) = 1./12.

        
        ai( 1 ) = 1.
        ai( 2 ) = 12.
        ai( 3 ) = 180.
        !ai( 4 ) = 2800.
    elseif( n_g .eq. 5 )then
        x_GL( 1 ) = -0.5
        x_GL( 2 ) = -0.5*sqrt( 3./7. )
        x_GL( 3 ) = 0.
        x_GL( 4 ) = -x_GL(2)
        x_GL( 5 ) = -x_GL(1)

        w_GL( 1 ) = 1./20.
        w_GL( 2 ) = 49./180.
        w_GL( 3 ) = 16./45.
        w_GL( 4 ) = w_GL(2)
        w_GL( 5 ) = w_GL(1)
     
        
        ai( 1 ) = 1.
        ai( 2 ) = 12.
        ai( 3 ) = 180.
        ai( 4 ) = 2800.
        !ai( 5 ) = 44100.0
    endif

    do i = 1 , nel_x
        do k = 1 , n_gl
            glx(k,:,i,: ) = x_GL( k )*dx + x(i)

        enddo
        do k = 1 , n_g
            gx(k,:,i,: ) = x_G( k )*dx + x(i)

        enddo
    enddo
           
    do j = 1 , nel_y
        do k = 1 , n_gl
            gly( :,k,:,j ) = x_GL( k )*dy + y(j)
        enddo
        do k = 1 , n_g
            gy( :,k,:,j ) = x_G( k )*dy + y(j)
        enddo
    enddo


    



    
    end subroutine grid_maker_GL