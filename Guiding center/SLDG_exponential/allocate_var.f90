    subroutine allocate_var
    implicit none
    integer :: k

    allocate( xgrid(1-nghost:nx+1+nghost) )
    allocate( ygrid(1-nghost:ny+1+nghost) )

    allocate( x(1-nghost:nx+nghost) )
    allocate( y(1-nghost:ny+nghost) )

    allocate( vertex(  1:(nx+1) , 1:(ny+1)  ) )
    allocate( vertex_star(  1:(nx+1) , 1:(ny+1)  ) )

    allocate( nodex( 1:nx,1:ny+1 ) )
    allocate( nodey( 1:nx+1,1:ny ) )
    allocate( nodec( 1:nx,1:ny ) )
    allocate( nodex_star( 1:nx,1:ny+1 ) )
    allocate( nodey_star( 1:nx+1,1:ny ) )
    allocate( nodec_star( 1:nx,1:ny ) )

    allocate( face_lr(1:nx,1:ny+1) )
    allocate( face_bt(1:nx+1,1:ny) )

    allocate( element(1-nghost:nx+nghost,1-nghost:ny+nghost,0:iexprk) )
    allocate( element_star(1:nx,1:ny,0:iexprk) )
    do i = 1-nghost , nx+nghost
        do j = 1-nghost , ny+nghost
            do k = 0,iexprk
                allocate( element(i,j,k)%umodal(1:n_moment) )
            enddo
        enddo
    enddo

    allocate( umod_t(1:nx,1:ny,1:n_moment) )

    allocate( com_mass(1:nx,1:ny) )


    !*******************************************
    !allocate( temp11(1:nx,1:ny,1:n_moment) )
    allocate( ele_dg(-kdg*nghost  :nx*kdg     + kdg*nghost ,0:iexprk ) )
    allocate( ee(0:nx,0:iexprk) )
    allocate( ee_c(1:nx,0:iexprk) )

    !*******************************************
    ! for incompressible
    allocate( phi_x_io(1-nghost:nx+nghost,1-nghost:ny+nghost,1:n_moment_LDG,0:iexprk) )
    allocate( phi_y_io(1-nghost:nx+nghost,1-nghost:ny+nghost,1:n_moment_LDG,0:iexprk) )

    allocate( phi_x_v(1:nx+1,1:ny+1,0:iexprk) )
    allocate( phi_x_c(1:nx,1:ny,0:iexprk) )
    allocate( phi_x_dex(1:nx,1:ny+1,0:iexprk) )
    allocate( phi_x_dey(1:nx+1,1:ny,0:iexprk) )

    allocate( phi_y_v(1:nx+1,1:ny+1,0:iexprk) )
    allocate( phi_y_c(1:nx,1:ny,0:iexprk) )
    allocate( phi_y_dex(1:nx,1:ny+1,0:iexprk) )
    allocate( phi_y_dey(1:nx+1,1:ny,0:iexprk) )

    end subroutine allocate_var
    !******************************************************
    subroutine deallocate_var
    implicit none
    integer :: k

    deallocate( xgrid )
    deallocate( ygrid )

    deallocate( x )
    deallocate( y )

    deallocate( vertex )
    deallocate( vertex_star )

    deallocate( nodex )
    deallocate( nodey )
    deallocate( nodec )
    deallocate( nodex_star )
    deallocate( nodey_star )
    deallocate( nodec_star )

    deallocate( face_lr )
    deallocate( face_bt )

    do i = 1-nghost , nx+nghost
        do j = 1-nghost , ny+nghost
            do k = 0,iexprk
                deallocate( element(i,j,k)%umodal )
            enddo
        enddo
    enddo

    deallocate( element )
    deallocate( element_star )


    deallocate( umod_t )

    deallocate( com_mass )
    !*******************************************
    !deallocate( temp11 )
    deallocate( ele_dg )
    deallocate( ee )
    deallocate( ee_c )

    !*******************************************
    ! for incompressible
    deallocate( phi_x_io )
    deallocate( phi_y_io )

    deallocate( phi_x_v )
    deallocate( phi_x_c )
    deallocate( phi_x_dex )
    deallocate( phi_x_dey )

    deallocate( phi_y_v )
    deallocate( phi_y_c )
    deallocate( phi_y_dex )
    deallocate( phi_y_dey )

    end subroutine deallocate_var
    !******************************************************