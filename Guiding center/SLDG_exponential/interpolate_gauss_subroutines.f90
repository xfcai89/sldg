    subroutine local_rk4_gauss
    implicit none

    type(type_point) :: cell_gau(1:3,1:3) ! only gauss points

    integer :: i,j
    integer :: ig,jg ! loop from 1:nk+1

    integer :: io1
    real :: phix(0:3),phiy(0:3),vel_x, vel_y
    integer :: idx,idy
    real :: vx1(1:2),vx2(1:2),vx3(1:2)

    !*******************************************************************
    ! programmed in a cell by cell way
    ! for getting the ``gauss points''
    ! which is from the gauss points in the Eulerian cells
    ! !
    ! we can use nodec for the center of each background cell
    ! and uniform cell length of side, dx,dy
    !******

    io1=io
    if(io==3) io1= io-1

    do i = 1,nx
        do j = 1,ny
            !*******
            do ig = 1,nk+1
                do jg = 1,nk+1
                    cell_gau(ig,jg)%coor(1)=nodec(i,j)%coor(1)+dx*gau(ig,1)
                    cell_gau(ig,jg)%coor(2)=nodec(i,j)%coor(2)+dy*gau(jg,1)

                    ! rk4 below
                    idx = ceiling( (cell_gau(ig,jg)%coor(1)-xleft)/dx )
                    idy = ceiling( (cell_gau(ig,jg)%coor(2)-ybottom)/dy )
                    call get_velocity_nodes(idx,idy,io1,cell_gau(ig,jg)%coor(1:2),phix,phiy )
                    call get_rk_stage_velocity(io1,phix(0:io1),phiy(0:io1),vel_x,vel_y)
                    !
                    vx1(1) = cell_gau(ig,jg)%coor(1) - 0.5*vel_x * dt
                    vx1(2) = cell_gau(ig,jg)%coor(2) - 0.5*vel_y * dt
                    !
                    !***********************************************
                    idx = ceiling( (vx1(1)-xleft)/dx )
                    idy = ceiling( (vx1(2)-ybottom)/dy )
                    !
                    call get_velocity_nodes(idx,idy,io1,vx1(1:2),phix,phiy )
                    call get_rk_stage_velocity(io1,phix(0:io1),phiy(0:io1),vel_x,vel_y)
                    !
                    vx2(1) = cell_gau(ig,jg)%coor(1) - 0.5*vel_x * dt
                    vx2(2) = cell_gau(ig,jg)%coor(2) - 0.5*vel_y * dt
                    !!********************************************************
                    idx = ceiling( (vx2(1)-xleft)/dx )
                    idy = ceiling( (vx2(2)-ybottom)/dy )
                    !
                    call get_velocity_nodes(idx,idy,io1,vx2(1:2),phix,phiy )
                    call get_rk_stage_velocity(io1,phix(0:io1),phiy(0:io1),vel_x,vel_y)
                    !
                    vx3(1) = cell_gau(ig,jg)%coor(1) - vel_x * dt
                    vx3(2) = cell_gau(ig,jg)%coor(2) - vel_y * dt
                    !*********
                    !*********
                    !********************************************************
                    idx = ceiling( (vx3(1)-xleft)/dx )
                    idy = ceiling( (vx3(2)-ybottom)/dy )
                    !
                    call get_velocity_nodes(idx,idy,io1,vx3(1:2),phix,phiy )
                    call get_rk_stage_velocity(io1,phix(0:io1),phiy(0:io1),vel_x,vel_y)
                    !
                    element_star(i,j,io)%gauss(ig,jg)%coor(1) =  1./3.*( -cell_gau(ig,jg)%coor(1)+vx1(1)+2.*vx2(1)+vx3(1) )  - 1./6.*vel_x*dt
                    element_star(i,j,io)%gauss(ig,jg)%coor(2) =  1./3.*( -cell_gau(ig,jg)%coor(2)+vx1(2)+2.*vx2(2)+vx3(2) )  - 1./6.*vel_y*dt
                enddo
            enddo
            !*******
        enddo
    enddo


    end subroutine local_rk4_gauss
    !*******************************************************************
    subroutine interpolate_gauss(nm,pe,pes,temp_L,temp_U,aa)
    implicit none
    integer,intent(in) :: nm
    type(type_element),pointer :: pe
    type(type_element_upstream),pointer :: pes
    real,intent(inout) :: temp_L(6,6),temp_U(6,6)
    real,intent(out) :: aa(1:n_moment)    
    
    real :: x_base,y_base
    real :: vert_temp(1:9,1:5)
    real :: A_temp(6,6)
    real :: b(1:6),psi(1:9)

    
    x_base = 0.25*( pes%vertex1%coor(1) + pes%vertex2%coor(1)  &
        + pes%vertex3%coor(1)  + pes%vertex4%coor(1)  )
    y_base = 0.25*( pes%vertex1%coor(2) + pes%vertex2%coor(2)  &
        + pes%vertex3%coor(2)  + pes%vertex4%coor(2)  )
    
    if(nm == 1)then
        aa(1) = 1.
        aa(2:n_moment) = 0.
    else
        if( nm ==2 )then

            if(n_moment<=3)then
                vert_temp(1,1) = ( pes%gauss(1,1)%coor(1) - x_base )/dx
                vert_temp(2,1) = ( pes%gauss(2,1)%coor(1) - x_base )/dx
                vert_temp(3,1) = ( pes%gauss(1,2)%coor(1) - x_base )/dx
                vert_temp(4,1) = ( pes%gauss(2,2)%coor(1) - x_base )/dx

                vert_temp(1,2) = ( pes%gauss(1,1)%coor(2) - y_base )/dy
                vert_temp(2,2) = ( pes%gauss(2,1)%coor(2) - y_base )/dy
                vert_temp(3,2) = ( pes%gauss(1,2)%coor(2) - y_base )/dy
                vert_temp(4,2) = ( pes%gauss(2,2)%coor(2) - y_base )/dy
                call get_matrix_a( vert_temp(1:4,1:2),  A_temp(1:3,1:3) )

                call doolittle(A_temp(1:3,1:3),temp_L(1:3,1:3),temp_U(1:3,1:3),3)
            elseif(n_moment==6)then
                vert_temp(1,1) = ( pes%gauss(1,1)%coor(1) - x_base )/dx
                vert_temp(2,1) = ( pes%gauss(2,1)%coor(1) - x_base )/dx
                vert_temp(3,1) = ( pes%gauss(3,1)%coor(1) - x_base )/dx
                vert_temp(4,1) = ( pes%gauss(1,2)%coor(1) - x_base )/dx
                vert_temp(5,1) = ( pes%gauss(2,2)%coor(1) - x_base )/dx
                vert_temp(6,1) = ( pes%gauss(3,2)%coor(1) - x_base )/dx
                vert_temp(7,1) = ( pes%gauss(1,3)%coor(1) - x_base )/dx
                vert_temp(8,1) = ( pes%gauss(2,3)%coor(1) - x_base )/dx
                vert_temp(9,1) = ( pes%gauss(3,3)%coor(1) - x_base )/dx

                vert_temp(1,2) = ( pes%gauss(1,1)%coor(2) - y_base )/dy
                vert_temp(2,2) = ( pes%gauss(2,1)%coor(2) - y_base )/dy
                vert_temp(3,2) = ( pes%gauss(3,1)%coor(2) - y_base )/dy
                vert_temp(4,2) = ( pes%gauss(1,2)%coor(2) - y_base )/dy
                vert_temp(5,2) = ( pes%gauss(2,2)%coor(2) - y_base )/dy
                vert_temp(6,2) = ( pes%gauss(3,2)%coor(2) - y_base )/dy
                vert_temp(7,2) = ( pes%gauss(1,3)%coor(2) - y_base )/dy
                vert_temp(8,2) = ( pes%gauss(2,3)%coor(2) - y_base )/dy
                vert_temp(9,2) = ( pes%gauss(3,3)%coor(2) - y_base )/dy

                vert_temp(1:9,3) = vert_temp(1:9,1)*vert_temp(1:9,1) - 1./12.
                vert_temp(1:9,4) = vert_temp(1:9,1)*vert_temp(1:9,2)
                vert_temp(1:9,5) = vert_temp(1:9,2)*vert_temp(1:9,2) - 1./12.

                call get_matrix2_a( vert_temp(1:9,1:5),  A_temp(1:6,1:6) )
                call doolittle(A_temp(1:6,1:6),temp_L(1:6,1:6),temp_U(1:6,1:6),6)
            endif
        endif

        if(n_moment<=3)then
            psi(1) = fphi( nm, gau(1,1), gau(1,1) ) ; ! ig=1,jg=1
            psi(2) = fphi( nm, gau(2,1), gau(1,1) ) ; ! ig=2,jg=1
            psi(3) = fphi( nm, gau(1,1), gau(2,1) ) ; ! ig=1,jg=2
            psi(4) = fphi( nm, gau(2,1), gau(2,1) ) ; ! ig=2,jg=2

            call get_vector_b( vert_temp(1:4,1:2),psi(1:4),b(1:3) )
            call solve17(temp_L(1:3,1:3),temp_U(1:3,1:3),b(1:3),aa(1:3),3)
        elseif(n_moment==6)then
            psi(1) = fphi( nm, gau(1,1), gau(1,1) ) ; ! ig=1,jg=1
            psi(2) = fphi( nm, gau(2,1), gau(1,1) ) ; ! ig=2,jg=1
            psi(3) = fphi( nm, gau(3,1), gau(1,1) ) ; ! ig=3,jg=1
            psi(4) = fphi( nm, gau(1,1), gau(2,1) ) ; ! ig=1,jg=2
            psi(5) = fphi( nm, gau(2,1), gau(2,1) ) ; ! ig=2,jg=2
            psi(6) = fphi( nm, gau(3,1), gau(2,1) ) ; ! ig=2,jg=2
            psi(7) = fphi( nm, gau(1,1), gau(3,1) ) ; ! ig=1,jg=3
            psi(8) = fphi( nm, gau(2,1), gau(3,1) ) ; ! ig=2,jg=3
            psi(9) = fphi( nm, gau(3,1), gau(3,1) ) ; ! ig=3,jg=3

            call get_vector2_b( vert_temp(1:9,1:5),psi(1:9),b(1:6) )
            call solve17(temp_L(1:6,1:6),temp_U(1:6,1:6),b(1:6),aa(1:6),6)
        endif

        !if(time>0.8)write(99,*) i,j,nm,aa(1:6)

    endif

    end subroutine interpolate_gauss
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    !*******************************************************************
    subroutine interpolate_gauss_gl(nm,pe,pes,temp_L,temp_U,aa)
    implicit none
    integer,intent(in) :: nm
    type(type_element),pointer :: pe
    type(type_element_upstream),pointer :: pes
    real,intent(inout) :: temp_L(6,6),temp_U(6,6)
    real,intent(out) :: aa(1:n_moment)    
    
    real :: x_base,y_base
    real :: vert_temp(1:9,1:5)
    real :: A_temp(6,6)
    real :: b(1:6),psi(1:9)

    
    x_base = 0.25*( pes%vertex1%coor(1) + pes%vertex2%coor(1)  &
        + pes%vertex3%coor(1)  + pes%vertex4%coor(1)  )
    y_base = 0.25*( pes%vertex1%coor(2) + pes%vertex2%coor(2)  &
        + pes%vertex3%coor(2)  + pes%vertex4%coor(2)  )
    
    if(nm == 1)then
        aa(1) = 1.
        aa(2:n_moment) = 0.
    else
        if( nm ==2 )then

            if(n_moment<=3)then
                vert_temp(1,1) = ( pes%vertex1%coor(1) - x_base )/dx
                vert_temp(2,1) = ( pes%vertex2%coor(1) - x_base )/dx
                vert_temp(3,1) = ( pes%vertex3%coor(1) - x_base )/dx
                vert_temp(4,1) = ( pes%vertex4%coor(1) - x_base )/dx

                vert_temp(1,2) = ( pes%vertex1%coor(2) - y_base )/dy
                vert_temp(2,2) = ( pes%vertex2%coor(2) - y_base )/dy
                vert_temp(3,2) = ( pes%vertex3%coor(2) - y_base )/dy
                vert_temp(4,2) = ( pes%vertex4%coor(2) - y_base )/dy
                call get_matrix_a( vert_temp(1:4,1:2),  A_temp(1:3,1:3) )

                call doolittle(A_temp(1:3,1:3),temp_L(1:3,1:3),temp_U(1:3,1:3),3)
            elseif(n_moment==6)then
                vert_temp(1,1) = ( pes%vertex1%coor(1) - x_base )/dx
                vert_temp(2,1) = ( pes%vertex2%coor(1) - x_base )/dx
                vert_temp(3,1) = ( pes%vertex3%coor(1) - x_base )/dx
                vert_temp(4,1) = ( pes%vertex4%coor(1) - x_base )/dx
                vert_temp(5,1) = ( pes%vertex5%coor(1) - x_base )/dx
                vert_temp(6,1) = ( pes%vertex6%coor(1) - x_base )/dx
                vert_temp(7,1) = ( pes%vertex7%coor(1) - x_base )/dx
                vert_temp(8,1) = ( pes%vertex8%coor(1) - x_base )/dx
                vert_temp(9,1) = ( pes%vertex9%coor(1) - x_base )/dx

                vert_temp(1,2) = ( pes%vertex1%coor(2) - y_base )/dy
                vert_temp(2,2) = ( pes%vertex2%coor(2) - y_base )/dy
                vert_temp(3,2) = ( pes%vertex3%coor(2) - y_base )/dy
                vert_temp(4,2) = ( pes%vertex4%coor(2) - y_base )/dy
                vert_temp(5,2) = ( pes%vertex5%coor(2) - y_base )/dy
                vert_temp(6,2) = ( pes%vertex6%coor(2) - y_base )/dy
                vert_temp(7,2) = ( pes%vertex7%coor(2) - y_base )/dy
                vert_temp(8,2) = ( pes%vertex8%coor(2) - y_base )/dy
                vert_temp(9,2) = ( pes%vertex9%coor(2) - y_base )/dy

                vert_temp(1:9,3) = vert_temp(1:9,1)*vert_temp(1:9,1) - 1./12.
                vert_temp(1:9,4) = vert_temp(1:9,1)*vert_temp(1:9,2)
                vert_temp(1:9,5) = vert_temp(1:9,2)*vert_temp(1:9,2) - 1./12.

                call get_matrix2_a( vert_temp(1:9,1:5),  A_temp(1:6,1:6) )
                call doolittle(A_temp(1:6,1:6),temp_L(1:6,1:6),temp_U(1:6,1:6),6)
            endif
        endif

        if(n_moment<=3)then
            psi(1) = fphi( nm, (pe%vertex1%coor(1)-x(i) )/dx,(pe%vertex1%coor(2)-y(j))/dy ) ;
            psi(2) = fphi( nm, (pe%vertex2%coor(1)-x(i) )/dx,(pe%vertex2%coor(2)-y(j))/dy ) ;
            psi(3) = fphi( nm, (pe%vertex3%coor(1)-x(i) )/dx,(pe%vertex3%coor(2)-y(j))/dy ) ;
            psi(4) = fphi( nm, (pe%vertex4%coor(1)-x(i) )/dx,(pe%vertex4%coor(2)-y(j))/dy ) ;

            call get_vector_b( vert_temp(1:4,1:2),psi(1:4),b(1:3) )
            call solve17(temp_L(1:3,1:3),temp_U(1:3,1:3),b(1:3),aa(1:3),3)
        elseif(n_moment==6)then
            psi(1) = fphi( nm, (pe%vertex1%coor(1)-x(i) )/dx,(pe%vertex1%coor(2)-y(j))/dy ) ;
            psi(2) = fphi( nm, (pe%vertex2%coor(1)-x(i) )/dx,(pe%vertex2%coor(2)-y(j))/dy ) ;
            psi(3) = fphi( nm, (pe%vertex3%coor(1)-x(i) )/dx,(pe%vertex3%coor(2)-y(j))/dy ) ;
            psi(4) = fphi( nm, (pe%vertex4%coor(1)-x(i) )/dx,(pe%vertex4%coor(2)-y(j))/dy ) ;
            psi(5) = fphi( nm, (pe%vertex5%coor(1)-x(i) )/dx,(pe%vertex5%coor(2)-y(j))/dy ) ;
            psi(6) = fphi( nm, (pe%vertex6%coor(1)-x(i) )/dx,(pe%vertex6%coor(2)-y(j))/dy ) ;
            psi(7) = fphi( nm, (pe%vertex7%coor(1)-x(i) )/dx,(pe%vertex7%coor(2)-y(j))/dy ) ;
            psi(8) = fphi( nm, (pe%vertex8%coor(1)-x(i) )/dx,(pe%vertex8%coor(2)-y(j))/dy ) ;
            psi(9) = fphi( nm, (pe%vertex9%coor(1)-x(i) )/dx,(pe%vertex9%coor(2)-y(j))/dy ) ;

            call get_vector2_b( vert_temp(1:9,1:5),psi(1:9),b(1:6) )
            call solve17(temp_L(1:6,1:6),temp_U(1:6,1:6),b(1:6),aa(1:6),6)
        endif


        !if(time>0.8)write(99,*) i,j,nm,aa(1:6)

    endif

    end subroutine interpolate_gauss_gl