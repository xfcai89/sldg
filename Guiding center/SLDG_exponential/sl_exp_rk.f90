    subroutine sl_exp_rk
    implicit none
    integer :: lx,ly,ic,jc
    real :: xrg,yrg

    integer :: k,ib

    !n_moment_2d = 3
    !n_moment_LDG = 3
    if(iadditive(io)==0 )then
        do i = 1 , nx
            do j = 1 , ny
                temp11(i,j,1:n_moment) = element(i,j,io)%umodal(1:n_moment)
            enddo
        enddo
        call Poisson_2D_LDG
    endif

    !open(125,file='temp_11.plt')
    !write(125,*)'zone    ','i=',nx*3,',    j=',ny*3
    !
    !DO  J=1,NY
    !    do jc = 1,3
    !        DO  I=1 ,NX
    !            do ic = 1,3
    !                xrg =  x(i)+ (ic-2.) *dx*0.25
    !                yrg =  y(j)+ (jc-2.) *dy*0.25
    !                WRITE(125,*) xrg,yrg,  polynomial2d( phi_x(i,j,1:n_moment),xrg,x(i),dx,yrg,y(j),dy,n_moment )
    !            enddo
    !        enddo
    !    enddo
    !enddo
    !
    !
    !CLOSE(125)
    !
    !pause
    !
    do i = 1,nx
        do j = 1,ny
            do k = 1,n_moment_LDG
                phi_x_io(i,j,k,io) =  phi_x(i,j,k)
                phi_y_io(i,j,k,io) =  phi_y(i,j,k)
            enddo
        enddo
    enddo

    ! boundary condition for velocity field
    do i = 1 , nx
        do ib = 1 , nghost
            phi_x_io(i,1-ib,1:n_moment_LDG,io) = phi_x_io(i,ny+1-ib,1:n_moment_LDG,io)
            phi_x_io(i,ny+ib,1:n_moment_LDG,io) = phi_x_io(i,0+ib,1:n_moment_LDG,io)

            phi_y_io(i,1-ib,1:n_moment_LDG,io) = phi_y_io(i,ny+1-ib,1:n_moment_LDG,io)
            phi_y_io(i,ny+ib,1:n_moment_LDG,io) = phi_y_io(i,0+ib,1:n_moment_LDG,io)
        enddo
    enddo

    do j = 1 - nghost , ny + nghost
        do ib = 1 , nghost
            phi_x_io(1-ib,j,1:n_moment_LDG,io) = phi_x_io(nx+1-ib,j,1:n_moment_LDG,io)
            phi_x_io(nx+ib,j,1:n_moment_LDG,io) = phi_x_io(0+ib,j,1:n_moment_LDG,io)

            phi_y_io(1-ib,j,1:n_moment_LDG,io) = phi_y_io(nx+1-ib,j,1:n_moment_LDG,io)
            phi_y_io(nx+ib,j,1:n_moment_LDG,io) = phi_y_io(0+ib,j,1:n_moment_LDG,io)
        enddo
    enddo
    !*************************


    if(io==0) call setdt

    !call RK_to_upstream

    call velocity_nodes

    !call local_rk1
    call local_rk4
    call local_rk4_gauss

    call get_intersections_outersegments
    call get_innersegments

    call get_integral

    end subroutine sl_exp_rk
    !********************************************
    subroutine local_rk1
    implicit none

    do i = 1,nx+1
        do j = 1,ny+1
            vertex_star( i,j )%coor(1) = vertex( i,j )%coor(1) +sign_rho* phi_y_v(i,j,io)* dt

            vertex_star( i,j )%coor(2) = vertex( i,j )%coor(2) -sign_rho* phi_x_v(i,j,io) * dt


            vertex_star( i,j )%id(1) = ceiling( (vertex_star( i,j )%coor(1)-xleft)/dx )
            vertex_star( i,j )%id(2) = ceiling( (vertex_star( i,j )%coor(2)-ybottom)/dy )
        enddo
    enddo

    do i = 1,nx
        do j = 1,ny+1
            !call RK2D( nodex( i,j )%coor(1:2),nodex_star( i,j )%coor(1:2),time,dt,-1. )
            nodex_star( i,j )%coor(1) = nodex( i,j )%coor(1) + sign_rho*phi_y_dex(i,j,io) * dt

            nodex_star( i,j )%coor(2) = nodex( i,j )%coor(2) - sign_rho*phi_x_dex(i,j,io) * dt
        enddo
    enddo

    do i = 1,nx+1
        do j = 1,ny
            !call RK2D( nodey( i,j )%coor(1:2),nodey_star( i,j )%coor(1:2),time,dt,-1. )
            nodey_star( i,j )%coor(1) = nodey( i,j )%coor(1) +sign_rho*phi_y_dey(i,j,io) * dt

            nodey_star( i,j )%coor(2) = nodey( i,j )%coor(2) -sign_rho*phi_x_dey(i,j,io) * dt
        enddo
    enddo

    do i = 1,nx
        do j = 1,ny
            !call RK2D( nodec( i,j )%coor(1:2),nodec_star( i,j )%coor(1:2),time,dt,-1. )
            nodec_star( i,j )%coor(1) = nodec( i,j )%coor(1)+sign_rho* phi_y_c(i,j,io) * dt

            nodec_star( i,j )%coor(2) = nodec( i,j )%coor(2) -sign_rho*phi_x_c(i,j,io) * dt
        enddo
    enddo

    ! get the vertexes of face_lr
    do j = 1 , ny+1
        do i = 1 , nx
            face_lr(i,j)%point_origin = vertex_star(i,j)
            face_lr(i,j)%point_end = vertex_star(i+1,j)
            ! quadratic-curved
            face_lr(i,j)%point_midt = nodex_star(i,j)
        enddo
    enddo

    ! get the vertexes of face_bt
    do i = 1 , nx+1
        do j = 1 , ny
            face_bt(i,j)%point_origin = vertex_star(i,j)
            face_bt(i,j)%point_end = vertex_star(i,j+1)
            ! quadratic-curved
            face_bt(i,j)%point_midt = nodey_star(i,j)
        enddo
    enddo


    end subroutine local_rk1
    !*************************************************************
    subroutine local_rk3
    implicit none
    real :: vel_x, vel_y
    real :: vx1(1:2),vx2(1:2)
    integer :: idx,idy
    real :: t1,t2,t3

    integer :: ii,io1

    real :: phix(0:3),phiy(0:3)

    do i = 1,nx+1
        do j = 1,ny+1
            if(io == 0)then
                vel_x = bt_con(2,1)*phi_y_v(i,j,0)
            elseif( io == 1 )then
                vel_x = bt_con(3,1)*phi_y_v(i,j,0)+bt_con(3,2)*phi_y_v(i,j,1)
            elseif(io==2)then
                vel_x = bt_con(4,1)*phi_y_v(i,j,0)+bt_con(4,2)*phi_y_v(i,j,1)+bt_con(4,3)*phi_y_v(i,j,2)
            elseif(io==3)then
                vel_x = bt_con(5,1)*phi_y_v(i,j,0)+bt_con(5,2)*phi_y_v(i,j,1)+bt_con(5,3)*phi_y_v(i,j,2)
            endif
            vel_x = -sign_rho*vel_x

            if(io == 0)then
                vel_y = bt_con(2,1)*phi_x_v(i,j,0)
            elseif(io==1)then
                vel_y = bt_con(3,1)*phi_x_v(i,j,0)+bt_con(3,2)*phi_x_v(i,j,1)
            elseif(io==2)then
                vel_y = bt_con(4,1)*phi_x_v(i,j,0)+bt_con(4,2)*phi_x_v(i,j,1)+bt_con(4,3)*phi_x_v(i,j,2)
            elseif( io==3 )then
                vel_y = bt_con(5,1)*phi_x_v(i,j,0)+bt_con(5,2)*phi_x_v(i,j,1)+bt_con(5,3)*phi_x_v(i,j,2)
            endif
            vel_y = sign_rho*vel_y

            vx1(1) = vertex( i,j )%coor(1) - vel_x * dt
            vx1(2) = vertex( i,j )%coor(2) - vel_y * dt
                      

            !***********************************************
            idx = ceiling( (vx1(1)-xleft)/dx )
            idy = ceiling( (vx1(2)-ybottom)/dy )

            io1=io
            if(io==3) io1= io-1
            call get_velocity_nodes(idx,idy,io1,vx1(1:2),phix,phiy )
            
            vel_x =0.
            do ii = 0,io1
                vel_x = vel_x + bt_con(io+2,ii+1)*phiy(ii)
            enddo
            vel_x = -sign_rho*vel_x

            vel_y =0.
            do ii = 0,io1
                vel_y = vel_y + bt_con(io+2,ii+1)*phix(ii)
            enddo
            vel_y = sign_rho*vel_y
            
            vx2(1) = 0.75*vertex( i,j )%coor(1) + 0.25*( vx1(1) - vel_x*dt  )
            vx2(2) = 0.75*vertex( i,j )%coor(2) + 0.25*( vx1(2) - vel_y * dt  )
            
            
            !********************************************************
            idx = ceiling( (vx2(1)-xleft)/dx )
            idy = ceiling( (vx2(2)-ybottom)/dy )
            
            io1=io
            if(io==3) io1= io-1
            call get_velocity_nodes(idx,idy,io1,vx2(1:2),phix,phiy )
            
            vel_x =0.
            do ii = 0,io1
                vel_x = vel_x + bt_con(io+2,ii+1)*phiy(ii)
            enddo
            vel_x = -sign_rho*vel_x

            vel_y =0.
            do ii = 0,io1
                vel_y = vel_y + bt_con(io+2,ii+1)*phix(ii)
            enddo
            vel_y = sign_rho*vel_y
            !
            vertex_star( i,j )%coor(1) = 1./3.*vertex( i,j )%coor(1) + 2./3.*( vx2(1) - vel_x*dt )
            vertex_star( i,j )%coor(2) = 1./3.*vertex( i,j )%coor(2) + 2./3.*( vx2(2) - vel_y * dt )

            vertex_star( i,j )%id(1) = ceiling( (vertex_star( i,j )%coor(1)-xleft)/dx )
            vertex_star( i,j )%id(2) = ceiling( (vertex_star( i,j )%coor(2)-ybottom)/dy )
        enddo
    enddo

    do i = 1,nx
        do j = 1,ny+1

            if(io == 0)then
                vel_x = bt_con(2,1)*phi_y_dex(i,j,0)
            elseif( io == 1 )then
                vel_x = bt_con(3,1)*phi_y_dex(i,j,0)+ bt_con(3,2)*phi_y_dex(i,j,1)
            elseif(io==2)then
                vel_x = bt_con(4,1)*phi_y_dex(i,j,0)+ bt_con(4,2)*phi_y_dex(i,j,1)+bt_con(4,3)*phi_y_dex(i,j,2)
            elseif(io==3)then
                vel_x = bt_con(5,1)*phi_y_dex(i,j,0)+ bt_con(5,2)*phi_y_dex(i,j,1)+bt_con(5,3)*phi_y_dex(i,j,2)
            endif
            vel_x = -sign_rho*vel_x

            vx1(1) = nodex( i,j )%coor(1) - vel_x * dt

            if(io == 0)then
                vel_y=bt_con(2,1)*phi_x_dex(i,j,0)
            elseif( io == 1 )then
                vel_y=bt_con(3,1)*phi_x_dex(i,j,0)+bt_con(3,2)*phi_x_dex(i,j,1)
            elseif(io==2)then
                vel_y=bt_con(4,1)*phi_x_dex(i,j,0)+bt_con(4,2)*phi_x_dex(i,j,1)+bt_con(4,3)*phi_x_dex(i,j,2)
            elseif(io==3)then
                vel_y=bt_con(5,1)*phi_x_dex(i,j,0)+bt_con(5,2)*phi_x_dex(i,j,1)+bt_con(5,3)*phi_x_dex(i,j,2)
            endif
            vel_y = sign_rho*vel_y

            vx1(2) = nodex( i,j )%coor(2) - vel_y * dt

            !***********************************************
            idx = ceiling( (vx1(1)-xleft)/dx )
            idy = ceiling( (vx1(2)-ybottom)/dy )

            io1=io
            if(io==3) io1= io-1
            call get_velocity_nodes(idx,idy,io1,vx1(1:2),phix,phiy )
            
            vel_x =0.
            do ii = 0,io1
                vel_x = vel_x + bt_con(io+2,ii+1)*phiy(ii)
            enddo
            vel_x = -sign_rho*vel_x

            vel_y =0.
            do ii = 0,io1
                vel_y = vel_y + bt_con(io+2,ii+1)*phix(ii)
            enddo
            vel_y = sign_rho*vel_y
            !
            vx2(1) = 0.75*nodex( i,j )%coor(1) + 0.25*( vx1(1) - vel_x*dt  )
            vx2(2) = 0.75*nodex( i,j )%coor(2) + 0.25*( vx1(2) - vel_y * dt  )
            !********************************************************
            idx = ceiling( (vx2(1)-xleft)/dx )
            idy = ceiling( (vx2(2)-ybottom)/dy )

            io1=io
            if(io==3) io1= io-1
            call get_velocity_nodes(idx,idy,io1,vx2(1:2),phix,phiy )
            
            vel_x =0.
            do ii = 0,io1
                vel_x = vel_x + bt_con(io+2,ii+1)*phiy(ii)
            enddo
            vel_x = -sign_rho*vel_x

            vel_y =0.
            do ii = 0,io1
                vel_y = vel_y + bt_con(io+2,ii+1)*phix(ii)
            enddo
            vel_y = sign_rho*vel_y

            nodex_star( i,j )%coor(1) = 1./3.*nodex( i,j )%coor(1) + 2./3.*( vx2(1) - vel_x*dt )
            nodex_star( i,j )%coor(2) = 1./3.*nodex( i,j )%coor(2) + 2./3.*( vx2(2) - vel_y * dt )

        enddo
    enddo

    do i = 1,nx+1
        do j = 1,ny
            if(io == 0)then
                vel_x = bt_con(2,1)*phi_y_dey(i,j,0)
            elseif( io == 1 )then
                vel_x = bt_con(3,1)*phi_y_dey(i,j,0)+bt_con(3,2)*phi_y_dey(i,j,1)
            elseif(io==2)then
                vel_x = bt_con(4,1)*phi_y_dey(i,j,0)+bt_con(4,2)*phi_y_dey(i,j,1)+bt_con(4,3)*phi_y_dey(i,j,2)
            elseif(io==3)then
                vel_x = bt_con(5,1)*phi_y_dey(i,j,0)+bt_con(5,2)*phi_y_dey(i,j,1)+bt_con(5,3)*phi_y_dey(i,j,2)
            endif
            vel_x = -sign_rho*vel_x

            vx1(1) = nodey( i,j )%coor(1) - vel_x * dt

            if(io == 0)then
                vel_y = bt_con(2,1)*phi_x_dey(i,j,0)
            elseif( io == 1 )then
                vel_y = bt_con(3,1)*phi_x_dey(i,j,0)+bt_con(3,2)*phi_x_dey(i,j,1)
            elseif(io==2)then
                vel_y = bt_con(4,1)*phi_x_dey(i,j,0)+bt_con(4,2)*phi_x_dey(i,j,1)+bt_con(4,3)*phi_x_dey(i,j,2)
            elseif(io==3)then
                vel_y = bt_con(5,1)*phi_x_dey(i,j,0)+bt_con(5,2)*phi_x_dey(i,j,1)+bt_con(5,3)*phi_x_dey(i,j,2)
            endif
            vel_y = sign_rho*vel_y

            vx1(2) = nodey( i,j )%coor(2) - vel_y * dt

            !***********************************************
            idx = ceiling( (vx1(1)-xleft)/dx )
            idy = ceiling( (vx1(2)-ybottom)/dy )
            io1=io
            if(io==3) io1= io-1
            call get_velocity_nodes(idx,idy,io1,vx1(1:2),phix,phiy )
            
            vel_x =0.
            do ii = 0,io1
                vel_x = vel_x + bt_con(io+2,ii+1)*phiy(ii)
            enddo
            vel_x = -sign_rho*vel_x

            vel_y =0.
            do ii = 0,io1
                vel_y = vel_y + bt_con(io+2,ii+1)*phix(ii)
            enddo
            vel_y = sign_rho*vel_y
            !
            vx2(1) = 0.75*nodey( i,j )%coor(1) + 0.25*( vx1(1) - vel_x*dt  )
            vx2(2) = 0.75*nodey( i,j )%coor(2) + 0.25*( vx1(2) - vel_y*dt  )
            !********************************************************
            idx = ceiling( (vx2(1)-xleft)/dx )
            idy = ceiling( (vx2(2)-ybottom)/dy )

            io1=io
            if(io==3) io1= io-1
            call get_velocity_nodes(idx,idy,io1,vx2(1:2),phix,phiy )
            
            vel_x =0.
            do ii = 0,io1
                vel_x = vel_x + bt_con(io+2,ii+1)*phiy(ii)
            enddo
            vel_x = -sign_rho*vel_x

            vel_y =0.
            do ii = 0,io1
                vel_y = vel_y + bt_con(io+2,ii+1)*phix(ii)
            enddo
            vel_y = sign_rho*vel_y
            !
            nodey_star( i,j )%coor(1) = 1./3.*nodey( i,j )%coor(1) + 2./3.*( vx2(1) - vel_x*dt )
            nodey_star( i,j )%coor(2) = 1./3.*nodey( i,j )%coor(2) + 2./3.*( vx2(2) - vel_y * dt )

        enddo
    enddo

    do i = 1,nx
        do j = 1,ny

            if(io == 0)then
                vel_x = bt_con(2,1)*phi_y_c(i,j,0)
            elseif( io == 1 )then
                vel_x = bt_con(3,1)*phi_y_c(i,j,0)+bt_con(3,2)*phi_y_c(i,j,1)
            elseif(io==2)then
                vel_x = bt_con(4,1)*phi_y_c(i,j,0)+bt_con(4,2)*phi_y_c(i,j,1)+bt_con(4,3)*phi_y_c(i,j,2)
            elseif(io==3)then
                vel_x = bt_con(5,1)*phi_y_c(i,j,0)+bt_con(5,2)*phi_y_c(i,j,1)+bt_con(5,3)*phi_y_c(i,j,2)
            endif
            vel_x = -sign_rho*vel_x

            vx1(1) = nodec( i,j )%coor(1) - vel_x * dt

            if(io == 0)then
                vel_y = bt_con(2,1)*phi_x_c(i,j,0)
            elseif( io == 1 )then
                vel_y = bt_con(3,1)*phi_x_c(i,j,0)+bt_con(3,2)*phi_x_c(i,j,1)
            elseif(io==2)then
                vel_y = bt_con(4,1)*phi_x_c(i,j,0)+bt_con(4,2)*phi_x_c(i,j,1)+bt_con(4,3)*phi_x_c(i,j,2)
            elseif(io==3)then
                vel_y = bt_con(5,1)*phi_x_c(i,j,0)+bt_con(5,2)*phi_x_c(i,j,1)+bt_con(5,3)*phi_x_c(i,j,2)
            endif
            vel_y = sign_rho*vel_y

            vx1(2) = nodec( i,j )%coor(2) - vel_y * dt

            !***********************************************
            idx = ceiling( (vx1(1)-xleft)/dx )
            idy = ceiling( (vx1(2)-ybottom)/dy )

            io1=io
            if(io==3) io1= io-1
            call get_velocity_nodes(idx,idy,io1,vx1(1:2),phix,phiy )
            
            vel_x =0.
            do ii = 0,io1
                vel_x = vel_x + bt_con(io+2,ii+1)*phiy(ii)
            enddo
            vel_x = -sign_rho*vel_x

            vel_y =0.
            do ii = 0,io1
                vel_y = vel_y + bt_con(io+2,ii+1)*phix(ii)
            enddo
            vel_y = sign_rho*vel_y 

            vx2(1) = 0.75*nodec( i,j )%coor(1) + 0.25*( vx1(1) - vel_x*dt  )
            vx2(2) = 0.75*nodec( i,j )%coor(2) + 0.25*( vx1(2) - vel_y * dt  )
            !********************************************************
            idx = ceiling( (vx2(1)-xleft)/dx )
            idy = ceiling( (vx2(2)-ybottom)/dy )

            io1=io
            if(io==3) io1= io-1
            call get_velocity_nodes(idx,idy,io1,vx2(1:2),phix,phiy )
            
            vel_x =0.
            do ii = 0,io1
                vel_x = vel_x + bt_con(io+2,ii+1)*phiy(ii)
            enddo
            vel_x = -sign_rho*vel_x

            vel_y =0.
            do ii = 0,io1
                vel_y = vel_y + bt_con(io+2,ii+1)*phix(ii)
            enddo
            vel_y = sign_rho*vel_y 
            
            nodec_star( i,j )%coor(1) = 1./3.*nodec( i,j )%coor(1) + 2./3.*( vx2(1) - vel_x*dt )
            nodec_star( i,j )%coor(2) = 1./3.*nodec( i,j )%coor(2) + 2./3.*( vx2(2) - vel_y * dt )

        enddo
    enddo

    ! get the vertexes of face_lr
    do j = 1 , ny+1
        do i = 1 , nx
            face_lr(i,j)%point_origin = vertex_star(i,j)
            face_lr(i,j)%point_end = vertex_star(i+1,j)
            ! quadratic-curved
            face_lr(i,j)%point_midt = nodex_star(i,j)
        enddo
    enddo

    ! get the vertexes of face_bt
    do i = 1 , nx+1
        do j = 1 , ny
            face_bt(i,j)%point_origin = vertex_star(i,j)
            face_bt(i,j)%point_end = vertex_star(i,j+1)
            ! quadratic-curved
            face_bt(i,j)%point_midt = nodey_star(i,j)
        enddo
    enddo


    end subroutine local_rk3

    !*******************************************************************
    subroutine local_rk4
    implicit none
    real :: vel_x, vel_y
    real :: vx1(1:2),vx2(1:2),vx3(1:2)
    integer :: idx,idy
    real :: t1,t2,t3

    integer :: ii,io1

    real :: phix(0:3),phiy(0:3)

    do i = 1,nx+1
        do j = 1,ny+1
            if(io == 0)then
                vel_x = bt_con(2,1)*phi_y_v(i,j,0)
            elseif( io == 1 )then
                vel_x = bt_con(3,1)*phi_y_v(i,j,0)+bt_con(3,2)*phi_y_v(i,j,1)
            elseif(io==2)then
                vel_x = bt_con(4,1)*phi_y_v(i,j,0)+bt_con(4,2)*phi_y_v(i,j,1)+bt_con(4,3)*phi_y_v(i,j,2)
            elseif(io==3)then
                vel_x = bt_con(5,1)*phi_y_v(i,j,0)+bt_con(5,2)*phi_y_v(i,j,1)+bt_con(5,3)*phi_y_v(i,j,2)
            endif
            vel_x = -sign_rho*vel_x

            if(io == 0)then
                vel_y = bt_con(2,1)*phi_x_v(i,j,0)
            elseif(io==1)then
                vel_y = bt_con(3,1)*phi_x_v(i,j,0)+bt_con(3,2)*phi_x_v(i,j,1)
            elseif(io==2)then
                vel_y = bt_con(4,1)*phi_x_v(i,j,0)+bt_con(4,2)*phi_x_v(i,j,1)+bt_con(4,3)*phi_x_v(i,j,2)
            elseif( io==3 )then
                vel_y = bt_con(5,1)*phi_x_v(i,j,0)+bt_con(5,2)*phi_x_v(i,j,1)+bt_con(5,3)*phi_x_v(i,j,2)
            endif
            vel_y = sign_rho*vel_y

            vx1(1) = vertex( i,j )%coor(1) - 0.5*vel_x * dt
            vx1(2) = vertex( i,j )%coor(2) - 0.5*vel_y * dt

            !***********************************************
            idx = ceiling( (vx1(1)-xleft)/dx )
            idy = ceiling( (vx1(2)-ybottom)/dy )

            io1=io
            if(io==3) io1= io-1
            call get_velocity_nodes(idx,idy,io1,vx1(1:2),phix,phiy )
            
            vel_x =0.
            do ii = 0,io1
                vel_x = vel_x + bt_con(io+2,ii+1)*phiy(ii)
            enddo
            vel_x = -sign_rho*vel_x

            vel_y =0.
            do ii = 0,io1
                vel_y = vel_y + bt_con(io+2,ii+1)*phix(ii)
            enddo
            vel_y = sign_rho*vel_y
            
            vx2(1) = vertex( i,j )%coor(1) - 0.5*vel_x * dt
            vx2(2) = vertex( i,j )%coor(2) - 0.5*vel_y * dt
            !********************************************************
            idx = ceiling( (vx2(1)-xleft)/dx )
            idy = ceiling( (vx2(2)-ybottom)/dy )
            
            io1=io
            if(io==3) io1= io-1
            call get_velocity_nodes(idx,idy,io1,vx2(1:2),phix,phiy )
            
            vel_x =0.
            do ii = 0,io1
                vel_x = vel_x + bt_con(io+2,ii+1)*phiy(ii)
            enddo
            vel_x = -sign_rho*vel_x

            vel_y =0.
            do ii = 0,io1
                vel_y = vel_y + bt_con(io+2,ii+1)*phix(ii)
            enddo
            vel_y = sign_rho*vel_y
            !
            vx3(1) = vertex( i,j )%coor(1) - vel_x * dt
            vx3(2) = vertex( i,j )%coor(2) - vel_y * dt
            !*********
            !*********
            !********************************************************
            idx = ceiling( (vx3(1)-xleft)/dx )
            idy = ceiling( (vx3(2)-ybottom)/dy )
            
            io1=io
            if(io==3) io1= io-1
            call get_velocity_nodes(idx,idy,io1,vx3(1:2),phix,phiy )
            
            vel_x =0.
            do ii = 0,io1
                vel_x = vel_x + bt_con(io+2,ii+1)*phiy(ii)
            enddo
            vel_x = -sign_rho*vel_x

            vel_y =0.
            do ii = 0,io1
                vel_y = vel_y + bt_con(io+2,ii+1)*phix(ii)
            enddo
            vel_y = sign_rho*vel_y            
            
            vertex_star( i,j )%coor(1) =  1./3.*( -vertex( i,j )%coor(1)+vx1(1)+2.*vx2(1)+vx3(1) )  - 1./6.*vel_x*dt 
            vertex_star( i,j )%coor(2) =  1./3.*( -vertex( i,j )%coor(2)+vx1(2)+2.*vx2(2)+vx3(2) )  - 1./6.*vel_y*dt 
 

            vertex_star( i,j )%id(1) = ceiling( (vertex_star( i,j )%coor(1)-xleft)/dx )
            vertex_star( i,j )%id(2) = ceiling( (vertex_star( i,j )%coor(2)-ybottom)/dy )
        enddo
    enddo

    do i = 1,nx
        do j = 1,ny+1

            if(io == 0)then
                vel_x = bt_con(2,1)*phi_y_dex(i,j,0)
            elseif( io == 1 )then
                vel_x = bt_con(3,1)*phi_y_dex(i,j,0)+ bt_con(3,2)*phi_y_dex(i,j,1)
            elseif(io==2)then
                vel_x = bt_con(4,1)*phi_y_dex(i,j,0)+ bt_con(4,2)*phi_y_dex(i,j,1)+bt_con(4,3)*phi_y_dex(i,j,2)
            elseif(io==3)then
                vel_x = bt_con(5,1)*phi_y_dex(i,j,0)+ bt_con(5,2)*phi_y_dex(i,j,1)+bt_con(5,3)*phi_y_dex(i,j,2)
            endif
            vel_x = -sign_rho*vel_x

            vx1(1) = nodex( i,j )%coor(1) - 0.5*vel_x * dt

            if(io == 0)then
                vel_y=bt_con(2,1)*phi_x_dex(i,j,0)
            elseif( io == 1 )then
                vel_y=bt_con(3,1)*phi_x_dex(i,j,0)+bt_con(3,2)*phi_x_dex(i,j,1)
            elseif(io==2)then
                vel_y=bt_con(4,1)*phi_x_dex(i,j,0)+bt_con(4,2)*phi_x_dex(i,j,1)+bt_con(4,3)*phi_x_dex(i,j,2)
            elseif(io==3)then
                vel_y=bt_con(5,1)*phi_x_dex(i,j,0)+bt_con(5,2)*phi_x_dex(i,j,1)+bt_con(5,3)*phi_x_dex(i,j,2)
            endif
            vel_y = sign_rho*vel_y

            vx1(2) = nodex( i,j )%coor(2) - 0.5*vel_y * dt

            !***********************************************
            idx = ceiling( (vx1(1)-xleft)/dx )
            idy = ceiling( (vx1(2)-ybottom)/dy )

            io1=io
            if(io==3) io1= io-1
            call get_velocity_nodes(idx,idy,io1,vx1(1:2),phix,phiy )
            
            vel_x =0.
            do ii = 0,io1
                vel_x = vel_x + bt_con(io+2,ii+1)*phiy(ii)
            enddo
            vel_x = -sign_rho*vel_x

            vel_y =0.
            do ii = 0,io1
                vel_y = vel_y + bt_con(io+2,ii+1)*phix(ii)
            enddo
            vel_y = sign_rho*vel_y
            
            vx2(1) = nodex( i,j )%coor(1) - 0.5*vel_x * dt
            vx2(2) = nodex( i,j )%coor(2) - 0.5*vel_y * dt
            !********************************************************
            idx = ceiling( (vx2(1)-xleft)/dx )
            idy = ceiling( (vx2(2)-ybottom)/dy )
            
            io1=io
            if(io==3) io1= io-1
            call get_velocity_nodes(idx,idy,io1,vx2(1:2),phix,phiy )
            
            vel_x =0.
            do ii = 0,io1
                vel_x = vel_x + bt_con(io+2,ii+1)*phiy(ii)
            enddo
            vel_x = -sign_rho*vel_x

            vel_y =0.
            do ii = 0,io1
                vel_y = vel_y + bt_con(io+2,ii+1)*phix(ii)
            enddo
            vel_y = sign_rho*vel_y
            !
            vx3(1) = nodex( i,j )%coor(1) - vel_x * dt
            vx3(2) = nodex( i,j )%coor(2) - vel_y * dt
            !*********
            !*********
            !********************************************************
            idx = ceiling( (vx3(1)-xleft)/dx )
            idy = ceiling( (vx3(2)-ybottom)/dy )
            
            io1=io
            if(io==3) io1= io-1
            call get_velocity_nodes(idx,idy,io1,vx3(1:2),phix,phiy )
            
            vel_x =0.
            do ii = 0,io1
                vel_x = vel_x + bt_con(io+2,ii+1)*phiy(ii)
            enddo
            vel_x = -sign_rho*vel_x

            vel_y =0.
            do ii = 0,io1
                vel_y = vel_y + bt_con(io+2,ii+1)*phix(ii)
            enddo
            vel_y = sign_rho*vel_y            
            
            nodex_star( i,j )%coor(1) =  1./3.*( -nodex( i,j )%coor(1)+vx1(1)+2.*vx2(1)+vx3(1) )  - 1./6.*vel_x*dt  
            nodex_star( i,j )%coor(2) =  1./3.*( -nodex( i,j )%coor(2)+vx1(2)+2.*vx2(2)+vx3(2) )  - 1./6.*vel_y*dt  
 

        enddo
    enddo

    do i = 1,nx+1
        do j = 1,ny
            if(io == 0)then
                vel_x = bt_con(2,1)*phi_y_dey(i,j,0)
            elseif( io == 1 )then
                vel_x = bt_con(3,1)*phi_y_dey(i,j,0)+bt_con(3,2)*phi_y_dey(i,j,1)
            elseif(io==2)then
                vel_x = bt_con(4,1)*phi_y_dey(i,j,0)+bt_con(4,2)*phi_y_dey(i,j,1)+bt_con(4,3)*phi_y_dey(i,j,2)
            elseif(io==3)then
                vel_x = bt_con(5,1)*phi_y_dey(i,j,0)+bt_con(5,2)*phi_y_dey(i,j,1)+bt_con(5,3)*phi_y_dey(i,j,2)
            endif
            vel_x = -sign_rho*vel_x

            vx1(1) = nodey( i,j )%coor(1) - 0.5*vel_x * dt

            if(io == 0)then
                vel_y = bt_con(2,1)*phi_x_dey(i,j,0)
            elseif( io == 1 )then
                vel_y = bt_con(3,1)*phi_x_dey(i,j,0)+bt_con(3,2)*phi_x_dey(i,j,1)
            elseif(io==2)then
                vel_y = bt_con(4,1)*phi_x_dey(i,j,0)+bt_con(4,2)*phi_x_dey(i,j,1)+bt_con(4,3)*phi_x_dey(i,j,2)
            elseif(io==3)then
                vel_y = bt_con(5,1)*phi_x_dey(i,j,0)+bt_con(5,2)*phi_x_dey(i,j,1)+bt_con(5,3)*phi_x_dey(i,j,2)
            endif
            vel_y = sign_rho*vel_y

            vx1(2) = nodey( i,j )%coor(2) - 0.5*vel_y * dt

            !***********************************************
            idx = ceiling( (vx1(1)-xleft)/dx )
            idy = ceiling( (vx1(2)-ybottom)/dy )

            io1=io
            if(io==3) io1= io-1
            call get_velocity_nodes(idx,idy,io1,vx1(1:2),phix,phiy )
            
            vel_x =0.
            do ii = 0,io1
                vel_x = vel_x + bt_con(io+2,ii+1)*phiy(ii)
            enddo
            vel_x = -sign_rho*vel_x

            vel_y =0.
            do ii = 0,io1
                vel_y = vel_y + bt_con(io+2,ii+1)*phix(ii)
            enddo
            vel_y = sign_rho*vel_y
            
            vx2(1) = nodey( i,j )%coor(1) - 0.5*vel_x * dt
            vx2(2) = nodey( i,j )%coor(2) - 0.5*vel_y * dt
            !********************************************************
            idx = ceiling( (vx2(1)-xleft)/dx )
            idy = ceiling( (vx2(2)-ybottom)/dy )
            
            io1=io
            if(io==3) io1= io-1
            call get_velocity_nodes(idx,idy,io1,vx2(1:2),phix,phiy )
            
            vel_x =0.
            do ii = 0,io1
                vel_x = vel_x + bt_con(io+2,ii+1)*phiy(ii)
            enddo
            vel_x = -sign_rho*vel_x

            vel_y =0.
            do ii = 0,io1
                vel_y = vel_y + bt_con(io+2,ii+1)*phix(ii)
            enddo
            vel_y = sign_rho*vel_y
            !
            vx3(1) = nodey( i,j )%coor(1) - vel_x * dt
            vx3(2) = nodey( i,j )%coor(2) - vel_y * dt
            !*********
            !*********
            !********************************************************
            idx = ceiling( (vx3(1)-xleft)/dx )
            idy = ceiling( (vx3(2)-ybottom)/dy )
            
            io1=io
            if(io==3) io1= io-1
            call get_velocity_nodes(idx,idy,io1,vx3(1:2),phix,phiy )
            
            vel_x =0.
            do ii = 0,io1
                vel_x = vel_x + bt_con(io+2,ii+1)*phiy(ii)
            enddo
            vel_x = -sign_rho*vel_x

            vel_y =0.
            do ii = 0,io1
                vel_y = vel_y + bt_con(io+2,ii+1)*phix(ii)
            enddo
            vel_y = sign_rho*vel_y            
            
            nodey_star( i,j )%coor(1) =  1./3.*( -nodey( i,j )%coor(1)+vx1(1)+2.*vx2(1)+vx3(1) )  - 1./6.*vel_x*dt  
            nodey_star( i,j )%coor(2) =  1./3.*( -nodey( i,j )%coor(2)+vx1(2)+2.*vx2(2)+vx3(2) )  - 1./6.*vel_y*dt  
 


        enddo
    enddo

    do i = 1,nx
        do j = 1,ny

            if(io == 0)then
                vel_x = bt_con(2,1)*phi_y_c(i,j,0)
            elseif( io == 1 )then
                vel_x = bt_con(3,1)*phi_y_c(i,j,0)+bt_con(3,2)*phi_y_c(i,j,1)
            elseif(io==2)then
                vel_x = bt_con(4,1)*phi_y_c(i,j,0)+bt_con(4,2)*phi_y_c(i,j,1)+bt_con(4,3)*phi_y_c(i,j,2)
            elseif(io==3)then
                vel_x = bt_con(5,1)*phi_y_c(i,j,0)+bt_con(5,2)*phi_y_c(i,j,1)+bt_con(5,3)*phi_y_c(i,j,2)
            endif
            vel_x = -sign_rho*vel_x

            vx1(1) = nodec( i,j )%coor(1) - 0.5*vel_x * dt

            if(io == 0)then
                vel_y = bt_con(2,1)*phi_x_c(i,j,0)
            elseif( io == 1 )then
                vel_y = bt_con(3,1)*phi_x_c(i,j,0)+bt_con(3,2)*phi_x_c(i,j,1)
            elseif(io==2)then
                vel_y = bt_con(4,1)*phi_x_c(i,j,0)+bt_con(4,2)*phi_x_c(i,j,1)+bt_con(4,3)*phi_x_c(i,j,2)
            elseif(io==3)then
                vel_y = bt_con(5,1)*phi_x_c(i,j,0)+bt_con(5,2)*phi_x_c(i,j,1)+bt_con(5,3)*phi_x_c(i,j,2)
            endif
            vel_y = sign_rho*vel_y

            vx1(2) = nodec( i,j )%coor(2) - 0.5*vel_y * dt

            !***********************************************
            idx = ceiling( (vx1(1)-xleft)/dx )
            idy = ceiling( (vx1(2)-ybottom)/dy )

            io1=io
            if(io==3) io1= io-1
            call get_velocity_nodes(idx,idy,io1,vx1(1:2),phix,phiy )
            
            vel_x =0.
            do ii = 0,io1
                vel_x = vel_x + bt_con(io+2,ii+1)*phiy(ii)
            enddo
            vel_x = -sign_rho*vel_x

            vel_y =0.
            do ii = 0,io1
                vel_y = vel_y + bt_con(io+2,ii+1)*phix(ii)
            enddo
            vel_y = sign_rho*vel_y
            
            vx2(1) = nodec( i,j )%coor(1) - 0.5*vel_x * dt
            vx2(2) = nodec( i,j )%coor(2) - 0.5*vel_y * dt
            !********************************************************
            idx = ceiling( (vx2(1)-xleft)/dx )
            idy = ceiling( (vx2(2)-ybottom)/dy )
            
            io1=io
            if(io==3) io1= io-1
            call get_velocity_nodes(idx,idy,io1,vx2(1:2),phix,phiy )
            
            vel_x =0.
            do ii = 0,io1
                vel_x = vel_x + bt_con(io+2,ii+1)*phiy(ii)
            enddo
            vel_x = -sign_rho*vel_x

            vel_y =0.
            do ii = 0,io1
                vel_y = vel_y + bt_con(io+2,ii+1)*phix(ii)
            enddo
            vel_y = sign_rho*vel_y
            !
            vx3(1) = nodec( i,j )%coor(1) - vel_x * dt
            vx3(2) = nodec( i,j )%coor(2) - vel_y * dt
            !*********
            !*********
            !********************************************************
            idx = ceiling( (vx3(1)-xleft)/dx )
            idy = ceiling( (vx3(2)-ybottom)/dy )
            
            io1=io
            if(io==3) io1= io-1
            call get_velocity_nodes(idx,idy,io1,vx3(1:2),phix,phiy )
            
            vel_x =0.
            do ii = 0,io1
                vel_x = vel_x + bt_con(io+2,ii+1)*phiy(ii)
            enddo
            vel_x = -sign_rho*vel_x

            vel_y =0.
            do ii = 0,io1
                vel_y = vel_y + bt_con(io+2,ii+1)*phix(ii)
            enddo
            vel_y = sign_rho*vel_y            
            
            nodec_star( i,j )%coor(1) =  1./3.*( -nodec( i,j )%coor(1)+vx1(1)+2.*vx2(1)+vx3(1) )  - 1./6.*vel_x*dt  
            nodec_star( i,j )%coor(2) =  1./3.*( -nodec( i,j )%coor(2)+vx1(2)+2.*vx2(2)+vx3(2) )  - 1./6.*vel_y*dt 

        enddo
    enddo

    ! get the vertexes of face_lr
    do j = 1 , ny+1
        do i = 1 , nx
            face_lr(i,j)%point_origin = vertex_star(i,j)
            face_lr(i,j)%point_end = vertex_star(i+1,j)
            ! quadratic-curved
            face_lr(i,j)%point_midt = nodex_star(i,j)
        enddo
    enddo

    ! get the vertexes of face_bt
    do i = 1 , nx+1
        do j = 1 , ny
            face_bt(i,j)%point_origin = vertex_star(i,j)
            face_bt(i,j)%point_end = vertex_star(i,j+1)
            ! quadratic-curved
            face_bt(i,j)%point_midt = nodey_star(i,j)
        enddo
    enddo


    end subroutine local_rk4