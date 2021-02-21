    subroutine velocity_nodes
    implicit none
    real :: t1,t2,t3,t4
    real :: xrg,yrg

    !*******************************************************************
    do i = 1,nx+1
        do j = 1 ,ny+1
            xrg = x(i-1)+0.5*dx
            yrg = y(j-1)+0.5*dy
            t1 = polynomial2d( phi_x(i-1,j-1,1:n_moment_LDG),xrg,x(i-1),dx,yrg,y(j-1),dy,n_moment_LDG )
            xrg = x(i)-0.5*dx
            yrg = y(j-1)+0.5*dy
            t2 = polynomial2d( phi_x(i,j-1,1:n_moment_LDG),xrg,x(i),dx,yrg,y(j-1),dy,n_moment_LDG )
            xrg = x(i-1)+0.5*dx
            yrg = y(j)-0.5*dy
            t3 = polynomial2d( phi_x(i-1,j,1:n_moment_LDG),xrg,x(i-1),dx,yrg,y(j),dy,n_moment_LDG )
            xrg = x(i)-0.5*dx
            yrg = y(j)-0.5*dy
            t4 = polynomial2d( phi_x(i,j,1:n_moment_LDG),xrg,x(i),dx,yrg,y(j),dy,n_moment_LDG )
            phi_x_v(i,j,io) = (t1+t2+t3+t4)*0.25
            !*************************************************
            !*************************************************
            xrg = x(i-1)+0.5*dx
            yrg = y(j-1)+0.5*dy
            t1 = polynomial2d( phi_y(i-1,j-1,1:n_moment_LDG),xrg,x(i-1),dx,yrg,y(j-1),dy,n_moment_LDG )
            xrg = x(i)-0.5*dx
            yrg = y(j-1)+0.5*dy
            t2 = polynomial2d( phi_y(i,j-1,1:n_moment_LDG),xrg,x(i),dx,yrg,y(j-1),dy,n_moment_LDG )
            xrg = x(i-1)+0.5*dx
            yrg = y(j)-0.5*dy
            t3 = polynomial2d( phi_y(i-1,j,1:n_moment_LDG),xrg,x(i-1),dx,yrg,y(j),dy,n_moment_LDG )
            xrg = x(i)-0.5*dx
            yrg = y(j)-0.5*dy
            t4 = polynomial2d( phi_y(i,j,1:n_moment_LDG),xrg,x(i),dx,yrg,y(j),dy,n_moment_LDG )
            phi_y_v(i,j,io) = (t1+t2+t3+t4)*0.25
            
    
            
        enddo
    enddo
    
    !*************************************
    !*************************************
    !*************************************
    do i = 1,nx
        do j = 1,ny+1
            xrg = x(i) 
            yrg = y(j-1)+0.5*dy
            t1 = polynomial2d( phi_x(i,j-1,1:n_moment_LDG),xrg,x(i),dx,yrg,y(j-1),dy,n_moment_LDG )
            xrg = x(i) 
            yrg = y(j)-0.5*dy
            t2 = polynomial2d( phi_x(i,j,1:n_moment_LDG),xrg,x(i),dx,yrg,y(j),dy,n_moment_LDG )
            phi_x_dex(i,j,io) = (t1+t2)*0.5
            !*************************************************
            !*************************************************
            xrg = x(i) 
            yrg = y(j-1)+0.5*dy
            t1 = polynomial2d( phi_y(i,j-1,1:n_moment_LDG),xrg,x(i),dx,yrg,y(j-1),dy,n_moment_LDG )
            xrg = x(i)
            yrg = y(j)-0.5*dy
            t2 = polynomial2d( phi_y(i,j,1:n_moment_LDG),xrg,x(i),dx,yrg,y(j),dy,n_moment_LDG )
            phi_y_dex(i,j,io) = (t1+t2)*0.5
        enddo
    enddo

    do i = 1,nx+1
        do j = 1,ny
            xrg = x(i-1)+0.5*dx
            yrg = y(j)
            t1 = polynomial2d( phi_x(i-1,j,1:n_moment_LDG),xrg,x(i-1),dx,yrg,y(j),dy,n_moment_LDG )
            xrg = x(i)-0.5*dx
            yrg = y(j)
            t2 = polynomial2d( phi_x(i,j,1:n_moment_LDG),xrg,x(i),dx,yrg,y(j),dy,n_moment_LDG )
            phi_x_dey(i,j,io) = (t1+t2)*0.5
            !*************************************************
            !*************************************************
            xrg = x(i-1)+0.5*dx
            yrg = y(j)
            t1 = polynomial2d( phi_y(i-1,j,1:n_moment_LDG),xrg,x(i-1),dx,yrg,y(j),dy,n_moment_LDG )
            xrg = x(i)-0.5*dx
            yrg = y(j)
            t2 = polynomial2d( phi_y(i,j,1:n_moment_LDG),xrg,x(i),dx,yrg,y(j),dy,n_moment_LDG )
            phi_y_dey(i,j,io) = (t1+t2)*0.5
        enddo
    enddo

    do i = 1,nx
        do j = 1,ny             
            phi_x_c(i,j,io) = polynomial2d( phi_x(i,j,1:n_moment_LDG),x(i),x(i),dx,y(j),y(j),dy,n_moment_LDG )
            phi_y_c(i,j,io) = polynomial2d( phi_y(i,j,1:n_moment_LDG),x(i),x(i),dx,y(j),y(j),dy,n_moment_LDG )
        enddo
    enddo

    end subroutine velocity_nodes