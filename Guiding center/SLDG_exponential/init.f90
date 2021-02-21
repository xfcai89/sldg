
    subroutine init
    implicit none
    real :: utemp
    integer :: kk,lx,ly
    integer :: k


    dx = (xright-xleft)/nx
    dy = (ytop-ybottom)/ny
    !******************************************************
    !subroutine mesh_generator
    !implicit none
    ! x(i) ----> center
    ! y(j) ----> center
    DO I = 1 - nghost, nx + nghost
        X(I) = XLEFT + (I-0.5) * DX
    ENDDO

    DO J = 1 - nghost , NY + nghost
        Y(J) = Ybottom + (J-0.5) * DY
    ENDDO

    DO I = 1 - nghost, nx + 1  + nghost
        Xgrid(I) = XLEFT + (I-1.) * DX
    ENDDO

    DO J =  1-nghost, NY+ 1 + nghost
        ygrid(J) = Ybottom + (J-1.) * DY
    ENDDO

    do i = 1,nx+1
        do j = 1,ny+1
            vertex(i,j)%coor(1) = xleft + (i-1) * dx
            vertex(i,j)%coor(2) = ybottom + (j-1)*dy
        enddo
    enddo

    do i = 1,nx
        do j = 1,ny+1
            nodex(i,j)%coor(1) = xleft + (i-0.5) * dx
            nodex(i,j)%coor(2) = ybottom + (j-1)*dy
        enddo
    enddo

    do i = 1,nx+1
        do j = 1,ny
            nodey(i,j)%coor(1) = xleft + (i-1) * dx
            nodey(i,j)%coor(2) = ybottom + (j-0.5)*dy
        enddo
    enddo

    do i = 1,nx
        do j = 1,ny
            nodec(i,j)%coor(1) = xleft + (i-0.5) * dx
            nodec(i,j)%coor(2) = ybottom + (j-0.5)*dy
        enddo
    enddo

    do i = 1 , nx
        do j = 1 ,ny
            do k = 0 , iexprk
                element(i,j,k)%vertex1 => vertex(i,j)
                element(i,j,k)%vertex2 => vertex(i+1,j)
                element(i,j,k)%vertex3 => vertex(i+1,j+1)
                element(i,j,k)%vertex4 => vertex(i,j+1)

                if(n_moment==6)then
                    element(i,j,k)%vertex5 => nodex(i,j)
                    element(i,j,k)%vertex6 => nodey(i+1,j)
                    element(i,j,k)%vertex7 => nodex(i,j+1)
                    element(i,j,k)%vertex8 => nodey(i,j)
                    element(i,j,k)%vertex9 => nodec(i,j)
                endif

            enddo

        enddo
    enddo
    !end subroutine mesh_generator
    !******************************************************

    do i = 1 , nx
        do j = 1 , ny
            do kk = 1 , n_moment
                utemp =0.0
                do lx=1,6
                    do ly=1,6
                        utemp = utemp +exact(x(i)+xg(lx)*dx ,y(j)+xg(ly)*dy,0.) *fphi(kk,xg(lx),xg(ly) )*wg(lx)*wg(ly)
                    enddo
                enddo
                element(i,j,0)%umodal(kk)= utemp * ai(kk)

            enddo
        enddo
    enddo

    end subroutine init
    !*******************************************************************   
    real function exact(x,y,t)
    implicit none
    real, intent(in) :: x,y,t

    exact = -2.* sin(x) * sin(y)


    return
    end function exact