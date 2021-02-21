    subroutine advection_dat
    implicit none

    if(icase .eq. 0)then
        call smooth_sin
    endif

    do j = 1 , nel_y
        do i = 1 , nel_x
            elem(i,j)%psi(1:n_g,1:n_g) = fel(1:n_g,1:n_g,i,j)
        enddo
    enddo

    end subroutine advection_dat
    !**************************************************************************************
    subroutine smooth_sin
    implicit none
    integer :: n1,n2

    do j = 1 , nel_y
        do i = 1 , nel_x
            !
            do n2 = 1 , n_g
                do n1 = 1 ,n_g
                    fel( n1,n2,i,j ) = exact(  gx(n1,n2,i,j),gy(n1,n2,i,j) ,0. )

                enddo
            enddo
        enddo
    enddo


    end subroutine smooth_sin
    !***********************************************
    real function exact(x,y,t)
    implicit none
    real, intent(in) :: x,y,t
    real :: rb,rb0

    exact = -2.* sin(x) * sin(y)

    return
    end function exact
    !***************************************************************************************
