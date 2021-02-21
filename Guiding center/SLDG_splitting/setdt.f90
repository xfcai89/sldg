    subroutine setdt
    implicit none
    real :: temp_max1,temp_max2,rrrr

    temp_max1 = 0.
    temp_max2 = 0.
    do i = 1 , nel_x
        do j = 1 ,nel_y
            temp_max1 = max( temp_max1, abs(phi_x(i,j,1)) )
            temp_max2 = max( temp_max2, abs(phi_y(i,j,1)) )
        enddo
    enddo
    
    rrrr = 1.
    dt = cfl/( temp_max2/(dx**rrrr) + temp_max1/(dy**rrrr) )
    !dt = cfl /( vx_max/dx + vy_max/dy )

    if(tn+dt >tprint )then
        dt = tprint - tn
    endif

    tn = tn + dt

    nt = nt + 1
    if(nt/2*2==nt) print *,tn,tn/tprint*100,"%"

    end subroutine setdt