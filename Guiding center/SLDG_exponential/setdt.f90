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
    if(time+dt>time_final) dt = time_final- time

    time = time +dt


    end subroutine setdt