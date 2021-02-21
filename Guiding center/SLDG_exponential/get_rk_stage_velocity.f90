    subroutine get_rk_stage_velocity(io1,phix,phiy,vel_x,vel_y)
    implicit none
    integer,intent(in) :: io1
    real,intent(in) :: phix(0:io1),phiy(0:io1)
    real,intent(out) :: vel_x,vel_y
    integer :: ii

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

    end subroutine get_rk_stage_velocity