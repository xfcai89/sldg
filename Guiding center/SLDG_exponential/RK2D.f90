    subroutine RK2D(vx,vx_s,t0,dt,fb )
    implicit none
    real,intent(in) :: vx(1:2),t0,dt,fb
    real,intent(out) :: vx_s(1:2)
    real :: v1(1:2),v2(1:2)
    real :: v3(1:2)
    real :: rk1(1:2),rk2(1:2),rk3(1:2),rk4(1:2),rk5(1:2),rk6(1:2)

    if(irk ==3)then
        v1(1) = vx(1) +fb* vel_x( vx(1), vx(2) ,t0 ) * dt
        v1(2) = vx(2) +fb* vel_y( vx(1), vx(2) ,t0 ) * dt
        v2(1) = 0.75*vx(1) + 0.25*(  v1(1) +fb* vel_x( v1(1),v1(2),t0+ fb*dt ) * dt )
        v2(2) = 0.75*vx(2) + 0.25*(  v1(2) +fb* vel_y( v1(1),v1(2),t0+ fb*dt ) * dt )
        vx_s(1) = vx(1)/3. + 2./3.*( v2(1) +fb* vel_x( v2(1),v2(2),t0+fb*dt*0.5 ) * dt  )
        vx_s(2) = vx(2)/3. + 2./3.*( v2(2) +fb* vel_y( v2(1),v2(2),t0+fb*dt*0.5 ) * dt  )

    elseif(irk ==4)then

        v1(1) = vx(1) +fb*0.5* vel_x( vx(1),vx(2),t0 ) * dt
        v1(2) = vx(2) +fb*0.5* vel_y( vx(1),vx(2),t0 ) * dt
        v2(1) = vx(1) +fb*0.5* vel_x( v1(1),v1(2),t0+fb*dt*0.5 ) * dt
        v2(2) = vx(2) +fb*0.5* vel_y( v1(1),v1(2),t0+fb*dt*0.5 ) * dt
        v3(1) = vx(1) +fb*vel_x( v2(1),v2(2),t0+fb*dt*0.5 ) *dt
        v3(2) = vx(2) +fb*vel_y( v2(1),v2(2),t0+fb*dt*0.5 ) *dt
        vx_s(1) = 1./3.*( -vx(1) + v1(1) +2.*v2(1)+v3(1) ) +fb* 1./6.*dt* vel_x( v3(1),v3(2),t0+fb*dt  )
        vx_s(2) = 1./3.*( -vx(2) + v1(2) +2.*v2(2)+v3(2) ) +fb* 1./6.*dt* vel_y( v3(1),v3(2),t0+fb*dt  )

    elseif(irk==5)then

        rk1(1) = vel_x( vx(1),vx(2),t0 )
        rk1(2) = vel_y( vx(1),vx(2),t0 )

        rk2(1) = vel_x( vx(1)+fb*0.25*rk1(1)*dt , vx(2)+fb*0.25*rk1(2)*dt   ,t0+fb*0.25*dt)
        rk2(2) = vel_y( vx(1)+fb*0.25*rk1(1)*dt , vx(2)+fb*0.25*rk1(2)*dt   ,t0+fb*0.25*dt)

        rk3(1) = vel_x( vx(1)+fb*rk1(1)*dt/8.+fb*rk2(1)*dt/8. , &
            vx(2)+fb*rk1(2)*dt/8.+fb*rk2(2)*dt/8. ,t0+fb*0.25*dt)
        rk3(2) = vel_y( vx(1)+fb*rk1(1)*dt/8.+fb*rk2(1)*dt/8. , &
            vx(2)+fb*rk1(2)*dt/8.+fb*rk2(2)*dt/8. ,t0+fb*0.25*dt)

        rk4(1) = vel_x(  vx(1)-fb*rk2(1)*dt/2.+fb*rk3(1)*dt , &
            vx(2)-fb*rk2(2)*dt/2.+fb*rk3(2)*dt ,t0+fb*0.5*dt)
        rk4(2) = vel_y( vx(1)-fb*rk2(1)*dt/2.+fb*rk3(1)*dt , &
            vx(2)-fb*rk2(2)*dt/2.+fb*rk3(2)*dt ,t0+fb*0.5*dt)

        rk5(1) = vel_x(  vx(1)+fb*rk1(1)*dt*3./16.+fb*rk4(1)*dt*9./16. , &
            vx(2)+fb*rk1(2)*dt*3./16.+fb*rk4(2)*dt*9./16. ,t0+fb*0.75*dt)
        rk5(2) = vel_y( vx(1)+fb*rk1(1)*dt*3./16.+fb*rk4(1)*dt*9./16. , &
            vx(2)+fb*rk1(2)*dt*3./16.+fb*rk4(2)*dt*9./16. ,t0+fb*0.75*dt)

        rk6(1) = vel_x( vx(1)-fb*rk1(1)*dt*3./7.+fb*rk2(1)*dt*2./7.+fb*rk3(1)*dt*12./7.-fb*rk4(1)*dt*12./7. &
            +fb*rk5(1)*dt*8./7. , &
            vx(2)-fb*rk1(2)*dt*3./7.+fb*rk2(2)*dt*2./7.+fb*rk3(2)*dt*12./7.-fb*rk4(2)*dt*12./7. &
            +fb*rk5(2)*dt*8./7. ,t0+fb*dt)
        rk6(2) = vel_y(  vx(1)-fb*rk1(1)*dt*3./7.+fb*rk2(1)*dt*2./7.+fb*rk3(1)*dt*12./7.-fb*rk4(1)*dt*12./7. &
            +fb*rk5(1)*dt*8./7. , &
            vx(2)-fb*rk1(2)*dt*3./7.+fb*rk2(2)*dt*2./7.+fb*rk3(2)*dt*12./7.-fb*rk4(2)*dt*12./7. &
            +fb*rk5(2)*dt*8./7. ,t0+fb*dt)

        vx_s(1) = vx(1) +fb* (7.*rk1(1)+32.*rk3(1)+12.*rk4(1)+32.*rk5(1)+7.*rk6(1) )*dt/90.
        vx_s(2) = vx(2) +fb* (7.*rk1(2)+32.*rk3(2)+12.*rk4(2)+32.*rk5(2)+7.*rk6(2) )*dt/90.

    endif

    end subroutine RK2D