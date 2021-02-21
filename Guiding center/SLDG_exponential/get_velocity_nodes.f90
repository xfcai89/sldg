    subroutine get_velocity_nodes(idx,idy,io1,v,phix,phiy )
    implicit none
    integer,intent(in) :: idx,idy,io1
    real,intent(in) :: v(1:2)
    real,intent(out) :: phix(0:io1),phiy(0:io1)
    integer :: ii
    real :: tp1,tp2,tp3,tp4
    real :: eps1
    
    eps1 = 1e-7*dx
    eps1 = eps
    
    do ii =0,io1
        phix(ii)=polynomial2d( phi_x_io(idx,idy,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy),dy,n_moment_LDG )
        phiy(ii)=polynomial2d( phi_y_io(idx,idy,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy),dy,n_moment_LDG )
    enddo
    
    
    if( abs( v(1)-(x(idx)-0.5*dx) )<=eps1 )then
    
        do ii = 0,io1
            tp1=polynomial2d( phi_x_io(idx,idy,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy),dy,n_moment_LDG )
            tp2=polynomial2d( phi_x_io(idx-1,idy,1:n_moment_LDG,ii),v(1),x(idx-1),dx,v(2),y(idy),dy,n_moment_LDG )
            phix(ii) = (tp1+tp2)/2.
            tp1=polynomial2d( phi_y_io(idx,idy,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy),dy,n_moment_LDG )
            tp2=polynomial2d( phi_y_io(idx-1,idy,1:n_moment_LDG,ii),v(1),x(idx-1),dx,v(2),y(idy),dy,n_moment_LDG )
            phiy(ii) = (tp1+tp2)/2.
        enddo
        !pause
    elseif( abs( v(1)-(x(idx)+0.5*dx) )<=eps1 )then
    
        do ii = 0,io1
            tp1=polynomial2d( phi_x_io(idx,idy,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy),dy,n_moment_LDG )
            tp2=polynomial2d( phi_x_io(idx+1,idy,1:n_moment_LDG,ii),v(1),x(idx+1),dx,v(2),y(idy),dy,n_moment_LDG )
            phix(ii) = (tp1+tp2)/2.
            tp1=polynomial2d( phi_y_io(idx,idy,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy),dy,n_moment_LDG )
            tp2=polynomial2d( phi_y_io(idx+1,idy,1:n_moment_LDG,ii),v(1),x(idx+1),dx,v(2),y(idy),dy,n_moment_LDG )
            phiy(ii) = (tp1+tp2)/2.
        enddo
        !pause
    elseif( abs( v(2)-(y(idy)-0.5*dy) )<=eps1 )then
    
        do ii = 0,io1
            tp1=polynomial2d( phi_x_io(idx,idy,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy),dy,n_moment_LDG )
            tp2=polynomial2d( phi_x_io(idx,idy-1,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy-1),dy,n_moment_LDG )
            phix(ii) = (tp1+tp2)/2.
            tp1=polynomial2d( phi_y_io(idx,idy,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy),dy,n_moment_LDG )
            tp2=polynomial2d( phi_y_io(idx,idy-1,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy-1),dy,n_moment_LDG )
            phiy(ii) = (tp1+tp2)/2.
        enddo
        !pause
    elseif( abs( v(2)-(y(idy)+0.5*dy) )<=eps1 )then
    
        do ii = 0,io1
            tp1=polynomial2d( phi_x_io(idx,idy,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy),dy,n_moment_LDG )
            tp2=polynomial2d( phi_x_io(idx,idy+1,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy+1),dy,n_moment_LDG )
            phix(ii) = (tp1+tp2)/2.
            tp1=polynomial2d( phi_y_io(idx,idy,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy),dy,n_moment_LDG )
            tp2=polynomial2d( phi_y_io(idx,idy+1,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy+1),dy,n_moment_LDG )
            phiy(ii) = (tp1+tp2)/2.
        enddo
        !pause
    else
        do ii =0,io1
            phix(ii)=polynomial2d( phi_x_io(idx,idy,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy),dy,n_moment_LDG )
            phiy(ii)=polynomial2d( phi_y_io(idx,idy,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy),dy,n_moment_LDG )
        enddo
    endif
    
    if( abs( v(1)-(x(idx)-0.5*dx) )<=eps1 .and. abs( v(2)-(y(idy)-0.5*dy) )<=eps1 )then
        do ii = 0,io1
            tp1=polynomial2d( phi_x_io(idx-1,idy-1,1:n_moment_LDG,ii),v(1),x(idx-1),dx,v(2),y(idy-1),dy,n_moment_LDG )
            tp2=polynomial2d( phi_x_io(idx,idy-1,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy-1),dy,n_moment_LDG )
            tp3=polynomial2d( phi_x_io(idx-1,idy,1:n_moment_LDG,ii),v(1),x(idx-1),dx,v(2),y(idy),dy,n_moment_LDG )
            tp4=polynomial2d( phi_x_io(idx,idy,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy),dy,n_moment_LDG )
            phix(ii) = (tp1+tp2+tp3+tp4)/4.
            tp1=polynomial2d( phi_y_io(idx-1,idy-1,1:n_moment_LDG,ii),v(1),x(idx-1),dx,v(2),y(idy-1),dy,n_moment_LDG )
            tp2=polynomial2d( phi_y_io(idx,idy-1,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy-1),dy,n_moment_LDG )
            tp3=polynomial2d( phi_y_io(idx-1,idy,1:n_moment_LDG,ii),v(1),x(idx-1),dx,v(2),y(idy),dy,n_moment_LDG )
            tp4=polynomial2d( phi_y_io(idx,idy,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy),dy,n_moment_LDG )
            phiy(ii) = (tp1+tp2+tp3+tp4)/4.
        enddo
        !pause
    elseif( abs( v(1)-(x(idx)+0.5*dx) )<=eps1 .and. abs( v(2)-(y(idy)-0.5*dy) )<=eps1 )then
        do ii = 0,io1
            tp1=polynomial2d( phi_x_io(idx,idy-1,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy-1),dy,n_moment_LDG )
            tp2=polynomial2d( phi_x_io(idx+1,idy-1,1:n_moment_LDG,ii),v(1),x(idx+1),dx,v(2),y(idy-1),dy,n_moment_LDG )
            tp3=polynomial2d( phi_x_io(idx,idy,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy),dy,n_moment_LDG )
            tp4=polynomial2d( phi_x_io(idx+1,idy,1:n_moment_LDG,ii),v(1),x(idx+1),dx,v(2),y(idy),dy,n_moment_LDG )
            phix(ii) = (tp1+tp2+tp3+tp4)/4.
            tp1=polynomial2d( phi_y_io(idx,idy-1,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy-1),dy,n_moment_LDG )
            tp2=polynomial2d( phi_y_io(idx+1,idy-1,1:n_moment_LDG,ii),v(1),x(idx+1),dx,v(2),y(idy-1),dy,n_moment_LDG )
            tp3=polynomial2d( phi_y_io(idx,idy,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy),dy,n_moment_LDG )
            tp4=polynomial2d( phi_y_io(idx+1,idy,1:n_moment_LDG,ii),v(1),x(idx+1),dx,v(2),y(idy),dy,n_moment_LDG )
            phiy(ii) = (tp1+tp2+tp3+tp4)/4.
        enddo
        !pause
    elseif( abs( v(1)-(x(idx)-0.5*dx) )<=eps1 .and. abs( v(2)-(y(idy)+0.5*dy) )<=eps1 )then
        do ii = 0,io1
            tp1=polynomial2d( phi_x_io(idx-1,idy,1:n_moment_LDG,ii),v(1),x(idx-1),dx,v(2),y(idy),dy,n_moment_LDG )
            tp2=polynomial2d( phi_x_io(idx,idy,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy),dy,n_moment_LDG )
            tp3=polynomial2d( phi_x_io(idx-1,idy+1,1:n_moment_LDG,ii),v(1),x(idx-1),dx,v(2),y(idy+1),dy,n_moment_LDG )
            tp4=polynomial2d( phi_x_io(idx,idy+1,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy+1),dy,n_moment_LDG )
            phix(ii) = (tp1+tp2+tp3+tp4)/4.
            tp1=polynomial2d( phi_y_io(idx-1,idy,1:n_moment_LDG,ii),v(1),x(idx-1),dx,v(2),y(idy),dy,n_moment_LDG )
            tp2=polynomial2d( phi_y_io(idx,idy,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy),dy,n_moment_LDG )
            tp3=polynomial2d( phi_y_io(idx-1,idy+1,1:n_moment_LDG,ii),v(1),x(idx-1),dx,v(2),y(idy+1),dy,n_moment_LDG )
            tp4=polynomial2d( phi_y_io(idx,idy+1,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy+1),dy,n_moment_LDG )
            phiy(ii) = (tp1+tp2+tp3+tp4)/4.
        enddo
        !pause
    elseif( abs( v(1)-(x(idx)+0.5*dx) )<=eps1 .and. abs( v(2)-(y(idy)+0.5*dy) )<=eps1 )then
        do ii = 0,io1
            tp1=polynomial2d( phi_x_io(idx,idy,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy),dy,n_moment_LDG )
            tp2=polynomial2d( phi_x_io(idx+1,idy,1:n_moment_LDG,ii),v(1),x(idx+1),dx,v(2),y(idy),dy,n_moment_LDG )
            tp3=polynomial2d( phi_x_io(idx,idy+1,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy+1),dy,n_moment_LDG )
            tp4=polynomial2d( phi_x_io(idx+1,idy+1,1:n_moment_LDG,ii),v(1),x(idx+1),dx,v(2),y(idy+1),dy,n_moment_LDG )
            phix(ii) = (tp1+tp2+tp3+tp4)/4.
            tp1=polynomial2d( phi_y_io(idx,idy,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy),dy,n_moment_LDG )
            tp2=polynomial2d( phi_y_io(idx+1,idy,1:n_moment_LDG,ii),v(1),x(idx+1),dx,v(2),y(idy),dy,n_moment_LDG )
            tp3=polynomial2d( phi_y_io(idx,idy+1,1:n_moment_LDG,ii),v(1),x(idx),dx,v(2),y(idy+1),dy,n_moment_LDG )
            tp4=polynomial2d( phi_y_io(idx+1,idy+1,1:n_moment_LDG,ii),v(1),x(idx+1),dx,v(2),y(idy+1),dy,n_moment_LDG )
            phiy(ii) = (tp1+tp2+tp3+tp4)/4.
        enddo
        !pause
    endif

    end subroutine get_velocity_nodes