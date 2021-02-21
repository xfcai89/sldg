    subroutine boundary
    implicit none
    integer :: ib

    do i = 1 , nx
        do ib = 1 , nghost
            element(i,1-ib,io)%umodal(1:n_moment) = element(i,ny+1-ib,io)%umodal(1:n_moment)
            element(i,ny+ib,io)%umodal(1:n_moment) = element(i,0+ib,io)%umodal(1:n_moment)
            !element(i,1-ib,io)%umodal(1:n_moment) = 0.
            !element(i,ny+ib,io)%umodal(1:n_moment) = 0.
        enddo
    enddo

    do j = 1 - nghost , ny + nghost
        do ib = 1 , nghost
            element(1-ib,j,io)%umodal(1:n_moment) = element(nx+1-ib,j,io)%umodal(1:n_moment)
            element(nx+ib,j,io)%umodal(1:n_moment) = element(0+ib,j,io)%umodal(1:n_moment)
        enddo
    enddo


    end subroutine boundary