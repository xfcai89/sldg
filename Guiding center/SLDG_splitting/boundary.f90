    subroutine bc_x
    implicit none
    integer :: ibc

    do ibc = 1 , ighost
        elem( 1-ibc,je)%split_modal( 1:n_moment ) = elem( nel_x+1 -ibc,je)%split_modal( 1:n_moment )
        elem( nel_x+ibc,je)%split_modal( 1:n_moment ) = elem(  ibc,je)%split_modal( 1:n_moment )
    enddo

    end subroutine bc_x
    !********************************
    subroutine bc_y
    implicit none
    integer :: ibc

    do ibc = 1 , ighost
        elem( ie, 1-ibc )%split_modal( 1:n_moment ) = elem( ie,nel_y+1 -ibc )%split_modal( 1:n_moment )
        elem( ie,nel_y+ibc )%split_modal( 1:n_moment ) = elem(  ie,ibc )%split_modal( 1:n_moment )
    enddo

    end subroutine bc_y