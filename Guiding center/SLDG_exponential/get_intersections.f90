    subroutine get_intersections(p)
    implicit none
    real :: vorigin(1:2),vend(1:2)
    integer :: mx,my
    integer :: ia,ib,ic,id
    integer :: kstart

    type(type_face),pointer :: p

    ia = p%point_origin%id(1)
    ib = p%point_end%id(1)

    ic = p%point_origin%id(2)
    id = p%point_end%id(2)
    vorigin(1:2) = p%point_origin%coor(1:2)
    vend(1:2) = p%point_end%coor(1:2)

    mx = abs(ib-ia)

    my = abs(id-ic)

    kstart = 0
    p%nsub_outer = 1+ mx+my
    p%point_inter( 0)%coor( 1 ) = vorigin(1)
    p%point_inter( 0)%coor( 2 ) = vorigin(2)
    p%point_inter( 0)%ixy_type = 0
    p%point_inter( 0)%idx2 = 2*ia
    p%point_inter( 0)%idy2 = 2*ic
    !******
    p%point_inter( 1+ mx+my)%coor( 1 ) = vend(1)
    p%point_inter( 1+ mx+my)%coor( 2 ) = vend(2)
    p%point_inter( 1+ mx+my)%ixy_type = 0
    p%point_inter( 1+ mx+my)%idx2 = 2*ib
    p%point_inter( 1+ mx+my)%idy2 = 2*id
    !**************************************************************************

    call get_mono_inters(ia,ib,ic,id,mx,my,vorigin(1:2),vend(1:2),p,kstart)


    end subroutine get_intersections
