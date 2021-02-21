    subroutine super2final_inner_y(pes,superpoint1,superpoint2)
    implicit none
    type(type_element_upstream),pointer :: pes
    type(type_inter_point) :: superpoint1,superpoint2
    integer :: ix0,nbx
    integer :: ii,jj
    real :: xmid_t
    
    ix0 = superpoint1%id_xy
    nbx = ix0 - superpoint2%id_xy
 
    point_inner(1) = superpoint1

    if( nbx>0 )then
        do ii = 1 , nbx
            point_inner(ii+1)%coor(1) = xgrid(ix0+1-ii)
            point_inner(ii+1)%coor(2) = point_inner(1)%coor(2)
            point_inner(ii+1)%ixy_type = 3
            point_inner(ii+1)%idx2 =  2*(ix0+1-ii)-1
            point_inner(ii+1)%idy2 = point_inner(1)%idy2

        enddo
    elseif( nbx<0 )then
        do ii = 1 , abs(nbx)
            point_inner(ii+1)%coor(1) = xgrid( ix0+ii )
            point_inner(ii+1)%coor(2) = point_inner(1)%coor(2)
            point_inner(ii+1)%ixy_type = 3
            point_inner(ii+1)%idx2 =  2*(ix0+ii)-1
            point_inner(ii+1)%idy2 = point_inner(1)%idy2
        enddo
    else

    endif

    point_inner(2+abs(nbx) ) = superpoint2

    do ii = 1 , 1 + abs(nbx)
        jj = pes%nsub_inner + ii
        pes%segment_inner(jj)%porigin = point_inner(ii)
        pes%segment_inner(jj)%pend = point_inner(ii+1)

        if( pes%segment_inner(jj)%porigin%coor(1)<pes%segment_inner(jj)%pend%coor(1) )then
            ! from left to right
            pes%segment_inner(jj)%id(2)=(pes%segment_inner(jj)%porigin%idy2+1)/2
        else
            pes%segment_inner(jj)%id(2)=(pes%segment_inner(jj)%porigin%idy2-1)/2
        endif

        if(pes%segment_inner(jj)%porigin%ixy_type==3 .and. pes%segment_inner(jj)%pend%ixy_type==3 )then
            xmid_t = (pes%segment_inner(jj)%porigin%coor(1)+pes%segment_inner(jj)%pend%coor(1) )/2.
            pes%segment_inner(jj)%id(1) = ceiling( (xmid_t-xleft)/dx )
        elseif( pes%segment_inner(jj)%porigin%ixy_type==2 )then
            pes%segment_inner(jj)%id(1) = pes%segment_inner(jj)%porigin%id_xy
        elseif( pes%segment_inner(jj)%porigin%ixy_type==3 .and. pes%segment_inner(jj)%pend%ixy_type==2 )then
            pes%segment_inner(jj)%id(1) = pes%segment_inner(jj)%pend%id_xy
        endif

    enddo

    pes%nsub_inner = pes%nsub_inner + abs(nbx) + 1

    do ii = 1 , 1 + abs(nbx)
        jj = pes%nsub_inner + ii
        
        pes%segment_inner(jj)%porigin = point_inner(ii+1)
        pes%segment_inner(jj)%pend = point_inner(ii)


        if( pes%segment_inner(jj)%porigin%coor(1)<pes%segment_inner(jj)%pend%coor(1) )then
            ! from left to right
            pes%segment_inner(jj)%id(2)=(pes%segment_inner(jj)%porigin%idy2+1)/2
        else
            pes%segment_inner(jj)%id(2)=(pes%segment_inner(jj)%porigin%idy2-1)/2
        endif

        if(pes%segment_inner(jj)%porigin%ixy_type==3 .and. pes%segment_inner(jj)%pend%ixy_type==3 )then
            xmid_t = (pes%segment_inner(jj)%porigin%coor(1)+pes%segment_inner(jj)%pend%coor(1) )/2.
            pes%segment_inner(jj)%id(1) = ceiling( (xmid_t-xleft)/dx )
        elseif( pes%segment_inner(jj)%porigin%ixy_type==2 )then
            pes%segment_inner(jj)%id(1) = pes%segment_inner(jj)%porigin%id_xy
        elseif( pes%segment_inner(jj)%porigin%ixy_type==3 .and. pes%segment_inner(jj)%pend%ixy_type==2 )then
            pes%segment_inner(jj)%id(1) = pes%segment_inner(jj)%pend%id_xy
        endif

    enddo
    pes%nsub_inner = pes%nsub_inner + abs(nbx) + 1

    end subroutine super2final_inner_y