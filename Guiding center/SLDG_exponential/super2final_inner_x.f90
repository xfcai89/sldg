    subroutine super2final_inner_x(pes,superpoint1,superpoint2)
    implicit none
    type(type_element_upstream),pointer :: pes
    type(type_inter_point) :: superpoint1,superpoint2    
    integer :: ii,jj
    integer :: nby,iy0
    real :: ymid_t


    iy0 = superpoint1%id_xy
    nby = iy0 - superpoint2%id_xy
 

    point_inner(1) =  superpoint1

    if( nby>0 )then
        do ii = 1,nby
            point_inner(ii+1)%coor(1) = point_inner(1)%coor(1)
            point_inner(ii+1)%coor(2) = ygrid(iy0+1-ii)
            point_inner(ii+1)%ixy_type = 3
            point_inner(ii+1)%idx2 =  point_inner(1)%idx2
            point_inner(ii+1)%idy2 = 2*(iy0+1-ii)-1
        enddo
    elseif( nby<0 )then
        do ii = 1 , abs(nby)
            point_inner(ii+1)%coor(1) = point_inner(1)%coor(1)
            point_inner(ii+1)%coor(2) = ygrid(iy0+ii)
            point_inner(ii+1)%ixy_type = 3
            point_inner(ii+1)%idx2 =  point_inner(1)%idx2
            point_inner(ii+1)%idy2 = 2*(iy0+ii)-1
        enddo
    else

    endif

    point_inner(2+abs(nby) ) = superpoint2

    ! final inner segments
    do ii = 1 , 1 + abs(nby)
        jj = pes%nsub_inner + ii
        pes%segment_inner(jj)%porigin = point_inner(ii)
        pes%segment_inner(jj)%pend = point_inner(ii+1)

        if( pes%segment_inner(jj)%porigin%coor(2)<pes%segment_inner(jj)%pend%coor(2) )then
            ! from bottom to top
            pes%segment_inner(jj)%id(1) = (pes%segment_inner(jj)%porigin%idx2-1)/2
        else
            pes%segment_inner(jj)%id(1) = (pes%segment_inner(jj)%porigin%idx2+1)/2
        endif

        if( pes%segment_inner(jj)%porigin%ixy_type == 3 .and. pes%segment_inner(jj)%pend%ixy_type==3 )then
            ymid_t = (pes%segment_inner(jj)%porigin%coor(2)+pes%segment_inner(jj)%pend%coor(2) )/2.
            pes%segment_inner(jj)%id(2) = ceiling( (ymid_t-ybottom)/dy )
        elseif( pes%segment_inner(jj)%porigin%ixy_type == 1  )then
            pes%segment_inner(jj)%id(2) = pes%segment_inner(jj)%porigin%id_xy
        elseif( pes%segment_inner(jj)%porigin%ixy_type == 3 .and. pes%segment_inner(jj)%pend%ixy_type==1  )then
            pes%segment_inner(jj)%id(2) = pes%segment_inner(jj)%pend%id_xy
        endif



    enddo

    pes%nsub_inner = pes%nsub_inner + abs(nby)+1

    do ii = 1 , 1+abs(nby)
        jj = pes%nsub_inner + ii
        pes%segment_inner(jj)%porigin = point_inner(ii+1)
        pes%segment_inner(jj)%pend = point_inner(ii)

        if( pes%segment_inner(jj)%porigin%coor(2)<pes%segment_inner(jj)%pend%coor(2) )then
            ! from bottom to top
            pes%segment_inner(jj)%id(1) = (pes%segment_inner(jj)%porigin%idx2-1)/2
        else
            pes%segment_inner(jj)%id(1) = (pes%segment_inner(jj)%porigin%idx2+1)/2
        endif

        if( pes%segment_inner(jj)%porigin%ixy_type == 3 .and. pes%segment_inner(jj)%pend%ixy_type==3 )then
            ymid_t = (pes%segment_inner(jj)%porigin%coor(2)+pes%segment_inner(jj)%pend%coor(2) )/2.
            pes%segment_inner(jj)%id(2) = ceiling( (ymid_t-ybottom)/dy )
        elseif( pes%segment_inner(jj)%porigin%ixy_type == 1  )then
            pes%segment_inner(jj)%id(2) = pes%segment_inner(jj)%porigin%id_xy
        elseif( pes%segment_inner(jj)%porigin%ixy_type == 3 .and. pes%segment_inner(jj)%pend%ixy_type==1  )then
            pes%segment_inner(jj)%id(2) = pes%segment_inner(jj)%pend%id_xy
        endif
    enddo

    pes%nsub_inner = pes%nsub_inner + abs(nby)+1


    end subroutine super2final_inner_x