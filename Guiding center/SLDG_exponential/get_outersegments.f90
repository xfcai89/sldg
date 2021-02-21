    subroutine get_outersegments(p)
    implicit none

    type(type_face),pointer :: p
    type(type_segment),pointer :: pouter
    integer :: k
    real :: xmid_t,ymid_t

    do k = 1 , p%nsub_outer
        p%segment_outer(k)%porigin = p%point_inter(k-1)
        p%segment_outer(k)%pend = p%point_inter(k)
        
        pouter => p%segment_outer(k)
        if( pouter%porigin%ixy_type==0 )then
            pouter%id(1) = pouter%porigin%idx2/2
            pouter%id(2) = pouter%porigin%idy2/2
        elseif( pouter%porigin%ixy_type==1 .and. pouter%pend%ixy_type==0 )then
            pouter%id(1) = pouter%pend%idx2/2
            pouter%id(2) = pouter%pend%idy2/2            
        elseif( pouter%porigin%ixy_type==1 .and. pouter%pend%ixy_type==1 )then
            xmid_t = ( pouter%porigin%coor(1)+pouter%pend%coor(1) )/2.
            pouter%id(1) =  ceiling( (xmid_t - xleft)/dx )
            pouter%id(2) = pouter%porigin%idy2/2
        elseif( pouter%porigin%ixy_type==1 .and. pouter%pend%ixy_type==2 )then
            pouter%id(1) = pouter%pend%idx2/2
            pouter%id(2) = pouter%porigin%idy2/2
        elseif( pouter%porigin%ixy_type==2 .and. pouter%pend%ixy_type==0 )then
            pouter%id(1) = pouter%pend%idx2/2
            pouter%id(2) = pouter%pend%idy2/2             
        elseif( pouter%porigin%ixy_type==2 .and. pouter%pend%ixy_type==1 )then
            pouter%id(1) = pouter%porigin%idx2/2
            pouter%id(2) = pouter%pend%idy2/2
        elseif( pouter%porigin%ixy_type==2 .and. pouter%pend%ixy_type==2 )then
            pouter%id(1) = pouter%porigin%idx2/2
            ymid_t = ( pouter%porigin%coor(2)+pouter%pend%coor(2) )/2.
            pouter%id(2) = ceiling( (ymid_t - ybottom)/dy )
        endif
    enddo

    end subroutine get_outersegments