    subroutine get_intersections_outersegments
    implicit none    

    type(type_face),pointer :: p

    do i = 1 , nx
        do j = 1 , ny + 1
            p=>face_lr(i,j)

            if( iqc==0 )then
                call get_intersections(p)
            elseif( iqc==1 )then
                call get_subfaces(p)
                call get_intersections_qc(p)
            endif
            call get_outersegments(p)
        enddo
    enddo

    do i = 1 , nx +1
        do j = 1 , ny
            p=>face_bt(i,j)
            if( iqc==0 )then
                call get_intersections(p)
            elseif( iqc==1 )then
                call get_subfaces(p)
                call get_intersections_qc(p)
            endif
            call get_outersegments(p)
        enddo
    enddo
    
    end subroutine get_intersections_outersegments