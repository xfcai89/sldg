    subroutine get_innersegments
    implicit none
    type(type_element_upstream),pointer :: pes
    type(type_face),pointer :: pf
    integer :: iy0,nby
    integer :: ix0,nbx
    integer :: ii,jj
    type(type_inter_point) :: superpoint1,superpoint2

    real :: xmid_t,ymid_t

    do i = 1 , nx
        do j = 1 , ny
            !Determine an upstream element's four vertexes and four faces.
            pes => element_star(i,j,io)
            pes%vertex1 => vertex_star(i,j)
            pes%vertex2 => vertex_star(i+1,j)
            pes%vertex3 => vertex_star(i+1,j+1)
            pes%vertex4 => vertex_star(i,j+1)

            pes%edge1_lr => face_lr(i,j)
            pes%edge2_bt => face_bt(i+1,j)
            pes%edge3_lr => face_lr(i,j+1)
            pes%edge4_bt => face_bt(i,j)

            !End determine an upstream element's four vertexes and four faces.

            if( iqc==0 )then
                iam = min(pes%vertex1%id(1),pes%vertex2%id(1),pes%vertex3%id(1),pes%vertex4%id(1) )
                ibm = max(pes%vertex1%id(1),pes%vertex2%id(1),pes%vertex3%id(1),pes%vertex4%id(1) )
                icm = min(pes%vertex1%id(2),pes%vertex2%id(2),pes%vertex3%id(2),pes%vertex4%id(2) )
                idm = max(pes%vertex1%id(2),pes%vertex2%id(2),pes%vertex3%id(2),pes%vertex4%id(2) )
            elseif( iqc==1 )then
                iam = min( pes%edge1_lr%point_origin%id(1),pes%edge1_lr%point_end%id(1) )
                ibm = max( pes%edge1_lr%point_origin%id(1),pes%edge1_lr%point_end%id(1) )
                icm = min( pes%edge1_lr%point_origin%id(2),pes%edge1_lr%point_end%id(2) )
                idm = max( pes%edge1_lr%point_origin%id(2),pes%edge1_lr%point_end%id(2) )
                if( pes%edge1_lr%nsubface==2 )then
                    iam = min( iam,pes%edge1_lr%point_mid1%id(1) )
                    ibm = max( ibm,pes%edge1_lr%point_mid1%id(1) )
                    icm = min( icm,pes%edge1_lr%point_mid1%id(2) )
                    idm = max( idm,pes%edge1_lr%point_mid1%id(2) )
                elseif( pes%edge1_lr%nsubface==3 )then
                    iam = min( iam,pes%edge1_lr%point_mid1%id(1),pes%edge1_lr%point_mid2%id(1) )
                    ibm = max( ibm,pes%edge1_lr%point_mid1%id(1),pes%edge1_lr%point_mid2%id(1) )
                    icm = min( icm,pes%edge1_lr%point_mid1%id(2),pes%edge1_lr%point_mid2%id(2) )
                    idm = max( idm,pes%edge1_lr%point_mid1%id(2),pes%edge1_lr%point_mid2%id(2) )
                endif
                !*******************
                iam = min( iam,pes%edge2_bt%point_origin%id(1),pes%edge2_bt%point_end%id(1) )
                ibm = max( ibm,pes%edge2_bt%point_origin%id(1),pes%edge2_bt%point_end%id(1) )
                icm = min( icm,pes%edge2_bt%point_origin%id(2),pes%edge2_bt%point_end%id(2) )
                idm = max( idm,pes%edge2_bt%point_origin%id(2),pes%edge2_bt%point_end%id(2) )
                if( pes%edge2_bt%nsubface==2 )then
                    iam = min( iam,pes%edge2_bt%point_mid1%id(1) )
                    ibm = max( ibm,pes%edge2_bt%point_mid1%id(1) )
                    icm = min( icm,pes%edge2_bt%point_mid1%id(2) )
                    idm = max( idm,pes%edge2_bt%point_mid1%id(2) )
                elseif( pes%edge2_bt%nsubface==3 )then
                    iam = min( iam,pes%edge2_bt%point_mid1%id(1),pes%edge2_bt%point_mid2%id(1) )
                    ibm = max( ibm,pes%edge2_bt%point_mid1%id(1),pes%edge2_bt%point_mid2%id(1) )
                    icm = min( icm,pes%edge2_bt%point_mid1%id(2),pes%edge2_bt%point_mid2%id(2) )
                    idm = max( idm,pes%edge2_bt%point_mid1%id(2),pes%edge2_bt%point_mid2%id(2) )
                endif
                !*******************
                iam = min( iam,pes%edge3_lr%point_origin%id(1),pes%edge3_lr%point_end%id(1) )
                ibm = max( ibm,pes%edge3_lr%point_origin%id(1),pes%edge3_lr%point_end%id(1) )
                icm = min( icm,pes%edge3_lr%point_origin%id(2),pes%edge3_lr%point_end%id(2) )
                idm = max( idm,pes%edge3_lr%point_origin%id(2),pes%edge3_lr%point_end%id(2) )
                if( pes%edge3_lr%nsubface==2 )then
                    iam = min( iam,pes%edge3_lr%point_mid1%id(1) )
                    ibm = max( ibm,pes%edge3_lr%point_mid1%id(1) )
                    icm = min( icm,pes%edge3_lr%point_mid1%id(2) )
                    idm = max( idm,pes%edge3_lr%point_mid1%id(2) )
                elseif( pes%edge3_lr%nsubface==3 )then
                    iam = min( iam,pes%edge3_lr%point_mid1%id(1),pes%edge3_lr%point_mid2%id(1) )
                    ibm = max( ibm,pes%edge3_lr%point_mid1%id(1),pes%edge3_lr%point_mid2%id(1) )
                    icm = min( icm,pes%edge3_lr%point_mid1%id(2),pes%edge3_lr%point_mid2%id(2) )
                    idm = max( idm,pes%edge3_lr%point_mid1%id(2),pes%edge3_lr%point_mid2%id(2) )
                endif
                !*******************
                iam = min( iam,pes%edge4_bt%point_origin%id(1),pes%edge4_bt%point_end%id(1) )
                ibm = max( ibm,pes%edge4_bt%point_origin%id(1),pes%edge4_bt%point_end%id(1) )
                icm = min( icm,pes%edge4_bt%point_origin%id(2),pes%edge4_bt%point_end%id(2) )
                idm = max( idm,pes%edge4_bt%point_origin%id(2),pes%edge4_bt%point_end%id(2) )
                if( pes%edge4_bt%nsubface==2 )then
                    iam = min( iam,pes%edge4_bt%point_mid1%id(1) )
                    ibm = max( ibm,pes%edge4_bt%point_mid1%id(1) )
                    icm = min( icm,pes%edge4_bt%point_mid1%id(2) )
                    idm = max( idm,pes%edge4_bt%point_mid1%id(2) )
                elseif( pes%edge4_bt%nsubface==3 )then
                    iam = min( iam,pes%edge4_bt%point_mid1%id(1),pes%edge4_bt%point_mid2%id(1) )
                    ibm = max( ibm,pes%edge4_bt%point_mid1%id(1),pes%edge4_bt%point_mid2%id(1) )
                    icm = min( icm,pes%edge4_bt%point_mid1%id(2),pes%edge4_bt%point_mid2%id(2) )
                    idm = max( idm,pes%edge4_bt%point_mid1%id(2),pes%edge4_bt%point_mid2%id(2) )
                endif
            endif

            isx(:) = 0
            isy(:) = 0
            pf =>pes%edge1_lr
            call super_inner(pf )
            pf =>pes%edge2_bt
            call super_inner(pf )
            pf =>pes%edge3_lr
            call super_inner(pf )
            pf =>pes%edge4_bt
            call super_inner(pf )




            pes%nsub_inner = 0
            ! get all inner segments!!
            do ix = 1,ibm-iam
                !**********sorting, if isx(ix)>2
                if( isx(ix)>2 )then
                    do ii = 1,3
                        do jj = 1, 4-ii
                            if( point_inner_x( ix,jj)%coor(2)>point_inner_x( ix,jj+1)%coor(2) )then
                                temp = point_inner_x( ix,jj+1)
                                point_inner_x( ix,jj+1) = point_inner_x( ix,jj)
                                point_inner_x( ix,jj) = temp
                            endif

                        enddo
                    enddo

                endif

                superpoint1 = point_inner_x( ix,1)
                superpoint2 = point_inner_x( ix,2)
                call super2final_inner_x(pes,superpoint1,superpoint2)
                if(isx(ix)>2 )then
                    superpoint1 = point_inner_x( ix,3)
                    superpoint2 = point_inner_x( ix,4)
                    call super2final_inner_x(pes,superpoint1,superpoint2)
                endif
            enddo

            !*********************************************************

            do iy = 1,idm-icm

                if( isy(iy)>2 )then
                    do ii = 1,3
                        do jj = 1, 4-ii
                            if( point_inner_y( iy,jj)%coor(1)>point_inner_y( iy,jj+1)%coor(1) )then
                                temp = point_inner_y( iy,jj+1)
                                point_inner_y( iy,jj+1) = point_inner_y( iy,jj)
                                point_inner_y( iy,jj) = temp
                            endif

                        enddo
                    enddo

                endif

                superpoint1 = point_inner_y( iy,1)
                superpoint2 = point_inner_y( iy,2)
                call super2final_inner_y(pes,superpoint1,superpoint2)
                if( isy(iy)>2 )then
                    superpoint1 = point_inner_y( iy,3)
                    superpoint2 = point_inner_y( iy,4)
                    call super2final_inner_y(pes,superpoint1,superpoint2)
                endif

            enddo

        enddo
    enddo

    end subroutine get_innersegments
    !**************************************************
    subroutine super_inner(pf )
    implicit none
    type(type_face),pointer :: pf

    integer :: k


    do k = 1, pf%nsub_outer-1
        if( pf%point_inter(k)%ixy_type==1 )then
            ix = pf%point_inter(k)%igrid - iam
            if( isx(ix) == 0 )then
                isx(ix) = 1
                point_inner_x( ix,1) = pf%point_inter(k)
            elseif( isx(ix)==1 )then
                isx(ix) = 2
                point_inner_x( ix,2) = pf%point_inter(k)
            elseif( isx(ix)==2 )then
                isx(ix) = 3
                point_inner_x( ix,3) = pf%point_inter(k)
            elseif( isx(ix)==3 )then
                isx(ix) = 4
                point_inner_x( ix,4) = pf%point_inter(k)
            endif
        elseif(pf%point_inter(k)%ixy_type==2)then
            iy = pf%point_inter(k)%igrid - icm
            if( isy(iy) == 0 )then
                isy(iy) = 1
                point_inner_y( iy,1 ) = pf%point_inter(k)
            elseif( isy(iy)==1 )then
                isy(iy) = 2
                point_inner_y( iy,2 ) = pf%point_inter(k)
            elseif( isy(iy)==2 )then
                isy(iy) = 3
                point_inner_y( iy,3 ) = pf%point_inter(k)
            elseif( isy(iy)==3 )then
                isy(iy) = 4
                point_inner_y( iy,4 ) = pf%point_inter(k)
            endif

        endif
    enddo


    end subroutine super_inner