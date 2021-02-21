    subroutine green(sum16)
    implicit none

    real,intent(out) :: sum16(1:n_moment)
    type(type_element_upstream),pointer :: pes
    type(type_element),pointer :: pe
    type(type_face),pointer :: pf
    real :: anticlock
    integer :: nm
    real :: sum
    !*********** least square ***********
    real :: x_base,y_base
    real :: aa(1:n_moment)
    real :: vert_temp(1:9,1:5)
    real :: A_temp(6,6)
    real :: temp_L(6,6),temp_U(6,6)
    real :: psi(1:9),b(1:6)


    pes => element_star(i,j,io)


    if(iadditive(io)==0 )then
        pe => element(i,j,0)
    elseif( iadditive(io)==1 )then
        pe => element(i,j,io)
    endif

    if(io==2 .and. iexp_case==7 )then
        pe => element(i,j,1)
    endif

    x_base = 0.25*( pes%vertex1%coor(1) + pes%vertex2%coor(1)  &
        + pes%vertex3%coor(1)  + pes%vertex4%coor(1)  )
    y_base = 0.25*( pes%vertex1%coor(2) + pes%vertex2%coor(2)  &
        + pes%vertex3%coor(2)  + pes%vertex4%coor(2)  )
    !x_base = 1./9.*( pes%vertex1%coor(1) + pes%vertex2%coor(1)  &
    !    + pes%vertex3%coor(1)  + pes%vertex4%coor(1) &
    !    +pes%vertex5%coor(1) + pes%vertex6%coor(1)  &
    !    + pes%vertex7%coor(1)  + pes%vertex8%coor(1) + pes%vertex9%coor(1) )
    !y_base = 1./9.*( pes%vertex1%coor(2) + pes%vertex2%coor(2)  &
    !    + pes%vertex3%coor(2)  + pes%vertex4%coor(2) &
    !    + pes%vertex5%coor(2) + pes%vertex6%coor(2)  &
    !    + pes%vertex7%coor(2)  + pes%vertex8%coor(2) + pes%vertex9%coor(2) )

    do nm = 1 ,n_moment
        sum = 0.
        !
        
        call interpolate_gauss(nm,pe,pes,temp_L,temp_U,aa)
        !if(nm == 1)then
        !    aa(1) = 1.
        !    aa(2:n_moment) = 0.
        !else
        !    if( nm ==2 )then
        !
        !        if(n_moment<=3)then
        !            vert_temp(1,1) = ( pes%vertex1%coor(1) - x_base )/dx
        !            vert_temp(2,1) = ( pes%vertex2%coor(1) - x_base )/dx
        !            vert_temp(3,1) = ( pes%vertex3%coor(1) - x_base )/dx
        !            vert_temp(4,1) = ( pes%vertex4%coor(1) - x_base )/dx
        !
        !            vert_temp(1,2) = ( pes%vertex1%coor(2) - y_base )/dy
        !            vert_temp(2,2) = ( pes%vertex2%coor(2) - y_base )/dy
        !            vert_temp(3,2) = ( pes%vertex3%coor(2) - y_base )/dy
        !            vert_temp(4,2) = ( pes%vertex4%coor(2) - y_base )/dy
        !            call get_matrix_a( vert_temp(1:4,1:2),  A_temp(1:3,1:3) )
        !
        !            call doolittle(A_temp(1:3,1:3),temp_L(1:3,1:3),temp_U(1:3,1:3),3)
        !        elseif(n_moment==6)then
        !            vert_temp(1,1) = ( pes%vertex1%coor(1) - x_base )/dx
        !            vert_temp(2,1) = ( pes%vertex2%coor(1) - x_base )/dx
        !            vert_temp(3,1) = ( pes%vertex3%coor(1) - x_base )/dx
        !            vert_temp(4,1) = ( pes%vertex4%coor(1) - x_base )/dx
        !            vert_temp(5,1) = ( pes%vertex5%coor(1) - x_base )/dx
        !            vert_temp(6,1) = ( pes%vertex6%coor(1) - x_base )/dx
        !            vert_temp(7,1) = ( pes%vertex7%coor(1) - x_base )/dx
        !            vert_temp(8,1) = ( pes%vertex8%coor(1) - x_base )/dx
        !            vert_temp(9,1) = ( pes%vertex9%coor(1) - x_base )/dx
        !
        !            vert_temp(1,2) = ( pes%vertex1%coor(2) - y_base )/dy
        !            vert_temp(2,2) = ( pes%vertex2%coor(2) - y_base )/dy
        !            vert_temp(3,2) = ( pes%vertex3%coor(2) - y_base )/dy
        !            vert_temp(4,2) = ( pes%vertex4%coor(2) - y_base )/dy
        !            vert_temp(5,2) = ( pes%vertex5%coor(2) - y_base )/dy
        !            vert_temp(6,2) = ( pes%vertex6%coor(2) - y_base )/dy
        !            vert_temp(7,2) = ( pes%vertex7%coor(2) - y_base )/dy
        !            vert_temp(8,2) = ( pes%vertex8%coor(2) - y_base )/dy
        !            vert_temp(9,2) = ( pes%vertex9%coor(2) - y_base )/dy
        !
        !            vert_temp(1:9,3) = vert_temp(1:9,1)*vert_temp(1:9,1) - 1./12.
        !            vert_temp(1:9,4) = vert_temp(1:9,1)*vert_temp(1:9,2)
        !            vert_temp(1:9,5) = vert_temp(1:9,2)*vert_temp(1:9,2) - 1./12.
        !
        !            call get_matrix2_a( vert_temp(1:9,1:5),  A_temp(1:6,1:6) )
        !            call doolittle(A_temp(1:6,1:6),temp_L(1:6,1:6),temp_U(1:6,1:6),6)
        !        endif
        !    endif
        !
        !    if(n_moment<=3)then
        !        psi(1) = fphi( nm, (pe%vertex1%coor(1)-x(i) )/dx,(pe%vertex1%coor(2)-y(j))/dy ) ;
        !        psi(2) = fphi( nm, (pe%vertex2%coor(1)-x(i) )/dx,(pe%vertex2%coor(2)-y(j))/dy ) ;
        !        psi(3) = fphi( nm, (pe%vertex3%coor(1)-x(i) )/dx,(pe%vertex3%coor(2)-y(j))/dy ) ;
        !        psi(4) = fphi( nm, (pe%vertex4%coor(1)-x(i) )/dx,(pe%vertex4%coor(2)-y(j))/dy ) ;
        !
        !        call get_vector_b( vert_temp(1:4,1:2),psi(1:4),b(1:3) )
        !        call solve17(temp_L(1:3,1:3),temp_U(1:3,1:3),b(1:3),aa(1:3),3)
        !    elseif(n_moment==6)then
        !        psi(1) = fphi( nm, (pe%vertex1%coor(1)-x(i) )/dx,(pe%vertex1%coor(2)-y(j))/dy ) ;
        !        psi(2) = fphi( nm, (pe%vertex2%coor(1)-x(i) )/dx,(pe%vertex2%coor(2)-y(j))/dy ) ;
        !        psi(3) = fphi( nm, (pe%vertex3%coor(1)-x(i) )/dx,(pe%vertex3%coor(2)-y(j))/dy ) ;
        !        psi(4) = fphi( nm, (pe%vertex4%coor(1)-x(i) )/dx,(pe%vertex4%coor(2)-y(j))/dy ) ;
        !        psi(5) = fphi( nm, (pe%vertex5%coor(1)-x(i) )/dx,(pe%vertex5%coor(2)-y(j))/dy ) ;
        !        psi(6) = fphi( nm, (pe%vertex6%coor(1)-x(i) )/dx,(pe%vertex6%coor(2)-y(j))/dy ) ;
        !        psi(7) = fphi( nm, (pe%vertex7%coor(1)-x(i) )/dx,(pe%vertex7%coor(2)-y(j))/dy ) ;
        !        psi(8) = fphi( nm, (pe%vertex8%coor(1)-x(i) )/dx,(pe%vertex8%coor(2)-y(j))/dy ) ;
        !        psi(9) = fphi( nm, (pe%vertex9%coor(1)-x(i) )/dx,(pe%vertex9%coor(2)-y(j))/dy ) ;
        !
        !        call get_vector2_b( vert_temp(1:9,1:5),psi(1:9),b(1:6) )
        !        call solve17(temp_L(1:6,1:6),temp_U(1:6,1:6),b(1:6),aa(1:6),6)
        !    endif
        !
        !
        !    !if(time>0.8)write(99,*) i,j,nm,aa(1:6)
        !
        !endif
        !******************************************************************************************

        if( nm == 1 )then
            pf => pes%edge1_lr
            call green_outersegments_face( pf,x_base,y_base  )
            pf => pes%edge2_bt
            call green_outersegments_face( pf,x_base,y_base  )
            pf => pes%edge3_lr
            call green_outersegments_face( pf,x_base,y_base  )
            pf => pes%edge4_bt
            call green_outersegments_face( pf,x_base,y_base  )

            call green_innersegments(pes,x_base,y_base)
        endif

        pf => pes%edge1_lr
        anticlock = 1.
        call sum_outersegments(nm, pf,anticlock,aa(1:n_moment),sum )
        pf => pes%edge2_bt
        anticlock = 1.
        call sum_outersegments(nm, pf,anticlock,aa(1:n_moment),sum )
        pf => pes%edge3_lr
        anticlock = -1.
        call sum_outersegments(nm, pf,anticlock,aa(1:n_moment),sum )
        pf => pes%edge4_bt
        anticlock = -1.
        call sum_outersegments(nm, pf,anticlock,aa(1:n_moment),sum )


        call sum_innersegments(nm, pes,aa(1:n_moment),sum )


        sum16(nm) = sum * ai(nm)

    enddo

    end subroutine green
    !**************************************
    subroutine sum_innersegments( nm, pes,aa,sum )
    implicit none
    integer,intent(in) :: nm
    type(type_element_upstream),pointer :: pes
    real,intent(inout) :: sum
    real,intent(in) :: aa(1:n_moment)
    integer :: ii,kk

    integer :: idx,idy,idx1,idy1
    real :: sum_temp

    do kk = 1 , pes%nsub_inner
        do ii = 1 ,n_moment
            sum = sum + aa(ii)*pes%segment_inner(kk)%c_ab(ii)
        enddo
        !*****************************************test mass
        !*****************************************important part
        if(nm==1)then
            idx = pes%segment_inner(kk)%id(1)
            idy = pes%segment_inner(kk)%id(2)
            if( idx > nx )then
                idx1 = idx -nx
            elseif(idx<1)then
                idx1 = nx+idx
            else
                idx1 = idx
            endif

            if( idy > ny )then
                idy1 = idy -ny
            elseif(idy<1)then
                idy1 = ny+idy
            else
                idy1 = idy
            endif

            sum_temp = 0.
            do ii = 1 ,n_moment
                sum_temp = sum_temp + aa(ii)*pes%segment_inner(kk)%c_ab(ii)
            enddo

            com_mass( idx1,idy1 ) = com_mass(idx1,idy1) +  sum_temp
            
            !if( time>0.2 .and. idx1==122 .and. idy1==56 .and. io==2 )then
            !    print *,pes%segment_inner(kk)%porigin%coor(1:2), pes%segment_inner(kk)%pend%coor(1:2)
            !    print *,i,j
            !    pause
            !endif
            
        endif
        !*****************************************test mass
        !*****************************************important part
        
        
        
    enddo

    end subroutine sum_innersegments
    !**************************************
    subroutine sum_outersegments( nm,pf,anticlock,aa,sum )
    implicit none
    integer,intent(in) :: nm
    type(type_face),pointer :: pf
    real,intent(inout) :: sum
    real,intent(in) :: aa(1:n_moment)
    real,intent(in) :: anticlock
    integer :: ii,kk

    integer :: idx,idy,idx1,idy1
    real :: sum_temp

    do kk = 1 , pf%nsub_outer
        do ii = 1 ,n_moment
            sum = sum + anticlock * aa(ii)*pf%segment_outer(kk)%c_ab(ii)
        enddo

        !*****************************************test mass
        !*****************************************important part
        if(nm==1)then
            idx = pf%segment_outer(kk)%id(1)
            idy = pf%segment_outer(kk)%id(2)
            if( idx > nx )then
                idx1 = idx -nx
            elseif(idx<1)then
                idx1 = nx+idx
            else
                idx1 = idx
            endif

            if( idy > ny )then
                idy1 = idy -ny
            elseif(idy<1)then
                idy1 = ny+idy
            else
                idy1 = idy
            endif

            sum_temp = 0.
            do ii = 1 ,n_moment
                sum_temp = sum_temp + anticlock* aa(ii)*pf%segment_outer(kk)%c_ab(ii)
            enddo

            com_mass( idx1,idy1 ) = com_mass(idx1,idy1) +  sum_temp
            
            !if( time>0.2 .and. idx1==122 .and. idy1==56 .and. io==2 )then
            !    print *,pf%segment_outer(kk)%porigin%coor(1:2),pf%segment_outer(kk)%pend%coor(1:2)
            !    print *,i,j,anticlock
            !    pause
            !endif
            
        endif
        !*****************************************test mass
        !*****************************************important part
    enddo

    end subroutine sum_outersegments
    !**************************************
    subroutine green_outersegments_face(pf,x_base,y_base)
    implicit none
    type(type_face),pointer :: pf
    integer :: k
    real :: v1(1:2),v2(1:2)
    real,intent(in) :: x_base,y_base
    integer :: idx,idy
    real :: s_ab(1:n_moment)
    real :: sol(1:n_moment)
    real :: dl(1:2),de(1:2),centrl(1:2)

    real :: v11(2),v22(2),v33(2)

    ! for qc
    v11(1:2) = pf%point_origin%coor(1:2)
    v22(1:2) = pf%point_midt%coor(1:2)
    v33(1:2) = pf%point_end%coor(1:2)
    ! end for qc

    do k = 1 , pf%nsub_outer
        v1(1:2) = pf%segment_outer(k)%porigin%coor(1:2)
        v2(1:2) = pf%segment_outer(k)%pend%coor(1:2)

        idx = pf%segment_outer(k)%id(1)
        idy = pf%segment_outer(k)%id(2)
        !call green_p1_gauss2_line_integral(v1,v2,x_base,y_base,idx,idy,s_ab)
        !call green_p2_gauss3_line_integral(v1,v2,x_base,y_base,idx,idy,s_ab)

        if(iadditive(io)==0 )then
            sol(1:n_moment) = element(idx,idy,0)%umodal(1:n_moment)
        elseif( iadditive(io)==1 )then
            sol(1:n_moment) = element(idx,idy,io)%umodal(1:n_moment)
        endif

        if(io==2 .and. iexp_case==7 )then
            sol(1:n_moment) = element(idx,idy,1)%umodal(1:n_moment)
        endif

        if( iqc==0 )then
            dl(1) = dx
            dl(2) = dy
            de(1) = dx
            de(2) = dy
            centrl(1) = x(idx)
            centrl(2) = y(idy)
            call segment_int(n_moment,v1(1:2),v2(1:2),x_base,y_base,idx,idy,centrl,dl,de,sol(1:n_moment),s_ab(1:n_moment))
        elseif( iqc==1 )then
            call green_p2qc_gauss3_outer( v1,v2,x_base,y_base,idx,idy,s_ab,v11,v22,v33,sol(1:n_moment) )
        endif
        pf%segment_outer(k)%c_ab(1:n_moment) = s_ab(1:n_moment)
    enddo

    end subroutine green_outersegments_face
    !************************************************
    subroutine green_innersegments(pes,x_base,y_base)
    implicit none
    type(type_element_upstream),pointer :: pes
    real :: x_base,y_base
    integer :: k
    integer :: idx,idy
    real :: s_ab(1:n_moment)
    real :: v1(1:2),v2(1:2)
    real :: sol(1:n_moment),dl(1:2),de(1:2),centrl(1:2)

    do k = 1 , pes%nsub_inner
        v1(1:2) = pes%segment_inner(k)%porigin%coor(1:2)
        v2(1:2) = pes%segment_inner(k)%pend%coor(1:2)

        idx = pes%segment_inner(k)%id(1)
        idy = pes%segment_inner(k)%id(2)

        !idebug = 1
        !call green_p1_gauss2_line_integral(v1,v2,x_base,y_base,idx,idy,s_ab)
        !call green_p2_gauss3_line_integral(v1,v2,x_base,y_base,idx,idy,s_ab)
        if(iadditive(io)==0 )then
            sol(1:n_moment) = element(idx,idy,0)%umodal(1:n_moment)
        elseif( iadditive(io)==1 )then
            sol(1:n_moment) = element(idx,idy,io)%umodal(1:n_moment)
        endif

        if(io==2 .and. iexp_case==7 )then
            sol(1:n_moment) = element(idx,idy,1)%umodal(1:n_moment)
        endif


        dl(1) = dx
        dl(2) = dy
        de(1) = dx
        de(2) = dy
        centrl(1) = x(idx)
        centrl(2) = y(idy)
        call segment_int(n_moment,v1(1:2),v2(1:2),x_base,y_base,idx,idy,centrl,dl,de,sol(1:n_moment),s_ab(1:n_moment))
        pes%segment_inner(k)%c_ab(1:n_moment) = s_ab(1:n_moment)
    enddo

    end subroutine green_innersegments