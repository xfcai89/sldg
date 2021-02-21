    subroutine get_mono_inters(ia,ib,ic,id,mx,my,vorigin,vend,p,kstart)
    implicit none
    integer,intent(in) :: ia,ib,ic,id,mx,my
    real,intent(in) :: vorigin(1:2),vend(1:2)
    integer,intent(in) :: kstart
    integer :: ii,jj
    integer :: k,inter,kk,ki
    integer :: ia2,ib2,ic2,id2
    type(type_face),pointer :: p
    
     
    !**************************************************************************
    ! Algorithm of getting the intersection points between face and gridlines
    ! find intersections between face and x=x_i
    
    if(mx .ne. 0)then

        do k = 1, mx
            if(ia<ib) then
                inter = ia + k
            else
                inter = ia + 1 - k
            endif
            
            kk = kstart + k
            p%point_inter( kk)%coor( 1 ) = xgrid( inter )
            p%point_inter( kk)%coor( 2 ) = yside( vorigin(1:2), vend(1:2),  xgrid(inter) )
            p%point_inter( kk)%ixy_type = 1
            p%point_inter( kk)%igrid = inter
            p%point_inter( kk)%idx2 = 2*inter-1
            p%point_inter( kk)%id_xy = ceiling( (p%point_inter( kk)%coor( 2 )-ybottom)/dy )
            p%point_inter( kk)%idy2 = 2*p%point_inter( kk)%id_xy
        enddo
    endif
    !find intersections between face and y=y_j
    
    if(my .ne. 0)then
 
        do k = 1 + mx, mx+my
            ki = k - mx
            if(ic<id)then
                inter = ic + ki
            else
                inter = ic  + 1 - ki
            endif
            
            kk = kstart + k
            p%point_inter( kk)%coor( 2 ) = ygrid( inter )
            p%point_inter( kk)%coor( 1 ) = xside( vorigin(1:2), vend(1:2),  ygrid(inter) )
            p%point_inter( kk)%ixy_type = 2
            p%point_inter( kk)%igrid = inter
            p%point_inter( kk)%idy2 = 2*inter-1
            p%point_inter( kk)%id_xy = ceiling( (p%point_inter( kk)%coor( 1 )-xleft)/dx )
            p%point_inter( kk)%idx2 = 2*p%point_inter( kk)%id_xy

        enddo
    endif
    !**************************************************************************
    ! Algorithm of sorting point_inter(1:mx+my) based on coordinates
    if(mx+my>1)then
        if(mx>=my)then
            if( vorigin(1)<vend(1) )then
                ! sort them in a increasing order
                do ii = 1 ,mx+my-1
                    do jj = kstart + 1, kstart+mx+my-ii
                        if( p%point_inter(jj)%coor(1)>p%point_inter(jj+1)%coor(1) )then
                            temp = p%point_inter(jj+1)
                            p%point_inter(jj+1) = p%point_inter(jj)
                            p%point_inter(jj) = temp
                        endif
                    enddo
                enddo
            else
                do ii = 1,mx+my-1
                    do jj = kstart+1 , kstart+mx+my-ii
                        if( p%point_inter(jj)%coor(1)<p%point_inter(jj+1)%coor(1) )then
                            temp = p%point_inter(jj+1)
                            p%point_inter(jj+1) = p%point_inter(jj)
                            p%point_inter(jj) = temp
                        endif
                    enddo
                enddo
            endif
        elseif(mx<my)then !mx>=my
            if( vorigin(2)<vend(2) )then
                do ii = 1, mx+my-1
                    do jj = kstart+1,kstart+mx+my-ii
                        if( p%point_inter(jj)%coor(2)>p%point_inter(jj+1)%coor(2) )then
                            temp = p%point_inter(jj+1)
                            p%point_inter(jj+1) = p%point_inter(jj)
                            p%point_inter(jj) = temp
                        endif
                    enddo
                enddo
            else
                do ii = 1,mx+my-1
                    do jj = kstart+1,kstart+mx+my-ii
                        if( p%point_inter(jj)%coor(2)<p%point_inter(jj+1)%coor(2) )then
                            temp = p%point_inter(jj+1)
                            p%point_inter(jj+1) = p%point_inter(jj)
                            p%point_inter(jj) = temp
                        endif
                    enddo
                enddo
            endif
        endif
    endif

    ia2 = 2*ia; ib2 = 2*ib; ic2 = 2*ic; id2 = 2*id;
    !if( mx+my>0 )then
    if( ia2<=ib2 )then
        do kk = kstart+1,kstart+mx+my
            if( p%point_inter(kk)%ixy_type == 2 )then
                if( p%point_inter(kk)%idx2<p%point_inter(kk-1)%idx2 )then
                    p%point_inter(kk)%idx2 = (p%point_inter(kk-1)%idx2 +1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                elseif( p%point_inter(kk)%idx2>p%point_inter(kk-1)%idx2 .and. p%point_inter(kk)%idx2>p%point_inter(kk+1)%idx2 )then
                    p%point_inter(kk)%idx2 = (p%point_inter(kk)%idx2-1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                endif

            endif

        enddo
        !******************reverse correction
        do kk = kstart+mx+my,kstart+1,-1
            if( p%point_inter(kk)%ixy_type == 2 )then
                if( p%point_inter(kk+1)%idx2<p%point_inter(kk)%idx2 )then
                    p%point_inter(kk)%idx2 = (p%point_inter(kk+1)%idx2 -1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                elseif( p%point_inter(kk+1)%idx2>p%point_inter(kk)%idx2 .and. p%point_inter(kk)%idx2<p%point_inter(kk-1)%idx2 )then
                    p%point_inter(kk)%idx2 = (p%point_inter(kk)%idx2+2)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                endif

            endif

        enddo            

    elseif( ia2>ib2 )then
        do kk = kstart+1, kstart+mx+my
            if( p%point_inter(kk)%ixy_type == 2 )then
                if( p%point_inter(kk)%idx2>p%point_inter(kk-1)%idx2 )then
                    p%point_inter(kk)%idx2 = (p%point_inter(kk-1)%idx2 -1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                elseif( p%point_inter(kk)%idx2<p%point_inter(kk-1)%idx2 .and. p%point_inter(kk)%idx2<p%point_inter(kk+1)%idx2 )then
                    p%point_inter(kk)%idx2 = (p%point_inter(kk)%idx2+2)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                endif
            endif

        enddo
        !******************reverse correction
        do kk = kstart+mx+my, kstart+1,-1
            if( p%point_inter(kk)%ixy_type == 2 )then
                if( p%point_inter(kk+1)%idx2>p%point_inter(kk)%idx2 )then
                    p%point_inter(kk)%idx2 = (p%point_inter(kk+1)%idx2 +1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                elseif( p%point_inter(kk+1)%idx2<p%point_inter(kk)%idx2 .and. p%point_inter(kk)%idx2>p%point_inter(kk-1)%idx2 )then
                    p%point_inter(kk)%idx2 = (p%point_inter(kk)%idx2-1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                endif
            endif

        enddo        
    endif

    !************************
    if( ic2<=id2 )then
        do kk = kstart+1, kstart+mx+my
            if( p%point_inter(kk)%ixy_type == 1 )then
                if( p%point_inter(kk)%idy2<p%point_inter(kk-1)%idy2 )then
                    p%point_inter(kk)%idy2 = (p%point_inter(kk-1)%idy2 +1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
                elseif( p%point_inter(kk)%idy2>p%point_inter(kk-1)%idy2 .and. p%point_inter(kk)%idy2>p%point_inter(kk+1)%idy2 )then
                    p%point_inter(kk)%idy2 = (p%point_inter(kk)%idy2-1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
                endif
            endif

        enddo
        !******************reverse correction
        do kk = kstart+mx+my,kstart+1,-1
            if( p%point_inter(kk)%ixy_type == 1 )then
                if( p%point_inter(kk+1)%idy2<p%point_inter(kk)%idy2 )then
                    p%point_inter(kk)%idy2 = (p%point_inter(kk+1)%idy2 -1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
                elseif( p%point_inter(kk+1)%idy2>p%point_inter(kk)%idy2 .and. p%point_inter(kk)%idy2<p%point_inter(kk-1)%idy2 )then
                    p%point_inter(kk)%idy2 = (p%point_inter(kk)%idy2+2)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
                endif

            endif

        enddo  
    elseif( ic2>id2 )then
        do kk = kstart+1,kstart+mx+my
            if( p%point_inter(kk)%ixy_type == 1 )then
                if( p%point_inter(kk)%idy2>p%point_inter(kk-1)%idy2 )then
                    p%point_inter(kk)%idy2 = (p%point_inter(kk-1)%idy2 -1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
                elseif( p%point_inter(kk)%idy2<p%point_inter(kk-1)%idy2 .and. p%point_inter(kk)%idy2<p%point_inter(kk+1)%idy2 )then
                    p%point_inter(kk)%idy2 = (p%point_inter(kk)%idy2+2)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
                endif
            endif

        enddo
        !******************reverse correction
        do kk = kstart+mx+my, kstart+1,-1
            if( p%point_inter(kk)%ixy_type == 1 )then
                if( p%point_inter(kk+1)%idy2>p%point_inter(kk)%idy2 )then
                    p%point_inter(kk)%idy2 = (p%point_inter(kk+1)%idy2 +1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
                elseif( p%point_inter(kk+1)%idy2<p%point_inter(kk)%idy2 .and. p%point_inter(kk)%idy2>p%point_inter(kk-1)%idy2 )then
                    p%point_inter(kk)%idy2 = (p%point_inter(kk)%idy2-1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
                endif
            endif

        enddo 
    endif
    !endif


    !*********************************************************************************
    ! check the monotonicity
    ia2 = 2*ia; ib2 = 2*ib; ic2 = 2*ic; id2 = 2*id;

    if( mx+my>0 )then
        if(ia2<=ib2)then
            if( p%point_inter(kstart+1)%idx2<ia2 )then
                write(911,*) 'lose monotone!'
                print *,'lose monotone!'
                pause
            endif
            do k = kstart+2,kstart+mx+my
                if( p%point_inter(k)%idx2<p%point_inter(k-1)%idx2 )then
                    write(911,*) 'lose monotone!'
                    print *,i,j
                    print *,vorigin(1:2),vend(1:2),mx,my
                    print *,ia,ib,ic,id
                    print *,p%point_inter(k)%idx2,p%point_inter(k-1)%idx2,'k=',k
                    print *,p%point_inter(k)%ixy_type,p%point_inter(k-1)%ixy_type,'k=',k
                    print *,'lose monotone!3'
                    pause
                endif
            enddo
            if( p%point_inter(kstart+mx+my)%idx2>ib2 )then
                write(911,*) 'lose monotone!'
            
                print *,'lose monotone!4'
                pause
            endif
        elseif( ia2>ib2 )then
            if( p%point_inter(kstart+1)%idx2>ia2 )then
                write(911,*) 'lose monotone!'
                print *,'lose monotone!'
                pause
            endif
            do k = kstart+2,kstart+mx+my
                if( p%point_inter(k)%idx2>p%point_inter(k-1)%idx2 )then
                    write(911,*) 'lose monotone!'
                    print *,i,j
                    print *,vorigin(1:2),vend(1:2),mx,my
                    print *,ia,ib,ic,id
                    print *,p%point_inter(k)%idx2,p%point_inter(k-1)%idx2,'k=',k
                    print *,p%point_inter(k)%ixy_type,p%point_inter(k-1)%ixy_type,'k=',k
                    print *,'lose monotone!2'
                    pause
                endif
            enddo
            if( p%point_inter(kstart+mx+my)%idx2<ib2 )then
                write(911,*) 'lose monotone!'
                print *,'lose monotone!'
                pause
            endif
        endif
        !*************************************************
        if(ic2<=id2)then
            if( p%point_inter(kstart+1)%idy2<ic2 )then
                write(911,*) 'lose monotone!'
                print *,'lose monotone!'
                pause
            endif
            do k = kstart+2,kstart+mx+my
                if( p%point_inter(k)%idy2<p%point_inter(k-1)%idy2 )then
                    write(911,*) 'lose monotone!'
                    print *,i,j
                    print *,vorigin(1:2),vend(1:2),mx,my
                    print *,ia,ib,ic,id
                    print *,p%point_inter(k)%idy2,p%point_inter(k-1)%idy2,'k=',k
                    print *,'lose monotone!1',p%point_inter(1)%idx2,ia2
                    print *,'lose monotone!1'
                    pause
                endif
            enddo
            if( p%point_inter(kstart+mx+my)%idy2>id2 )then
                write(911,*) 'lose monotone!'
                print *,'lose monotone!'
                pause
            endif
        elseif( ic2>id2 )then
            if( p%point_inter(kstart+1)%idy2>ic2 )then
                write(911,*) 'lose monotone!'
                print *,'lose monotone!'
                pause
            endif
            do k = kstart+2,kstart+mx+my
                if( p%point_inter(k)%idy2>p%point_inter(k-1)%idy2 )then
                    write(911,*) 'lose monotone!'
                    print *,'lose monotone!'
                    pause
                endif
            enddo
            if( p%point_inter(kstart+mx+my)%idy2<id2 )then
                write(911,*) 'lose monotone!'
                print *,'lose monotone!'
                pause
            endif
        endif
    endif
    !*********************************************************************************

    end subroutine get_mono_inters
    
    
    
    !*********************************************************************************
    real function xside( v1, v2, y)
    implicit none
    real,intent(in) :: v1(2),v2(2)
    real,intent(in) :: y
    real :: r

    r =  ( y- v1(2) )/(v2(2)-v1(2) )

    if( r>1. )then

        xside = v2(1)

    elseif(r<0.)then
        xside = v1(1)
    else
        xside = ( v2(1)-v1(1) )*r + v1(1)
    endif

    end function xside
    !*********************************************************************************
    real function yside( v1, v2, x)
    implicit none
    real,intent(in) :: v1(2),v2(2)
    real,intent(in) :: x
    real :: r

    r = ( x- v1(1) )/(v2(1)-v1(1) )

    if(r>1)then
        yside = v2(2)
    elseif(r<0)then
        yside = v1(2)
    else
        yside = ( v2(2)-v1(2) )*r + v1(2)
    endif

    end function yside