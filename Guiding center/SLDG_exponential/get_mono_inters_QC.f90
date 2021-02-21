    subroutine get_mono_inters_QC(ia,ib,ic,id,mx,my,vorigin,vend,xi_o,xi_e,p,kstart,v11,v33)
    implicit none
    integer,intent(in) :: ia,ib,ic,id,mx,my
    real,intent(in) :: vorigin(1:2),vend(1:2),xi_o,xi_e
    real,intent(in) :: v11(1:2),v33(1:2)
    integer,intent(inout) :: kstart
    integer :: ii,jj
    integer :: k,inter,kk,ki
    integer :: ia2,ib2,ic2,id2
    type(type_face),pointer :: p
    integer :: kend

    real :: cxi
    real :: t,tc
    !debug
    real :: td

    integer :: istraight

    real :: t11,t22

    !**************************************************************************
    ! Algorithm of getting the intersection points between face and gridlines
    ! find intersections between face and x=x_i


    if(mx .ne. 0)then

        do k = 1 , mx
            if(ia<ib) then
                inter = ia + k
            else
                inter = ia + 1 - k
            endif

            kk = kstart +k
            p%point_inter( kk)%coor( 1 ) = xgrid( inter )
            call get_kexi_xgrid( v11(1:2), v33(1:2),p%et2(1:2), xgrid( inter ),xi_o,xi_e,cxi,istraight  )


            t = quadratic_function( p%ya,p%yb,p%yc,cxi )
            !t11 = quadratic_function( p%ya,p%yb,p%yc,cxi )
            if( istraight==1 .or. p%ist==1 )then
                t11 = t
                call get_kexi_xgrid_s(v11(1:2), v33(1:2),p%et2(1:2), xgrid( inter ),xi_o,xi_e,cxi)
                !t = reverse_to_y( vorigin(1:2), vend(1:2),cxi,p%et2(2) )
                t22 = reverse_to_y( v11(1:2), v33(1:2),cxi,p%et2(2) )

                if( abs(t11-t22)<0.001*dy   )then
                    t = t11
                        !print *,abs(t11-t22),0.1*dx**2
                        !pause
                else
                    t = t22
                endif
                 t = t22
                t= t11
                 !pause
            endif



            call confine( vorigin(2), vend(2),t,tc )
            p%point_inter( kk)%coor( 2 ) = tc
            !p%point_inter( k)%coor( 2 ) = yside( vorigin(1:2), vend(1:2),  xgrid(inter) )

            td = yside( vorigin(1:2), vend(1:2),  xgrid(inter) )

            !if( abs(tc-td)>6d-2 )then
            !    print *,tc,td,abs(tc-td)
            !    pause
            !
            !endif
            p%point_inter( kk)%xxii = cxi
            p%point_inter( kk)%ixy_type = 1
            p%point_inter( kk)%igrid = inter
            p%point_inter( kk)%idx2 = 2*inter-1
            p%point_inter( kk)%id_xy = ceiling( (p%point_inter( kk)%coor( 2 )-ybottom)/dy )
            p%point_inter( kk)%idy2 = 2*p%point_inter( kk)%id_xy
            
            !if(time>0.0755 .and. i==155 .and. j==108 .and. io==3 )then
            !    print *,p%point_inter( kk)
            !    pause
            !endif
            
        enddo
    endif
    !find intersections between face and y=y_j

    if(my .ne. 0)then
        !print *,vorigin(1:2)
        !print *,vend(1:2)
        !pause
        do k = 1 + mx, mx+my
            ki = k - mx
            if(ic<id)then
                inter = ic + ki
            else
                inter = ic  + 1 - ki
            endif

            kk = kstart +k
            p%point_inter( kk)%coor( 2 ) = ygrid( inter )
            call get_kexi_ygrid( v11(1:2), v33(1:2),p%et2(1:2), ygrid( inter ),xi_o,xi_e,cxi,istraight  )
            t = quadratic_function( p%xa,p%xb,p%xc,cxi )
            !t11 = quadratic_function( p%xa,p%xb,p%xc,cxi )
            if( istraight==1 .or. p%ist==1 )then
                t11 =t
                call get_kexi_ygrid_s(v11(1:2), v33(1:2),p%et2(1:2), ygrid( inter ),xi_o,xi_e,cxi)
                !t = reverse_to_x( vorigin(1:2), vend(1:2),cxi,p%et2(2) )
                t22 = reverse_to_x( v11(1:2), v33(1:2),cxi,p%et2(2) )

                if( abs(t11-t22)<0.001*dx  )then
                    t = t11
                    !print *,abs(t11-t22),0.01*dy**3
                    !pause
                else
                    t = t22
                endif
                t = t22
                t =t11
                !pause
            endif




            call confine( vorigin(1), vend(1),t,tc )
            p%point_inter( kk)%coor( 1 ) = tc
            !p%point_inter( k)%coor( 1 ) = xside( vorigin(1:2), vend(1:2),  ygrid(inter) )
            td = xside( vorigin(1:2), vend(1:2),  ygrid(inter) )

            !if( abs(tc-td)>6d-2 )then
            !    print *,tc,td,abs(tc-td),1,i,j
            !    pause
            !endif

            p%point_inter( kk)%xxii = cxi
            p%point_inter( kk)%ixy_type = 2
            p%point_inter( kk)%igrid = inter
            p%point_inter( kk)%idy2 = 2*inter-1
            p%point_inter( kk)%id_xy = ceiling( (p%point_inter( kk)%coor( 1 )-xleft)/dx )
            p%point_inter( kk)%idx2 = 2*p%point_inter( kk)%id_xy

            !if(time>0.0755 .and. i==155 .and. j==108 .and. io==3 )then
            !    print *,p%point_inter( kk)
            !    pause
            !endif


        enddo
    endif
    !**************************************************************************
    !p%nsub_outer = 1+ mx+my
    !p%point_inter( 0)%coor( 1 ) = vorigin(1)
    !p%point_inter( 0)%coor( 2 ) = vorigin(2)
    !p%point_inter( 0)%ixy_type = 0
    !p%point_inter( 0)%idx2 = 2*ia
    !p%point_inter( 0)%idy2 = 2*ic
    !!******
    !p%point_inter( 1+ mx+my)%coor( 1 ) = vend(1)
    !p%point_inter( 1+ mx+my)%coor( 2 ) = vend(2)
    !p%point_inter( 1+ mx+my)%ixy_type = 0
    !p%point_inter( 1+ mx+my)%idx2 = 2*ib
    !p%point_inter( 1+ mx+my)%idy2 = 2*id
    !!**************************************************************************
    ! Algorithm of sorting point_inter(1:mx+my) based on coordinates
    !if(mx+my>1)then
    !    if(mx>=my)then
    !        if( vorigin(1)<vend(1) )then
    !            ! sort them in a increasing order
    !            do ii = 1 ,mx+my-1
    !                do jj = kstart+1 , kstart+mx+my-ii
    !                    if( p%point_inter(jj)%coor(1)>p%point_inter(jj+1)%coor(1) )then
    !                        temp = p%point_inter(jj+1)
    !                        p%point_inter(jj+1) = p%point_inter(jj)
    !                        p%point_inter(jj) = temp
    !                    endif
    !                enddo
    !            enddo
    !        else
    !            do ii = 1,mx+my-1
    !                do jj = kstart+1 , kstart+mx+my-ii
    !                    if( p%point_inter(jj)%coor(1)<p%point_inter(jj+1)%coor(1) )then
    !                        temp = p%point_inter(jj+1)
    !                        p%point_inter(jj+1) = p%point_inter(jj)
    !                        p%point_inter(jj) = temp
    !                    endif
    !                enddo
    !            enddo
    !        endif
    !    elseif(mx<my)then !mx>=my
    !        if( vorigin(2)<vend(2) )then
    !            do ii = 1, mx+my-1
    !                do jj = kstart+1 , kstart+mx+my-ii
    !                    if( p%point_inter(jj)%coor(2)>p%point_inter(jj+1)%coor(2) )then
    !                        temp = p%point_inter(jj+1)
    !                        p%point_inter(jj+1) = p%point_inter(jj)
    !                        p%point_inter(jj) = temp
    !                    endif
    !                enddo
    !            enddo
    !        else
    !            do ii = 1,mx+my-1
    !                do jj = kstart+1 , kstart+mx+my-ii
    !                    if( p%point_inter(jj)%coor(2)<p%point_inter(jj+1)%coor(2) )then
    !                        temp = p%point_inter(jj+1)
    !                        p%point_inter(jj+1) = p%point_inter(jj)
    !                        p%point_inter(jj) = temp
    !                    endif
    !                enddo
    !            enddo
    !        endif
    !    endif
    !endif

    !
    if(mx+my>1)then
         !sort them in a increasing order
        do ii = 1 ,mx+my-1
            do jj = kstart+1 , kstart+mx+my-ii
                if( p%point_inter(jj)%xxii>p%point_inter(jj+1)%xxii )then
                    temp = p%point_inter(jj+1)
                    p%point_inter(jj+1) = p%point_inter(jj)
                    p%point_inter(jj) = temp
                endif
            enddo
        enddo
    endif

    ia2 = 2*ia; ib2 = 2*ib; ic2 = 2*ic; id2 = 2*id;
    !if( mx+my>0 )then
    if( ia2<=ib2 )then
        do kk = kstart+1,kstart+mx+my
            if( p%point_inter(kk)%ixy_type == 2 )then
                if( p%point_inter(kk)%idx2<p%point_inter(kk-1)%idx2 )then
                    p%point_inter(kk)%idx2 = p%point_inter(kk)%idx2 +2 
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                    
                    !if( p%point_inter(kk)%coor(1)<p%point_inter(kk-1)%coor(1) )then
                    !
                    !    pause
                    !endif
                    
                elseif( p%point_inter(kk)%idx2>p%point_inter(kk-1)%idx2 .and. p%point_inter(kk)%idx2>p%point_inter(kk+1)%idx2 )then
                    p%point_inter(kk)%idx2 = p%point_inter(kk)%idx2- 2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                    
                    !if( p%point_inter(kk)%coor(1)>p%point_inter(kk+1)%coor(1) )then
                    !    print *,mx,my
                    !    print *,1
                    !    pause
                    !endif
                    
                endif

            endif

        enddo
        !******************reverse correction
        do kk = kstart+mx+my,kstart+1,-1
            if( p%point_inter(kk)%ixy_type == 2 )then
                if( p%point_inter(kk+1)%idx2<p%point_inter(kk)%idx2 )then
                    p%point_inter(kk)%idx2 =  p%point_inter(kk)%idx2 -2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                    
                    !if( p%point_inter(kk+1)%coor(1)<p%point_inter(kk)%coor(1) )then
                    !    print *,mx,my
                    !    print *,2
                    !    pause
                    !endif
                    
                elseif( p%point_inter(kk+1)%idx2>p%point_inter(kk)%idx2 .and. p%point_inter(kk)%idx2<p%point_inter(kk-1)%idx2 )then
                    p%point_inter(kk)%idx2 =  p%point_inter(kk)%idx2+2 
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                    
                    !if( p%point_inter(kk)%coor(1)<p%point_inter(kk-1)%coor(1) )then
                    !    print *,mx,my
                    !    print *,3
                    !    pause
                    !endif
                    
                endif

            endif

        enddo

    elseif( ia2>ib2 )then
        do kk = kstart+1, kstart+mx+my
            if( p%point_inter(kk)%ixy_type == 2 )then
                if( p%point_inter(kk)%idx2>p%point_inter(kk-1)%idx2 )then
                    p%point_inter(kk)%idx2 =  p%point_inter(kk)%idx2 -2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                    
                    !if( p%point_inter(kk)%coor(1)>p%point_inter(kk-1)%coor(1) )then
                    !    print *,mx,my
                    !    print *,4
                    !    pause
                    !endif                    
                    
                elseif( p%point_inter(kk)%idx2<p%point_inter(kk-1)%idx2 .and. p%point_inter(kk)%idx2<p%point_inter(kk+1)%idx2 )then
                    p%point_inter(kk)%idx2 = p%point_inter(kk)%idx2+2 
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                    
                    !if( p%point_inter(kk)%coor(1)<p%point_inter(kk+1)%coor(1) )then
                    !    print *,mx,my
                    !    print *,5
                    !    pause
                    !endif                     
                    
                endif
            endif

        enddo
        !******************reverse correction
        do kk = kstart+mx+my, kstart+1,-1
            if( p%point_inter(kk)%ixy_type == 2 )then
                if( p%point_inter(kk+1)%idx2>p%point_inter(kk)%idx2 )then
                    p%point_inter(kk)%idx2 = p%point_inter(kk)%idx2 +2 
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                    
                    !if( p%point_inter(kk+1)%coor(1)<p%point_inter(kk)%coor(1) )then
                    !    print *,mx,my
                    !    print *,6
                    !    pause
                    !endif 
                    
                elseif( p%point_inter(kk+1)%idx2<p%point_inter(kk)%idx2 .and. p%point_inter(kk)%idx2>p%point_inter(kk-1)%idx2 )then
                    p%point_inter(kk)%idx2 = p%point_inter(kk)%idx2-2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                    
                    !if( p%point_inter(kk)%coor(1)>p%point_inter(kk-1)%coor(1) )then
                    !    print *,mx,my
                    !    print *,7
                    !    pause
                    !endif 
                    
                endif
            endif

        enddo
    endif

    !************************
    if( ic2<=id2 )then
        do kk = kstart+1, kstart+mx+my
            if( p%point_inter(kk)%ixy_type == 1 )then
                if( p%point_inter(kk)%idy2<p%point_inter(kk-1)%idy2 )then
                    p%point_inter(kk)%idy2 = p%point_inter(kk)%idy2 +2 
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
                    
                    !if( p%point_inter(kk)%coor(2)<p%point_inter(kk-1)%coor(2) )then
                    !    print *,mx,my
                    !    print *,8
                    !    pause
                    !endif                     
                    
                elseif( p%point_inter(kk)%idy2>p%point_inter(kk-1)%idy2 .and. p%point_inter(kk)%idy2>p%point_inter(kk+1)%idy2 )then
                    p%point_inter(kk)%idy2 = p%point_inter(kk)%idy2-2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
                    
                    !if( p%point_inter(kk)%coor(2)>p%point_inter(kk+1)%coor(2) )then
                    !    print *,mx,my,p%point_inter(kk)%coor(2),p%point_inter(kk+1)%coor(2)
                    !    print *,abs(p%point_inter(kk)%coor(2)-p%point_inter(kk+1)%coor(2))
                    !    print *,9
                    !    pause
                    !endif                     
                    
                endif
            endif

        enddo
        !******************reverse correction
        do kk = kstart+mx+my,kstart+1,-1
            if( p%point_inter(kk)%ixy_type == 1 )then
                if( p%point_inter(kk+1)%idy2<p%point_inter(kk)%idy2 )then
                    p%point_inter(kk)%idy2 = p%point_inter(kk)%idy2 -2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
                    
                    !if( p%point_inter(kk+1)%coor(2)<p%point_inter(kk)%coor(2) )then
                    !    print *,mx,my
                    !    print *,10
                    !    pause
                    !endif                     
                    
                elseif( p%point_inter(kk+1)%idy2>p%point_inter(kk)%idy2 .and. p%point_inter(kk)%idy2<p%point_inter(kk-1)%idy2 )then
                    p%point_inter(kk)%idy2 = p%point_inter(kk)%idy2+2 
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
                    
                    !if( p%point_inter(kk)%coor(2)<p%point_inter(kk-1)%coor(2) )then
                    !    print *,mx,my
                    !    print *,11
                    !    pause
                    !endif 
                    
                endif

            endif

        enddo
    elseif( ic2>id2 )then
        do kk = kstart+1,kstart+mx+my
            if( p%point_inter(kk)%ixy_type == 1 )then
                if( p%point_inter(kk)%idy2>p%point_inter(kk-1)%idy2 )then
                    p%point_inter(kk)%idy2 = p%point_inter(kk)%idy2 -2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
                    
                    !if( p%point_inter(kk)%coor(2)>p%point_inter(kk-1)%coor(2) )then
                    !    print *,mx,my
                    !    print *,12
                    !    pause
                    !endif 
                    
                elseif( p%point_inter(kk)%idy2<p%point_inter(kk-1)%idy2 .and. p%point_inter(kk)%idy2<p%point_inter(kk+1)%idy2 )then
                    p%point_inter(kk)%idy2 = p%point_inter(kk)%idy2+2 
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
                    
                    !if( p%point_inter(kk)%coor(2)<p%point_inter(kk+1)%coor(2) )then
                    !    print *,mx,my
                    !    print *,13
                    !    pause
                    !endif                     
                    
                endif
            endif

        enddo
        !******************reverse correction
        do kk = kstart+mx+my, kstart+1,-1
            if( p%point_inter(kk)%ixy_type == 1 )then
                if( p%point_inter(kk+1)%idy2>p%point_inter(kk)%idy2 )then
                    p%point_inter(kk)%idy2 = p%point_inter(kk)%idy2 +2 
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
                    
                    !if( p%point_inter(kk+1)%coor(2)>p%point_inter(kk)%coor(2) )then
                    !    print *,mx,my
                    !    print *,14
                    !    pause
                    !endif                   
                    
                elseif( p%point_inter(kk+1)%idy2<p%point_inter(kk)%idy2 .and. p%point_inter(kk)%idy2>p%point_inter(kk-1)%idy2 )then
                    p%point_inter(kk)%idy2 = p%point_inter(kk)%idy2-2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
                    
                    !if( p%point_inter(kk)%coor(2)>p%point_inter(kk-1)%coor(2) )then
                    !    print *,mx,my
                    !    print *,15
                    !    pause
                    !endif                     
                    
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
                print *,'lose monotone!5'
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
                print *,i,j
                print *,vorigin(1:2),vend(1:2),mx,my
                print *,ia,ib,ic,id
                print *,p%point_inter(kend)%idx2,p%point_inter(kend-1)%idx2,'k=',kend,kstart
                print *,p%point_inter(kend)%ixy_type,p%point_inter(kend-1)%ixy_type,'k=',kend
                print *,'lose monotone!4'
                pause
            endif
        elseif( ia2>ib2 )then
            if( p%point_inter(kstart+1)%idx2>ia2 )then
                write(911,*) 'lose monotone!'
                print *,'lose monotone!6'
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
                print *,'lose monotone!7'
                pause
            endif
        endif
        !*************************************************
        if(ic2<=id2)then
            if( p%point_inter(kstart+1)%idy2<ic2 )then
                write(911,*) 'lose monotone!'
                    print *,i,j
                    print *,vorigin(1:2),vend(1:2),mx,my
                    print *,ia,ib,ic,id
                    print *,p%point_inter(kstart+1)%idy2,ic2,'k=',kstart+1,io
                    print *,p%point_inter(kstart+1)
                print *,'lose monotone!8'
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
                pause
            endif
        elseif( ic2>id2 )then
            if( p%point_inter(kstart+1)%idy2>ic2 )then
                write(911,*) 'lose monotone!'
                pause
            endif
            do k = kstart+2,kstart+mx+my
                if( p%point_inter(k)%idy2>p%point_inter(k-1)%idy2 )then
                    write(911,*) 'lose monotone!'
                    pause
                endif
            enddo
            if( p%point_inter(kstart+mx+my)%idy2<id2 )then
                write(911,*) 'lose monotone!'
                pause
            endif
        endif
    endif
    !*********************************************************************************

    end subroutine get_mono_inters_QC

    !*********************************************************************************

    subroutine get_kexi_xgrid( v1,v3,e2, x, xorigin,xend, cxi,istraight )
    implicit none
    real,intent(in) :: v1(2),v3(2),e2(2),x,xorigin,xend
    real,intent(out) :: cxi
    integer,intent(out) :: istraight

    real :: at,bt,ct
    real :: t1,t2,t3
    real :: xi(2),eta(1:2)
    !real :: xorigin,xend

    if( abs( v3(1)-v1(1) )<=abs( v3(2)-v1(2) ) )then
        at = e2(2)/( e2(1)**2-1. )
        bt = ( v3(1)-v1(1) )/( v3(2)-v1(2) )
        ct =  2.*( ( v1(1)+v3(1) )/2.-x)/( v3(2)-v1(2) )-at
        
        if( abs(at)>eps )then 
        call get_quadratic_roots( at,bt,ct,xi(1:2) )
        elseif( abs(at)<=eps .and. abs(bt)>eps )then
       
            xi(1:2) = -ct/bt 
        else
            pause
            xi(1:2) = xorigin*0.5 + xend*0.5
        endif
        !pause
    else
        t1 = e2(2)/( e2(1)**2-1. )
        t2 = ( v3(2)-v1(2) )/( v3(1)-v1(1) )
        t3 = ( x-0.5*v1(1)-0.5*v3(1) )/( v3(1)-v1(1) )
        at = t1*t2**2
        bt = -1. - 4.*t3*t2*t1
        ct = t1*( 4.*t3**2-1. )
        
        if( abs(at)>eps )then
        call get_quadratic_roots( at,bt,ct,eta(1:2) )
        xi(1) = 2*t3-t2*eta(1)
        xi(2) = 2*t3-t2*eta(2)
        elseif( abs(at)<=eps .and. abs(bt)>eps )then
            eta(1:2) = -ct/bt
            xi(1:2) = 2.*t3-t2*eta(1:2)
        else
            pause
            xi(1:2) = xorigin*0.5 + xend*0.5
        endif        

    endif

    ! pick one root!
    !xorigin = -1.
    !xend = 1.
    call pick_1root( xi(1:2),xorigin,xend,cxi,istraight )
    
    if( bt**2 - 4.*at*ct <-eps )then
        istraight = 1
        print *,221, bt**2 - 4.*at*ct,-eps, at
        pause
    endif    
    
    if(istraight==1)then
        print *,v1(1:2),v3(1:2),e2(1:2),x,xorigin,xend
        print *,xi(1:2),xorigin,xend,cxi,istraight
        print *,abs( v3(1)-v1(1) ),abs( v3(2)-v1(2) )
        print *,eta(1:2)
        print *,2222,t3,t2,eta(1)
        pause
    endif
    end subroutine get_kexi_xgrid
    !*****************************************
    subroutine get_kexi_ygrid( v1,v3,e2, y,xorigin,xend,cxi,istraight )
    implicit none
    real,intent(in) :: v1(2),v3(2),e2(2),y,xorigin,xend
    real,intent(out) :: cxi
    integer,intent(out) :: istraight

    real :: at,bt,ct
    real :: t1,t2,t3
    real :: xi(2),eta(1:2)
    !real :: xorigin,xend

    if( abs( v3(1)-v1(1) )<=abs( v3(2)-v1(2) ) )then
        t1 = e2(2)/( e2(1)**2-1. )
        t2 = ( v3(1)-v1(1) )/( v3(2)-v1(2) )
        t3 = ( y-0.5*v1(2)-0.5*v3(2) )/( v3(2)-v1(2) )
        at = t1*t2**2
        bt = -1.+4.*t3*t2*t1
        ct =  4.*t1*t3**2  -t1
        
        if( abs(at)>eps )then
        call get_quadratic_roots( at,bt,ct,eta(1:2) )
        xi(1) = 2.*t3+t2*eta(1)
        xi(2) = 2.*t3+t2*eta(2)
        elseif( abs(at)<=eps .and. abs(bt)>eps )then
            eta(1:2) = -ct/bt
            xi(1:2) = 2.*t3+t2*eta(1:2)
        else
            pause
            xi(1:2) = xorigin*0.5 + xend*0.5
        endif

    else

        at = - e2(2)/( e2(1)**2-1. )
        bt = ( v3(2)-v1(2) )/(v3(1)-v1(1) )
        ct = 2*( ( v1(2)+v3(2) )/2.-y)/(v3(1)-v1(1) )-at
        if( abs(at)>eps )then 
        call get_quadratic_roots( at,bt,ct,xi(1:2) )
        elseif( abs(at)<=eps .and. abs(bt)>eps )then
       
            xi(1:2) = -ct/bt 
        else
            pause
            xi(1:2) = xorigin*0.5 + xend*0.5
        endif
        !pause
    endif

    ! pick one root!
    !xorigin = -1.
    !xend = 1.
    
    call pick_1root( xi(1:2),xorigin,xend,cxi,istraight )
    

    
    
    !if(istraight==1)then
    !    print *,xi(1:2),xorigin,xend,cxi,istraight
    !    pause
    !endif
    
    
    if( bt**2 - 4.*at*ct <-eps )then
        istraight = 1
        print *,abs( v3(1)-v1(1) ),abs( v3(2)-v1(2) )
        print *,223, bt**2 - 4.*at*ct,-eps, at
        pause
    endif
    
    end subroutine get_kexi_ygrid
    !*****************************************
    subroutine pick_1root( xi,xorigin,xend,cxi,istraight )
    implicit none
    real,intent(in) :: xi(1:2),xorigin,xend
    real,intent(out) :: cxi
    integer :: istraight

    istraight = 0
    if( ( xi(1) - xorigin )*( xi(1) - xend ) <=0. .and. ( xi(2) - xorigin )*( xi(2) - xend )>0.  )then
        cxi = xi(1)
    elseif( ( xi(1) - xorigin )*( xi(1) - xend ) >0. .and. ( xi(2) - xorigin )*( xi(2) - xend )<=0.  )then
        cxi = xi(2)
    elseif( ( xi(1) - xorigin )*( xi(1) - xend ) <=0. .and. ( xi(2) - xorigin )*( xi(2) - xend )<=0. )then
        cxi = xi(1)*0.5 + xi(2)*0.5
        
        !if( abs(xi(1)-xi(2))>10*eps )then
        !print *,2233,xi(1:2),xorigin,xend
        !pause
        !endif
    elseif( xi(1)>xend .and. xi(2)>xend )then
        cxi = xend
    elseif( xi(1)<xorigin .and. xi(2)<xorigin )then
        cxi = xorigin
        ! special cases
    elseif(  xi(1)<xorigin .and. xi(2)>xend )then
        if( abs( xi(1)-xorigin )< abs( xi(2)-xend) )then
        cxi = xorigin
        elseif( abs( xi(1)-xorigin )>= abs( xi(2)-xend) )then
            cxi = xend
        endif
    elseif(  xi(2)<xorigin .and. xi(1)>xend)then
        if( abs( xi(2)-xorigin )< abs( xi(1)-xend) )then
        cxi = xorigin
        elseif( abs( xi(2)-xorigin )>= abs( xi(1)-xend) )then
            cxi = xend
        endif
 
    !else
    !    pause
    !    cxi = xorigin*0.5 + xend*0.5
    !    istraight = 1
        
        !print *,xi(1:2),xorigin,xend
        
        !pause
    endif

    end subroutine pick_1root
    !******************************************
    subroutine confine( t1,t2,t,tc )
    implicit none
    real,intent(in) :: t1,t2,t
    real,intent(out) :: tc
    real :: tmin,tmax

    tmin = min( t1,t2 )
    tmax = max( t1,t2 )
    if( t<tmin )then
        tc = tmin
    elseif( t>tmax )then
        tc = tmax
    else
        tc = t
    endif

    end subroutine confine
    !******************************************
    subroutine confine3( t1,t2,t3,t,tc )
    implicit none
    real,intent(in) :: t1,t2,t3,t
    real,intent(out) :: tc
    real :: tmin,tmax

    tmin = min( t1,t2,t3 )
    tmax = max( t1,t2,t3 )
    if( t<tmin )then
        tc = tmin
    elseif( t>tmax )then
        tc = tmax
    else
        tc = t
    endif

    end subroutine confine3
    !******************************************
    subroutine get_quadratic_roots( a,b,c, xi )
    implicit none
    real,intent(in) :: a,b,c
    real,intent(out) :: xi(2)

    real :: gamma,delta

    delta = abs(b**2-4.*a*c)

    if(b>=0.)then
        gamma = 1.
    else
        gamma = -1.
    endif

    xi(1) = 2.*c/( -b - gamma *sqrt(delta) )
    xi(2) = ( -b-gamma*sqrt(delta)  )/a*0.5

    end subroutine get_quadratic_roots




    !****************************************************
    subroutine get_kexi_xgrid_s( v1,v3,e2, x, xorigin,xend, cxi  )
    implicit none
    real,intent(in) :: v1(1:2),v3(1:2),e2(1:2),x,xorigin,xend
    real,intent(out) :: cxi


    real :: xi

    xi = 2.*( x-(v1(1)+v3(1))/2.-(v3(2)-v1(2))/2.*e2(2)   )/( v3(1)-v1(1) )

    call confine(  xorigin,xend,xi,cxi  )

    end subroutine get_kexi_xgrid_s
    !****************************************************
    subroutine get_kexi_ygrid_s( v1,v3,e2, y, xorigin,xend, cxi  )
    implicit none
    real,intent(in) :: v1(1:2),v3(1:2),e2(1:2),y,xorigin,xend
    real,intent(out) :: cxi


    real :: xi

    xi = 2.*( y-(v1(2)+v3(2))/2.+(v3(1)-v1(1))/2.*e2(2)   )/( v3(2)-v1(2) )

    call confine(  xorigin,xend,xi,cxi  )

    end subroutine get_kexi_ygrid_s

    !**********************************************************
    real function reverse_to_x(v1,v3, exi20,eta20)
    implicit none
    real,intent(in) ::  v1(1:2),v3(1:2), exi20,eta20

    reverse_to_x = ( v3(1)-v1(1) )/2.*exi20 + ( v3(2)-v1(2) )/2.*eta20 + ( v1(1)+v3(1) )/2.


    end function reverse_to_x
    !**********************************************************
    real function reverse_to_y( v1,v3, exi20,eta20 )
    implicit none
    real,intent(in) :: v1(1:2),v3(1:2), exi20,eta20

    reverse_to_y = ( v3(2)-v1(2) )/2.*exi20 - (  v3(1)-v1(1) )/2.*eta20 + (  v1(2)+v3(2) )/2.

    end function reverse_to_y
