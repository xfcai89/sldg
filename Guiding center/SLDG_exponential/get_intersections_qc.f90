    subroutine get_intersections_qc(p)
    implicit none
    real :: vorigin(1:2),vend(1:2),v11(1:2),v33(1:2)
    integer :: mx,my
    integer :: ia,ib,ic,id
    integer :: kstart,kend
    
    real :: xi_o,xi_e

    type(type_face),pointer :: p
    
    
    v11(1:2) = p%point_origin%coor(1:2)
    v33(1:2) = p%point_end%coor(1:2)

    p%nsub_outer = 0
    kstart = p%nsub_outer
    if( p%nsubface==1 )then
        ia = p%point_origin%id(1)
        ib = p%point_end%id(1)

        ic = p%point_origin%id(2)
        id = p%point_end%id(2)
        vorigin(1:2) = p%point_origin%coor(1:2)
        vend(1:2) = p%point_end%coor(1:2)

        mx = abs(ib-ia)
        my = abs(id-ic)
        
        kstart = p%nsub_outer

        p%nsub_outer = 1+ mx+my
        p%point_inter( 0)%coor( 1 ) = vorigin(1)
        p%point_inter( 0)%coor( 2 ) = vorigin(2)
        p%point_inter( 0)%ixy_type = 0
        p%point_inter( 0)%idx2 = 2*ia
        p%point_inter( 0)%idy2 = 2*ic
        !******
        !print *,vorigin(1:2),vend(1:2)
        p%point_inter( 1+ mx+my)%coor( 1 ) = vend(1)
        p%point_inter( 1+ mx+my)%coor( 2 ) = vend(2)
        p%point_inter( 1+ mx+my)%ixy_type = 0
        p%point_inter( 1+ mx+my)%idx2 = 2*ib
        p%point_inter( 1+ mx+my)%idy2 = 2*id
        !**************************************************************************
        ! we will got kstart+mx+my-1 - kstart +1 =mx+my intersections
        ! for example, kstart = 1 to kend = mx+my = kstart + mx + my - 1
         
        !call get_mono_inters(ia,ib,ic,id,mx,my,vorigin(1:2),vend(1:2),p,kstart)
        
        xi_o = -1.
        xi_e = 1.
        call get_mono_inters_QC(ia,ib,ic,id,mx,my,vorigin(1:2),vend(1:2),xi_o,xi_e,p,kstart,v11,v33)
        
        
    elseif( p%nsubface==2 )then
        ia = p%point_origin%id(1)
        ib = p%point_mid1%id(1)

        ic = p%point_origin%id(2)
        id = p%point_mid1%id(2)
        vorigin(1:2) = p%point_origin%coor(1:2)
        vend(1:2) = p%point_mid1%coor(1:2)

        mx = abs(ib-ia)
        my = abs(id-ic)
        
        kstart = p%nsub_outer
        
        p%point_inter( p%nsub_outer)%coor( 1 ) = vorigin(1)
        p%point_inter( p%nsub_outer)%coor( 2 ) = vorigin(2)
        p%point_inter( p%nsub_outer)%ixy_type = 0
        p%point_inter( p%nsub_outer)%idx2 = 2*ia
        p%point_inter( p%nsub_outer)%idy2 = 2*ic
        !******
        p%nsub_outer = p%nsub_outer+ 1+ mx+my
        !******
        p%point_inter( p%nsub_outer)%coor( 1 ) = vend(1)
        p%point_inter( p%nsub_outer)%coor( 2 ) = vend(2)
        p%point_inter( p%nsub_outer)%ixy_type = 0
        p%point_inter( p%nsub_outer)%idx2 = 2*ib
        p%point_inter( p%nsub_outer)%idy2 = 2*id
        !**************************************************************************
        ! we will got kstart+mx+my-1 - kstart +1 =mx+my intersections
        ! for example, kstart = 1 to kend = mx+my = kstart + mx + my - 1
        !call get_mono_inters(ia,ib,ic,id,mx,my,vorigin(1:2),vend(1:2),p,kstart)
        
        xi_o = -1.
        xi_e = p%point_mid1%xxii
        call get_mono_inters_QC(ia,ib,ic,id,mx,my,vorigin(1:2),vend(1:2),xi_o,xi_e,p,kstart,v11,v33)
        
         
        
        !***************************************************************************
        !***************************************************************************
        ia = p%point_mid1%id(1)
        ib = p%point_end%id(1)

        ic = p%point_mid1%id(2)
        id = p%point_end%id(2)
        vorigin(1:2) = p%point_mid1%coor(1:2)
        vend(1:2) = p%point_end%coor(1:2)

        mx = abs(ib-ia)
        my = abs(id-ic)
        
        kstart = p%nsub_outer
        p%point_inter( p%nsub_outer)%coor( 1 ) = vorigin(1)
        p%point_inter( p%nsub_outer)%coor( 2 ) = vorigin(2)
        p%point_inter( p%nsub_outer)%ixy_type = 0
        p%point_inter( p%nsub_outer)%idx2 = 2*ia
        p%point_inter( p%nsub_outer)%idy2 = 2*ic
        !******
        p%nsub_outer = p%nsub_outer+ 1+ mx+my
        !******
        p%point_inter( p%nsub_outer)%coor( 1 ) = vend(1)
        p%point_inter( p%nsub_outer)%coor( 2 ) = vend(2)
        p%point_inter( p%nsub_outer)%ixy_type = 0
        p%point_inter( p%nsub_outer)%idx2 = 2*ib
        p%point_inter( p%nsub_outer)%idy2 = 2*id
        !**************************************************************************
        ! we will got kstart+mx+my-1 - kstart +1 =mx+my intersections
        ! for example, kstart = 1 to kend = mx+my = kstart + mx + my - 1
 
        !call get_mono_inters(ia,ib,ic,id,mx,my,vorigin(1:2),vend(1:2),p,kstart)
        
        xi_o = p%point_mid1%xxii
        xi_e = 1.
        call get_mono_inters_QC(ia,ib,ic,id,mx,my,vorigin(1:2),vend(1:2),xi_o,xi_e,p,kstart,v11,v33)
        
 
        
    elseif( p%nsubface==3 )then
        ia = p%point_origin%id(1)
        ib = p%point_mid1%id(1)

        ic = p%point_origin%id(2)
        id = p%point_mid1%id(2)
        vorigin(1:2) = p%point_origin%coor(1:2)
        vend(1:2) = p%point_mid1%coor(1:2)

        mx = abs(ib-ia)
        my = abs(id-ic)
        
        kstart = p%nsub_outer
        p%point_inter( p%nsub_outer)%coor( 1 ) = vorigin(1)
        p%point_inter( p%nsub_outer)%coor( 2 ) = vorigin(2)
        p%point_inter( p%nsub_outer)%ixy_type = 0
        p%point_inter( p%nsub_outer)%idx2 = 2*ia
        p%point_inter( p%nsub_outer)%idy2 = 2*ic
        !******
        p%nsub_outer = p%nsub_outer+ 1+ mx+my
        !******
        p%point_inter( p%nsub_outer)%coor( 1 ) = vend(1)
        p%point_inter( p%nsub_outer)%coor( 2 ) = vend(2)
        p%point_inter( p%nsub_outer)%ixy_type = 0
        p%point_inter( p%nsub_outer)%idx2 = 2*ib
        p%point_inter( p%nsub_outer)%idy2 = 2*id
        !**************************************************************************
        ! we will got kstart+mx+my-1 - kstart +1 =mx+my intersections
        ! for example, kstart = 1 to kend = mx+my = kstart + mx + my - 1
 
        !call get_mono_inters(ia,ib,ic,id,mx,my,vorigin(1:2),vend(1:2),p,kstart)
        
        xi_o = -1.
        xi_e = p%point_mid1%xxii
        call get_mono_inters_QC(ia,ib,ic,id,mx,my,vorigin(1:2),vend(1:2),xi_o,xi_e,p,kstart,v11,v33)
        
 
        !***************************************************************************
        !***************************************************************************
        ia = p%point_mid1%id(1)
        ib = p%point_mid2%id(1)

        ic = p%point_mid1%id(2)
        id = p%point_mid2%id(2)
        vorigin(1:2) = p%point_mid1%coor(1:2)
        vend(1:2) = p%point_mid2%coor(1:2)

        mx = abs(ib-ia)
        my = abs(id-ic)
        
       kstart = p%nsub_outer 
        p%point_inter( p%nsub_outer)%coor( 1 ) = vorigin(1)
        p%point_inter( p%nsub_outer)%coor( 2 ) = vorigin(2)
        p%point_inter( p%nsub_outer)%ixy_type = 0
        p%point_inter( p%nsub_outer)%idx2 = 2*ia
        p%point_inter( p%nsub_outer)%idy2 = 2*ic
        !******
        p%nsub_outer = p%nsub_outer+ 1+ mx+my
        !******
        p%point_inter( p%nsub_outer)%coor( 1 ) = vend(1)
        p%point_inter( p%nsub_outer)%coor( 2 ) = vend(2)
        p%point_inter( p%nsub_outer)%ixy_type = 0
        p%point_inter( p%nsub_outer)%idx2 = 2*ib
        p%point_inter( p%nsub_outer)%idy2 = 2*id
        !**************************************************************************
 
        !call get_mono_inters(ia,ib,ic,id,mx,my,vorigin(1:2),vend(1:2),p,kstart)
        
        xi_o = p%point_mid1%xxii
        xi_e = p%point_mid2%xxii
        call get_mono_inters_QC(ia,ib,ic,id,mx,my,vorigin(1:2),vend(1:2),xi_o,xi_e,p,kstart,v11,v33)
 
        !***************************************************************************
        !***************************************************************************
        ia = p%point_mid2%id(1)
        ib = p%point_end%id(1)

        ic = p%point_mid2%id(2)
        id = p%point_end%id(2)
        vorigin(1:2) = p%point_mid2%coor(1:2)
        vend(1:2) = p%point_end%coor(1:2)

        mx = abs(ib-ia)
        my = abs(id-ic)
        
        kstart = p%nsub_outer
        p%point_inter( p%nsub_outer)%coor( 1 ) = vorigin(1)
        p%point_inter( p%nsub_outer)%coor( 2 ) = vorigin(2)
        p%point_inter( p%nsub_outer)%ixy_type = 0
        p%point_inter( p%nsub_outer)%idx2 = 2*ia
        p%point_inter( p%nsub_outer)%idy2 = 2*ic
        !******
        p%nsub_outer = p%nsub_outer+ 1+ mx+my
        !******
        p%point_inter( p%nsub_outer)%coor( 1 ) = vend(1)
        p%point_inter( p%nsub_outer)%coor( 2 ) = vend(2)
        p%point_inter( p%nsub_outer)%ixy_type = 0
        p%point_inter( p%nsub_outer)%idx2 = 2*ib
        p%point_inter( p%nsub_outer)%idy2 = 2*id
        !**************************************************************************
        ! we will got kstart+mx+my-1 - kstart +1 =mx+my intersections
        ! for example, kstart = 1 to kend = mx+my = kstart + mx + my - 1
 
        !call get_mono_inters(ia,ib,ic,id,mx,my,vorigin(1:2),vend(1:2),p,kstart)
        
        xi_o = p%point_mid2%xxii
        xi_e = 1.
        call get_mono_inters_QC(ia,ib,ic,id,mx,my,vorigin(1:2),vend(1:2),xi_o,xi_e,p,kstart,v11,v33)
 
    endif

    end subroutine get_intersections_qc
