    subroutine SLDG1D_QiuGuoCai(    ixy,dxy, tnum,dt17,euler_gl, euler_gauss, xy_split,umod_temp )
    implicit none
    integer,intent(in) :: ixy
    real,intent(in) :: dxy,dt17,tnum
    real,intent(in) :: euler_gl(1:n_gl)
    real,intent(in) :: euler_gauss(1:n_g)
    real,intent(in) :: xy_split
    real :: upstream_gl(1:n_gl), upstream_gauss(1:n_g)
    real,intent(out) :: umod_temp(1:n_moment )
    real :: xl_star,xr_star
    integer :: ia,ib
    integer :: mx
    integer :: kc,inter
    integer :: isub , nsub
    integer :: idx
    real,allocatable :: zz(:)

    integer :: ik
    real :: temp_e,xrg,yrg
    real :: temp_e1,temp_e2

    ! STEP1. Locate the foots of trajectory x_l^{n,star}, x_r^{n,star}.
    ! We numerically solve the following final-value problem (trajectory equation):
    !             d x(t) / dt = a( x(t) , t )
    ! with the final-value x( t^{n+1} ) = x_l^n, x_r^n by a high order numerical integrator such as a classical fourth-order
    ! numerical integrator such as a classical fourth-order Runge-Kutta method.


    ! STEP2. locate Gauss-lobatto points

    if( ixy == 1 )then
        xrg = x(ie) + x_GL( 1 )*dx
        temp_e1 = polynomial2d(  phi_y( ie-1, je, 1:n_moment_LDG ),xrg,x(ie-1),dx,xy_split,y(je),dy,n_moment_LDG )
        temp_e2 = polynomial2d(  phi_y( ie, je, 1:n_moment_LDG ),xrg,x(ie),dx,xy_split,y(je),dy,n_moment_LDG )
        upstream_gl(1) = euler_gl(1) + sign_rho* (temp_e1+temp_e2)*0.5 * dt17

        xrg = x(ie) + x_GL( 2 )*dx
        temp_e1 = polynomial2d(  phi_y( ie, je, 1:n_moment_LDG ),xrg,x(ie),dx,xy_split,y(je),dy,n_moment_LDG )
        temp_e2 = polynomial2d(  phi_y( ie+1, je, 1:n_moment_LDG ),xrg,x(ie+1),dx,xy_split,y(je),dy,n_moment_LDG )
        upstream_gl(2) = euler_gl(2) + sign_rho* (temp_e1+temp_e2)*0.5 * dt17

        do ik = 1,n_g
            xrg = x(ie) + x_G( ik )*dx
            temp_e = polynomial2d(  phi_y( ie, je, 1:n_moment_LDG ),xrg,x(ie),dx,xy_split,y(je),dy,n_moment_LDG )
            upstream_gauss(ik) = euler_gauss(ik) + sign_rho* temp_e * dt17
        enddo

    elseif(ixy == 2)then
        yrg = y(je) + x_GL( 1 )*dy
        temp_e1 = polynomial2d(  phi_x( ie, je-1, 1:n_moment_LDG ),xy_split,x(ie),dx,yrg,y(je-1),dy,n_moment_LDG )
        temp_e2 = polynomial2d(  phi_x( ie, je, 1:n_moment_LDG ),xy_split,x(ie),dx,yrg,y(je),dy,n_moment_LDG )
        upstream_gl(1) = euler_gl(1) - sign_rho* (temp_e1+temp_e2)*0.5 * dt17

        yrg = y(je) + x_GL( 2 )*dy
        temp_e1 = polynomial2d(  phi_x( ie, je, 1:n_moment_LDG ),xy_split,x(ie),dx,yrg,y(je),dy,n_moment_LDG )
        temp_e2 = polynomial2d(  phi_x( ie, je+1, 1:n_moment_LDG ),xy_split,x(ie),dx,yrg,y(je+1),dy,n_moment_LDG )
        upstream_gl(2) = euler_gl(2) - sign_rho* (temp_e1+temp_e2)*0.5 * dt17

        do ik = 1,n_g
            yrg = y(je) + x_G( ik )*dy
            temp_e = polynomial2d(  phi_x( ie, je, 1:n_moment_LDG ),xy_split,x(ie),dx,yrg,y(je),dy,n_moment_LDG )
            upstream_gauss(ik) = euler_gauss(ik) - sign_rho* temp_e * dt17
        enddo
    endif


    !
    xl_star =  upstream_gl( 1 )
    xr_star =  upstream_gl( n_gl )


    call id_get17( xl_star,dxy,ia,ixy )
    call id_get17( xr_star,dxy,ib,ixy )


    mx = ib - ia

    allocate( zz(1:2+mx) )

    zz(1) = xl_star
    zz(2+mx ) = xr_star

    if(mx .ne. 0)then
        do kc = 1 , mx
            inter = ia + kc
            if( ixy == 1 )then
                zz( 1+kc ) = xgrid( inter )
            elseif(ixy == 2)then
                zz( 1+kc ) = ygrid(inter)
            endif
        enddo
    endif

    do isub = 1 , 1 + mx
        Dn_star_sub(isub)%vpoint(1) = zz(isub)
        Dn_star_sub(isub)%vpoint(2) = zz(isub+1)

        call id_get17( (zz(isub)+zz(isub+1) )/2.,dxy,idx ,ixy )
       
        Dn_star_sub(isub)%id =  idx

    enddo

    deallocate(zz)
    nsub = 1 + mx

    if(nod==1)then
        call Green2( upstream_gl(1:n_gl ), upstream_gauss(1:n_g),umod_temp(1:n_moment ) , nsub, dxy,ixy)
    elseif(nod==2)then
        !call Green3( umod_temp(1:n_gl ) )
        call Green3( upstream_gl(1:n_gl ), upstream_gauss(1:n_g), umod_temp(1:n_moment ) , nsub, dxy,ixy)
    elseif(nod==3)then
        !call Green4( umod_temp(1:n_gl ) )
    endif


    end subroutine SLDG1D_QiuGuoCai