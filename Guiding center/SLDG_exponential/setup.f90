    subroutine setup
    implicit none
    real :: alpha,beta,gamma
    real :: phi,sigma

    !*******************************************************************
    ! set up the number of cells
    !
    nx = 16*2**(kkkk-1)
    !nx = 4

    nx = 20*kkkk

    !nx = 160
    ny = nx
    ! we will set up the number of ghost cells later,
    ! since it is relevant to CFL.
    !***********************************

    nk = 2
    n_moment = (nk+1)*(nk+2)/2
    n_moment = 6

    n_moment_2d = 6
    n_moment_LDG = 10
    iqc = 1

    time_final = 1.

    if(kkkk==0)then
        cfl = 1.
    else
        cfl = (kkkk)*5.
    endif
    !cfl = 1.
    cfl = 1.

    nghost = int(cfl) +5

    !iexprk =1

    !iexample = 2
    irk = 3
    irelative = 1

    idebug = 0

    !if(iexample==1)then
    !    xleft = -pi
    !    xright = pi
    !    ybottom = -pi
    !    ytop = pi
    !elseif(iexample==2)then
    !    xleft = -pi
    !    xright = pi
    !    ybottom = -pi
    !    ytop = pi
    !elseif(iexample ==3)then
    !    xleft = -pi*2.
    !    xright = pi*2.
    !    ybottom = -pi*2.
    !    ytop = pi*2.
    !endif

    !*****************************************************
    !VP setup

    iexp_case = 7

    if( iexp_case == 1 )then
        bt_con(2,1) = 1.
        iadditive(0) = 0
        iexprk = 1
    elseif(iexp_case == 2)then
        bt_con(2,1) = 0.5
        iadditive(0) = 0
        bt_con(3,1) = 0. ;        bt_con(3,2) = 1.
        iadditive(1) = 0
        iexprk = 2
    elseif( iexp_case == 3 )then
        bt_con(2,1) = (2.-sqrt(2.) )/2.
        iadditive(0) = 0
        bt_con(3,1) = (-2.*sqrt(2.) )/3. ;      bt_con(3,2) = 1 - bt_con(3,1)
        iadditive(1) = 0
        bt_con(4,1) =   0. ;                    bt_con(4,2) = 1- bt_con(2,1) ; bt_con(4,3) = bt_con(2,1)
        iadditive(2) = 0
        iexprk = 3
    elseif( iexp_case == 4 )then
        bt_con(2,1) = 1./2.
        iadditive(0) = 0
        bt_con(3,1) = -1.     ;      bt_con(3,2) = 2.;
        iadditive(1) = 0
        bt_con(4,1) = 1./12.  ;      bt_con(4,2) = 1./3. ;   bt_con(4,3) = -1./4.;
        iadditive(2) = 0
        bt_con(5,1) = 1./12.  ;      bt_con(5,2) = 1./3. ;   bt_con(5,3) = 5./12.;
        iadditive(3) = 1
        iexprk = 4
    elseif( iexp_case == 5 )then
        gamma = (3.+sqrt(3.) )/6.
        phi = 1./(6.*(2.*gamma-1.) )
        bt_con(2,1) = gamma
        iadditive(0) = 0
        bt_con(3,1) = gamma-1.     ;      bt_con(3,2) = 2.*(1.-gamma) ;
        iadditive(1) = 0
        bt_con(4,1) = 0.  ;      bt_con(4,2) = 1./2.-phi ;   bt_con(4,3) = 1./2.+phi ;
        iadditive(2) = 0
        bt_con(5,1) = 0.  ;      bt_con(5,2) = phi ;   bt_con(5,3) = -phi;
        iadditive(3) = 1
        iexprk = 4
    elseif( iexp_case == 6 )then
        alpha = 0.5
        beta = 1./6.
        gamma = ( 3.+sqrt(3.) )/6.
        sigma = (  alpha+beta*(1-2.*gamma)-1./3.    )/( 1-2.*gamma )
        bt_con(2,1) = gamma
        iadditive(0) = 0
        bt_con(3,1) = gamma-1.     ;      bt_con(3,2) = 2.*(1.-gamma) ;
        iadditive(1) = 0
        bt_con(4,1) = alpha  ;      bt_con(4,2) = beta ;   bt_con(4,3) = sigma ;
        iadditive(2) = 0
        bt_con(5,1) = -alpha  ;      bt_con(5,2) = 0.5-beta ;   bt_con(5,3) = 0.5-sigma ;
        iadditive(3) = 1
        iexprk = 4
    elseif(iexp_case == 7)then
        bt_con(2,1) = 1./3.
        iadditive(0) = 0
        bt_con(3,1) = 0. ;        bt_con(3,2) = 2./3.
        iadditive(1) = 0
        bt_con(4,1) = -1./12. ;   bt_con(4,2) = 0. ;      bt_con(4,3) = 3./4.
        iadditive(2) = 0
        iexprk = 3
    endif

    kdg = 3
    call SLDG_poisson1d_paras

    ! write(*,*) 'i_case: '
    i_case = 1
    !*******************************************************************
    ! currently, we only test the code for
    ! a state and smooth steady incompressible Euler equations.

    !*******************************************************************
    ! compatational domain
    xleft = 0.
    xright = 2.*pi
    ybottom = 0.
    ytop = 2.*pi


    !***************************************************
    if( n_moment==3 )then
        gau(1,1)=-sqrt(3.0)/6.0
        gau(2,1)=sqrt(3.0)/6.0
        gau(1,2)=0.5
        gau(2,2)=0.5
    elseif( n_moment==6 )then
        gau(1,1) = -sqrt(15.)/10.
        gau(2,1) = 0.
        gau(3,1) = sqrt(15.)/10.
        gau(1,2) = 5./18.
        gau(2,2) =  4./9.
        gau(3,2) = 5./18.
    endif

    end subroutine setup
    !*******************************************************************
    ! velocity fields for linear transport problem
    real function  vel_x( x,y,t )
    implicit none
    real,intent(in) :: x,y,t

    if( iexample == 1 )then
        vel_x = 1.
    elseif( iexample == 2 )then
        vel_x = -cos(x/2)**2 * sin(y) * cos(pi*t/time_final )*pi
    elseif( iexample == 3 )then
        vel_x = -y
    endif

    end function vel_x

    real function  vel_y( x,y,t )
    implicit none

    real,intent(in) :: x,y,t

    if( iexample == 1 )then
        vel_y = 1.
    elseif( iexample == 2 )then
        vel_y = sin(x)*cos(y/2.)**2 * cos(pi*t/time_final )*pi
    elseif(iexample == 3)then
        vel_y = x
    endif

    end function vel_y