
    module poisson_periodic_LDG_mod
    !************************************************************
    !  The first version at March 28th, at UD
    !  This version solve 
    !      Laplace psi = sign_rho * rho
    !  with periodic boundary condition,
    !  where sign_rho is sign for positive or negative
    !  
    !  If sign_rho = 1, for solving incompressible Euler equations.
    !  If sign_rho = -1, for solving Vlasov Guiding center.
    !
    !  n_moment_2d and n_moment_LDG is for  LDG with various 
    ! polynomials.
    ! if n_moment_2d =3 and n_moment_LDG =6,
    ! it means input P1 solution, and output p2 solution.
    !*************************************************************
    ! please set sign_rho, n_moment_2d and n_moment_LDG in
    ! subroutine poisson_parameters
    !************************************************************
    use polynomial_mod
    use poisson_system_solver
    use variable , only : nel_x,nel_y,dx,dy,nt
    use poisson_element_mod
    implicit none
    ! for LDG solver of 2D Poisson equation
    !*************************************************************
    ! temp11(:,:,:) is for
    !  INPUTing the DG solution for the right-hand side
    !   temp11( :     ,    :   ,    :   )
    !           |          |        |
    !    index of x        |        |
    !                  index of y   |
    !                           DG moments of element (i,j)
    real,allocatable,public :: temp11(:,:,:)
    !*************************************************************
    ! phi_x(:,:,:) is for
    !  OUTPUTing the DG solution for the velocity feild
    !   phi_x( :     ,     :   ,    :   )
    !           |          |        |
    !    index of x        |        |
    !                  index of y   |
    !                           DG moments of element (i,j)
    real,allocatable,public :: phi_x(:,:,:),phi_y(:,:,:)
    !**************************************************************
    integer,public :: n_moment_2d  ! for inputing temp11
    integer,public :: n_moment_LDG ! for outputing phi_x,phi_y
    !**************************************************
    real,public :: sign_rho
    !**************************************************************
    real,public :: ai2d(10) ! for the transfer from nodal to modal
    !*********************************************************************
    !set them in subroutine parameters
    real,private :: gauss(5,2)
    real,private :: c(10)
    !*************************************************************
    real,private :: block_a1(10,10),block_a1w(10,10)
    real,private :: block_a2(10,10),block_a2s(10,10)

    real,private :: block_b1(10,10),block_b1e(10,10)
    real,private :: block_b2(10,10),block_b2n(10,10)

    real,private :: block_c(10,10)
    real,private :: block_ce(10,10),block_cw(10,10),block_cn(10,10),block_cs(10,10)

    real,private :: block_dp( 10,10 )
    real,private :: block_dn( 10,10 )
    real,private :: block_de( 10,10 )
    real,private :: block_dw( 10,10 )
    real,private :: block_ds( 10,10 )
    ! Compressed Sparse Row formate, used in calling pmgmres
    integer,private :: na !
    real,allocatable,private :: aa(:)
    integer,allocatable,private :: ia(:),ja(:)
    ! the row and column indices of the matrix values, used in calling pmgmres
    real,allocatable,private :: bb(:),ss(:),ss_phi(:)
    integer,private :: nz_num
    ! the number of nonzero matrix values, used in calling pmgmres
    integer,private :: iele

    integer,private :: iter
    real,private :: err
    ! Compressed Sparse Row formate
    !*********************************************************************
    real,allocatable,private :: aa1(:)
    integer,allocatable,private :: ia1(:),ja1(:)
    real,allocatable,private :: aa2(:)
    integer,allocatable,private :: ia2(:),ja2(:)
    integer,private :: nz_num1
    !*********************************************************************
    integer,allocatable,private :: inx(:),iny(:)
    ! for LDG solver of 2D Poisson equation
    !*********************************************************************

    CONTAINS

    subroutine allocate_var_poisson
    implicit none
    !************************************************************************
    ! for calling LDG
    allocate( temp11(1-1:nel_x+1,1-1:nel_y+1,1:n_moment_2d) )
    allocate( phi_x(1-1:nel_x+1,1-1 :nel_y+1 ,1:n_moment_LDG) )
    allocate( phi_y(1-1 :nel_x+1 ,1-1 :nel_y+1 ,1:n_moment_LDG) )
    allocate( element_poisson(1:nel_x*nel_y) )
    !***************used in subroutine poisson_rhs_vector.f90
    ! compute the right hand side vector of the linear system
    allocate( bb( 1:nel_x*nel_y*n_moment_LDG ) , ss( 1:nel_x*nel_y*n_moment_LDG ),ss_phi(1:nel_x*nel_y*n_moment_LDG) )
    allocate( inx(1:nel_x*nel_y) ,iny(1:nel_x*nel_y) )

    !***************used in subroutine poisson_rhs_vector.f90
    !************************************************************************
    end subroutine allocate_var_poisson

    subroutine deallocate_var_poisson
    implicit none
    !************************************************************************
    ! for calling LDG
    deallocate( temp11 )
    deallocate( phi_x )
    deallocate( phi_y )
    deallocate( element_poisson  )
    !***************used in subroutine poisson_rhs_vector.f90
    ! compute the right hand side vector of the linear system
    deallocate( bb , ss ,ss_phi  )
    deallocate( inx  ,iny  )

    !*******************************************
    deallocate( aa )
    deallocate( ja )
    deallocate( ia )
    ! compute the right hand side vector of the linear system

    !*******************
    deallocate( aa1 )
    deallocate( ja1 )
    deallocate( ia1 )
    !*******************
    deallocate( aa2 )
    deallocate( ja2 )
    deallocate( ia2 )
    !*******************************************
    !***************used in subroutine poisson_rhs_vector.f90
    !************************************************************************
    end subroutine deallocate_var_poisson


    subroutine poisson_parameters
    implicit none

    n_moment_2d = 3
    n_moment_LDG = 6

    ! the points of 10th order Gauss quadrature
    gauss(1,1) = -sqrt( 5.+2.*sqrt(10./7.) )/6.
    gauss(2,1) = -sqrt( 5.-2.*sqrt(10./7.) )/6.
    gauss(3,1) = 0.
    gauss(4,1) = sqrt( 5.-2.*sqrt(10./7.) )/6.
    gauss(5,1) = sqrt( 5.+2.*sqrt(10./7.) )/6.
    ! coefficients of 10th order Gauss quadrature
    gauss(1,2) =  ( 322.-13.*sqrt(70.) )/900. *0.5
    gauss(2,2) =  ( 322.+13.*sqrt(70.) )/900. *0.5
    gauss(3,2) =  128./225.*0.5
    gauss(4,2) =  ( 322.+13.*sqrt(70.) )/900. *0.5
    gauss(5,2) =  ( 322.-13.*sqrt(70.) )/900. *0.5

    !*************************
    !LDG
    if(n_moment_LDG ==1)then
        c(1) = 1.
    elseif(n_moment_LDG ==3)then
        c(1) = 1.
        c(2) = 1./12.
        c(3) = 1./12.
    elseif(n_moment_LDG == 6)then
        c(1) = 1.
        c(2) = 1./12.
        c(3) = 1./12.
        c(4) = 1./180.
        c(5) = 1./144.
        c(6) = 1./180.
    elseif(n_moment_LDG == 10)then
        c(1) = 1.
        c(2) = 1./12.
        c(3) = 1./12.
        c(4) = 1./180.
        c(5) = 1./144.
        c(6) = 1./180.
        c(7) = 1./2800.
        c(8) = 1./2160.
        c(9) = 1./2160.
        c(10) = 1./2800.
    endif

    ai2d(1) = 1.
    ai2d(2) = 12.
    ai2d(3) = 12.
    ai2d(4) = 180.
    ai2d(5) = 144.
    ai2d(6) = 180.    
    ai2d(7) = 2800.
    ai2d(8) = 2160.
    ai2d(9) = 2160.
    ai2d(10) = 2800.    
    
    end subroutine poisson_parameters



    subroutine poisson_coefficients
    implicit none
    integer :: l,k,mx,my,ii,jj
    real :: temp
    real :: block_temp(1:n_moment_LDG,1:n_moment_LDG)
    real :: block_b1a1(1:n_moment_LDG,1:n_moment_LDG)
    real :: block_b1ea1w(1:n_moment_LDG,1:n_moment_LDG)
    real :: block_b2a2(1:n_moment_LDG,1:n_moment_LDG)
    real :: block_b2na2s(1:n_moment_LDG,1:n_moment_LDG)
    real :: block_b2na2(1:n_moment_LDG,1:n_moment_LDG)
    real :: block_b1ea1(1:n_moment_LDG,1:n_moment_LDG)
    real :: block_b1a1w(1:n_moment_LDG,1:n_moment_LDG)
    real :: block_b2a2s(1:n_moment_LDG,1:n_moment_LDG)


    do l= 1,n_moment_LDG
        do k= 1,n_moment_LDG
            temp=0.0
            do my=1,5
                temp=temp+gauss(my,2)*fle_comp( k-1 ,0.5,gauss(my,1))*fle_comp( L-1 ,0.5,gauss(my,1))
            enddo
            do mx = 1,5
                do my = 1,5
                    temp = temp - gauss(mx,2)*gauss(my,2)*fle_comp( k-1 ,gauss(mx,1),gauss(my,1))*flex( L-1 ,gauss(mx,1),gauss(my,1))
                enddo
            enddo

            block_a1(L,k)=temp/dx/c(L)
            !***********************************************************************************************
            temp=0.0
            do my=1,5
                temp=temp+gauss(my,2)*fle_comp( k-1 ,0.5,gauss(my,1)  )*fle_comp( L-1 ,-0.5,gauss(my,1))
            enddo
            block_a1w(L,k) = -temp/dx/c(L)
            !***********************************************************************************************
            temp=0.0
            do mx=1,5
                temp=temp+gauss(mx,2)*fle_comp( k-1 ,gauss(mx,1),0.5 )*fle_comp( l-1 , gauss(mx,1),0.5 )
            enddo
            do mx = 1,5
                do my = 1,5
                    temp = temp - gauss(mx,2)*gauss(my,2)*fle_comp( k-1 ,gauss(mx,1),gauss(my,1))*fley( L-1 ,gauss(mx,1),gauss(my,1))
                enddo
            enddo

            block_a2(l,k)=temp/dy/c(l)
            !***********************************************************************************************
            temp=0.0
            do mx=1,5
                temp=temp+gauss(mx,2)*fle_comp( k-1 , gauss(mx,1),0.5 )*fle_comp( l-1 ,gauss(mx,1) ,-0.5 )
            enddo
            block_a2s(l,k) = -temp/dy/c(L)
            !***********************************************************************************************
            !***********************************************************************************************
            temp = 0.
            do my = 1 , 5
                temp = temp + gauss(my,2)*fle_comp(k-1,-0.5,gauss(my,1) )*fle_comp(L-1,-0.5,gauss(my,1) )
            enddo
            do mx = 1, 5
                do my = 1,5
                    temp = temp + gauss(mx,2)*gauss(my,2)*fle_comp(k-1,gauss(mx,1),gauss(my,1))*flex(L-1,gauss(mx,1),gauss(my,1))
                enddo
            enddo

            block_b1(L,k) = - temp/dx
            !***********************************************************************************************
            temp = 0.
            do my = 1 , 5
                temp = temp + gauss(my,2)*fle_comp(k-1,-0.5,gauss(my,1) )*fle_comp(L-1,0.5,gauss(my,1) )
            enddo
            block_b1e(L,k) = temp/dx
            !***********************************************************************************************
            temp = 0.
            do mx = 1 , 5
                temp = temp + gauss(mx,2)*fle_comp(k-1,gauss(mx,1),-0.5 )*fle_comp(L-1,gauss(mx,1),-0.5 )
            enddo
            do mx = 1 , 5
                do my = 1 , 5
                    temp = temp + gauss(mx,2)*gauss(my,2)*fle_comp(k-1,gauss(mx,1),gauss(my,1) )*fley( L-1,gauss(mx,1),gauss(my,1) )
                enddo
            enddo
            block_b2(L,k) = -temp/dy
            !***********************************************************************************************
            temp = 0.
            do mx = 1 , 5
                temp = temp + gauss(mx,2)*fle_comp(k-1,gauss(mx,1) ,-0.5)*fle_comp(L-1,gauss(mx,1),0.5)
            enddo
            block_b2n(L,k) = temp/dy
            !***********************************************************************************************
            !***********************************************************************************************
            temp = 0.0
            do my = 1 , 5
                temp = temp + gauss(my,2) * fle_comp( k-1,0.5 ,gauss(my,1) ) * fle_comp( l-1,0.5,gauss(my,1) ) &
                    + gauss(my,2) * fle_comp( k-1,-0.5 ,gauss(my,1) ) * fle_comp( l-1,-0.5,gauss(my,1) )
            enddo
            temp = temp * dy/dx
            do mx = 1 , 5
                temp = temp + gauss(mx,2) * fle_comp( k-1,gauss(mx,1) ,0.5 ) * fle_comp( L-1,gauss(mx,1),0.5 ) &
                    + gauss(mx,2) * fle_comp( k-1,gauss(mx,1) ,-0.5 ) * fle_comp( L-1,gauss(mx,1),-0.5  )

            enddo
            block_c(L,k) = -temp/dy


            !************************************************************************************************
            temp = 0.0
            do my = 1 , 5
                temp = temp + gauss(my,2)*fle_comp( k-1,-0.5,gauss(my,1)  )*fle_comp( L-1,0.5,gauss(my,1) )
            enddo
            block_ce(l,k) = temp/dx
            !************************************************************************************************
            temp = 0.0
            do my = 1 , 5
                temp = temp + gauss(my,2)*fle_comp( k-1,0.5,gauss(my,1)  )*fle_comp(L-1,-0.5,gauss(my,1) )
            enddo
            block_cw(l,k) = temp/dx
            !************************************************************************************************
            temp = 0.
            do mx = 1 , 5
                temp = temp + gauss(mx,2)*fle_comp(k-1,gauss(mx,1) ,-0.5)*fle_comp(L-1,gauss(mx,1),0.5 )
            enddo
            block_cn(l,k) = temp/dy
            !************************************************************************************************
            temp = 0.0
            do mx = 1,5
                temp = temp + gauss(mx,2)*fle_comp(k-1,gauss(mx,1) ,0.5 )*fle_comp(L-1,gauss(mx,1),-0.5 )
            enddo
            block_cs(L,k) = temp/dy

        enddo
    enddo

    !#########################################################################################################
    block_dp(1:n_moment_LDG,1:n_moment_LDG) = block_c(1:n_moment_LDG,1:n_moment_LDG)
    call blocks_multiply( n_moment_LDG,block_b1(1:n_moment_LDG,1:n_moment_LDG),block_a1(1:n_moment_LDG,1:n_moment_LDG),block_b1a1(1:n_moment_LDG,1:n_moment_LDG) )
    call blocks_multiply( n_moment_LDG,block_b1e(1:n_moment_LDG,1:n_moment_LDG),block_a1w(1:n_moment_LDG,1:n_moment_LDG),block_b1ea1w(1:n_moment_LDG,1:n_moment_LDG) )
    call blocks_multiply( n_moment_LDG,block_b2(1:n_moment_LDG,1:n_moment_LDG),block_a2(1:n_moment_LDG,1:n_moment_LDG),block_b2a2(1:n_moment_LDG,1:n_moment_LDG) )
    call blocks_multiply( n_moment_LDG,block_b2n(1:n_moment_LDG,1:n_moment_LDG),block_a2s(1:n_moment_LDG,1:n_moment_LDG),block_b2na2s(1:n_moment_LDG,1:n_moment_LDG) )

    call blocks_multiply( n_moment_LDG,block_b2n(1:n_moment_LDG,1:n_moment_LDG),block_a2(1:n_moment_LDG,1:n_moment_LDG),block_b2na2(1:n_moment_LDG,1:n_moment_LDG) )
    call blocks_multiply( n_moment_LDG,block_b1e(1:n_moment_LDG,1:n_moment_LDG),block_a1(1:n_moment_LDG,1:n_moment_LDG),block_b1ea1(1:n_moment_LDG,1:n_moment_LDG) )
    call blocks_multiply( n_moment_LDG,block_b1(1:n_moment_LDG,1:n_moment_LDG),block_a1w(1:n_moment_LDG,1:n_moment_LDG),block_b1a1w(1:n_moment_LDG,1:n_moment_LDG) )
    call blocks_multiply( n_moment_LDG,block_b2(1:n_moment_LDG,1:n_moment_LDG),block_a2s(1:n_moment_LDG,1:n_moment_LDG),block_b2a2s(1:n_moment_LDG,1:n_moment_LDG) )

    do ii = 1 ,n_moment_LDG
        do jj = 1 ,n_moment_LDG
            !******************************************************************
            ! for phi_xx + phi_yy = rho
            ! phi_x = -E_1
            ! phi_y = -E_2
            ! (E_1)_x + (E_2)_x = rho
            !block_dp(ii,jj) = (block_c(ii,jj)+( block_b1a1(ii,jj)+block_b1ea1w(ii,jj)+block_b2a2(ii,jj)+block_b2na2s(ii,jj) ) )*dx*dy
            !
            !block_dn(ii,jj) = (block_cn(ii,jj) + block_b2na2(ii,jj) )*dx*dy
            !block_de(ii,jj) = (block_ce(ii,jj) + block_b1ea1(ii,jj) )*dx*dy
            !block_dw(ii,jj) = (block_cw(ii,jj) + block_b1a1w(ii,jj) )*dx*dy
            !block_ds(ii,jj) = (block_cs(ii,jj) + block_b2a2s(ii,jj) )*dx*dy
            !*******************************************************************
            ! for phi_xx + phi_yy = - rho
            ! phi_x = -E_1
            ! phi_y = -E_2
            ! (E_1)_x + (E_2)_x = rho
            block_dp(ii,jj) = (block_c(ii,jj) + sign_rho*( block_b1a1(ii,jj)+block_b1ea1w(ii,jj)+block_b2a2(ii,jj)+block_b2na2s(ii,jj) ) )*dx*dy


            block_dn(ii,jj) = (block_cn(ii,jj) + sign_rho* block_b2na2(ii,jj) )*dx*dy
            block_de(ii,jj) = (block_ce(ii,jj) + sign_rho* block_b1ea1(ii,jj) )*dx*dy
            block_dw(ii,jj) = (block_cw(ii,jj) + sign_rho* block_b1a1w(ii,jj) )*dx*dy
            block_ds(ii,jj) = (block_cs(ii,jj) + sign_rho* block_b2a2s(ii,jj) )*dx*dy

        enddo
    enddo

    end subroutine poisson_coefficients
    !***************************************************************************************************
    subroutine blocks_multiply( nn,ann,bnn,cnn )
    implicit none
    integer, intent(in) :: nn
    real,intent(in) :: ann(nn,nn),bnn( nn,nn )
    real,intent(out) :: cnn(nn,nn)
    integer :: ii,jj,kk
    real :: temp

    do ii = 1, nn
        do jj = 1 ,nn
            temp = 0.
            do kk = 1,nn
                temp = temp + ann(ii,kk)*bnn(kk,jj)

            enddo
            cnn(ii,jj) = temp
        enddo
    enddo


    end subroutine blocks_multiply
    !*****************************************************
    subroutine blocks_minus( nn,ann,bnn,cnn )
    implicit none
    integer, intent(in) :: nn
    real,intent(in) :: ann(nn,nn),bnn( nn,nn )
    real,intent(out) :: cnn(nn,nn)
    integer :: ii,jj



    do ii = 1, nn
        do jj = 1 ,nn
            cnn(ii,jj) = ann(ii,jj) - bnn(ii,jj)
        enddo
    enddo


    end subroutine blocks_minus



    subroutine poisson_CSR_init
    implicit none
    ! to get element_poisson(iele) % inside_or_boundary

    do iele = 1 , nel_x*nel_y
        if( iele == 1 )then
            element_poisson(iele)%inside_or_boundary = 5
        elseif(iele == nel_x)then
            element_poisson(iele)%inside_or_boundary = 6
        elseif(iele == nel_x*nel_y)then
            element_poisson(iele)%inside_or_boundary = 7
        elseif(iele == nel_x*(nel_y-1)+1  )then
            element_poisson(iele)%inside_or_boundary = 8
        elseif( iele>1 .and. iele<nel_x )then
            element_poisson(iele)%inside_or_boundary = 1
        elseif( iele>nel_x*(nel_y-1)+1 .and. iele<nel_x*nel_y )then
            element_poisson(iele)%inside_or_boundary = 3
        elseif( iele/nel_x*nel_x== iele .and. iele > nel_x .and. iele <nel_x*nel_y )then
            element_poisson(iele)%inside_or_boundary = 2
        elseif( iele/nel_x*nel_x+1 ==iele .and. iele >1 .and. iele<nel_x*(nel_y-1) + 1 )then
            element_poisson(iele)%inside_or_boundary = 4
        else
            element_poisson(iele)%inside_or_boundary = 0
        endif
    enddo
    !******************************
    do iele = 1 , nel_x*nel_y
        select case(element_poisson(iele)%inside_or_boundary )
        case(0)
            element_poisson(iele)%pivot%irank = 3
            element_poisson(iele)%east%irank = 4
            element_poisson(iele)%west%irank = 2
            element_poisson(iele)%north%irank = 5
            element_poisson(iele)%south%irank = 1

            element_poisson(iele)%pivot%icolumn = iele
            element_poisson(iele)%east%icolumn = iele + 1
            element_poisson(iele)%west%icolumn = iele - 1
            element_poisson(iele)%north%icolumn = iele + nel_x
            element_poisson(iele)%south%icolumn = iele - nel_x
        case(1)
            element_poisson(iele)%pivot%irank = 2
            element_poisson(iele)%east%irank = 3
            element_poisson(iele)%west%irank = 1
            element_poisson(iele)%north%irank = 4
            element_poisson(iele)%south%irank = 5

            element_poisson(iele)%pivot%icolumn = iele
            element_poisson(iele)%east%icolumn = iele + 1
            element_poisson(iele)%west%icolumn = iele - 1
            element_poisson(iele)%north%icolumn = iele + nel_x
            element_poisson(iele)%south%icolumn = iele + nel_x*(nel_y-1)
        case(2)
            element_poisson(iele)%pivot%irank = 4
            element_poisson(iele)%east%irank = 2
            element_poisson(iele)%west%irank = 3
            element_poisson(iele)%north%irank = 5
            element_poisson(iele)%south%irank = 1

            element_poisson(iele)%pivot%icolumn = iele
            element_poisson(iele)%east%icolumn = nel_x*( iele/nel_x - 1 ) + 1
            element_poisson(iele)%west%icolumn = iele - 1
            element_poisson(iele)%north%icolumn = iele + nel_x
            element_poisson(iele)%south%icolumn = iele - nel_x
        case(3)
            element_poisson(iele)%pivot%irank = 4
            element_poisson(iele)%east%irank = 5
            element_poisson(iele)%west%irank = 3
            element_poisson(iele)%north%irank = 1
            element_poisson(iele)%south%irank = 2

            element_poisson(iele)%pivot%icolumn = iele
            element_poisson(iele)%east%icolumn = iele + 1
            element_poisson(iele)%west%icolumn = iele - 1
            element_poisson(iele)%north%icolumn = iele - nel_x*(nel_y-1)
            element_poisson(iele)%south%icolumn = iele - nel_x
        case(4)
            element_poisson(iele)%pivot%irank = 2
            element_poisson(iele)%east%irank = 3
            element_poisson(iele)%west%irank = 4
            element_poisson(iele)%north%irank = 5
            element_poisson(iele)%south%irank = 1

            element_poisson(iele)%pivot%icolumn = iele
            element_poisson(iele)%east%icolumn = iele + 1
            element_poisson(iele)%west%icolumn = ( iele/nel_x + 1 )*nel_x
            element_poisson(iele)%north%icolumn = iele + nel_x
            element_poisson(iele)%south%icolumn = iele - nel_x
        case(5)
            element_poisson(iele)%pivot%irank = 1
            element_poisson(iele)%east%irank = 2
            element_poisson(iele)%west%irank = 3
            element_poisson(iele)%north%irank = 4
            element_poisson(iele)%south%irank = 5

            element_poisson(iele)%pivot%icolumn = iele
            element_poisson(iele)%east%icolumn = iele + 1
            element_poisson(iele)%west%icolumn = nel_x
            element_poisson(iele)%north%icolumn = iele + nel_x
            element_poisson(iele)%south%icolumn = iele + nel_x*(nel_y-1)
        case(6)
            element_poisson(iele)%pivot%irank = 3
            element_poisson(iele)%east%irank = 1
            element_poisson(iele)%west%irank = 2
            element_poisson(iele)%north%irank = 4
            element_poisson(iele)%south%irank = 5

            element_poisson(iele)%pivot%icolumn = iele
            element_poisson(iele)%east%icolumn = 1
            element_poisson(iele)%west%icolumn = iele - 1
            element_poisson(iele)%north%icolumn = 2*nel_x
            element_poisson(iele)%south%icolumn = nel_x*nel_y
        case(7)
            element_poisson(iele)%pivot%irank = 5
            element_poisson(iele)%east%irank = 3
            element_poisson(iele)%west%irank = 4
            element_poisson(iele)%north%irank = 1
            element_poisson(iele)%south%irank = 2

            element_poisson(iele)%pivot%icolumn = iele
            element_poisson(iele)%east%icolumn = nel_x*(nel_y-1) + 1
            element_poisson(iele)%west%icolumn = iele - 1
            element_poisson(iele)%north%icolumn = nel_x
            element_poisson(iele)%south%icolumn = nel_x*(nel_y-1)
        case(8)
            element_poisson(iele)%pivot%irank = 3
            element_poisson(iele)%east%irank = 4
            element_poisson(iele)%west%irank = 5
            element_poisson(iele)%north%irank = 1
            element_poisson(iele)%south%irank = 2

            element_poisson(iele)%pivot%icolumn = iele
            element_poisson(iele)%east%icolumn = iele + 1
            element_poisson(iele)%west%icolumn = nel_x*(nel_y)
            element_poisson(iele)%north%icolumn = 1
            element_poisson(iele)%south%icolumn = iele - nel_x
        end select
    enddo

    end subroutine poisson_CSR_init


    subroutine poisson_CSR_format
    implicit none
    integer :: ii,kk
    integer :: nz_num_subrow1,nz_num_row1
    integer :: ip_start,ie_start,iw_start,in_start,is_start

    real :: ai1
    ! store in the Compressed Sparse Row format

    ai1 = dx*dy



    na = nel_x*nel_y*n_moment_LDG
    nz_num = n_moment_LDG*n_moment_LDG * 5*nel_x*nel_y + nel_x*nel_y -5
    allocate( aa(nz_num),ja(nz_num) )
    allocate( ia(na+1) )

    nz_num_subrow1 = n_moment_LDG*5 + nel_x*nel_y - 5


    ia(1) = 1
    ia(2) = 1 + nz_num_subrow1

    do ii = 2,na
        ia(ii+1) = ia(ii) + 5*n_moment_LDG
    enddo

    aa( 2*n_moment_LDG+1:2*n_moment_LDG+nel_x-3 ) = ai1


    do ii = 1 , nel_x - 3
        ja( 2*n_moment_LDG+ii ) = (ii+1) * n_moment_LDG + 1
    enddo

    aa( 4*n_moment_LDG+nel_x-3  +1 : 4*n_moment_LDG+nel_x-3 + nel_x*(nel_y-2)-1 ) = ai1



    do ii = 1 , nel_x*(nel_y-2)-1
        ja( 4*n_moment_LDG+nel_x-3 +ii ) = (nel_x+ii) * n_moment_LDG + 1
    enddo

    aa( 5*n_moment_LDG + nel_x*(nel_y-1)-4 + 1 : nz_num_subrow1 ) = ai1


    do ii = 1 , nel_x - 1
        ja( 5*n_moment_LDG+nel_x*(nel_y-1)-4 + ii ) = ( nel_x*(nel_y-1) + ii ) * n_moment_LDG + 1
    enddo

    ! compute element_poisson(1)

    aa(1) = block_dp(1,1) + ai1

    aa(2:n_moment_LDG) = block_dp(1,2:n_moment_LDG)

    aa(1+n_moment_LDG) = block_de(1,1) + ai1
    aa(2+n_moment_LDG:2*n_moment_LDG ) = block_de(1,2:n_moment_LDG)

    do ii = 1 , 2*n_moment_LDG
        ja( ii) = ii
    enddo

    iw_start = 2*n_moment_LDG + nel_x - 3
    aa(iw_start + 1) = block_dw(1,1) + ai1
    aa(iw_start + 2: iw_start+n_moment_LDG ) = block_dw(1,2:n_moment_LDG)

    aa(iw_start + n_moment_LDG +1 ) = block_dn(1,1) + ai1
    aa(iw_start + n_moment_LDG +2 : iw_start+2*n_moment_LDG ) = block_dn(1,2:n_moment_LDG)

    do ii = 1 , 2*n_moment_LDG
        ja(iw_start+ii) = (nel_x-1)*n_moment_LDG + ii
    enddo

    is_start = 4*n_moment_LDG + nel_x*(nel_y-1) - 4

    aa(is_start+1) = block_ds(1,1) + ai1
    aa(is_start+2 : is_start+n_moment_LDG ) = block_ds(1,2:n_moment_LDG)

    do ii = 1,n_moment_LDG
        ja(is_start+ii) = nel_x*(nel_y-1)*n_moment_LDG + ii
    enddo

    !*******************
    do kk = 2,n_moment_LDG
        ip_start = nz_num_subrow1 + (kk-2)*5*n_moment_LDG &
            +(element_poisson(1)%pivot%irank-1 )*n_moment_LDG
        ie_start = nz_num_subrow1 + (kk-2)*5*n_moment_LDG &
            +(element_poisson(1)%east%irank-1 )*n_moment_LDG
        iw_start = nz_num_subrow1 + (kk-2)*5*n_moment_LDG &
            +(element_poisson(1)%west%irank-1 )*n_moment_LDG
        in_start = nz_num_subrow1 + (kk-2)*5*n_moment_LDG &
            +(element_poisson(1)%north%irank-1 )*n_moment_LDG
        is_start = nz_num_subrow1 + (kk-2)*5*n_moment_LDG &
            +(element_poisson(1)%south%irank-1 )*n_moment_LDG
        !********
        aa( ip_start+1:ip_start+n_moment_LDG ) = block_dp(kk,1:n_moment_LDG)
        aa( ie_start+1:ie_start+n_moment_LDG ) = block_de(kk,1:n_moment_LDG)
        aa( iw_start+1:iw_start+n_moment_LDG ) = block_dw(kk,1:n_moment_LDG)
        aa( in_start+1:in_start+n_moment_LDG ) = block_dn(kk,1:n_moment_LDG)
        aa( is_start+1:is_start+n_moment_LDG ) = block_ds(kk,1:n_moment_LDG)
        !********
        do ii = 1 , n_moment_LDG
            ja(ip_start+ii) = ( element_poisson(1)%pivot%icolumn-1 )*n_moment_LDG + ii
            ja(ie_start+ii) = ( element_poisson(1)%east%icolumn-1 )*n_moment_LDG + ii
            ja(iw_start+ii) = ( element_poisson(1)%west%icolumn-1 )*n_moment_LDG + ii
            ja(in_start+ii) = ( element_poisson(1)%north%icolumn-1 )*n_moment_LDG + ii
            ja(is_start+ii) = ( element_poisson(1)%south%icolumn-1 )*n_moment_LDG + ii
        enddo

    enddo

    nz_num_row1 = nz_num_subrow1 + ( n_moment_LDG - 1 )*5*n_moment_LDG

    !*************************************************************************************
    !*************************************************************************************
    do iele = 2,nel_x*nel_y
        do kk = 1 , n_moment_LDG
            ip_start = nz_num_row1 + (kk-1)*5*n_moment_LDG + (iele-2)*5*n_moment_LDG*n_moment_LDG &
                +( element_poisson(iele)%pivot%irank-1 )*n_moment_LDG
            ie_start = nz_num_row1 + (kk-1)*5*n_moment_LDG + (iele-2)*5*n_moment_LDG*n_moment_LDG &
                +( element_poisson(iele)%east%irank-1 )*n_moment_LDG
            iw_start = nz_num_row1 + (kk-1)*5*n_moment_LDG + (iele-2)*5*n_moment_LDG*n_moment_LDG &
                +( element_poisson(iele)%west%irank-1 )*n_moment_LDG
            in_start = nz_num_row1 + (kk-1)*5*n_moment_LDG + (iele-2)*5*n_moment_LDG*n_moment_LDG &
                +( element_poisson(iele)%north%irank-1 )*n_moment_LDG
            is_start = nz_num_row1 + (kk-1)*5*n_moment_LDG + (iele-2)*5*n_moment_LDG*n_moment_LDG &
                +( element_poisson(iele)%south%irank-1 )*n_moment_LDG
            !*******
            aa(ip_start+1:ip_start+n_moment_LDG) = block_dp(kk,1:n_moment_LDG)
            aa(ie_start+1:ie_start+n_moment_LDG) = block_de(kk,1:n_moment_LDG)
            aa(iw_start+1:iw_start+n_moment_LDG) = block_dw(kk,1:n_moment_LDG)
            aa(in_start+1:in_start+n_moment_LDG) = block_dn(kk,1:n_moment_LDG)
            aa(is_start+1:is_start+n_moment_LDG) = block_ds(kk,1:n_moment_LDG)


            !*******
            do ii = 1 , n_moment_LDG
                ja( ip_start + ii ) = ( element_poisson(iele)%pivot%icolumn-1 )*n_moment_LDG + ii
                ja( ie_start + ii ) = ( element_poisson(iele)%east%icolumn-1 )*n_moment_LDG + ii
                ja( iw_start + ii ) = ( element_poisson(iele)%west%icolumn-1 )*n_moment_LDG + ii
                ja( in_start + ii ) = ( element_poisson(iele)%north%icolumn-1 )*n_moment_LDG + ii
                ja( is_start + ii ) = ( element_poisson(iele)%south%icolumn-1 )*n_moment_LDG + ii
            enddo
        enddo
    enddo


    end subroutine poisson_CSR_format

    subroutine poisson_CSR_aa1
    implicit none
    integer :: ip_start, iw_start
    integer :: n_rank
    integer :: ii,kk

    ! store in the Compressed Sparse Row format

    na = nel_x*nel_y*n_moment_LDG
    nz_num1 = n_moment_LDG*n_moment_LDG * 2*nel_x*nel_y
    allocate( aa1(nz_num1),ja1(nz_num1) )
    allocate( ia1(na+1) )

    !!*************
    !bb(1:na) = ss_phi(1:na)
    !ss = 0.
    !!*************

    ia1(1) = 1
    do ii = 1,na
        ia1(ii+1) = ia1(ii) + 2*n_moment_LDG
    enddo

    !*************************************************************************************
    do iele = 1,nel_x*nel_y
        do kk = 1 , n_moment_LDG
            n_rank = ( sign(1,element_poisson(iele)%pivot%irank - element_poisson(iele)%west%irank) + 1)/2

            ip_start =  (kk-1)*2*n_moment_LDG + (iele-1)*2*n_moment_LDG*n_moment_LDG &
                + n_rank *n_moment_LDG

            n_rank = ( sign(1,element_poisson(iele)%west%irank - element_poisson(iele)%pivot%irank) + 1)/2

            iw_start = (kk-1)*2*n_moment_LDG + (iele-1)*2*n_moment_LDG*n_moment_LDG &
                + n_rank *n_moment_LDG

            !*******
            aa1(ip_start+1:ip_start+n_moment_LDG) = block_a1(kk,1:n_moment_LDG)

            aa1(iw_start+1:iw_start+n_moment_LDG) = block_a1w(kk,1:n_moment_LDG)
            !*******
            do ii = 1 , n_moment_LDG
                ja1( ip_start + ii ) = ( element_poisson(iele)%pivot%icolumn-1 )*n_moment_LDG + ii
                ja1( iw_start + ii ) = ( element_poisson(iele)%west%icolumn-1 )*n_moment_LDG + ii
            enddo
            ! dot product
            !do ii = 1 , n_moment_LDG
            !    ss( (iele-1)*n_moment_LDG +kk ) = ss( (iele-1)*n_moment_LDG +kk ) &
            !        + aa1( ip_start +ii )*bb( ja1(ip_start + ii) ) + aa1( iw_start + ii )*bb( ja1(iw_start+ii) )
            !
            !enddo

        enddo
    enddo

    !ss(1:na) = -ss(1:na)



    end subroutine poisson_CSR_aa1


    subroutine poisson_CSR_aa2
    implicit none
    integer :: ip_start, is_start
    integer :: n_rank
    integer :: ii,kk
    ! store in the Compressed Sparse Row format

    allocate( aa2(nz_num1),ja2(nz_num1) )
    allocate( ia2(na+1) )
    !!*************
    !bb(1:na) = ss_phi(1:na)
    !ss = 0.
    !!*************

    ia2(1) = 1
    do ii = 1,na
        ia2(ii+1) = ia2(ii) + 2*n_moment_LDG
    enddo

    !*************************************************************************************
    do iele = 1,nel_x*nel_y
        do kk = 1 , n_moment_LDG
            n_rank = ( sign(1,element_poisson(iele)%pivot%irank - element_poisson(iele)%south%irank) + 1)/2
            ip_start =  (kk-1)*2*n_moment_LDG + (iele-1)*2*n_moment_LDG*n_moment_LDG &
                + n_rank *n_moment_LDG

            n_rank = ( sign(1,element_poisson(iele)%south%irank - element_poisson(iele)%pivot%irank) + 1)/2
            is_start = (kk-1)*2*n_moment_LDG + (iele-1)*2*n_moment_LDG*n_moment_LDG &
                + n_rank *n_moment_LDG

            !*******
            aa2(ip_start+1:ip_start+n_moment_LDG) = block_a2(kk,1:n_moment_LDG)

            aa2(is_start+1:is_start+n_moment_LDG) = block_a2s(kk,1:n_moment_LDG)
            !*******
            do ii = 1 , n_moment_LDG
                ja2( ip_start + ii ) = ( element_poisson(iele)%pivot%icolumn-1 )*n_moment_LDG + ii
                ja2( is_start + ii ) = ( element_poisson(iele)%south%icolumn-1 )*n_moment_LDG + ii
            enddo

        enddo
    enddo

    !ss(1:na) = -ss(1:na)



    end subroutine poisson_CSR_aa2


    subroutine Poisson_2D_LDG
    implicit none
    integer :: irow
    integer :: ii
    integer :: i,j

    call poisson_rhs_vector

    if(nt==0)then
        ss = 0.
    else
        ss(1:nel_x*nel_y*n_moment_LDG) = ss_phi(1:nel_x*nel_y*n_moment_LDG)
    endif

    call pmgmres_ilu_cr(na,nz_num,ia,ja,aa,ss,bb,na/40+1,40,1.0e-10,iter,err) !1d-6 enough for p1

    ss_phi(1:nel_x*nel_y*n_moment_LDG) = ss(1:nel_x*nel_y*n_moment_LDG)

    call poisson_phi_x
    !**************************transfer to x-y coordinate
    irow = 0
    do iele = 1 , nel_x*nel_y
        phi_x( inx(iele),iny(iele),1:n_moment_LDG ) = ss( irow + 1 : irow + n_moment_LDG )
        irow = irow + n_moment_LDG
    enddo
    !**********************************************************************************
    call poisson_phi_y
    !**************************transfer to x-y coordinate
    irow = 0
    do iele = 1 , nel_x*nel_y
        phi_y( inx(iele),iny(iele),1:n_moment_LDG ) = ss( irow + 1 : irow + n_moment_LDG )
        irow = irow + n_moment_LDG
    enddo
    !**********************************************************************************

    do i = 1 , nel_x
        do j = 1 ,1
            phi_x(i,1-j,:) = phi_x(i,nel_y+1 -j,:)
            phi_x(i,nel_y +j,:) = phi_x(i,0+j,:)

            phi_y(i,1-j,:) = phi_y(i,nel_y+1 -j,:)
            phi_y(i,nel_y +j,:) = phi_y(i,0+j,:)
        enddo
    enddo

    do j = 1 -1    , nel_y +1
        do i = 1 , 1
            phi_x(1-i,j,:) = phi_x(nel_x+1 -i,j,:)
            phi_x(nel_x +i,j,:) = phi_x(0+i,j,:)

            phi_y(1-i,j,:) = phi_y(nel_x+1 -i,j,:)
            phi_y(nel_x +i,j,:) = phi_y(0+i,j,:)
        enddo
    enddo

    end subroutine Poisson_2D_LDG


    subroutine poisson_rhs_vector
    implicit none
    integer :: irow
    integer :: m,n

    integer :: ii,jj



    bb = 0.0
    irow = 0

    do iele = 1,nel_x*nel_y
        !***********************
        if( iele/nel_x*nel_x == iele )then
            iny(iele) = iele/nel_x
            inx(iele) = iele - nel_x*(iny(iele)-1)
        else
            iny(iele) = iele/nel_x + 1
            inx(iele) = iele -  iele/nel_x*nel_x
        endif
        !***********************

        do ii = 1 , n_moment_2d
            bb( irow + ii ) =  temp11( inx(iele) ,iny(iele), ii ) * c(ii) * dx * dy
        enddo

        irow = irow + n_moment_LDG
    enddo



    end subroutine poisson_rhs_vector


    subroutine poisson_phi_x
    implicit none
    integer :: n_rank,kk,ii,ip_start,iw_start

    !*************
    bb(1:na) = ss_phi(1:na)
    ss = 0.
    !*************
    do iele = 1,nel_x*nel_y
        do kk = 1 , n_moment_LDG
            n_rank = ( sign(1,element_poisson(iele)%pivot%irank - element_poisson(iele)%west%irank) + 1)/2

            ip_start =  (kk-1)*2*n_moment_LDG + (iele-1)*2*n_moment_LDG*n_moment_LDG &
                + n_rank *n_moment_LDG

            n_rank = ( sign(1,element_poisson(iele)%west%irank - element_poisson(iele)%pivot%irank) + 1)/2

            iw_start = (kk-1)*2*n_moment_LDG + (iele-1)*2*n_moment_LDG*n_moment_LDG &
                + n_rank *n_moment_LDG

            !*******
            ! dot product
            do ii = 1 , n_moment_LDG
                ss( (iele-1)*n_moment_LDG +kk ) = ss( (iele-1)*n_moment_LDG +kk ) &
                    + aa1( ip_start +ii )*bb( ja1(ip_start + ii) ) + aa1( iw_start + ii )*bb( ja1(iw_start+ii) )

            enddo

        enddo
    enddo

    ss(1:na) = sign_rho * ss(1:na)

    end subroutine poisson_phi_x


    subroutine poisson_phi_y
    implicit none
    integer :: ip_start, is_start
    integer :: n_rank
    integer :: ii,kk

    !*************
    bb(1:na) = ss_phi(1:na)
    ss = 0.
    !*************

    !*************************************************************************************
    do iele = 1,nel_x*nel_y
        do kk = 1 , n_moment_LDG
            n_rank = ( sign(1,element_poisson(iele)%pivot%irank - element_poisson(iele)%south%irank) + 1)/2
            ip_start =  (kk-1)*2*n_moment_LDG + (iele-1)*2*n_moment_LDG*n_moment_LDG &
                + n_rank *n_moment_LDG

            n_rank = ( sign(1,element_poisson(iele)%south%irank - element_poisson(iele)%pivot%irank) + 1)/2
            is_start = (kk-1)*2*n_moment_LDG + (iele-1)*2*n_moment_LDG*n_moment_LDG &
                + n_rank *n_moment_LDG

            ! dot product
            do ii = 1 , n_moment_LDG
                ss( (iele-1)*n_moment_LDG +kk ) = ss( (iele-1)*n_moment_LDG +kk ) &
                    + aa2( ip_start +ii )*bb( ja2(ip_start + ii) ) + aa2( is_start + ii )*bb( ja2(is_start+ii) )
            enddo
        enddo
    enddo

    ss(1:na) = sign_rho * ss(1:na)

    end subroutine poisson_phi_y


    end module poisson_periodic_LDG_mod