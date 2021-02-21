    module globals2d

    integer :: i,j
    integer :: nx,ny
    
    integer :: nel_x,nel_y
    
    real :: pi,eps
    real :: xleft,xright,ybottom,ytop
    real :: dx,dy
    real,allocatable :: xgrid(:),ygrid(:)
    real,allocatable :: x(:),y(:)
    integer :: nghost
    integer :: iexample,irk
    integer :: nk,n_moment,iqc
    integer :: nt
    integer :: irelative
    real :: r_norm1
    
    real :: time_final,time,dt,cfl
    
    real,allocatable :: umod_t(:,:,:)
    
    !********************************************
    integer :: iam,ibm,icm,idm
    integer :: isx(6),isy(6),ix,iy
    !********************************************
    
    !********************************************
    ! parameters
    real :: gau2(1:2,1:2)
    real :: gau3(1:3,1:2)
    real :: ai(1:6)
    real :: xg(1:6),wg(1:6)
    !***********
    real :: gau(1:3,1:2)
    !********************************************
    
    !********************************************
    ! order
    real :: er11,er22,er33
    integer :: kkkk
    !********************************************
    
    integer :: idebug
    
    real,allocatable :: com_mass(:,:)
    
    !*********************************************
    ! RK
    integer :: io,iexprk,iexp_case
    real :: bt_con(5,3)
    integer :: iadditive(0:5)
    
    !*********************************************
    ! variables for 1d poisson
    real :: pint
    integer :: kdg
    real :: vm(0:3),vp(0:3)
    real :: wq(4),vgau(4,0:3), vxgau(4,0:3)
    real :: ain( 4,4 ),BL(4,4)
    !****************************
    ! for calling poisson
    !real,allocatable :: temp11(:,:,:)
    real,allocatable :: ele_dg(:,:)
    real,allocatable :: ee(:,:),ee_c(:,:)
    
    !******************************
    ! for Vlasov-Poisson
    integer :: i_case
    integer :: mmp !??????????
    
    real :: er111,er222,er333
    real,allocatable :: ref(:,:,:,:)
    
    !********************************
    ! for incompressible 
    real,allocatable :: phi_x_io(:,:,:,:),phi_y_io(:,:,:,:)
    integer :: n_moment_2d,n_moment_LDG
    real,allocatable :: phi_x_v(:,:,:),phi_x_c(:,:,:),phi_x_dex(:,:,:),phi_x_dey(:,:,:)
    real,allocatable :: phi_y_v(:,:,:),phi_y_c(:,:,:),phi_y_dex(:,:,:),phi_y_dey(:,:,:)
    
    end module globals2d