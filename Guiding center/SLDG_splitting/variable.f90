    module variable
    implicit none

    !
    real,parameter :: pi = 4.*atan(1.)
    !

    integer :: kkkk
    integer :: nod,n_moment
    !integer,parameter :: nod = 2
    integer :: i,j,k,kk
    integer :: ig,jg
    integer :: ie,je
    integer :: icase
    integer :: irk
    ! nod is the degree of the polynomial
    integer :: n_gl,nel_x,nel_y,ighost
    integer :: nel
    real :: dx,dy
    real :: xleft,xright,yleft,yright

    real :: tprint,tn,dt,cfl
    integer :: nt

    real :: x_GL(5),w_GL(5)
    real :: x_g(5),w_g(5)
    real :: xg(6),wg(6)
    integer :: n_g
    !real :: gauss(4,1)
    real :: ai(5)
    real, allocatable :: glx(:,:,:,:),gly(:,:,:,:),fel(:,:,:,:)
    real, allocatable :: gx(:,:,:,:),gy(:,:,:,:)
    real, allocatable :: x(:),y(:)
    real, allocatable :: xgrid(:),ygrid(:)

    ! test order
    real :: er11,er22,er33
    integer :: norder(0:30)
    real :: begin_time,end_time

    !***********
    real :: ftemp

    end module variable