    module LU
    !--------------------------------------------------module coment
    !  Description :  LU matrix factorization
    !
    !---------------------------------------------------------------

    contains

    subroutine solve(A,b,xx,N)
    implicit none
    integer,intent(in) :: N
    real,intent(in) :: A(N,N),b(N)
    real,intent(out) :: xx(N)

    real :: L(N,N),U(N,N)

    real :: yy(N)

    call doolittle(A,L,U,N)

    call  downtri(L,b,yy,N)

    call uptri(U,yy,xx,N)

    end subroutine solve

    subroutine solve17(L,U,b,xx,N)
    implicit none
    integer,intent(in) :: N
    real,intent(in) :: L(N,N),U(N,N),b(N)
    real,intent(out) :: xx(N)

    real :: yy(N)

    call  downtri(L,b,yy,N)

    call uptri(U,yy,xx,N)

    end subroutine solve17

    subroutine doolittle(A,L,U,N)
    !-------------------------------------------subroutine  comment
    !  Purpose   :  Doolittle function
    !              A=LU
    !--------------------------------------------------------------
    !  Input  parameters  :
    !       1.    A  n*n matrix
    !       2.    N  order
    !  Output parameters  :
    !       1.   L
    !       2.   U
    !---------------------------------------------------------------
    implicit none
    integer,intent(in) :: n
    real,intent(in) :: A(N,N)
    real,intent(out) ::  L(N,N),U(N,N)
    integer :: k,i,j,m
    real :: s

    ! the first row of U

    U(1,:)=A(1,:)

    ! the first column of L
    L(:,1)=a(:,1)/U(1,1)

    do k=2,N

        l(k,k)=1

        do j=k,n
            s=0
            do m=1,k-1
                s=s+l(k,m)*u(m,j)
            end do
            u(k,j)=a(k,j)-s
        end do


        do i=k+1,n
            s=0
            do m=1,k-1
                s=s+l(i,m)*u(m,k)
            end do
            l(i,k)=(a(i,k)-s)/u(k,k)

        end do


    end do

    end subroutine doolittle


    subroutine uptri(A,b,xx,N)
    !--------------------------------------------subroutine  comment
    !  Purpose   :  backward substitution of upper triangular matrix
    !                 Ax=b
    !---------------------------------------------------------------
    !  Input  parameters  :
    !       1.   A(N,N)
    !       2.   b(N)
    !       3.   N
    !  Output parameters  :
    !       1.  xx
    !
    !---------------------------------------------------------------
    implicit none

    integer,intent(in) :: n

    real,intent(in) :: A(N,N),b(N)
    real,intent(out) :: xx(N)
    integer :: i,j

    xx(N)=b(N)/A(N,N)

    ! backward substitution
    do i=n-1,1,-1

        xx(i)=b(i)
        do j=i+1,N
            xx(i)=xx(i)-a(i,j)*xx(j)
        end do
        xx(i)=xx(i)/A(i,i)

    end do

    end subroutine uptri


    subroutine downtri(A,b,xx,N)
    !--------------------------------------------subroutine  comment
    !  Purpose   :  backward substitution of lower triangular matrix
    !                 Ax=b
    !----------------------------------------------------------------
    !  Input  parameters  :
    !       1.   A(N,N)
    !       2.   b(N)
    !       3.   N
    !  Output parameters  :
    !       1.  xx
    !
    !---------------------------------------------------------------
    implicit none
    integer,intent(in) :: N
    real,intent(in) :: A(N,N),b(N)
    real,intent(out) :: xx(N)
    integer :: k,i,j

    xx(1)=b(1)/a(1,1)

    do k=2,N
        xx(k)=b(k)
        do i=1,k-1
            xx(k)=xx(k)-a(k,i)*xx(i)
        end do
        xx(k)=xx(k)/a(k,k)

    end do

    end subroutine downtri

    end module LU