    !***************************************************************************
    ! semi-Lagrangian Discontinuous Galerkin (SLDG) with operator splitting
    ! for incompressible Euler equations
    !
    ! 2D Cartesian plane and a Nodal DG formulation
    !                  code by
    !                         Jingmei Qiu, jingqiu@udel.edu
    !                         Xiaofeng Cai, xfcai@udel.edu
    !                         Mingchang Ding,
    !                         March 22th, 2018, University of Delaware
    !                         March 28th, 2018, University of Delaware
    !***************************************************************************
    Program sldg2d_split
    use polynomial_mod
    use variable
    use element_mod

    use poisson_system_solver
    use poisson_element_mod
    use poisson_periodic_LDG_mod

    implicit none

    do kkkk = 1,5
        !nel_x = 10*2**(kkkk-1)
        nel_x = 20*kkkk
        nel_y = nel_x

        norder(kkkk) = nel_x
        nel = max(nel_x,nel_y)

        nod = 1 ! nod stands for p^nod 1d polynomial
        n_moment = nod +1

        sign_rho = 1. ! if sign_rho = 1, the code is for incompressible Euler equations.

        n_gl = 2
        n_g = nod + 1
        icase = 0
        irk = 5
        tprint = 1.
        cfl = 1.
        ighost = int(cfl) + 3

        call allocate_variable
        call CPU_TIME(begin_time)


        ! grid, GL points & Initial (exact) data
        call grid_maker_GL
        call advection_dat

        call poisson_parameters
        call allocate_var_poisson
        call poisson_coefficients
        call poisson_CSR_init
        call poisson_CSR_format
        call poisson_CSR_aa1
        call poisson_CSR_aa2

        ! strang splitting
        tn = 0.0d0
        nt = 0
        do while(tn<tprint)         !Time Loop
            call splitting
        enddo ! do while
        call CPU_TIME(end_time)
        write(2016,*) 'total_time',end_time-begin_time
        call order_DG
        call deallocate_variable
        call deallocate_var_poisson
    enddo

    contains
    include "grid_maker_GL.f90"
    include "allocate_variable.f90"
    include "advection_dat.f90"
    include "setdt.f90"
    
    include "SLDG1D_QiuGuoCai.f90"
    include "id_get.f90"
    include "green2.f90"
    include "green3.f90"

    include "splitting.f90"
    include "boundary.f90"
    include "order.f90"
    include "output.f90"

    end program sldg2d_split