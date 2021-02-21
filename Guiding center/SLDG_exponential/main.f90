    !*******************************************************************
    ! A well-defined SLDG code for
    ! the incompressible Euler equations with
    ! the periodic boundary conditions in two directions.
    !             Author:
    !                    Xiaofeng Cai
    !*******************************************************************
    program SLDG2D
    use module_data2d
    use globals2d
    use LU

    use poisson_element_mod
    use poisson_periodic_LDG_mod
    implicit none

    call parameters

    allocate( ref(1:100,1:100,1:6 ,1:6) )
    do kkkk = 1,5

        call setup
        call allocate_var
        call init

        !********************************************************************
        call poisson_parameters
        nel_x = nx
        nel_y = ny
        sign_rho = 1. ! if sign_rho = 1, the code is for incompressible Euler equations.
        call allocate_var_poisson
        call poisson_coefficients
        call poisson_CSR_init
        call poisson_CSR_format
        call poisson_CSR_aa1
        call poisson_CSR_aa2
        !********************************************************************
 

        time = 0.
        nt = 0
        !******************** BEGIN TIME EVOLUTION ***************************
        do while(time<time_final  )
            !call get_norm
            do io = 0,iexprk-1
                call boundary
                call sl_exp_rk
            enddo

            !update solution!
            do i = 1 ,nx
                do j = 1,ny
                    element(i,j,0)%umodal(1:n_moment) =  element(i,j,iexprk)%umodal(1:n_moment)
                enddo
            enddo

            nt = nt + 1
            if(nt/1*1==nt) print *,time,time/time_final*100,"%"
        enddo

        !call get_norm
        call order_DG

        !call order_time

        !call order
        !call output
        !pause

        call deallocate_var

        call deallocate_var_poisson
    enddo

    deallocate( ref )

    contains

    include "sl_exp_rk.f90"

    include "RK2D.f90"

    include "get_norm.f90"

    include "setdt.f90"
    include "setup.f90"
    include "order_DG.f90"

    include "parameters.f90"

    include "init.f90"
    include "boundary.f90"

    include "get_integral.f90"
    include "green.f90"
 
    include "green_gauss_line_integral.f90"

    include "get_matrix_vector.f90"
    include "get_intersections_outersegments.f90"
    include "get_intersections.f90"
    include "get_intersections_qc.f90"
    include "get_outersegments.f90"
    include "get_innersegments.f90"

    include "polynomials.f90"

    include "output.f90"

    include "get_mono_inters.f90"

    include "super2final_inner_x.f90"
    include "super2final_inner_y.f90"


    include "get_subfaces.f90"
    include "get_mono_inters_QC.f90"

    include "green_p2qc_gauss3_outer.f90"

    !include "achive_get_intersections.f90"
    include "SLDG_poisson1d.f90"


    include "order_time.f90"

    include "data_save.f90"


    include "velocity_nodes.f90"
    include "get_velocity_nodes.f90"

    include "order.f90"
    
    include "allocate_var.f90"
    
    ! the following interpolate_gauss_subroutines.f90 is added at
    ! 07/19/2019 for 
    ! interpolating the test function by gauss points.
    include "interpolate_gauss_subroutines.f90"
    include "get_rk_stage_velocity.f90"
    !*******************************************************************
    subroutine RK_to_upstream
    implicit none


    do i = 1,nx+1
        do j = 1,ny+1

            call RK2D( vertex( i,j )%coor(1:2),vertex_star( i,j )%coor(1:2),time,dt,-1. )

            vertex_star( i,j )%id(1) = ceiling( (vertex_star( i,j )%coor(1)-xleft)/dx )
            vertex_star( i,j )%id(2) = ceiling( (vertex_star( i,j )%coor(2)-ybottom)/dy )
        enddo
    enddo

    do i = 1,nx
        do j = 1,ny+1
            call RK2D( nodex( i,j )%coor(1:2),nodex_star( i,j )%coor(1:2),time,dt,-1. )
        enddo
    enddo

    do i = 1,nx+1
        do j = 1,ny
            call RK2D( nodey( i,j )%coor(1:2),nodey_star( i,j )%coor(1:2),time,dt,-1. )
        enddo
    enddo

    do i = 1,nx
        do j = 1,ny
            call RK2D( nodec( i,j )%coor(1:2),nodec_star( i,j )%coor(1:2),time,dt,-1. )
        enddo
    enddo

    ! get the vertexes of face_lr
    do j = 1 , ny+1
        do i = 1 , nx
            face_lr(i,j)%point_origin = vertex_star(i,j)
            face_lr(i,j)%point_end = vertex_star(i+1,j)
            ! quadratic-curved
            face_lr(i,j)%point_midt = nodex_star(i,j)
        enddo
    enddo

    ! get the vertexes of face_bt
    do i = 1 , nx+1
        do j = 1 , ny
            face_bt(i,j)%point_origin = vertex_star(i,j)
            face_bt(i,j)%point_end = vertex_star(i,j+1)
            ! quadratic-curved
            face_bt(i,j)%point_midt = nodey_star(i,j)
        enddo
    enddo
    
    


    end subroutine RK_to_upstream


    end program SLDG2D