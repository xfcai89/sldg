    module module_data2d
    !********************************************
    ! data structure of vertexes
    type, public :: type_point
        sequence
        real coor(1:2)
    end type

    type(type_point),allocatable,target, public :: vertex( :,: )
    type(type_point),allocatable,target, public :: nodex( :,: )
    type(type_point),allocatable,target, public :: nodey( :,: )
    type(type_point),allocatable,target, public :: nodec( :,: )
    !********************************************
    type, public :: type_point_upstream
        sequence
        real coor(1:2)
        integer id(1:2)
        ! qc
        real :: xxii
    end type
    type(type_point_upstream),allocatable,target, public :: vertex_star( :,: )
    type(type_point_upstream),allocatable,target, public :: nodex_star( :,: )
    type(type_point_upstream),allocatable,target, public :: nodey_star( :,: )
    type(type_point_upstream),allocatable,target, public :: nodec_star( :,: )
    !********************************************

    !********************************************
    type, public :: type_inter_point
        sequence
        real coor(1:2)
        integer ixy_type
        integer idx2
        integer idy2
        ! idx2 is for checking monotone
        ! note that on grid, index is i+1/2, so idx2=2i+1
        integer igrid
        integer id_xy
        ! for qc
        real xxii
    end type

    type(type_inter_point) :: temp
    type(type_inter_point) :: point_inner_x(10,4),point_inner_y(10,4)
    type(type_inter_point) :: point_inner(10)
    !*************************************************
    ! segment
    type, public :: type_segment
        sequence
        type(type_inter_point) :: porigin,pend
        integer :: id(1:2)
        real :: c_ab(1:6)
    end type
    !*************************************************
    ! data structure of faces
    type, public :: type_face
        sequence
        type(type_point_upstream) :: point_origin,point_end

        type(type_inter_point) :: point_inter( 0:10 )
        integer :: nsub_outer
        type(type_segment) :: segment_outer( 10 )
        !***********************************************
        ! for QC
        type(type_point_upstream) :: point_midt,point_mid1,point_mid2
        integer :: nsubface 
        integer :: ist
        real :: xa,xb,xc,ya,yb,yc
        real :: et2(1:2)
    end type

    type(type_face),allocatable,target, public :: face_lr(:,:)
    type(type_face),allocatable,target, public :: face_bt(:,:)

    !********************************************
    ! data structure of Eulerian elements
    ! attribute:
    !            modals
    type, public :: type_element
        sequence
        type(type_point),pointer :: vertex1,vertex2,vertex3,vertex4
        type(type_point),pointer :: vertex5,vertex6,vertex7,vertex8
        type(type_point),pointer :: vertex9
        real,allocatable :: umodal(:)
    end type

    type(type_element),allocatable,target, public :: element(:,:,:)
    !*******************************************************************
    ! data structure of upstream elements
    type, public :: type_element_upstream
        sequence
        type(type_point_upstream),pointer :: vertex1,vertex2,vertex3,vertex4
        type(type_point_upstream),pointer :: vertex5,vertex6,vertex7,vertex8
        type(type_point_upstream),pointer :: vertex9
        
        type(type_point) :: gauss(1:3,1:3) 
        ! 3 means the gauss points can allow for P^2
        
        type(type_face),pointer :: edge1_lr,edge2_bt,edge3_lr,edge4_bt
        integer :: nsub_inner
        type(type_segment) :: segment_inner( 30 )
    end type

    type(type_element_upstream),allocatable,target, public :: element_star(:,:,:)
    ! arguments 
    ! the first is for
    ! the second is for 
    ! the third is for

    end module module_data2d