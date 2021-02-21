
    module poisson_element_mod

    type,public :: block_type
        sequence
        integer :: irank
        integer :: icolumn
    end type
    type,public :: element_poisson_type
        sequence
        integer :: inside_or_boundary
        type(block_type) pivot
        type(block_type) east
        type(block_type) west
        type(block_type) north
        type(block_type) south
    end type

    type(element_poisson_type),allocatable, public :: element_poisson(:)


    end module poisson_element_mod
    