    module element_mod

    type, public :: elements
        sequence
        real :: psi(5,5)
        real :: split_modal(5)
    end type

    type(elements),allocatable, public :: elem( :,: )

    type,public :: sub_elements
        sequence
        real :: vpoint(2)
        integer :: id
    end type

    type(sub_elements),allocatable,public :: Dn_star_sub(:)

    end module element_mod