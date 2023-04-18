module plate_mod

    use json_mod
    use json_xtnsn_mod
    use element_mod

    implicit none

    type :: plate

        real :: a, b, h, G, E, v
        real,dimension(:,:),allocatable :: A_mat, D_mat, G_stiff
        integer :: Nx, Ny, N_elements, N_nodes, N_eq, N_dof, N_hbw
        type(element),dimension(:,:),allocatable :: elements
        type(node),dimension(:,:),allocatable :: nodes
    
    contains

        procedure :: init => plate_init
        procedure :: calc_D => plate_calc_D
        procedure :: discretize => plate_discretize
        procedure :: solve => plate_solve
        procedure :: assemble_global_stiffness => plate_assemble_global_stiffness
    
    end type plate

    
contains


    subroutine plate_init(this, settings)
        ! Initializes the plate

        implicit none
        
        class(plate),intent(inout) :: this
        type(json_value),pointer :: settings

        ! Store
        call json_xtnsn_get(settings, 'geometry.a', this%a)
        call json_xtnsn_get(settings, 'geometry.b', this%b)
        call json_xtnsn_get(settings, 'geometry.Nx', this%Nx)
        call json_xtnsn_get(settings, 'geometry.Ny', this%Ny)
        call json_xtnsn_get(settings, 'geometry.h', this%h)
        call json_xtnsn_get(settings, 'material.E', this%E)
        call json_xtnsn_get(settings, 'material.nu', this%v)

        ! Calculate shear modulus
        this%G = this%E/(2.*(1.+this%v))

        ! Initialize a few things
        call this%calc_D()
        call this%discretize()
        
    end subroutine plate_init


    subroutine plate_calc_D(this)
        ! Calculates the [D_mat] matrix

        implicit none
        
        class(plate),intent(inout) :: this

        ! Initialize storage
        allocate(this%A_mat(4:5,4:5), source=0.)
        allocate(this%D_mat(1:6,1:6), source=0.)

        ! Calculate A_mat matrix
        this%A_mat(4,4) = this%G*this%h*5./6.
        this%A_mat(5,5) = this%G*this%h*5./6.

        ! Calculate D_mat matrix
        this%D_mat(1,1) = this%E*this%h**3/(12.*(1-this%v**2))
        this%D_mat(1,2) = this%v*this%E*this%h**3/(12.*(1-this%v**2))
        this%D_mat(2,1) = this%v*this%E*this%h**3/(12.*(1-this%v**2))
        this%D_mat(2,2) = this%E*this%h**3/(12.*(1-this%v**2))
        this%D_mat(6,6) = this%G*this%h**3/12.
        
    end subroutine plate_calc_D


    subroutine plate_discretize(this)
        ! Generates the nodes and elements for the plate

        implicit none
        
        class(plate), intent(inout) :: this

        integer :: i, j, ind
        real :: dx, dy

        ! Determine numbers of elements and nodes
        this%N_elements = this%Nx*this%Ny
        this%N_nodes = (this%Nx+1)*(this%Ny+1)
        
        ! Set up nodes
        allocate(this%nodes(this%Nx+1,this%Ny+1))
        dx = this%a/this%Nx
        dy = this%b/this%Ny
        do i=1,this%Nx+1
            do j=1,this%Ny+1
                ind = j + (i-1)*(this%Ny+1)
                this%nodes(i,j)%index = ind
                this%nodes(i,j)%loc(1) = dx*(i-1)
                this%nodes(i,j)%loc(2) = dy*(j-1)
            end do
        end do

        ! Set up elements
        allocate(this%elements(this%Nx,this%Ny))
        do i=1,this%Nx
            do j=1,this%Ny
                call this%elements(i,j)%init(this%nodes(i,j)%loc, &
                                             this%nodes(i,j+1)%loc, &
                                             this%nodes(i+1,j+1)%loc, &
                                             this%nodes(i,j+1)%loc, &
                                             this%nodes(i,j)%index, &
                                             this%nodes(i,j+1)%index, &
                                             this%nodes(i+1,j+1)%index, &
                                             this%nodes(i,j+1)%index)
            end do
        end do

    end subroutine plate_discretize


    subroutine plate_solve(this)
        ! Solves the FEM over the plate
        
        implicit none
        
        class(plate), intent(inout) :: this

        ! Set up global stiffness matrix
        call this%assemble_global_stiffness()
        
    end subroutine plate_solve


    subroutine plate_assemble_global_stiffness(this)
        ! Assembles the global stiffness matrix

        implicit none
        
        class(plate),intent(inout) :: this

        integer :: i, j, k, l
        real,dimension(:,:,:,:,:,:),allocatable :: E

        ! Get elemental stiffness matrices
        allocate(E(3, 3, 4, 4, this%Nx, this%Ny))
        do k=1,this%Nx
            do l=1,this%Ny
                do i=1,4
                    do j=1,4
                        E(:,:,i,j,k,l) = this%elements(k,l)%get_elemental_stiffness(i, j, this%A_mat, this%D_mat)
                    end do
                end do
            end do
        end do

        ! Group elemental stiffnesses by node

        ! Determine number of equations
        this%N_dof = 3
        this%N_eq = this%N_nodes*this%N_dof
        
        ! Determine half bandwidth
        this%N_hbw = 0

        ! Loop through elements
        do i=1,this%Nx
            do j=1,this%Ny

                ! Loop through nodes on element
                do k=1,4
                    do l=1,4
                        this%N_hbw = max(this%N_hbw, &
                            abs(this%elements(i,j)%global_node_indices(k) - this%elements(i,j)%global_node_indices(l)))
                    end do
                end do
            end do
        end do

        ! Initialize global stiffness matrix
        allocate(this%G_stiff(this%N_eq, this%N_hbw), source=0.)
    
    end subroutine plate_assemble_global_stiffness

    
end module plate_mod