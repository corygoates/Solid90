module element_mod

    use gauss_integration_mod

    implicit none

    type node

        real,dimension(2) :: loc
        integer :: index

    end type node

    type master_element
        ! A linear element

        real,dimension(2,4) :: node_locs

        contains

            procedure :: init => master_element_init
            procedure :: psi => master_element_psi
            procedure :: d_psi_d_xi => master_element_d_psi_d_xi
            procedure :: d_psi_d_eta => master_element_d_psi_d_eta

    end type master_element


    type element

        real,dimension(2,4) :: node_locs
        integer,dimension(4) :: global_node_indices
        real,dimension(2,2) :: J, J_star
        real :: det_J, A
        type(master_element) :: master
        real,dimension(2) :: centr

        contains

            procedure :: init => element_init
            procedure :: calc_J => element_calc_J
            procedure :: global_to_local => element_global_to_local
            procedure :: local_to_global => element_local_to_global
            procedure :: psi => element_psi
            procedure :: d_psi_d_x => element_d_psi_d_x
            procedure :: d_psi_d_y => element_d_psi_d_y
            procedure :: get_K11_ij => element_get_K11_ij
            procedure :: get_K12_ij => element_get_K12_ij
            procedure :: get_K13_ij => element_get_K13_ij
            procedure :: get_K21_ij => element_get_K21_ij
            procedure :: get_K22_ij => element_get_K22_ij
            procedure :: get_K23_ij => element_get_K23_ij
            procedure :: get_K31_ij => element_get_K31_ij
            procedure :: get_K32_ij => element_get_K32_ij
            procedure :: get_K33_ij => element_get_K33_ij
            procedure :: get_elemental_stiffness => element_get_elemental_stiffness

    end type element
    
contains

    subroutine master_element_init(this)
        ! Initializes the master element

        implicit none
        
        class(master_element), intent(inout) :: this
    
        ! Node locations
        this%node_locs(:,1) = (/-1., -1./)
        this%node_locs(:,2) = (/1., -1./)
        this%node_locs(:,3) = (/1., 1./)
        this%node_locs(:,4) = (/-1., 1./)
        
    end subroutine master_element_init


    function master_element_psi(this, xi, eta, j) result(P)
        ! Shape function

        implicit none
        
        class(master_element),intent(in) :: this
        real,intent(in) :: xi, eta
        integer :: j

        real :: P
    
        select case (j)

        case (1)
            P = 0.25*(1.-xi)*(1.-eta)

        case (2)
            P = 0.25*(1.+xi)*(1.-eta)

        case (3)
            P = 0.25*(1.+xi)*(1.+eta)

        case (4)
            P = 0.25*(1.-xi)*(1.+eta)

        end select
        
    end function master_element_psi


    function master_element_d_psi_d_xi(this, xi, eta, j) result(dP)
        ! xi derivative of the shape function

        implicit none
        
        class(master_element),intent(in) :: this
        real,intent(in) :: xi, eta
        integer :: j

        real :: dP
    
        select case (j)

        case (1)
            dP = -0.25*(1.-eta)

        case (2)
            dP = 0.25*(1.-eta)

        case (3)
            dP = 0.25*(1.+eta)

        case (4)
            dP = -0.25*(1.+eta)

        end select
        
    end function master_element_d_psi_d_xi


    function master_element_d_psi_d_eta(this, xi, eta, j) result(dP)
        ! eta derivative of the shape function

        implicit none
        
        class(master_element),intent(in) :: this
        real,intent(in) :: xi, eta
        integer :: j

        real :: dP
    
        select case (j)

        case (1)
            dP = -0.25*(1.-xi)

        case (2)
            dP = -0.25*(1.+xi)

        case (3)
            dP = 0.25*(1.+xi)

        case (4)
            dP = 0.25*(1.-xi)

        end select
        
    end function master_element_d_psi_d_eta


    subroutine element_init(this, loc1, loc2, loc3, loc4, i1, i2, i3, i4)
        ! Initializes the element, assumed to be rectangular

        implicit none
        
        class(element),intent(inout) :: this
        real,dimension(2),intent(in) :: loc1, loc2, loc3, loc4
        integer, intent(in) :: i1, i2, i3, i4

        ! Store
        this%node_locs(:,1) = loc1
        this%node_locs(:,2) = loc2
        this%node_locs(:,3) = loc3
        this%node_locs(:,4) = loc4
        this%global_node_indices = (/ i1, i2, i3, i4/)

        ! Calculate center and area
        this%centr(1) = 0.5*(loc2(1) + loc1(1))
        this%centr(2) = 0.5*(loc3(2) + loc2(2))
        this%A = (loc2(1) - loc1(1))*(loc3(2) - loc2(2))

        ! Initialize some other things
        call this%master%init()
        call this%calc_J()
        
    end subroutine element_init


    subroutine element_calc_J(this)
        ! Calculates the Jacobian matrices

        implicit none
        
        class(element), intent(inout) :: this
    
        ! Caclulate Jacobian
        this%J(1,1) = 0.5*(this%node_locs(1,2) - this%node_locs(1,1))
        this%J(2,2) = 0.5*(this%node_locs(2,3) - this%node_locs(2,2))
        this%J(1,2) = 0.5*(this%node_locs(2,2) - this%node_locs(2,1))
        this%J(2,1) = 0.5*(this%node_locs(1,3) - this%node_locs(1,2))

        ! Calculate determinant
        this%det_J = this%J(1,1)*this%J(2,2) - this%J(1,2)*this%J(2,1)

        ! Calculate inverse
        this%J_star(1,1) = this%J(2,2)/this%det_J
        this%J_star(2,2) = this%J(1,1)/this%det_J
        this%J_star(1,2) = this%J(1,2)/this%det_J
        this%J_star(2,1) = this%J(2,1)/this%det_J
        
    end subroutine element_calc_J


    function element_global_to_local(this, x_y) result(xi_eta)
        ! Converts the global location to local coordinates

        implicit none
        
        class(element),intent(in) :: this
        real,dimension(2),intent(in) :: x_y

        real,dimension(2) :: xi_eta

        xi_eta = matmul(this%J_star, x_y - this%centr)
        
    end function element_global_to_local


    function element_local_to_global(this, xi_eta) result(x_y)
        ! Converts the local location to global coordinates

        implicit none
        
        class(element),intent(in) :: this
        real,dimension(2),intent(in) :: xi_eta

        real,dimension(2) :: x_y

        x_y = matmul(this%J, xi_eta) + this%centr
        
    end function element_local_to_global


    function element_psi(this, x, y, j) result(P)
        ! Calculates the shape function for this element

        implicit none
        
        class(element),intent(in) :: this
        real,intent(in) :: x, y
        integer :: j

        real :: P

        real,dimension(2) :: xi_eta

        ! Transform to local
        xi_eta = this%global_to_local((/x, y/))

        ! Get psi
        P = this%master%psi(xi_eta(1), xi_eta(2), j)

    end function element_psi


    function element_d_psi_d_x(this, xi, eta, j) result(dP)
        ! Calculates the derivative of the shape function for this element wrt x

        implicit none
        
        class(element),intent(in) :: this
        real,intent(in) :: xi, eta
        integer :: j

        real :: dP

        real :: dp_dxi, dp_deta

        ! Get derivatives
        dp_dxi = this%master%d_psi_d_xi(xi, eta, j)
        dp_deta = this%master%d_psi_d_eta(xi, eta, j)

        dP = this%J_star(1,1)*dp_dxi + this%J_star(1,2)*dp_deta

    end function element_d_psi_d_x


    function element_d_psi_d_y(this, xi, eta, j) result(dP)
        ! Calculates the derivative of the shape function for this element wrt y

        implicit none
        
        class(element),intent(in) :: this
        real,intent(in) :: xi, eta
        integer :: j

        real :: dP

        real :: dp_dxi, dp_deta

        ! Get derivatives
        dp_dxi = this%master%d_psi_d_xi(xi, eta, j)
        dp_deta = this%master%d_psi_d_eta(xi, eta, j)

        dP = this%J_star(2,1)*dp_dxi + this%J_star(2,2)*dp_deta

    end function element_d_psi_d_y


    function element_get_K11_ij(this, i, j, A44, A55) result(K11)
        ! Returns the K11 stiffness coefficient for the given node j and the shape function i

        implicit none
        
        class(element),intent(in) :: this
        integer,intent(in) :: i, j
        real,intent(in) :: A44, A55

        real :: K11

        ! Integrate
        K11 = gauss_quad_2D(integrand, 1) ! Reduced to prevent shear locking

        contains

            function integrand(xi, eta) result(int)
                ! Integrand of the K11 integral

                implicit none
                
                real,intent(in) :: xi, eta

                real :: int

                int = (A55*this%d_psi_d_x(xi, eta, j)*this%d_psi_d_x(xi, eta, i) &
                    + A44*this%d_psi_d_y(xi, eta, j)*this%d_psi_d_y(xi, eta, i))*this%det_J

            end function integrand
        
    end function element_get_K11_ij


    function element_get_K12_ij(this, i, j, A55) result(K12)
        ! Returns the K12 stiffness coefficient for the given node j and the shape function i

        implicit none
        
        class(element),intent(in) :: this
        integer,intent(in) :: i, j
        real,intent(in) :: A55

        real :: K12

        ! Integrate
        K12 = gauss_quad_2D(integrand, 1) ! Reduced to prevent shear locking

        contains

            function integrand(xi, eta) result(int)
                ! Integrand of the K12 integral

                implicit none
                
                real,intent(in) :: xi, eta

                real :: int

                int = A55*this%master%psi(xi, eta, j)*this%d_psi_d_x(xi, eta, i)*this%det_J

            end function integrand
    
    end function element_get_K12_ij


    function element_get_K13_ij(this, i, j, A44) result(K13)
        ! Returns the K13 stiffness coefficient for the given node j and the shape function i

        implicit none
        
        class(element),intent(in) :: this
        integer,intent(in) :: i, j
        real,intent(in) :: A44

        real :: K13

        ! Integrate
        K13 = gauss_quad_2D(integrand, 1) ! Reduced to prevent shear locking

        contains

            function integrand(xi, eta) result(int)
                ! Integrand of the K13 integral

                implicit none
                
                real,intent(in) :: xi, eta

                real :: int

                int = A44*this%master%psi(xi, eta, j)*this%d_psi_d_y(xi, eta, i)*this%det_J

            end function integrand
    
        
    end function element_get_K13_ij


    function element_get_K21_ij(this, i, j, A55) result(K21)
        ! Returns the K21 stiffness coefficient for the given node j and the shape function i

        implicit none
        
        class(element),intent(in) :: this
        integer,intent(in) :: i, j
        real,intent(in) :: A55

        real :: K21

        ! Cheat
        K21 = this%get_K12_ij(j, i, A55)
    
    end function element_get_K21_ij


    function element_get_K22_ij(this, i, j, A55, D11, D66) result(K22)
        ! Returns the K22 stiffness coefficient for the given node j and the shape function i

        implicit none
        
        class(element),intent(in) :: this
        integer,intent(in) :: i, j
        real,intent(in) :: A55, D11, D66

        real :: K22

        ! Integrate
        K22 = gauss_quad_2D(integrand, 2) + gauss_quad_2D(integrand_red, 1)

        contains

            function integrand(xi, eta) result(int)
                ! Integrand of the K22 integral

                implicit none
                
                real,intent(in) :: xi, eta

                real :: int

                int = (D11*this%d_psi_d_x(xi, eta, i)*this%d_psi_d_x(xi, eta, j) &
                    + D66*this%d_psi_d_y(xi, eta, i)*this%d_psi_d_y(xi, eta, j) &
                    )*this%det_J

            end function integrand

            function integrand_red(xi, eta) result(int)
                ! Reduced integration integrand of the K22 integral

                implicit none
                
                real,intent(in) :: xi, eta

                real :: int

                int = A55*this%master%psi(xi, eta, j)*this%master%psi(xi, eta, i)*this%det_J

            end function integrand_red
    
    end function element_get_K22_ij


    function element_get_K23_ij(this, i, j, D12, D66) result(K23)
        ! Returns the K23 stiffness coefficient for the given node j and the shape function i

        implicit none
        
        class(element),intent(in) :: this
        integer,intent(in) :: i, j
        real,intent(in) :: D12, D66

        real :: K23

        ! Integrate
        K23 = gauss_quad_2D(integrand, 2)

        contains

            function integrand(xi, eta) result(int)
                ! Integrand of the K23 integral

                implicit none
                
                real,intent(in) :: xi, eta

                real :: int

                int = (D12*this%d_psi_d_x(xi, eta, i)*this%d_psi_d_y(xi, eta, j) &
                    + D66*this%d_psi_d_y(xi, eta, i)*this%d_psi_d_x(xi, eta, j) &
                    )*this%det_J

            end function integrand
    
    end function element_get_K23_ij


    function element_get_K31_ij(this, i, j, A44) result(K31)
        ! Returns the K31 stiffness coefficient for the given node j and the shape function i

        implicit none
        
        class(element),intent(in) :: this
        integer,intent(in) :: i, j
        real,intent(in) :: A44

        real :: K31

        ! Cheat
        K31 = this%get_K13_ij(j, i, A44)
    
    end function element_get_K31_ij


    function element_get_K32_ij(this, i, j, D12, D66) result(K32)
        ! Returns the K32 stiffness coefficient for the given node j and the shape function i

        implicit none
        
        class(element),intent(in) :: this
        integer,intent(in) :: i, j
        real,intent(in) :: D12, D66

        real :: K32

        ! Cheat
        K32 = this%get_K23_ij(j, i, D12, D66)
    
    end function element_get_K32_ij


    function element_get_K33_ij(this, i, j, A44, D22, D66) result(K33)
        ! Returns the K33 stiffness coefficient for the given node j and the shape function i

        implicit none
        
        class(element),intent(in) :: this
        integer,intent(in) :: i, j
        real,intent(in) :: A44, D22, D66

        real :: K33

        ! Integrate
        K33 = gauss_quad_2D(integrand, 2) + gauss_quad_2D(integrand_red, 1)

        contains

            function integrand(xi, eta) result(int)
                ! Integrand of the K23 integral

                implicit none
                
                real,intent(in) :: xi, eta

                real :: int

                int = (D22*this%d_psi_d_y(xi, eta, i)*this%d_psi_d_y(xi, eta, j) &
                    + D66*this%d_psi_d_x(xi, eta, i)*this%d_psi_d_x(xi, eta, j) &
                    )*this%det_J

            end function integrand

            function integrand_red(xi, eta) result(int)
                ! Reduced integration integrand of the K23 integral

                implicit none
                
                real,intent(in) :: xi, eta

                real :: int

                int = A44*this%master%psi(xi, eta,i)*this%master%psi(xi, eta, j)*this%det_J

            end function integrand_red
    
    end function element_get_K33_ij


    function element_get_elemental_stiffness(this, i, j, A, D) result(K)
        ! Returns the i,j elemental stiffness matrix

        implicit none
        
        class(element),intent(in) :: this
        integer,intent(in) :: i, j
        real,dimension(:,:),allocatable :: A, D

        real,dimension(3,3) :: K

        ! Get elements
        K(1,1) = this%get_K11_ij(i, j, A(4,4), A(5,5))
        K(1,2) = this%get_K12_ij(i, j, A(5,5))
        K(1,3) = this%get_K13_ij(i, j, A(4,4))
        K(2,1) = this%get_K21_ij(i, j, A(5,5))
        K(2,2) = this%get_K22_ij(i, j, A(5,5), D(1,1), D(6,6))
        K(2,3) = this%get_K23_ij(i, j, D(1,2), D(6,6))
        K(3,1) = this%get_K31_ij(i, j, A(4,4))
        K(3,2) = this%get_K32_ij(i, j, D(1,2), D(6,6))
        K(3,3) = this%get_K33_ij(i, j, A(4,4), D(2,2), D(6,6))
        
    end function element_get_elemental_stiffness

end module element_mod