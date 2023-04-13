module element_mod

    use gauss_integration_mod

    implicit none

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
        this%J(1,2) = 0.
        this%J(2,1) = 0.

        ! Calculate determinant
        this%det_J = this%J(1,1)*this%J(2,2)

        ! Calculate inverse
        this%J_star(1,1) = this%J(2,2)/this%det_J
        this%J_star(2,2) = this%J(1,1)/this%det_J
        this%J_star(1,2) = 0.
        this%J_star(2,1) = 0.
        
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
        real :: dp_dxi, dp_deta

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
        K11 = gauss_quad_2D(integrand, 3)

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

    
end module element_mod