module gauss_integration_mod

    implicit none

    real,dimension(1),parameter :: locs1 = (/ 0. /)
    real,dimension(2),parameter :: locs2 = (/ -0.57735027, 0.57735027 /)
    real,dimension(3),parameter :: locs3 = (/ -0.77459667, 0., 0.77459667 /)
    real,dimension(4),parameter :: locs4 = (/ -0.86113631, -0.33998104, 0.33998104, 0.86113631 /)
    real,dimension(1),parameter :: weights1 = (/ 2. /)
    real,dimension(2),parameter :: weights2 = (/ 1., 1. /)
    real,dimension(3),parameter :: weights3 = (/ 5./9., 8./9., 5./9. /)
    real,dimension(4),parameter :: weights4 = (/ 0.34785485, 0.65214515, 0.65214515, 0.34785485 /)
    
contains


    function gauss_quad_1D(f, N) result(I)
        ! Uses N-point Gaussian quadrature to integrate f from -1 to 1

        implicit none

        interface
            real function f(x)
                real,intent(in) :: x
            end function f
        end interface
        integer :: N

        real :: I

        real,dimension(:),allocatable :: locs, weights
        integer :: j

        ! Get locations and weights
        select case (N)

        case (1)
            locs = locs1
            weights = weights1

        case (2)
            locs = locs2
            weights = weights2

        case (3)
            locs = locs3
            weights = weights3

        case (4)
            locs = locs4
            weights = weights4

        end select

        ! Integrate
        I = 0.
        do j=1,N
            I = I + weights(j)*f(locs(j))
        end do
        
    end function gauss_quad_1D


    function gauss_quad_2D(f, N) result(I)
        ! Uses N-point Gaussian quadrature to integrate f over [-1,1]x[-1,1]

        implicit none

        interface
            real function f(x, y)
                real,intent(in) :: x, y
            end function f
        end interface
        integer :: N

        real :: I

        real,dimension(:),allocatable :: locs, weights
        integer :: j, k

        ! Get locations and weights
        select case (N)

        case (1)
            locs = locs1
            weights = weights1

        case (2)
            locs = locs2
            weights = weights2

        case (3)
            locs = locs3
            weights = weights3

        case (4)
            locs = locs4
            weights = weights4

        end select

        ! Integrate
        I = 0.
        do j=1,N
            do k=1,N
                I = I + weights(j)*weights(k)*f(locs(j), locs(k))
            end do
        end do
        
    end function gauss_quad_2D

    
end module gauss_integration_mod