program main

    use element_mod
    use gauss_integration_mod

    implicit none

    type(element) :: test_element
    real :: coef

    ! Initialize
    call test_element%init((/0., 0./), &
                           (/1., 0./), &
                           (/1., 1./), &
                           (/0., 1./), &
                           1, 2, 3, 4)

    coef = test_element%get_K11_ij(1, 1, 200., 300.)
    write(*,*) coef

end program main