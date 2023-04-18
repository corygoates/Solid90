program main

    use element_mod
    use gauss_integration_mod
    use json_mod
    use json_xtnsn_mod

    implicit none

    type(element) :: test_element
    real :: coef
    character(100) :: input_file
    type(json_file) :: input_json
    type(json_value),pointer :: plate_settings
    logical :: found

    ! Initialize
    call test_element%init((/0., 0./), &
                           (/1., 0./), &
                           (/1., 1./), &
                           (/0., 1./), &
                           1, 2, 3, 4)

    coef = test_element%get_K11_ij(1, 1, 200., 300.)
    write(*,*) coef

    ! Set up
    call json_initialize()

    ! Get input file
    call getarg(1, input_file)
    input_file = trim(input_file)

    ! Load JSON
    call input_json%load_file(filename=input_file)
    call json_check()

    ! Initialize plate
    call input_json%get('plate', plate_settings, found)

end program main