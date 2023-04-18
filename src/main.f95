program main

    use json_mod
    use json_xtnsn_mod
    use plate_mod

    implicit none

    type(plate) :: my_plate
    character(100) :: input_file
    type(json_file) :: input_json
    type(json_value),pointer :: plate_settings
    logical :: found

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
    call my_plate%init(plate_settings)

    ! Solve
    call my_plate%solve()

end program main