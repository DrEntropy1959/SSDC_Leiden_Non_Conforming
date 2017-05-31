! This module contains subroutine and functions used for the mesh converter. 

module gmsh_converter_utilities

  ! Public subroutine and functions
  public get_file_unit
  public set_n_nodes_entities

  contains

  !============================================================================
  
  !============================================================================

  function get_file_unit(max_unit)
    ! This function returns a unit number that is not in use

    ! Nothing is defined implicitly
    implicit none

    integer, intent(in) :: max_unit
    integer :: i_unit
    integer :: j
    integer :: iostat
    logical :: opened
    integer :: get_file_unit

    continue

    ! Check which unit is free in a descending order
    j = max_unit

    do i_unit = j,1,-1
      inquire(unit=i_unit,opened=opened,iostat=iostat)
      if (iostat .ne. 0) then 
        cycle
        if (.not. opened) then 
          exit
        endif
      endif
    enddo

    ! Assign free unit to get_file_unit
    get_file_unit = i_unit

    return
  end function get_file_unit

  !============================================================================
  
  !============================================================================
  
  function set_n_nodes_entities()
    
    ! Nothing is implicitly defined
    implicit none

    integer, dimension(93) :: set_n_nodes_entities

    continue

    ! Set number of vertices for each gmsh entity
    set_n_nodes_entities = 0

    set_n_nodes_entities(1)  = 2
    set_n_nodes_entities(2)  = 3 
    set_n_nodes_entities(3)  = 4
    set_n_nodes_entities(4)  = 4
    set_n_nodes_entities(5)  = 8
    set_n_nodes_entities(6)  = 6
    set_n_nodes_entities(7)  = 5
    set_n_nodes_entities(8)  = 3
    set_n_nodes_entities(9)  = 6
    set_n_nodes_entities(10) = 9
    set_n_nodes_entities(11) = 10
    set_n_nodes_entities(12) = 27
    set_n_nodes_entities(13) = 18
    set_n_nodes_entities(14) = 14
    set_n_nodes_entities(15) = 1
    set_n_nodes_entities(16) = 8
    set_n_nodes_entities(17) = 20
    set_n_nodes_entities(18) = 15
    set_n_nodes_entities(19) = 13
    set_n_nodes_entities(20) = 9
    set_n_nodes_entities(21) = 10
    set_n_nodes_entities(22) = 12
    set_n_nodes_entities(23) = 15
    set_n_nodes_entities(24) = 15
    set_n_nodes_entities(25) = 21
    set_n_nodes_entities(26) = 4
    set_n_nodes_entities(27) = 5
    set_n_nodes_entities(28) = 6
    set_n_nodes_entities(29) = 20
    set_n_nodes_entities(30) = 35
    set_n_nodes_entities(31) = 56
    set_n_nodes_entities(92) = 64
    set_n_nodes_entities(93) = 125

    return
  end function set_n_nodes_entities

  !============================================================================
  
  !============================================================================

end module gmsh_converter_utilities

