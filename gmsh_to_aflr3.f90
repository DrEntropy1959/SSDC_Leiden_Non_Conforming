! This program converts the grid constructed in gmsh and saved as .msh to the
! AFLR3 format.
!
! IMPORTANT NOTE: Only grids formed by quadrilateral elements defined by four  
! vertices and hexahedrons elements defined by eight vertices can be read and
! converted correctly to the AFLR3 format.

program gmsh_to_aflr3

  ! Load modules
  !use precision_vars
  use gmsh_converter_utilities

  ! Nothing is implicitly defined
  implicit none

  integer, parameter :: wp = 8
  character(120)     :: file_in
  character(120)     :: header
  integer            :: length
  integer            :: i_err, max_unit, i_unit
  integer            :: file_type_gmsh
  integer            :: data_size_gmsh
  real(wp)           :: version_gmsh

  integer                :: trash_int
  integer                :: n_nodes_gmsh
  integer                :: i, j
  integer                :: n_entities_gmsh
  integer                :: entity_id_gmsh, entity_type_gmsh, n_entity_tags_gmsh
  integer                :: cnt_hexas_gmsh, cnt_quads_gmsh
  integer, dimension(93) :: n_nodes_entities
  
  integer,  allocatable, dimension(:)    :: nodes_id_gmsh
  integer,  allocatable, dimension(:)    :: entities_info_gmsh
  integer,  allocatable, dimension(:,:)  :: quads_info_gmsh, hexas_info_gmsh
  real(wp), allocatable, dimension(:,:)  :: nodes_coords_gmsh
  
  character(120) :: output_file_name
  integer        :: n_nodes, n_surf_trias, n_surf_quads, n_vol_tets, &
                  & n_vol_pents_5, n_vol_pents_6, n_vol_hexs

  integer :: tmp

  continue

  ! Set the number of nodes for each "gmsh" geometrical entity
  n_nodes_entities = set_n_nodes_entities()

  ! Parse input file name
  if (command_argument_count().gt.0) then
    call get_command_argument(1,file_in,length,i_err)
  end if

  ! Write at screen the name of the mesh that will be converted
  write(*,*)
  write(*,*) "Conversion of ", trim(file_in), " to AFLR3 format in", &
    & " progress..."
  write(*,*)

  ! Ensure that the casefile exists
  max_unit = 98
  i_unit = get_file_unit(max_unit)

  open(unit=i_unit,file=trim(file_in)//'.msh',status='old',iostat=i_err)
  if(i_err /= 0) then
    write(*,*) "The .msh grid ", trim(file_in)," doesn't exist."
    stop
  end if
  close(i_unit)

  ! If the mesh exists, open the mesh file and read it
  open(unit=i_unit,file=trim(file_in)//'.msh',status='old',form='formatted', &
    & iostat=i_err)

  ! Read the string at the beginning of the section: $MeshFormat
  read(i_unit,*) header
    
  ! Read version number, file_type and data_size
  read(i_unit,*) version_gmsh, file_type_gmsh, data_size_gmsh
    
  ! Read the string at the end of the section: $EndMeshFormat
  read(i_unit,*) header

  ! Read the string at the beginning of the section: $Nodes
  read(i_unit,*) header

  ! Read the total number of nodes
  read(i_unit,*) n_nodes_gmsh

  ! Allocate memory for the nodes ID. Gmsh assumes that the node-number must be 
  ! a postive (non-zero) integer BUT it does not assume that node-numbers 
  ! have to form necessarily a dense nor an ordered sequence. 
  allocate(nodes_id_gmsh(n_nodes_gmsh))
  nodes_id_gmsh = 0

  ! Allocate memory for the vertex coordinates
  allocate(nodes_coords_gmsh(1:3,1:n_nodes_gmsh))
  nodes_coords_gmsh = 0.0_wp

  ! Read the ID of the node and its coordinates
  read(i_unit,*) (nodes_id_gmsh(i),nodes_coords_gmsh(1,i), &
    & nodes_coords_gmsh(2,i),nodes_coords_gmsh(3,i),i = 1,n_nodes_gmsh)
    
  ! Read the string at the end of the section: $EndNodes
  read(i_unit,*) header

  ! Read the string at the beginning of the section: $Elements
  read(i_unit,*) header

  ! Read the number of entities.
  ! NOTE: In Gmsh every single geometrical element is an entity: points, lines, 
  ! quads, hexahedron, etc.
  read(i_unit,*) n_entities_gmsh

  ! Read the entities information. This is the first pass just to count the
  ! number of quadrilateral elements with 4 vertices and hexahedrons with 8
  ! vertices.
  cnt_hexas_gmsh = 0
  cnt_quads_gmsh = 0

  do i = 1, n_entities_gmsh
    read(i_unit,*) entity_id_gmsh, entity_type_gmsh, n_entity_tags_gmsh
!    write(*,*) 'entity_id_gmsh', entity_id_gmsh
!    write(*,*) 'entity_type_gmsh', entity_type_gmsh
!    write(*,*) 'n_entity_tags_gmsh', n_entity_tags_gmsh

    ! Entity type = 3 means quadrilateral defined by 4 vertices
    if (entity_type_gmsh == 3) then
      cnt_quads_gmsh = cnt_quads_gmsh + 1

      backspace(i_unit)
      
      ! Trash information
      read(i_unit,*) (trash_int,j = 1, 3 + n_entity_tags_gmsh + &
        & n_nodes_entities(entity_type_gmsh))
      
    ! Entity type = 5 means hexahedrons defined by 8 vertices
    else if (entity_type_gmsh == 5) then
      cnt_hexas_gmsh = cnt_hexas_gmsh + 1

      backspace(i_unit)
      
      ! Trash information
      read(i_unit,*) (trash_int,j = 1, 3 + n_entity_tags_gmsh + &
        & n_nodes_entities(entity_type_gmsh))

    else
      backspace(i_unit)
      
      ! Trash information
      read(i_unit,*) (trash_int,j = 1, 3 + n_entity_tags_gmsh + &
        & n_nodes_entities(entity_type_gmsh))
    end if

  end do
  
  ! Close i_unit 
  close(i_unit)
 
  ! Write at screen the number of quadrilateral faces with 4 vertices and 
  ! hexahedron with 8 vertices
  write(*,*)
  write(*,*) 'Number of hexahedrons with 8 vertices: ', cnt_hexas_gmsh
  write(*,*) 'Number of quadrilateral faces with 4 vertices: ', cnt_quads_gmsh
  write(*,*)

  ! Allocate memory for reading the information of quadrilateral faces defined
  ! by 4 vertices and hexahedron defined by 8 vertices
  allocate(quads_info_gmsh(3+2+n_nodes_entities(3),cnt_quads_gmsh))
  allocate(hexas_info_gmsh(3+2+n_nodes_entities(5),cnt_hexas_gmsh))

  ! Re-open the file and read the information needed. This is the second pass.
  open(unit=i_unit,file=trim(file_in)//'.msh',status='old',form='formatted', &
    & iostat=i_err)

  ! Read the string at the beginning of the section: $MeshFormat
  read(i_unit,*) header
    
  ! Read version number, file_type and data_size
  read(i_unit,*) version_gmsh, file_type_gmsh, data_size_gmsh
    
  ! Read the string at the end of the section: $EndMeshFormat
  read(i_unit,*) header

  ! Read the string at the beginning of the section: $Nodes
  read(i_unit,*) header

  ! Read the total number of nodes (redundant but necessary operation)
  read(i_unit,*) n_nodes_gmsh
  nodes_id_gmsh = 0
  nodes_coords_gmsh = 0.0_wp

  ! Read the ID of the node and its coordinates
  read(i_unit,*) (nodes_id_gmsh(i),nodes_coords_gmsh(1,i), &
    & nodes_coords_gmsh(2,i),nodes_coords_gmsh(3,i),i = 1,n_nodes_gmsh)
    
  ! Read the string at the end of the section: $EndNodes
  read(i_unit,*) header

  ! Read the string at the beginning of the section: $Elements
  read(i_unit,*) header

  ! Read the number of geometrical entities
  read(i_unit,*) n_entities_gmsh

  ! Read the information related to the quadrilateral faces defined by 4 
  ! vertices and hexahedrons defined by 8 vertices
  cnt_hexas_gmsh = 0
  cnt_quads_gmsh = 0

  do i = 1, n_entities_gmsh
    read(i_unit,*) entity_id_gmsh, entity_type_gmsh, n_entity_tags_gmsh

    ! Entity type = 3 means quadrilateral defined by 4 vertices
    if (entity_type_gmsh == 3) then
      cnt_quads_gmsh = cnt_quads_gmsh + 1

      backspace(i_unit)
      
      ! Trash information
      read(i_unit,*) (quads_info_gmsh(j,cnt_quads_gmsh),j = 1, 3 + n_entity_tags_gmsh + &
        & n_nodes_entities(entity_type_gmsh))
      
    ! Entity type = 5 means hexahedrons defined by 8 vertices
    else if (entity_type_gmsh == 5) then
      cnt_hexas_gmsh = cnt_hexas_gmsh + 1
      backspace(i_unit)
      
      ! Trash information
      read(i_unit,*) (hexas_info_gmsh(j,cnt_hexas_gmsh),j = 1, 3 + n_entity_tags_gmsh + &
        & n_nodes_entities(entity_type_gmsh))

    else
      backspace(i_unit)
      
      ! Trash information
      read(i_unit,*) (trash_int,j = 1, n_entity_tags_gmsh + &
        & n_nodes_entities(entity_type_gmsh))
    end if

  end do
  
  ! Close i_unit
  close(i_unit)

  ! Re-oder nodes of the boundary quadrilaterla faces defined by 4 vertices so
  ! that they follow the righ-hand-side rule.

  !do i = 1, cnt_quads_gmsh
  !  if (quads_info_gmsh(4,i) == 300 .or. quads_info_gmsh(4,i) == 301 .or. quads_info_gmsh(4,i) == 302) then
  !    quads_info_gmsh(6,i) = quads_info_gmsh(6,i)
  !    tmp = quads_info_gmsh(7,i)
  !    quads_info_gmsh(7,i) = quads_info_gmsh(9,i)
  !    quads_info_gmsh(8,i) = quads_info_gmsh(8,i)
  !    quads_info_gmsh(9,i) = tmp
  !  end if
  !end do

  ! Define quatities to be written in the output file
  n_nodes = n_nodes_gmsh
  n_surf_trias = 0
  n_surf_quads = cnt_quads_gmsh 
  n_vol_tets = 0
  n_vol_pents_5 = 0
  n_vol_pents_6 = 0
  n_vol_hexs = cnt_hexas_gmsh
  i_unit = get_file_unit(max_unit)
  
  ! Output file name
  output_file_name = trim(file_in) // ".b8.ugrid"
  
  ! Open file
  max_unit = 98
  i_unit = get_file_unit(max_unit)
  open(i_unit,file=output_file_name, status="unknown", access="stream")

  ! Write file, i.e. grid in AFLR3 format.
  write(i_unit) n_nodes
  write(i_unit) n_surf_trias
  write(i_unit) n_surf_quads
  write(i_unit) n_vol_tets
  write(i_unit) n_vol_pents_5
  write(i_unit) n_vol_pents_6
  write(i_unit) n_vol_hexs

  write(i_unit) (nodes_coords_gmsh(1,i),nodes_coords_gmsh(2,i), &
    & nodes_coords_gmsh(3,i),i=1,n_nodes)
  
  write(i_unit) (quads_info_gmsh(6,i),quads_info_gmsh(7,i), &
    & quads_info_gmsh(8,i),quads_info_gmsh(9,i),i=1,n_surf_quads)
  
  write(i_unit) (quads_info_gmsh(4,i),i=1,n_surf_quads)
  
  write(i_unit) (hexas_info_gmsh(6,i),hexas_info_gmsh(7,i), &
    & hexas_info_gmsh(8,i),hexas_info_gmsh(9,i),hexas_info_gmsh(10,i), &
    & hexas_info_gmsh(11,i),hexas_info_gmsh(12,i),hexas_info_gmsh(13,i), &
    & i=1,n_vol_hexs)
  

  ! Close i_unit
  close(i_unit)


  ! Re-oder nodes of the boundary quadrilaterla faces defined by 4 vertices so
  ! that they follow the righ-hand-side rule.

  !do i = 1, cnt_quads_gmsh
  !  if (quads_info_gmsh(4,i) == 300 .or. quads_info_gmsh(4,i) == 301 .or. quads_info_gmsh(4,i) == 302) then
  !    quads_info_gmsh(6,i) = quads_info_gmsh(6,i)
  !    tmp = quads_info_gmsh(7,i)
  !    quads_info_gmsh(7,i) = quads_info_gmsh(9,i)
  !    quads_info_gmsh(8,i) = quads_info_gmsh(8,i)
  !    quads_info_gmsh(9,i) = tmp
  !  end if
  !end do

  ! Define quatities to be written in the output file
  n_nodes = n_nodes_gmsh
  n_surf_trias = 0
  n_surf_quads = cnt_quads_gmsh 
  n_vol_tets = 0
  n_vol_pents_5 = 0
  n_vol_pents_6 = 0
  n_vol_hexs = cnt_hexas_gmsh
  i_unit = get_file_unit(max_unit)
  
  ! Output file name
  output_file_name = trim(file_in) // ".ascii"
  
  ! Open file
  max_unit = 98
  i_unit = get_file_unit(max_unit)
  open(i_unit,file=output_file_name, status="unknown")

  ! Write file, i.e. grid in AFLR3 format.
  write(i_unit,*) n_nodes
  write(i_unit,*) n_surf_trias
  write(i_unit,*) n_surf_quads
  write(i_unit,*) n_vol_tets
  write(i_unit,*) n_vol_pents_5
  write(i_unit,*) n_vol_pents_6
  write(i_unit,*) n_vol_hexs

  write(i_unit,*) (nodes_coords_gmsh(1,i),nodes_coords_gmsh(2,i), &
    & nodes_coords_gmsh(3,i),i=1,n_nodes)
  
  write(i_unit,*) (quads_info_gmsh(6,i),quads_info_gmsh(7,i), &
    & quads_info_gmsh(8,i),quads_info_gmsh(9,i),i=1,n_surf_quads)
  
  write(i_unit,*) (quads_info_gmsh(4,i),i=1,n_surf_quads)
  
  write(i_unit,*) (hexas_info_gmsh(6,i),hexas_info_gmsh(7,i), &
    & hexas_info_gmsh(8,i),hexas_info_gmsh(9,i),hexas_info_gmsh(10,i), &
    & hexas_info_gmsh(11,i),hexas_info_gmsh(12,i),hexas_info_gmsh(13,i), &
    & i=1,n_vol_hexs)
  

  ! Close i_unit
  close(i_unit)



  return
end program gmsh_to_aflr3




