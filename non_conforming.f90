module non_conforming

  use precision_vars
  use initcollocation

  !-- Nothing is implicitely defined
  implicit none

  !-- private subroutines functions etc
  private

  !-- public subroutines functions etc
  public Lagrange_interpolant_basis_1D
  public Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo
  public Vandermonde_1D_monomial
  public Rotate_xione_2_xitwo_and_back
  public Rotate_GL_2_G_and_back_I
  public h_refine_boundary
contains
!==================================================================================================
!
!
! Purpose: This subroutine h refines the mesh at the boundaries
! 
! Comments: Only setup for 3D. Numbering for the splitting: we only consider a equal subvdivision
! of the hex in each direction. The nod numbering convention used here is given
! in the following figure:
!
!                        8-------15--------7
!                       /.                /|
!                      / .               / |
!                     /  .              /  |
!                    /                 /   |
!                   16   .            14   |
!                  /     20          /    19
!                 /      .          /      |
!                /       .         /       |
!               5---------13------6        |       
!               |        4.....11.|........3
!               |       .         |       /
!               |      .          |      /
!               17     .          18    /  zeta ^     
!               |    12           |    10       |    / eta
!               |   .             |   /         |   /
!               |  .              |  /          |  /
!               | .               | /           | /
!               |.                |/            |/
!               1-------9---------2             ----------> xi
!               5---------13------6    6---------14------7   7---------15------8     
!               |        |        |    |        |        |   |        |        |  
!               |        |        |    |        |        |   |        |        |  
!               |        |        |    |        |        |   |        |        |  
!               |        |        |    |        |        |   |        |        |
!               17       21       18   18       22       19  19       23       20
!               |        |        |    |        |        |   |        |        |
!               |        |        |    |        |        |   |        |        |
!               |        |        |    |        |        |   |        |        |
!               |        |        |    |        |        |   |        |        |
!               1-------9---------2    4-------10--------3   3-------11--------4   
!               8---------16------5    4---------11------3   8---------15------7     
!               |        |        |    |        |        |   |        |        |  
!               |        |        |    |        |        |   |        |        |  
!               |        |        |    |        |        |   |        |        |  
!               |        |        |    |        |        |   |        |        |
!               20       24       17   12       25       10  16       26       14
!               |        |        |    |        |        |   |        |        |
!               |        |        |    |        |        |   |        |        |
!               |        |        |    |        |        |   |        |        |
!               |        |        |    |        |        |   |        |        |
!               4-------12---------1   1-------9--------2    5-------13--------6   


!               the node in the middle of the volume is numbered 27
!
! Outputs:
!         
!
!===================================================================================================
   subroutine h_refine_boundary()

     use referencevariables, only : nelems, nfacesperelem, ndim
     use variables, only : ef2e, e_edge2e
     use variables, only : vx_master, ic2nh

     !-- local variables
     integer :: ielem, iface, iedge, ipartner, max_partners, nvertex,&
                element_count, vertex_count
     integer :: nelems_old, nvertex_old, nfaces_old
     integer, allocatable, dimension(:,:,:) :: ef2e_temp                                            !-- (7 ,nfaceperelem, nelements)
     integer, allocatable, dimension(:,:,:,:,:) :: e_edge2e_temp                                    !-- (3,number_of_edges_per_face,max_partners,nfaceperelem,nelems)
     integer, allocatable, dimension(:) :: e_old2e
     integer, allocatable, dimension(:,:) :: ic2nh_temp

     integer :: v1, v2, v3, v4, v5, v6, v7, v8
     integer :: v1n, v2n, v3n, v4n, v5n, v6n, v7n, v8n,&
                v9n, v10n, v11n, v12n, v13n, v14n, v15n, v16n,&
                v17n, v18n, v19n, v20n, v21n, v22n, v23n, v24n,&
                v25n, v26n, v27n
     integer :: element_number_adjoining, faceID_adjoining

     real(wp), dimension (3) :: xyz1, xyz2, xyz3, xyz4,& 
                                xyz5, xyz6, xyz7, xyz8,&
                                xyz9, xyz10, xyz11, xyz12,&
                                xyz13, xyz14, xyz15, xyz16,&
                                xyz17, xyz18, xyz19, xyz20,&
                                xyz21, xyz22, xyz23, xyz24,&
                                xyz25, xyz26, xyz27
     real(wp), allocatable, dimension(:,:) :: vx_master_temp                                         !-- (3,number of vertices)

     logical :: bcelement = .false.

     !-- update the number of faces per element (this is a maximum for one level refinment)
     nfaces_old = nfacesperelem
     nfacesperelem = 4*6

     !-- first we loop over all elements to determine how many boundary elements we have
     !-- and construct a map between the old element numbers and new element numbers
     allocate(e_old2e(nelems))
     element_count = 1
     e_loop_bc1 : do ielem = 1,nelems

       bcelement = .false.
       f_loop_bc1: do iface = 1,nfaces_old

         !-- check if the element has a face with a boundary condition
         f_if_bc1: if (ef2e(1,iface,ielem) < 0) then
           bcelement = .true.
           !-- exit f_loop_bc2
           exit
         endif f_if_bc1

       end do f_loop_bc1
     
       if(bcelement)then
         e_old2e(ielem) = -1000
         !-- update counter
         element_count = element_count+8
       else
         e_old2e(ielem) = element_count

         !-- update counter
         element_count = element_count+1
       endif

     enddo e_loop_bc1
     
     !-- allocate the necessary temporary arrays
     nelems_old = nelems
     nelems = element_count 
     nvertex_old = size(vx_master(1,:)) 
     nvertex = nelems*8
     max_partners = size(e_edge2e(1,1,:,1,1))
                                     
     allocate(vx_master_temp(3,nvertex))
     allocate(ef2e_temp(7,nfacesperelem,nelems))     
     allocate(e_edge2e_temp(3,2**ndim,max_partners,nfacesperelem,nelems))
     allocate(ic2nh_temp(8,nelems))

     !-- store unsplit elements
     e_loop_bc2: do ielem = 1,nelems_old
       bcelement = .false.
       f_loop_bc2: do iface = 1,nfaces_old

         !-- check if the element has a face with a boundary condition
         f_if_bc2: if (ef2e(1,iface,ielem) < 0) then
           !-- split the element and append the information of the new elements at the bottom of the arrays
           bcelement = .true.
           !-- exit f_loop_bc2
           exit
         endif f_if_bc2
       
       end do f_loop_bc2
       if(bcelement)then
         !-- do nothing this is an element that will be split
       else
         !-- first set it equal to original ef2e and e_edge2e
         ef2e_temp(:,:,e_old2e(ielem)) = ef2e(:,:,ielem)
         e_edge2e_temp(:,:,:,:,e_old2e(ielem)) = e_edge2e(:,:,:,:,ielem)
         !-- update the connectivity
         iface_loop: do iface = 1,nfaces_old
           ef2e_temp(2,iface,e_old2e(ielem)) = e_old2e(ef2e(2,iface,ielem))
           iedge_loop: do iedge = 1,4
             ipartner_loop: do ipartner = 1,max_partners
               e_edge2e_temp(2,iedge,ipartner,iface,e_old2e(ielem)) = e_old2e(e_edge2e(2,iedge,ipartner,iface,ielem))
             enddo ipartner_loop
           enddo iedge_loop
         enddo iface_loop
       endif
     enddo e_loop_bc2
    
     !-- loop over the element, split boundary elements, update face conectivity
     !of adjoining elements and populate vx_master_temp, and icn2h_temp
     element_count = 1
     vertex_count = 1
     e_loop_bc3 : do ielem = 1,nelems_old

       bcelement = .false.
       f_loop_bc3: do iface = 1,nfaces_old

         !-- check if the element has a face with a boundary condition
         f_if_bc3: if (ef2e(1,iface,ielem) < 0) then
           !-- split the element and append the information of the new elements at the bottom of the arrays
           bcelement = .true.
           !-- exit f_loop_bc2
           exit
         endif f_if_bc3

       end do f_loop_bc3
       write(*,*)'ielem = ',ielem
       !-- vertex number of the original 8 vertices
       v1 = ic2nh(1,ielem);v2 = ic2nh(2,ielem); v3 = ic2nh(3,ielem); v4 = ic2nh(4,ielem);
       v5 = ic2nh(5,ielem);v6 = ic2nh(6,ielem); v7 = ic2nh(7,ielem); v8 = ic2nh(8,ielem);
write(*,*)'here1'
       !-- spatial locations of the original 8 vertices
       xyz1 = Vx_master(1:3,v1); xyz2 = Vx_master(1:3,v2); xyz3 = Vx_master(1:3,v3); xyz4 = Vx_master(1:3,v4)
       xyz5 = Vx_master(1:3,v5); xyz6 = Vx_master(1:3,v6); xyz7 = Vx_master(1:3,v7); xyz8 = Vx_master(1:3,v8)
write(*,*)'here2'
       !-- set the new vertex numbers
       v1n = vertex_count; v2n = vertex_count+1; v3n = vertex_count+2; v4n = vertex_count+3; 
       v5n =vertex_count+4 ; v6n = vertex_count+5; v7n = vertex_count+6; v8n = vertex_count+7;
write(*,*)'here3'
       !-- if the element is a boundary element split otherwise copy relevant information
       if(bcelement)then
         !-- split the element


         !-- construct new vertices (see figure at top)
         xyz9 = 0.5_wp*(xyz1+xyz2); xyz10 = 0.5_wp*(xyz2+xyz3); xyz11 = 0.5_wp*(xyz3+xyz4); xyz12 = 0.5_wp*(xyz4+xyz1) 
         xyz13 = 0.5_wp*(xyz5+xyz6); xyz14 = 0.5_wp*(xyz6+xyz7); xyz15 = 0.5_wp*(xyz7+xyz8); xyz16 = 0.5_wp*(xyz8+xyz5)
         xyz17 = 0.5_wp*(xyz1+xyz5); xyz18 = 0.5_wp*(xyz2+xyz6); xyz19 = 0.5_wp*(xyz3+xyz7); xyz20 = 0.5_wp*(xyz4+xyz8)
         xyz21 = 0.5_wp*(xyz9+xyz13); xyz22 = 0.5_wp*(xyz10+xyz14); xyz23 = 0.5_wp*(xyz11+xyz15); xyz24 = 0.5_wp*(xyz12+xyz16)
         xyz25 = 0.5_wp*(xyz9+xyz11); xyz26 = 0.5_wp*(xyz13+xyz15);xyz27 = 0.5_wp*(xyz21+xyz23)
write(*,*)'here4'
         !-- set the new vertex numbers
         v9n = vertex_count+8; v10n = vertex_count+9; v11n = vertex_count+10; v12n = vertex_count+11;
         v13n = vertex_count+12; v14n = vertex_count+13; v15n = vertex_count+14; v16n = vertex_count+15;
         v17n = vertex_count+16; v18n = vertex_count+17; v19n = vertex_count+18; v20n = vertex_count+19;
         v21n = vertex_count+20; v22n = vertex_count+21; v23n = vertex_count+22; v24n = vertex_count+23;
         v25n = vertex_count+24; v26n = vertex_count+25; v27n = vertex_count+26;
write(*,*)'here5'
         !-- store the vertices
         vx_master_temp(1:3,v1n) = xyz1; vx_master_temp(1:3,v2n) = xyz2; vx_master_temp(1:3,v3n) = xyz3; 
         vx_master_temp(1:3,v4n) = xyz4;
         vx_master_temp(1:3,v5n) = xyz5; vx_master_temp(1:3,v6n) = xyz6; vx_master_temp(1:3,v7n) = xyz7; 
         vx_master_temp(1:3,v8n) = xyz8;
         vx_master_temp(1:3,v9n) = xyz9;vx_master_temp(1:3,v10n) = xyz10;vx_master_temp(1:3,v11n) = xyz11;
         vx_master_temp(1:3,v12n) = xyz12;  
         vx_master_temp(1:3,v13n) = xyz13;vx_master_temp(1:3,v14n) = xyz14;vx_master_temp(1:3,v15n) = xyz15;
         vx_master_temp(1:3,v16n) = xyz16; 
         vx_master_temp(1:3,v17n) = xyz17;vx_master_temp(1:3,v18n) = xyz18;vx_master_temp(1:3,v19n) = xyz19;
         vx_master_temp(1:3,v20n) = xyz20; 
         vx_master_temp(1:3,v21n) = xyz21;vx_master_temp(1:3,v22n) = xyz22;vx_master_temp(1:3,v23n) = xyz23;
         vx_master_temp(1:3,v24n) = xyz24;
         vx_master_temp(1:3,v25n) = xyz25;vx_master_temp(1:3,v26n) = xyz26;vx_master_temp(1:3,v27n) = xyz27;

write(*,*)'here6'
         !-- store the vertex to element connectivity
         ic2nh_temp(1:8,element_count) = (/v1n,v9n,v25n,v12n,v17n,v21n,v27n,v24n/)
         ic2nh_temp(1:8,element_count+1) = (/V9n,V2n,V10n,V25n,v21n,v18n,v22n,v27n/)
         ic2nh_temp(1:8, element_count+2) = (/v25n,v10n,v3n,v11n,v27n,v22n,v19n,v23n/)
         ic2nh_temp(1:8, element_count+3) = (/v12n,v25n,v11n,v4n,v24n,v27n,v23n,v20n/)
         ic2nh_temp(1:8, element_count+4) = (/v17n,v21n,v27n,v24n,v5n,v13n,v26n,v16n/)
         ic2nh_temp(1:8, element_count+5) = (/v21n,v18n,v22n,v27n,v13n,v6n,v24n,v26n/)
         ic2nh_temp(1:8, element_count+6) = (/v27n,v22n,v19n,v23n,v26n,v14n,v7n,v15n/)
         ic2nh_temp(1:8, element_count+7) = (/v24n,v27n,v23n,v20n,v16n,v26n,v15n,v8n/)
if(ielem.EQ.1)then
write(*,*)'xyz1 = [',xyz1,'];'
write(*,*)'xyz9 = [',xyz9,'];'
write(*,*)'xyz25 = [',xyz25,'];'
write(*,*)'xyz12 = [',xyz12,'];'
write(*,*)'xyz17 = [',xyz17,'];'
write(*,*)'xyz21 = [',xyz21,'];'
write(*,*)'xyz27 = [',xyz27,'];'
write(*,*)'xyz24 = [',xyz24,'];'
write(*,*)'1 = ',vx_master_temp(1:3,ic2nh_temp(1,element_count))
write(*,*)'2 = ',vx_master_temp(1:3,ic2nh_temp(2,element_count))
write(*,*)'3 = ',vx_master_temp(1:3,ic2nh_temp(3,element_count))
write(*,*)'4 = ',vx_master_temp(1:3,ic2nh_temp(4,element_count))
write(*,*)'5 = ',vx_master_temp(1:3,ic2nh_temp(5,element_count))
write(*,*)'6 = ',vx_master_temp(1:3,ic2nh_temp(6,element_count))
write(*,*)'7 = ',vx_master_temp(1:3,ic2nh_temp(7,element_count))
write(*,*)'8 = ',vx_master_temp(1:3,ic2nh_temp(8,element_count))
write(*,*)'element_count = ',element_count
endif

write(*,*)'here7'
if(.false.)then
         !-- update ef2e
         !-- new element 1
         ef2e_temp(:,:,element_count) = ef2e(:,:,ielem)
         !-- new element 1 face 1
         ef2e_temp(1,1,element_count) = ef2e(1,1,ielem)        !--Adjoining element face ID
         ef2e_temp(2,1,element_count) = e_old2e(ef2e(2,1,ielem))!--Adjoining element ID
         !-- new element 1 face 2
         ef2e_temp(1,2,element_count) = ef2e(1,2,ielem)        !--Adjoining element face ID
         ef2e_temp(2,2,element_count) = e_old2e(ef2e(2,1,ielem))!--Adjoining element ID
         !-- new element 1 face 3
         ef2e_temp(1,3,element_count) = 5                      !--Adjoining element face ID
         ef2e_temp(2,3,element_count) = element_count+1        !--Adjoining element ID 
         !-- new element 1 face 4
         ef2e_temp(1,4,element_count) = 2                      !--Adjoining element face ID
         ef2e_temp(2,4,element_count) = element_count+3        !--Adjoining element ID 
         !-- new element 1 face 5
         ef2e_temp(1,5,element_count) = ef2e(1,5,ielem)        !--Adjoining element face ID
         ef2e_temp(2,5,element_count) = e_old2e(ef2e(2,5,ielem))!--Adjoining element ID
         !-- new element 1 face 6
         ef2e_temp(1,6,element_count) = 1                      !--Adjoining element face ID
         ef2e_temp(2,6,element_count) = element_count+4        !--Adjoining element ID
write(*,*)'here8'
         !-- new element 2
         ef2e_temp(:,:,element_count+1) = ef2e(:,:,ielem)
         !-- new element 2 face 1
         ef2e_temp(1,1,element_count+1) = ef2e(1,1,ielem)        !--Adjoining element face ID
         ef2e_temp(2,1,element_count+1) = e_old2e(ef2e(2,1,ielem))!--Adjoining element ID
         !-- new element 2 face 2
         ef2e_temp(1,2,element_count+1) = ef2e(1,2,ielem)        !--Adjoining element face ID
         ef2e_temp(2,2,element_count+1) = e_old2e(ef2e(2,2,ielem))!--Adjoining element ID
         !-- new element 2 face 3
         ef2e_temp(1,3,element_count+1) = ef2e(1,3,ielem)        !--Adjoining element face ID
         ef2e_temp(2,3,element_count+1) = e_old2e(ef2e(2,3,ielem))!--Adjoining element ID 
         !-- new element 2 face 4
         ef2e_temp(1,4,element_count+1) = 2                      !--Adjoining element face ID
         ef2e_temp(2,4,element_count+1) = element_count+2        !--Adjoining element ID 
         !-- new element 2 face 5
         ef2e_temp(1,5,element_count+1) = 3                      !--Adjoining element face ID
         ef2e_temp(2,5,element_count+1) = element_count         !--Adjoining element ID
         !-- new element 2 face 6
         ef2e_temp(1,6,element_count+1) = 1                      !--Adjoining element face ID
         ef2e_temp(2,6,element_count+1) = element_count+5        !--Adjoining element ID
write(*,*)'here9'
         !-- new element 3
         ef2e_temp(:,:,element_count+2) = ef2e(:,:,ielem)
         !-- new element 3 face 1
         ef2e_temp(1,1,element_count+2) = ef2e(1,1,ielem)        !--Adjoining element face ID
         ef2e_temp(2,1,element_count+2) = e_old2e(ef2e(2,1,ielem))!--Adjoining element ID
         !-- new element 3 face 2
         ef2e_temp(1,2,element_count+2) = 4                      !--Adjoining element face ID
         ef2e_temp(2,2,element_count+2) = element_count+1        !--Adjoining element ID
         !-- new element 3 face 3
         ef2e_temp(1,3,element_count+2) = ef2e(1,3,ielem)         !--Adjoining element face ID
         ef2e_temp(2,3,element_count+2) = e_old2e(ef2e(2,3,ielem))!--Adjoining element ID 
         !-- new element 3 face 4
         ef2e_temp(1,4,element_count+2) = ef2e(1,4,ielem)         !--Adjoining element face ID
         ef2e_temp(2,4,element_count+2) = e_old2e(ef2e(2,4,ielem))!--Adjoining element ID
         !-- new element 3 face 5
         ef2e_temp(1,5,element_count+2) = 3                      !--Adjoining element face ID
         ef2e_temp(2,5,element_count+2) = element_count+3        !--Adjoining element ID
         !-- new element 3 face 6
         ef2e_temp(1,6,element_count+2) = 1                      !--Adjoining element face ID
         ef2e_temp(2,6,element_count+2) = element_count+6        !--Adjoining element ID

         !-- new element 4
         ef2e_temp(:,:,element_count+3) = ef2e(:,:,ielem)
         !-- new element 4 face 1
         ef2e_temp(1,1,element_count+3) = ef2e(1,1,ielem)        !--Adjoining element face ID
         ef2e_temp(2,1,element_count+3) = e_old2e(ef2e(2,1,ielem))!--Adjoining element ID
         !-- new element 4 face 2
         ef2e_temp(1,2,element_count+3) = 4                      !--Adjoining element face ID
         ef2e_temp(2,2,element_count+3) = element_count          !--Adjoining element ID
         !-- new element 4 face 3
         ef2e_temp(1,3,element_count+3) = ef2e(1,3,ielem)         !--Adjoining element face ID
         ef2e_temp(2,3,element_count+3) = e_old2e(ef2e(2,3,ielem))!--Adjoining element ID 
         !-- new element 4 face 4
         ef2e_temp(1,4,element_count+3) = ef2e(1,4,ielem)         !--Adjoining element face ID
         ef2e_temp(2,4,element_count+3) = e_old2e(ef2e(2,4,ielem))!--Adjoining element ID 
         !-- new element 4 face 5
         ef2e_temp(1,5,element_count+3) = ef2e(1,5,ielem)         !--Adjoining element face ID
         ef2e_temp(2,5,element_count+3) = e_old2e(ef2e(2,5,ielem))!--Adjoining element ID
         !-- new element 4 face 6
         ef2e_temp(1,6,element_count+3) = 1                      !--Adjoining element face ID
         ef2e_temp(2,6,element_count+3) = element_count+7        !--Adjoining element ID

         !-- new element 5
         ef2e_temp(:,:,element_count+4) = ef2e(:,:,ielem)
         !-- new element 5 face 1
         ef2e_temp(1,1,element_count+4) = 6                       !--Adjoining element face ID
         ef2e_temp(2,1,element_count+4) = element_count           !--Adjoining element ID
         !-- new element 5 face 2
         ef2e_temp(1,2,element_count+4) = ef2e(1,2,ielem)         !--Adjoining element face ID
         ef2e_temp(2,2,element_count+4) = e_old2e(ef2e(2,2,ielem))!--Adjoining element ID
         !-- new element 5 face 3
         ef2e_temp(1,3,element_count+4) = 5                       !--Adjoining element face ID
         ef2e_temp(2,3,element_count+4) = element_count+5         !--Adjoining element ID 
         !-- new element 5 face 4
         ef2e_temp(1,4,element_count+4) = 2                       !--Adjoining element face ID
         ef2e_temp(2,4,element_count+4) = element_count+7         !--Adjoining element ID 
         !-- new element 5 face 5
         ef2e_temp(1,5,element_count+4) = ef2e(1,5,ielem)         !--Adjoining element face ID
         ef2e_temp(2,5,element_count+4) = e_old2e(ef2e(2,5,ielem))!--Adjoining element ID
         !-- new element 5 face 6
         ef2e_temp(1,6,element_count+4) = ef2e(1,6,ielem)         !--Adjoining element face ID
         ef2e_temp(2,6,element_count+4) = e_old2e(ef2e(2,6,ielem))!--Adjoining element ID

         !-- new element 6
         ef2e_temp(:,:,element_count+5) = ef2e(:,:,ielem)
         !-- new element 6 face 1
         ef2e_temp(1,1,element_count+5) = 6                       !--Adjoining element face ID
         ef2e_temp(2,1,element_count+5) = element_count+1         !--Adjoining element ID
         !-- new element 6 face 2
         ef2e_temp(1,2,element_count+5) = ef2e(1,2,ielem)         !--Adjoining element face ID
         ef2e_temp(2,2,element_count+5) = e_old2e(ef2e(2,2,ielem))!--Adjoining element ID
         !-- new element 6 face 3
         ef2e_temp(1,3,element_count+5) = ef2e(1,3,ielem)         !--Adjoining element face ID
         ef2e_temp(2,3,element_count+5) = e_old2e(ef2e(2,3,ielem))!--Adjoining element ID 
         !-- new element 6 face 4
         ef2e_temp(1,4,element_count+5) = 2                       !--Adjoining element face ID
         ef2e_temp(2,4,element_count+5) = element_count+6         !--Adjoining element ID 
         !-- new element 6 face 5
         ef2e_temp(1,5,element_count+5) = 3                       !--Adjoining element face ID
         ef2e_temp(2,5,element_count+5) = element_count+4         !--Adjoining element ID
         !-- new element 6 face 6
         ef2e_temp(1,6,element_count+5) = ef2e(1,6,ielem)         !--Adjoining element face ID
         ef2e_temp(2,6,element_count+5) = e_old2e(ef2e(2,6,ielem))!--Adjoining element ID

         !-- new element 7
         ef2e_temp(:,:,element_count+6) = ef2e(:,:,ielem)
         !-- new element 7 face 1
         ef2e_temp(1,1,element_count+6) = 6                       !--Adjoining element face ID
         ef2e_temp(2,1,element_count+6) = element_count+2         !--Adjoining element ID
         !-- new element 7 face 2
         ef2e_temp(1,2,element_count+6) = 4                       !--Adjoining element face ID
         ef2e_temp(2,2,element_count+6) = element_count+5         !--Adjoining element ID
         !-- new element 7 face 3
         ef2e_temp(1,3,element_count+6) = ef2e(1,3,ielem)         !--Adjoining element face ID
         ef2e_temp(2,3,element_count+6) = e_old2e(ef2e(2,3,ielem))!--Adjoining element ID
         !-- new element 7 face 4
         ef2e_temp(1,4,element_count+6) = ef2e(1,4,ielem)         !--Adjoining element face ID
         ef2e_temp(2,4,element_count+6) = e_old2e(ef2e(2,4,ielem))!--Adjoining element ID
         !-- new element 7 face 5
         ef2e_temp(1,5,element_count+6) = 3                       !--Adjoining element face ID
         ef2e_temp(2,5,element_count+6) = element_count+7         !--Adjoining element ID
         !-- new element 7 face 6
         ef2e_temp(1,6,element_count+6) = ef2e(1,6,ielem)         !--Adjoining element face ID
         ef2e_temp(2,6,element_count+6) = e_old2e(ef2e(2,6,ielem))!--Adjoining element ID

         !-- new element 8
         ef2e_temp(:,:,element_count+7) = ef2e(:,:,ielem)
         !-- new element 8 face 1
         ef2e_temp(1,1,element_count+7) = 6                       !--Adjoining element face ID
         ef2e_temp(2,1,element_count+7) = element_count+3         !--Adjoining element ID
         !-- new element 8 face 2
         ef2e_temp(1,2,element_count+7) = 4                       !--Adjoining element face ID
         ef2e_temp(2,2,element_count+7) = element_count+4         !--Adjoining element ID
         !-- new element 8 face 3
         ef2e_temp(1,3,element_count+7) = 5                       !--Adjoining element face ID
         ef2e_temp(2,3,element_count+7) = element_count+6         !--Adjoining element ID
         !-- new element 8 face 4
         ef2e_temp(1,4,element_count+7) = ef2e(1,4,ielem)         !--Adjoining element face ID
         ef2e_temp(2,4,element_count+7) = e_old2e(ef2e(2,4,ielem))!--Adjoining element ID
         !-- new element 8 face 5
         ef2e_temp(1,5,element_count+7) = ef2e(1,5,ielem)         !--Adjoining element face ID
         ef2e_temp(2,5,element_count+7) = e_old2e(ef2e(2,5,ielem))!--Adjoining element ID
         !-- new element 8 face 6
         ef2e_temp(1,6,element_count+7) = ef2e(1,6,ielem)         !--Adjoining element face ID
         ef2e_temp(2,6,element_count+7) = e_old2e(ef2e(2,6,ielem))!--Adjoining element ID

         !-- update the ef2e of elements touching the new elements 
         !-- old face 1
         element_number_adjoining = e_old2e(ef2e(2,1,ielem))
         faceID_adjoining = ef2e(1,1,ielem)
         ef2e_temp(1,4*(faceID_adjoining-1)+1:4*faceID_adjoining,element_number_adjoining) = 1
         ef2e_temp(2,4*(faceID_adjoining-1)+1:4*faceID_adjoining,element_number_adjoining) = (/1,2,3,4/)

         !-- old face 2
         element_number_adjoining = e_old2e(ef2e(2,2,ielem))
         faceID_adjoining = ef2e(1,2,ielem)
         ef2e_temp(1,4*(faceID_adjoining-1)+1:4*faceID_adjoining,element_number_adjoining) = 2
         ef2e_temp(2,4*(faceID_adjoining-1)+1:4*faceID_adjoining,element_number_adjoining) = (/1,2,6,5/)

         !-- old face 3
         element_number_adjoining = e_old2e(ef2e(2,3,ielem))
         faceID_adjoining = ef2e(1,3,ielem)
         ef2e_temp(1,4*(faceID_adjoining-1)+1:4*faceID_adjoining,element_number_adjoining) = 3
         ef2e_temp(2,4*(faceID_adjoining-1)+1:4*faceID_adjoining,element_number_adjoining) = (/2,3,7,6/)

         !-- old face 4
         element_number_adjoining = e_old2e(ef2e(2,4,ielem))
         faceID_adjoining = ef2e(1,4,ielem)
         ef2e_temp(1,4*(faceID_adjoining-1)+1:4*faceID_adjoining,element_number_adjoining) = 4
         ef2e_temp(2,4*(faceID_adjoining-1)+1:4*faceID_adjoining,element_number_adjoining) = (/4,3,7,6/)

         !-- old face 5
         element_number_adjoining = e_old2e(ef2e(2,5,ielem))
         faceID_adjoining = ef2e(1,5,ielem)
         ef2e_temp(1,4*(faceID_adjoining-1)+1:4*faceID_adjoining,element_number_adjoining) = 5
         ef2e_temp(2,4*(faceID_adjoining-1)+1:4*faceID_adjoining,element_number_adjoining) = (/1,4,8,5/)

         !-- old face 6
         element_number_adjoining = e_old2e(ef2e(2,6,ielem))
         faceID_adjoining = ef2e(1,6,ielem)
         ef2e_temp(1,4*(faceID_adjoining-1)+1:4*faceID_adjoining,element_number_adjoining) = 6
         ef2e_temp(2,4*(faceID_adjoining-1)+1:4*faceID_adjoining,element_number_adjoining) = (/5,6,7,8/)

         !-- update e_edge2e of elements touching the new elements
endif         
         !-- update counters
         vertex_count = vertex_count+26
         element_count = element_count+7
       else
         !-- store old values
         vx_master_temp(1:3,v1n) = xyz1; vx_master_temp(1:3,v2n) = xyz2; vx_master_temp(1:3,v3n) = xyz3; 
         vx_master_temp(1:3,v4n) = xyz4;
         vx_master_temp(1:3,v5n) = xyz5; vx_master_temp(1:3,v6n) = xyz6; vx_master_temp(1:3,v7n) = xyz7; 
         vx_master_temp(1:3,v8n) = xyz8;

         ic2nh_temp(1:8,element_count) = (/v1n,v2n,v3n,v4n,v5n,v6n,v7n,v8n/)

         !-- update counters
         vertex_count = vertex_count+7
         element_count = element_count+1
       endif

     enddo e_loop_bc3
write(*,*)'1 = ',vx_master_temp(1:3,ic2nh_temp(1,1))
write(*,*)'2 = ',vx_master_temp(1:3,ic2nh_temp(2,1))
write(*,*)'3 = ',vx_master_temp(1:3,ic2nh_temp(3,1))
write(*,*)'4 = ',vx_master_temp(1:3,ic2nh_temp(4,1))
write(*,*)'5 = ',vx_master_temp(1:3,ic2nh_temp(5,1))
write(*,*)'6 = ',vx_master_temp(1:3,ic2nh_temp(6,1))
write(*,*)'7 = ',vx_master_temp(1:3,ic2nh_temp(7,1))
write(*,*)' = ',vx_master_temp(1:3,ic2nh_temp(8,1))

write(*,*)'finishes spliting'
     !-- assigne temp arrays to arrays used in main code
     deallocate(ef2e); allocate(ef2e(7,nfacesperelem,nelems))
     ef2e(:,:,:) = ef2e_temp(:,:,:)
write(*,*)'after ef2e'
     deallocate(e_edge2e);allocate(e_edge2e(3,2**ndim,max_partners,nfacesperelem,nelems))
     e_edge2e(:,:,:,:,:) = e_edge2e_temp(:,:,:,:,:)
write(*,*)'after e_edgef2e'
     deallocate(vx_master); allocate(vx_master(3,nelems*8))
write(*,*) 'shape(vx_master) = ',shape(vx_master),'shape(vx_master_temp) = ',shape(vx_master_temp)
     vx_master(:,:) = vx_master_temp(:,:)
write(*,*)'after vx_master'
     deallocate(ic2nh); 
allocate(ic2nh(8,nelems))
write(*,*)'after deallocate allocate'
     ic2nh(:,:) = ic2nh_temp(:,:)
write(*,*)'after ic2nh'
     
     !-- deallocate statements
     deallocate(ef2e_temp)
     deallocate(e_edge2e_temp)
     deallocate(e_old2e)
     deallocate(vx_master_temp)
     deallocate(ic2nh_temp)
   end subroutine h_refine_boundary
!==================================================================================================
!
! Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo()
!
! Purpose: Constructs the Vandermonde matrix from the Lagrange basis functions on the nodal set 
!          xione evaluated at the nodes of xitwo
! 
! Comments: The computation of the interpolation basis functions is based on 
!           "Barycentric Lagrange Interpolation", SIAM REview, Vol 46, No 3 pp. 501-517 
!           (see docs/barycentric.pdf). The notation in this subroutine is consistent with that 
!           in the paper and specifically with the section 3 "An improved Lagrange Formula".
!
!           the barycentric weights, wj, are computed as wj = 1/(Pi_{k=0^{n-1,K!=j}(x_{j}-x_{k})}
!           the Lagrange basis functions are constructed as l_{j} =l*w_{j}/(x-x_{j}) where 
!           l = Pi_{k=0}^{n-1}(x-x_{j})
! Additional documentation: docs/barycentric.pdf
!
! Unit tests: unit_tests/Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo_program 
!             compile using sh compile_Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo
!
! Inputs: nxione (int): number of nodes in the nodal set xione
!         nxitwo (int): number of nodes in the nodal set xitwo
!         xione (real(wp)) size (nxioneX1): nodal set that is used to construct the Lagrange basis functions
!         xitwo (real(wp)) size (nxitwoX1): nodal set at which the Lagrange basis functions will be evaluated at
!
! Outputs:
!         Vandermonde (real(wp)) size (nxitwoXnxione): the above mentioned Vandermonde matrix
!
!===================================================================================================
   subroutine Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo(nxione,nxitwo,xione,xitwo,Vandermonde)

     !-- Nothing is implicitely defined
     implicit none

     !-- input variables
     integer, intent(in)                         :: nxione, nxitwo
     real(wp), intent(in)                        :: xione(nxione), xitwo(nxitwo)
     real(wp), intent(inout)                     :: Vandermonde(nxitwo,nxione)
     
     !-- local variables
     integer                                     :: i, j, k
     real(wp)                                    :: wj(nxione),l
 
     !-- construct the bary centric weights, wj,
     do j = 0,nxione-1 
       wj(j+1) = 1.0_wp
       do k = 0,nxione-1
         if (k.NE.j) then
           wj(j+1) = wj(j+1)*(xione(j+1)-xione(k+1))
         end if
       end do
       wj(j+1) = 1.0_wp/wj(j+1)
     end do
     !-- construct the Vandermonde matrix V(j,k) = l_k(xitwo(j)), where l_k is the Lagrange 
     !-- interpolation basis construted from xione
     do j = 1,nxione
       do k = 1,nxitwo
         !-- construct l
         l = 1.0_wp
         do i = 1, nxione
           l = l*(xitwo(k)-xione(i))
         end do
         !-- check to see if the nodal location on the second set of nodes matches the node 
         !-- location associated with the kth Lagrange basis function 
         if (abs(xitwo(k)-xione(j)) <=2.0_wp*epsilon(1.0_wp)) then
           Vandermonde(k,j) = 1.0_wp
         else
           Vandermonde(k,j)  = l*wj(j)/(xitwo(k)-xione(j))
         end if
       end do
     end do

   end subroutine Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo

!==================================================================================================
!
! Vandermonde_1D_monomial()
!
! Purpose: Constructs the Vandermonde matrix from the monomials on the node set xi 
! 
! Comments: 
!
! Additional documentation:
!
! Unit tests: unit_tests/Vandermonde_1D_monomial.f90 
!             compile using sh unit_tests/comiplie_Vandermonde_1D_monomial.sh
!
! Inputs: n (int): number of nodes in the nodal set xi
!         p (int): the highest degree monomial that will be evaluated
!         xi (real(wp)) size (nX1): nodal on which the monomials will be evaluated
!
! Outputs:
!         Vandermonde (real(wp)) size (nX(p+1)): the above mentioned Vandermonde matrix
!
!===================================================================================================
   subroutine Vandermonde_1D_monomial(n,p,xi,Vandermonde)

     !-- Nothing is implicitely defined
     implicit none

     !-- input variables
     integer, intent(in)                         :: n, p
     real(wp), intent(in)                        :: xi(n)
     real(wp), intent(inout)                     :: Vandermonde(n,(p+1))
     
     !-- local variables
     integer                                     :: i, j
 
     !-- construct the Vandermonde matrix
     do i = 1,n
       do j = 1,p+1
         Vandermonde(i,j) = xi(i)**(j-1)
       end do
     end do
   end subroutine Vandermonde_1D_monomial

!==================================================================================================
!
! Rotate_xione_2_xitwo_and_back()
!
! Purpose: This constructs the interpolants from xione to xitwo and back.
!
! Comments: This is done by constructing the Vandermonde matrix of the Lagrange basis functions, 
!           constructed from xione, evaluated on xitwo, Ixione2xitwo. The interpolant back 
!           is given by Pxione^-1*Ixione2xitwo^T*Pxitwo, where Pxione and Pxitwo are the 
!           diagonal norm matrices on xione and xitwo, respectively. 
!
! Additional documentation
!
! Inputs: nxione (int): number of nodes in the one-dimensional nodal distribution xione
!         nxitwo (int): number of nodes in the one-dimensional nodal distribution xitwo
!         xione (real(wp)) size (nxioneX1): first nodal distribution
!         xitwo (real(wp)) size (nxitwoX1): second nodal distribution
!         Bxione (real(wp) size (nxioneX1): quadrature weights on the first nodal distribution
!         Bxitwo (real(wp)) size (nxitwoX1): quadrature weights on the second nodal distribution
!
! Outputs:
!         xione2xitwo (real(wp)) size(nxioneXnxitwo): interpolation matrix from xione to xitwo
!         xitwo2xione (real(wp)) size(nxitwoXnxione): interpolation matrix from xitwo to xione
!
!===================================================================================================
   subroutine Rotate_xione_2_xitwo_and_back(nxione,nxitwo,xione,xitwo,Bxione,Bxitwo,xione2xitwo,xitwo2xione)

     !-- Nothing is implicitely defined
     implicit none

     !-- input variables
     integer, intent(in)                         :: nxione,nxitwo
     real(wp), intent(in)                        :: xione(nxione),xitwo(nxitwo),Bxione(nxione),Bxitwo(nxitwo)
     real(wp), intent(inout)                     :: xione2xitwo(nxitwo,nxione), xitwo2xione(nxione,nxitwo)

     !-- local variables
     integer                                     :: i,j
     !-- construct the Vandermonde matrix from the Lagrange basis funcitons on xione evaluated at
     !-- the nodal locations of xitwo

     call Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo(nxione,nxitwo,xione,xitwo,xione2xitwo)
    
     !-- construct the interpolant from xitwo to xione
     do i = 1, nxione
       do j = 1,nxitwo
         xitwo2xione(i,j) = 1.0_wp/Bxione(i)*(xione2xitwo(j,i)*Bxitwo(j))
       end do 
     end do 
   end subroutine Rotate_xione_2_xitwo_and_back

!==================================================================================================
!
! Rotate_GL_2_G_and_back_I()
!
! Purpose: This gives back the interpolant from GL to G and back
!
! Comments: The interpolants were solved for in the Maple script maple/Construct_1D_interpolants_Identity.mw
!
! Additional documentation
!
! Inputs: nGL (int): number of nodes in the one-dimensional Gauss Lobatto (GL) nodal distribution
!         nG (int): number of nodes in the one-dimensional Gauss (G) nodal distribution
!         xGL (real(wp)) size (nxioneX1): GL nodal distribution
!         xG (real(wp)) size (nxitwoX1): G nodal distribution
!         BxiGL (real(wp) size (nxioneX1): quadrature weights on the GL nodal distribution
!         BxiG (real(wp)) size (nxitwoX1): quadrature weights on the G nodal distribution
!
! Outputs:
!         IGL2G (real(wp)) size(nGXnGL): interpolation matrix from the GL nodes to the G nodes
!         IG2GL (real(wp)) size(nGLXnG): interpolation matrix from the G nodes to GL nodes
!
!===================================================================================================
   subroutine Rotate_GL_2_G_and_back_I(nGL,nG,xiGL,xiG,BxiGL,BxiG,IGL2G,IG2GL)

     !-- Nothing is implicitely defined
     implicit none

     !-- input variables
     integer, intent(in)                         :: nGL,nG
     real(wp), intent(in)                        :: xiGL(nGL),xiG(nG),BxiGL(nGL),BxiG(nG)
     real(wp), intent(inout)                     :: IGL2G(nG,nGL), IG2GL(nGL,nG)

     !-- local variables
     

    if (nGL==2.AND.nG==2) then
      IGL2G(1,1) = 1.000000000000000000000000000000_wp
      IGL2G(1,2) = 0.000000000000000000000000000000_wp
      IGL2G(2,1) = 0.000000000000000000000000000000_wp
      IGL2G(2,2) = 1.000000000000000000000000000000_wp
      
      IG2GL(1,1) = 1.000000000000000000000000000000_wp
      IG2GL(1,2) = 0.000000000000000000000000000000_wp
      IG2GL(2,1) = 0.000000000000000000000000000000_wp
      IG2GL(2,2) = 1.000000000000000000000000000000_wp
    else if (nGL==3.AND.nG==3) then
      IGL2G(1,1) = 0.343146490609516399718000303683_wp
      IGL2G(1,2) = 1.088303688022450577599852472590_wp
      IGL2G(1,3) = -0.431450178631966977317852776273_wp
      IGL2G(2,1) = 0.430189805014031610999907795369_wp
      IGL2G(2,2) = 0.139620389971936778000184409261_wp
      IGL2G(2,3) = 0.430189805014031610999907795369_wp
      IGL2G(3,1) = -0.431450178631966977317852776273_wp
      IGL2G(3,2) = 1.088303688022450577599852472590_wp
      IGL2G(3,3) = 0.343146490609516399718000303683_wp
      
      IG2GL(1,1) = 0.571910817682527332863333839473_wp
      IG2GL(1,2) = 1.147172813370750962666420787650_wp
      IG2GL(1,3) = -0.719083631053278295529754627123_wp
      IG2GL(2,1) = 0.453459870009354407333271863580_wp
      IG2GL(2,2) = 0.093080259981291185333456272841_wp
      IG2GL(2,3) = 0.453459870009354407333271863580_wp
      IG2GL(3,1) = -0.719083631053278295529754627123_wp
      IG2GL(3,2) = 1.147172813370750962666420787650_wp
      IG2GL(3,3) = 0.571910817682527332863333839473_wp      
    else if (nGL==4.AND.nG==4) then 
      IGL2G(1,1) = 0.670133597061877309375642466273_wp
      IGL2G(1,2) = 0.382690211445752349043404473625_wp
      IGL2G(1,3) = -0.059634895378013858655525747055_wp
      IGL2G(1,4) = 0.006811086870384200236478807157_wp
      IGL2G(2,1) = -0.124994022375115610291702237717_wp
      IGL2G(2,2) = 1.094392949182378686886803453950_wp
      IGL2G(2,3) = 0.011123163321311394153889248044_wp
      IGL2G(2,4) = 0.019477909871425529251009535715_wp
      IGL2G(3,1) = 0.019477909871425529251009535715_wp
      IGL2G(3,2) = 0.011123163321311394153889248052_wp
      IGL2G(3,3) = 1.094392949182378686886803453960_wp
      IGL2G(3,4) = -0.124994022375115610291702237717_wp
      IGL2G(4,1) = 0.006811086870384200236478807157_wp
      IGL2G(4,2) = -0.059634895378013858655525747055_wp
      IGL2G(4,3) = 0.382690211445752349043404473625_wp
      IGL2G(4,4) = 0.670133597061877309375642466273_wp
      
      IG2GL(1,1) = 1.398655311764185408399936481650_wp
      IG2GL(1,2) = -0.489085476472273963755402275747_wp
      IG2GL(1,3) = 0.076214547296997107993785971074_wp
      IG2GL(1,4) = 0.014215617411091447361679823021_wp
      IG2GL(2,1) = 0.159744773085697986032413787639_wp
      IG2GL(2,2) = 0.856443671190025127617787284033_wp
      IG2GL(2,3) = 0.008704700480085614686859237835_wp
      IG2GL(2,4) = -0.024893144755808728337060309506_wp
      IG2GL(3,1) = -0.024893144755808728337060309506_wp
      IG2GL(3,2) = 0.008704700480085614686859237829_wp
      IG2GL(3,3) = 0.856443671190025127617787284040_wp
      IG2GL(3,4) = 0.159744773085697986032413787639_wp
      IG2GL(4,1) = 0.014215617411091447361679823021_wp
      IG2GL(4,2) = 0.076214547296997107993785971074_wp
      IG2GL(4,3) = -0.489085476472273963755402275747_wp
      IG2GL(4,4) = 1.398655311764185408399936481650_wp
    else if (nGL==5.AND.nG==5) then
      IGL2G(1,1) = 0.470502966472490507460316236133_wp
      IGL2G(1,2) = 0.803126567593681665261854299464_wp
      IGL2G(1,3) = -0.491470918189824019737687496142_wp
      IGL2G(1,4) = 0.369914191164023953735669915576_wp
      IGL2G(1,5) = -0.152072807040372106720152955027_wp
      IGL2G(2,1) = 0.071767876791755725608492888402_wp
      IGL2G(2,2) = 0.529448395544771508275461087673_wp
      IGL2G(2,3) = 0.689001782387354883935218360338_wp
      IGL2G(2,4) = -0.492612611092600584063108759494_wp
      IGL2G(2,5) = 0.202394556368718466243936423086_wp
      IGL2G(3,1) = -0.187499999999999999999999999998_wp
      IGL2G(3,2) = 0.437500000000000000000000000001_wp
      IGL2G(3,3) = 0.500000000000000000000000000006_wp
      IGL2G(3,4) = 0.437500000000000000000000000001_wp
      IGL2G(3,5) = -0.187499999999999999999999999998_wp
      IGL2G(4,1) = 0.202394556368718466243936423086_wp
      IGL2G(4,2) = -0.492612611092600584063108759494_wp
      IGL2G(4,3) = 0.689001782387354883935218360338_wp
      IGL2G(4,4) = 0.529448395544771508275461087673_wp
      IGL2G(4,5) = 0.071767876791755725608492888402_wp
      IGL2G(5,1) = -0.152072807040372106720152955027_wp
      IGL2G(5,2) = 0.369914191164023953735669915576_wp
      IGL2G(5,3) = -0.491470918189824019737687496142_wp
      IGL2G(5,4) = 0.803126567593681665261854299464_wp
      IGL2G(5,5) = 0.470502966472490507460316236133_wp
      
      IG2GL(1,1) = 1.114748022560237464834855958820_wp
      IG2GL(1,2) = 0.343501634534003810651995244455_wp
      IG2GL(1,3) = -1.066666666666666666666666666660_wp
      IG2GL(1,4) = 0.968718374310688038674062097476_wp
      IG2GL(1,5) = -0.360301364738262647494246634073_wp
      IG2GL(2,1) = 0.349498057896440617267939081403_wp
      IG2GL(2,2) = 0.465445435697663304036304550688_wp
      IG2GL(2,3) = 0.457142857142857142857142857140_wp
      IG2GL(2,4) = -0.433062586135970603458639415191_wp
      IG2GL(2,5) = 0.160976235399009539297252925968_wp
      IG2GL(3,1) = -0.163747509950278330491814509669_wp
      IG2GL(3,2) = 0.463747509950278330491814509666_wp
      IG2GL(3,3) = 0.400000000000000000000000000005_wp
      IG2GL(3,4) = 0.463747509950278330491814509666_wp
      IG2GL(3,5) = -0.163747509950278330491814509669_wp
      IG2GL(4,1) = 0.160976235399009539297252925968_wp
      IG2GL(4,2) = -0.433062586135970603458639415191_wp
      IG2GL(4,3) = 0.457142857142857142857142857140_wp
      IG2GL(4,4) = 0.465445435697663304036304550688_wp
      IG2GL(4,5) = 0.349498057896440617267939081403_wp
      IG2GL(5,1) = -0.360301364738262647494246634073_wp
      IG2GL(5,2) = 0.968718374310688038674062097476_wp
      IG2GL(5,3) = -1.066666666666666666666666666660_wp
      IG2GL(5,4) = 0.343501634534003810651995244455_wp
      IG2GL(5,5) = 1.114748022560237464834855958820_wp     
    else if (nGL==6.AND.nG==6) then

    else if (nGL==7.AND.nG==7) then
    else if (nGL==8.AND.nG==8) then
    else if (nGL==9.AND.nG==9) then
    else if (nGL==10.AND.nG==10) then
    else if (nGL==11.AND.nG==11) then
    else if (nGL==12.AND.nG==12) then
    else if (nGL==13.AND.nG==13) then
    else if (nGL==15.AND.nG==15) then
    else if (nGL==16.AND.nG==16) then
    else if (nGL==17.AND.nG==17) then
    end if 
   end subroutine Rotate_GL_2_G_and_back_I

!==================================================================================================
!
! Lagrange_interpolant_basis_1D()
!
! Purpose: Evaluates the jth Lagrange basis function, constructed from the nodes xivec at the 
!          point xi.
!
! Additional documentation
!
! Inputs: 
!       n (int): number of nodes in the one-dimensional nodal distribution xivec
!       xivec (real(wp) vector of size [nX1]): one-dimensional nodal distribution
!       xi_eval (real(wp)): location at which the jth Lagrange basis function is evaluated
!       j (int): which Lagrange basis (numbered 0 through n-1)
!
! Outputs:
! 
!        lj_eval (real(wp)): evaluation of the jth Lagrange basis at xi_eval
!
!===================================================================================================
   subroutine Lagrange_interpolant_basis_1D(n,xivec,xi_eval,j,lj_eval)

     !-- arguments
     integer, intent(in)                         :: n,j
     real(wp), intent(in)                        :: xivec(1:n)
     real(wp), intent(in)                        :: xi_eval
     real(wp), intent(inout)                     :: lj_eval
    
     !-- local variables
     integer                                     :: k
     real(wp)                                    :: wj,l,xij

     !-- construct wj
     wj = 1.0_wp
     xij = xivec(j+1)
     do k = 0,n-1
       if (k.NE.j) then
         wj = wj*(xij-xivec(k+1))
       end if
     end do
     wj = 1.0_wp/wj

     !-- construct l
     l = 0.0_wp
     do k = 0,n-1
       l = l*(xi_eval)
     end do
     lj_eval = 1.0_wp
   end subroutine Lagrange_interpolant_basis_1D

!==================================================================================================
!
! Derivative_Lagrange_interpolant_basis_1D()
!
! Purpose:
!
! Additional documentation
!
! Inputs:
!
! Outputs:
!
!===================================================================================================
   !subroutine Derivative_Lagrange_interpolant_basis_1D()
   !end subroutine Derivative_Lagrange_interpolant_basis_1D
end module non_conforming
