module controlvariables
  use precision_vars
  use dataTypes
  implicit none

  character(120) :: casefile
  character(120) :: grid_format
  character(120) :: dbfile
  character(120) :: statfile
  character(120) :: runcasename
  character(120) :: solutionfile
  character(120) :: outputfile

  character(120) :: physics                  = 'Nonlinear Navier-Stokes'
  character(120) :: RK_Method                = 'Williamson_Low_Storage_45'
  character(120) :: IMEX_Element             = 'explicit'
  character(120) :: IMEX_penalty             = 'explicit'
  character(120) :: discretization           = 'SSDC'
  character(120) :: flux_entropy_correction  = 'normal'
  character(120) :: Riemann_Diss             = 'Roe'
  character(120) :: Riemann_Diss_BC          = 'Roe'

  character(120) :: Grid_Topology            = 'linear'  !  'cylinder'
  real(wp), dimension(3) :: cylinder_x0, cylinder_x1

  integer   :: itimestep
  integer   :: restart_time_steps = 0
  real(wp)  :: timelocal, timeglobal
  real(wp)  :: timestep, timemaximum
  real(wp)  :: CFL = 1.0_wp
  real(wp)  :: WENO_Bias = 0.75_wp

  real(wp)  :: err_time_L2  = 0.0_wp
  real(wp)  :: err_time_Lf  = 0.0_wp
  real(wp)  :: err_space_L2 = 0.0_wp
  real(wp)  :: err_space_Lf = 0.0_wp
  
  real(wp)  :: err_time_tm1 = 0.0_wp
  real(wp)  :: err_time_tm2 = 0.0_wp

  real(wp)  :: c_1_aero_coeff = 0.0_wp
  real(wp)  :: c_2_aero_coeff = 0.0_wp
  real(wp)  :: c_3_aero_coeff = 0.0_wp
  real(wp), dimension(3)  :: area_sum = 0.0_wp

  real(wp) :: l2_error_bc_no_slip_wall_u_1 = 0.0_wp
  real(wp) :: l2_error_bc_no_slip_wall_u_2 = 0.0_wp
  real(wp) :: l2_error_bc_no_slip_wall_u_3 = 0.0_wp

  real(wp) :: linf_error_bc_no_slip_wall_u_1 = 0.0_wp
  real(wp) :: linf_error_bc_no_slip_wall_u_2 = 0.0_wp
  real(wp) :: linf_error_bc_no_slip_wall_u_3 = 0.0_wp
  
  real(wp) :: linf_error_heat_entropy_flow_wall_bc = 0.0_wp

  real(wp) :: heat_entropy_flow_wall_bc = 0.0_wp

  logical          :: new                            = .true.

  logical          :: time_averaging                 = .true.
  logical          :: filter_solution                = .false.
  logical          :: Entropy_Correction             = .false.
  logical          :: variable_viscosity             = .false.
  logical          :: Dt_by_CFL                      = .false.

  logical          :: write_restart                  = .true.
  character(120)   :: write_restart_dir              = 'restart'
  character(120)   :: write_restart_common_name      = 'restart_p'
  logical          :: write_restart_formatted        = .true.
  integer          :: write_restart_frequency        = 10000
  
  character(120)   :: read_restart_dir               = 'restart_'
  character(120)   :: read_restart_common_name       = 'restart_p'
  character(120)   :: read_restart_time              = 't0.00000000'
  logical          :: read_restart_formatted         = .true.
  logical          :: read_restart_time_averaging    = .false.
  
  logical          :: write_solution                 = .true.
  character(120)   :: write_solution_dir             = 'solution'
  character(120)   :: write_solution_common_name     = 'solution'
  logical          :: write_solution_formatted       = .true.
  integer          :: write_solution_frequency       = 1000
  logical          :: write_aero_coeffs              = .false.
  logical          :: write_dke_dt                   = .false.
  logical          :: write_errors_wall              = .false.
  
  logical          :: verbose                        = .false.
  logical          :: variabletimestep               = .false.
  logical          :: calcjacobian                   = .false.

  logical          :: blprofileread                  = .false.
  integer          :: nblprofiles                    = 1
  
  logical          :: perturb_internal_vertices      = .false.

end module controlvariables
