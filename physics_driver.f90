module physics_driver
  use precision_vars
  implicit none

contains

  subroutine physics_init()
    use controlvariables
    use navierstokes
    use initialize_CSR
    implicit none

    select case (physics)
    case ('Nonlinear Navier-Stokes')
      call Navier_Stokes_init() ! (timeinteg)
      select case (RK_Method)
      case ('Williamson_Low_Storage_45')
      case ('Kraaij_LS_RK_35')
      case ('heun_method')
      case ('IMEX_RK_46')
        !           write(*,*)'before CSR_Get_Pointers'
        !           call CSR_Get_Pointers()
        !           call CSR_Get_Pointers()
        !           write(*,*)'after CSR_Get_Pointers'
        !           call CSR_Combine_Pointers()
        !           write(*,*)'after CSR_Combine_Pointers'
      end select
    case ('linearized Navier-Stokes')
    case ('linearized Euler')
    case default
      write(*,*) 'physics_init'
      write(*,*)'Not a valid set of physical equations'
      write(*,*)'Check physics in setup file'
      write(*,*)'Stopping'
      stop
    end select

  end subroutine physics_init

  subroutine physics_timederivative()
    use controlvariables
    use navierstokes
    implicit none

    select case (physics)
    case ('Nonlinear Navier-Stokes')
      select case (RK_Method)
      case ('Williamson_Low_Storage_45')
        call nse_calc_dudt_LSRK(timelocal)

      case ('Kraaij_LS_RK_35')
        call nse_calc_dudt_LSRK(timelocal)
      
      case ('heun_method') 
        call nse_calc_dudt_LSRK(timelocal)
      
      case ('IMEX_RK_46')
      case default
        write(*,*) 'Invalid time stepping scheme'
        stop
      end select
    case ('linearized Navier-Stokes')
    case ('linearized Euler')
    case default
      write(*,*) 'physics_timederivative'
      write(*,*)'Not a valid set of physical equations'
      write(*,*)'Check physics in setup file'
      write(*,*)'Stopping'
      stop
    end select
  end subroutine physics_timederivative

end module physics_driver
