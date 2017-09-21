program ssdcsolver_ad
	use precision_vars
	use fileio
	use referencevariables
	use polyinit
	use initgrid
	use mpimod
	use advectiondiffusion
	use timeinteg
	use variables, only: xg, ug, uhat, uold, dudt
	use controlvariables
	implicit none

	integer :: ierr
	! character(120) :: casefile
	integer :: irkstep, itimestep
	! real(wp) :: timemaximum,timestep
	real(wp) :: timelocal, timeglobal

	! PDE indepdent bookkeeping
	! =========================
	! set dimension
	ndim = 2
	! set polynomial order
	npoly = 4
	! initialize 
	ierr = rmapInit(npoly,ndim)
	! read grid file
	call cgnsBCInit()
	casefile = 'test2.uns.4x4'
	call cgnsReadUnstructuredGrid(casefile)
	call init_quad_4()
	! initialize connectivity
	call E2EConnectivity_cgns()
	! initialize nodes
	call calcnodes()
	! calculate metrics
	call calcmetrics()
	! calculate connections
	call calcfacenodeconnectivity()
	! calculate normals
	call calcfacenormals()
	! initialize low storage RK scheme
	call lsrkinit()

	timestep = 1.0e-03_wp

	timemaximum = 5.0e-00_wp

	timeglobal = 0.0_wp
	! =============================
	! End PDE indepdent bookkeeping

	! Wave Equation Solver
	! ====================
	! calculate initial condition
	call ad_calcinitialcondition()
	! allocate memory for semidiscretization
	call ad_initializesemidiscretization()

	itimestep = 0
	do 
		if (timeglobal >= timemaximum) exit
		timestep = min(timestep,timemaximum-timeglobal)
		itimestep = itimestep+1
		uold = 0.0_wp
		do irkstep = 1, rksteps
			! update local time
			timelocal = timeglobal + crk(irkstep)*timestep
			! calculate time derivative approximation
			call ad_calctimederivative(timelocal)
			! integrate in time
			call LSRK(ug, uhat, uold, dudt, irkstep, timestep, &
				nequations, nodesperelem, nelems)
		end do
		timeglobal = timeglobal + timestep
		write(*,100) itimestep, timeglobal
	end do

	100 format('step: ',I10.1,1X,'solution time: ',ES12.5)

	call ad_outputsolution()
	call ad_calcerror(timeglobal)

	! ========================
	! End Wave Equation Solver

end program ssdcsolver_ad
