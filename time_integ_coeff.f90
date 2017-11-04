! This module contains the coefficients of all Runge-Kutta schemes that can be
! use to solve the system of ODEs arsing from the spatial discretization.

module time_integ_coeff

  ! Load modules
  use precision_vars

  ! Nothing is implicitly defined
  implicit none

  ! Global variables of the modules
  real(wp), allocatable, dimension(:,:) :: arkexp, arkimp
  real(wp), allocatable, dimension(:) :: brk, brkh, crk, alsrk
  real(wp), allocatable, dimension(:,:) :: svp
  integer :: narksteps, rksteps
  integer :: current_stage_imex = 1 

contains

  !============================================================================
  
  !============================================================================
  ! ls_rk_initialization - Sets the coefficients of the explicit low-storage
  !                        Runge-Kutta scheme.

  subroutine ls_TVDrk_initialization()

    ! Nothing is implicitly defined
    implicit none

    real(wp) :: b1, b2, b3, b4, b5

    ! Williamson
    ! ==========

    ! Number of stages
    rksteps = 5

    ! Allocate memory for coefficients and initialize them
    allocate(alsrk(rksteps),brk(rksteps),brkh(rksteps),crk(rksteps))
    alsrk = zero
    brk   = zero
    brkh  = zero
    crk   = zero


    alsrk(1) = +0.00000000000000_wp               
    alsrk(2) = -2.60810978953486_wp
    alsrk(3) = -0.08977353434746_wp
    alsrk(4) = -0.60081019321053_wp
    alsrk(5) = -0.72939715170280_wp

    brk(1) = +0.67892607116139_wp
    brk(2) = +0.20654657933371_wp
    brk(3) = +0.27959340290485_wp
    brk(4) = +0.31738259840613_wp
    brk(5) = +0.30319904778284_wp

    ! Weights of the embedded method for the error estimation
    b1 = +0.307294750679102_wp
    b2 = +0.00505716120539947_wp
    b3 = -0.101135470481044_wp
    b4 = +0.731007770927835_wp
    b5 = +0.0577757876687071_wp

    brkh(5) = b5
    brkh(4) = b4 - alsrk(5)*b5
    brkh(3) = b3 - alsrk(4)*b4
    brkh(2) = b2 - alsrk(3)*b3
    brkh(1) = b1 - alsrk(2)*b2

    ! Nodes
    crk(1) = +0.00000000000000_wp
    crk(2) = +0.67892607116139_wp
    crk(3) = +0.34677649493991_wp
    crk(4) = +0.66673359500982_wp
    crk(5) = +0.76590087429032_wp

    return

  end subroutine ls_TVDrk_initialization

  !============================================================================
  
  !============================================================================
  ! ls_rk_initialization - Sets the coefficients of the explicit low-storage
  !                        Runge-Kutta scheme.

  subroutine ls_rk_initialization()

    ! Nothing is implicitly defined
    implicit none

    real(wp) :: b1, b2, b3, b4, b5!, b6
    real(wp) :: AA1, AA2, AA3, AA4, AA5, AA6 
    real(wp) :: BB1, BB2, BB3, BB4, BB5, BB6 
    real(wp) :: BH1, BH2, BH3, BH4, BH5, BH6 

    real(wp), dimension(6) :: b, bh, db

    integer, parameter  :: method = 0

    ! Williamson
    ! ==========

    select case(method)

      case(0)

        rksteps = 5                ! Number of stages

        ! Allocate memory for coefficients and initialize them
        allocate(alsrk(rksteps),brk(rksteps),brkh(rksteps),crk(rksteps))
        alsrk = zero
        brk   = zero
        brkh  = zero
        crk   = zero

        ! Stage coefficients
        alsrk(1) =  0.0_wp
        alsrk(2) = -0.4801594388478_wp
        alsrk(3) = -1.4042471952_wp
        alsrk(4) = -2.016477077503_wp
        alsrk(5) = -1.056444269767_wp

        ! Weights
        brk(1) = 0.1028639988105_wp
        brk(2) = 0.7408540575767_wp
        brk(3) = 0.7426530946684_wp
        brk(4) = 0.4694937902358_wp
        brk(5) = 0.1881733382888_wp

        ! Helping variables
        b1 = 1150667._wp/35155867._wp
        b2 = 12453995._wp/50501979._wp
        b3 = 6259373._wp/22243429._wp
        b4 = 18446371._wp/66784475._wp
        b5 = 9499292._wp/58258297._wp
  
        ! Weights of the embedded method for the error estimation
        brkh(5) = b5
        brkh(4) = b4 - alsrk(5)*b5
        brkh(3) = b3 - alsrk(4)*b4
        brkh(2) = b2 - alsrk(3)*b3
        brkh(1) = b1 - alsrk(2)*b2

        ! Nodes
        crk(1) = 0._wp
        crk(2) = 0.1028639988105_wp
        crk(3) = 0.487989987833_wp
        crk(4) = 0.6885177231562_wp
        crk(5) = 0.9023816453077_wp

      case(1)

        rksteps = 6                ! Number of stages

        ! Allocate memory for coefficients and initialize them
        allocate(alsrk(rksteps),brk(rksteps),brkh(rksteps),crk(rksteps))
        alsrk = zero
        brk   = zero
        brkh  = zero
        crk   = zero

        ! Stage coefficients
        alsrk(2) = -56022528.0_wp/ 150050353.0_wp
        alsrk(3) = -96646469.0_wp/ 120630563.0_wp
        alsrk(4) = -144992129.0_wp/ 98109555.0_wp
        alsrk(5) = -79666811.0_wp/ 37868774.0_wp
        alsrk(6) = -85387295.0_wp/ 44036756.0_wp
        AA1 = alsrk(1) ; 
        AA2 = alsrk(2) ; 
        AA3 = alsrk(3) ; 
        AA4 = alsrk(4) ; 
        AA5 = alsrk(5) ; 
        AA6 = alsrk(6) ; 

        ! Weights of main mathod
        brk(1)  = 10764601.0_wp/ 113944427.0_wp
        brk(2)  = 36150543.0_wp/ 139452415.0_wp
        brk(3)  = 182008983.0_wp/ 504724679.0_wp
        brk(4)  = 92067757.0_wp/ 111839021.0_wp
        brk(5)  = 46859183.0_wp/ 63300472.0_wp
        brk(6)  = 47126472.0_wp/ 380443417.0_wp
        BB1 = brk(1) ; 
        BB2 = brk(2) ; 
        BB3 = brk(3) ; 
        BB4 = brk(4) ; 
        BB5 = brk(5) ; 
        BB6 = brk(6) ; 

        !  Butcher weights for b_j
        b(1)  = BB1 + AA2*BB2 + AA2*AA3*BB3 + AA2*AA3*AA4*BB4 + AA2*AA3*AA4*AA5*BB5 + AA2*AA3*AA4*AA5*AA6*BB6;
        b(2)  = BB2 + AA3*BB3 + AA3*AA4*BB4 + AA3*AA4*AA5*BB5 + AA3*AA4*AA5*AA6*BB6;
        b(3)  = BB3 + AA4*BB4 + AA4*AA5*BB5 + AA4*AA5*AA6*BB6;
        b(4)  = BB4 + AA5*BB5 + AA5*AA6*BB6;
        b(5)  = BB5 + AA6*BB6;
        b(6)  = BB6;

        ! db_j = b_j - bh_j
        db(1) = -5062502208571.0_wp/10831071574082.0_wp
        db(2) =  7350684609029.0_wp/8562057514121.0_wp
        db(3) = -18221243228277.0_wp/40554070092353.0_wp
        db(4) =  129818982301.0_wp/4807233560427.0_wp
        db(5) =  253026190555.0_wp/10259350328147.0_wp
        db(6) =  15875538881.0_wp/2432482339802.0_wp

        ! Weights of embedded mathod
        bh(:) = b(:) - db(:)

        crk(1)  = 0.0_wp
        crk(2)  = 10764601.0_wp/113944427.0_wp
        crk(3)  = 1122593674877.0_wp/4369462009356.0_wp
        crk(4)  = 4914898956286.0_wp/11260214860537.0_wp
        crk(5)  = 2210250771543.0_wp/3380123205067.0_wp
        crk(6)  = 34536188621667.0_wp/ 35138252455445.0_wp

        BH6  = bh(6);
        BH5  = bh(5) - AA6*BH6;
        BH4  = bh(4) - AA5*BH5 - AA5*AA6*BH6;
        BH3  = bh(3) - AA4*BH4 - AA4*AA5*BH5 - AA4*AA5*AA6*BH6;
        BH2  = bh(2) - AA3*BH3 - AA3*AA4*BH4 - AA3*AA4*AA5*BH5 - AA3*AA4*AA5*AA6*BH6;
        BH1  = bh(1) - AA2*BH2 - AA2*AA3*BH3 - AA2*AA3*AA4*BH4 - AA2*AA3*AA4*AA5*BH5 - AA2*AA3*AA4*AA5*AA6*BH6;

        brkh(1)  = BH1
        brkh(2)  = BH2
        brkh(3)  = BH3
        brkh(4)  = BH4
        brkh(5)  = BH5
        brkh(6)  = BH6

    end select

    return
  end subroutine ls_rk_initialization

  !============================================================================
  
  !============================================================================
  ! kraaijemanger_initialization - Sets the coefficients of the explicit 
  !                                Kraaijemanger 54 Runge-Kutta scheme.
  ! Hans Kraaijemanger  RK54
  ! Max stability coefficient "r"  r <= 1.5081800491898379228

  subroutine kraaijemanger_initialization

    ! Nothing is implicitly defined
    implicit none

    !real(wp) :: al32
    !real(wp) :: al42, al43
    !real(wp) :: al52, al53, al54
    !real(wp) :: al62, al63, al64, al65

    ! Number of stages
    narksteps = 5

    ! Allocate memory for the coefficients and initialize them
    allocate(arkexp(narksteps,narksteps),brk(narksteps),brkh(narksteps), &
      & crk(narksteps))
    arkexp = zero
    brk = zero
    brkh = zero
    crk = zero

    ! Stage coefficients
    arkexp(2,1) = 0.3917522265718890583279931517854343_wp
    arkexp(3,1) = 0.21766909626116921035766572386831129_wp
    arkexp(3,2) = 0.36841059305037202075016687846678052_wp
    arkexp(4,1) = 0.082692086657810754405152309996115072_wp
    arkexp(4,2) = 0.13995850219189573937901744272469101_wp
    arkexp(4,3) = 0.25189177427169263983636985780844332_wp
    arkexp(5,1) = 0.067966283637114963235790587691024752_wp
    arkexp(5,2) = 0.11503469850463199466939672337939896_wp
    arkexp(5,3) = 0.20703489859738471850986385308584489_wp
    arkexp(5,4) = 0.54497475022851992203743357849557077_wp
    
    ! Weights
    brk(1) = 0.1468118760847864495633257353492341_wp
    brk(2) = 0.24848290944497614757118961041719405_wp
    brk(3) = 0.10425883033198029566679642356126373_wp
    brk(4) = 0.27443890090134945680513992422870496_wp
    brk(5) = 0.22600748323690765039354833410343016_wp

    ! Weights of the embedded method for the error estimation
    ! (3rd order embedded scheme)
    brkh(1) = 0.20482128441985996296373209258976
    brkh(2) = 0
    brkh(3) = 0.478448536408634339291743883212572
    brkh(4) = 0.1662547879827892615039004382002
    brkh(5) = 0.15047539118871643624062358599751
    
    ! Nodes
    crk(2) = 0.3917522265718890583279931517854343_wp
    crk(3) = 0.5860796893115412311078326023350918_wp
    crk(4) = 0.47454236312139913362053961052924941_wp
    crk(5) = 0.93501063096765159845248474265183937_wp

    return 
  end subroutine kraaijemanger_initialization

  !============================================================================
  
  !============================================================================
  ! rk_imex_46_initialization - Sets the coefficients of the implicit-explicit 
  !                             46 Runge-Kutta scheme.

  subroutine rk_imex_46_initialization

    ! Nothing is implicitly defined
    real(wp) :: al32
    real(wp) :: al42, al43
    real(wp) :: al52, al53, al54
    real(wp) :: al62, al63, al64, al65

    ! Number of stages
    narksteps = 6

    allocate(arkexp(narksteps,narksteps),arkimp(narksteps,narksteps), &
      & brk(narksteps),brkh(narksteps),crk(narksteps),svp(narksteps,narksteps))
    arkexp = zero
    arkimp = zero
    brk = zero
    brkh = zero
    crk = zero
    svp = zero

    ! Stage coefficients of the explicit portion of the RK method
    arkexp(2,1) = 1._wp/2._wp
    arkexp(3,1) = 13861._wp/62500._wp
    arkexp(3,2) = 6889._wp/62500._wp
    arkexp(4,1) =-116923316275._wp/2393684061468._wp
    arkexp(4,2) =-2731218467317._wp/15368042101831._wp
    arkexp(4,3) = 9408046702089._wp/11113171139209._wp
    arkexp(5,1) =-451086348788._wp/2902428689909._wp
    arkexp(5,2) =-2682348792572._wp/7519795681897._wp
    arkexp(5,3) = 12662868775082._wp/11960479115383._wp
    arkexp(5,4) = 3355817975965._wp/11060851509271._wp
    arkexp(6,1) = 647845179188._wp/3216320057751._wp
    arkexp(6,2) = 73281519250._wp/8382639484533._wp
    arkexp(6,3) = 552539513391._wp/3454668386233._wp
    arkexp(6,4) = 3354512671639._wp/8306763924573._wp
    arkexp(6,5) = 4040._wp/17871._wp

    ! Stage coefficients of the implicit portion of the RK method
    arkimp(2,1) = 1._wp/4._wp
    arkimp(2,2) = 1._wp/4._wp
    arkimp(3,1) = 8611._wp/62500._wp
    arkimp(3,2) =-1743._wp/31250._wp
    arkimp(3,3) = 1._wp/4._wp
    arkimp(4,1) = 5012029._wp/34652500._wp
    arkimp(4,2) =-654441._wp/2922500._wp
    arkimp(4,3) = 174375._wp/388108._wp
    arkimp(4,4) = 1._wp/4._wp
    arkimp(5,1) = 15267082809._wp/155376265600._wp
    arkimp(5,2) =-71443401._wp/120774400._wp
    arkimp(5,3) = 730878875._wp/902184768._wp
    arkimp(5,4) = 2285395._wp/8070912._wp
    arkimp(5,5) = 1._wp/4._wp
    arkimp(6,1) = 82889._wp/524892._wp
    arkimp(6,2) = 0._wp
    arkimp(6,3) = 15625._wp/83664._wp
    arkimp(6,4) = 69875._wp/102672._wp
    arkimp(6,5) =-2260._wp/8211._wp
    arkimp(6,6) = 1._wp/4._wp

    ! Weigths 
    brk(1) = 82889._wp/524892._wp
    brk(2) = 0._wp
    brk(3) = 15625._wp/83664._wp
    brk(4) = 69875._wp/102672._wp
    brk(5) =-2260._wp/8211._wp
    brk(6) = 1._wp/4._wp

    ! Weights of the embedded method for the error estimation
    brkh(1) = 4586570599._wp/29645900160._wp
    brkh(2) = 0._wp
    brkh(3) = 178811875._wp/945068544._wp
    brkh(4) = 814220225._wp/1159782912._wp
    brkh(5) = -3700637._wp/11593932._wp
    brkh(6) = 61727._wp/225920._wp

    ! Nodes
    crk(1) = 0._wp
    crk(2) = 1._wp/2._wp
    crk(3) = 83._wp/250._wp
    crk(4) = 31._wp/50._wp
    crk(5) = 17._wp/20._wp
    crk(6) = 1._wp

    ! Coefficients of the stage value predictor (SVP)
    ! These are the ``Hairer'' coefficients for this scheme.
    ! I call these Hairer coefficients because the idea come from HWII.
    ! However, they are derived for this specific scheme.
    al32 = 83._wp/125._wp           
    al42 = 372._wp/175._wp          
    al43 = -775._wp/581._wp         
    al52 = 52003._wp/16272._wp
    al53 = -65639125._wp/16206912._wp
    al54 = 5825645._wp/6053184._wp
    al62 = -155._wp/567._wp
    al63 = -113125._wp/141183._wp
    al64 = 43300._wp/173259._wp
    al65 = 36160._wp/24633._wp

    svp(3,1) = (arkimp(2,1)*al32) !*Fn
    svp(3,2) = (arkimp(2,2)*al32) !*F2

    svp(4,1) = (arkimp(2,1)*al42 + arkimp(3,1)*al43) !*Fn
    svp(4,2) = (arkimp(2,2)*al42 + arkimp(3,2)*al43) !*F2
    svp(4,3) = (arkimp(3,3)*al43) !*F3

    svp(5,1) = (arkimp(2,1)*al52 + arkimp(3,1)*al53 + arkimp(4,1)*al54) !*Fn
    svp(5,2) = (arkimp(2,2)*al52 + arkimp(3,2)*al53 + arkimp(4,2)*al54) !*F2
    svp(5,3) = (arkimp(3,3)*al53 + arkimp(4,3)*al54) !*F3
    svp(5,4) = (arkimp(4,4)*al54) !*F4

    svp(6,1) = (arkimp(2,1)*al62 + arkimp(3,1)*al63 + arkimp(4,1)*al64 + arkimp(5,1)*al65) !*Fn
    svp(6,2) = (arkimp(2,2)*al62 + arkimp(3,2)*al63 + arkimp(4,2)*al64 + arkimp(5,2)*al65) !*F2
    svp(6,3) = (arkimp(3,3)*al63 + arkimp(4,3)*al64 + arkimp(5,3)*al65) !*F3
    svp(6,4) = (arkimp(4,4)*al64 + arkimp(5,4)*al65) !*F4
    svp(6,5) = (arkimp(5,5)*al65) !*F5

    return
  end subroutine rk_imex_46_initialization

  !============================================================================

end module time_integ_coeff
