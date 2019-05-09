    !!!!!module hod!!!!!
    module hod
        use settings
        use libnpcf
        implicit none
        type, bind(C)::float3
                real(kind=4)::x,y,z
        end type
        type(npcf) :: corr
        type(float3), dimension(:), allocatable :: pos_g,pos_h,pos_s,triangles,temp
        real(mcp), dimension(:), allocatable :: mass
        real(mcp), dimension(:), allocatable :: twopcf_data,threepcf_data
        real(mcp), dimension(:), allocatable :: num_sat, ind_sat
        real(mcp), dimension(:,:), allocatable :: cov2,invcov2,cov3,invcov3
        integer, dimension(:), allocatable :: work2, work3
        integer :: timesRans,numShells, numTriangles,dummy
        integer :: nhalos,nsates
        real(mcp) :: v_r, r_max, r_min, V_box
        real(mcp), dimension(:), allocatable :: threePoint,twoPoint,shells
        logical :: use_2pcf,use_3pcf
    contains
        subroutine init_random_seed()

        INTEGER :: i, n, clock
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed

        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))

        CALL SYSTEM_CLOCK(COUNT=clock)

        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        CALL RANDOM_SEED(PUT = seed)

        DEALLOCATE(seed)
        end
        SUBROUTINE DVERT(V,LV,N,W)
      INTEGER LV
      REAL*8 V(LV,1),S,T
      INTEGER W(1),I,J,K,L,M,N,P
      IF ( N .EQ. 1 ) GOTO 110
      L = 0
      M = 1
10    IF ( L .EQ. N ) GOTO 90
      K = L
      L = M
      M = M + 1
!C     ---------------------------------------
!C     |*** FIND PIVOT AND START ROW SWAP ***|
!C     ---------------------------------------
      P = L
      IF ( M .GT. N ) GOTO 30
      S = ABS(V(L,L))
      DO 20 I = M,N
           T = ABS(V(I,L))
           IF ( T .LE. S ) GOTO 20
           P = I
           S = T
20    CONTINUE
      W(L) = P
30    S = V(P,L)
      V(P,L) = V(L,L)
      IF ( S .EQ. 0. ) GOTO 120
!C     -----------------------------
!C     |*** COMPUTE MULTIPLIERS ***|
!C     -----------------------------
      V(L,L) = -1.
      S = 1./S
      DO 40 I = 1,N
40         V(I,L) = -S*V(I,L)
      J = L
50    J = J + 1
      IF ( J .GT. N ) J = 1
      IF ( J .EQ. L ) GOTO 10
      T = V(P,J)
      V(P,J) = V(L,J)
      V(L,J) = T
      IF ( T .EQ. 0. ) GOTO 50
!C     ------------------------------
!C     |*** ELIMINATE BY COLUMNS ***|
!C     ------------------------------
      IF ( K .EQ. 0 ) GOTO 70
      DO 60 I = 1,K
60         V(I,J) = V(I,J) + T*V(I,L)
70    V(L,J) = S*T
      IF ( M .GT. N ) GOTO 50
      DO 80 I = M,N
80         V(I,J) = V(I,J) + T*V(I,L)
      GOTO 50
!C     -----------------------
!C     |*** PIVOT COLUMNS ***|
!C     -----------------------
90    L = W(K)
      DO 100 I = 1,N
           T = V(I,L)
           V(I,L) = V(I,K)
100        V(I,K) = T
      K = K - 1
      IF ( K .GT. 0 ) GOTO 90
      RETURN
110   IF ( V(1,1) .EQ. 0. ) GOTO 120
      V(1,1) = 1./V(1,1)
      RETURN
120   WRITE(6,*) 'ERROR: MATRIX HAS NO INVERSE'
      STOP
      END
      
        function poisson_inphase(mean)
        real(mcp) :: mean,s,p,u
        integer :: poisson_inphase,x
        x=0
        s=dexp(-mean)
        p=dexp(-mean)
        call RANDOM_NUMBER(u)
        do while (u>s)
                x=x+1
                p=p*mean/x
                s=s+p
        end do
        poisson_inphase=x 
        
        end function

    end module hod
    

    !!!!!module hod end!
    module CalcLike
    use DataLikelihoodList
    use GeneralTypes
    use BaseParameters
    use MatrixUtils
    use MiscUtils
    implicit none
    private

    Type, extends(TConfigClass) :: TLikeCalculator
        real(mcp) :: Temperature =1
        logical :: test_likelihood= .false.
        logical :: timing = .false.
        real(mcp), allocatable :: test_cov_matrix(:,:)
    contains
    procedure :: AddLike
    procedure :: AddLikeTemp
    procedure :: GetLogLikeBounds
    procedure :: GetLogPriors
    procedure :: GetLogLike
    procedure :: GetLogLikeMain
    procedure :: ReadParams => TLikeCalculator_ReadParams
    procedure :: TestLikelihoodFunction => TLikeCalculator_TestLikelihoodFunction
    procedure :: WritePerformanceStats => TLikeCalculator_WritePerformanceStats
    procedure :: WriteParamsHumanText_unit => TLikeCalculator_WriteParamsHumanText
    procedure :: WriteParamsHumanText_file => TLikeCalculator_WriteParamsHumanTextFile
    procedure :: WriteParamPointTextData => TLikeCalculator_WriteParamPointTextData
    generic :: WriteParamsHumanText => WriteParamsHumanText_unit, WriteParamsHumanText_file
    end type TLikeCalculator

    type, extends(TLikeCalculator) :: TGenericLikeCalculator
    contains
    procedure :: GetLogLikeMain => Generic_GetLogLikeMain
    end type TGenericLikeCalculator

    type, extends(TLikeCalculator) :: TTheoryLikeCalculator
        class(TTheoryParams), allocatable :: TheoryParams
        class(TCalculationAtParamPoint), pointer :: Params
        logical changeMask(max_num_params)
        logical SlowChanged, SemiSlowchanged
    contains
    procedure :: GetLogLikeMain => TheoryLike_GetLogLikeMain
    procedure :: GetLogLikePost => TheoryLike_GetLogLikePost
    procedure :: SetTheoryParams => TheoryLike_SetTheoryParams
    procedure :: CheckPriorCuts => TheoryLike_CheckPriorCuts
    procedure :: CalculateRequiredTheoryChanges =>TheoryLike_CalculateRequiredTheoryChanges
    procedure :: GetTheoryForLike=>TheoryLike_GetTheoryForLike
    procedure :: GetTheoryForImportance=>TheoryLike_GetTheoryForImportance
    procedure :: GetLogLikeWithTheorySet => TheoryLike_LogLikeWithTheorySet
    procedure :: UpdateTheoryForLikelihoods => TheoryLike_UpdateTheoryForLikelihoods
    procedure :: SetNewTheoryResults => TheoryLike_SetNewTheoryResults
    procedure :: TestLikelihoodFunction => TheoryLike_TestLikelihoodFunction
    procedure :: WriteParamsHumanText_unit => TTheoryLike_WriteParamsHumanText
    procedure :: WriteParamPointTextData => TTheoryLike_WriteParamPointTextData
    end type TTheoryLikeCalculator

    type, extends(TCheckpointable) :: TLikelihoodUser
        class(TLikeCalculator), pointer :: LikeCalculator => null()
    end type

    type, extends(TCheckpointable) :: TTheoryLikelihoodUser
        class(TTheoryLikeCalculator), pointer :: LikeCalculator => null()
    end type

    public TLikeCalculator, TGenericLikeCalculator, TTheoryLikeCalculator, TLikelihoodUser, TTheoryLikelihoodUser
    contains

    subroutine AddLike(this, CurrentLike, LikeToAdd)
    class(TLikeCalculator) :: this
    real(mcp), intent(in) :: LikeToAdd
    real(mcp) CurrentLike

    if (CurrentLike/=LogZero) then
        if (LikeToAdd == logZero) then
            CurrentLike = LogZero
        else
            CurrentLike = CurrentLike + LikeToAdd
        end if
    end if
    end subroutine AddLike

    subroutine AddLikeTemp(this, CurrentLike, LikeToAdd)
    class(TLikeCalculator) :: this
    real(mcp), intent(in) :: LikeToAdd
    real(mcp) CurrentLike

    if (CurrentLike/=LogZero) then
        if (LikeToAdd == logZero) then
            CurrentLike = LogZero
        else
            CurrentLike = CurrentLike + LikeToAdd/this%Temperature
        end if
    end if
    end subroutine AddLikeTemp


    function GetLogLikeBounds(this,Params)
    class(TLikeCalculator) :: this
    class(TCalculationAtParamPoint):: Params
    real(mcp) :: GetLogLikeBounds

    if (any(Params%P(:num_params) > BaseParams%PMax(:num_params)) .or. &
        & any(Params%P(:num_params) < BaseParams%PMin(:num_params))) then
    GetLogLikeBounds = logZero
    else
        GetLogLikeBounds=0
    end if

    end function GetLogLikeBounds

    function GetLogPriors(this, P) result(logLike)
    class(TLikeCalculator) :: this
    integer i
    real(mcp), intent(in) :: P(num_params)
    real(mcp) logLike

    logLike=0
    do i=1,num_params
        if ((BaseParams%varying(i) .or. BaseParams%include_fixed_parameter_priors) &
            .and. BaseParams%GaussPriors%std(i)/=0) then
        logLike = logLike + ((P(i)-BaseParams%GaussPriors%mean(i))/BaseParams%GaussPriors%std(i))**2
        end if
    end do

    do i= 1, size(BaseParams%LinearCombinations)
        associate(Comb => BaseParams%LinearCombinations(i))
            if (Comb%std/=0) then
                logLike = logLike + ((dot_product(Comb%Combination,P) -Comb%mean)/Comb%std)**2
            end if
        end associate
    end do
    logLike=logLike/2

    end function GetLogPriors

    function GetLogLike(this, Params) !Get -Ln(Likelihood) for chains
    class(TLikeCalculator) :: this
    class(TCalculationAtParamPoint) :: Params
    real(mcp) GetLogLike

    GetLogLike = this%GetLogLikeBounds(Params)
    if (GetLogLike==LogZero) return
    if (this%test_likelihood) then
        call this%AddLikeTemp(GetLogLike,this%TestLikelihoodFunction(Params))
    else
        call this%AddLikeTemp(GetLogLike,this%GetLogLikeMain(Params))
    end if
    if (GetLogLike==LogZero) return
    call this%AddLikeTemp(GetLogLike,this%getLogPriors(Params%P))

    end function GetLogLike

    function GetLogLikeMain(this, Params) result(LogLike)!Get -Ln(Likelihood) for chains
    class(TLikeCalculator) :: this
    class(TCalculationAtParamPoint) :: Params
    real(mcp) :: LogLike

    LogLike = LogZero !error - must be overridden

    end function GetLogLikeMain

    subroutine TLikeCalculator_ReadParams(this, Ini)
    class(TLikeCalculator) :: this
    class(TSettingIni) :: Ini
    character(LEN=:), pointer :: covMatrix

    this%test_likelihood = Ini%Read_Logical('test_likelihood', .false.)
    if (this%test_likelihood) then
        print *,'** Using test Gaussian likelihood from covariance + hard priors **'
        covMatrix => Ini%Read_String('test_covariance')
        if (covMatrix/='') then
            allocate(this%test_cov_matrix(num_params_used, num_params_used))
            call BaseParams%ReadSetCovMatrix(covMatrix, this%test_cov_matrix)
        end if
    end if
    call Ini%Read('temperature',this%Temperature)

    end subroutine TLikeCalculator_ReadParams

    function TLikeCalculator_TestLikelihoodFunction(this,Params) result(LogLike)
    class(TLikeCalculator) :: this
    class(TCalculationAtParamPoint) Params
    real(mcp) :: LogLike
    real(mcp), allocatable, save :: covInv(:,:)
    real(mcp) X(num_params_used)

    if (.not. allocated(covInv)) then
        allocate(covInv(num_params_used,num_params_used))
        if (.not. allocated(this%test_cov_matrix)) then
            covInv = BaseParams%covariance_estimate
        else
            covInv = this%test_cov_matrix
        endif
        call Matrix_Inverse(covInv)
    end if
    X = Params%P(params_used) - BaseParams%Center(params_used)
    LogLike = dot_product(X, matmul(covInv, X))/2

    end function TLikeCalculator_TestLikelihoodFunction

    subroutine TLikeCalculator_WritePerformanceStats(this, unit)
    class(TLikeCalculator) :: this
    integer, intent(in) :: unit

    end subroutine TLikeCalculator_WritePerformanceStats


    subroutine TLikeCalculator_WriteParamsHumanText(this, aunit, P, LogLike, weight)
    class(TLikeCalculator) :: this
    class(TCalculationAtParamPoint) P
    real(mcp), intent(in), optional :: LogLike,weight
    integer, intent(in) :: aunit
    integer isused,i

    if (present(weight)) then
        write (aunit,*) ' weight    = ',weight
    end if

    if (present(LogLike)) then
        write (aunit,*) '-log(Like) = ',LogLike
        write (aunit,*) ' chi-sq    = ',LogLike*2
        write (aunit,*) ''
    end if

    do isused = 0,1
        do i=1, num_params
            if (isused==0 .and. BaseParams%varying(i) .or. isused==1 .and. .not. BaseParams%varying(i)) then
                write(aunit,'(1I5,1E15.7,"   ",1A22)', advance='NO') &
                    i, P%P(i), BaseParams%NameMapping%name(i)
                write (aunit,'(a)') trim(BaseParams%NameMapping%label(i))
            end if
        end do
        write (aunit,*) ''
    end do

    end subroutine TLikeCalculator_WriteParamsHumanText


    subroutine TLikeCalculator_WriteParamsHumanTextFile(this, fname, P, LogLike, weight)
    class(TLikeCalculator) :: this
    character(LEN=*), intent(in) :: fname
    class(TCalculationAtParamPoint) P
    real(mcp), intent(in), optional :: LogLike, weight
    Type(TTextFile) F

    call F%CreateFile(fname)
    call this%WriteParamsHumanText(F%unit, P, LogLike, weight)
    call F%Close()

    end subroutine TLikeCalculator_WriteParamsHumanTextFile


    subroutine TLikeCalculator_WriteParamPointTextData(this, output_root, Params)
    class(TLikeCalculator) :: this
    character(LEN=*), intent(in) :: output_root
    class(TCalculationAtParamPoint) Params
    end subroutine TLikeCalculator_WriteParamPointTextData


    function Generic_GetLogLikeMain(this, Params) result(LogLike)!Get -Ln(Likelihood) for chains
    use hod
    class(TGenericLikeCalculator) :: this
    class(TCalculationAtParamPoint) :: Params
    real(mcp) LogLike
    !!!!!!!!!hod!!!!!!!!!
    real(mcp) :: loglike_2pcf,loglike_3pcf
    real(mcp) :: m_cut,sigma,kappa,m1,alpha
    real(mcp) :: N_cent,N_sat
    real(mcp) :: rand
    integer :: count_cent,count_sate,count_total,tmp_nsate
    integer :: i,j
    m_cut=10.d0**Params%P(1)
    m1=10.d0**Params%P(2)
    sigma=Params%P(3)
    kappa=Params%P(4)
    alpha=Params%P(5)
    count_cent=0
    count_sate=0
    count_total=0
    allocate(temp(nhalos))
    do i=1,nhalos
!        write(*,*) i
!        if(mass(i).gt.m_cut) then
        N_cent = 1.d0/2.d0*erfc(dlog10(m_cut/mass(i))/(dsqrt(2.d0)*sigma))
        N_sat = ((mass(i)-kappa*m_cut)/m1)**alpha
        call RANDOM_NUMBER(rand)
        if(rand.lt.N_cent) then
            if(pos_h(i)%x.le.2500.d0.and.pos_h(i)%y.le.2500.d0.and.pos_h(i)%z.le.2500.d0) then
            count_total=count_total+1
            count_cent=count_cent+1
            temp(count_total)%x=pos_h(i)%x
            temp(count_total)%y=pos_h(i)%y
            temp(count_total)%z=pos_h(i)%z
            end if
        end if
!        write(*,*) N_sat,'breakpoint2'
        if(mass(i).le.kappa*m_cut) then
        tmp_nsate=0
        else
        tmp_nsate=poisson_inphase(N_sat)
        end if
!        write(*,*) tmp_nsate,'breakpoint3'
        if(tmp_nsate.gt.num_sat(i)) tmp_nsate=num_sat(i)
        do j=1,tmp_nsate
            count_sate=count_sate+1
            count_total=count_total+1
            temp(count_total)=pos_s(ind_sat(i)+j-1)
        end do
!        end if
    end do
    write(*,*) count_total, count_sate
    allocate(pos_g(count_total))
    pos_g(1:count_total)=temp(1:count_total)
    
    !!!!npcf calculation from here!!!!
    dummy=corr%setNumParticles(count_total)
        !!!!npcf calculation end here!!!!
    !!!!loglike calculation from here!!!!
    loglike_2pcf=0.d0
    loglike_3pcf=0.d0
    if(use_2pcf) then
    dummy=corr%calculate2pt(pos_g)
    dummy=corr%get2pt(twoPoint)
        do i=1,numShells
            do j=1,numShells
!                if(i==1.and.j==1) write(*,*) twopoint(i),twopcf_data(i)
                loglike_2pcf=loglike_2pcf+&
                (twopoint(i)-twopcf_data(i))*invcov2(i,j)*&
                (twopoint(j)-twopcf_data(j))    
!                write(*,*) i,j,loglike_2pcf
            end do
        end do
    end if
    if(use_3pcf) then
    dummy=corr%calculateCorrelations(pos_g)
    dummy=corr%get3pt(threePoint)
        do i=1,numTriangles
            do j=1,numTriangles
                loglike_3pcf=loglike_3pcf+&
                (threepoint(i)-threepcf_data(i))*invcov3(i,j)*&
                (threepoint(j)-threepcf_data(j))    
            end do
        end do
    end if
    deallocate(pos_g)
    deallocate(temp)


!    Loglike=loglike_2pcf+loglike_3pcf
    Loglike=loglike_2pcf
    write(*,*) Loglike   


    !!!!!!!!!hod end!!!!!
    !Used when you want to plug in your own CMB-independent likelihood function:
    !Parameter array is Params%P, so e.g. 2D unit Gaussian would be
    !LogLike = (Params%P(1)**2+Params%P(2)**2)/2
    !LogLike = LogZero
    !call MpiStop('Generic_GetLogLikeMain: need to write this function!')

    end function Generic_GetLogLikeMain


    subroutine TheoryLike_SetTheoryParams(this, Params)
    class(TTheoryLikeCalculator) :: this
    class(TCalculationAtParamPoint), target :: Params

    if (.not. allocated(this%TheoryParams)) call this%Config%Parameterization%NewTheoryParams(this%TheoryParams)
    call this%Config%Parameterization%ParamArrayToTheoryParams(Params%P,this%TheoryParams)
    this%Params => Params

    end subroutine TheoryLike_SetTheoryParams

    subroutine TheoryLike_SetNewTheoryResults(this, Params)
    class(TTheoryLikeCalculator) :: this
    class(TCalculationAtParamPoint) :: Params

    if (.not. allocated(Params%Theory)) call this%Config%NewTheory(Params%Theory)

    end subroutine TheoryLike_SetNewTheoryResults


    function TheoryLike_GetLogLikeMain(this, Params) result(LogLike)
    !Get -Ln(Likelihood), not accounting for temperature
    class(TTheoryLikeCalculator) :: this
    class(TCalculationAtParamPoint) :: Params
    real(mcp) LogLike

    call this%SetTheoryParams(Params)
    LogLike = this%Config%Parameterization%NonBaseParameterPriors(this%TheoryParams)
    if (LogLike == logZero) return
    if (.not. Params%validInfo) then
        this%changeMask(1:num_params) = .true.
    else
        this%changeMask(1:num_params) = Params%lastParamArray(1:num_params)/=Params%P(1:num_params)
    end if
    call this%SetNewTheoryResults(Params)
    if (this%CalculateRequiredTheoryChanges()) then
        call this%AddLike(LogLike, this%GetLogLikeWithTheorySet())
        Params%lastParamArray(1:num_params) = Params%P(1:num_params)
    else
        LogLike = logZero
    end if
    if (LogLike==logZero) return

    if (Feedback>2) call DataLikelihoods%WriteLikelihoodContribs(stdout, Params%likelihoods)

    end function TheoryLike_GetLogLikeMain

    function TheoryLike_CheckPriorCuts(this, Params) result(checkPriorCuts)
    class(TTheoryLikeCalculator) :: this
    real(mcp)  CheckPriorCuts
    class(TCalculationAtParamPoint) :: Params

    CheckPriorCuts = this%GetLogLikeBounds(Params)
    if (CheckPriorCuts==LogZero) return

    call this%SetTheoryParams(Params)
    CheckPriorCuts = this%Config%Parameterization%NonBaseParameterPriors(this%TheoryParams)

    end function TheoryLike_CheckPriorCuts


    function TheoryLike_GetLogLikePost(this,Params, do_like) result(LogLike)
    !for importance sampling where theory may be pre-stored
    class(TTheoryLikeCalculator) :: this
    real(mcp)  LogLike
    class(TCalculationAtParamPoint) :: Params
    logical, optional, intent(in) :: do_like(DataLikelihoods%count)

    LogLike = this%GetLogLikeBounds(Params)
    if (LogLike==LogZero) return

    this%SlowChanged = .false.
    this%ChangeMask = .true.

    call this%SetTheoryParams(Params)
    call this%AddLikeTemp(LogLike,this%Config%Parameterization%NonBaseParameterPriors(this%TheoryParams))
    if (LogLike == logZero) return
    call this%AddLikeTemp(LogLike, this%GetLogLikeWithTheorySet(do_like))
    if (LogLike == logZero) return
    call this%AddLikeTemp(LogLike, this%GetLogPriors(Params%P))

    end function TheoryLike_GetLogLikePost


    function TheoryLike_LogLikeWithTheorySet(this, likelihood_mask) result(logLike)
    class(TTheoryLikeCalculator) :: this
    logical, intent(in), optional :: likelihood_mask(DataLikelihoods%count)
    real(mcp) logLike
    real(mcp) itemLike
    class(TDataLikelihood), pointer :: like => null()
    integer i
    logical :: do_like(DataLikelihoods%count)
    Type(TTimer) Timer

    if (present(likelihood_Mask)) then
        do_like = likelihood_mask
    else
        do_like = .true.
    end if
    logLike = logZero
    call this%GetTheoryForLike(null()) !chance to initalize
    do i= 1, DataLikelihoods%count
        if (do_like(i)) then
            like => DataLikelihoods%Item(i)
            if (any(like%dependent_params(1:num_params) .and. this%changeMask(1:num_params) )) then
                call this%GetTheoryForLike(like)
                if (this%timing) call Timer%Start
                itemLike = like%GetLogLike(this%TheoryParams, this%Params%Theory, this%Params%P(like%nuisance_indices))
                if (this%timing) call Timer%WriteTime('Time for '//trim(like%name))
                if (itemLike == logZero) return
                this%Params%Likelihoods(i) = itemLike
            end if
        end if
    end do
    logLike = sum(this%Params%likelihoods(1:DataLikelihoods%Count))

    end function TheoryLike_LogLikeWithTheorySet


    logical function TheoryLike_CalculateRequiredTheoryChanges(this)
    class(TTheoryLikeCalculator) :: this

    !Set this%Params theory entries (this%Params%Theory) as required for current this%changeMask and this%TheoryParams
    this%Params%validInfo = .true.
    TheoryLike_CalculateRequiredTheoryChanges = .true.

    end function TheoryLike_CalculateRequiredTheoryChanges

    subroutine TheoryLike_GetTheoryForLike(this,Like)
    class(TTheoryLikeCalculator) :: this
    class(TDataLikelihood), pointer :: like

    !If needed, likelihood specific calculation/initalization; like=null for first initial call
    end subroutine TheoryLike_GetTheoryForLike

    subroutine TheoryLike_GetTheoryForImportance(this,Params, error)
    class(TTheoryLikeCalculator) :: this
    class(TCalculationAtParamPoint), target :: Params
    integer error

    call this%SetTheoryParams(Params)
    call this%Config%Calculator%GetTheoryForImportance(this%TheoryParams, Params%Theory, error)

    end subroutine TheoryLike_GetTheoryForImportance


    subroutine TheoryLike_UpdateTheoryForLikelihoods(this, Params)
    class(TTheoryLikeCalculator) :: this
    class(TCalculationAtParamPoint) :: Params

    end subroutine TheoryLike_UpdateTheoryForLikelihoods

    function TheoryLike_TestLikelihoodFunction(this,Params) result(LogLike)
    class(TTheoryLikeCalculator) :: this
    class(TCalculationAtParamPoint) Params
    real(mcp) :: LogLike

    call this%SetNewTheoryResults(Params) !so derived parameters can be output OK
    LogLike = this%TLikeCalculator%TestLikelihoodFunction(Params)

    end function TheoryLike_TestLikelihoodFunction


    subroutine TTheoryLike_WriteParamsHumanText(this, aunit, P, LogLike, weight)
    class(TTheoryLikeCalculator) :: this
    class(TCalculationAtParamPoint) P
    real(mcp), intent(in), optional :: LogLike, weight
    integer, intent(in) :: aunit
    real(mcp), allocatable :: derived(:)
    integer :: numderived = 0
    integer i

    call this%TLikeCalculator%WriteParamsHumanText(aunit, P, LogLike, weight)

    call this%Config%Parameterization%CalcDerivedParams(P%P,P%Theory, derived)
    call DataLikelihoods%addLikelihoodDerivedParams(P%P, P%Theory, derived, P%Likelihoods, LogLike)
    if (allocated(derived)) numderived = size(derived)
    do i=1, numderived
        write(aunit,'(1I5,1E15.7,"   ",1A22)', advance='NO') &
            num_params+i, derived(i), BaseParams%NameMapping%name(num_params + i )
        write (aunit,'(a)') trim(BaseParams%NameMapping%label(num_params+i))
    end do

    if (present(LogLike)) then
        write(aunit,*) ''
        write(aunit,*) '-log(Like)     chi-sq   data'
        call DataLikelihoods%WriteLikelihoodContribs(aunit, P%likelihoods)
    end if

    end subroutine TTheoryLike_WriteParamsHumanText

    subroutine TTheoryLike_WriteParamPointTextData(this, output_root, Params)
    class(TTheoryLikeCalculator) :: this
    character(LEN=*), intent(in) :: output_root
    class(TCalculationAtParamPoint) Params

    if (allocated(Params%Theory)) then
        call DataLikelihoods%WriteDataForLikelihoods(Params%P, Params%Theory, output_root)
        call Params%Theory%WriteTextData(output_root)
    end if

    end subroutine TTheoryLike_WriteParamPointTextData


    end module CalcLike
