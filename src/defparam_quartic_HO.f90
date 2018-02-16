MODULE quartic_HO
  !
  USE nrtype
  !
#if _OPENMP
  USE OMP_LIB
#endif
  !
  IMPLICIT NONE
  !
  ! Parameters
  INTEGER(KIND = I4B) :: ierr ! Error flag
  !
  INTEGER(KIND = I4B) :: HO_Dimension ! Truncated HO basis dimension v = 0, ... dim-1 
  !
  INTEGER(KIND = I4B), PARAMETER :: nHamoper = 3_I4B ! Number of Hamiltonian operators
  !
  REAL(KIND = DP) ::  s_value, A_value, B_value, C_value, D_value ! Scale and Hamiltonian parameters 
  ! 
  LOGICAL :: eigenvec ! If .T. compute eigenvalues and eigenvectors
  LOGICAL :: excitation ! If .T. compute excitation energy with respect to L = 0 g.s.
  !
  ! Arrays
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: Ham_quartic_mat ! Hamiltonian matrix
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: eigenval_vec ! Hamiltonian Eigenvalues
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: emean_vec ! Mean energy value of each two eigenstates
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: omegaeff_vec ! Energy increment for adjancent levels
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: ipr_vec ! Inverse Participation Ratio
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: Husimi_ipr_vec ! Husimi inverse Participation Ratio
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: expected_n_vec ! n expectation value
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: expected_x_vec ! x expectation value
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: expected_x2_vec ! x^2 expectation value
  !
  ! Benchmarking variables
  LOGICAL :: benchmark ! If .T. carry out benchmark calculations
  INTEGER(KIND = I4B) :: benchmark_total_iter = 100, benchmark_iter
  REAL(KIND = SP) :: time_check, time_check_ref
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: benchmark_MFLOPS
  REAL(KIND = DP) :: mean_Mflops, std_Mflops
  !
  !Convergence variables
  LOGICAL :: convergence ! If .T. carry out convergence checks
  INTEGER(KIND = I4B) :: dim_step ! HO_dimension step
  REAL(KIND = DP) :: delta_energy ! Energy error tolerance
  !
  !
  REAL(KIND = DP) :: GS_energy ! Ground state energy 
  !
  !
  REAL(KIND = DP), PARAMETER :: Zero_Parameter = 1.0E-20_DP ! Zero limit to discard a parameter evaluation
  !
  !
  ! Control verbosity
  INTEGER(KIND = I4B) :: Iprint = 1_I4B
  !
  ! Flag to control storing results
  LOGICAL :: save_output
  CHARACTER(LEN=75) :: file_name
  !
  ! OpenMP Definitions
  INTEGER(KIND = I4B) :: N_threads = 1
  !
CONTAINS
  !
#ifndef  __GFORTRAN__
  SUBROUTINE Allocate_arrays_sub(HO_Dimension, output_states, eigenvec, save_output)
#else
  SUBROUTINE Allocate_arrays_sub(HO_Dimension, output_states, save_output)    
#endif
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: HO_Dimension ! Dimension of the truncated HO basis
    INTEGER(KIND = I4B), INTENT(IN) :: output_states ! States saved and considered after diagonalization
    !
#ifndef  __GFORTRAN__
    LOGICAL, INTENT(IN) :: eigenvec
#endif
    LOGICAL, INTENT(IN) :: save_output
    !
    ! ALLOCATE HAMILTONIAN
     ALLOCATE(Ham_quartic_mat(1:HO_Dimension,1:HO_Dimension), STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "Ham_quartic_mat allocation request denied."
        STOP
     ENDIF
     !
     Ham_quartic_mat = 0.0_DP
     !
     !
     ! ALLOCATE EIGENVALUES VECTOR
     ALLOCATE(eigenval_vec(1:HO_Dimension), STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "eigenval_vec allocation request denied."
        STOP
     ENDIF
     !
     eigenval_vec = 0.0_DP
     !
#ifndef  __GFORTRAN__
     IF (eigenvec) THEN
        ! ifort
        ALLOCATE(Z(1:HO_Dimension, 1:HO_Dimension), STAT = IERR)    
        IF (IERR /= 0) THEN
           WRITE(UNIT = *, FMT = *) "Z allocation request denied."
           STOP
        ENDIF
        !
        Z = 0.0_DP     
     ENDIF
#endif
     !
     ! mean energy value and energy steps
     IF (save_output) THEN
        ALLOCATE(emean_vec(1:output_states-1), STAT = IERR)    
        IF (IERR /= 0) THEN
           WRITE(UNIT = *, FMT = *) "emean_vec allocation request denied."
           STOP
        ENDIF
        emean_vec = 0.0_DP     
        !
        ALLOCATE(omegaeff_vec(1:output_states-1), STAT = IERR)    
        IF (IERR /= 0) THEN
           WRITE(UNIT = *, FMT = *) "omegaeff_vec allocation request denied."
           STOP
        ENDIF
        omegaeff_vec = 0.0_DP     
        !
        ALLOCATE(ipr_vec(1:output_states), STAT = IERR)    
        IF (IERR /= 0) THEN
           WRITE(UNIT = *, FMT = *) "ipr_vec allocation request denied."
           STOP
        ENDIF
        ipr_vec = 0.0_DP
        !
        ALLOCATE(Husimi_ipr_vec(1:output_states), STAT = IERR)    
        IF (IERR /= 0) THEN
           WRITE(UNIT = *, FMT = *) "Husimi_ipr_vec allocation request denied."
           STOP
        ENDIF
        Husimi_ipr_vec = 0.0_DP
        !
        ALLOCATE(expected_n_vec(1:output_states), STAT = IERR)    
        IF (IERR /= 0) THEN
           WRITE(UNIT = *, FMT = *) "expected_n_vec allocation request denied."
           STOP
        ENDIF
        expected_n_vec = 0.0_DP     
        !
        ALLOCATE(expected_x_vec(1:output_states), STAT = IERR)    
        IF (IERR /= 0) THEN
           WRITE(UNIT = *, FMT = *) "expected_x_vec allocation request denied."
           STOP
        ENDIF
        expected_x_vec = 0.0_DP     
        !
        ALLOCATE(expected_x2_vec(1:output_states), STAT = IERR)    
        IF (IERR /= 0) THEN
           WRITE(UNIT = *, FMT = *) "expected_x2_vec allocation request denied."
           STOP
        ENDIF
        expected_x2_vec = 0.0_DP     
        !
     ENDIF
     !
   END SUBROUTINE Allocate_arrays_sub
   !
#ifndef  __GFORTRAN__
   SUBROUTINE Deallocate_arrays_sub(eigenvec, save_output)
#else
   SUBROUTINE Deallocate_arrays_sub(save_output)
#endif
    !
    IMPLICIT NONE
    !
#ifndef  __GFORTRAN__
    LOGICAL, INTENT(IN) :: eigenvec
#endif
    LOGICAL, INTENT(IN) :: save_output
    !
    ! DEALLOCATE HAMILTONIAN
     DEALLOCATE(Ham_quartic_mat, STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "Ham_quartic_mat deallocation request denied."
        STOP
     ENDIF
     !
     ! DEALLOCATE EIGENVALUES VECTOR
     DEALLOCATE(eigenval_vec, STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "eigenval_vec deallocation request denied."
        STOP
     ENDIF
     !
#ifndef  __GFORTRAN__
     IF (eigenvec) THEN
        ! ifort
        DEALLOCATE(Z, STAT = IERR)    
        IF (IERR /= 0) THEN
           WRITE(UNIT = *, FMT = *) "Z deallocation request denied."
           STOP
        ENDIF
        !
     ENDIF
#endif
     !
     ! mean energy value and energy steps
     IF (save_output) THEN
        DEALLOCATE(emean_vec, STAT = IERR)    
        IF (IERR /= 0) THEN
           WRITE(UNIT = *, FMT = *) "emean_vec deallocation request denied."
           STOP
        ENDIF
        !
        DEALLOCATE(omegaeff_vec, STAT = IERR)    
        IF (IERR /= 0) THEN
           WRITE(UNIT = *, FMT = *) "omegaeff_vec deallocation request denied."
           STOP
        ENDIF
        !
        DEALLOCATE(ipr_vec, STAT = IERR)    
        IF (IERR /= 0) THEN
           WRITE(UNIT = *, FMT = *) "ipr_vec deallocation request denied."
           STOP
        ENDIF
        !
        DEALLOCATE(Husimi_ipr_vec, STAT = IERR)    
        IF (IERR /= 0) THEN
           WRITE(UNIT = *, FMT = *) "Husimi_ipr_vec deallocation request denied."
           STOP
        ENDIF
        !
        DEALLOCATE(expected_n_vec, STAT = IERR)    
        IF (IERR /= 0) THEN
           WRITE(UNIT = *, FMT = *) "expected_n_vec deallocation request denied."
           STOP
        ENDIF
        !
        DEALLOCATE(expected_x_vec, STAT = IERR)    
        IF (IERR /= 0) THEN
           WRITE(UNIT = *, FMT = *) "expected_x_vec deallocation request denied."
           STOP
        ENDIF
        !
        DEALLOCATE(expected_x2_vec, STAT = IERR)    
        IF (IERR /= 0) THEN
           WRITE(UNIT = *, FMT = *) "expected_x2_vec deallocation request denied."
           STOP
        ENDIF
        !
     ENDIF
     !
   END SUBROUTINE Deallocate_arrays_sub
   !
   !
   SUBROUTINE QUARTIC_HAMILTONIAN_HO(N_val, Ham_quartic_mat, alpha_value, beta_value, gamma_value, delta_value) 
    !
    ! Subroutine to build the Hamiltonian for a Quartic potential
    ! Potential:
    ! V_C(x) = A\, x^4 + B\, x^2 + C\,x + D\,x^3~,
    !
    ! n = 1 case
    !\begin{align}
    ! \hat{\cal H}^{(n=1)} =& \epsilon_0 + \left( 3\alpha + \beta + \frac{1}{2}\right)\hat n
    ! + \frac{\alpha}{4}({a^\dagger}^4 + a^4) +  \frac{3\alpha}{2} {a^\dagger}^2 a^2 + \alpha
    ! ({a^\dagger} a^3 +{a^\dagger}^3 a)  \label{qhamn1}\\
    ! &+ \frac{6\alpha +2\beta -1}{4}({a^\dagger}^2 +  a^2)+ \frac{\gamma}{\sqrt{2}} (a^\dagger + a) ~~,\nonumber
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! Dimension of the truncated HO
    !
    REAL(KIND = DP), DIMENSION(:,:), INTENT(OUT) :: Ham_quartic_mat ! Hamiltonian matrix
    !
    REAL(KIND = DP) ::  alpha_value, beta_value, gamma_value, delta_value ! Scaled Hamiltonian parameters
    !
    !
    INTEGER(KIND = I4B) :: state_index
    REAL(KIND = DP) :: nval
    !
    !
    !
    ! Build Hamiltonian
    !
    !
    !
    !
    !$OMP PARALLEL SECTIONS DEFAULT(NONE) PRIVATE(state_index, nval) &
    !$OMP & SHARED(N_val, alpha_value, beta_value, gamma_value, delta_value, Ham_quartic_mat, Iprint)
    !
    !$OMP SECTION
#if _OPENMP
    IF (Iprint > 1) WRITE(*,*) "Operators number 1, 3 :: a^d a,  and zero energy, Thread Number ", &
         OMP_GET_THREAD_NUM() 
#else
    IF (Iprint > 1) WRITE(*,*) "Operators number 1, 3 :: a^d a, (a^d a)**2, and zero energy"
#endif
    !
    DO state_index = 1, N_val
       !
       nval = REAL(state_index-1,DP)
       !
       Ham_quartic_mat(state_index, state_index) = Ham_quartic_mat(state_index, state_index) + &
            !
            (3.0_DP*alpha_value + beta_value + 0.5_DP)*nval + &  ! number operator n
            !
            1.5_DP*alpha_value*(nval-1.0_DP)*nval + &          ! (a^d a)**2 ( state_index = 1, 2 terms are zero)
            !
            0.75_DP*alpha_value + 0.5_DP*beta_value + 0.25_DP    ! Zero energy
       !
    ENDDO
    !
    !
    !$OMP SECTION
#if _OPENMP
    IF (Iprint > 1) WRITE(*,*) "Operator number 2 :: (a^d)**4, Thread Number ", OMP_GET_THREAD_NUM() 
#else
    IF (Iprint > 1) WRITE(*,*) "Operator number 2 :: (a^d)**4"
#endif
    !
    DO state_index = 1, N_val - 4
       !
       nval = REAL(state_index-1,DP)
       Ham_quartic_mat(state_index, state_index + 4) = Ham_quartic_mat(state_index, state_index + 4 ) + &
            0.25_DP*alpha_value*SQRT((nval+1.0_DP)*(nval+2.0_DP)*(nval+3.0_DP)*(nval+4.0_DP))
       !
    ENDDO
    !
    !
    !$OMP SECTION
#if _OPENMP
    IF (Iprint > 1) WRITE(*,*) "Operator number 4 and 5 :: (a^d)**3 a and (a^d)**2, Thread Number ", &
         OMP_GET_THREAD_NUM() 
#else
    IF (Iprint > 1)  WRITE(*,*) "Operators number 4 and 5 :: (a^d)**3 a and (a^d)**2"
#endif
    !
    DO state_index = 1, N_val - 2
       !
       nval = REAL(state_index-1,DP)
       Ham_quartic_mat(state_index, state_index + 2) = Ham_quartic_mat(state_index, state_index + 2) + &
            !
            alpha_value*nval*SQRT((nval+1.0_DP)*(nval + 2.0_DP)) + & ! (a^d)**3 a (state_index = 1 is zero)
            !
            (1.5_DP*alpha_value + 0.5_DP*beta_value - 0.25_DP)*SQRT((nval+1.0_DP)*(nval + 2.0_DP)) ! (a^d)**2
       !
    ENDDO
    !
    !
    !$OMP SECTION
#if _OPENMP
    IF (Iprint > 1) WRITE(*,*) "Operators number 6 and 8 :: a^d and (a^d)**2a, Thread Number ", &
         OMP_GET_THREAD_NUM() 
#else
    IF (Iprint > 1)  WRITE(*,*) "Operators number 6 and 8 :: a^d and (a^d)**2a"
#endif
    !
    DO state_index = 1, N_val - 1
       !
       nval = REAL(state_index-1,DP)
       Ham_quartic_mat(state_index, state_index + 1) = Ham_quartic_mat(state_index, state_index + 1) + &
            ((gamma_value*SQRT(0.5_DP)+3.0_DP*delta_value*SQRT(0.125_DP)) + &
            3.0_DP*delta_value*nval*SQRT(0.125_DP))*SQRT(nval+1.0_DP)
       !
    ENDDO
    !
    !$OMP SECTION
#if _OPENMP
    IF (Iprint > 1) WRITE(*,*) "Operator number 7 :: (a^d)**3, Thread Number ", OMP_GET_THREAD_NUM() 
#else
    IF (Iprint > 1)  WRITE(*,*) "Operator number 7 :: (a^d)**3"
#endif
    !
    DO state_index = 1, N_val - 3
       !
       nval = REAL(state_index-1,DP)
       Ham_quartic_mat(state_index, state_index + 3) = Ham_quartic_mat(state_index, state_index + 3) + &
            delta_value*SQRT(0.125_DP)*SQRT((nval+1.0_DP)*(nval+2.0_DP)*(nval+3.0_DP))
       !
    ENDDO
    !
    !$OMP END PARALLEL SECTIONS
    !
  END SUBROUTINE QUARTIC_HAMILTONIAN_HO
  !
  SUBROUTINE Convergence_check(dim, dim_step, delta_energy, num_states, &
       alpha_value, beta_value, gamma_value, delta_value)
    !
    ! Check enegy convergence. Criterion: |E_calc(dim1) - E_calc(dim1+delta_energy)| < delta_energy for the first num_states states
    !
    ! Arguments
    !
    ! dim  INOUT :: initial and final dimension of the HO basis
    ! dim_step     IN :: increments in basis dimension for iteration
    ! delta_energy IN :: Maximum allowed energy difference to state convergence attained
    ! num_states   IN :: Number of eigenstates to chek
    ! alpha_value, beta_value, gamma_value IN :: Hamiltonian parameters
    ! Ham_quartic_mat  OUT :: System converged eigenstates
    ! Ham_quartic_aval OUT :: System converged eigenvalues
    !
    !
    ! Lapack 95
#ifdef  __GFORTRAN__
    ! gfortran
    USE LA_PRECISION, ONLY: WP => DP
    USE F95_LAPACK, ONLY: LA_SYEVR
#else
    !ifort
    USE F95_PRECISION, ONLY: WP => DP
    USE LAPACK95, ONLY: SYEVR
#endif
    !
    !
    IMPLICIT NONE
    !
    INTEGER(KIND=I4B), INTENT(INOUT) :: dim
    INTEGER(KIND=I4B), INTENT(IN) :: dim_step, num_states
    REAL(KIND = DP), INTENT(IN) :: delta_energy
    REAL(KIND = DP), INTENT(IN) :: alpha_value, beta_value, gamma_value, delta_value
    !
    ! Local Variables
    REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: Ham_mat_sub
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: Ham_aval_sub
    !
    LOGICAL :: conv_search 
    !
    REAL(KIND = DP), DIMENSION(1:num_states) :: Ham_aval_0
    INTEGER(KIND=I4B)  ::  dim_val
    !
    ! Initial diagonalization
    conv_search = .TRUE.
    ! 
    ! ALLOCATE HAMILTONIAN
    ALLOCATE(Ham_mat_sub(1:dim,1:dim), STAT = IERR)    
    IF (IERR /= 0) THEN
       WRITE(UNIT = *, FMT = *) "Ham_mat_sub allocation request denied."
       STOP
    ENDIF
    !
    Ham_mat_sub = 0.0_DP
    !
    !
    ! ALLOCATE EIGENVALUES VECTOR
    ALLOCATE(Ham_aval_sub(1:dim), STAT = IERR)    
    IF (IERR /= 0) THEN
       WRITE(UNIT = *, FMT = *) "Ham_aval_sub allocation request denied."
       STOP
    ENDIF
    !
    Ham_aval_sub = 0.0_DP
    !
    ! BUILD HAMILTONIAN
    CALL QUARTIC_HAMILTONIAN_HO(dim, Ham_mat_sub, alpha_value, beta_value, gamma_value, delta_value)
    !      
    ! Diagonalize Hamiltonian matrix (LAPACK95)
#ifdef  __GFORTRAN__
    ! gfortran
    CALL LA_SYEVR(A=Ham_mat_sub, W=Ham_aval_sub, JOBZ='N', UPLO='U')
#else
    !ifort
    CALL SYEVR(A=Ham_mat_sub, W=Ham_aval_sub, UPLO='U')
#endif           
    !
    Ham_aval_0 = Ham_aval_sub(1:num_states)
    !
    IF (Iprint > 0) THEN
       WRITE(UNIT = *, FMT = *) "Convergence search started."
       WRITE(UNIT = *, FMT = *) "Initial dimension = ", dim, "; num_states = ", num_states
       WRITE(UNIT = *, FMT = *) "Initial eigenvalues: ", Ham_aval_0
    ENDIF
    !
    ! DEALLOCATE HAMILTONIAN
    DEALLOCATE(Ham_mat_sub, STAT = IERR)    
    IF (IERR /= 0) THEN
       WRITE(UNIT = *, FMT = *) "Ham_mat_sub deallocation request denied."
       STOP
    ENDIF
    !
    ! DEALLOCATE EIGENVALUES VECTOR
    DEALLOCATE(Ham_aval_sub, STAT = IERR)    
    IF (IERR /= 0) THEN
       WRITE(UNIT = *, FMT = *) "Ham_aval_sub deallocation request denied."
       STOP
    ENDIF
    !
    ! New dim value
    dim_val = dim + dim_step
    !
    ! ALLOCATE HAMILTONIAN
    ALLOCATE(Ham_mat_sub(1:dim_val,1:dim_val), STAT = IERR)    
    IF (IERR /= 0) THEN
       WRITE(UNIT = *, FMT = *) "Ham_mat_sub allocation request denied."
       STOP
    ENDIF
    !
    Ham_mat_sub = 0.0_DP
    !
    !
    ! ALLOCATE EIGENVALUES VECTOR
    ALLOCATE(Ham_aval_sub(1:dim_val), STAT = IERR)    
    IF (IERR /= 0) THEN
       WRITE(UNIT = *, FMT = *) "Ham_aval_sub allocation request denied."
       STOP
    ENDIF
    !
    Ham_aval_sub = 0.0_DP
    !
    DO WHILE (conv_search)
       !
       IF (Iprint > 0) THEN
          !
          WRITE(UNIT = *, FMT = *) "Building Hamiltonian for dimension = ", dim_val
       ENDIF
       !
       ! BUILD HAMILTONIAN
       CALL QUARTIC_HAMILTONIAN_HO(dim_val, Ham_mat_sub, alpha_value, beta_value, gamma_value, delta_value)
       ! Diagonalize Hamiltonian matrix (LAPACK95)
#ifdef  __GFORTRAN__
       ! gfortran
       CALL LA_SYEVR(A=Ham_mat_sub, W=Ham_aval_sub, JOBZ='N', UPLO='U')
#else
       !ifort
       CALL SYEVR(A=Ham_mat_sub, W=Ham_aval_sub, UPLO='U')
#endif           
       !
       IF (Iprint > 1) THEN
          WRITE(UNIT = *, FMT = *) "Dimension = ", dim_val
          WRITE(UNIT = *, FMT = *) ABS(Ham_aval_0 - Ham_aval_sub(1:num_states))
       ENDIF
       !
       IF (ANY(ABS(Ham_aval_0 - Ham_aval_sub(1:num_states)) >  delta_energy)) THEN
          ! one step more
          Ham_aval_0 = Ham_aval_sub(1:num_states)
          !
          ! DEALLOCATE HAMILTONIAN
          DEALLOCATE(Ham_mat_sub, STAT = IERR)    
          IF (IERR /= 0) THEN
             WRITE(UNIT = *, FMT = *) "Ham_mat_sub deallocation request denied."
             STOP
          ENDIF
          !
          ! DEALLOCATE EIGENVALUES VECTOR
          DEALLOCATE(Ham_aval_sub, STAT = IERR)    
          IF (IERR /= 0) THEN
             WRITE(UNIT = *, FMT = *) "Ham_aval_sub deallocation request denied."
             STOP
          ENDIF
          !
          ! New dim value
          dim_val = dim_val + dim_step
          !
          ! ALLOCATE HAMILTONIAN
          ALLOCATE(Ham_mat_sub(1:dim_val,1:dim_val), STAT = IERR)    
          IF (IERR /= 0) THEN
             WRITE(UNIT = *, FMT = *) "Ham_mat_sub allocation request denied."
             STOP
          ENDIF
          !
          Ham_mat_sub = 0.0_DP
          !
          !
          ! ALLOCATE EIGENVALUES VECTOR
          ALLOCATE(Ham_aval_sub(1:dim_val), STAT = IERR)    
          IF (IERR /= 0) THEN
             WRITE(UNIT = *, FMT = *) "Ham_aval_sub allocation request denied."
             STOP
          ENDIF
          !
          Ham_aval_sub = 0.0_DP
       ELSE
          ! convergence reached
          dim = dim_val
          conv_search = .FALSE.
          IF (Iprint > 0) WRITE(UNIT = *, FMT = *) "Convergence attained for dimension = ", dim
       ENDIF
       !
    END DO
    !  
  END SUBROUTINE Convergence_check
  ! 
  SUBROUTINE Calculations_Sc(dim, n_states, eigenval_vec, eigenvec_arr, &
       emean_vec, omegaeff_vec, ipr_vec, Husimi_ipr_vec, expected_n_vec, &
       expected_x_vec, expected_x2_vec)
    ! 
    ! Function to compute the mean value of eigenvalues and the omega_eff as
    !
    ! Eigenvalues E_1, E_2, ... , E_N
    !
    ! emean_vec_i = (E_i+1 + E_i)/2
    ! omegaeff_vec_i = (E_i+1 - E_i)/\Delta_n with \Delta_n = 1
    !
    ! dim is the total basis dimension and n_states is the number of states under study
    !
    IMPLICIT NONE
    !
    INTEGER(KIND=I4B), INTENT(IN) :: dim, n_states
    !
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: eigenval_vec
    REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: eigenvec_arr
    REAL(KIND = DP), DIMENSION(:), INTENT(OUT) :: emean_vec, omegaeff_vec, ipr_vec, Husimi_ipr_vec, &
         expected_n_vec, expected_x_vec, expected_x2_vec
    !
    ! Local variables
    INTEGER(KIND = I4B) state_index
    !
    ! emean omega_eff
    DO state_index = 1, n_states-1
       emean_vec(state_index) = (eigenval_vec(state_index) + eigenval_vec(state_index+1))/2.0_DP
       omegaeff_vec(state_index) = eigenval_vec(state_index+1) - eigenval_vec(state_index)
    ENDDO
    !
    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(state_index) &
    !$OMP & SHARED(dim, n_states, eigenvec_arr, ipr_vec, Husimi_ipr_vec, expected_n_vec, expected_x_vec, expected_x2_vec, Iprint)
    !
    ! ipr, Husimi, <n>, <x>, and <x2>
    DO state_index = 1, n_states
#if _OPENMP
       IF (Iprint > 1) WRITE(*,*) "ipr, Husimi, <n>, <x>, and <x2>, state ", state_index, &
            " Thread Number  ", OMP_GET_THREAD_NUM() 
#else
       IF (Iprint > 1) WRITE(*,*) "ipr, Husimi, <n>, <x>, and <x2>, state ", state_index
#endif
       !
       ipr_vec(state_index) = Inv_Part_Ratio(eigenvec_arr(:,state_index))
       !       Husimi_ipr_vec(state_index) = IPR_Husimi(dim, eigenvec_arr(:,state_index))
       Husimi_ipr_vec(state_index) = 0.0_DP
       expected_n_vec(state_index) = Exp_Value_n(dim, eigenvec_arr(:,state_index))
       expected_x_vec(state_index) = Exp_Value_x(dim, eigenvec_arr(:,state_index))
       expected_x2_vec(state_index) = Exp_Value_x2(dim, eigenvec_arr(:,state_index))
    ENDDO
    !
    !$OMP END PARALLEL DO 
    !
  END SUBROUTINE Calculations_Sc
  !
  !
  FUNCTION IPR_Husimi(dim, eigenvec)
    ! 
    ! Function to compute the IPR for cusp eigenvectors in a HO basis
    !
    IMPLICIT NONE
    !
    INTEGER(KIND=I4B), INTENT(IN) :: dim
    REAL(KIND = DP), DIMENSION(:) :: eigenvec
    !
    REAL(KIND=DP) :: IPR_Husimi
    !
    ! Local variables
    INTEGER(KIND = I4B) n_1, n_2, n_3, n_4
    !
    IPR_Husimi = 0.0_DP
    !
    ! TODO Can be optimized n4 = n1+n3-n2, remove a cycle and split the operations in the cycle
    DO n_1 = 0, dim-1
       DO n_2 = 0, dim-1
          DO n_3 = 0, dim-1
             DO n_4 = 0, dim-1
                !
                IF (n_1+n_3 /= n_2+n_4) CYCLE
                !
                IPR_Husimi = IPR_Husimi + &
                     eigenvec(n_1+1)*eigenvec(n_2+1)*eigenvec(n_3+1)*eigenvec(n_4+1) * &
                     gamma_quot_ipr((/n_1,n_2,n_3,n_4/))
                !
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !
  END FUNCTION IPR_HUSIMI
  !
  !
  FUNCTION gamma_quot_ipr(nvals)
    !
    ! 2^-(n_1+n_3+1) * Gamma(n_1 + n_3 + 1) / \prod_{k = 1}^4 SQRT(Gamma(n_k +1))
    !
    ! TODO can be optimized calculating and saving factor values. Change arguments structure (n1+n3, 
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), DIMENSION(1:4), INTENT(IN) :: nvals
    !
    REAL(KIND = DP) :: gamma_quot_ipr
    !
    ! Local variables
    !
    INTEGER(KIND = I4B) :: index
    REAL(KIND = DP), DIMENSION(1:4) :: z
    REAL(KIND = DP) :: Term_1!, Term_2
    !
    z = REAL(nvals, DP)
    !
    Term_1 = -(z(1) + z(3) + 1.0_DP)*LOG(2.0_DP)
    !
    Term_1 = Term_1 + loggamma(z(1) + z(3) + 1.0_DP)
    !
    !
    !Term_2 = 0.0_DP
    !
    DO index = 1, 4
       !
       Term_1 = Term_1 - loggamma(z(index) + 1.0_DP)/2.0_DP
       !
    ENDDO
    !
    gamma_quot_ipr = EXP(Term_1)
    !
  END FUNCTION gamma_quot_ipr
  !
  !
  FUNCTION loggamma(z)
    !
    ! Value of ln[Gamma(z)]. For z > gamma_limit compute
    ! the asymptotic expansion of Ln(Gamma(z)) Abram. 6.1.41
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), intent(in) :: z
    !
    REAL(KIND = DP) :: loggamma
    !
    REAL(KIND = DP), PARAMETER :: gamma_limit = 30.0_DP
    !
    IF (z < gamma_limit) THEN
       !
       loggamma = LOG(GAMMA(z))
    ELSE
       !
       loggamma = -z + LOG(z)*(z-0.5_DP) + 0.5_DP*LOG(2.0_DP*PI_D)  &
            + 1.0_DP/(12.0_DP*z) &
            - 1.0_DP/(360.0_DP*z**3) 
       !
    ENDIF
  END FUNCTION loggamma
  !
  FUNCTION Inv_Part_Ratio(eigenvector)
    !
    ! Subroutine to compute the IPR for a Quantum Cusp eigenstate
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: eigenvector ! Eigenvector components (HO basis)
    !
    REAL(KIND = DP) :: Inv_Part_Ratio
    !
    !
    Inv_Part_Ratio = 1.0_DP/SUM(eigenvector**4)
    !
  END FUNCTION Inv_Part_Ratio
  !
  FUNCTION Exp_Value_X(dim_block, eigenvector)
    !
    ! Subroutine to compute the < X > for a Quantum Cusp eigenstate
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_block
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: eigenvector ! Eigenvector components (HO basis)
    !
    REAL(KIND = DP) :: Exp_Value_X
    !
    ! Local variables
    INTEGER(KIND = I4B) :: basis_index
    !
    Exp_Value_X = 0.0_DP
    !
    DO basis_index = 2, dim_block ! This is the Fortran index. n = basis_index - 1
       Exp_Value_X = Exp_Value_X + eigenvector(basis_index)*eigenvector(basis_index-1)*SQRT(REAL(basis_index - 1,DP))
       !
    ENDDO
    !
    !
    Exp_Value_X = Exp_Value_X*SQRT(2.0_DP)
    !
  END FUNCTION Exp_Value_X
  !
  FUNCTION Exp_Value_X2(dim_block, eigenvector)
    !
    ! Subroutine to compute the < X^2 > for a Quantum Cusp eigenstate
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_block
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: eigenvector ! Eigenvector components (HO basis)
    !
    REAL(KIND = DP) :: Exp_Value_X2
    !
    ! Local variables
    INTEGER(KIND = I4B) :: basis_index
    !
    Exp_Value_X2 = 0.0_DP
    !
    ! Diagonal part ( 1 + 2 a^\dagger a)/2
    DO basis_index = 1, dim_block ! This is the Fortran index. n = basis_index - 1
       !
       Exp_Value_X2 = Exp_Value_X2 + (REAL(basis_index-1,DP) + 0.5_DP)*eigenvector(basis_index)**2
       !
    ENDDO
    !
    ! Non-Diagonal part (a^\dagger^2 + a^2)/2
    DO basis_index = 3, dim_block ! This is the Fortran index. n = basis_index - 1
       !
       Exp_Value_X2 = Exp_Value_X2 + eigenvector(basis_index)*eigenvector(basis_index-2)*&
            SQRT(REAL(basis_index - 1,DP)*REAL(basis_index - 2,DP))
       !
    ENDDO
    !
  END FUNCTION Exp_Value_X2
  !
  FUNCTION Exp_Value_n(dim_block, eigenvector)
    !
    ! Subroutine to compute the < n = a^\dagger a > for a Quantum Cusp eigenstate
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_block
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: eigenvector ! Eigenvector components (HO basis)
    !
    REAL(KIND = DP) :: Exp_Value_n
    !
    ! Local variables
    INTEGER(KIND = I4B) :: basis_index
    !
    Exp_Value_n = 0.0_DP
    !
    DO basis_index = 1, dim_block ! This is the Fortran index. n = basis_index - 1
       !
       Exp_Value_n = Exp_Value_n + REAL(basis_index - 1,DP)*eigenvector(basis_index)**2
       !
    ENDDO
    !
    !
  END FUNCTION Exp_Value_n
  !
  !
  FUNCTION output_filename(dim_value, filename_basis)
    !
    ! Function to prepare the a filename string according with the rule
    !
    ! output_filename = filename_basis_<dim_value>.dat
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_value ! Hamiltonian block dimension
    !
    CHARACTER(LEN=*), INTENT(IN) :: filename_basis
    !
    CHARACTER(LEN=75) :: output_filename
    !
    !
    ! Build filename
    IF ( dim_value < 10) THEN !to avoid spaces
       WRITE(output_filename, '(A, I1, ".dat")') TRIM(filename_basis), dim_value
    ELSE IF ( dim_value < 100) THEN 
       WRITE(output_filename, '(A, I2, ".dat")') TRIM(filename_basis), dim_value
    ELSE IF ( dim_value < 1000) THEN 
       WRITE(output_filename, '(A, I3, ".dat")') TRIM(filename_basis), dim_value
    ELSE IF ( dim_value < 10000) THEN 
       WRITE(output_filename, '(A, I4, ".dat")') TRIM(filename_basis), dim_value
    ELSE IF ( dim_value < 100000) THEN 
       WRITE(output_filename, '(A, I5, ".dat")') TRIM(filename_basis), dim_value
    ELSE
       WRITE(output_filename, '(A, I6, ".dat")') TRIM(filename_basis), dim_value
    ENDIF
    !
  END FUNCTION output_filename
  !
  SUBROUTINE Save_Output_Vector(filename, X_vector, Y_vector, dim_X)
    !
    ! Arguments (All INTENT(IN))
    !
    ! filename            :: Character, name of the file. By default the file is opened with STATUS = "UNKNOWN"
    ! X_vector (OPTIONAL) :: Abscissa data
    ! Y_vector            :: Ordinate (Vector) data
    ! dim_X               :: Vector dimension
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: filename
    !
    REAL(KIND = DP), DIMENSION(:), OPTIONAL, INTENT(IN) :: X_vector
    !
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: Y_vector
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_X
    !
    ! Local variables
    INTEGER(KIND = I4B) :: X_index
    !
    ! Open file
    OPEN (UNIT = 1000, FILE = filename, STATUS = "UNKNOWN")
    !
    ! Save results
    IF (PRESENT(X_vector)) THEN
       !
       DO X_index = 1, dim_X
          !
          WRITE(UNIT = 1000, FMT = *) X_vector(X_index), Y_vector(X_index)
          !
       ENDDO
       !
    ELSE
       !
       DO X_index = 1, dim_X
          !
          WRITE(UNIT = 1000, FMT = *) X_index, Y_vector(X_index)
          !
       ENDDO
       !
    ENDIF
    !
    CLOSE(UNIT = 1000)
    !
  END SUBROUTINE Save_Output_Vector
  !
  SUBROUTINE Save_Output_Array(filename, Y_array, dim_X, dim_Y, X_vector, column_major_mode)
    !
    ! Arguments (All INTENT(IN))
    !
    ! filename            :: Character, name of the file. By default the file is opened with STATUS = "UNKNOWN"
    ! Y_array             :: Array data (dim_X x dim_Y)
    ! dim_X               :: Array number of rows
    ! dim_Y    (OPTIONAL) :: Array number of columns. Default value dim_Y = dim_X.
    ! X_vector (OPTIONAL) :: Abscissa data
    ! column_major_mode (OPTIONAL) :: Logical
    !
    !    If column_major_mode = .TRUE. (DEFAULT)
    !
    !    Output file:
    !
    !    [X_1       X_2       ....   X_dim_Y]
    !     Y_11      Y_12      ....   Y_1dim_Y
    !     Y_21      Y_22      ....   Y_2dim_Y
    !     .         .                .
    !     .         .                .
    !     .         .                .
    !     Y_dim_X1  Y_dim_X2  ....   Y_dim_Xdim_Y
    !
    !    If column_major_mode = .FALSE.
    !
    !    Output file:
    !
    !     [X_1]     Y_11      Y_21      ....   Y_dim_X1
    !     [X_2]     Y_12      Y_22      ....   Y_dim_X2
    !     .         .         .
    !     .         .         .
    !     .         .         .
    !     [X_dim_Y] Y_1dim_Y  Y_2dim_Y  ....   Y_dim_Xdim_Y
    !
    CHARACTER(LEN=64), INTENT(IN) :: filename
    !
    REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: Y_array
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_X
    INTEGER(KIND = I4B), OPTIONAL, INTENT(IN) :: dim_Y
    !
    REAL(KIND = DP), DIMENSION(:), OPTIONAL, INTENT(IN) :: X_vector
    !
    LOGICAL, OPTIONAL, INTENT(IN) :: column_major_mode
    !
    ! Local variables
    INTEGER(KIND = I4B) :: X_index, Y_index, Y_dimension
    !
    ! Open file
    OPEN (UNIT = 1000, FILE = filename, STATUS = "UNKNOWN")
    !
    ! Define ordinate dimension
    IF (PRESENT(dim_Y)) THEN
       Y_dimension = dim_Y
    ELSE
       Y_dimension = dim_X
    ENDIF
    !
    ! Save results
    IF (PRESENT(X_vector)) THEN
       !
       IF (PRESENT(column_major_mode)) THEN
          !
          IF (column_major_mode) THEN
             WRITE(UNIT = 1000, FMT = *) X_vector(1:Y_dimension)
             DO X_index = 1, dim_X
                ! Column-major ordering
                WRITE(UNIT = 1000, FMT = *) Y_array(X_index, 1:Y_dimension)
             ENDDO
          ELSE  ! Row-major ordering
             DO Y_index = 1, Y_dimension
                WRITE(UNIT = 1000, FMT = *) X_vector(Y_index),  Y_array(1:dim_X, Y_index)
             ENDDO
          ENDIF
       ELSE ! Default case: column-major ordering
          WRITE(UNIT = 1000, FMT = *) X_vector(1:Y_dimension)
          DO X_index = 1, dim_X
             ! Column-major ordering
             WRITE(UNIT = 1000, FMT = *) Y_array(X_index, 1:Y_dimension)
          ENDDO
       ENDIF
       !
       ! No X_vector
    ELSE
       !
       IF (PRESENT(column_major_mode)) THEN
          !
          IF (column_major_mode) THEN
             !
             DO X_index = 1, dim_X
                ! Column-major ordering
                WRITE(UNIT = 1000, FMT = *) Y_array(X_index, 1:Y_dimension)
             ENDDO
          ELSE
             DO Y_index = 1, Y_dimension
                WRITE(UNIT = 1000, FMT = *) Y_array(1:dim_X, Y_index)
             ENDDO
          ENDIF
       ELSE ! Default case: column-major ordering
          DO X_index = 1, dim_X
             ! Column-major ordering
             WRITE(UNIT = 1000, FMT = *) Y_array(X_index, 1:Y_dimension)
          ENDDO
       ENDIF
       !
    ENDIF
    !
    CLOSE(UNIT = 1000)
    !
  END SUBROUTINE SAVE_OUTPUT_ARRAY
  !
  !
  !
END MODULE quartic_HO
