PROGRAM quartic_eigensystem_HO
  !
  ! Eigenvalues and eigenstates of a Hamiltonian with a quartic potential in a HO truncated basis
  !
  ! Potential:
  ! V_C(x) = A\, x^4 + B\, x^2 + C\,x^n~,
  !
  !\begin{align}
  ! \hat{\cal H}^{(n=1)} =& \epsilon_0 + \left( 3\alpha + \beta + \frac{1}{2}\right)\hat n
  ! + \frac{\alpha}{4}({a^\dagger}^4 + a^4) +  \frac{3\alpha}{2} {a^\dagger}^2 a^2 + \alpha
  ! ({a^\dagger} a^3 +{a^\dagger}^3 a)  \label{qhamn1}\\
  ! &+ \frac{6\alpha +2\beta -1}{4}({a^\dagger}^2 +  a^2)+ \frac{\gamma}{\sqrt{2}} (a^\dagger + a) ~~,\nonumber
  !  \hat{\cal H}^{(n=3)} =& \epsilon_3 + ( 3\alpha + 2\beta)\hat n
  !  + \frac{\alpha}{4}({a^\dagger}^4 + a^4) +  \frac{3\alpha}{2} {a^\dagger}^2 a^2 + \alpha
  !  ({a^\dagger} a^3 +{a^\dagger}^3 a) + \frac{3\alpha +2\beta -1}{2}{a^\dagger}^2 a^2 \label{qhamn3}\\
  !  &+ \frac{\gamma}{2\sqrt{2}} ({a^\dagger}^3 + a^3)+
  !  \frac{3\gamma}{2\sqrt{2}} ({a^\dagger}a^2+{a^\dagger}^2a) + \frac{3\gamma}{\sqrt{2}} (a^\dagger + a) ~~.\nonumber
  ! \end{align}
  ! \noindent where $\epsilon_0 = \frac{3\alpha + 2\beta +1}{4}$ and $\alpha = s^2A/2$, $\beta = B/2$, and $\gamma = s^{n-2}C/2$.
  !
  !
  USE nrtype
  !
  USE quartic_HO
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
#if _OPENMP
  USE OMP_LIB
#endif
  !
  IMPLICIT NONE
  !
  !
  REAL(KIND = DP) ::  alpha_value, beta_value, gamma_value ! Scaled Hamiltonian parameters
  !
  INTEGER(KIND = I4B) :: output_states ! Number of states saved (If 0 then info about all the states is saved)
  
  INTEGER(KIND = I4B) :: state_index, np 
  !
#ifndef  __GFORTRAN__
  !ifort Why is Z necessary???
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: Z
#endif 
  !
  !
  ! NAMELISTS
  NAMELIST/par_aux/ Iprint, eigenvec, excitation, save_output, output_states, convergence, dim_step, delta_energy
  NAMELIST/par_0/ HO_Dimension, benchmark, benchmark_total_iter
  NAMELIST/par_1/ s_value, A_value, B_value, C_value
  !
  !
  ! Read parameters
  !
  READ(UNIT = *, NML = par_aux)
  !
  IF (Iprint > 1) PRINT 10
  READ(UNIT = *, NML = par_0)
  !
  IF (Iprint > 1) PRINT 20
  READ(UNIT = *, NML = par_1)
  !
  !
  IF (Iprint > 1) THEN
     WRITE(UNIT = *, FMT = 5) Iprint, eigenvec, excitation, save_output, output_states
     IF (convergence) WRITE(UNIT = *, FMT = 6) convergence, dim_step, delta_energy
     WRITE(UNIT = *, FMT = 15) HO_Dimension, benchmark, benchmark_total_iter
     WRITE(UNIT = *, FMT = 25) s_value, A_value, B_value, C_value
  ENDIF
  !
  IF (save_output .AND. output_states == 0) output_states =  HO_Dimension
  !
  IF (benchmark) THEN
     Iprint = 0 ! Turn off display output
     ALLOCATE(benchmark_MFLOPS(1:benchmark_total_iter), STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "benchmark_MFLOPS allocation request denied."
        STOP
     ENDIF
  ELSE
     benchmark_total_iter = 1
  ENDIF
  !
  benchmarking: DO benchmark_iter = 1, benchmark_total_iter ! Benchmarking loop
     !
     ! Initialize time
     IF (benchmark) THEN
#if _OPENMP
        time_check_ref = OMP_GET_WTIME()
#else
        CALL CPU_TIME(time_check_ref)
#endif
     ENDIF
     !
     !
     ! Scale Hamiltonian parameters
     alpha_value =  0.5_DP*A_value*s_value**2
     beta_value = 0.5_DP*B_value
     gamma_value = 0.5_DP*C_value/s_value
     !
     IF (convergence) THEN
        ! Check eigenvalues convergence and fix new HO_dimension
        CALL Convergence_check(HO_Dimension, dim_step, delta_energy, output_states, &
             alpha_value, beta_value, gamma_value)
     ENDIF
     ! Allocate arrays
#ifndef  __GFORTRAN__
     CALL Allocate_arrays_sub(HO_Dimension, output_states, eigenvec, save_output)
#else
     CALL Allocate_arrays_sub(HO_Dimension, output_states, save_output)    
#endif
     ! BUILD HAMILTONIAN
     CALL QUARTIC_HAMILTONIAN_HO(HO_Dimension, Ham_quartic_mat, alpha_value, beta_value, gamma_value)
     !      
     ! Diagonalize Hamiltonian matrix (LAPACK95)
     IF (eigenvec) THEN
#ifdef  __GFORTRAN__
        ! gfortran
        CALL LA_SYEVR(A=Ham_quartic_mat, W=eigenval_vec, JOBZ='V', UPLO='U')
#else
        !ifort
        CALL SYEVR(A=Ham_quartic_mat, W=eigenval_vec, Z=Z, UPLO='U', IL=1, IU=HO_Dimension)
        !
        Ham_quartic_mat = Z
#endif           
     ELSE
#ifdef  __GFORTRAN__
        ! gfortran
        CALL LA_SYEVR(A=Ham_quartic_mat, W=eigenval_vec, JOBZ='N', UPLO='U')
#else
        !ifort
        CALL SYEVR(A=Ham_quartic_mat, W=eigenval_vec, UPLO='U')
#endif           
     ENDIF
     !
     IF (excitation) THEN
        !
        GS_energy = eigenval_vec(1)
        !
        eigenval_vec = eigenval_vec - GS_energy
        !
     ENDIF
     !
     !
     !  Display Output
     IF (Iprint > 0 .AND. .NOT. benchmark) THEN
        DO state_index = 1, output_states
           !
           np = state_index - 1_I4B
           !
           WRITE(UNIT = *, FMT = *) np, eigenval_vec(state_index)
           !
        ENDDO
     ENDIF
     !
     ! Save Output
     IF (save_output .AND. .NOT. benchmark) THEN
        ! Save eigenvalues
        file_name = output_filename(HO_dimension, "energies_")
        CALL Save_Output_Vector(filename = file_name, Y_vector = eigenval_vec, dim_X = output_states)
        !
        IF (eigenvec) THEN
           !
           ! Calculations
           CALL Calculations_Sc(HO_dimension, output_states, eigenval_vec, Ham_quartic_mat, &
                emean_vec, omegaeff_vec, ipr_vec, Husimi_ipr_vec, expected_n_vec, &
                expected_x_vec, expected_x2_vec)
           ! Save calculation results
           file_name = output_filename(HO_dimension, "emean_omegaeff_")
           CALL Save_Output_Vector(filename = file_name, X_vector = emean_vec, Y_vector = omegaeff_vec, &
                dim_X = output_states - 1)              ! Save mean w^eff vs energy
           file_name = output_filename(HO_dimension, "ipr_") 
           CALL Save_Output_Vector(filename = file_name, X_vector = eigenval_vec, Y_vector = ipr_vec, &
                dim_X = output_states)        ! Save ipr  vs energy
           file_name = output_filename(HO_dimension, "ipr_Husimi_") 
           CALL Save_Output_Vector(filename = file_name, X_vector = eigenval_vec, Y_vector = Husimi_ipr_vec, &
                dim_X = output_states)        ! Save Husimi ipr vs energy
           file_name = output_filename(HO_dimension, "expected_n_") 
           CALL Save_Output_Vector(filename = file_name, X_vector = eigenval_vec, Y_vector = expected_n_vec, &
                dim_X = output_states)        ! Save expectation value of n vs energy
           file_name = output_filename(HO_dimension, "expected_x_") 
           CALL Save_Output_Vector(filename = file_name, X_vector = eigenval_vec, Y_vector = expected_x_vec, &
                dim_X = output_states)        ! Save  expectation value of x vs energy
           file_name = output_filename(HO_dimension, "expected_x2_") 
           CALL Save_Output_Vector(filename = file_name, X_vector = eigenval_vec, Y_vector = expected_x2_vec, &
                dim_X = output_states)        ! Save expectation value of x2  vs enery
           ! Save eigenvectors
           ! eigenvectors saved in row-major order as
           ! E_i psi_i,1 psi_i,2 psi_i,3 ... psi_i,HO_dimension
           file_name = output_filename(HO_dimension, "eigenstates_")
           CALL Save_Output_Array(filename = file_name, Y_array = Ham_quartic_mat, &
                dim_X = HO_dimension, dim_Y = output_states, X_vector = eigenval_vec, &
                column_major_mode = .TRUE.)           
        ENDIF
        !
     ENDIF
     !
     !
     ! Deallocate arrays
#ifndef  __GFORTRAN__
     CALL Deallocate_arrays_sub(eigenvec, save_output)
#else
     CALL Deallocate_arrays_sub(save_output)
#endif
     !
     ! Check time
     IF (benchmark) THEN
#if _OPENMP
        time_check = OMP_GET_WTIME()
#else
        CALL CPU_TIME(time_check)
#endif
        IF (iprint > 1) &
             WRITE(UNIT = *, FMT = *) " Time :: Hamilt. building and diagonal. ", time_check - time_check_ref
        ! 5+12+5+8+12+8+5 = 55*dim flop hamiltonian building
        ! (4/3)*dim**3 flop Householder reflections to tridiagonal form
        ! 6*dim**3  QR algorithm 
        benchmark_MFLOPS(benchmark_iter) = &
             1.0D-6*(55.0_DP*HO_Dimension + (4.0_DP/3.0_DP + 6.0)*(1.0_DP*HO_Dimension)**3)/(time_check-time_check_ref)
        IF (Iprint > 0) WRITE(UNIT = *, FMT = *) time_check-time_check_ref, benchmark_MFLOPS(benchmark_iter)
     ENDIF
     !
  ENDDO benchmarking
  !
  IF (benchmark) THEN
     !
     ! Mean value
     mean_Mflops = SUM(benchmark_MFLOPS)/REAL(benchmark_total_iter, DP)
     !
     ! Standard deviation
     std_Mflops = SQRT(SUM((benchmark_MFLOPS - mean_Mflops)**2)/REAL(benchmark_total_iter-1, DP))
     !
     WRITE(unit = *, FMT = *) mean_Mflops, " +/- ", std_Mflops  
     !
     DEALLOCATE(benchmark_MFLOPS, STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "benchmark_MFLOPS deallocation request denied."
        STOP
     ENDIF
     !
     !
  ENDIF
  !
5 FORMAT(1X, " Iprint = ", I2, "; eigenvec = ", L2, "; excitation = ", L2, "; save_output = ", L2, &
       "; output_states = ", I6)
6 FORMAT(1X, " Convergence search on. ", " convergence = ", L2, "; dim_step = ", I4, "; delta_energy = ", ES14.7)
10 FORMAT(1X, "Reading  HO_Dimension, benchmark, benchmark_total_iter")
15 FORMAT(1X, "HO_Dimension = ", I6, " benchmark = ", L2, " benchmark_total_iter = ", I6)
20 FORMAT(1X, "Reading  kval, aval, bval")
25 FORMAT(1X, "s_value = ", ES14.7, "; A_value = ", ES14.7, /, &
        "B_value = ", ES14.7, "; C_value = ", ES14.7)
  !
  !
  !
END PROGRAM quartic_eigensystem_HO
