MODULE ho_mod_f
  !
  USE nrtype
  !
  IMPLICIT NONE
  !
  ! Eigenvectors X dependence calculation
  LOGICAL :: avec_X_calc ! If .T. compute x dependence of a set of eigenvectors
  !
  ! NDX -> XGRID DIMENSION
  INTEGER(KIND = I4B) :: dim_X
  !
  ! Number of eigenstates whose X dependence is saved.
  INTEGER(KIND = I4B) :: avec_X_states
  !
  ! Xmin Xmax and delta_X step
  REAL(KIND = DP) :: X_min, X_max, Delta_X
  ! XGRID
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: X_Grid      
  ! HARMONIC BASIS
  REAL(KIND = DP),  DIMENSION(:,:), ALLOCATABLE :: Har_Bas_X, Avec_Har_X
  !
CONTAINS
  !
  SUBROUTINE EIGENVEC_X(dim, s_value, alpha_value, beta_value, gamma_value, delta_value, &
       eigenval_vec, eigenvec_mat, Iprint)
    !
    ! Subroutine to build and save the space dependence of the first avec_X_states  eigenvectors
    !
    IMPLICIT NONE
    !
    ! Dimension of the Harmonic Basis
    INTEGER(KIND = I4B), INTENT(IN) :: dim
    ! Length Scale and scaled Hamiltonian parameters of the problem
    REAL(KIND = DP), INTENT(IN) ::  s_value, alpha_value, beta_value, gamma_value, delta_value
    ! Eigenvalues matrix
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: eigenval_vec
    ! Eigenvectors matrix
    REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: eigenvec_mat
    ! Verbosity control Iprint > 1 ==> normalization test
    INTEGER(KIND = I4B), INTENT(IN) :: Iprint
    !
    ! Local variables
    INTEGER(KIND = I4B) :: Ierr, state_index, X_index
    REAL(KIND = DP), PARAMETER ::wf_scaling = 5.0_DP
    !
    dim_X = dim_X + 2 ! TO ACCOMODATE X_min AND X_max
    !
    ! DEFINE X GRID (unitless length)
    ALLOCATE(X_grid(1:dim_X), STAT = Ierr)    
    IF (Ierr /= 0) THEN
       PRINT*, "X_grid allocation request denied."
       STOP
    ENDIF
    !
    !
    Delta_X = (X_max-X_min)/REAL(dim_X - 1, DP)  
    !
    X_grid = X_min + Delta_X*REAL((/ (X_index, X_index = 0, dim_X-1) /),DP)
    !
    !
    ! Allocate and define Harmonic Basis
    ALLOCATE(Har_Bas_X(1:dim_X, 1:dim), STAT = Ierr)    
    IF (Ierr /= 0) THEN
       PRINT*, "Har_Bas_X allocation request denied."
       STOP
    ENDIF
    !
    CALL HO_1D_BASIS(s_value, Dim, Iprint) ! Oscillator Length = 1.0
    !
    ! Allocate and define Spatial dependence of eigenvectors
    ALLOCATE(Avec_Har_X(1:dim_X, 1:avec_X_states), STAT = Ierr)    
    IF (Ierr /= 0) THEN
       PRINT*, "Avec_Har_X allocation request denied."
       STOP
    ENDIF
    !
    DO state_index = 1, avec_X_states
       DO X_index = 1, dim_X
          avec_Har_X(X_index,state_index) = DOT_PRODUCT(eigenvec_mat(1:Dim,state_index), Har_Bas_X(X_index,1:Dim))
       ENDDO
    ENDDO
    !
    ! Save eigenvectors (dirty hack)
    DO X_index = 1, dim_X
       !
       write(30,*) X_grid(X_index), (avec_Har_X(X_index, state_index), state_index = 1, avec_X_states)
       !
       write(35,*) X_grid(X_index), &
            alpha_value*X_grid(X_index)**4 + beta_value*X_grid(X_index)**2 + &
            gamma_value*X_grid(X_index) + delta_value*X_grid(X_index)**3, &
            (eigenval_vec(state_index) + &
            wf_scaling*avec_Har_X(X_index, state_index), state_index = 1, avec_X_states)
       !
       write(40,*) X_grid(X_index), &
            alpha_value*X_grid(X_index)**4 + beta_value*X_grid(X_index)**2 + &
            gamma_value*X_grid(X_index) + delta_value*X_grid(X_index)**3, &
            (eigenval_vec(state_index) + &
            wf_scaling*avec_Har_X(X_index, state_index)**2, state_index = 1, avec_X_states)
       !
    ENDDO
    !
    ! Deallocate Harmonic basis and spatial eigenvectors
    DEALLOCATE(Har_Bas_X, STAT = Ierr)    
    IF (Ierr /= 0) THEN
       PRINT*, "Har_Bas_X deallocation request denied."
       STOP
    ENDIF
    !
    DEALLOCATE(Avec_Har_X, STAT = Ierr)    
    IF (Ierr /= 0) THEN
       PRINT*, "Avec_Har_X deallocation request denied."
       STOP
    ENDIF
    !
  END SUBROUTINE EIGENVEC_X
  !
  SUBROUTINE HO_1D_BASIS(apar, NdimH, Iprint)
    !     
    !     COMPUTES A 1D HARMONIC BASIS
    !
    !     INPUT  :: X_GRID  --> VECTOR WITH X GRID
    !               apar    --> LENGTH SCALE OF THE PROBLEM
    !               NdimH   --> DIMENSION OF THE HARMONIC BASIS
    !
    !
    !     OUTPUT :: HAR_BAS_X  --> MATRIX WITH HARMONIC BASIS
    !
    !     FORMAT :: Iprint  --> VERBOSITY CONTROL
    !
    !     Taken from B2UCoupling with Pseudostates source code
    !     ~/PR_Fortran/B2UCoupling/Padova/1D/One_Body/Pseudostates_Calc/src
    !
    !     by Currix TM.
    !
    !
    IMPLICIT NONE
    !
    ! ARGUMENTS
    INTEGER(KIND = I4B), INTENT(IN) :: NdimH, Iprint
    REAL(KIND = DP), INTENT(IN) :: apar 
    !
    !
    ! OTHER VARIABLES
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: HO_norm_test
    REAL(KIND = DP) :: PI14, apar2, ainteg
    INTEGER(KIND = I4B) :: kx, Ierr
    !
    IF (Iprint > 2) PRINT*, "BUILDING HARMONIC BASIS"
    !
    PI14 = SQRT(SQRT(PI_D))
    apar2 = apar*apar
    !
    !
    !     HO n = 0
    HAR_BAS_X(:,1) = SQRT(apar)/(PI14)*EXP(-apar2*X_GRID*X_GRID/2.0_DP)
    !
    !     HO n = 1
    IF (NdimH > 1) HAR_BAS_X(:,2) = (SQRT(apar)/(PI14))*SQRT(2.0_DP)*apar*X_GRID*EXP(-apar2*X_GRID*X_GRID/2.0_DP)
    !
    !    RECURRENCE RELATION (WATCH OUT THE INDEXES :: MATRIX START AT 1, NOT 0)
    DO kx = 2, NdimH-1
       HAR_BAS_X(:,kx+1) = &
            SQRT(2.0_DP/(1.0_DP*kx))*apar*X_GRID*HAR_BAS_X(:,kx) - &
            SQRT((1.0_DP*(kx-1))/(1.0_DP*kx))*HAR_BAS_X(:,kx-1)
    ENDDO
    !
    !     TESTING NORMALIZATION 
    IF (Iprint >= 1) THEN
       !
       ! DEFINE HO_norm_test
       ALLOCATE(HO_norm_test(1:dim_X), STAT = Ierr)    
       IF (Ierr /= 0) THEN
          PRINT*, "HO_norm_test allocation request denied."
          STOP
       ENDIF
       !
       DO kx = 1, NdimH
          HO_norm_test = HAR_BAS_X(:,kx)**2
          ! 
          ! Integrate with simpsn intlib subroutine
          CALL simpsn(dim_X, Delta_X, HO_norm_test, ainteg)
          !
          WRITE(*,*) "HARMONIC FUNCTION ", kx, " NORMALIZATION", ainteg
       ENDDO
       !
       DEALLOCATE(HO_norm_test, STAT = Ierr)    
       IF (Ierr /= 0) THEN
          WRITE(*,*) "HO_norm_test deallocation request denied."
          STOP
       ENDIF
       !
       IF (Iprint > 2) WRITE(*,*) "DONE"
       !
    ENDIF
    !     
    !     
  END SUBROUTINE HO_1D_BASIS
  !
END MODULE ho_mod_f
