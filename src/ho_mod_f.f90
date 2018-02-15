MODULE ho_mod_f
  !
  USE nrtype
  !
  IMPLICIT NONE
  !
  !
  ! NDX -> XGRID DIMENSION
  INTEGER(KIND = I4B) :: dim_X
  !
  ! Xmin Xmax and delta_X step
  REAL(KIND = DP) :: X_min, X_max, Delta_X
  ! XGRID
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: X_Grid      
  ! HARMONIC BASIS
  REAL(KIND = DP),  DIMENSION(:,:), ALLOCATABLE :: Har_Bas_X, Avec_Har_X
  !
END MODULE ho_mod_f
