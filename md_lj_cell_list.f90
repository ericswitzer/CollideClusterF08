MODULE md_lj_cell_list_module
    USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit, error_unit
    
    IMPLICIT NONE
    PRIVATE
    
    PUBLIC :: allocate_arrays, deallocate_arrays, integrate, force
    
    TYPE, PUBLIC :: potential_type ! A composite variable for interactions comprising
        REAL    :: cut ! the potential energy cut (but not shifted) at r_cut and
        REAL    :: pot ! the potential energy cut-and-shifted at r_cut and
        REAL    :: vir ! the virial and
        REAL    :: lap ! the Laplacian and
    CONTAINS
        PROCEDURE :: add_potential_type
        GENERIC   :: OPERATOR(+) => add_potential_type
    END TYPE potential_type        
    
    CONTAINS    
    
        FUNCTION add_potential_type ( a, b ) RESULT (c)
            IMPLICIT NONE
            TYPE(potential_type)              :: c    ! Result is the sum of
            CLASS(potential_type), INTENT(in) :: a, b ! the two inputs
            c%cut = a%cut  +   b%cut
            c%pot = a%pot  +   b%pot
            c%vir = a%vir  +   b%vir
            c%lap = a%lap  +   b%lap
        END FUNCTION add_potential_type        
        
        SUBROUTINE allocate_arrays(n,r_cut,r,v,f,r_start)
            USE verlet_list_module, ONLY : initialize_list
            IMPLICIT NONE
            INTEGER*8, INTENT(in) :: n ! Number of particles
            REAL*8, INTENT(in) :: r_cut ! Potential cuttoff distance
            REAL*8, DIMENSION(:,:), ALLOCATABLE, INTENT(out) :: r, v, f ! Position, velocity, and force arrays
            REAL*8, DIMENSION(:,:), ALLOCATABLE, OPTIONAL, INTENT(out) :: r_start
            
            IF(PRESENT(r_start)) THEN
                ALLOCATE(r(3,n),v(3,n),f(3,n),r_start(3,n))
            ELSE
                ALLOCATE(r(3,n),v(3,n),f(3,n))
            END IF
            
            CALL initialize_list(n,r_cut)

        END SUBROUTINE allocate_arrays

        SUBROUTINE deallocate_arrays(r,v,f)
            IMPLICIT NONE
            REAL*8, DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: r, v, f
        
            DEALLOCATE(r,v,f)
        END SUBROUTINE deallocate_arrays
        
        SUBROUTINE integrate(n,dt,r_cut,r,v,f,integration_type,simulation_type,total,nc1,nc2,cohesive_flag,cohesive,r_start,k_spring,lambda,energyloss_visc,pbc_box)
            IMPLICIT NONE
            INTEGER*8, INTENT(in) :: n
            REAL*8, INTENT(in) :: r_cut
            REAL*8, INTENT(in):: dt        ! Time step size
            CHARACTER*15, INTENT(in) :: integration_type
            CHARACTER*9, INTENT(in) :: simulation_type
            TYPE(potential_type), INTENT(inout) :: total
            REAL*8, DIMENSION(:,:), ALLOCATABLE, OPTIONAL, INTENT(in) :: r_start
            REAL*8, DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: r, v, f
            INTEGER*8, OPTIONAL, INTENT(in) :: nc1
            INTEGER*8, OPTIONAL, INTENT(in) :: nc2
            REAL*8, OPTIONAL, INTENT(in) :: cohesive, k_spring, lambda, pbc_box
            LOGICAL, OPTIONAL, INTENT(in) :: cohesive_flag
            REAL*8, OPTIONAL, INTENT(inout) :: energyloss_visc
            REAL*8 :: dt_sq      ! dt^2
            REAL*8 :: dt_sq_half    ! 0.5*(dt^2)
            REAL*8 :: dt_half       ! 0.5*dt
            
            dt_half = 0.5d0*dt
            dt_sq = dt*dt
            dt_sq_half = 0.5d0*dt_sq

            
            IF (integration_type == 'Velocity Verlet') THEN
                IF (cohesive_flag == .TRUE. .and. simulation_type == 'clst-clst') THEN
                    r(:,:) = r(:,:) + dt * v(:,:) + dt_sq_half * f(:,:)                 ! Drift step
                    v(:,:) = v(:,:) + dt_half * f(:,:)                                  ! Kick half-step
                    CALL force(n,r_cut,r,f,total,nc1,nc2,cohesive_flag,cohesive)        ! Force evaluation
                    v(:,:) = v(:,:) + dt_half * f(:,:)                                  ! Kick half-step
                ELSEIF (cohesive_flag == .TRUE. .and. simulation_type == 'clst-wall') THEN
                    r(:,:) = r(:,:) + dt * v(:,:) + dt_sq_half * f(:,:)                    ! Drift step
                    r(:,nc1+1:n) = r(:,nc1+1:n) - pbc_box*DNINT(r(:,nc1+1:n)/pbc_box)      ! PBC conditions for the wall
                    v(:,:) = v(:,:) + dt_half * f(:,:)                                     ! Kick half-step
                    CALL force(n,r_cut,r,f,total,nc1,nc2,cohesive_flag,cohesive,r_start,v,k_spring,lambda,pbc_box) ! Force evaluation
                    v(:,:) = v(:,:) + dt_half * f(:,:)                                     ! Kick half-step
                    energyloss_visc = energyloss_visc + dt*SUM(v(:,nc1+1:n)**2.0d0)*lambda
                ELSE
                    r(:,:) = r(:,:) + dt * v(:,:) + dt_sq_half * f(:,:)                 ! Drift step
                    v(:,:) = v(:,:) + dt_half * f(:,:)                                  ! Kick half-step
                    CALL force(n,r_cut,r,f,total)                                       ! Force evaluation
                    v(:,:) = v(:,:) + dt_half * f(:,:)                                  ! Kick half-step
                END IF
            END IF
            
        END SUBROUTINE integrate
        
        SUBROUTINE force(n,r_cut,r,f,total,nc1,nc2,cohesive_flag,cohesive,r_start,v,k_spring,lambda,pbc_box,list_flag,check_flag,eps_strength,eps_array)
            USE verlet_list_module, ONLY : point, list, make_list
            IMPLICIT NONE
            INTEGER*8, INTENT(in) :: n
            REAL*8, INTENT(in) :: r_cut
            REAL*8, DIMENSION(:,:), ALLOCATABLE, INTENT(in) :: r
            REAL*8, DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: f
            REAL*8, DIMENSION(:,:), ALLOCATABLE, OPTIONAL, INTENT(in) :: r_start, v
            TYPE(potential_type), INTENT(out) :: total
            INTEGER*8, OPTIONAL, INTENT(in) :: nc1, nc2
            REAL*8, OPTIONAL, INTENT(in) :: cohesive, k_spring, lambda, pbc_box
            LOGICAL, OPTIONAL, INTENT(in) :: cohesive_flag,check_flag,list_flag
            REAL*8, DIMENSION(3,n) :: f_spring, f_visc
            REAL*8, OPTIONAL :: eps_strength
            REAL*8, DIMENSION(:), OPTIONAL :: eps_array
            
            
            INTEGER*8 :: i, j, k, m, index
            REAL*8 :: rij_sq, r_cut_sq, dri_sq
            REAL*8 :: sr2, sr6, sr12, pot_cut
            REAL*8, DIMENSION(3) :: rij, fij, dri
            REAL*8 :: pot_spring
            TYPE(potential_type) :: pair
            
            ! total%pot is the nonbonded cut-and-shifted potential energy for whole system
            ! total%cut is the nonbonded cut (but not shifted) potential energy for whole system
            ! total%vir is the corresponding virial
            ! total%lap is the corresponding Laplacian
            
            ! Calculate potential at cutoff
            sr2     = 1.0d0 / r_cut**2.0d0
            sr6     = sr2 ** 3.0d0
            sr12    = sr6 ** 2.0d0
            pot_cut = sr12 - sr6 ! Without numerical factor 4
        
            ! Calculate r cutoff squared
            r_cut_sq = r_cut**2.0d0
            
            ! Make Verlet Neighbor List
            IF(PRESENT(list_flag)) THEN
                WRITE(*,*)'Intentional reset of Verlet List'
                CALL make_list(n,r,.TRUE.)
            ELSE
                CALL make_list(n,r)
            END IF
            
            f = 0.0d0
            f_spring = 0.0d0
            f_visc = 0.0d0
            pot_spring = 0.0d0
            
            total = potential_type(pot=0.0d0,cut=0.0d0,vir=0.0d0,lap=0.0d0)
            
            IF(PRESENT(check_flag)) THEN
                OPEN(unit=10,file='.\Output\check.dat',status='old',position='append')
                WRITE(unit=10,fmt='(a)')'Checking pointer array...'
                DO m=1,n
                    WRITE(10,*)'point ',m,'value ',point(m)
                END DO
                WRITE(10,*)'End of Check'
                WRITE(10,*)''
                CLOSE(10)
            END IF
            
            IF(cohesive_flag == .TRUE. .and. PRESENT(k_spring)) THEN
                DO i=1,n-1 ! Begin outer loop over atoms
                    DO k=point(i),point(i+1)-1 ! Begin inner loop over neighbor atoms (if any)
                        
                        j = list(k)                            
                        rij(:) = r(:,i) - r(:,j) ! Seperation vector
                        
                        IF(i >= nc1+1 .and. j >= nc1+1) THEN
                            rij(:) = rij(:) - pbc_box*DNINT(rij(:)/pbc_box) ! Mirror Image on wall particles only
                        END IF
                        
                        rij_sq = SUM(rij**2.0d0) ! Squared seperation
                    
                        IF(rij_sq < r_cut_sq) THEN ! Check within cutoff
                            sr2 = 1.0d0 / rij_sq
                            sr6 = sr2 ** 3.0d0
                            sr12 = sr6 ** 2.0d0
                            IF(i <= nc1 .and. j > nc1) THEN
                                pair%cut = sr12 - (cohesive*sr6) ! LJ pair potential (cut but not shifted)
                            ELSE
                                pair%cut = sr12 - sr6
                            END IF
                            pair%vir = pair%cut + sr12 ! LJ pair virial
                            pair%pot = pair%cut - pot_cut ! LJ pair potential (cut-and-shifted)
                            pair%lap = (22.0d0 * sr12 - 5.0d0 * sr6) * sr2 ! LJ pair Laplacian
                            fij = rij * pair%vir * sr2 ! LJ pair forces
                            total = total + pair
                            f(:,i) = f(:,i) + fij
                            f(:,j) = f(:,j) - fij
                        END IF
                        
                    END DO
                                                
                    ! Wall's spring force and energy relaxation (cfr Saitoh and Hayakawa 2009).
                    ! For the potential, we only calculate the potential added by the spring
                    ! force. I have not accounted for the viscocity potential loss, nor the impact
                    ! on the Laplacian or Virial
                    IF(i >= nc1+1) THEN
                        dri(:) = r(:,i) - r_start(:,i) 
                        dri_sq = SUM(dri**2.0d0)
                        f_spring(:,i) = f_spring(:,i) - (k_spring*dri(:))
                        f_visc(:,i) = f_visc(:,i) - (lambda*v(:,i))
                        pot_spring = pot_spring + 0.5d0*k_spring*dri_sq
                    END IF
                    
                END DO
                    
            ELSEIF (cohesive_flag == .TRUE.) THEN
                DO i=1,n-1 ! Begin outer loop over atoms
                    DO k=point(i),point(i+1)-1 ! Begin inner loop over neighbor atoms (if any)
                        j = list(k)
                        rij(:) = r(:,i) - r(:,j) ! Seperation vector
                        rij_sq = SUM(rij**2.0d0) ! Squared seperation
                    
                        IF(rij_sq < r_cut_sq) THEN ! Check within cutoff
                            sr2 = 1.0d0 / rij_sq
                            sr6 = sr2 ** 3.0d0
                            sr12 = sr6 ** 2.0d0
                            IF(i <= nc1 .and. j > nc1) THEN
                                pair%cut = sr12 - (cohesive*sr6) ! LJ pair potential (cut but not shifted)
                            ELSE
                                pair%cut = sr12 - sr6
                            END IF
                            pair%vir = pair%cut + sr12 ! LJ pair virial
                            pair%pot = pair%cut - pot_cut ! LJ pair potential (cut-and-shifted)
                            pair%lap = (22.0d0 * sr12 - 5.0d0 * sr6) * sr2 ! LJ pair Laplacian
                            fij = rij * pair%vir * sr2 ! LJ pair forces
                            total = total + pair
                            f(:,i) = f(:,i) + fij
                            f(:,j) = f(:,j) - fij
                            
                            IF
                            DO i=1,n
                                IF(eps_array(i) \=0) THEN
                                    index = eps_array(i)
                                    f(:,index) = eps_strength*f(:,i)
                                ELSE IF
                                    
                            END DO
                        END IF
                    END DO
                END DO
            ELSE
                DO i=1,n-1 ! Begin outer loop over atoms
                    DO k=point(i),point(i+1)-1 ! Begin inner loop over neighbor atoms (if any)
                        j = list(k)                   
                        rij(:) = r(:,i) - r(:,j) ! Seperation vector
                        rij_sq = SUM(rij**2.0d0) ! Squared seperation
                    
                        IF(rij_sq < r_cut_sq) THEN ! Check within cutoff
                            sr2 = 1.0d0 / rij_sq
                            sr6 = sr2 ** 3.0d0
                            sr12 = sr6 ** 2.0d0
                            pair%cut = sr12 - sr6 ! LJ pair potential (cut but not shifted)
                            pair%vir = pair%cut + sr12 ! LJ pair virial
                            pair%pot = pair%cut - pot_cut ! LJ pair potential (cut-and-shifted)
                            pair%lap = (22.0d0 * sr12 - 5.0d0 * sr6) * sr2 ! LJ pair Laplacian
                            fij = rij * pair%vir * sr2 ! LJ pair forces
                            total = total + pair
                            f(:,i) = f(:,i) + fij
                            f(:,j) = f(:,j) - fij
                        END IF
                    END DO
                END DO
                    
            END IF
            
            ! If Epsilon modifier exists, apply the appropriate force terms
            
            ! Multiply results by numerical factors
            f         = f         * 24.0d0       ! 24*epsilon
            f         = f + f_spring + f_visc  ! Force of spring and viscocity
            total%cut = total%cut * 4.0d0        ! 4*epsilon
            total%cut = total%cut + pot_spring
            total%pot = total%pot * 4.0d0        ! 4*epsilon
            total%pot = total%pot + pot_spring
            total%vir = total%vir * 24.0d0 / 3.0d0 ! 24*epsilon and divide virial by 3
            total%lap = total%lap * 24.0d0 * 2.0d0 ! 24*epsilon and factor 2 for ij and ji
        
        END SUBROUTINE force
    
END MODULE md_lj_cell_list_module