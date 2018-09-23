MODULE collide_cluster_module
    USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
    USE md_lj_cell_list_module, ONLY : allocate_arrays, deallocate_arrays, integrate, potential_type, force
    USE config_io_module, ONLY : read_config_atoms, snap_positions_vmd, time_stamp
    USE verlet_list_module, ONLY : finalize_list
    USE cluster_count_module, ONLY: cluster_count
    USE cluster_adhere_module, ONLY: cluster_adhere
    
    IMPLICIT NONE
    PRIVATE
    
    PUBLIC :: collide_cluster
    
    CONTAINS
    
        SUBROUTINE collide_cluster(r_cut,r_sep,v_com,dt,nsteps,integration_type,energy_snapshot,position_snapshot,&
            & cohesive_flag_num,cohesive,simulation_type,k_spring,lambda,clst_str,epsilon_modifier_flag,epsilon_modifier_strength)

            IMPLICIT NONE
  
            INTEGER*8                              :: nc1, nc2, nwall       ! Number of particles for each cluster or wall
            INTEGER*8                              :: n                     ! Total sum of particles
            INTEGER*8                              :: nsteps                ! Total number of steps
            INTEGER*8                              :: energy_snapshot       ! Snapshot size for energy measurement
            INTEGER*8                              :: position_snapshot     ! Snapshot size for energy measurement
            INTEGER*8                              :: i,j                     ! Dummy index
            INTEGER*8                              :: rnmeas = 0            ! Number of energy measurements made
            INTEGER*8                              :: cohesive_flag_num     ! Flag for cohesive in 1 or 0
            REAL*8, INTENT(in)                     :: r_cut                 ! Potential cut off distance
            REAL*8                                 :: dt                    ! Time step size

            REAL*8                                 :: energy = 0.0d0, energy_sq = 0.0d0
            REAL*8                                 :: sum_energy = 0.0d0, sum_energy_sq = 0.0d0
            REAL*8                                 :: sigma_energy = 0.0d0
    
            REAL*8                                 :: energy_kin = 0.0d0, energy_kin_sq = 0.0d0
            REAL*8                                 :: sum_energy_kin = 0.0d0, sum_energy_kin_sq = 0.0d0
            REAL*8                                 :: sigma_energy_kin = 0.0d0
    
            REAL*8                                 :: sum_energy_pot = 0.0d0, sum_energy_pot_sq = 0.0d0
            REAL*8                                 :: sigma_energy_pot = 0.0d0
    
            REAL*8                                 :: sum_energy_pot_int = 0.0d0, sum_energy_pot_int_sq = 0.0d0
            REAL*8                                 :: sigma_energy_pot_int = 0.0d0
    
            CHARACTER*15                           :: integration_type      ! Type of integration performed
            CHARACTER*9                            :: simulation_type       ! Simulation type: cluster to wall or cluster to cluster
            REAL*8,    DIMENSION(3)                :: r_com                 ! Dummy center of mass
            REAL*8,    DIMENSION(3), INTENT(in)    :: r_sep                 ! Start seperation between c.o.m. of clusters
            REAL*8,    DIMENSION(3)                :: r_mod                 ! Modified value of r_mod
            REAL*8,    DIMENSION(3)                :: v_com                 ! Start velocities between c.o.m. of clusters
            REAL*8,    DIMENSION(:,:), ALLOCATABLE :: r, v, f, r_start      ! Postion, Velocity, Force, and Displacement arrays
    
            REAL*8                                 :: cohesive              ! Cohesive LJ modification factor
            LOGICAL                                :: cohesive_flag         ! Cohesive flag

            REAL*8                                 :: energyloss_visc       ! Energy loss due to viscocity
            REAL*8                                 :: pbc_box               ! Size of pbc box for wall calculation
            REAL*8                                 :: k_spring              ! Spring constant for wall
            REAL*8                                 :: lambda                ! Viscocity coefficient for wall
    
            CHARACTER*33                           :: clst_str              ! Name of cluster
            CHARACTER*9                            :: clst_inp              ! Shortened name of cluster
            INTEGER*8                              :: perc, perc_flag       ! Percentage Complete counter
            
            INTEGER*8                              :: num_clust             ! Number of clusters after collision
            INTEGER*8, DIMENSION(:), ALLOCATABLE   :: npart_clust           ! Total number of particles in each cluster
            INTEGER*8, DIMENSION(:), ALLOCATABLE   :: clust_point           ! Pointer of cluster(i)
            REAL*8,    DIMENSION(:,:), ALLOCATABLE :: r_1                   ! Subset of r - particles belonging to initial cluster
            INTEGER*8                              :: biggest_cluster_index          ! Location of biggest cluster in pointer array
            INTEGER*8                              :: npart_adsorb       
            INTEGER*8                              :: biggest_cluster_npart      
            INTEGER*8                              :: biggest_cluster_start       
            INTEGER*8                              :: biggest_cluster_end
            
            REAL*8                                 :: r_adsorb
            REAL*8                                 :: r_adsorb_sq
            INTEGER*8                              :: r_temp
            REAL*8, DIMENSION(3)                   :: r_adsorb_i
            REAL*8, DIMENSION(3)                   :: r_adsorb_ij
            REAL*8                                 :: r_adsorb_ij_sq
            
            REAL*8                                 :: energy_interact,energy_kin_c1,energy_kin_c2,energy_kin_tot
            INTEGER*8                              :: adh_flag
            
            INTEGER*8                              :: epsilon_modifier_flag
            LOGICAL                                :: epsilon_flag
            REAL*8, OPTIONAL                       :: epsilon_modifier_strength
            REAL*8, DIMENSION(:), ALLOCATABLE      :: epsilon_mantle_array
            
    
            TYPE(potential_type) :: pot_total, pot_total_sq
            

            ! Check for the cohesive parameter modification
            IF (cohesive_flag_num == 1) THEN
                cohesive_flag = .TRUE.
            ELSE
                cohesive_flag = .FALSE.
            END IF
            
            ! See if the simulation is a cluster to cluster collision, or if it's something rlse
            IF(simulation_type == 'clst-clst') THEN
                
                WRITE(clst_inp,'(a9)') clst_str
                
                ! Read in initial configuration and allocate necessary arrays for simulation
                CALL read_config_atoms('.\Input\'//clst_inp,nc1)     ! Get n for cluster 1
                CALL read_config_atoms('.\Input\target',nc2)     ! Get n for cluster 2
                n = nc1 + nc2                                      ! Calculate total n 
                CALL allocate_arrays(n,r_cut,r,v,f)                ! Allocate combo array
                CALL read_config_atoms('.\Input\'//clst_inp,nc1,r,v)     ! Get positions and velocities of cluster 1
                CALL read_config_atoms('.\Input\target',nc2,r(:,nc1+1:n),v(:,nc1+1:n))     ! Get positions and velocities of cluster 2
    
                ! Ensure the c.o.m. of both clusters is zero
                r_com(:) = SUM(r(:,:nc1),dim=2)/DBLE(nc1) ! Cluster 1 c.o.m.
                r(:,:nc1) = r(:,:nc1) - SPREAD(r_com(:),dim=2,ncopies=nc1) ! Move cluster 1
                r_com(:) = SUM(r(:,nc1+1:n),dim=2)/DBLE(nc2) ! Cluster 2 c.o.m.
                r(:,nc1+1:n) = r(:,nc1+1:n) - SPREAD(r_com(:),dim=2,ncopies=nc2) ! Move cluster 2
    
                ! Move the c.o.m. of clusters to the chosen seperation distance
                r(:,:nc1) = r(:,:nc1) - SPREAD(r_sep(:)/2.0d0,dim=2,ncopies=nc1) ! Move cluster 1
                r(:,nc1+1:n) = r(:,nc1+1:n) + SPREAD(r_sep(:)/2.0d0,dim=2,ncopies=nc2) ! Move cluster 2
    
                ! Set translational c.o.m. of each cluster
                v(:,:nc1) = v(:,:nc1) + SPREAD(v_com(:)/2.0d0,dim=2,ncopies=nc1) ! Set cluster 1 c.o.m. velocity
                v(:,nc1+1:n) = v(:,nc1+1:n) - SPREAD(v_com(:)/2.0d0,dim=2,ncopies=nc2) ! Set cluster 2 c.o.m. velocity
        
                ! Prep the output XYZ files and Energy file
                OPEN(unit=10,file='.\Output\collision_'//clst_str//'.xyz',status='replace')
                CLOSE(unit=10)
                OPEN(unit=10,file='.\Output\collision_target_'//clst_str//'.xyz',status='replace')
                CLOSE(unit=10)
                OPEN(unit=10,file='.\Output\energy_'//clst_str//'.dat',status='replace')
                CLOSE(unit=10)
        
            ELSEIF(simulation_type == 'clst-wall') THEN
                
                WRITE(clst_inp,'(a9)') clst_str
                
                ! Read in initial configuration and allocate necessary arrays for simulation
                CALL read_config_atoms('.\Input\'//clst_inp,nc1)     ! Get n for cluster 1
                CALL read_config_atoms('.\Input\target',nwall)       ! Get n for the wall
                n = nc1 + nwall                                    ! Calculate total n 
                CALL allocate_arrays(n,r_cut,r,v,f,r_start)        ! Allocate combo array
                CALL read_config_atoms('.\Input\'//clst_inp,nc1,r,v) ! Get positions and velocities of cluster 1
                CALL read_config_atoms('.\Input\target',nwall,r(:,nc1+1:n),v(:,nc1+1:n))   ! Get positions and velocities of wall
        
                ! Ensure the c.o.m. of the cluster and the wall is zero
                r_com(:) = SUM(r(:,:nc1),dim=2)/DBLE(nc1) ! Cluster 1 c.o.m.
                r(:,:nc1) = r(:,:nc1) - SPREAD(r_com(:),dim=2,ncopies=nc1) ! Move cluster 1
                r_com(:) = SUM(r(:,nc1+1:n),dim=2)/DBLE(nwall) ! Wall c.o.m.
                r(:,nc1+1:n) = r(:,nc1+1:n) - SPREAD(r_com(:),dim=2,ncopies=nwall) ! Move wall

                ! Move the bottom middle of cluster 1 to the chosen seperation distance
                r_mod(3) = r_sep(3) - MINVAL(r(3,:nc1))
                r(:,:nc1) = r(:,:nc1) + SPREAD(r_mod(:),dim=2,ncopies=nc1) ! Move cluster 1

                ! Set translational c.o.m. of cluster 1
                v(:,:nc1) = v(:,:nc1) + SPREAD(v_com(:),dim=2,ncopies=nc1) ! Set cluster 1 c.o.m. velocity
        
                ! Set r_start array values
                r_start = r
        
                nc2 = nwall
        
                pbc_box = MAXVAL(r(1,nc1+1:n)) - MINVAL(r(1,nc1+1:n)) + (2.0d0**(1.0d0/6.0d0))
        
                IF (pbc_box < 2.0d0*r_cut) THEN
                    WRITE(unit=output_unit,fmt='(a,f4.1,a,f3.1,a)')'Error: r_cut x2 (',2.0d0*r_cut,') is larger than the PBC box (',pbc_box,').'
                    WRITE(unit=output_unit,fmt='(a)')'Press any key to exit.'
                    READ(unit=input_unit,fmt=*)
                    STOP
                END IF
        
                ! Prep the output XYZ files and Energy File
                OPEN(unit=10,file='.\Output\collision_'//clst_str//'.xyz',status='replace')
                CLOSE(unit=10)
                OPEN(unit=10,file='.\Output\collision_target_'//clst_str//'.xyz',status='replace')
                CLOSE(unit=10)
                OPEN(unit=10,file='.\Output\energy_'//clst_str//'.dat',status='replace')
                CLOSE(unit=10)
        
            END IF
    
    
            ! Perform initial checks
            WRITE(unit=output_unit,fmt='(a)')'Start collision time information'
            CALL time_stamp()
            WRITE(unit=output_unit,fmt='(a)')    
            
            ! Setup modified mantle potential if applicable
            IF (epsilon_modifier_flag == 1) THEN
                epsilon_flag = .TRUE.
                OPEN(unit=10,file='.\Input\cluster01mantle.dat',status='old')
                DO i=1,nc1
                    READ(unit=10,fmt=*) epsilon_mantle_array(i)
                END DO
                CLOSE(unit=10)
                IF (PRESENT(epsilon_modifier_strength) == .FALSE) THEN
                    epsilon_modifier_strength = 1000.0d0
                END IF
            ELSE
                epsilon_flag = .FALSE.
            END IF
            
            ! Call force routine to initialize make_list and setup initial energy values
            IF (epsilon_flag == .TRUE.) THEN
                CALL force(n,r_cut,r,f,pot_total,nc1,nc2,cohesive_flag,cohesive,r_start,v,k_spring,lambda,pbc_box,.TRUE.,.FALSE.,epsilon_modifier_strength,&
                & epsilon_mantle_array)
            ELSE
                CALL force(n,r_cut,r,f,pot_total,nc1,nc2,cohesive_flag,cohesive,r_start,v,k_spring,lambda,pbc_box,.TRUE.)
            END IF
            
            ! Record initial energy
            energy_kin = 0.5d0*SUM(v**2.0d0)
            OPEN(unit=10,file='.\Output\energy_'//clst_str//'.dat',status='old',position='append')
            WRITE(unit=10,fmt='(a,f15.8)')'Initial kinetic energy   = ',energy_kin
            WRITE(unit=10,fmt='(a,f15.8)')'Initial potential energy   = ',pot_total%pot
            WRITE(unit=10,fmt='(a,f15.8)')'Initial internal potential energy   = ',pot_total%cut
            WRITE(unit=10,fmt='(a,f15.8)')'Initial total energy   = ',energy_kin + pot_total%pot
            WRITE(unit=10,fmt=*)
            CLOSE(unit=10)
            energy_kin = 0.d0
    
            f = 0.0d0
            energyloss_visc = 0.0d0
            
            perc_flag = nsteps/100
    
            ! Run Simulation
            DO i=1,nsteps
        
                IF(simulation_type == 'clst-wall') THEN
                    CALL integrate(n,dt,r_cut,r,v,f,integration_type,simulation_type,pot_total,nc1,nc2,cohesive_flag,cohesive,r_start,k_spring,lambda,energyloss_visc,pbc_box)
                ELSEIF(cohesive_flag == .TRUE.) THEN
                    CALL integrate(n,dt,r_cut,r,v,f,integration_type,simulation_type,pot_total,nc1,nc2,cohesive_flag,cohesive)
                ELSE
                    CALL integrate(n,dt,r_cut,r,v,f,integration_type,simulation_type,pot_total)
                END IF
                

                IF (MOD(i,perc_flag) == 0) THEN
                    perc = i/perc_flag
                    WRITE(unit=output_unit,fmt='(i3,a25,a)') perc, '% complete on simulation ',clst_str
                END IF
                
                IF (MOD(i,energy_snapshot) == 0) THEN
                    rnmeas = rnmeas + 1.0d0
                    pot_total_sq%cut = (pot_total%cut)**2.0d0
                    pot_total_sq%pot = (pot_total%pot)**2.0d0
                    pot_total_sq%vir = (pot_total%vir)**2.0d0
                    pot_total_sq%lap = (pot_total%lap)**2.0d0
                    energy_kin = 0.5d0*SUM(v**2.0d0)
                    energy_kin_sq = energy_kin ** 2.0d0
                    energy = pot_total%pot + energy_kin
                    energy_sq = energy ** 2.0d0
            
                    ! Total Potential Energy
                    sum_energy_pot = sum_energy_pot + pot_total%pot
                    sum_energy_pot_sq = sum_energy_pot_sq + pot_total_sq%pot
            
                    ! Total Internal Potential Energy
                    sum_energy_pot_int = sum_energy_pot_int + pot_total%cut
                    sum_energy_pot_int_sq = sum_energy_pot_int_sq + pot_total_sq%cut
            
                    ! Total Kinetic Energy
                    sum_energy_kin = sum_energy_kin + energy_kin
                    sum_energy_kin_sq = sum_energy_kin_sq + energy_kin_sq
            
                    ! Total Energy
                    sum_energy = sum_energy + energy
                    sum_energy_sq = sum_energy_sq + energy_sq         
                END IF
        
                IF (MOD(i,position_snapshot) == 0 .and. simulation_type == 'clst-wall') THEN
                    WRITE(unit=output_unit,fmt='(a,i8)')'Position snap on step ',i
                    CALL snap_positions_vmd(nc1,nwall,n,r,clst_str,clst_str)    
                ELSEIF (MOD(i,position_snapshot) == 0) THEN
                    WRITE(unit=output_unit,fmt='(a,i8)')'Position snap on step ',i
                    CALL snap_positions_vmd(nc1,nc2,n,r,clst_str)    
                END IF
        
            END DO
            WRITE(unit=output_unit,fmt='(a)')
            WRITE(unit=output_unit,fmt='(a)')'Stop collision time information'
            CALL time_stamp()
            WRITE(unit=output_unit,fmt='(a)') 'Collision Complete!'
            WRITE(unit=output_unit,fmt='(a)')
    
            ! Record Final energies
            energy_kin = 0.5d0*SUM(v**2.0d0)
            CALL force(n,r_cut,r,f,pot_total,nc1,nc2,cohesive_flag,cohesive,r_start,v,k_spring,lambda,pbc_box)
            OPEN(unit=10,file='.\Output\energy_'//clst_str//'.dat',status='old',position='append')
            WRITE(unit=10,fmt='(a,f15.8)')'Final kinetic energy   = ',energy_kin
            WRITE(unit=10,fmt='(a,f15.8)')'Final potential energy   = ',pot_total%pot
            WRITE(unit=10,fmt='(a,f15.8)')'Final internal potential energy   = ',pot_total%cut
            WRITE(unit=10,fmt='(a,f15.8)')'Total energy loss due to viscocity   = ',energyloss_visc
            WRITE(unit=10,fmt='(a,f15.8)')'Final total energy   = ',energy_kin + pot_total%pot + energyloss_visc
            WRITE(unit=10,fmt=*)
            CLOSE(unit=10)
    
            ! Record averages
            OPEN(unit=10,file='.\Output\energy_'//clst_str//'.dat',status='old',position='append')
                sum_energy_kin = sum_energy_kin/rnmeas
                sum_energy_kin_sq = sum_energy_kin_sq/rnmeas
                sigma_energy_kin = DSQRT(sum_energy_kin_sq - (sum_energy_kin**2.0d0)) 
                sum_energy_pot = sum_energy_pot/rnmeas
                sum_energy_pot_sq = sum_energy_pot_sq/rnmeas
                sigma_energy_pot = DSQRT(sum_energy_pot_sq - (sum_energy_pot**2.0d0)) 
                sum_energy_pot_int = sum_energy_pot_int/rnmeas
                sum_energy_pot_int_sq = sum_energy_pot_int_sq/rnmeas
                sigma_energy_pot_int = DSQRT(sum_energy_pot_int_sq - (sum_energy_pot_int**2.0d0)) 
                sum_energy = sum_energy/rnmeas
                sum_energy_sq = sum_energy_sq/rnmeas
                sigma_energy = DSQRT(sum_energy_sq - (sum_energy**2.0d0))
                energyloss_visc = energyloss_visc/rnmeas
                WRITE(unit=10,fmt='(a,i15)')'# of measurements =   ',rnmeas
                WRITE(unit=10,fmt='(a,f15.8,a3,f15.8)')'Average kinetic energy               = ',sum_energy_kin,' +/- ',sigma_energy_kin
                WRITE(unit=10,fmt='(a,f15.8,a3,f15.8)')'Average potential energy             = ',sum_energy_pot,' +/- ',sigma_energy_pot
                WRITE(unit=10,fmt='(a,f15.8,a3,f15.8)')'Average internal potential energy    = ',sum_energy_pot_int,' +/- ',sigma_energy_pot_int
                IF(simulation_type == 'clst-wall') THEN
                    WRITE(unit=10,fmt='(a,f15.8)')     'Average energy loss due to viscocity = ',energyloss_visc
                END IF
                WRITE(unit=10,fmt='(a,f15.8,a3,f15.8)')'Average total energy                 = ',sum_energy,' +/- ',sigma_energy
        
                
            CLOSE(unit=10)
            
            IF (simulation_type == 'clst-wall') THEN
                ALLOCATE (npart_clust(nc1),clust_point(nc1))
                ALLOCATE (r_1(3,nc1))
                r_1 = r(:,:nc1)
                npart_adsorb = 0
                r_adsorb = 1.6
                r_adsorb_sq = r_adsorb**2.0d0
                CALL cluster_count(nc1,r_cut,r_1,num_clust,npart_clust,clust_point)
                OPEN(unit=10,file='.\Output\cluster_stats.dat',status='old',position='append')
                WRITE(unit=10,fmt='(a,f15.8,f15.8,f15.8,a,i16,a,i16)')'Distance: ',&
                    &r_sep(1),r_sep(2),r_sep(3),'       Total number of clusters found: ',&
                    &num_clust
                DO i=1,num_clust          
                    WRITE(unit=10,fmt='(a,i2,a,i8,a)')'Cluster ',i,' has ',npart_clust(i),' particles'
                END DO
                CLOSE(unit=10)
                
                biggest_cluster_index = MAXLOC(npart_clust,DIM=1)
                biggest_cluster_npart = npart_clust(biggest_cluster_index)
                biggest_cluster_end = 0
                
                DO i=1,biggest_cluster_index
                    biggest_cluster_end = biggest_cluster_end + npart_clust(i)
                END DO
                
                biggest_cluster_start = biggest_cluster_end - biggest_cluster_npart + 1
                
                IF (biggest_cluster_start == 0) THEN
                    biggest_cluster_start = 1
                END IF
            
                DO i=biggest_cluster_start,biggest_cluster_end,1
                    r_temp = clust_point(i)
                    r_adsorb_i(:) = r(:,r_temp)
                    
                    DO j=1,nc2
                        r_adsorb_ij(:) = r_adsorb_i(:) - r(:,nc1+j)
                        r_adsorb_ij_sq = SUM(r_adsorb_ij**2.0d0)
                
                        IF (r_adsorb_ij_sq < r_adsorb_sq) THEN
                            npart_adsorb = npart_adsorb + 1
                            EXIT
                        END IF
                        
                    END DO
                        
                END DO
                
                
                OPEN(unit=10,file='.\Output\cluster_data_saitoh_hayakawa.dat',status='old',position='append')
                WRITE(unit=10,fmt='(f15.8,i8,i8,i8)')DSQRT(SUM(v_com**2.0d0)),nc1,MAXVAL(npart_clust),npart_adsorb
                CLOSE(unit=10)
                DEALLOCATE (npart_clust,clust_point,r_1)
                
            END IF
            
            OPEN(unit=10,file='.\Output\cluster_data_adhere.dat',status='old',position='append')
            CALL cluster_adhere(nc1,nc2,r_sep(3),r_cut,cohesive,r,v,energy_interact,energy_kin_c1,energy_kin_c2,energy_kin_tot,adh_flag)
            WRITE(unit=10,fmt='(f15.8,f15.8,f15.8,f15.8,f15.8,i15)')DSQRT(SUM(v_com**2.0d0)),energy_interact,energy_kin_c1,energy_kin_c2,energy_kin_tot,adh_flag
            CLOSE(unit=10)
    
            CALL deallocate_arrays(r,v,f)
            CALL finalize_list
            RETURN
        
        END SUBROUTINE collide_cluster

END MODULE collide_cluster_module
