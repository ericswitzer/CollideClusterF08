MODULE cluster_adhere_module
    USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, error_unit, output_unit, iostat_end, iostat_eor
    IMPLICIT NONE
    PRIVATE

    PUBLIC :: cluster_adhere
    
CONTAINS
                           
    SUBROUTINE cluster_adhere(npart_cluster,npart_target,r_sepz,r_cut,cohesive,r,v,energy_interact,energy_kinetic_cluster,energy_kinetic_target,energy_kinetic_total,adhere_flag)
	
	IMPLICIT NONE
    
    INTEGER*8, INTENT(in)               :: npart_cluster,npart_target
    REAL*8, INTENT(in)                  :: r_cut, cohesive
    REAL*8, DIMENSION(:,:), INTENT(in)  :: r, v
    INTEGER*8, INTENT(out)              :: adhere_flag
    REAL*8, INTENT(out)                 :: energy_interact, energy_kinetic_cluster, energy_kinetic_target, energy_kinetic_total
    
    INTEGER*8                           :: i, j
    REAL*8                              :: rij_sq, r_cut_sq, sr2, sr6, sr12, energy_pair
    REAL*8, DIMENSION(3)                :: rij
    
    REAL*8                              :: src2, src6, src12, pot_cut
    
    REAL*8                              :: r_comz, r_sepz
    
    r_cut_sq = r_cut**2.0d0
    
    energy_interact = 0.0d0
    energy_kinetic_cluster = 0.0d0
    energy_kinetic_target = 0.0d0
    energy_kinetic_total = 0.0d0
    adhere_flag = 0
    
    ! Calculate potential at cutoff
    src2     = 1.0d0 / r_cut**2.0d0
    src6     = src2 ** 3.0d0
    src12    = src6 ** 2.0d0
    pot_cut = src12 - src6 ! Without numerical factor 4
    
    ! Calculate interaction energy
    DO i=1,npart_cluster
        DO j=1,npart_target
            rij(:) = r(:,i) - r(:,npart_cluster+j)
            rij_sq = SUM(rij**2.0d0)
            IF(rij_sq < r_cut_sq) THEN ! Check within cutoff
                sr2 = 1.0d0 / rij_sq
                sr6 = sr2 ** 3.0d0
                sr12 = sr6 ** 2.0d0
                energy_pair = sr12 - (cohesive*sr6)
                energy_pair = energy_pair - pot_cut
                energy_interact = energy_interact + energy_pair
            END IF
        END DO
    END DO
    energy_interact = energy_interact * 4.0d0
    
    ! Calculate Kinetic energy of cluster and target
    energy_kinetic_cluster = 0.5d0*SUM(v(:,:npart_cluster)**2.0d0)
    energy_kinetic_target = 0.5d0*SUM(v(:,npart_cluster+1:)**2.0d0)
    energy_kinetic_total = energy_kinetic_cluster + energy_kinetic_target
    
    ! Set adherance flag
    r_comz = SUM(r(3,:npart_cluster))/DBLE(npart_cluster) ! Cluster's z-axis c.o.m.
    IF (r_comz <= 0.8d0*r_sepz) THEN
        adhere_flag = 1
    ELSE
        adhere_flag = 0
    END IF
    
    
    RETURN
    
    END SUBROUTINE cluster_adhere

END MODULE cluster_adhere_module