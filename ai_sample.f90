program sample
    use pcg_module 
    use civilization_mod
    implicit none 

    real :: start_time, end_time 
    integer :: i, n, j, tgt 
    integer, parameter :: N_pt = 5000, No_of_sims = 1
    real(8), parameter :: H_fixed = 0.001
    integer :: occupied
    integer, dimension(1) :: min_init_distance_id
    integer, dimension(1) :: min_dist_id
    integer, dimension(1) :: l_min_distances_id 
    real(8) :: min_distance, l_min_distance
    real(8), dimension(N_pt) :: init_distance
    real(8), allocatable :: minimum_distances(:)
    real(8) :: minimum_distance
    real(8), allocatable :: dist_n(:) 
    type(civil), allocatable :: civil_xyz(:)
    real(8) :: ti, tj, clock, H
    real(8) :: ti_matter, tj_matter, k
    real(8) :: random_number
    integer :: failed_counter
    integer(8) :: no_of_sims_counter
    real(8) :: cpath, t_total, std_t_statistic
    integer(8) :: seed
    integer :: possible_targets 

    ! Random number array allocation and index tracking
    real(8), allocatable :: random_number_array(:)
    integer :: rna_size, rna_counter, rn_checker
    real(8) :: r

    call cpu_time(start_time)

    cpath = 0.0_8
    t_total = 0.0_8 
    std_t_statistic = 0.0_8

    print*, "Starting simulation..."

    ! Allocate memory
    allocate(civil_xyz(N_pt))
    allocate(dist_n(N_pt))
    allocate(minimum_distances(N_pt))
    allocate(random_number_array(500000)) ! Large buffer to prevent overflow

    open(unit = 10, file = 'sims=1,k=vary', action = 'write', status = 'replace')

    do no_of_sims_counter = 1, No_of_sims   

        ! Reset parameters for each simulation
        seed = 48
        occupied = 0
        failed_counter = 0
        ti = 0.0_8
        tj = 0.0_8
        ti_matter = 0.0_8
        tj_matter = 0.0_8
        clock = 0.0_8
        n = 1
        H = 0.06_8
        k = 0.80_8
        min_init_distance_id(:) = 0
        minimum_distance = 0.0_8

        call pcg_init(seed, seed)
        call initialize_universe(civil_xyz, N_pt, init_distance)

        ! 1. Generate random numbers until N_pt successful values (>= k) exist
        rna_size = 0
        rn_checker = 0
        do while (rn_checker < N_pt)
            rna_size = rna_size + 1
            if (rna_size > size(random_number_array)) then
                print*, "Error: Exceeded random array allocation size."
                stop
            end if
            r = pcg_random_real()
            random_number_array(rna_size) = r 
            if (r >= k) rn_checker = rn_checker + 1
        end do
        
        rna_counter = 1 ! Fortran 1-based indexing
        print*, "Generated random numbers. Total size:", rna_size, "Successes (>= k):", rn_checker

        ! 2. Locate starting civilization closest to center
        min_init_distance_id = minloc(init_distance(1:N_pt))
        civil_xyz(min_init_distance_id(1))%status = 1
        occupied = 1 

        ! Main simulation loop
        do while (occupied < N_pt)
            
            minimum_distances = huge(1.0_8)
            possible_targets = count(civil_xyz(:)%status == 0)

            ! Step 1: Scan active colonies (status 1) and assign targets
            do i = 1, N_pt
                if (civil_xyz(i)%status == 1 .and. possible_targets > 0) then 
                    call find_neighbor(i, N_pt, civil_xyz, min_dist_id, min_distance)

                    ! Deadlock recovery: Clear origin history if all local neighbors are blocked
                    if (min_dist_id(1) == 0 .or. min_distance >= huge(1.0_8) / 2.0_8) then
                        do j = 1, N_pt
                            if (civil_xyz(j)%status == 0) civil_xyz(j)%origin_history = 0
                        end do
                        call find_neighbor(i, N_pt, civil_xyz, min_dist_id, min_distance)
                    end if

                    ! Reserve target if valid neighbor found
                    if (min_dist_id(1) > 0 .and. min_distance < huge(1.0_8) / 2.0_8) then
                        tgt = min_dist_id(1)
                        civil_xyz(i)%target = tgt
                        civil_xyz(i)%status = 2        ! Travelling
                        civil_xyz(tgt)%status = 3     ! Reserved target
                        civil_xyz(i)%total_distance = civil_xyz(i)%total_distance + min_distance
                        civil_xyz(i)%remaining_distance = min_distance
                        minimum_distances(i) = min_distance
                    end if

                else if (civil_xyz(i)%status == 2) then 
                    minimum_distances(i) = civil_xyz(i)%remaining_distance
                end if 
            end do 
            
            ! Step 2: Determine global minimum distance step
            l_min_distances_id = minloc(minimum_distances(1:N_pt))
            l_min_distance = minimum_distances(l_min_distances_id(1))

            ! Exit guard only when no ships are in flight and no moves remain
            if (l_min_distance >= huge(1.0_8) / 2.0_8) exit

            ! Step 3: Advance simulation clock
            clock = clock + l_min_distance
            tj = (-1.0_8/H)*log(exp(-H*ti) - H*l_min_distance)
            tj_matter = (2.0_8/(3.0_8*H))*((1.0_8 + (3.0_8/2.0_8)*H*ti_matter)**(1.0_8/3.0_8) + H*l_min_distance/2.0_8)**3.0_8 - (2.0_8/(3.0_8*H))
            cpath = cpath + l_min_distance

            ! Step 4: Advance travelling ships
            do i = 1, N_pt 
                if (civil_xyz(i)%status == 2) then 
                    civil_xyz(i)%remaining_distance = civil_xyz(i)%remaining_distance - l_min_distance
                    civil_xyz(i)%total_distance = civil_xyz(i)%total_distance + l_min_distance

                    if (civil_xyz(i)%remaining_distance <= 1.0e-10_8) then 
                        random_number = random_number_array(rna_counter)
                        rna_counter = rna_counter + 1 
                        civil_xyz(i)%status = 1
                        minimum_distances(i) = huge(1.0_8) 
                        
                        tgt = civil_xyz(i)%target

                        if (random_number >= k) then
                            civil_xyz(tgt)%status = 1 
                            occupied = occupied + 1 
                        else
                            civil_xyz(tgt)%status = 0
                            failed_counter = failed_counter + 1
                            civil_xyz(i)%counter = civil_xyz(i)%counter + 1 
                            
                            ! Record origin i into target tgt's history
                            civil_xyz(tgt)%origin_history(1) = i
                        end if 
                        
                        civil_xyz(i)%target = 0 ! Clear active target reference
                    else
                        minimum_distances(i) = civil_xyz(i)%remaining_distance
                    end if 
                end if 
            end do 
          
            ti = tj 
            ti_matter = tj_matter
             
        end do 
        
        write(10,*) k, clock
        print*, k, clock

    end do 

    call cpu_time(end_time)
    
    print*, "DONE"
    print*, "k = ", k, " Success occupation: ", occupied, ' Failed occupation: ', failed_counter
    print*, "Number of unoccupied sites: ", count(civil_xyz(:)%status == 0 )
    print*, "Number of occupied sites: ", count(civil_xyz(:)%status == 1 )
    print*, "Number of travelling: ", count(civil_xyz(:)%status == 2 )
    print*, "Number of target: ", count(civil_xyz(:)%status == 3 )
    print*, 'Time elapsed: ', clock

    close(10)
end program sample