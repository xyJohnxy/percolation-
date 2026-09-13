
program sample
    use pcg_module 
    use civilization_mod
    implicit none 

    real :: start_time, end_time 
    integer :: i, n, j 
    integer, parameter :: N_pt = 5000, No_of_sims = 1
    real(8), parameter :: H_fixed = 0.001 !, k = 0.6 
    integer :: occupied
    integer, dimension(1) :: min_init_distance_id
    integer, dimension(1) :: min_dist_id
    integer, dimension(1) :: l_min_distances_id 
    real(8) :: min_distance, l_min_distance
    real(8), dimension(N_pt) :: init_distance
    real(8), allocatable :: minimum_distances(:)
    ! real(8), dimension(No_of_sims) :: t_terminate_static  !array that contains the end times for each simulation
    real(8) :: minimum_distance
    real(8), allocatable :: dist_n(:) 
    type(civil), allocatable :: civil_xyz(:)
    real(8) :: ti, tj, clock, H
    real(8) :: ti_matter, tj_matter, k
    real(8) :: random_number !death rate, random_number
    integer :: failed_counter
    integer(8) :: no_of_sims_counter
    real(8) :: cpath, t_total, std_t_statistic
    integer(8) :: seed, possible_targets 

    !generate a dynamic array of containing the random numbers 
    real(8), allocatable :: random_number_array(:)
    integer :: rna_size = 0, rna_counter = 1, rn_checker = 0 
    real(8) :: r


    call cpu_time(start_time)

    cpath = 0
    t_total = 0  
    std_t_statistic = 0

    print*, "working"
    ! Allocate array sizes 
    allocate(civil_xyz(N_pt))
    allocate(dist_n(N_pt))
    allocate(minimum_distances(N_pt))
    allocate(random_number_array(100000))

    random_number_array = 0

    minimum_distances = huge(1.0_8) 

    open(unit = 10, file = 'sims=1,k=vary', action = 'write', status = 'replace')

    do no_of_sims_counter = 1, No_of_sims   

    ! do j = 0, 18
    
    !every simulation, reset the values of the parameters
    seed = 48
    occupied = 0
    failed_counter = 0
    ti = 0.0_8
    tj = 0.0_8
    ti_matter = 0.0_8
    tj_matter = 0.0_8
    clock = 0.0_8
    n = 1
    H = 0.06
    k = 0.80
    min_init_distance_id(:) = 0
    minimum_distance = 0
    j = 1 


    call pcg_init(seed, seed)
    call initialize_universe(civil_xyz, N_pt, init_distance)  !distribute the planets in the universe

    !generate the elements of the random_number_array. Don't stop until count(random_number_array >= k) == 5000
        
    do while(rn_checker < N_pt)
        r = pcg_random_real()
        rna_size = rna_size + 1
        random_number_array(rna_size) = r 
        rn_checker = count(random_number_array > k)
    end do

    print*, "Generated random numbers info : ", (count(random_number_array >= k)), rna_size

    ! Locate starting civilization closest to center
    min_init_distance_id = minloc(init_distance(1:N_pt))
    minimum_distance = sqrt(init_distance(min_init_distance_id(1)))

    civil_xyz(min_init_distance_id(1))%status = 1
    occupied = 1 

    ! Main simulation loop
    do while (occupied < N_pt) !this simulation will not stop until all the points are occupied
        ! allocate(random_number_array(rna_size))
        ! print*, "Successfull generated random number array. Size = ", rna_size, count(random_number_array >= k)

        !we will ony assign targets if there are still unoccupied planets. 
        ! Step 1: Scan for active colonies (status 1) and assign targets
        possible_targets = count(civil_xyz(:)%status == 0)
        do i = 1, N_pt
            if (civil_xyz(i)%status == 1 .and. possible_targets > 0) then 
                call find_neighbor(i, N_pt, civil_xyz, min_dist_id, min_distance)
                ! Inside find_neighbor or right after calling it:
                if (min_distance >= huge(1.0_8) / 2.0_8) then
                    ! Reset origin history to allow retries if trapped
                    civil_xyz(min_dist_id(1))%origin_history = 0
                end if

                ! Ensure a valid unoccupied neighbor was found
                if (min_distance < huge(1.0_8) / 2.0_8) then
                    civil_xyz(i)%target = min_dist_id(1)
                    civil_xyz(i)%status = 2  ! Now travelling
                    civil_xyz(min_dist_id(1))%status = 3  ! Target reserved

                    civil_xyz(i)%total_distance = civil_xyz(i)%total_distance + min_distance
                    civil_xyz(i)%remaining_distance = min_distance
                    minimum_distances(i) = min_distance   !stores the distance to the array of minimum distances
                end if
            else if (civil_xyz(i)%status == 1 .and. possible_targets == 0) then 
                !do nothing. All the occupied planets do not have any possible target.
            else if (civil_xyz(i)%status == 2) then 
                minimum_distances(i) = civil_xyz(i)%remaining_distance
            end if 

        end do !this loop will stop only when all the occupied points find its target 
        
        ! Step 2: Determine global minimum remaining distance step
        l_min_distances_id = minloc(minimum_distances(1:N_pt))
        l_min_distance = minimum_distances(l_min_distances_id(1))

        ! Guard against infinite loop if no targets remain
        if (l_min_distance >= huge(1.0_8) / 2.0_8) then
            print*, "Loop ended because of this line"
        exit
        end if 

        ! Step 3: Advance global simulation clock
        clock = clock + l_min_distance
        tj = (-1/H)*log(exp(-H*ti) - H*l_min_distance) !for dark energy
        tj_matter=  (2.0/(3.0*H))*((1 + (3.0/2.0)*H*ti_matter)**(1.0/3.0) + H*l_min_distance/2.0)**3.0 - (2.0/(3.0*H))

        cpath = cpath + l_min_distance
        ! Step 4: Move all travelling ships forward by l_min_distance
        do i = 1, N_pt  !this loop will stop when all planets of status 2 moved a distance of l_min_distance 
            if (civil_xyz(i)%status == 2) then 
                civil_xyz(i)%remaining_distance = civil_xyz(i)%remaining_distance - l_min_distance
                civil_xyz(i)%total_distance = civil_xyz(i)%total_distance + l_min_distance !we add the total distance regardless of outcome

                if (civil_xyz(i)%remaining_distance <= 1.0e-10_8) then !if it reaches its target 
                    random_number = random_number_array(rna_counter)
                    rna_counter = rna_counter + 1 
                    civil_xyz(i)%status = 1
                    minimum_distances(i) = huge(1.0_8)  ! Reset distance tracker for completed trip
                    
                    
                    if(random_number >= k) then   !the civilisation survives
                        civil_xyz(civil_xyz(i)%target)%status = 1 
                        occupied = occupied + 1 
                        

                    else  !the civilisation does not survive
                        civil_xyz(civil_xyz(i)%target)%status = 0
                        failed_counter = failed_counter + 1
                        civil_xyz(i)%counter = civil_xyz(i)%counter + 1 
                        ! print*, "number of possible targets : ", possible_targets, "Generated num : ", random_number

                        !the target planet will remember the failed occupations and prevent the origin planet to occupy it for 10 more unsuccessful attempts 
                        if(mod(civil_xyz(i)%counter, 1) == 0) then
                            civil_xyz(civil_xyz(i)%target)%origin_history(1) = i
                        else
                            civil_xyz(civil_xyz(i)%target)%origin_history(mod(civil_xyz(i)%counter, 1)) = i
                        end if 
                    
                    end if 
                    civil_xyz(i)%target = 0 !reset the target 
                    
                    ! print*, "Occupied = ", occupied, "failed = ", failed_counter, 'Static Time = ', clock, &
                    !         'Dark energy dominated = ', tj, "Matter dominated energy =", tj_matter

                else
                    minimum_distances(i) = civil_xyz(i)%remaining_distance
                end if 
            end if 
        end do !this loop will stop when all planets of status 1 moved a distance of l_min_distance 
      
        ti = tj 
        ti_matter = tj_matter
         
    end do !this loop will stop when all the planets are occupied 
    
    
    write(10,*) k, clock
    print*, k, clock

! end do 

end do !this loop will stop after 500 simulations 

    call cpu_time(end_time)
    
    print*, "DONE"
    print*, "k = ", k, "Success occupation: ",occupied ,'Failed occupation', failed_counter
    print*, "Number of unoccupied sites: ", count(civil_xyz(:)%status == 0 )
    print*, "Number of occupied sites: ", count(civil_xyz(:)%status == 1 )
    print*, "Number of travelling: ", count(civil_xyz(:)%status == 2 )
    print*, "Number of target: ", count(civil_xyz(:)%status == 3 )
    print*, 'Time elapsed: ', clock

    close(10)
end program sample



