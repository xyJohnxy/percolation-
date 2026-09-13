
program sample
    use pcg_module 
    use civilization_mod
    implicit none 

    real :: start_time, end_time 
    integer :: i, n
    integer, parameter :: N_pt = 5000, No_of_sims = 1
    real(8), parameter :: k_fixed = 0.80, H_fixed = 0.01
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
    real(8) :: ti_matter, tj_matter
    real(8) :: k, random_number !death rate, random_number
    integer :: failed_counter
    integer(8) :: no_of_sims_counter
    integer(8) :: seed

    

    call cpu_time(start_time)
    occupied = 0
    failed_counter = 0
    ti = 0.0_8
    tj = 0.0_8
    ti_matter = 0.0_8
    tj_matter = 0.0_8
    clock = 0.0_8
    n = 1
    H = 0.02
    seed = 23
    k = 0.90_8


    print*, "Simulating..."
    ! Allocate array sizes 
    allocate(civil_xyz(N_pt))
    allocate(dist_n(N_pt))
    allocate(minimum_distances(N_pt))

    

    ! open(unit = 10, file = 'initial_points.dat', action = 'write', status = 'replace')
    ! open(unit = 11, file = 'n(t)_vs_t_k=0.6.dat', action = 'write', status = 'replace')
    ! open(unit = 12, file = 'n(t)_vs_t_k=0.8,h=0.1.dat', action = 'write', status = 'replace')
    ! open(unit = 13, file = "H vs time_completed.dat", action = "write", status = "replace")
    ! open(unit = 14, file = "500sims_average n(t) vs time", action = "write", status = "replace")
    open(unit = 15, file = 'sims=1,k=vary2.dat', action ='write',status='replace')

    ! print*, "Occupied = ", occupied, "failed = ", failed_counter, 'Static Time = ', clock, 'Dark energy dominated = ', tj
    ! write(11,*) occupied, clock, tj
    ! write(12,*) occupied, clock, tj, tj_matter
    ! Generate random civilizations



    print*, "Simulation 1"

    ! do while(k <= 0.90)   !this loop runs the 500 simulations 
    do no_of_sims_counter = 1, No_of_sims 
        call pcg_init(seed, seed)
        minimum_distances = huge(1.0_8) 

        !initialize the parameters per simulation 
        occupied = 0
        failed_counter = 0
        ti = 0.0_8
        tj = 0.0_8
        ti_matter = 0.0_8
        tj_matter = 0.0_8
        clock = 0.0_8
        n = 1

        ! call initialize_universe(civil_xyz, N_pt, init_distance)
  
    do while (n <= N_pt)  !this loop generates the planets
        civil_xyz(n)%rx = 10.0_8 * pcg_random_real() - 5.0_8
        civil_xyz(n)%ry = 10.0_8 * pcg_random_real() - 5.0_8
        civil_xyz(n)%rz = 10.0_8 * pcg_random_real() - 5.0_8

        !initialize the planets objects 
        civil_xyz(n)%status = 0
        civil_xyz(n)%target = 0
        civil_xyz(n)%remaining_distance = 0 
        civil_xyz(n)%total_distance = 0
        civil_xyz(n)%counter = 0 
        civil_xyz(n)%origin_history = 0 

        init_distance(n) = civil_xyz(n)%rx**2 + civil_xyz(n)%ry**2 + civil_xyz(n)%rz**2

        if (init_distance(n) > radius**2) cycle
        n = n + 1 
    end do !this loop will terminate if all the planets have been generated 

    ! Locate starting civilization closest to center
    min_init_distance_id = minloc(init_distance(1:N_pt))
    minimum_distance = sqrt(init_distance(min_init_distance_id(1)))

    civil_xyz(min_init_distance_id(1))%status = 1
    occupied = 1 

    ! print*, "Occupied = ", occupied, "failed = ", failed_counter, 'Static Time = ', clock, 'Dark energy dominated = ', tj
    ! write(11,*) occupied, clock, tj
    ! write(12,*) occupied,clock, tj, tj_matter
    ! Main simulation loop

    do while (occupied < N_pt) !this loop will not stop untill all the planets are occupied
        ! Step 1: Scan for active colonies (status 1) and assign targets
        do i = 1, N_pt
            if (civil_xyz(i)%status == 1) then 
                
                call find_neighbor(i, N_pt, civil_xyz, min_dist_id, min_distance)
                ! Ensure a valid unoccupied neighbor was found
                if (min_distance < huge(1.0_8) / 2.0_8) then
                    civil_xyz(i)%target = min_dist_id(1)
                    civil_xyz(i)%status = 2  ! Now travelling
                    civil_xyz(min_dist_id(1))%status = 3  ! Target reserved

                    civil_xyz(i)%total_distance = civil_xyz(i)%total_distance + min_distance
                    civil_xyz(i)%remaining_distance = min_distance
                    minimum_distances(i) = min_distance   !stores the distance to the array of minimum distances
                end if
            else if (civil_xyz(i)%status == 2) then 
                minimum_distances(i) = civil_xyz(i)%remaining_distance
            end if 
        end do   !this loop will terminate if all the occupied points find its target

        
        ! Step 2: Determine global minimum remaining distance step
        l_min_distances_id = minloc(minimum_distances(1:N_pt))
        l_min_distance = minimum_distances(l_min_distances_id(1))

        ! Guard against infinite loop if no targets remain
        if (l_min_distance >= huge(1.0_8) / 2.0_8) exit

        ! Step 3: Advance global simulation clock
        clock = clock + l_min_distance
        tj = (-1/H)*log(exp(-H*ti) - H*l_min_distance) !for dark energy
        tj_matter=  (2.0/(3.0*H))*((1 + (3.0/2.0)*H*ti_matter)**(1.0/3.0) + H*l_min_distance/2.0)**3.0 - (2.0/(3.0*H))

        ! Step 4: Move all travelling planets forward by l_min_distance
        do i = 1, N_pt 
            if (civil_xyz(i)%status == 2) then 
                civil_xyz(i)%remaining_distance = civil_xyz(i)%remaining_distance - l_min_distance

                !if it reaches its target 
                if (civil_xyz(i)%remaining_distance <= 1.0e-10_8) then 
                    random_number = pcg_random_real()
                    civil_xyz(i)%status = 1  !the origin planet automatically resets to 0 regardless of the outcome 
                    minimum_distances(i) = huge(1.0_8)  ! Reset distance tracker for completed trip

                    ! stores the planets last 5 ventures
                    civil_xyz(i)%counter = mod(civil_xyz(i)%counter, 30) + 1
                    civil_xyz(i)%origin_history(civil_xyz(i)%counter) = civil_xyz(i)%target
                    
                    if(random_number > k) then   !the civilisation survives
                        civil_xyz(civil_xyz(i)%target)%status = 1 
                        occupied = occupied + 1 

                    else  !the civilisation does not survive
                        civil_xyz(civil_xyz(i)%target)%status = 0
                        failed_counter = failed_counter + 1
                    
                    end if 
                    ! print*, "Occupied = ", occupied, "failed = ", failed_counter, 'Static Time = ', clock, &
                    !         'Dark energy dominated = ', tj, "Matter dominated energy =", tj_matter

                    civil_xyz(i)%target = 0 !resets the target to 0 

                else
                    minimum_distances(i) = civil_xyz(i)%remaining_distance
                end if 
            end if 

        ! print*, "Planet ID: starting, histories: ", civil_xyz(min_init_distance_id(1))%origin_history, "Status: ", civil_xyz(min_init_distance_id(1))%status
        end do   !this loop will terminate if all the planets move by l_min_distance 

        ! write(11,*) clock, occupied
        ! print*, "Number of occupied sites = ", occupied , " | Clock = ", 
        ! write(11,*) occupied, clock, tj
        ! write(12,*) occupied, clock, tj, tj_matter
        ti = tj 
        ti_matter = tj_matter
        print*,occupied, clock
    end do    !this loop will terminate if all the points are occupied

    ! close(11)
    ! close(12)
    ! print*, "K = ", k , "Occupation time : ", clock
    ! write(15,*)k ,clock

    ! print*, "occupied: " , occupied, "| End time (s, D.E., M.D.): ", clock, tj, tj_matter 
    ! print*, H, tj, tj_matter
    ! write(13, *) H, tj, tj_matter
    ! write(12,*) occupied, clock, tj, tj_matter

    ! H = H + 0.005
    
! end do
end do !this loop will terminate if there have been 500 simulations
! k = k+0.05_08
 print*, "(10 Last 5 targets) : ", civil_xyz(10)%origin_history
    print*, "Number of unoccupied sites: ", count(civil_xyz(:)%status == 0 )
    print*, "Number of occupied sites: ", count(civil_xyz(:)%status == 1 )
    print*, "Number of travelling: ", count(civil_xyz(:)%status == 2 )
    print*, "Number of target: ", count(civil_xyz(:)%status == 3 )


    print*, "K = ", k , "Occupation time : ", clock
! end do 
   
    call cpu_time(end_time)
    print*, 'Time elapsed: ', end_time - start_time
    ! close(13)
    close(15)


    

end program sample


