
program sample
    use pcg_module 
    use civilization_mod
    implicit none 

    real :: start_time, end_time 
    integer :: i, n
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
    integer(8) :: seed

    call cpu_time(start_time)
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
    k = 0.0
    cpath = 0
    t_total = 0  
    std_t_statistic = 0

    print*, "working"
    ! Allocate array sizes 
    allocate(civil_xyz(N_pt))
    allocate(dist_n(N_pt))
    allocate(minimum_distances(N_pt))

    minimum_distances = huge(1.0_8) 

    open(unit = 10, file = 'sims=1,k=vary', action = 'write', status = 'replace')
    ! open(unit = 11, file = 'n(t)_vs_t_k=0.5.dat', action = 'write', status = 'replace')
    ! open(unit = 12, file = 'n(t)_vs_t_k=0.8,h=0.1.dat', action = 'write', status = 'replace')
    ! open(unit = 13, file = "H vs time_completed.dat", action = "write", status = "replace")
    ! open(unit = 14, file = "500sims_average n(t) vs time", action = "write", status = "replace")
    ! open(unit = 15, file = "n(t)_vs_t_vary_k=0.0.dat", action = "write", status = 'replace')
    ! open(unit = 16, file = "n(t)_vs_t_vary_k.dat", action = "write", status = 'replace')
    ! print*, "Occupied = ", occupied, "failed = ", failed_counter, 'Static Time = ', clock, 'Dark energy dominated = ', tj
    ! write(11,*) occupied, clock, tj
    ! write(12,*) occupied, clock, tj, tj_matter
    ! Generate random civilizations
    ! do while(H <= 0.08)


    do no_of_sims_counter = 1, No_of_sims 
    ! do while(k <= 0.8)
    call pcg_init(seed+no_of_sims_counter,seed+no_of_sims_counter)
    occupied = 0
    failed_counter = 0
    ti = 0.0_8
    tj = 0.0_8
    ti_matter = 0.0_8
    tj_matter = 0.0_8
    clock = 0.0_8
    n = 1
    ! H = H + 0.005
    ! k = k + 0.1

    call initialize_universe(civil_xyz, N_pt, init_distance)

    ! Locate starting civilization closest to center
    min_init_distance_id = minloc(init_distance(1:N_pt))
    minimum_distance = sqrt(init_distance(min_init_distance_id(1)))

    civil_xyz(min_init_distance_id(1))%status = 1
    occupied = 1 
    ! cpath = cpath + minimum_distance
    ! print*,"(#occupied, #status = 1):", occupied, count(civil_xyz(:)%status == 1)

    ! print*,"Starting distane : ", minimum_distance

    ! print*, "Occupied = ", occupied, "failed = ", failed_counter, 'Static Time = ', clock, 'Dark energy dominated = ', tj
    write(10,*) occupied, clock, tj, tj_matter
    ! write(12,*) occupied,clock, tj, tj_matter
    ! write(15,*) occupied, clock 
    
    ! Main simulation loop
    do while (occupied < N_pt)
        
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
        end do
        
        ! Step 2: Determine global minimum remaining distance step
        l_min_distances_id = minloc(minimum_distances(1:N_pt))
        l_min_distance = minimum_distances(l_min_distances_id(1))

        ! Guard against infinite loop if no targets remain
        if (l_min_distance >= huge(1.0_8) / 2.0_8) exit

        ! Step 3: Advance global simulation clock
        clock = clock + l_min_distance
        tj = (-1/H)*log(exp(-H*ti) - H*l_min_distance) !for dark energy
        tj_matter=  (2.0/(3.0*H))*((1 + (3.0/2.0)*H*ti_matter)**(1.0/3.0) + H*l_min_distance/2.0)**3.0 - (2.0/(3.0*H))

        cpath = cpath + l_min_distance
        ! Step 4: Move all travelling ships forward by l_min_distance
        do i = 1, N_pt 
            if (civil_xyz(i)%status == 2) then 
                civil_xyz(i)%remaining_distance = civil_xyz(i)%remaining_distance - l_min_distance

                if (civil_xyz(i)%remaining_distance <= 1.0e-10_8) then 
                    random_number = pcg_random_real()
                    civil_xyz(i)%status = 1
                    minimum_distances(i) = huge(1.0_8)  ! Reset distance tracker for completed trip
                    
                    if(random_number > k) then   !the civilisation survives
                        civil_xyz(civil_xyz(i)%target)%status = 1 
                        occupied = occupied + 1 
                        

                    else  !the civilisation does not survive
                        civil_xyz(civil_xyz(i)%target)%status = 0
                        failed_counter = failed_counter + 1
                    
                    end if 

                    
                    ! print*, "Occupied = ", occupied, "failed = ", failed_counter, 'Static Time = ', clock, &
                    !         'Dark energy dominated = ', tj, "Matter dominated energy =", tj_matter

                else
                    minimum_distances(i) = civil_xyz(i)%remaining_distance
                end if 
            end if 
        end do

        !let's print out the number of status 1 planets. It should be 2 (the origin and the target)
        
        ! print*,clock,"(#occupied, #status = 1):", occupied, count(civil_xyz(:)%status == 1)
        ! write(11,*) clock, occupied
        ! write(15,*) occupied, clock 
        write(10,*) occupied, clock, tj, tj_matter
        ! print*, "Number of occupied sites = ", occupied , " | Clock = ", clock
        ! write(11,*) occupied, clock, tj
        ! write(11,*) occupied, clock
        ! write(12,*) occupied, clock, tj, tj_matter
        ti = tj 
        ti_matter = tj_matter
         
    end do 

    ! close(11)
    ! close(12)
    ! close(15)

    ! print*, "Simulation Number: " , no_of_sims_counter, "| End time (s, D.E., M.D.): ", clock, tj, tj_matter 
    ! print*, "Clock : ", clock 
    ! print*, "Clock - Cpath : ", clock - cpath 
    ! write(13, *) H, tj, tj_matter
    ! write(12,*) occupied, clock, tj, tj_matter
    ! print*, no_of_sims_counter
    ! t_total = t_total + clock    !the total end times of all 500 simulations. Used for calculating the mean end times.
    ! t_terminate_static(no_of_sims_counter) = clock  !store the end time to an array 
    ! write(16,*) k, clock
    ! print*, k, clock

    !Reset values 
    print*, 'Time elapsed: ', clock
   
end do
! close(16)
    
    
! end do


! !calculate the statistics 
! mean_time = t_total/No_of_sims

! do i = 1,No_of_sims 
!     print*,i, t_terminate_static(i)
!     std_t_statistic = std_t_statistic + (t_terminate_static(i) - mean_time)**2
! end do 

! std_t_statistic = sqrt(std_t_statistic/(No_of_sims-1))


    
    call cpu_time(end_time)
    
    ! print*, "Average elapsed time (500 sims): ", mean_time 
!     print*, "Standard deviation (500 sims): ", std_t_statistic
!     print*, "Std/mean_time*100 : ", std_t_statistic/mean_time  
    print*, "Number of unoccupied sites: ", count(civil_xyz(:)%status == 0 )
    print*, "Number of occupied sites: ", count(civil_xyz(:)%status == 1 )
    print*, "Number of travelling: ", count(civil_xyz(:)%status == 2 )
    print*, "Number of target: ", count(civil_xyz(:)%status == 3 )
    print*, 'Time elapsed: ', clock
    ! close(13)
    ! close(11)



    

end program sample



