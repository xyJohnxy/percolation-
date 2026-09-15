module civilization_mod
    use pcg_module
    implicit none

    ! Global parameters
    real(8), parameter :: radius = 5.0_8
    real(8), parameter :: diameter = 10.0_8

    ! Derived type definition accessible to all routines using this module
    type :: civil 
        real(8) :: rx, ry, rz
        integer :: status             ! 0 (unoccupied), 1 (active colony), 2 (travelling), 3 (targeted)
        integer :: target             ! index of the destination 
        real(8) :: total_distance = 0.0_8
        real(8) :: remaining_distance
    end type civil 

contains

    subroutine find_neighbor(i, N_pt, civil_xyz, minimum_distance_id, minimum_distance)
        integer, intent(in) :: i, N_pt
        type(civil), intent(in) :: civil_xyz(:)
        
        ! Outputs must use intent(out)
        integer, intent(out) :: minimum_distance_id(1) 
        real(8), intent(out) :: minimum_distance

        integer :: j
        real(8) :: dist_min, dist_dummy

        ! Use huge() so non-matching sites are never picked by minloc
        dist_min = huge(1.0_8)
        
        
        do j = 1, N_pt
            if (j /= i .and. civil_xyz(j)%status == 0) then 
                dist_dummy = (civil_xyz(i)%rx - civil_xyz(j)%rx)**2 + &
                                 (civil_xyz(i)%ry - civil_xyz(j)%ry)**2 + &
                                 (civil_xyz(i)%rz - civil_xyz(j)%rz)**2

                if(dist_dummy < dist_min) then
                    dist_min = dist_dummy 
                    minimum_distance_id(1) = j
                end if 
            end if 
        end do 

        ! minimum_distance_id = minloc(dist_n(1:N_pt))
        minimum_distance = sqrt(dist_min)

    end subroutine find_neighbor

end module civilization_mod

!------------------------------------------------------------------------------------------------------------

program sample
    use pcg_module 
    use civilization_mod
    implicit none 

    real :: start_time, end_time 
    integer :: i, n
    integer, parameter :: N_pt = 5000, No_of_sims = 500
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
    real(8) :: cpath, t_total
    integer(8) :: seed

    !statistics variables 
    real(8), dimension(No_of_sims) :: occupation_t_static 
    real(8) :: mean_time_s = 0, std_static = 0

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
    k = 0.75_8
    cpath = 0
    t_total = 0  


    print*, "working"
    ! Allocate array sizes 
    allocate(civil_xyz(N_pt))
    allocate(dist_n(N_pt))
    allocate(minimum_distances(N_pt))

    minimum_distances = huge(1.0_8) 

    open(unit = 10 , file = "sims=500,k=vary.dat", position = "APPEND", action = "write")
    
    do no_of_sims_counter = 1, No_of_sims 
    
    call pcg_init(seed,no_of_sims_counter)
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

    do while (n <= N_pt)
        civil_xyz(n)%rx = 10.0_8 * pcg_random_real() - 5.0_8
        civil_xyz(n)%ry = 10.0_8 * pcg_random_real() - 5.0_8
        civil_xyz(n)%rz = 10.0_8 * pcg_random_real() - 5.0_8

        civil_xyz(n)%status = 0

        init_distance(n) = civil_xyz(n)%rx**2 + civil_xyz(n)%ry**2 + civil_xyz(n)%rz**2

        if (init_distance(n) > radius**2) cycle
        n = n + 1 
    end do !this loop will end when there have been 5000 generated points 
    ! close(10)

    ! Locate starting civilization closest to center
    min_init_distance_id = minloc(init_distance(1:N_pt))
    minimum_distance = sqrt(init_distance(min_init_distance_id(1)))

    civil_xyz(min_init_distance_id(1))%status = 1
    occupied = 1 
    
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
        end do  !this loop will end when all status 1 planets find neighbors
        
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
        end do !this loop will end when all travelling planets moved 
        ti = tj 
        ti_matter = tj_matter
         
    end do !this loop will end when all the planets are occupied

    
    ! print*, "k: ", k , "Time elapsed: ", clock
    occupation_t_static(no_of_sims_counter) = clock
    if(no_of_sims_counter == 1) then  
        print*, "Done 1 sims"
    else if(no_of_sims_counter == 10) then 
        print*, "Done 10 sims"
    else if(no_of_sims_counter == 20) then 
        print*, "Done 20 sims"
    else if(no_of_sims_counter == 50) then 
        print*, "Done 50 sims"
    else if(no_of_sims_counter == 100) then 
        print*, "Done 100 sims"
    else if(no_of_sims_counter == 250) then 
        print*, "Done 250 sims"
    else if(no_of_sims_counter == 500) then 
        print*, "Done 500 sims"
    end if 

end do !this loop will end when all 500 simulations are finished 
! close(16)



    !calculate the mean and std
    !mean
    print*, "Calculating the statistics "
    do i=1,No_of_sims
        mean_time_s= mean_time_s + occupation_t_static(i)
    end do 
    mean_time_s = mean_time_s/No_of_sims !-> mean time of the after 5000 measurements 
    
    !std 
    do i = 1,No_of_sims
        std_static = std_static + (occupation_t_static(i)-mean_time_s)**2 
    end do 
    std_static = sqrt(std_static/N_pt)




    print*, "Number of simulations: ", No_of_sims
    print*, "mortality factor: ", k, "mean: ", mean_time_s, "std: ", std_static

    write(10, *) k, mean_time_s, std_static

    call cpu_time(end_time)
     
    ! print*, "Number of unoccupied sites: ", count(civil_xyz(:)%status == 0 )
    ! print*, "Number of occupied sites: ", count(civil_xyz(:)%status == 1 )
    ! print*, "Number of travelling: ", count(civil_xyz(:)%status == 2 )
    ! print*, "Number of target: ", count(civil_xyz(:)%status == 3 )
    ! print*, 'Time elapsed: ', clock
    ! close(13)
    ! close(11)



    

end program sample



