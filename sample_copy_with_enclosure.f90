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
        logical, dimension(8) :: surrounded = .false.
    end type civil 

contains

    function distance(i, j, civil_xyz) result(dist)
        integer, intent(in) :: i,j
        type(civil), intent(in) :: civil_xyz(:)

        real(8):: dist

        dist = (civil_xyz(i)%rx - civil_xyz(j)%rx)**2 + &
                                    (civil_xyz(i)%ry - civil_xyz(j)%ry)**2 + &
                                    (civil_xyz(i)%rz - civil_xyz(j)%rz)**2
    end function distance

    subroutine find_neighbor(i, N_pt, civil_xyz, minimum_distance_id, minimum_distance)
        integer, intent(in) :: i, N_pt
        type(civil), intent(inout) :: civil_xyz(:)
        
        ! Outputs must use intent(out)
        integer, intent(out) :: minimum_distance_id(1) 
        real(8), intent(out) :: minimum_distance

        integer :: j
        real(8) :: dist_min, dist_dummy, max_range 

        max_range = 0.9
        ! Use huge() so non-matching sites are never picked by minloc
        dist_min = huge(1.0_8)
        minimum_distance_id(1) =0 !-> default value 
        
        
        do j = 1, N_pt
            if (j /= i .and. civil_xyz(j)%status == 0) then 
                dist_dummy = distance(i,j,civil_xyz)

                if(dist_dummy < dist_min) then
                    dist_min = dist_dummy 
                    minimum_distance_id(1) = j
                end if 
            end if 
        end do 

        ! minimum_distance_id = minloc(dist_n(1:N_pt))
        minimum_distance = sqrt(dist_min)

        if(minimum_distance > max_range) then 
            civil_xyz(i)%status = 4 !retired
            minimum_distance = huge(1.0_8)
        end if 
            
    end subroutine find_neighbor


end module civilization_mod

!------------------------------------------------------------------------------------------------------------

program sample
    use pcg_module 
    use civilization_mod
    implicit none 

    real :: start_time, end_time 
    integer :: i, n
    integer, parameter :: N_pt = 5000, No_of_sims = 1
    real(8), parameter :: k = 0.0, H_fixed = 0.001

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
    real(8) :: random_number !death rate, random_number
    integer :: failed_counter, retired 
    integer(8) :: seed = 48, no_of_sims_counter = 1 

    !statistics variables 
    real(8), dimension(No_of_sims) :: occupation_t_static, occupation_t_de, occupation_t_m
    real(8) :: mean_time_s = 0, std_static = 0
    real(8) :: mean_time_de = 0, std_de = 0
    real(8) :: mean_time_m = 0, std_m = 0

    integer :: totally_occupied = 0 ! checks for the number of simulation that successfully occupied all the planets

    call cpu_time(start_time)
    occupied = 0
    retired = 0
    failed_counter = 0
    ti = 0.0_8
    tj = 0.0_8
    ti_matter = 0.0_8
    tj_matter = 0.0_8
    clock = 0.0_8
    n = 1
    H = 0.04

    print*, "working"
    ! Allocate array sizes 
    allocate(civil_xyz(N_pt))
    allocate(dist_n(N_pt))
    allocate(minimum_distances(N_pt))

    minimum_distances = huge(1.0_8) 

    ! open(unit = 10 , file = "sims=10,k=vary, type=octant_enclosed.dat", position = "APPEND", action = "write")
    ! open(unit = 20 , file = "3D_with_enclosure_static,plot,k=vary,sims=500,hb=0.6.dat", action = "write", status = "replace")
    ! open(unit = 10 , file = "3D_with_enclosure_static,curve,k=vary,sims=500,hb=3.0.dat", position = "APPEND", action = "write")
    open(unit = 11 , file = "n(t)_vs_t,sims=1,k=0.0,max_range=0.9.dat", status = "replace" , action = "write")
    ! open(unit = 12 , file = "n(t)_vs_t,dark_energy,max_range=0.9,sims=1,H=0.12.dat", status = "replace" , action = "write")

do no_of_sims_counter = 1, No_of_sims 
        ! 1. Re-initialize scalars
        occupied = 0
        retired = 0
        failed_counter = 0
        ti = 0.0_8
        tj = 0.0_8
        ti_matter = 0.0_8
        tj_matter = 0.0_8
        clock = 0.0_8
        n = 1
        H = 0.02

        ! 2. Reset array tracking (CRITICAL FIX)
        minimum_distances = huge(1.0_8) 

        call pcg_init(seed, no_of_sims_counter)

        ! 3. Generate points and reset ALL struct fields (CRITICAL FIX)
        do while (n <= N_pt)
            civil_xyz(n)%rx = 10.0_8 * pcg_random_real() - 5.0_8
            civil_xyz(n)%ry = 10.0_8 * pcg_random_real() - 5.0_8
            civil_xyz(n)%rz = 10.0_8 * pcg_random_real() - 5.0_8

            civil_xyz(n)%status = 0
            civil_xyz(n)%target = 0
            civil_xyz(n)%total_distance = 0.0_8
            civil_xyz(n)%remaining_distance = 0.0_8

            init_distance(n) = civil_xyz(n)%rx**2 + civil_xyz(n)%ry**2 + civil_xyz(n)%rz**2

            if (init_distance(n) > radius**2) cycle
            n = n + 1 
        end do
    

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
                !after scanning for occupied planet, we have to check first if it is already surrounded. 
                ! If it is, update the status to 4 (retired) and move on to the next planet

                call find_neighbor(i, N_pt, civil_xyz, min_dist_id, min_distance)
                if(civil_xyz(i)%status == 4) then
                    cycle
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

                else
                    minimum_distances(i) = civil_xyz(i)%remaining_distance
                end if 
            end if 
        end do

        write(11,*) occupied, clock, tj,tj_matter 
        
        if(mod(occupied,100)==0) then 
        print*, occupied, clock, tj, tj_matter
        end if 

        ti = tj 
        ti_matter = tj_matter
    end do !-> this loop will end when all the planets are occupied 

    if(occupied == N_pt) then 
        totally_occupied = totally_occupied+1
    end if 
    
    occupation_t_static(no_of_sims_counter) = clock
    occupation_t_de(no_of_sims_counter) = tj 
    occupation_t_m(no_of_sims_counter) = tj_matter
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
    else if(no_of_sims_counter == 400) then 
        print*, "Done 400 sims"
    else if(no_of_sims_counter == 500) then 
        print*, "Done 500 sims"
    end if 

end do !-> this loop will terminate after 500 simulations 

!this loop is to visually see how the planets occupy their neighbors
    ! do i = 1, N_pt 
    !     write(20,*) civil_xyz(i)%rx, civil_xyz(i)%ry, civil_xyz(i)%rz, civil_xyz(i)%status
    ! end do 



    !calculate the mean and std
    !mean
    print*, "Calculating the statistics "
    do i=1,No_of_sims
        mean_time_s= mean_time_s + occupation_t_static(i)
        mean_time_de= mean_time_de + occupation_t_de(i)
        mean_time_m= mean_time_m + occupation_t_m(i)
    end do 
    mean_time_s = mean_time_s/No_of_sims !-> mean time after 500 simulations 
    mean_time_de = mean_time_de/No_of_sims
    mean_time_m = mean_time_m/No_of_sims
    
    !std 
    do i = 1,No_of_sims
        std_static = std_static + (occupation_t_static(i)-mean_time_s)**2 
        std_de = std_de + (occupation_t_de(i)-mean_time_de)**2 
        std_m = std_m + (occupation_t_m(i)-mean_time_m)**2 
    end do 

    std_static = sqrt(std_static/N_pt)
    std_de= sqrt(std_de/N_pt)
    std_m = sqrt(std_m/N_pt)

    print*, "Number of simulations: ", No_of_sims
    print*, "mortality factor: ", k
    print*, "Static: mean: ", mean_time_s, "std: ", std_static
    print*, "Dark energy: mean: ", mean_time_de, "std: ", std_de
    print*, "Matter: mean: ", mean_time_m, "std: ", std_m

    ! write(10, *) k, mean_time_s, mean_time_de, mean_time_m
!-------------------------------------------------------------------------------------------------------

    print*, "Results : ", occupied,clock, tj, tj_matter
    call cpu_time(end_time)
    print*, 'Time elapsed: ', end_time - start_time
    print*, "Number of simulations that did not occupy all the planets", No_of_sims - totally_occupied

    print*, "Number of unoccupied sites: ", count(civil_xyz(:)%status == 0 )
    print*, "Number of occupied sites: ", count(civil_xyz(:)%status == 1 )
    print*, "Number of travelling: ", count(civil_xyz(:)%status == 2 )
    print*, "Number of target: ", count(civil_xyz(:)%status == 3 )
    print*, "Number of retired: ", count(civil_xyz(:)%status == 4 )
    close(10)
    ! close(20)
    close(12)
    close(11)


    

end program sample



