
module civilization_mod
    use pcg_module
    implicit none

    ! Global parameters
    real(8), parameter :: radius = 5.0_8
    real(8), parameter :: diameter = 10.0_8

    ! Derived type definition accessible to all routines using this module
    type :: civil 
        real(8) :: rx, ry
        integer :: status             ! 0 (unoccupied), 1 (active colony), 2 (travelling), 3 (targeted)
        integer :: target             ! index of the destination 
        real(8) :: total_distance = 0.0_8
        real(8) :: remaining_distance
        integer, dimension(5000) :: neighbors 
        integer :: neighbors_counter 
        
    end type civil 

contains

    function distance(i, j, civil_xyz) result(dist)
        integer, intent(in) :: i,j
        type(civil), intent(in) :: civil_xyz(:)

        real(8):: dist

        dist = (civil_xyz(i)%rx - civil_xyz(j)%rx)**2 + &
                                    (civil_xyz(i)%ry - civil_xyz(j)%ry)**2 
    end function distance

    subroutine find_neighbor(i, N_pt, civil_xyz, minimum_distance_id, minimum_distance)
        integer, intent(in) :: i, N_pt
        type(civil), intent(inout) :: civil_xyz(:)
        
        ! Outputs must use intent(out)
        integer, intent(out) :: minimum_distance_id(1) 
        real(8), intent(out) :: minimum_distance

        integer :: j
        real(8) :: dist_min, dist_dummy
        real(8) :: max_range = huge(1.0_8)
        

        ! Use huge() so non-matching sites are never picked by minloc
        dist_min = huge(1.0_8)
        
        
        do j = 1, N_pt
            if (j /= i .and. civil_xyz(j)%status == 0) then 
                dist_dummy = (civil_xyz(i)%rx - civil_xyz(j)%rx)**2 + &
                                 (civil_xyz(i)%ry - civil_xyz(j)%ry)**2 
                                 
                if(dist_dummy < dist_min) then
                    dist_min = dist_dummy 
                    minimum_distance_id(1) = j
                end if 
            end if 
        end do 

        ! minimum_distance_id = minloc(dist_n(1:N_pt))
        minimum_distance = sqrt(dist_min)
        if(minimum_distance > max_range) then 
            civil_xyz(i)%status = 4 !we convert it to retired if the nearest neighbor is outside its maximum range" 
            minimum_distance = huge(1.0_8)
        end if 

    end subroutine find_neighbor

    function area_of_triangle(a,b,c,civil_xyz) result(area)
        integer, intent(in) :: a,b,c 
        type(civil),intent(inout) :: civil_xyz(:)
        real(8) :: area 

        real(8) :: x1, x2, x3, y1,y2,y3 

        x1 = civil_xyz(a)%rx
        x2 = civil_xyz(b)%rx 
        x3 = civil_xyz(c)%rx 

        y1 = civil_xyz(a)%ry
        y2 = civil_xyz(b)%ry
        y3 = civil_xyz(c)%ry

        area = 0.5*abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))

    end function area_of_triangle

    subroutine triangulation(i, civil_xyz, retired)
        integer, intent(in) :: i
        type(civil), intent(inout) :: civil_xyz(:)
        integer, intent(inout) :: retired

        !local variables
        integer :: q,r,s
        integer :: q_id, r_id, s_id 
        real(8) :: A, A1,A2,A3

        do q = 1, civil_xyz(i)%neighbors_counter
            q_id =  civil_xyz(i)%neighbors(q)
            if(civil_xyz(q_id)%status == 0 .or. civil_xyz(q_id)%status == 3 ) cycle 
            if(civil_xyz(i)%status  == 4 ) exit

            do r = q+1,civil_xyz(i)%neighbors_counter
                r_id =  civil_xyz(i)%neighbors(r)
                if(civil_xyz(r_id)%status == 0 .or. civil_xyz(r_id)%status == 3 ) cycle 
                if(civil_xyz(i)%status  == 4 ) exit 

                do s = r+1,civil_xyz(i)%neighbors_counter
                s_id =  civil_xyz(i)%neighbors(s)
                if(civil_xyz(s_id)%status == 0 .or. civil_xyz(s_id)%status == 3 ) cycle

                A = area_of_triangle(q_id,r_id,s_id,civil_xyz)
                A1 = area_of_triangle(i,r_id,s_id,civil_xyz)
                A2 = area_of_triangle(i,q_id,s_id,civil_xyz)
                A3 = area_of_triangle(i,r_id,q_id,civil_xyz)

                if(abs(A - (A1+A2+A3)) <= 0.0_8) then
                    !enclosed
                    civil_xyz(i)%status = 4 !enclosed na sya 
                    retired = retired + 1 
                    print*, "Successful triangulation of planet ", i
                    exit
                end if 

                end do 
            end do 
        end do 
    end subroutine triangulation 


end module civilization_mod

!------------------------------------------------------------------------------------------------------------

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
    real(8) :: cpath, t_total, hubble_sphere, dum_dist
    integer(8) :: seed
    integer :: retired


    call cpu_time(start_time)
    seed = 48
    occupied = 0
    failed_counter = 0
    ti = 0.0_8
    tj = 0.0_8
    ti_matter = 0.0_8
    tj_matter = 0.0_8
    clock = 0.0_8
    hubble_sphere = 1.0


    n = 1
    H = 0.06
    k = 0.0
    cpath = 0
    t_total = 0  


    print*, "working"
    ! Allocate array sizes 
    allocate(civil_xyz(N_pt))
    allocate(dist_n(N_pt))
    allocate(minimum_distances(N_pt))

    minimum_distances = huge(1.0_8) 

    open(unit = 20 , file = "2D_status_checker_triangulation,N=vary,N_pt=5000.dat", action = "write", status = "replace")
    
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
    retired = 0
    ! H = H + 0.005
    ! k = k + 0.1

    do while (n <= N_pt)
        civil_xyz(n)%rx = 10.0_8 * pcg_random_real() - 5.0_8
        civil_xyz(n)%ry = 10.0_8 * pcg_random_real() - 5.0_8

        civil_xyz(n)%status = 0
        civil_xyz(n)%neighbors = 0
        civil_xyz(n)%neighbors_counter = 0

        init_distance(n) = civil_xyz(n)%rx**2 + civil_xyz(n)%ry**2

        if (init_distance(n) > radius**2) cycle
        n = n + 1 
    end do !this loop will end when there have been 5000 generated points 
    
    !finding neighbors to a planet
    do i=1,N_pt
        do j=1, N_pt 
            if(i==j) cycle 
            dum_dist = distance(i,j,civil_xyz)
            ! print*,dum_dist
            if( dum_dist<= hubble_sphere**2) then 
                civil_xyz(i)%neighbors_counter = civil_xyz(i)%neighbors_counter + 1
                civil_xyz(i)%neighbors(civil_xyz(i)%neighbors_counter) = j
                ! print*, "Planet ", i , "found a neighbor"
        end if 
    end do 
    end do 


    ! Locate starting civilization closest to center
    min_init_distance_id = minloc(init_distance(1:N_pt))
    minimum_distance = sqrt(init_distance(min_init_distance_id(1)))

    civil_xyz(min_init_distance_id(1))%status = 1
    occupied = 1 
    
    ! Main simulation loop
    do while (occupied < 5000)
        
        ! Step 1: Scan for active colonies (status 1) and assign targets
        do i = 1, N_pt
            if (civil_xyz(i)%status == 1) then 
                ! call occupied_neighbor(i, N_pt, civil_xyz, retired)  !we update the civil%surrounded array
                
                call triangulation(i,civil_xyz, retired)
                if(civil_xyz(i)%status == 4) then
                    cycle
                end if 

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

    !this loop is to visually see how the planets occupy their neighbors
    do i = 1, N_pt 
        write(20,*) civil_xyz(i)%rx, civil_xyz(i)%ry, civil_xyz(i)%status
    end do 

    
    ! print*, "k: ", k , "Time elapsed: ", clock
    
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
    

    call cpu_time(end_time)
     
    print*, "Number of unoccupied sites: ", count(civil_xyz(:)%status == 0 )
    print*, "Number of occupied sites: ", count(civil_xyz(:)%status == 1 )
    print*, "Number of travelling: ", count(civil_xyz(:)%status == 2 )
    print*, "Number of target: ", count(civil_xyz(:)%status == 3 )
    print*, "Number of retired: ", count(civil_xyz(:)%status == 4 )
    print*, "Failed occupation: ", failed_counter

    
    ! close(13)
    ! close(11)
    ! close(10)
    close(20)



    

end program sample




