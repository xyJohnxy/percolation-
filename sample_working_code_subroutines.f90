module civilization_mod
    use pcg_module
    implicit none

    ! Global parameters
    real(8), parameter :: radius = 5.0_8
    real(8), parameter :: diameter = 10.0_8

    ! Derived type definition accessible to all routines using this module
    type :: civil 
        !initialize the values of the objects of civil
        real(8) :: rx, ry, rz = 0
        integer :: status = 0             ! 0 (unoccupied), 1 (active colony), 2 (travelling), 3 (targeted)
        integer :: target = 0             ! index of the destination 
        real(8) :: total_distance = 0.0_8  !total distance travelled by the planet including the failed attempts 
        real(8) :: remaining_distance = 0.0_8
        integer :: counter = 0   !counts the number of failed attempts 
        integer, dimension(1) :: origin_history 
        
    end type civil 

contains

    subroutine find_neighbor(i, N_pt, civil_xyz, minimum_distance_id, minimum_distance)
        integer, intent(in) :: i, N_pt
        type(civil), dimension(N_pt), intent(in) :: civil_xyz(:)
        
        ! Outputs must use intent(out)
        integer, intent(out) :: minimum_distance_id(1) !index of the nearest neighbor 
        real(8), intent(out) :: minimum_distance  !distance of the nearest neighbor

        integer :: j
        real(8) :: dist_min, dist_dummy

        minimum_distance_id(1) = 0 
        minimum_distance = huge(1.0_8)
        
        ! we initialize the minimum distance to a large number so that it would immediately update the time it finds a closest neighbor
        dist_min = huge(1.0_8)
        
        !i -> index of the occupied planet
        !j -> index of the neighbors 
        !scans each planets 
        do j = 1, N_pt
            !searches only planets of status = 0 and bypass search the last 5 planets 
            
            if (j /= i .and. civil_xyz(j)%status == 0) then
                if(any(civil_xyz(j)%origin_history == i)) cycle !-> if j is an element of the origin history, we do not search it.
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


    subroutine initialize_universe(planets,N_pt,init_distance)
        implicit none 
        !parameters
        integer, intent(in) :: N_pt
        type(civil), dimension(N_pt) :: planets
        real(8), dimension(N_pt) :: init_distance

        !local parameters
        integer :: n = 1 

        do while(n<= N_pt)
            planets(n)%rx = 10.0_8 * pcg_random_real() - 5.0_8
            planets(n)%ry = 10.0_8 * pcg_random_real() - 5.0_8
            planets(n)%rz = 10.0_8 * pcg_random_real() - 5.0_8

            !every planet created, we initialize its objects to 0
            planets(n)%status = 0
            planets(n)%target = 0
            planets(n)%remaining_distance = 0 
            planets(n)%total_distance = 0
            planets(n)%origin_history(:) =0
            planets(n)%counter = 0 

            init_distance(n) = planets(n)%rx**2 + planets(n)%ry**2 + planets(n)%rz**2

            if (init_distance(n) > radius**2) cycle
            n = n + 1 
        end do 

        print*, "Successfully generated the planets"


    end subroutine initialize_universe

    subroutine planet_status(planet,index, N_pt)
        !this subroutine prints out the status of the planet 
        integer, intent(in) :: N_pt, index
        type(civil), dimension(N_pt) :: planet
        
        print*, "Planet position : ", planet(index)%rx, planet(index)%ry, planet(index)%rz 
        print*, "Status: ",planet(index)%status,  "Target: ", planet(index)%target
        ! print*, "Last targets: ", planet(index)%origin_history

    end subroutine 

end module civilization_mod

!------------------------------------------------------------------------------------------------------------