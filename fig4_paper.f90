!this code replicates figure 4 in the second paper of dr alinea. I want to show that the occupation time to reach the residual planets at the outer planets contribute signitificantly
!to the total occupation time 


program fig4
implicit none 

integer :: io_stat
real :: col1,col2, last_col2, col3, col4

! open(unit = 10, file = "original,average_n_vs_t,sims=500,k=0.0,H=0.02.dat" , status = "old", action = "read" )
! open(unit = 11, file = "fig4_data.dat" , status = "replace", action = "write" )

open(unit = 10, file = "enclosed_average_n_vs_t,sims=500,k=0.60,H=0.02.dat" , status = "old", action = "read" )
open(unit = 11, file = "fig4_enclosed_data.dat" , status = "replace", action = "write" )

do 
    read(10, *, iostat = io_stat) col1, col2
    if (io_stat /= 0) exit 

    last_col2 = col2 
end do 

rewind(10)

do 
    read(10, *, iostat = io_stat) col1, col2
    if (io_stat /= 0) exit 

    col3 = col1/5000.0 * 100.0
    col4 = col2/last_col2 * 100.0 
    print*, col3, col4

    write(11,*) col3, col4
end do 
print*, "done"
close(10)
close(11)
end program fig4