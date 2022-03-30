program isingmodel

! JS: NOTE:  There is a book on reserve in the science library,
! JS:  David Chandler, "Introduction to Modern Statistical Mechanics"
! JS:  that describes both the Ising model and the Monte Carlo procedure
! JS:  in some depth.  It might be a useful resource if you get stuck.

! array, defining variables

implicit none
real :: kT, E, E2, M, C, paccept,ptrial, Estart,Eflip

integer :: col, row, east,west,north,south
integer :: nMC
integer, parameter :: J = -1    ! define a ferromagnetic system
integer, dimension(16,16) :: spin 

! random number generator
call random_seed()

! initialized our spin array
do col=1,16
   do row=1,16
      call random_number(paccept)
      if (paccept .lt. 0.5) then
         spin(col,row) = +1
      else 
         spin(col,row) = -1
      end if
   end do
end do

! program structure:

! defining the loops (and what do we mean by loops)

! loop over temperature (kT in units of the coupling, J)

do kT=4.0, 0.1, -0.05      !JS: Start at kT=4, and go to kT=0.1 in steps of -0.05

! JS: First, run a set of MC steps to equilibrate our system. 
! JS: In other words, so that the configuration of the spins is
! JS: in an equilibrium configuration consistent with Boltzman's laws
! JS: for the particular temperature 

! loop MC trials
   do nMC=1,10000

! spatial loop in the 2D plane (incorporate pbc)

      do col=1,16
         do row=1,16

            east = col+1
            west = col-1
            north = row-1
            south = row+1

            if (col.eq.1) west = 16
            if (col.eq.16) east = 1
            if (row.eq.1) north = 16
            if (row.eq.16) south = 1

! try flip a spin at position (col,row)
! calculate energy (function?)

            Estart = J*(spin(col,row)*spin(east,row) &
                 + spin(col,row)*spin(west,row) &
                 + spin(col,row)*spin(col,north) &
                 + spin(col,row)*spin(col,south) )

            Eflip =  J*((-spin(col,row))*spin(east,row) &
                 + (-spin(col,row))*spin(west,row) &
                 + (-spin(col,row))*spin(col,north) &
                 + (-spin(col,row))*spin(col,south) )

! acceptance/rejectance for MC algorithm (define as subroutine?)

            if (Eflip .lt. Estart) then 
               spin(col,row) = -spin(col,row)
            else  
               call random_number(paccept)
               ptrial = EXP(-(Eflip-Estart)/kT)

               if (paccept .lt. ptrial) then
                  spin(col,row) = -spin(col,row)
               end if
            end if

         end do
      end do

   end do  ! loop over nMC

! JS: Now that the above has finished, and we are in equilibrium,
! JS: run another MC cycle again, this time to collect statistics
! JS:  In other words, another do-loop over nMC, row, col, etc.
! JS:  flipping spins as appropriate, etc.

! JS: BUT now, every time (mod(nMC, 50) .eq. 0), collect E, E^2, M to calculate
! JS: averages.  This E will have to be the total energy of the entire system
! JS: calculated over all the spins.  You might want to create functions to 
! JS: give you the value of E, M, E^2, or you could just do it inline.

! JS: When the datacollection is over (you have completed this second run of
! JS: nMC steps), it will be time to output the averages.

!OUTPUT:  Graphs; (Averages) Heat Capacity (<E^2>-<E>^2, ;derivative of <E(T)>, 
! Temperature, Energy <E> (to get Heat Capacity), Magnetization <M>
! JS: (You will want to report these as an intensive property,  
! JS:  i.e., <E>/N, C/N, <M>/N, where N is the total number of spins (16*16)

end do     !JS: loop over kT values (slowly cooling the system by lowering kT)

end program isingmodel

a
