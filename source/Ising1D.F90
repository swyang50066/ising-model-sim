Program Ising1D

! 1D-IsingModel with periodic boundary condition 

! Declare variables and arrays
Implicit None

Integer :: T, n
Real :: random, trial, Einitial, Eflip, dE , E
Real :: Eavg, E2avg, Mavg, M2avg, Cv, X! Evaluated expected physical quantities

Integer :: order, left, right
Integer :: nMonteCarlo ! Number of loop lf Monte-Carlo simulation
Real, Parameter :: J = 1.0    ! Define a feromagentic system.
Real, Parameter :: kB = 1.0 ! Boltzmann constant in simulation unit
Integer, Dimension(50) :: spin

! **** Random number generator ****
Call Random_Seed()

! Initialized spinor array
Do order = 1,50
        Call Random_Number(random)
      
        If (random .lt. 0.5) Then
                spin(order) = +1
        Else 
                spin(order) = -1
        End If
End Do

! **** Evaluation loop ****

! Loop over temperature
Open (Unit=1, File="Results.txt", Action="write", Status="replace") 
do T = 1,501,5      !along T = 1 ~ 500

! Reset the variables
Eavg = 0
E2avg = 0
Mavg = 0
M2avg = 0
Cv = 0
X = 0
        
        ! Do Monte-Carlo Simulation
        Do nMonteCarlo = 1,10000 ! 10000 trial

                ! Compare inital and fliped states
                Do order = 1,50

                        left = order-1
                        right = order-1

                        If (order .eq. 1) left = 16
                        If (order .eq. 16) right = 1 ! Declare boundary order as periodically

                        ! Try flip order's spinor and calculate initial & filped state energy
                        Einitial = J*(spin(order)*spin(left) &
                                + spin(order)*spin(right))

                        Eflip =  J*((-spin(order))*spin(order) &
                                + (-spin(order)*spin(left)))

                        dE = Eflip - Einitial

                        ! acceptance/rejectance for spin state
                        If (de .lt. 0.0) Then 
                                spin(order) = -spin(order)
                        Else  
                                Call Random_Number(random)
                                trial = EXP(-dE/(kB*T))

                                If (random .lt. trial) Then
                                spin(order) = -spin(order)
                                Einitial = Eflip
                                End If
                        End If
                
                        ! Sum interesting correction values
                        Eavg = Eavg - Einitial
                        E2avg = E2avg + Eavg**2
                        Mavg = Mavg + spin(order)
                        M2avg = M2avg + spin(order)**2       
                End Do
        End Do

! Find Monte-Carlo simulation average value
Eavg = Eavg/10000.0
E2avg = E2avg/10000.0
Mavg = Mavg/10000.0
M2avg = M2avg/10000.0

! Heat capacity, susceptibility
Cv = E2avg - Eavg**2
X = T*(M2avg - Mavg**2)

Write (1,*) T, " ", Eavg, " ", Mavg, " ", Cv, " ", X 
End Do
Close (1)

! JS: Now that the above has finished, and we are in equilibrium,
! JS: run another MC cycle again, this time to collect statistics
! JS:  In other words, another do-loop over nMC, row, col, etc.
! JS:  flipping spins as appropriate, etc.

! JS: BUT now, every time (mod(nMC, 50) .eq. 0), collect E, E^2, M to
! calculate
! JS: averages.  This E will have to be the total energy of the entire
! system
! JS: calculated over all the spins.  You might want to create functions
! to 
! JS: give you the value of E, M, E^2, or you could just do it inline.

! JS: When the datacollection is over (you have completed this second
! run of
! JS: nMC steps), it will be time to output the averages.

!OUTPUT:  Graphs; (Averages) Heat Capacity (<E^2>-<E>^2, ;derivative of
!<E(T)>, 
! Temperature, Energy <E> (to get Heat Capacity), Magnetization <M>
! JS: (You will want to report these as an intensive property,  
! JS:  i.e., <E>/N, C/N, <M>/N, where N is the total number of spins
! (16*16) 

End Program Ising1D


