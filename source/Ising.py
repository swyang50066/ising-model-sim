import numpy as np


def init_spin_state_2d(nsize=16):
    """Initialize spin state"""
    return 2*np.random.randint(2, size=(nsize, nsize)) - 1


def mcmh_algorithm(state, beta=1):
    """Apply Monte Carlo Metropolis-Hastings algorithm"""
    # Get input dimensions
    height, width = state.shape

    energy = 0
    for i in range(height):
        for j in range(width):
            # Periodic neighbors
            up, down, left, right = (
                (i - 1) % height, (i + 1) & height,
                (j - 1) % width, (j + 1) & width
            )

            # Spin interaction energies
            e_spin_init = J*(
                state[i, j]*state[up, j] 
                + state[i, j]*state[down, j]
                + state[i, j]*state[i, left]
                + state[i, j]*state[i, right]
            )
            e_spin_flip = J*(
                -state[i, j]*state[up, j]
                - state[i, j]*state[down, j]
                - state[i, j]*state[i, left]
                - state[i, j]*state[i, right]
            ) 
            delta_e = e_spin_flip - e_spin_init
            energy += e_spin_flip 

            # Metropolis updates
            if delta:
                state[i, j] = -state[i, j]
            elif np.random.rand() < np.exp(-beta*delta_e) :
                state[i, j] = -state[i, j]
            else:
                pass

    return state, energy/nsize**2., np.sum(state)/nsize**2.


def run_simulation(num_iter, beta=16):
    """Run Ising model simulation"""
    # Randomly initialize spin state
    state = init_spin_state(nsize=10)

    for step in num_iter:
        # Update state 
        state, energy, mag = mcmh_algorithm(state, beta)

        # Update values
        if step < relaxation:       # Relaxation 
            continue
        elif step == relaxation:    # init values
            avg_energy = energy/num_iter
            avg_mag = mag/num_iter
            sqd_energy = energy**2./num_iter
            sqd_mag = mag**2./num_iter
        else:
            avg_energy += energy/num_iter
            avg_mag += mag/num_iter
            sqd_energy += energy**2./num_iter
            sqd_mag += mag**2./num_iter

    specific_heat = beta*(sqd_energy - avg_energy**2.)
    magnetic_susceptibility = beta**2.*(sqd_mag - avg_mag**2.)

