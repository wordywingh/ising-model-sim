import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def weighted_die(num_steps):
    """
    Returns the simulated amount of money earned after playing the game given the number of simulation steps
    """
    #state the probability of each die face
    probability_map = {
        1:3,
        2:3,
        3:1,
        4:1,
        5:1,
        6:1
    }
    #initialize randomly the state of the dice
    #physically this would be the dice starts in any given position from 1 - 6
    current_state = np.random.randint(1,7)
    #initialize the money we have
    money = 0
    
    for x in range(num_steps):
        #propose a new state
        proposed_state = np.random.randint(1,7)
        #if the current state is the same as the proposed state
        if current_state == proposed_state:
            #accept the move aka don't do anything because the new state is the same as the old state
            pass
        else:
            p_accept = min(1,probability_map[proposed_state]/probability_map[current_state])
            if np.random.uniform(0,1) < p_accept:
                #accept the move with probabiliy p_accept
                current_state = proposed_state
            else:
                #if we don't accept the move we keep the current state
                pass
        if current_state == 1 or current_state == 2:
            money += 1
        else:
            money -= 1
            
    return money/num_steps

def two_dim_ising(L, temp, num_steps):
    """
    Returns the resulting lattice of a MCMC simulation performed on a randomly initialized lattice after num_steps
    and given the temp
    """
    
    lattice = random_lattice(L)
    
    return two_dim_ising_step(lattice,temp,num_steps)

def two_dim_ising_step(lattice,temp,num_steps):
    """
    Returns the resulting lattice of a MCMC simulation performed on a given lattice after num_steps
    and given the temp
    """
    
    L = len(lattice)
    
    for x in range(num_steps):
        #pick a random location in the lattice
        i = np.random.randint(L)
        j = np.random.randint(L)
        #calculate energy by using the formula and minding neighbor states. mod will take care of edge case looping around
        #terms for neighbors are top, bottom, left, and right
        delta_E = 2*lattice[i,j]*(lattice[(i-1)%L,j]+lattice[(i+1)%L,j]+lattice[i,(j-1)%L]+lattice[i,(j+1)%L])
        if delta_E <= 0:
            #accept the move if energy is less than or equal to zero
            lattice[i,j] *= -1
        else:
            #else the energy is greater than zero, so accept the move with probability exp(-delta_E/T)
            if np.random.uniform(0,1) < np.exp(-delta_E/temp):
                lattice[i,j] *= -1
            else:
                pass
    
    #after we have done num_steps, return the resulting lattice
    return lattice

def random_lattice(L):
    """
    Generates a random square lattice of side length L and returns it
    """
    
    #initialize a square lattice with random spins -1 or 1
    lattice = 2*np.random.randint(2,size = (L,L)) - 1
    
    return lattice

def aligned_lattice(L):
    """
    Generates a fully allgined lattice of spin +1 with side length L and returns it
    """
    
    lattice = np.ones((L,L))
    return lattice

def calculate_E(lattice):
    """
    Calculates the observable E
    """
    energy = 0
    L = len(lattice)
    for i in range(len(lattice)):
        for j in range(len(lattice)):
            energy += -(lattice[(i-1)%L,j]+lattice[(i+1)%L,j]+lattice[i,(j-1)%L]+lattice[i,(j+1)%L])*lattice[i,j]
    return energy/2

def calculate_S(lattice):
    """
    Calculates the observable S
    """
    return abs(np.sum(lattice))

def Onsager(x):
    """
    Onsager's exact solution for the plots I'm going to make
    """
    if x < 2.2692:
        return ( 1 - ( np.sinh(2/(x)) ) ** (-4) ) ** (1/8)
    else:
        return 0
