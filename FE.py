import matplotlib.pyplot as plt
import numpy as np
import math

# Function to calculate Level Symmetric Quadrature
def LQ(N):
    # Initialize the weights and ordinates
    N_pos = int(N/2)

    mu = np.zeros(N_pos)
    
    # Calculate the weights and ordinates
    mu[0] = 0.1389568
    for m in range(0, N_pos):
        print (m)
        mu[m] = math.sqrt(mu[0]**2 + 2*(m)*(1 - 3*mu[0]**2) / (N-2))
        
    return mu

mu = LQ(16)

print (mu)