import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint



# Parameters
alpha = 0.5    # Maximum rate of transcription
alpha_0 = 0.0005   # Baseline rate of transcription due to leakiness
n = 2         # Hill coefficient
beta = 20     # Degradation rate ratio (protein to mRNA)
Km = 40   # number of repressors that half-maximally repress a promoter

# ODE system (as defined in the article)
def repressilator(y,t):
    
    m_lacI, p_lacI, m_tetR, p_tetR, m_cI, p_cI = y
    
    dm_lacI_dt = -np.log(2)*m_lacI/2 + (alpha / (1 + (p_cI/Km)**n)) + alpha_0
    dp_lacI_dt = beta*m_lacI - np.log(2)*p_lacI/10
    
    dm_tetR_dt = -np.log(2)*m_tetR/2 + (alpha / (1 + (p_lacI/Km)**n)) + alpha_0
    dp_tetR_dt = beta*m_tetR - np.log(2)*p_tetR/10
    
    dm_cI_dt = -np.log(2)*m_cI/2 + (alpha / (1 + (p_tetR/Km)**n)) + alpha_0
    dp_cI_dt = beta*m_cI - np.log(2)*p_cI/10
    
    return [dm_lacI_dt, dp_lacI_dt, dm_tetR_dt, dp_tetR_dt, dm_cI_dt, dp_cI_dt]

# Initial conditions
m_lacI_0 = 2
p_lacI_0 = 0
m_tetR_0 = 2
p_tetR_0 = 0
m_cI_0 = 2
p_cI_0 = 0

y0 = [m_lacI_0, p_lacI_0, m_tetR_0, p_tetR_0, m_cI_0, p_cI_0]

t = np.linspace(0, 1000, 1000)

# Solving the ODE
solution = odeint(repressilator, y0, t)

# Plotting for protein concentrations with specified colors
plt.figure(figsize=(12, 6))

plt.plot(t, solution[:, 1], label='Protein LacI', color='blue')
plt.plot(t, solution[:, 3], label='Protein TetR', color='green')
plt.plot(t, solution[:, 5], label='Protein cI', color='red')

plt.xlabel('Time')
plt.ylabel('Protein Concentration')
plt.title('Repressilator Protein Dynamics')
plt.legend()
plt.show()

# New plotting for mRNA concentrations with specified colors
plt.figure(figsize=(12, 6))

plt.plot(t, solution[:, 0], label='mRNA LacI', color='darkblue')
plt.plot(t, solution[:, 2], label='mRNA TetR', color='orange')
plt.plot(t, solution[:, 4], label='mRNA cI', color='purple')

plt.xlabel('Time')
plt.ylabel('mRNA Concentration')
plt.title('Repressilator mRNA Dynamics')
plt.legend()
plt.show()


