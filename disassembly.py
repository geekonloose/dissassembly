import numpy as np
import math
from matplotlib import pyplot as plt
from temperature_profile import temp_profile
number = int(1E6) # number of iteration




# %% Initial temperature profile

temperature = temp_profile()
T = np.zeros(number)
T[2] = 2000 #initial temperature
T[1] = 1600




# %% Values for constants

alpha = 0.0035*1 # check number Reactivity insertion rate
beta = 0.0035 # Delay Neutron Fraction
rho0 = 0.0035 # initial reactivity insertion
sigmatr = 1.0
D = 1/3*sigmatr # value to be found
sigma = 1.0 # value to be found
KD = 1 #value of doppler coefficient is unknown
l = 1E-7  #  prompt generation time
dt = 1E-7  #  time steps
lmbda = 0.08
rhouo2 = 10970 # UO_2 density
rhopuo2 = 11460 # PuO_2 density
p = 0.1515 # weight fraction of plutonium in fuel
rho273 = (rhouo2*rhopuo2)/(p*rhouo2+(1-p)*rhopuo2) #density of fuel at 273
hfg = 2000
gas_const = 1.380648813131313131E-23 #boltzman gas constant
#gas_const = 8.31/238.0

a = []

# %% constant value of cp for UO2 (all value in j/gramm)
c1u = -78.4303E-3
c2u = 193.238E-3
c3u = 162.8647E-3
c4u = -104.0014E-3
c5u = 29.2056E-3
c6u = -1.9507E-3
c7u = 2.6441E-3


# %% constant value of cp for PUO2 (all value in j/gramm)

c1p = -118.2062E-3
c2p = 311.7866E-3
c3p = 19.629E-3
c4p = -0.752E-3
c5p = 0
c6p = 0
c7p = 7.0131E-3


# %% Variable initialisation

Q = 1000E6*np.ones(number)  #  reactor power
c = np.zeros(number)
y = np.zeros(number) # precursor concentration
rho = np.zeros(number) # total reactivity insertion
rho_in = np.zeros(number) # reactivity insertion
rho_dop = np.zeros(number)  # doppler reactivity
rho_dis = np.zeros(number) # reactivity addition due to fuel material displacement
E0 = np.zeros(number) #decay energy
k = np.zeros(number) # multiplication factor
k_exp = np.zeros(number) # k due to fuel expansion
k_dop = np.zeros(number)
R =  np.ones(number)

R[2] = float(((9.0/4.0)*(1.0/3.14))**(1.0/3.0)) #initial value of the bubble radius

E =np.ones(number) # additional energy added by power exvursion



LMBDA = np.zeros(number) # Mean Generation time
P = np.zeros(number) # pressure
P[2] = 2E5 # pressure in Pa
P[1] = 1E5
density_of_fuel = np.zeros(number)
vr = np.zeros(number) # velocity of the fuel displacement
dR = np.zeros(number) # change in the radius due to displacement
u = np.zeros(number)
cpuo2 = np.zeros(number)
cppuo2 = np.zeros(number)
cp = np.zeros(number)
m = np.zeros(number)
k_expansion = np.zeros(number)
V = np.zeros(number)
rho_expansion = np.zeros(number)







# %% Main Program Starts from here






for i in range(2,number-1):

    print('i',i)

#%% Total reactivity calculation

    k_expansion[i] = 2.5 * (1/(1+(3.14**2 / R[i]**2) + (D/sigma)))

    rho_expansion[i] = 1-(1/k_expansion[i])

    rho_in[i] = alpha * i * dt # reactivity insertion

#    KD = (T[i])**(-1/3.0) # doppler coefficient

    KD = 1

    E0[i] = 0.0631*1100E6*(i*dt)**(-0.1322)

    k_dop[i]  =  KD*(1 - (E0[i]**(1.0/2)/(E0[i]+E[i])**(1.0/2))) # multiplication coefficient due to doppler effect

    rho_dop[i] = 1-(1.0/k_dop[i]) # reactivity due to doppler effect

    # kd = (k[i] - k[i-1])/(T[i]- T[i-1])
    # rho_dop[i] = kd*np.log(T[i]/T[i-1])


    rho[i] = rho0 + rho_in[i] + rho_dop[i] + rho_dis[i] + 4*beta

    k[i] = 1/(1-rho[i])

#    LMBDA[i] = l/k[i]
    LMBDA[i] = 1E-6

#%% Point Kineticss equation

    Q[i+1] = Q[i] + dt*((rho[i]-beta)/LMBDA[i]  +  lmbda*c[i])

    c[i+1] = c[i] + dt*(beta/LMBDA[i] - lmbda*c[i])


#%% Pressure

    P[i] = 1E6*math.exp(-78847/T[i] + 53.152 - 4.208 * math.log(T[i]))


#    P[i] = P[1] * math.exp(hfg/(gas_const*T[i]) + hfg/(gas_const*T[1]))

#%% Density calculation

    Lt_divide_L273 =  0.99672 + 1.179E-5*T[i] -2.429E-9*T[i]**2 + 1.219E-12*T[i]**3

#    density_of_fuel[i] =  (rho273 * Lt_divide_L273**3) # density at T
    density_of_fuel[i] = 11000E3 # in gm/m3
    vr[i+1] = vr[i] + dt*(-1/(density_of_fuel[i]) * (P[i]-P[i-1])/dt * (dt/(R[i]-R[i-1]))) ###### check LUNCH

    dR[i+1] = dR[i] + dt*vr[i]

    R[i+1] = R[i] + dR[i]


    V[i] = 4.0/3 * 3.14 * R[i]**3

    m[i] = 11000E3*V[i]  # in gm/m3

    # %% calculation of cp

    y = 0.8
    t = T[i]/1000


    cpuo2[i] = c2u + 2*c3u*t + 3*c4u*t**2 + 4*c5u*t**3 + 5*c6u*t**4 - c7u*t**(-2)

    cppuo2[i] = c2p + 2*c3p*t + 3*c4p*t**2 + 4*c5p*t**3 + 5*c6p*t**4 - c7p*t**(-2)

    cp[i] = y * cpuo2[i] + (1-y) * cppuo2[i]

    a.append(Q[i]/cp[i]*m[i])


    # %% calculation of temperature

    T[i+1] = T[i] + Q[i]/cp[i]*m[i]

#    E[i] = Q[i]



# %% Plot

plt.plot([dt*i for i in range(0,i+2)],Q)







