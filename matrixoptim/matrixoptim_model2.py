#code written by Walker Law in July 2025
import numpy as np
import pandas as pd
import math
import time

starttime = time.time()

mp = .9382720813  # proton mass in GeV

df = pd.read_csv('//wsl.localhost/Ubuntu/home/walkerlaw/reu/matrixoptim/data_model2.csv', index_col=0, header=0, delimiter=',')

dW   = df.loc[:, 'dW'].to_numpy()
dEm  = df.loc[:, 'dEm'].to_numpy()
dPmx = df.loc[:, 'dPmx'].to_numpy()
dPmy = df.loc[:, 'dPmy'].to_numpy()
dPmz = df.loc[:, 'dPmz'].to_numpy()

Eb      = df.loc[:, 'Eb'].to_numpy()
Ef      = df.loc[:, 'Ef'].to_numpy()
thetae  = df.loc[:, 'thetae'].to_numpy()
thetap  = df.loc[:, 'thetap'].to_numpy()
phie    = df.loc[:, 'phie'].to_numpy()
phip    = df.loc[:, 'phip'].to_numpy()
Pp      = df.loc[:, 'Pp'].to_numpy()


M = np.empty((0, 7))

# columns are partials w.r.t. Eb, Ef, thetae, thetap, phie, phip, Pp
# rows are W, Em, Pmx, Pmy, Pmz in order from top to bottom for each run
for i in range(len(df)):
     W2 = mp**2 + 2*mp*(Eb[i]-Ef[i]) - 4*Eb[i]*Ef[i]*(math.sin(thetae[i]/2)**2)           #abs /#relative values
     A = [[Eb[i]*(2*mp-4*Ef[i]*(math.sin(thetae[i]/2)**2))/(2*math.sqrt(W2)),                   #dW/dEb
           Ef[i]*(-2*mp-4*Eb[i]*(math.sin(thetae[i]/2)**2))/(2*math.sqrt(W2)),                  #dW/dEf
           (-2*Eb[i]*Ef[i]*math.sin(thetae[i]/2)*math.cos(thetae[i]/2))/(math.sqrt(W2)),  #dW/dthetae
           0,                                                                             #dW/dthetap  
           0,                                                                             #dW/dphie
           Pp[i]*0,                                                                             #dW/dphip
           0],                                                                            #dW/dPp

          [Eb[i]*1,                                                                             #dEm/dEb                            
           Ef[i]*-1,                                                                            #dEm/dEf
           0,                                                                             #dEm/dthetae
           0,                                                                             #dEm/dthetap
           0,                                                                             #dEm/dphie
           Pp[i]*0,                                                                             #dEm/dphip
           -(Pp[i])/math.sqrt(Pp[i]**2+mp**2)],                                           #dEm/dPp

          [Eb[i]*0,                                                                             #dPmx/dEb             
           Ef[i]*-math.sin(thetae[i])*math.cos(phie[i]),                                        #dPmx/dEf
           -Ef[i]*math.cos(thetae[i])*math.cos(phie[i]),                                  #dPmx/dthetae
           -Pp[i]*math.cos(thetap[i])*math.cos(phip[i]),                                  #dPmx/dthetap
           Ef[i]*math.sin(thetae[i])*math.sin(phie[i]),                                   #dPmx/dphie
           Pp[i]*Pp[i]*math.sin(thetap[i])*math.sin(phip[i]),                                   #dPmx/dphip
           -math.sin(thetap[i])*math.cos(phip[i])],                                       #dPmx/dPp

          [Eb[i]*0,                                                                             #dPmy/dEb
           Ef[i]*-math.sin(phie[i]),                                                            #dPmy/dEf
           0,                                                                             #dPmy/dthetae
           0,                                                                             #dPmy/dthetap
           -Ef[i]*math.cos(phie[i]),                                                      #dPmy/dphie
           Pp[i]*-Pp[i]*math.cos(phip[i]),                                                      #dPmy/dphip
           -math.sin(phip[i])],                                                           #dPmy/dPp

          [Eb[i]*1,                                                                             #dPmz/dEb
           Ef[i]*-math.cos(thetae[i])*math.cos(phie[i]),                                        #dPmz/dEf
           Ef[i]*math.sin(thetae[i])*math.cos(phie[i]),                                   #dPmz/dthetae
           Pp[i]*math.sin(thetap[i])*math.cos(phip[i]),                                   #dPmz/dthetap
           Ef[i]*math.cos(thetae[i])*math.sin(phie[i]),                                   #dPmz/dphie
           Pp[i]*Pp[i]*math.cos(thetap[i])*math.sin(phip[i]),                                   #dPmz/dphip
           -math.cos(thetap[i])*math.cos(phip[i])]]                                       #dPmz/dPp
     M = np.vstack((M, A))

MT = np.transpose(M)

# Extract the uncertainty columns from the DataFrame
sigmas = []

for i in range(len(df)):
    sigmas.append(df.iloc[i]['sigdW'])
    sigmas.append(df.iloc[i]['sigdEm'])
    sigmas.append(df.iloc[i]['sigdPmx'])
    sigmas.append(df.iloc[i]['sigdPmy'])
    sigmas.append(df.iloc[i]['sigdPmz'])

# Convert to numpy array and square the values (since N = diag(sigma^2))
sigma_array = np.array(sigmas)

N = np.diag(sigma_array ** 2)
Ninv = np.linalg.inv(N)

y = np.empty((0, 1))

for i in range(len(df)):
     a = [[dW[i]],
          [dEm[i]],
          [dPmx[i]],
          [dPmy[i]],
          [dPmz[i]]]
     y = np.vstack((y, a))

# showing initial chi2 value, uncorrected spectrometer offsets (all zeros)
x = np.zeros((7, 1))

#initial chi2 calculation
initchi2 = (np.transpose(y-M@x))@Ninv@(y-M@x)
print(f'Initial chi²: {initchi2}, so initial chi2/dof value is {initchi2 / (len(df)*4-4)} for n*4 data points and 4 estimated parameters (since n*4 - 4 degrees of freedom).')

# Solve for the optimal spectrometer offsets
x = np.linalg.inv(MT @ Ninv @ M) @ (MT @ Ninv @ y)

# Calculate the final chi-squared value
finalchi2 = np.transpose(y-M@x)@Ninv@(y-M@x)

#Covariance matrix calculation (for spectrometer offset uncertainties)
Cov = np.linalg.inv(MT @ Ninv @ M)

endtime = time.time()

print(f'Final chi²: {finalchi2}, so final chi²/dof value is {finalchi2 / (len(df)*4-4)} for n*4 data points and 4 estimated parameters (since n*4 - 4 degrees of freedom).') 
print('Chi²/dof values less than 1.25 are generally considered acceptable.')
print(f'Spectrometer offsets: \ndelta E_b / E_b (factor - 1): {x[0]} ± {math.sqrt(Cov[0,0])} \ndelta E_f / E_f (factor - 1): {x[1]} ± {math.sqrt(Cov[1,1])} \ndelta theta_e (in radians): {x[2]} ± {math.sqrt(Cov[2,2])} \ndelta theta_p (factor - 1): {x[3]} ± {math.sqrt(Cov[3,3])} \ndelta phi_e (in radians): {x[4]} ± {math.sqrt(Cov[4,4])} \ndelta phi_p (in radians): {x[5]} ± {math.sqrt(Cov[5,5])} \ndelta P_p / P_p (factor - 1): {x[6]} ± {math.sqrt(Cov[6,6])}')
print(f'Time elapsed for calculation: {(endtime - starttime):.3} seconds')