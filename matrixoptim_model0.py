#pretty much the same as matrixoptim2.py, but with some changes to the derivatives in the matrices
import numpy as np
import pandas as pd
import math
import time

starttime = time.time()

# declare variables
mp = .9382720813  # proton mass in GeV

df = pd.read_csv('//wsl.localhost/Ubuntu/home/walkerlaw/reu/matrixoptim/data_raw.csv', index_col=0, header=0, delimiter=',')

dW   = df.loc[:, 'dW'].to_numpy()
dEm  = df.loc[:, 'dEm'].to_numpy()
dPmx = df.loc[:, 'dPmx'].to_numpy()
dPmz = df.loc[:, 'dPmz'].to_numpy()
Ef     = df.loc[:, 'Ef'].to_numpy()
Pp     = df.loc[:, 'Pp'].to_numpy()
thetae = np.deg2rad(df.loc[:, 'thetae'].to_numpy())
thetap = np.deg2rad(df.loc[:, 'thetap'].to_numpy())
Eb     = df.loc[:, 'Eb'].to_numpy()

M = np.empty((0, 4))

# columns are partials w.r.t. Eb, Ef, thetae, thetap
# rows are W, Em, Pmx, Pmz in order from top to bottom for each run
for i in range(len(df)):
     W2 = mp**2 + 2*mp*(Eb[i]-Ef[i]) - 4*Eb[i]*Ef[i]*(math.sin(thetae[i]/2)**2)
     print(math.sqrt(W2))
     A = [[Eb[i]*((mp-2*Ef[i]*(math.sin(thetae[i]/2)**2))/(math.sqrt(W2))), Ef[i]*((-mp-2*Eb[i]*(math.sin(thetae[i]/2)**2))/(math.sqrt(W2))), (-2*Eb[i]*Ef[i]*math.sin(thetae[i]/2)*math.cos(thetae[i]/2))/(math.sqrt(W2)), 0],
          [Eb[i]*1, Ef[i]*-1, 0, 0],
          [0, Ef[i]*-math.sin(thetae[i]), -Ef[i]*math.cos(thetae[i]), -Pp[i]*math.cos(thetap[i])],
          [Eb[i]*1, Ef[i]*(-math.cos(thetae[i])), Ef[i]*math.sin(thetae[i]), Pp[i]*math.sin(thetap[i])]]
     M = np.vstack((M, A))


MT = np.transpose(M)

# #Extract the uncertainty columns from the DataFrame
# sigmas = []
# 
# for i in range(len(df)):  # 6 runs
#     sigmas.append(df.iloc[i]['sigdW'])
#     sigmas.append(df.iloc[i]['sigdEm'])
#     sigmas.append(df.iloc[i]['sigdPmx'])
#     sigmas.append(df.iloc[i]['sigdPmz'])
# 
# #Convert to numpy array and square the values (since N = diag(sigma^2))
# sigma_array = np.array(sigmas)
# 
# N = np.diag(sigma_array ** 2)
# Ninv = np.linalg.inv(N)

#uncertainties for the spectrometer offsets (Eb, Ef, thetae, thetap)
sigx = np.array([0.0001, 0.0001, 0.0005, 0.0005])  # example uncertainties

# Square the sigx for later use
sigx2 = sigx**2

# Propagate spectrometer uncertainties to kinematic residual uncertainties (sigy2)
# sigy2_i = sum_j (d y_i / d x_j)^2 * sigma_x_j^2
sigy2 = np.sum((M**2) * sigx2[np.newaxis, :], axis=1)

# Construct the inverse covariance matrix Ninv
N = np.diag(sigy2)
Ninv = np.linalg.inv(N)

y = np.empty((0, 1))

for i in range(len(df)):
     a = [[dW[i]],
          [dEm[i]],
          [dPmx[i]],
          [dPmz[i]]]
     y = np.vstack((y, a))

# showing initial chi2 value, uncorrected spectrometer offsets (all zeros)
x = np.zeros((4, 1))

initchi2 = (np.transpose(y-M@x))@Ninv@(y-M@x)

print(f'Initial chi²: {initchi2}, so initial chi2/dof value is {initchi2 / ((4-1)*(len(df)-1))} for 4 kinematic settings and {len(df)} runs (since (4-1)*({len(df)}-1) = {4*(len(df)-1)} degrees of freedom).')

# Solve for the spectrometer offsets
x = np.linalg.inv(MT @ Ninv @ M) @ (MT @ Ninv @ y)

# Calculate the final chi-squared value
finalchi2 = np.transpose(y-M@x)@Ninv@(y-M@x)

endtime = time.time()

print(f'Final chi²: {finalchi2}, so final chi²/dof value is {finalchi2 / ((4-1)*(len(df)-1))} for 4 kinematic settings and {len(df)} runs (since (4-1)*({len(df)}-1) = {4*(len(df)-1)} degrees of freedom).') 
print('Chi²/dof values less than 1.25 are generally considered acceptable.')
print(f'Spectrometer offsets: \ndelta E_b / E_b (factor - 1): {x[0]} \ndelta E_f / E_f (factor - 1): {x[1]} \ndelta theta_e (in radians): {x[2]} \ndelta theta_p (in radians): {x[3]}')
print(f'Time elapsed for calculation: {(endtime - starttime):.3} seconds')