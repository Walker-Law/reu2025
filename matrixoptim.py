import numpy as np
import pandas as pd
import math
import time

starttime = time.time()

# declare variables
mp = .9382720813  # proton mass in GeV
Eb = 10.542

df = pd.read_csv('//wsl.localhost/Ubuntu/home/walkerlaw/reu/matrixoptim/data.csv', index_col=0, header=0, delimiter=',')

dW1 = df.iloc[0, 0]
dEm1 = df.iloc[0, 1]
dPmx1 = df.iloc[0, 2]
dPmz1 = df.iloc[0, 3]

dW2 = df.iloc[1, 0]
dEm2 = df.iloc[1, 1]
dPmx2 = df.iloc[1, 2]
dPmz2 = df.iloc[1, 3]

dW3 = df.iloc[2, 0]
dEm3 = df.iloc[2, 1]
dPmx3 = df.iloc[2, 2]
dPmz3 = df.iloc[2, 3]

dW4 = df.iloc[3, 0]
dEm4 = df.iloc[3, 1]
dPmx4 = df.iloc[3, 2]
dPmz4 = df.iloc[3, 3]

dW5 = df.iloc[4, 0]
dEm5 = df.iloc[4, 1]
dPmx5 = df.iloc[4, 2]
dPmz5 = df.iloc[4, 3]

dW6 = df.iloc[5, 0]
dEm6 = df.iloc[5, 1]
dPmx6 = df.iloc[5, 2]
dPmz6 = df.iloc[5, 3]

Ef1 = df.iloc[0, 4]
pp1 = df.iloc[0, 5]
thetae1 = np.radians(df.iloc[0, 6])
thetap1 = np.radians(df.iloc[0, 7])
Ep1 = df.iloc[0, 8]

Ef2 = df.iloc[1, 4]
pp2 = df.iloc[1, 5]
thetae2 = np.radians(df.iloc[1, 6])
thetap2 = np.radians(df.iloc[1, 7])
Ep2 = df.iloc[1, 8]

Ef3 = df.iloc[2, 4]
pp3 = df.iloc[2, 5]
thetae3 = np.radians(df.iloc[2, 6])
thetap3 = np.radians(df.iloc[2, 7])
Ep3 = df.iloc[2, 8]

Ef4 = df.iloc[3, 4]
pp4 = df.iloc[3, 5]
thetae4 = np.radians(df.iloc[3, 6])
thetap4 = np.radians(df.iloc[3, 7])
Ep4 = df.iloc[3, 8]

Ef5 = df.iloc[4, 4]
pp5 = df.iloc[4, 5]
thetae5 = np.radians(df.iloc[4, 6])
thetap5 = np.radians(df.iloc[4, 7])
Ep5 = df.iloc[4, 8]

Ef6 = df.iloc[5, 4]
pp6 = df.iloc[5, 5]
thetae6 = np.radians(df.iloc[5, 6])
thetap6 = np.radians(df.iloc[5, 7])
Ep6 = df.iloc[5, 8]


M = [[-(Eb*Ef1*math.sin(thetae1))/mp, 0, -Eb, 0],
     [0, 0, -Ef1, -(pp1**2)/Ep1],
     [-Ef1*math.sin(thetae1), -pp1*math.sin(thetap1), -Ef1*math.cos(thetae1), -pp1*math.cos(thetap1)],
     [-Ef1*math.cos(thetae1), -pp1*math.cos(thetap1),  Ef1*math.sin(thetae1),  pp1*math.sin(thetap1)],

     [-(Eb*Ef2*math.sin(thetae2))/mp, 0, -Eb, 0],
     [0, 0, -Ef2, -(pp2**2)/Ep2],
     [-Ef2*math.sin(thetae2), -pp2*math.sin(thetap2), -Ef2*math.cos(thetae2), -pp2*math.cos(thetap2)],
     [-Ef2*math.cos(thetae2), -pp2*math.cos(thetap2),  Ef2*math.sin(thetae2),  pp2*math.sin(thetap2)],

     [-(Eb*Ef3*math.sin(thetae3))/mp, 0, -Eb, 0],
     [0, 0, -Ef3, -(pp3**2)/Ep3],
     [-Ef3*math.sin(thetae3), -pp3*math.sin(thetap3), -Ef3*math.cos(thetae3), -pp3*math.cos(thetap3)],
     [-Ef3*math.cos(thetae3), -pp3*math.cos(thetap3),  Ef3*math.sin(thetae3),  pp3*math.sin(thetap3)],

     [-(Eb*Ef4*math.sin(thetae4))/mp, 0, -Eb, 0],
     [0, 0, -Ef4, -(pp4**2)/Ep4],
     [-Ef4*math.sin(thetae4), -pp4*math.sin(thetap4), -Ef4*math.cos(thetae4), -pp4*math.cos(thetap4)],
     [-Ef4*math.cos(thetae4), -pp4*math.cos(thetap4),  Ef4*math.sin(thetae4),  pp4*math.sin(thetap4)],

     [-(Eb*Ef5*math.sin(thetae5))/mp, 0, -Eb, 0],
     [0, 0, -Ef5, -(pp5**2)/Ep5],
     [-Ef5*math.sin(thetae5), -pp5*math.sin(thetap5), -Ef5*math.cos(thetae5), -pp5*math.cos(thetap5)],
     [-Ef5*math.cos(thetae5), -pp5*math.cos(thetap5),  Ef5*math.sin(thetae5),  pp5*math.sin(thetap5)],

     [-(Eb*Ef6*math.sin(thetae6))/mp, 0, -Eb, 0],
     [0, 0, -Ef6, -(pp6**2)/Ep6],
     [-Ef6*math.sin(thetae6), -pp6*math.sin(thetap6), -Ef6*math.cos(thetae6), -pp6*math.cos(thetap6)],
     [-Ef6*math.cos(thetae6), -pp6*math.cos(thetap6),  Ef6*math.sin(thetae6),  pp6*math.sin(thetap6)]]

MT = np.transpose(M)

# Extract the uncertainty columns from the DataFrame
sigmas = []

for i in range(6):  # 6 runs
    sigmas.append(df.iloc[i]['sigdW'])
    sigmas.append(df.iloc[i]['sigdEm'])
    sigmas.append(df.iloc[i]['sigdPmx'])
    sigmas.append(df.iloc[i]['sigdPmz'])

# Convert to numpy array and square the values (since N = diag(sigma^2))
sigma_array = np.array(sigmas)

N = np.diag(sigma_array ** 2)
Ninv = np.linalg.inv(N)

Ninv = np.linalg.inv(N)

y = [[dW1],
     [dEm1],
     [dPmx1],
     [dPmz1],

     [dW2],
     [dEm2],
     [dPmx2],
     [dPmz2],

     [dW3],
     [dEm3],
     [dPmx3],
     [dPmz3],
     
     [dW4],
     [dEm4],
     [dPmx4],
     [dPmz4],

     [dW5],
     [dEm5],
     [dPmx5],
     [dPmz5],

     [dW6],
     [dEm6],
     [dPmx6],
     [dPmz6]]

#showing initial chi2 value, uncorrected spectrometer offsets (all zeros)
x = np.zeros((4, 1))

initchi2 = np.transpose(y-M@x)@Ninv@(y-M@x)

print(f'Initial chi2: {initchi2}, so initial chi2/dof value is {initchi2 / (len(df)*4-4)} for 24 data points and 4 estimated parameters (since 24 - 4 = 20 degrees of freedom).')


# Solve for the spectrometer offsets
x = np.linalg.inv(MT @ Ninv @ M) @ (MT @ Ninv @ y)

# Calculate the final chi-squared value
finalchi2 = np.transpose(y-M@x)@Ninv@(y-M@x)

print(f'Final chi²: {finalchi2}, so final chi²/dof value is {finalchi2 / (len(df)*4-4)} for n*4 data points and 4 estimated parameters (since n*4 - 4 degrees of freedom).\nChi²/dof values less than 1.25 are generally considered acceptable.')

print(f'Spectrometer offsets: \n{x}')

endtime = time.time()

print(f'Time taken for calculation: {endtime - starttime:.4} seconds')