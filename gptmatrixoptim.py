import numpy as np
import pandas as pd

# === Load CSV file ===
df = pd.read_csv("//wsl.localhost/Ubuntu/home/walkerlaw/reu/matrixoptim/data_raw.csv")  # Replace with your actual filename

# === Define input variables ===
# Residuals: [dW, dEm, dPmx, dPmy, dPmz]
Y = df[["dW", "dEm", "dPmx", "dPmy", "dPmz"]].to_numpy().T  # shape (5, N)
Y = Y.flatten()  # shape (5*N,)

# Independent variables: [Ef, Pp, thetae, thetap, phie, phip]
X = df[["Ef", "Pp", "thetae", "thetap", "phie", "phip"]].to_numpy().T  # shape (6, N)

# Nominal Eb
Eb = df["Eb"].values[0]
mp = 0.938272

# Spectrometer offset sigmas (same for all rows)
sigx = df[["sigEf", "sigPp", "sigthetae", "sigthetap", "sigphie", "sigphip"]].iloc[0].to_numpy()

# === Compute matrix M (5 equations × 6 parameters × N events) ===
N = X.shape[1]
M = np.zeros((5 * N, 6))

for i in range(N):
    Ef, Pp, thetae, thetap, phie, phip = X[:, i]
    pe = np.sqrt(Ef**2 - 0.000511**2)  # Assume electron mass
    pp = Pp
    Epp = np.sqrt(Pp**2 + mp**2)

    # Derivatives of W, Em, Pmx, Pmy, Pmz wrt Ef, Pp, angles
    row = i * 5

    # ∂W / ∂params
    Q2 = 4 * Ef * (Eb - Ef) * np.sin(thetae / 2)**2
    nu = Eb - Ef
    dW_dEf = (2 / (2*mp + 2*nu)) * ((2*mp * (1 - Eb/Ef)) + Q2/Ef)
    M[row+0, 0] = dW_dEf  # dW/dEf
    M[row+0, 2] = 4*Eb*Ef*np.sin(thetae/2)*np.cos(thetae/2)/mp  # dW/dthetae

    # ∂Em / ∂params
    M[row+1, 0] = -1
    M[row+1, 1] = -Pp / Epp
    M[row+1, 3] = Pp * np.sin(thetap) / Epp

    # ∂Pmx / ∂params
    M[row+2, 0] = -pe * np.sin(thetae)
    M[row+2, 1] = -pp * np.sin(thetap)
    M[row+2, 2] = -pe * np.cos(thetae)
    M[row+2, 3] = -pp * np.cos(thetap)

    # ∂Pmy / ∂phie, phip only
    M[row+3, 4] = -pe
    M[row+3, 5] = pp

    # ∂Pmz / ∂params
    M[row+4, 0] = pe * np.cos(thetae)
    M[row+4, 1] = -pp * np.cos(thetap)
    M[row+4, 2] = -pe * np.sin(thetae)
    M[row+4, 3] = pp * np.sin(thetap)

# === Build N matrix: diagonal with sigy² = sum_j (∂y_i/∂x_j)^2 * σ_xj^2 ===
sigx2 = sigx**2
sigy2 = np.sum((M**2) * sigx2[np.newaxis, :], axis=1)
Ninv = np.diag(1 / sigy2)

# === Compute least-squares solution ===
MT_Ninv = M.T @ Ninv
cov = np.linalg.inv(MT_Ninv @ M)        # Covariance matrix of x
x_best = cov @ (MT_Ninv @ Y)            # Best-fit spectrometer offsets
x_uncert = np.sqrt(np.diag(cov))        # Uncertainties

# === Compute chi-square ===
residuals = Y - M @ x_best
chi2 = residuals.T @ Ninv @ residuals
ndf = len(Y) - len(x_best)

# === Print results ===
params = ["Ef", "Pp", "thetae", "thetap", "phie", "phip"]
print("Fitted Spectrometer Offsets:")
for name, val, err in zip(params, x_best, x_uncert):
    print(f"{name:8s}: {val:+.6f} ± {err:.6f}")

print(f"\nChi-squared / dof = {chi2:.3f} / {ndf} = {chi2/ndf:.3f}")
