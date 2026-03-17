import pandas as pd
import numpy as np
from scipy.signal import find_peaks, butter, filtfilt, hilbert
from scipy.fft import rfft, rfftfreq
import matplotlib.pyplot as plt

# --- 1️⃣ Načtení dat z XLSB ---
file_path = "poskoky.xlsb"
lower_time = 10
upper_time = 15
sheet_name = "List1"
citlivost = 1   #pokud je upravena citlivost a není změněna konstanta u měřicího zařízení

# 5 signálů
columns = ["1", "2", "3", "4", "5", "6"]  # uprav podle svých dat
df = pd.read_excel(file_path, sheet_name=sheet_name, engine="pyxlsb",
                   skiprows=49, names=columns)



Fs = 1200  # vzorkovací frekvence
threshold = 100  # mm/s^2 ořezání extrémních hodnot - nahradí nulou



## Čas dopočítávaný ze vzorkovací frekvence
N = len(df)
t = np.arange(N) / Fs  # čas v sekundách


# --- 3️⃣ Cyklus přes všechny signály ---
clean_data = pd.DataFrame()

for col in columns:
    signal = citlivost * df[col].to_numpy()
    
    # Odstranění odlehlých hodnot (ořez a nahrazení nulami)
    signal_clean = np.where(np.abs(signal) > threshold, 0, signal)

    # Uložení
    clean_data[col] = signal_clean
    

### --- 5️Vizualizace signálu  ---
##for col in columns:
##    plt.figure(figsize=(12,5))
##    plt.plot(t, clean_data[col], alpha=0.5, label="original data")
##    plt.xlabel("čas [s]")
##    plt.ylabel(f"zrychlení [mm/s^2]")
##    plt.legend()
##    plt.title(f"Snímač {col}")
##
### --- 6️Vizualizace všech signálů v jednom obrázku ---
##fig, axes = plt.subplots(len(columns), 1, figsize=(8.3, 11.7), sharex=True)
##
##for i, col in enumerate(columns):
##    ax = axes[i]
##    ax.plot(t, clean_data[col], alpha=0.5, label="originál")
##    ax.set_title(f"Snímač {col}")
##    ax.legend(fontsize=8, loc="upper right")
##
##axes[-1].set_xlabel("čas [s]")  # jen u spodního grafu
##plt.tight_layout()
##
##plt.show()

# --- 7 Omezení jen na signály, které dokmitávají
# Přidáme čas do DataFrame
clean_data["t"] = t

# Omezíme na interval 60–80 s a jen sloupce C a E
filtered = clean_data.loc[(clean_data["t"] >= lower_time) & (clean_data["t"] <= upper_time), ["t", "1", "2", "3", "4", "5", "6"]]

##plt.figure(figsize=(12,5))
##plt.plot(filtered["t"], filtered["1"], label="1")
##plt.plot(filtered["t"], filtered["2"], label="2")
##plt.xlabel("čas [s]")
##plt.ylabel("zrychlení [mm/s^2]")
##plt.legend()
##plt.title(f"Dokmitávání {lower_time}–{upper_time} s (1 a 2)")
##plt.show()

# --- 8 Kód pro stanovení tlumení

def bandpass(x, fs, f0, bw=1.0, order=3):
    # bw = pásmo ± v Hz kolem f0
    low = max(0.001, f0 - bw/2)
    high = min(fs/2 - 0.001, f0 + bw/2)
    b, a = butter(order, [low/(fs/2), high/(fs/2)], btype='band')
    return filtfilt(b, a, x)

def estimate_dominant_freq(x, fs):
    N = len(x)
    X = np.abs(rfft(x - np.mean(x)))
    freqs = rfftfreq(N, 1/fs)
    idx = np.argmax(X)
    return freqs[idx]

def log_decrement_from_peaks(t, x, fs, prominence=None, distance=None):
    # find peaks
    peaks, props = find_peaks(x, prominence=prominence, distance=distance)
    if len(peaks) < 2:
        return None
    A = np.abs(x[peaks])
    t_peaks = t[peaks]
    # odřízneme první 3 a poslední 3 vrcholy
    A = A[3:-3]
    t_peaks = t_peaks[3:-3]
    # compute per-cycle deltas
    deltas = np.log(A[:-1] / A[1:])
    deltas = deltas[deltas > 0]  # pojistka, kdyby se objevily nesmysly
    if len(deltas) == 0:
        return None
    Td_est = np.mean(np.diff(t_peaks))
    delta_mean = np.mean(deltas)
    delta_std = np.std(deltas, ddof=1) if len(deltas) > 1 else 0.0
    return {
        'peaks': peaks, 't_peaks': t_peaks, 'A': A,
        'deltas': deltas, 'Delta_mean': delta_mean,
        'Delta_std': delta_std, 'T_d': Td_est
    }

def log_decrement_from_envelope(t, x, T_d, trim=0.1):
    env = np.abs(hilbert(x - np.mean(x)))
    mask = env > 0

    ln_env = np.log(env[mask])
    t_env = t[mask]

    # 1️vynecháme část začátku a konce (např. 10 % - promenna "trim")
    n = len(t_env)
    i1 = int(trim * n)
    i2 = int((1 - trim) * n)
    t_env = t_env[i1:i2]
    ln_env = ln_env[i1:i2]

    if len(t_env) < 2:
        return None

    # 2️lineární fit ln(env) = ln(A0) - λ t
    p = np.polyfit(t_env, ln_env, 1)
    lam = -p[0]  # decay rate
    Delta = lam * T_d

    return {
        'lambda': lam,
        'Delta': Delta,
        'poly': p,
        't_env': t_env,
        'env': env
    }


def delta_to_zeta(Delta):
    return Delta / np.sqrt(4 * np.pi**2 + Delta**2)

# === Example: zpracovat oba kanaly (C a E) v DataFrame 'filtered' ===
channels = ["1", "2", "3", "4", "5", "6"]
results = {}
for ch in channels:
    t = filtered['t'].to_numpy()
    x = filtered[ch].to_numpy()
    x = x - np.mean(x)  # odečíst DC

    # 1) odhad dom. frekvence a bandpass (upravit bw podle potřeby, např. 2 Hz)
    f0 = estimate_dominant_freq(x, Fs)
    x_filt = bandpass(x, Fs, f0, bw=2.0, order=3)

    # 2) peaks method
    # prominence a distance lze upravit: distance ~= Fs / f0 * 0.5 (min. rozestup mezi vrcholy)
    distance = int(0.5 * Fs / f0) if f0>0 else None
    pk = log_decrement_from_peaks(t, x_filt, Fs, prominence=np.std(x_filt)*0.3, distance=distance)

    # 3) envelope method using T_d from peaks (pokud nejsou vrcholy, spočteme T_d = 1/f0)
    T_d = pk['T_d'] if (pk is not None and 'T_d' in pk and pk['T_d']>0) else 1.0/f0
    env_res = log_decrement_from_envelope(t, x_filt, T_d)

    # 4) vypočteme zeta
    Delta_peaks = pk['Delta_mean'] if pk is not None else None
    zeta_peaks = delta_to_zeta(Delta_peaks) if Delta_peaks is not None else None
    Delta_env = env_res['Delta']
    zeta_env = delta_to_zeta(Delta_env)

    results[ch] = {
        'f0': f0,
        'peak_res': pk,
        'env_res': env_res,
        'Delta_peaks': Delta_peaks,
        'Delta_peaks_std': pk['Delta_std'] if pk is not None else None,
        'zeta_peaks': zeta_peaks,
        'Delta_env': Delta_env,
        'zeta_env': zeta_env
    }

    # quick plots
    plt.figure(figsize=(10,5))
    plt.plot(t, x, alpha=0.3, label='raw')
    plt.plot(t, x_filt, label='bandpassed')
    ##plt.plot(t, np.abs(hilbert(x_filt)), '--', label='envelope')
    plt.plot(env_res['t_env'], np.exp(np.polyval(env_res['poly'], env_res['t_env'])),
         "--", lw=2, label="fitted envelope")
    if pk is not None:
        plt.plot(pk['t_peaks'], pk['A'], 'ro', label='used peaks')
    plt.title(rf"Channel {ch} — $f_0 = {f0:.3f}$ Hz,  $\delta_{{peaks}} = {Delta_peaks:.4f}$, $\zeta_{{peaks}} = {zeta_peaks:.4e}$, $\delta_{{env}} = {Delta_env:.4f}$, $\zeta_{{env}} = {zeta_env:.4e}$")
    plt.xlabel('t [s]')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{file_path}_{ch}_dokmit.png", dpi=300, bbox_inches="tight")

plt.show()

# tisk výsledků
for ch, r in results.items():
    print(f"=== {ch} ===")
    print(f"dominant f0 = {r['f0']:.3f} Hz")
    if r['Delta_peaks'] is not None:
        print(f"Delta (peaks) = {r['Delta_peaks']:.5f} ± {r['Delta_peaks_std']:.5f}")
        print(f"zeta (peaks)  = {r['zeta_peaks']:.5e}")
    print(f"Delta (envelope fit) = {r['Delta_env']:.5f}")
    print(f"zeta (envelope)      = {r['zeta_env']:.5e}")


