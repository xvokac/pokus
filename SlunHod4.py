import numpy as np
import matplotlib.pyplot as plt

phi = np.deg2rad(50)

# gnómon (osa Země nahoru)
G = np.array([0, -np.cos(phi), -np.sin(phi)])

# normála stěny (jih)
n = np.array([0, -1, 0])

def sun_vector(delta):
    # poledne!
    return np.array([0, np.cos(delta), np.sin(delta)])

def shadow_point(delta):
    S = sun_vector(delta)

    # průsečík
    t = -np.dot(n, G) / np.dot(n, S)
    P = G + t * (S)

    return P

# =========================
# VÝPOČET
# =========================
cases = [(0, "Horizont"),
         (-90+np.rad2deg(phi)-23.5, "Zima"),
         (-90+np.rad2deg(phi)+0, "Rovnodennost"),
         (-90+np.rad2deg(phi)+23.5, "Léto")]

plt.figure(figsize=(6,6))

for d, name in cases:
    P = shadow_point(np.deg2rad(d))

    # bokorys:
    x = np.dot(P, n)        # vzdálenost od stěny
    z = P[2]                # výška

    print(f"{name}: z = {z:.4f}")

    plt.plot(x, z, 'o', label=name)

# gnómon
g_x = np.dot(G, n)
g_z = G[2]

plt.plot([0, g_x], [0, g_z], 'r-', label="gnómon")

# stěna
plt.plot([0,0], [-0.5,1.5], 'k-', linewidth=2)

plt.axhline(0)
plt.axvline(0)

plt.gca().set_aspect('equal')
plt.legend()
plt.title("Bokorys – poledne")
plt.xlabel("od stěny")
plt.ylabel("výška")

plt.show()
