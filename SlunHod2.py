import numpy as np
import matplotlib.pyplot as plt

# =========================
# PARAMETRY
# =========================
phi = np.deg2rad(50.0)   # zeměpisná šířka
D   = np.deg2rad(0.0)   # azimut stěny (0 = jih)

L = 1.0  # délka gnómonu

# =========================
# ZÁKLADNÍ VEKTORY
# =========================
n = np.array([np.sin(D), -np.cos(D), 0])   # normála stěny
G = np.array([0, np.cos(phi), np.sin(phi)])  # gnómon

# =========================
# SOUSTAVA NA STĚNĚ
# =========================
z = np.array([0,0,1])

v = z - np.dot(z, n)*n
v /= np.linalg.norm(v)

u = np.cross(v, n)
u /= np.linalg.norm(u)

def to2D(P):
    x = np.dot(P, u)
    y = -np.dot(P, v)
    return np.array([x, y])

# =========================
# SLUNCE
# =========================
def sun_vector(H, delta):
    # Lokální soustava: x=východ, y=sever, z=nahoru.
    # H ... hodinový úhel (kladný na západ), delta ... deklinace Slunce.
    return np.array([
        np.cos(delta) * np.sin(H),
        np.cos(phi) * np.sin(delta) - np.sin(phi) * np.cos(delta) * np.cos(H),
        np.sin(phi) * np.sin(delta) + np.cos(phi) * np.cos(delta) * np.cos(H)
    ])



def shadow_point(H, delta):
    S = sun_vector(H, delta)

    denom = np.dot(n, S)
    if abs(denom) < 1e-12:
        return None

    # Průsečík přímky G + λS se stěnou n·P = 0.
    # (S je směr od pozorovatele ke Slunci.)
    t = -np.dot(n, G) / denom
    if t <= 0:
        return None

    return G + t * S



# =========================
# FIGURY
# =========================
fig, axs = plt.subplots(1, 3, figsize=(15,5))

# =========================
# 1️⃣ SLUNEČNÍ HODINY
# =========================
ax = axs[0]

roman_hours = {
    6: "VI", 7: "VII", 8: "VIII", 9: "IX", 10: "X", 11: "XI",
    12: "XII", 13: "XIII", 14: "XIV", 15: "XV", 16: "XVI",
    17: "XVII", 18: "XVIII"
}

for hour in range(6, 19):
    H = np.deg2rad(15 * (hour - 12))

    nh = np.array([
        np.sin(H),
        np.cos(H) * np.sin(phi),
        np.cos(H) * np.cos(phi)
    ])

    d = np.cross(n, nh)
    p = to2D(d)
    p /= np.linalg.norm(p)

    # Vykreslit jen spodní část hodinových čar (y <= 0).
    if p[1] > 0:
        p = -p

    # Protáhnout čáry až na hranici obdélníku: -3 < x < 3, -3 < y < 0.
    tx = np.inf if abs(p[0]) < 1e-12 else 3.0 / abs(p[0])
    ty = np.inf if abs(p[1]) < 1e-12 else 3.0 / abs(p[1])
    Lline = min(tx, ty)

    ax.plot([0, p[0] * Lline], [0, p[1] * Lline], 'k-')

    label_scale = min(Lline * 1.05, Lline + 0.2)
    ax.text(p[0] * label_scale, p[1] * label_scale, roman_hours[hour],
            ha='center', va='center', fontsize=9)


# křivky (pro fixní deklinaci δ)
all_front_points = []
for d, name in [(-23.44, "Zima"),
                (0.0, "Rovnodennost"),
                (23.44, "Léto")]:
    pts = []
    for h in np.linspace(-90, 90, 400):
        P = shadow_point(np.deg2rad(h), np.deg2rad(d))
        if P is None:
            continue
        # Slunce musí být nad horizontem
        if sun_vector(np.deg2rad(h), np.deg2rad(d))[2] <= 0:
            continue
        pts.append(to2D(P))

    if len(pts) > 2:
        pts = np.array(pts)
        ax.plot(pts[:,0], pts[:,1], label=name)
        all_front_points.append(pts)

# horizont: body pro východ/západ (výška Slunce = 0°) a reálné deklinace.
horizon_pts = []
for d in np.linspace(-23.44, 23.44, 800):
    delta = np.deg2rad(d)
    arg = -np.tan(phi) * np.tan(delta)
    if arg < -1 or arg > 1:
        continue

    H0 = np.arccos(arg)
    for H in (-H0, H0):
        P = shadow_point(H, delta)
        if P is None:
            continue
        horizon_pts.append(to2D(P))

if len(horizon_pts) > 2:
    horizon_pts = np.array(horizon_pts)
    ax.plot(horizon_pts[:,0], horizon_pts[:,1], label="Horizont")
    all_front_points.append(horizon_pts)

# gnómon
g2 = to2D(G)
ax.plot([0, g2[0]], [0, g2[1]], 'r-', linewidth=2)

ax.set_title("Stěna (front view)")
ax.set_aspect('equal')
ax.set_xlim(-4, 4)
ax.set_ylim(-4, 1)
ax.axhline(0, linewidth=0.5)
ax.axvline(0, linewidth=0.5)
ax.legend()

# =========================
# 2️⃣ PŮDORYS (TOP VIEW)
# =========================
ax = axs[1]

# projekce do XY roviny
Gn = np.array([-G[0], -G[1]])
nn = np.array([n[0], n[1]])

# normála stěny
ax.plot([0, nn[0]], [0, nn[1]], 'k-', label="normála stěny")

# gnómon
ax.plot([0, Gn[0]], [0, Gn[1]], 'r-', label="gnómon")

# stěna = kolmice na normálu
wall_dir = np.array([-nn[1], nn[0]])
ax.plot([-wall_dir[0], wall_dir[0]],
        [-wall_dir[1], wall_dir[1]], 'k--', label="stěna")

ax.set_title("Půdorys")
ax.set_aspect('equal')
ax.legend()


# =========================
# 3️⃣ BOKORYS (SIDE VIEW)
# =========================
ax = axs[2]

# vzdálenost od stěny (projekce do normály)
g_normal = np.dot(G, n)

# výška
g_height = G[2]

# svislice
ax.plot([0,0], [0,1], 'k--', label="svislice")

# gnómon
ax.plot([0, -g_normal], [0, -g_height], 'r-', label="gnómon")

# stěna
ax.plot([0,0], [-1,1], 'k-', linewidth=2)

ax.set_title("Bokorys")
ax.set_aspect('equal')
ax.legend()


plt.tight_layout()
plt.show()
