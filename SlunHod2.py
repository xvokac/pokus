import sys

import numpy as np
import matplotlib.pyplot as plt

# =========================
# PARAMETRY
# =========================


def parse_args(argv):
    default_lat = 50.0
    default_azimuth = 0.0

    if len(argv) == 3:
        try:
            lat = float(argv[1])
            azimuth = float(argv[2])
            return lat, azimuth
        except ValueError:
            print("Neplatné argumenty – očekávám čísla.")
    elif len(argv) != 1:
        print("Špatný počet argumentů.")

    print("Použití: python SlunHod2.py <zeměpisná_šířka_N> <azimut_stěny>")
    print("  <zeměpisná_šířka_N>: severní zeměpisná šířka ve stupních (např. 49.5)")
    print("  <azimut_stěny>: azimut stěny ve stupních, 0 = jih,")
    print("                  kladné hodnoty směrem k západu, záporné k východu")
    print(f"Spouštím s implicitními hodnotami: {default_lat} a {default_azimuth}.\n")

    return default_lat, default_azimuth


lat_deg, azimuth_deg = parse_args(sys.argv)

phi = np.deg2rad(lat_deg)   # zeměpisná šířka
D   = np.deg2rad(azimuth_deg)   # azimut stěny (0 = jih)

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

def to_roman(number):
    vals = [
        (1000, "M"), (900, "CM"), (500, "D"), (400, "CD"),
        (100, "C"), (90, "XC"), (50, "L"), (40, "XL"),
        (10, "X"), (9, "IX"), (5, "V"), (4, "IV"), (1, "I")
    ]
    result = []
    n = number
    for value, symbol in vals:
        while n >= value:
            result.append(symbol)
            n -= value
    return "".join(result)


def is_hour_visible(hour):
    """Vrací True, pokud se hodinová čára může během roku objevit na stěně."""
    H = np.deg2rad(15 * (hour - 12))
    for delta_deg in np.linspace(-23.44, 23.44, 361):
        delta = np.deg2rad(delta_deg)
        S = sun_vector(H, delta)
        # Slunce musí být nad horizontem.
        if S[2] <= 0:
            continue
        if shadow_point(H, delta) is not None:
            return True
    return False

for hour in range(1, 24):
    H = np.deg2rad(15 * (hour - 12))

    # Hodinová rovina je určena osou Země (gnómonem) a směrem Slunce
    # pro daný hodinový úhel. Pro výpočet normály hodinové roviny stačí
    # použít libovolnou deklinaci; zde volíme rovnodennost (delta = 0).
    s_eq = sun_vector(H, 0.0)
    hour_plane_normal = np.cross(G, s_eq)

    # Směr hodinové čáry na stěně je průsečnice stěny a hodinové roviny.
    d = np.cross(n, hour_plane_normal)
    p = to2D(d)
    p /= np.linalg.norm(p)

    # Hodinové čáry kreslíme jako polopřímky od paty gnómonu:
    # I–XI vlevo (x < 0), XII dolů (y < 0), XIII–XXIII vpravo (x > 0).
    extent = 3.0
    if hour <= 11:
        sign = -1.0 if p[0] > 0 else 1.0
    elif hour == 12:
        sign = -1.0 if p[1] > 0 else 1.0
    else:
        sign = 1.0 if p[0] > 0 else -1.0

    ray_dir = sign * p
    # Vykreslujeme jen čáry směřující do dolní poloviny ciferníku (y < 0)
    # a jen ty, které mohou být během roku reálně osvětlené.
    if ray_dir[1] >= 0:
        continue
    if not is_hour_visible(hour):
        continue

    t_limits = []
    if abs(ray_dir[0]) > 1e-12:
        t_limits.append(extent / abs(ray_dir[0]))
    if abs(ray_dir[1]) > 1e-12:
        t_limits.append(extent / abs(ray_dir[1]))
    if not t_limits:
        continue

    t_end = min(t_limits)
    endpoint = ray_dir * t_end

    ax.plot([0.0, endpoint[0]], [0.0, endpoint[1]], 'k-')

    label_pos = endpoint * 1.05
    ax.text(label_pos[0], label_pos[1], to_roman(hour),
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

ax.set_title("Stěna (pohled)")
ax.set_aspect('equal')
ax.set_xlim(-4, 4)
ax.set_ylim(-4, 4)
ax.axhline(0, linewidth=0.5)
ax.axvline(0, linewidth=0.5)
ax.legend()

# =========================
# 2️⃣ PŮDORYS (TOP VIEW)
# =========================
ax = axs[1]

# projekce do XY roviny
Gn = np.array([G[0], -G[1]])
nn = np.array([-n[0], n[1]])

# normála stěny
ax.plot([0, nn[0]], [0, nn[1]], 'k--', label="normála stěny")

# gnómon
ax.plot([0, Gn[0]], [0, Gn[1]], 'r-', label="gnómon")

# stěna = kolmice na normálu
wall_dir = np.array([-nn[1], nn[0]])
ax.plot([-wall_dir[0], wall_dir[0]],
        [-wall_dir[1], wall_dir[1]], 'k-', label="stěna")

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
ax.plot([0,0], [0,1], 'k--')

# gnómon
ax.plot([0, -g_normal], [0, -g_height], 'r-', label="gnómon")

# stěna
ax.plot([0,0], [-1,1], 'k-', linewidth=2, label="stěna")

ax.set_title("Bokorys")
ax.set_aspect('equal')
ax.legend()


plt.tight_layout()
plt.show()
