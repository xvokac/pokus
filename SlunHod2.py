import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# =========================
# PARAMETRY
# =========================


def parse_args(argv):
    default_lat = 50.0
    default_azimuth = 0.0
    default_gnomon_length = 1.0
    export_dxf_path = None

    args = argv[1:]
    if "--export-dxf" in args:
        idx = args.index("--export-dxf")
        if idx + 1 >= len(args):
            print("Chybí cesta k výstupnímu DXF souboru za --export-dxf.")
        else:
            export_dxf_path = args[idx + 1]
            del args[idx:idx + 2]

    if len(args) in (2, 3):
        try:
            lat = float(args[0])
            azimuth = float(args[1])
            gnomon_length = float(args[2]) if len(args) == 3 else default_gnomon_length
            return lat, azimuth, gnomon_length, export_dxf_path
        except ValueError:
            print("Neplatné argumenty – očekávám čísla.")
    elif len(args) != 0:
        print("Špatný počet argumentů.")

    print(
        "Použití: python SlunHod2.py <zeměpisná_šířka_N> <azimut_stěny> "
        "[délka_gnómu] [--export-dxf output_file.dxf]"
    )
    print("  <zeměpisná_šířka_N>: severní zeměpisná šířka ve stupních (např. 49.5)")
    print("  <azimut_stěny>: azimut stěny ve stupních, 0 = jih,")
    print("                  kladné hodnoty směrem k západu, záporné k východu")
    print("  [délka_gnómu]: délka gnómu (volitelné, implicitně 1)")
    print("  [--export-dxf output_file.dxf]: export pohledu na stěnu do DXF")
    print(
        "Spouštím s implicitními hodnotami: "
        f"{default_lat}, {default_azimuth} a {default_gnomon_length}.\n"
    )

    return default_lat, default_azimuth, default_gnomon_length, export_dxf_path


lat_deg, azimuth_deg, gnomon_length, export_dxf_path = parse_args(sys.argv)

phi = np.deg2rad(lat_deg)   # zeměpisná šířka
D   = np.deg2rad(azimuth_deg)   # azimut stěny (0 = jih)

L = gnomon_length  # délka gnómonu


def wall_equator_horizon_angle(phi_rad, azimuth_rad):
    """Úhel mezi rovníkem a horizontem přímo na stěně.

    Geometricky jde o úhel mezi průsečnicí stěny s rovinou rovníku
    a průsečnicí stěny s horizontální rovinou.
    """
    sin_phi = np.sin(phi_rad)
    cos_phi = np.cos(phi_rad)
    sin_D = np.sin(azimuth_rad)

    numerator = abs(sin_D * cos_phi)
    denominator = abs(sin_phi)
    return np.arctan2(numerator, denominator)


angle_exact_deg = np.rad2deg(wall_equator_horizon_angle(phi, D))
approx_coeff = np.inf if np.isclose(np.tan(phi), 0.0) else 1.0 / np.tan(phi)
angle_approx_deg = abs(azimuth_deg) * approx_coeff

print("Kontrola vztahu pro nenulový azimut:")
print(f"  Zeměpisná šířka φ = {lat_deg:.3f}°")
print(f"  Azimut stěny D = {azimuth_deg:.3f}°")
print(f"  Přesný úhel (rovník–horizont na stěně) = {angle_exact_deg:.3f}°")
if np.isfinite(angle_approx_deg):
    print(
        "  Maloúhlová aproximace: "
        f"α ≈ |D|·cot(φ) = {approx_coeff:.3f}·|D| = {angle_approx_deg:.3f}°"
    )
else:
    print("  Maloúhlová aproximace α ≈ |D|·cot(φ) zde není použitelná (cot(φ) → ∞).")
print()

# =========================
# ZÁKLADNÍ VEKTORY
# =========================
n = np.array([np.sin(D), -np.cos(D), 0])   # normála stěny
G = L * np.array([0, np.cos(phi), np.sin(phi)])  # gnómon

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


def write_dxf_wall_view(path, line_entities, text_entities):
    """Minimalistický ASCII DXF export (R12) pro 2D čáry a text."""
    lines = [
        "0", "SECTION",
        "2", "HEADER",
        "0", "ENDSEC",
        "0", "SECTION",
        "2", "TABLES",
        "0", "TABLE",
        "2", "LAYER",
        "70", "1",
        "0", "LAYER",
        "2", "0",
        "70", "0",
        "62", "7",
        "6", "CONTINUOUS",
        "0", "ENDTAB",
        "0", "ENDSEC",
        "0", "SECTION",
        "2", "ENTITIES",
    ]

    for x1, y1, x2, y2 in line_entities:
        lines.extend([
            "0", "LINE",
            "8", "0",
            "10", f"{x1:.10f}",
            "20", f"{y1:.10f}",
            "30", "0.0",
            "11", f"{x2:.10f}",
            "21", f"{y2:.10f}",
            "31", "0.0",
        ])

    for x, y, height, value in text_entities:
        lines.extend([
            "0", "TEXT",
            "8", "0",
            "10", f"{x:.10f}",
            "20", f"{y:.10f}",
            "30", "0.0",
            "40", f"{height:.10f}",
            "1", value,
        ])

    lines.extend(["0", "ENDSEC", "0", "EOF"])
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")

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


def solar_declination_and_eot(day_of_year):
    """Vrátí (deklinace, časová rovnice) pro daný den v roce.

    day_of_year: 1..365
    deklinace v radiánech, časová rovnice v minutách.
    Použitá je standardní aproximace NOAA.
    """
    gamma = 2.0 * np.pi * (day_of_year - 1) / 365.0

    delta = (
        0.006918
        - 0.399912 * np.cos(gamma)
        + 0.070257 * np.sin(gamma)
        - 0.006758 * np.cos(2.0 * gamma)
        + 0.000907 * np.sin(2.0 * gamma)
        - 0.002697 * np.cos(3.0 * gamma)
        + 0.00148 * np.sin(3.0 * gamma)
    )

    eot_minutes = 229.18 * (
        0.000075
        + 0.001868 * np.cos(gamma)
        - 0.032077 * np.sin(gamma)
        - 0.014615 * np.cos(2.0 * gamma)
        - 0.040849 * np.sin(2.0 * gamma)
    )
    return delta, eot_minutes


def analemma_color(day_of_year):
    """Barva osmičky časové rovnice podle ročního období."""
    # Přibližné dny v roce (nepřestupný rok):
    # jarní rovnodennost ≈ 79, letní slunovrat ≈ 172,
    # podzimní rovnodennost ≈ 266, zimní slunovrat ≈ 355.
    if day_of_year >= 355 or day_of_year <= 79:
        return "tab:blue"    # zima: zimní slunovrat -> jarní rovnodennost
    if 172 <= day_of_year <= 266:
        return "tab:red"     # léto: letní slunovrat -> podzimní rovnodennost
    return "gold"            # přechodné části: neutrální žlutá



# =========================
# FIGURY
# =========================
fig_wall, ax_wall = plt.subplots(figsize=(6, 6))
fig_top, ax_top = plt.subplots(figsize=(6, 6))
fig_side, ax_side = plt.subplots(figsize=(6, 6))

# =========================
# 1️⃣ SLUNEČNÍ HODINY
# =========================
ax = ax_wall
displayed_hour_angles_deg = []
dxf_lines = []
dxf_texts = []

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
    """Vrací True, pokud se hodinová čára může během roku objevit na stěně.

    Poznámka:
    Pro VI/XVIII na jižní stěně (D = 0) je stín při rovnodennosti v nekonečnu
    (paprsek je rovnoběžný se stěnou). V grafu je ale chceme ponechat jako
    hraniční hodinové čáry, proto sledujeme i tento limitní případ.
    """
    H = np.deg2rad(15 * (hour - 12))
    has_grazing_case = False
    for delta_deg in np.linspace(-23.44, 23.44, 361):
        delta = np.deg2rad(delta_deg)
        S = sun_vector(H, delta)

        # Slunce musí být alespoň na horizontu (pro limitní případ).
        if S[2] < 0:
            continue

        denom = np.dot(n, S)
        if abs(denom) < 1e-12:
            has_grazing_case = True
            continue

        # Průsečík je v přední polorovině stěny.
        t = -np.dot(n, G) / denom
        if t > 0:
            return True

    # Limitní viditelnost: typicky VI/XVIII pro D≈0.
    if has_grazing_case and np.isclose(np.sin(D), 0.0, atol=1e-12):
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
    # Vykreslujeme jen čáry směřující převážně do dolní poloviny ciferníku.
    # Malá tolerance (y < 0.1) pomůže zobrazit i téměř vodorovné čáry
    # kolem VI/XVIII pro azimut 0.
    if ray_dir[1] >= 0.1:
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
    dxf_lines.append((0.0, 0.0, endpoint[0], endpoint[1]))
    displayed_hour_angles_deg.append(np.rad2deg(H))

    label_pos = endpoint * 1.05
    ax.text(label_pos[0], label_pos[1], to_roman(hour),
            ha='center', va='center', fontsize=9)
    dxf_texts.append((label_pos[0], label_pos[1], 0.14, to_roman(hour)))


# Rozsah hodin pro doplňkové křivky odvodíme z hodin skutečně zobrazených
# v ciferníku (ne pevně VI–XVIII).
if displayed_hour_angles_deg:
    min_hour_angle_deg = min(displayed_hour_angles_deg)
    max_hour_angle_deg = max(displayed_hour_angles_deg)
else:
    min_hour_angle_deg = -90.0
    max_hour_angle_deg = 90.0

# křivky (pro fixní deklinaci δ)
all_front_points = []
for d, name in [(-23.44, "Zima"),
                (0.0, "Rovnodennost"),
                (23.44, "Léto")]:
    pts = []
    for h in np.linspace(min_hour_angle_deg, max_hour_angle_deg, 400):
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
        for i in range(len(pts) - 1):
            dxf_lines.append((pts[i, 0], pts[i, 1], pts[i + 1, 0], pts[i + 1, 1]))
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
        H_deg = np.rad2deg(H)
        if H_deg < min_hour_angle_deg or H_deg > max_hour_angle_deg:
            continue
        P = shadow_point(H, delta)
        if P is None:
            continue
        horizon_pts.append(to2D(P))

if len(horizon_pts) > 2:
    horizon_pts = np.array(horizon_pts)
    ax.plot(horizon_pts[:,0], horizon_pts[:,1], label="Horizont")
    for i in range(len(horizon_pts) - 1):
        dxf_lines.append((
            horizon_pts[i, 0], horizon_pts[i, 1],
            horizon_pts[i + 1, 0], horizon_pts[i + 1, 1]
        ))
    all_front_points.append(horizon_pts)

# osmička časové rovnice (analema) na hodinové čáře XII
analemma_day_pts = []
for day in range(1, 366):
    delta, eot_min = solar_declination_and_eot(day)
    # V 12:00 středního slunečního času je hodinový úhel
    # roven odchylce časové rovnice.
    H = np.deg2rad(0.25 * eot_min)
    S = sun_vector(H, delta)
    if S[2] <= 0:
        continue
    P = shadow_point(H, delta)
    if P is None:
        continue
    analemma_day_pts.append((day, to2D(P)))

if len(analemma_day_pts) > 2:
    segment_points = []
    for i in range(len(analemma_day_pts) - 1):
        day_a, p_a = analemma_day_pts[i]
        day_b, p_b = analemma_day_pts[i + 1]
        if day_b - day_a != 1:
            continue
        color = analemma_color(day_a)
        ax.plot(
            [p_a[0], p_b[0]],
            [p_a[1], p_b[1]],
            color=color,
            linewidth=1.8,
        )
        dxf_lines.append((p_a[0], p_a[1], p_b[0], p_b[1]))
        segment_points.extend([p_a, p_b])

    if segment_points:
        all_front_points.append(np.array(segment_points))

    # položky legendy
    ax.plot([], [], color="tab:blue", linewidth=1.8, label="Čas. rovnice: zima")
    ax.plot([], [], color="tab:red", linewidth=1.8, label="Čas. rovnice: léto")
    ax.plot([], [], color="gold", linewidth=1.8, label="Čas. rovnice: jaro/podzim")

# gnómon
g2 = to2D(G)
ax.plot([0, g2[0]], [0, g2[1]], 'r-', linewidth=2)
dxf_lines.append((0.0, 0.0, g2[0], g2[1]))

ax.set_title("Stěna (pohled)")
ax.set_aspect('equal')
ax.set_xlim(-4, 4)
ax.set_ylim(-4, 4)
ax.axhline(0, linewidth=0.5)
ax.axvline(0, linewidth=0.5)
ax.legend()

if export_dxf_path:
    if dxf_lines:
        min_x = min(min(x1, x2) for x1, _, x2, _ in dxf_lines)
        max_x = max(max(x1, x2) for x1, _, x2, _ in dxf_lines)
        min_y = min(min(y1, y2) for _, y1, _, y2 in dxf_lines)
    else:
        min_x, max_x, min_y = -1.0, 1.0, -1.0
    width = max(1.0, max_x - min_x)
    base_x = min_x
    base_y = min_y - 0.35
    line_gap = 0.22

    metadata = [
        f"Vstup: zemepisna sirka = {lat_deg:.6f} deg",
        f"Vstup: azimut steny = {azimuth_deg:.6f} deg",
        f"Vstup: delka gnomu = {gnomon_length:.6f}",
        (
            "Vrchol gnomu vzhledem k pate (x,y,z) = "
            f"({G[0]:.6f}, {G[1]:.6f}, {G[2]:.6f})"
        ),
    ]
    for i, text in enumerate(metadata):
        dxf_texts.append((base_x, base_y - i * line_gap, max(0.10, 0.03 * width), text))

    write_dxf_wall_view(export_dxf_path, dxf_lines, dxf_texts)
    print(f"DXF export uložen: {export_dxf_path}")

# =========================
# 2️⃣ PŮDORYS (TOP VIEW)
# =========================
ax = ax_top

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

ax.set_title("Půdorys (gnómon směřuje k jihu)")
ax.set_aspect('equal')
ax.legend()


# =========================
# 3️⃣ BOKORYS (SIDE VIEW)
# =========================
ax = ax_side

def to_side(P):
    """2D souřadnice v bokorysu: x = vzdálenost od stěny, y = výška."""
    return np.array([-np.dot(P, n), -P[2]])

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

# body na polední (XII) poloze
# H = 0 odpovídá slunečnímu poledni
delta_horizon = phi - np.pi / 2
side_markers = [
    ("Horizont", delta_horizon, "tab:gray"),
    ("Rovník", 0.0, "tab:blue"),
    ("Zima", np.deg2rad(-23.44), "tab:cyan"),
    ("Léto", np.deg2rad(23.44), "tab:orange"),
]

for label, delta, color in side_markers:
    P = shadow_point(0.0, delta)
    if P is None:
        continue
    P2 = to_side(P)
    ax.plot(P2[0], P2[1], marker='o', color=color, markersize=5)
    ax.text(P2[0] + 0.05, P2[1], label, color=color, fontsize=9, va='center')

ax.set_title("Bokorys")
ax.set_aspect('equal')
ax.legend()


for fig in (fig_wall, fig_top, fig_side):
    fig.tight_layout()
plt.show()
