"""
1D gas-side heat flux along chamber -> convergent -> throat -> bell nozzle
using the Bartz correlation with CEA-derived combustion properties.

Inputs below are taken straight from the design sheet (5 kN, 25 bar, O/F 1.331).
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.optimize import brentq

# ------------------------------------------------------------------ CEA inputs
Pc      = 25.0e5      # Pa, chamber (stagnation) pressure
Tc      = 2823.7      # K
gamma   = 1.2245
cp_g    = 2280.9      # J/kg-K
Pr_g    = 0.6111
mu_g    = 8.979e-5    # Pa-s
cstar_th = 1676.3     # m/s
eta_cs  = 0.92
cstar   = cstar_th * eta_cs   # delivered c*, consistent with mdot = Pc*At/c*

# ------------------------------------------------------------------ geometry (m)
D_t     = 0.042077          # throat diameter from sheet
D_c     = 3.0 * 0.0254      # lined chamber ID (3.25" tube - 2 x 0.125" ablative)
L_cyl   = 0.304             # cylindrical section length (12" tube)
theta_c = np.radians(30.0)  # convergent half-angle
R_u     = 1.5 * D_t / 2     # upstream throat arc radius (Bartz "R" uses this)
R_d     = 0.382 * D_t / 2   # downstream throat arc radius
eps     = 4.067             # exit area ratio from sheet
# Simplified Rao bell (parabolic approximation): arc R_d to theta_n, then a quadratic
# Bezier from theta_n to theta_e, length = bell_frac x 15deg-cone length.
bell_frac = 0.80
theta_n   = np.radians(26.0)   # initial parabola angle (Rao chart, eps~4, 80% bell)
theta_e   = np.radians(11.0)   # exit angle

# ------------------------------------------------------------------ wall temps (K)
# Bartz needs a wall temperature. Ablative surface runs near its pyrolysis temp;
# copper insert runs much colder early in the burn. Rough but defensible.
T_w_ablative = 1000.0
T_w_copper   = 600.0
T_w_steel    = 1000.0                 # 4130 expansion section
# Copper insert spans the whole convergent (from the end of the barrel) through the
# throat to x_copper_end downstream. Everything past that is 4130.
x_copper_end = 0.010                  # m downstream of throat plane

# ------------------------------------------------------------------ burn / ablative model
t_burn      = 5.0      # s
n_t         = 101
mult_start  = 0.12     # fraction of Bartz q reaching the wall through fresh ablative
mult_end    = 1.00     # fraction at burnout (liner consumed)
# Ablative multiplier applies ONLY to the cylindrical (lined) barrel. Convergent,
# throat (copper) and divergent (steel) get multiplier 1 at all times.

# ------------------------------------------------------------------ build contour
def build_contour(n=800):
    R_t = D_t / 2
    R_c = D_c / 2
    R_e = R_t * np.sqrt(eps)
    # upstream arc: from angle -theta_c to 0 on a circle of radius R_u centred at (0, R_t + R_u)
    x_arc_start = -R_u * np.sin(theta_c)
    r_arc_start = R_t + R_u * (1 - np.cos(theta_c))
    L_cone = (R_c - r_arc_start) / np.tan(theta_c)
    x_cone_start = x_arc_start - L_cone
    x_cyl_start = x_cone_start - L_cyl

    xs, rs = [], []
    x = np.linspace(x_cyl_start, x_cone_start, n // 4, endpoint=False)
    xs.append(x); rs.append(np.full_like(x, R_c))
    x = np.linspace(x_cone_start, x_arc_start, n // 4, endpoint=False)
    xs.append(x); rs.append(R_c - (x - x_cone_start) * np.tan(theta_c))
    th = np.linspace(-theta_c, 0.0, n // 4, endpoint=False)
    xs.append(R_u * np.sin(th)); rs.append(R_t + R_u * (1 - np.cos(th)))

    # --- bell: downstream arc to theta_n
    th = np.linspace(0.0, theta_n, n // 8, endpoint=False)
    xs.append(R_d * np.sin(th)); rs.append(R_t + R_d * (1 - np.cos(th)))
    # --- parabolic section as quadratic Bezier N -> Q -> E
    N = np.array([R_d * np.sin(theta_n), R_t + R_d * (1 - np.cos(theta_n))])
    L_bell = bell_frac * (R_e - R_t) / np.tan(np.radians(15.0))   # Rao length def.
    E = np.array([L_bell, R_e])
    # control point = intersection of tangent lines from N (slope tan theta_n) and E (slope tan theta_e)
    m1, m2 = np.tan(theta_n), np.tan(theta_e)
    xq = (E[1] - N[1] + m1 * N[0] - m2 * E[0]) / (m1 - m2)
    Q = np.array([xq, N[1] + m1 * (xq - N[0])])
    s = np.linspace(0.0, 1.0, n // 8)[:, None]
    B = (1 - s) ** 2 * N + 2 * (1 - s) * s * Q + s ** 2 * E
    xs.append(B[:, 0]); rs.append(B[:, 1])
    return np.concatenate(xs), np.concatenate(rs), x_cone_start

# ------------------------------------------------------------------ Mach from area ratio
def area_ratio(M):
    g = gamma
    return (1.0 / M) * ((2.0 / (g + 1)) * (1 + (g - 1) / 2 * M**2)) ** ((g + 1) / (2 * (g - 1)))

def mach_from_area(A_At, supersonic):
    if A_At <= 1.0 + 1e-9:
        return 1.0
    f = lambda M: area_ratio(M) - A_At
    return brentq(f, 1.0 + 1e-6, 20.0) if supersonic else brentq(f, 1e-6, 1.0 - 1e-6)

# ------------------------------------------------------------------ Bartz
def bartz(x, r):
    A_t = np.pi * (D_t / 2) ** 2
    A = np.pi * r**2
    M = np.array([mach_from_area(a / A_t, xi > 0) for a, xi in zip(A, x)])

    T_w = np.where(x >= x_barrel_end, T_w_copper, T_w_ablative)   # convergent = copper
    T_w = np.where(x > x_copper_end, T_w_steel, T_w)               # past insert = 4130

    g = gamma
    fac = 1 + (g - 1) / 2 * M**2
    sigma = 1.0 / ((0.5 * T_w / Tc * fac + 0.5) ** 0.68 * fac**0.12)

    h_g = (0.026 / D_t**0.2) * (mu_g**0.2 * cp_g / Pr_g**0.6) \
          * (Pc / cstar) ** 0.8 * (D_t / R_u) ** 0.1 * (A_t / A) ** 0.9 * sigma

    r_rec = Pr_g ** (1.0 / 3.0)                       # turbulent recovery factor
    T_aw = Tc * (1 + r_rec * (g - 1) / 2 * M**2) / fac
    q = h_g * (T_aw - T_w)
    return M, h_g, T_aw, T_w, q

# ------------------------------------------------------------------ run
x, r, x_barrel_end = build_contour()
M, h_g, T_aw, T_w, q_bartz = bartz(x, r)          # steady Bartz flux, no ablative credit

def ablative_multiplier(x, t):
    """Linear ramp mult_start -> mult_end over t_burn on the barrel; 1 elsewhere."""
    f_t = mult_start + (mult_end - mult_start) * np.clip(t / t_burn, 0.0, 1.0)
    return np.where(x < x_barrel_end, f_t, 1.0)

t = np.linspace(0.0, t_burn, n_t)
Q = np.array([q_bartz * ablative_multiplier(x, ti) for ti in t])   # [n_t, n_x]  W/m2
E_wall = np.trapezoid(Q, t, axis=0)                                # J/m2 delivered over burn

i_t = np.argmin(np.abs(x))
i_b = np.argmax(x >= x_barrel_end) - 1
print(f"Bell: length {x.max()*1e3:.1f} mm, exit D {2*r[-1]*1e3:.1f} mm")
print(f"Throat (no ablative):  q = {q_bartz[i_t]/1e6:.2f} MW/m2 constant, "
      f"E = {E_wall[i_t]/1e6:.1f} MJ/m2 over {t_burn:.0f} s")
print(f"Barrel: Bartz q = {q_bartz[i_b]/1e6:.2f} MW/m2 -> wall q "
      f"{Q[0,i_b]/1e6:.2f} (t=0) .. {Q[-1,i_b]/1e6:.2f} MW/m2 (t={t_burn:.0f} s), "
      f"E = {E_wall[i_b]/1e6:.1f} MJ/m2 (mean mult {(mult_start+mult_end)/2:.2f})")

# ------------------------------------------------------------------ plot helpers
def plot_contour(ax, fig):
    q_end = Q[-1]
    pts = np.array([x, r]).T.reshape(-1, 1, 2)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    norm = plt.Normalize(Q.min() / 1e6, Q.max() / 1e6)
    for sign in (1, -1):
        s = segs.copy(); s[:, :, 1] *= sign
        lc = LineCollection(s, cmap="inferno", norm=norm, linewidth=6)
        lc.set_array(q_end[:-1] / 1e6)
        ax.add_collection(lc)
    ax.fill_between(x, -r, r, color="0.92", zorder=0)
    ax.axvspan(x.min(), x_barrel_end, color="tab:brown", alpha=0.08, lw=0)
    ax.axvline(0, color="k", ls="--", lw=0.8)
    ax.set_xlim(x.min(), x.max()); ax.set_ylim(-1.15 * r.max(), 1.15 * r.max())
    ax.set_ylabel("radius [m]"); ax.set_xlabel("axial position from throat [m]")
    ax.set_title(f"Wall heat flux at burnout (t = {t_burn:.0f} s)  |  shaded = ablative-lined barrel")
    fig.colorbar(lc, ax=ax, pad=0.01).set_label("q [MW/m²]")

def plot_xt(ax, fig):
    X, T = np.meshgrid(x, t)
    cf = ax.contourf(X, T, Q / 1e6, levels=40, cmap="inferno")
    ax.contour(X, T, Q / 1e6, levels=[2, 4, 6, 10, 15, 20, 25], colors="w", linewidths=0.5, alpha=0.6)
    ax.axvline(x_barrel_end, color="c", ls=":", lw=1); ax.axvline(0, color="w", ls="--", lw=0.8)
    ax.set_ylabel("time [s]"); ax.set_xlabel("axial position from throat [m]")
    ax.set_xlim(x.min(), x.max())
    ax.set_title("Wall heat flux q(x, t) — ablative multiplier ramps "
                 f"{mult_start:.2f} → {mult_end:.2f} on barrel only")
    fig.colorbar(cf, ax=ax, pad=0.01).set_label("q [MW/m²]")
    return cf

def plot_profiles(ax, fig, align_cbar=None):
    ax.plot(x, q_bartz / 1e6, color="0.5", lw=1.2, ls="--", label="Bartz (no ablative)")
    ax.plot(x, Q[0] / 1e6, color="C0", lw=1.8, label="wall q, t = 0")
    ax.plot(x, Q[-1] / 1e6, color="C3", lw=1.8, label=f"wall q, t = {t_burn:.0f} s")
    ax.set_ylabel("q [MW/m²]"); ax.grid(alpha=0.3)
    axb = ax.twinx()
    axb.plot(x, E_wall / 1e6, color="C2", lw=1.5, ls=":", label="∫q dt over burn")
    axb.set_ylabel("E [MJ/m²]", color="C2")
    ax.axvline(0, color="k", ls="--", lw=0.8)
    ax.set_xlabel("axial position from throat [m]"); ax.set_xlim(x.min(), x.max())
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = axb.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper left")
    if align_cbar is not None:
        fig.colorbar(align_cbar, ax=ax, pad=0.01).ax.set_visible(False)

# ------------------------------------------------------------------ combined figure
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 11), layout="constrained",
                                    gridspec_kw={"height_ratios": [1.0, 1.3, 1.0]})
plot_contour(ax1, fig)
cf = plot_xt(ax2, fig)
plot_profiles(ax3, fig, align_cbar=cf)
fig.savefig("bartz_heat_flux.png", dpi=160)

# ------------------------------------------------------------------ individual figures
fig1, ax = plt.subplots(figsize=(12, 4), layout="constrained")
plot_contour(ax, fig1); fig1.savefig("bartz_contour.png", dpi=160)

fig2, ax = plt.subplots(figsize=(12, 5), layout="constrained")
plot_xt(ax, fig2); fig2.savefig("bartz_xt.png", dpi=160)

fig3, ax = plt.subplots(figsize=(12, 4.5), layout="constrained")
plot_profiles(ax, fig3); fig3.savefig("bartz_profiles.png", dpi=160)

plt.show()