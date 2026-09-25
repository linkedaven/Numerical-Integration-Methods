from functools import lru_cache

import matplotlib.pyplot as plt
import mpmath as mp
from matplotlib.ticker import ScalarFormatter, NullFormatter
from matplotlib.widgets import Slider, Button

mp.mp.dps = 60  # working precision, in decimal digits

# ----------------------------------------------------------------------
# Dark theme (same palette as the double pendulum / projectile motion apps)
# ----------------------------------------------------------------------
plt.style.use('dark_background')

BG_COLOR = '#1e1e1e'
PANEL_COLOR = '#2b2b2b'
TEXT_COLOR = '#e0e0e0'
GRID_COLOR = '#444444'
BUTTON_COLOR = '#3a3a3a'
BUTTON_HOVER = '#505050'
METHOD_COLORS = {
    'Trapezoidal':   '#ff9f43',
    'Simpson 1/3':   '#54a0ff',
    'Simpson 3/8':   '#1dd1a1',
    'Gaussian quad': '#feca57',
}

# ----------------------------------------------------------------------
# Change the function here to your choice
# ----------------------------------------------------------------------
def f(x):
    return mp.sin(x)

a, b = mp.mpf(0), mp.pi
I_exact = mp.mpf(2)

# ========================================================================
# 1) TRAPEZOIDAL RULE  (composite, degree-1 Newton-Cotes)
# ========================================================================
def trapezoidal_rule(f, a, b, n):
    h = (b - a) / n
    total = (f(a) + f(b)) / 2
    for i in range(1, n):
        total += f(a + i * h)
    return h * total

# ========================================================================
# 2) SIMPSON'S 1/3 RULE  (composite, degree-2 Newton-Cotes, n must be even)
# ========================================================================
def simpsons_1_3(f, a, b, n):
    if n % 2 != 0:
        n += 1
    h = (b - a) / n
    total = f(a) + f(b)
    for i in range(1, n):
        total += (4 if i % 2 else 2) * f(a + i * h)
    return h / 3 * total

# ========================================================================
# 3) SIMPSON'S 3/8 RULE  (composite, degree-3 Newton-Cotes, n multiple of 3)
# ========================================================================
def simpsons_3_8(f, a, b, n):
    if n % 3 != 0:
        n += 3 - (n % 3)
    h = (b - a) / n
    total = f(a) + f(b)
    for i in range(1, n):
        total += (2 if i % 3 == 0 else 3) * f(a + i * h)
    return 3 * h / 8 * total

# ========================================================================
# 4) GAUSSIAN QUADRATURE  (Gauss-Legendre, n nodes)
# ========================================================================
@lru_cache(maxsize=None)
def gauss_legendre_nodes_weights(n):
    J = mp.zeros(n, n)
    for k in range(1, n):
        beta = mp.mpf(k) / mp.sqrt(4 * mp.mpf(k) ** 2 - 1)
        J[k - 1, k] = beta
        J[k, k - 1] = beta
    E, EV = mp.eigsy(J)
    nodes = tuple(E[i] for i in range(n))
    weights = tuple(2 * EV[0, i] ** 2 for i in range(n))
    return nodes, weights

def gaussian_quadrature(f, a, b, n):
    nodes, weights = gauss_legendre_nodes_weights(n)
    total = mp.mpf(0)
    for x, w in zip(nodes, weights):
        xi = 0.5 * (b - a) * x + 0.5 * (b + a)
        total += w * f(xi)
    return 0.5 * (b - a) * total

# ========================================================================
# Figure & widgets
# ========================================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6.3),
                                gridspec_kw={'width_ratios': [1.3, 1]})
fig.patch.set_facecolor(BG_COLOR)
plt.subplots_adjust(bottom=0.22)

results_table = None  # created/replaced each time Run is clicked

def compute_results(n_points):
    """Single-n comparison, used for the results table."""
    return {
        'Trapezoidal':                 trapezoidal_rule(f, a, b, n_points),
        "Simpson's 1/3":               simpsons_1_3(f, a, b, n_points),
        "Simpson's 3/8":               simpsons_3_8(f, a, b, n_points),
        'Gaussian quadrature':         gaussian_quadrature(f, a, b, n_points),
    }

methods_conv = {
    'Trapezoidal':   lambda n: trapezoidal_rule(f, a, b, n),
    "Simpson 1/3":   lambda n: simpsons_1_3(f, a, b, n),
    "Simpson 3/8":   lambda n: simpsons_3_8(f, a, b, n),
    'Gaussian quad': lambda n: gaussian_quadrature(f, a, b, n),
}

def recompute_and_draw(n_points):
    global results_table

    results = compute_results(n_points)

    # ---- convergence plot: every n from 1 up to n_points ----
    ns = list(range(1, n_points + 1))
    ax1.clear()
    ax1.set_facecolor(BG_COLOR)
    for name, fn in methods_conv.items():
        errs = [float(abs(fn(n) - I_exact)) for n in ns]
        ax1.loglog(ns, errs, 'o-', color=METHOD_COLORS[name], label=name)
    ax1.set_xlabel('n (points / subintervals)', color=TEXT_COLOR)
    ax1.set_ylabel('|error| (log scale)', color=TEXT_COLOR)
    ax1.set_title(f'Convergence (mpmath, {mp.mp.dps}-digit precision):\n'
                  r'$\int_0^{\pi} \sin(x)\,dx = 2$', color=TEXT_COLOR)
    leg = ax1.legend(facecolor=PANEL_COLOR, edgecolor=GRID_COLOR)
    for txt in leg.get_texts():
        txt.set_color(TEXT_COLOR)
    ax1.grid(True, which='both', color=GRID_COLOR, alpha=0.4)
    ax1.tick_params(colors=TEXT_COLOR)
    for spine in ax1.spines.values():
        spine.set_color(GRID_COLOR)

    tick_step = max(1, len(ns) // 10)
    ax1.set_xticks(ns[::tick_step])
    ax1.xaxis.set_major_formatter(ScalarFormatter())
    ax1.xaxis.set_minor_formatter(NullFormatter())
    ax1.tick_params(axis='x', labelsize=8, colors=TEXT_COLOR)

    # ---- results table ----
    ax2.clear()
    ax2.set_facecolor(BG_COLOR)
    ax2.axis('off')
    ax2.set_title(f'Results at n = {n_points}', fontsize=11, color=TEXT_COLOR)
    rows = [["Method", "Approx.", "Error"]]
    for name, val in results.items():
        err = abs(val - I_exact)
        rows.append([name, mp.nstr(val, 10), mp.nstr(err, 3)])

    results_table = ax2.table(cellText=rows, cellLoc='center', loc='center',
                               colWidths=[0.5, 0.28, 0.22],
                               bbox=[0.0, 0.15, 1.0, 0.7])
    results_table.auto_set_font_size(False)
    results_table.set_fontsize(9)
    for (r, c), cell in results_table.get_celld().items():
        cell.set_edgecolor(GRID_COLOR)
        if r == 0:
            cell.set_facecolor(BUTTON_COLOR)
            cell.set_text_props(color=TEXT_COLOR, weight='bold')
        else:
            cell.set_facecolor(PANEL_COLOR)
            cell.set_text_props(color=TEXT_COLOR)

    fig.canvas.draw_idle()

# ---- initial draw ----
recompute_and_draw(10)

# ---- n slider + Run button, in the window itself ----
ax_slider = fig.add_axes([0.30, 0.07, 0.35, 0.05])
ax_slider.set_facecolor(PANEL_COLOR)
slider_n = Slider(ax_slider, 'n', valmin=1, valmax=60, valinit=10, valstep=1,
                   color='#00e5ff',   # cyan -- the filled/active part of the track
                   track_color=PANEL_COLOR,
                   initcolor=None)
slider_n.label.set_color(TEXT_COLOR)
slider_n.valtext.set_color(TEXT_COLOR)
for spine in ax_slider.spines.values():
    spine.set_color(GRID_COLOR)

ax_run = fig.add_axes([0.70, 0.06, 0.12, 0.06])
ax_run.set_facecolor(BUTTON_COLOR)
btn_run = Button(ax_run, 'Run', color=BUTTON_COLOR, hovercolor=BUTTON_HOVER)

btn_run.color = BUTTON_COLOR
btn_run.hovercolor = BUTTON_HOVER
btn_run.ax.set_facecolor(BUTTON_COLOR)
btn_run.label.set_color(TEXT_COLOR)
btn_run.label.set_fontsize(10)
for spine in ax_run.spines.values():
    spine.set_color(GRID_COLOR)

def on_run(event):
    try:
        n_points = int(round(slider_n.val))
        if n_points < 1:
            print("n must be at least 1.")
            return
        recompute_and_draw(n_points)
    except ValueError:
        print("Please enter a valid integer for n.")

btn_run.on_clicked(on_run)

plt.show()