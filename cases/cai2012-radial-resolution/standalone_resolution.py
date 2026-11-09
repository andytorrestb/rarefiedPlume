#!/usr/bin/env python3
"""Standalone cell-size and spherical continuum-breakdown estimate.

Copy this file anywhere; needs only numpy, scipy and matplotlib. Run:
    python standalone_resolution.py

Edit the inputs below. The physics is extracted/adapted from this repository's
plumetools/cai2012ext/interface.py (equivalent_source, area_mach,
supersonic_mach, mean_free_path, bird_parameter),
plumetools/markelov1999/sourceflow.py (angular normalization), and
plumetools/radialstudy/physics.py (predict_breakdown).

This retains the Cai gross outward-flux mapping, spherical isentropic
expansion and VHS convention. The model starts at an equivalent-source anchor,
not the physical exit disk. Its axis breakdown radius is interpreted as a
spherical surrogate; it is not a measured or angle-dependent plume surface.
At an actual crossing, lambda_b is recalculated from local n and T, and
Kn_breakdown_D_e = lambda_b / D_e uses the same length as the exit Kn.
Cell sizing is piecewise: use lambda_b where a breakdown crossing exists;
otherwise use lambda_e for cases already noncontinuum at the model anchor.
If the search ends before breakdown, inlet sizing is left unavailable.
Repository license: GPL-3.0-or-later. No project imports or input files needed.
"""
import csv
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq


# ------------------------------- edit inputs --------------------------------
D_e = 0.2                         # Physical exit diameter, m
Kn_e = np.array([0.0001, 0.001, 0.01, 0.1, 1, 10, 100])
T_e = 300.0                       # Exit temperature, K
S_e = 2.0                         # Exit speed ratio U_e / sqrt(2 R T_e)
gamma = 5.0 / 3.0                 # Monatomic argon
mass = 6.63e-26                   # Mass of one argon atom, kg
diameter = 4.17e-10                # VHS reference collision diameter, m
T_ref = 273.0                     # VHS reference temperature, K
omega = 0.74                     # VHS viscosity exponent
B_THRESHOLD = 0.05
BIRD_DEFINITION = 'general'       # Or 'radial_approx': Stewart Eq. (3)
MAX_RADIUS_OVER_D = 25.0          # Search bound, not a substitute for a root
# Cell-size band from the user's script; the study's stricter target is 0.8 lambda.
DX_MIN_OVER_LAMBDA = 0.33
DX_MAX_OVER_LAMBDA = 1.0
OUTPUT_DIR = Path('standalone-resolution-results')  # Relative to working folder
SHOW_PLOTS = True                 # False saves figures without opening windows
# ---------------------------------------------------------------------------


def area_mach(mach):
    """Spherical area ratio F(M), Stewart & Lumpkin Eq. (2)."""
    m = np.asarray(mach, dtype=float)
    return ((2/(gamma+1))*(1+(gamma-1)*m*m/2))**((gamma+1)/(2*(gamma-1))) / m


def supersonic_mach(area_ratio, lower):
    """Extracted log-space bisection on the supersonic branch."""
    target = np.asarray(area_ratio, dtype=float)
    low = np.full(target.shape, math.log(max(lower, 1.0)))
    high = np.full(target.shape, math.log(2*max(lower, 1.0)+1))
    for _ in range(200):
        short = area_mach(np.exp(high)) < target
        if not np.any(short):
            break
        high = np.where(short, high+math.log(2), high)
    else:
        raise ValueError('Cannot bracket supersonic Mach number')
    for _ in range(200):
        mid = (low+high)/2
        below = area_mach(np.exp(mid)) < target
        low, high = np.where(below, mid, low), np.where(below, high, mid)
        if np.all(high-low < 1e-15):
            break
    return np.exp((low+high)/2)


def mean_free_path(n, temperature):
    """Argon VHS mean free path, m; same convention as the study."""
    return ((np.asarray(temperature)/T_ref)**(omega-0.5)
            / (math.sqrt(2)*math.pi*diameter**2*np.asarray(n)))


def make_source(kn):
    """Minimal extraction of the existing exit-to-equivalent-source mapping."""
    if not np.isfinite(kn) or kn <= 0:
        raise ValueError('Exit Kn must be positive and finite')
    if (not np.isfinite([D_e, T_e, S_e, gamma, mass, diameter, T_ref, omega]).all()
            or min(D_e, T_e, S_e, mass, diameter, T_ref) <= 0 or gamma <= 1):
        raise ValueError('Invalid exit or gas inputs')
    R = 1.380649e-23/mass
    cp = gamma*R/(gamma-1)
    a = math.sqrt(2*R*T_e)
    U_e = S_e*a
    # Gross outward drifting-Maxwellian moments; copied from one_way_moments.
    H = math.erfc(-S_e)/2
    E = math.exp(-S_e*S_e)/(2*math.sqrt(math.pi))
    J1 = U_e*H+a*E
    J3 = (U_e**3+1.5*U_e*a*a)*H+a*(U_e*U_e+a*a)*E
    enthalpy = (J3+2*R*T_e*J1)/(2*J1)
    T_a = (enthalpy-U_e*U_e/2)/cp
    if T_a <= 0:
        raise ValueError('Equivalent-source anchor temperature must be positive')
    M_a = U_e/math.sqrt(gamma*R*T_a)
    if M_a <= 1:
        raise ValueError('Equivalent-source anchor must be supersonic')
    theta_L = math.pi/2*(math.sqrt((gamma+1)/(gamma-1))-1)
    exponent = (gamma+0.41)/(gamma-1)
    I = quad(lambda th: math.cos(math.pi*th/(2*theta_L))**exponent*math.sin(th),
             0, theta_L, limit=200)[0]
    area = math.pi*D_e**2/4
    r_a = math.sqrt(area*J1/(2*math.pi*I*U_e))
    n_e = (T_e/T_ref)**(omega-0.5)/(math.sqrt(2)*math.pi*diameter**2*kn*D_e)
    return dict(R=R, T_t=enthalpy/cp, M_a=M_a, r_a=r_a, n_e=n_e, U_e=U_e)


def state(radius, source):
    """Isentropic radial streamtube state on the plume axis."""
    r = np.asarray(radius, dtype=float)
    if np.any(~np.isfinite(r)) or np.any(r < source['r_a']*(1-1e-12)):
        raise ValueError('Radius must be at or downstream of the model anchor')
    M = supersonic_mach(area_mach(source['M_a'])*(r/source['r_a'])**2, source['M_a'])
    T = source['T_t']/(1+(gamma-1)*M*M/2)
    U = M*np.sqrt(gamma*source['R']*T)
    n = source['n_e']*source['U_e']*source['r_a']**2/(U*r*r)
    return dict(M=M, T=T, U=U, n=n, **{'lambda': mean_free_path(n, T)},
                c_bar=np.sqrt(8*source['R']*T/math.pi))


def bird(radius, source, definition=BIRD_DEFINITION):
    """Bird B = -lambda U dot grad(rho)/(rho c_bar)."""
    st = state(radius, source)
    simple = 2*st['lambda']*st['U']/(np.asarray(radius)*st['c_bar'])
    if definition == 'radial_approx':
        return simple  # Stewart Eq. (3), without the acceleration correction.
    if definition == 'general':
        return simple*st['M']**2/(st['M']**2-1)
    raise ValueError('Bird definition must be general or radial_approx')


def breakdown(kn, *, max_radius_over_D=MAX_RADIUS_OVER_D,
              threshold=B_THRESHOLD, definition=BIRD_DEFINITION):
    """Find the first outward B=threshold crossing, or report its absence."""
    source = make_source(kn)
    if not np.isfinite(threshold) or threshold <= 0:
        raise ValueError('Bird threshold must be positive and finite')
    outer = max_radius_over_D*D_e
    if not np.isfinite(outer) or outer <= source['r_a']:
        raise ValueError('Search bound must exceed the equivalent-source anchor')
    radius = np.geomspace(source['r_a'], outer, 257)
    values = bird(radius, source, definition)
    result = dict(Kn=float(kn), lambda_e_m=float(kn*D_e),
                  anchor_over_D=source['r_a']/D_e, B_anchor=float(values[0]),
                  r_m=None, r_over_D=None, lambda_breakdown_m=None,
                  Kn_breakdown_D_e=None)
    if values[0] >= threshold:
        return dict(result, status='ALREADY_NONCONTINUUM_AT_INLET')
    crossings = np.flatnonzero((values[:-1] < threshold) & (values[1:] >= threshold))
    if not len(crossings):
        return dict(result, status='NOT_REACHED_IN_DOMAIN')
    i = crossings[0]
    root = brentq(lambda r: float(bird(r, source, definition))-threshold,
                  radius[i], radius[i+1], xtol=1e-13)
    # Re-evaluate the local VHS MFP; do not reuse the prescribed exit MFP.
    lambda_b = float(state(root, source)['lambda'])
    return dict(result, status='CROSSING_IN_DOMAIN', r_m=root, r_over_D=root/D_e,
                lambda_breakdown_m=lambda_b, Kn_breakdown_D_e=lambda_b/D_e)


def main():
    from matplotlib import pyplot as plt
    output = Path(OUTPUT_DIR)
    output.mkdir(parents=True, exist_ok=True)
    results = [breakdown(kn) for kn in Kn_e]
    for row in results:
        if row['status'] == 'CROSSING_IN_DOMAIN':
            inlet_lambda, basis = row['lambda_breakdown_m'], 'breakdown'
        elif row['status'] == 'ALREADY_NONCONTINUUM_AT_INLET':
            inlet_lambda, basis = row['lambda_e_m'], 'exit'
        else:
            inlet_lambda, basis = None, 'unavailable'
        row.update(sizing_basis=basis, lambda_inlet_m=inlet_lambda,
                   dx_min_m=None if inlet_lambda is None else DX_MIN_OVER_LAMBDA*inlet_lambda,
                   dx_max_m=None if inlet_lambda is None else DX_MAX_OVER_LAMBDA*inlet_lambda,
                   dx_exit_min_m=DX_MIN_OVER_LAMBDA*row['lambda_e_m'],
                   dx_exit_max_m=DX_MAX_OVER_LAMBDA*row['lambda_e_m'])
        radius = '---' if row['r_over_D'] is None else f"{row['r_over_D']:.6g}"
        lambda_b = '---' if row['lambda_breakdown_m'] is None else f"{row['lambda_breakdown_m']:.6g}"
        kn_b = '---' if row['Kn_breakdown_D_e'] is None else f"{row['Kn_breakdown_D_e']:.6g}"
        band = 'unavailable' if inlet_lambda is None else f"{row['dx_min_m']:.6g} to {row['dx_max_m']:.6g} m"
        print(f"Kn_e={row['Kn']:8g}  r_breakdown/D_e={radius:>10}  "
              f"lambda_b(m)={lambda_b:>10}  Kn_b(D_e)={kn_b:>10}  {row['status']}")
        print(f"    Inlet sizing at {basis}: dx = {band}")
    with (output/'resolution.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)
    from matplotlib.lines import Line2D
    dx_min = np.array([np.nan if r['dx_min_m'] is None else r['dx_min_m'] for r in results])
    dx_max = np.array([np.nan if r['dx_max_m'] is None else r['dx_max_m'] for r in results])
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.loglog(Kn_e, dx_min, '-', color='tab:blue', label=f'dx_min = {DX_MIN_OVER_LAMBDA:g} lambda_inlet')
    ax.loglog(Kn_e, dx_max, '-', color='tab:orange', label=f'dx_max = {DX_MAX_OVER_LAMBDA:g} lambda_inlet')
    ax.fill_between(Kn_e, dx_min, dx_max, color='grey', alpha=.1)
    for basis, marker in [('breakdown', 'D'), ('exit', 'o')]:
        selected = np.array([r['sizing_basis'] == basis for r in results])
        ax.scatter(Kn_e[selected], dx_min[selected], marker=marker, color='tab:blue', zorder=3)
        ax.scatter(Kn_e[selected], dx_max[selected], marker=marker, color='tab:orange', zorder=3)
    ax.set(xlabel='Exit Knudsen number, Kn_e', ylabel='Inlet cell edge, dx (m)',
           title='Piecewise DSMC inlet cell-size requirement')
    ax.set_xticks(Kn_e, [f'{kn:g}' for kn in Kn_e])
    ax.minorticks_off()
    ax.grid(True, which='major', alpha=.25)
    handles, labels = ax.get_legend_handles_labels()
    handles += [Line2D([], [], color='grey', marker='D', linestyle='None'),
                Line2D([], [], color='grey', marker='o', linestyle='None')]
    labels += ['Continuum case: inlet at breakdown', 'Noncontinuum case: inlet at exit']
    ax.legend(handles, labels, fontsize=9)
    fig.subplots_adjust(left=.12, right=.97, top=.91, bottom=.25)
    fig.text(.12, .05,
             'Continuum: use local lambda at the predicted breakdown crossing.\n'
             'Already noncontinuum at model anchor: use prescribed exit lambda.\n'
             'Lines join sampled cases; they do not locate the regime transition.',
             fontsize=9, linespacing=1.5)
    fig.savefig(output/'cell_sizes.png', dpi=160)
    fig, ax = plt.subplots(figsize=(8, 5.5))
    # Zero is the requested plotting convention for an already-noncontinuum inlet.
    # A crossing beyond the search bound remains unknown, rather than zero.
    radii = [0.0 if r['status'] == 'ALREADY_NONCONTINUUM_AT_INLET'
             else np.nan if r['r_m'] is None else r['r_m'] for r in results]
    ax.semilogx(Kn_e, radii, 'o-', label=rf'Bird $B = {B_THRESHOLD:g}$ ({BIRD_DEFINITION})')
    at_inlet = np.array([r['status'] == 'ALREADY_NONCONTINUUM_AT_INLET' for r in results])
    ax.scatter(Kn_e[at_inlet], np.zeros(np.count_nonzero(at_inlet)), color='tab:orange',
               label='Noncontinuum at inlet: shown as 0', zorder=3)
    for row in results:
        if row['Kn_breakdown_D_e'] is not None:
            ax.annotate(rf"$\mathbf{{Kn}}_{{\mathbf{{b}}}} = \mathbf{{{row['Kn_breakdown_D_e']:.4g}}}$",
                        (row['Kn'], row['r_m']), xytext=(14, 10),
                        textcoords='offset points', fontsize=10, fontweight='bold',
                        ha='left', va='bottom')
    ax.set(xlabel=r'Exit Knudsen number, $\mathrm{Kn}_e$', ylabel=r'Breakdown radius, $r_b$ (m)',
           title='Predicted spherical continuum breakdown')
    ax.set_xticks(Kn_e, [f'{kn:g}' for kn in Kn_e])
    ax.minorticks_off()
    ax.axhline(0, color='grey', linewidth=.8)
    ax.set_ylim(top=ax.get_ylim()[1]*1.15)  # Room above the highest point's label.
    ax.grid(True, which='major', alpha=.25)
    ax.legend()
    fig.subplots_adjust(left=.12, right=.97, top=.9, bottom=.29)
    already = [f"{r['Kn']:g}" for r in results if r['status'] == 'ALREADY_NONCONTINUUM_AT_INLET']
    missing = [f"{r['Kn']:g}" for r in results if r['status'] == 'NOT_REACHED_IN_DOMAIN']
    note = [f'Radius from equivalent-source centre, in metres; search ends at {MAX_RADIUS_OVER_D*D_e:g} m.']
    if already:
        note.append('Already noncontinuum at anchor: Kn = '+', '.join(already)+'.')
    if missing:
        note.append('No crossing inside search interval: Kn = '+', '.join(missing)+'.')
    note.append(r'$\mathrm{Kn}_b = \lambda_b / D_e$, using local VHS mean free path at the crossing.')
    note.append('Already noncontinuum at inlet is plotted as 0; an unreached crossing remains unavailable.')
    fig.text(.12, .045, '\n'.join(note), fontsize=9, linespacing=1.5)
    fig.savefig(output/'breakdown_radius.png', dpi=160)
    print(f'Figures and CSV: {output.resolve()}')
    if SHOW_PLOTS:
        plt.show()
    else:
        plt.close('all')


if __name__ == '__main__':
    main()
