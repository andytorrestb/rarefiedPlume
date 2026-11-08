"""Axis wrappers over the existing equivalent source and Bird implementation."""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy.optimize import brentq

from plumetools.cai2012ext.interface import equivalent_source
from .config import digest


@dataclass(frozen=True)
class RadialModel:
    config: object
    reference: object
    source: object
    r_in: float
    r_out: float
    physical_hash: str
    inlet_state: dict | None = None
    sampling_record: dict | None = None

    def state(self, r):
        if self.inlet_state is not None and np.ndim(r)==0 and float(r)==self.r_in:
            return dict(self.inlet_state)
        return self.source.state(r, 0.)

    def identity(self, *, angle_rad):
        return {'source_domain_hash': self.physical_hash, 'r_in_m': self.r_in,
                'r_out_m': self.r_out, 'angle_rad': angle_rad}


def source_for_kn(config, kn: float) -> RadialModel:
    if not np.isfinite(kn) or kn <= 0:
        raise ValueError('Kn must be positive and finite')
    ref = replace(config.reference, exit=replace(config.reference.exit, knudsen=float(kn)))
    src = equivalent_source(ref)
    inlet = src.r_a if config.domain.inlet_over_D is None else config.domain.inlet_over_D * src.diameter_m
    outer = config.domain.outer_over_D * src.diameter_m
    if inlet < src.r_a * (1 - 1e-12) or outer <= inlet:
        raise ValueError('domain requires anchor <= inlet < outlet')
    threshold = config.breakdown.threshold
    if float(src.bird(src.r_a, 0, config.breakdown.definition)) < threshold <= float(src.bird(inlet, 0, config.breakdown.definition)):
        raise ValueError('inlet must remain upstream of the available breakdown crossing')
    return RadialModel(config, ref, src, inlet, outer, _physical_hash(config,kn,inlet,outer))


def _physical_hash(config,kn,inlet,outer):
    physical = {**config.physical_inputs, 'exit': {**config.physical_inputs['exit'], 'knudsen': float(kn)}}
    identity = {'reference': physical, 'r_in_m': inlet, 'r_out_m': outer,
                'mapping': 'cai2012ext_markelov_streamtube', 'axis': [1, 0, 0]}
    return digest(identity)


def _restore_source(template,record):
    """Restore recorded coefficients; reuse the existing source's equations."""
    e,t,m,c = (record[k] for k in ('exit','targets','markelov','closure'))
    targets = replace(template.targets,**{key:t[name] for key,name in (
        ('area_m2','area_m2'),('number_flux_per_m2_s','gross_number_flux_per_m2_s'),
        ('energy_flux_W_per_m2','gross_energy_flux_W_per_m2'),('mdot_kg_per_s','mdot_target_kg_per_s'),
        ('edot_W','edot_target_W'),('h_J_per_kg','h_target_J_per_kg'),
        ('mdot_net_kg_per_s','mdot_continuum_net_kg_per_s'),('h_net_J_per_kg','h_continuum_net_J_per_kg'))})
    return replace(template,gamma=record['gamma'],cp=record['cp_J_per_kg_K'],
        gas_constant=record['gas_constant_J_per_kg_K'],convention=record['mean_free_path_convention'],
        diameter_m=e['D_m'],knudsen=e['Kn'],n_e=e['n_e_per_m3'],T_e=e['T_e_K'],
        U_e=e['U_e_m_per_s'],rho_e=e['rho_e_kg_per_m3'],targets=targets,
        theta_L=m['theta_L_rad'],k_theta=m['k_theta'],I=m['I'],A_M=m['A_M'],
        exponent_offset=m['exponent_offset'],p_t_r_o2=m['p_t_r_o2_Pa_m2'],
        T_t=c['T_t_K'],T_a=c['T_a_K'],U_a=c['U_a_m_per_s'],M_a=c['M_a'],
        r_a=c['r_a_m'],r_star=c['r_star_m'],V_inf=c['V_inf_m_per_s'],
        source_point=tuple(record['source_point_m']),exit_origin=tuple(record['exit_origin_m']),axis=tuple(record['axis']))


def require_same_physics(parent, child):
    if parent != child:
        raise ValueError('numerical revision changes physical source/domain/inlet identity')


def predict_breakdown(model, *, definition='general', threshold=.05) -> dict:
    result = {'definition': definition, 'threshold': threshold, 'basis': 'predicted_axis_surrogate'}
    try:
        if definition not in ('general', 'radial_approx') or not np.isfinite(threshold) or threshold <= 0:
            raise ValueError('invalid Bird definition or threshold')
        if not np.isfinite(model.r_in + model.r_out) or not 0 < model.r_in < model.r_out:
            raise ValueError('invalid domain')
        r = np.geomspace(model.r_in, model.r_out, 257)
        st = model.state(r)
        if any(np.any(~np.isfinite(st[k])) | np.any(np.asarray(st[k]) <= 0)
               for k in ('n', 'T', 'U', 'lambda')):
            raise ValueError('nonpositive or nonfinite source state')
        b = np.asarray(model.source.bird(r, 0., definition))
        if not np.all(np.isfinite(b)):
            raise ValueError('nonfinite Bird parameter')
        result.update(B_in=float(b[0]), B_out=float(b[-1]))
        if b[0] >= threshold:
            return {**result, 'status': 'ALREADY_NONCONTINUUM_AT_INLET'}
        crossing = np.flatnonzero((b[:-1] < threshold) & (b[1:] >= threshold))
        if not len(crossing):
            return {**result, 'status': 'NOT_REACHED_IN_DOMAIN'}
        i = int(crossing[0])
        root = brentq(lambda x: float(model.source.bird(x, 0., definition)) - threshold, r[i], r[i+1], xtol=1e-13)
        return {**result, 'status': 'CROSSING_IN_DOMAIN', 'r_m': root,
                'r_over_D': root / model.source.diameter_m,
                'state': {k: float(v) for k, v in model.state(root).items()}}
    except (ValueError, TypeError, FloatingPointError) as exc:
        return {**result, 'status': 'INVALID_STATE', 'reason': str(exc)}
