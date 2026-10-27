"""Raw moment algebra; no clipping or covariance repairs."""
from __future__ import annotations
import numpy as np
from plumetools.cai2012.constants import BOLTZMANN_J_PER_K as KB


def derive(M0, M1, M2, mass_kg, *, relative_roundoff=1e-10):
    M0, M1, M2 = (np.asarray(a, dtype=float) for a in (M0, M1, M2))
    with np.errstate(divide='ignore', invalid='ignore'):
        U = M1/M0[..., None]
        tensor = mass_kg/KB*(M2/M0[..., None, None]-U[..., :, None]*U[..., None, :])
    finite = np.isfinite(tensor).all(axis=(-2, -1)) & np.isfinite(M0)
    eig = np.full(M0.shape+(3,), np.nan)
    if np.any(finite):
        eig[finite] = np.linalg.eigvalsh(tensor[finite])
    scale = np.maximum(np.max(np.abs(tensor), axis=(-2, -1)), 1.)
    symmetric = np.max(np.abs(tensor-np.swapaxes(tensor, -1, -2)), axis=(-2, -1)) <= relative_roundoff*scale
    valid = finite & (M0 > 0) & symmetric & (np.min(eig, axis=-1) >= -relative_roundoff*scale)
    return dict(M0=M0, M1=M1, M2=M2, n=M0, U=U, T_tensor=tensor,
                T_trans=np.trace(tensor, axis1=-2, axis2=-1)/3,
                P_tensor=M0[..., None, None]*KB*tensor, valid=valid,
                tensor_eigenvalues_K=eig, tensor_symmetric=symmetric,
                relative_roundoff_tolerance=relative_roundoff)


def aggregate(M0, M1, M2, weights, mass_kg):
    """Volume-weighted raw moments, then central moments (includes dispersion)."""
    w = np.asarray(weights, dtype=float)
    if np.any(w < 0) or not np.isfinite(w).all() or w.sum() <= 0:
        raise ValueError('aggregation weights must be finite, nonnegative, nonempty')
    w = w/w.sum()
    return derive(np.einsum('i,i->', w, M0), np.einsum('i,ij->j', w, M1),
                  np.einsum('i,ijk->jk', w, M2), mass_kg)
