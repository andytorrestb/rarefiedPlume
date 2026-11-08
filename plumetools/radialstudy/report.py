"""CSV/JSON/Markdown always; optional plots require matplotlib only on demand."""
import csv
import json
from pathlib import Path

import numpy as np


def _json(path,value):
    path.write_text(json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n')


def write_case_report(root,result,*,plots=True):
    root = Path(root)
    root.mkdir(parents=True,exist_ok=True)
    _json(root/'summary.json',result)
    profile = result.get('profile')
    if profile:
        names = ['r_m','n','Ux','T','lambda_predicted_m','lambda_measured_m','max_edge_m',
                 'mean_parcels','max_collision_ratio','max_transport_ratio']
        with (root/'radial.csv').open('w',newline='') as stream:
            writer = csv.writer(stream); writer.writerow(names)
            writer.writerows(zip(*(profile[k] for k in names)))
    if plots and profile:
        try:
            import matplotlib
            matplotlib.use('Agg')
            from matplotlib import pyplot as plt
        except ImportError:
            return
        d = result['diameter_m']
        r = np.asarray(profile['r_m'])/d
        fig,axes = plt.subplots(3,2,figsize=(11,11),layout='constrained')
        for key,ax,pred,norm in [('n',axes[0,0],'n_predicted','n_over_ne'),('Ux',axes[0,1],'U_predicted','Ux_over_Ue'),('T',axes[1,0],'T_predicted','T_over_Te')]:
            values = np.asarray(profile[key],dtype=float)
            error = np.asarray(profile.get('se',{}).get(key,np.full(len(r),np.nan)),dtype=float)
            ax.errorbar(r,values,yerr=error,label='measured ± block SE',linewidth=1)
            ax.plot(r,profile[pred],label='isentropic reference')
            ax.set_ylabel(key+' (SI)'); ax.legend(fontsize=8)
        ax = axes[1,1]
        ax.semilogy(r,profile['lambda_predicted_m'],label='predicted lambda')
        ax.semilogy(r,profile['lambda_measured_m'],label='measured lambda')
        ax.set_ylabel('mean free path (m)'); ax.legend(fontsize=8)
        ax = axes[2,0]
        for key in ('max_edge_over_measured_mfp','max_edge_over_predicted_mfp'):
            ax.plot(r,profile[key],label=key)
        for y in (.8,1,1/3):
            ax.axhline(y,color='grey',linestyle='--',linewidth=.8)
        ax.set_ylabel('largest edge / lambda'); ax.legend(fontsize=8)
        ax = axes[2,1]
        ax.plot(r,profile['mean_parcels'],label='mean parcels/cell')
        ax.axhline(result['assessment']['mean_parcels_target'],color='grey',linestyle='--')
        ax.set_ylabel('mean parcels/cell'); ax.legend(fontsize=8)
        crossing = result['predicted_breakdown'].get('r_m')
        for ax in axes.flat:
            ax.set_xlabel('nominal r / D_e')
            if crossing:
                ax.axvline(crossing/d,color='orange',linestyle=':')
        fig.savefig(root/'profiles.png',dpi=140); plt.close(fig)
        fig,ax = plt.subplots(layout='constrained')
        for key,label in [('max_collision_ratio','deltaT / tau'),('max_transport_ratio','transport ratio')]:
            ax.plot(r,profile[key],label=label)
        for key in ('collision_limit','transport_limit'):
            ax.axhline(result['timestep'][key],color='grey',linestyle='--')
        ax.set_xlabel('nominal r / D_e'); ax.set_ylabel('timestep ratio'); ax.legend()
        fig.savefig(root/'timestep.png',dpi=140); plt.close(fig)


def write_study_report(root,analysis,*,plots=True):
    root = Path(root)
    root.mkdir(parents=True,exist_ok=True)
    _json(root/'summary.json',analysis)
    fields = ['path','Kn','status','n_cells','particle_weight','planned_steps','max_edge_over_mfp','min_mean_parcels','breakdown_status']
    with (root/'summary.csv').open('w',newline='') as stream:
        writer = csv.DictWriter(stream,fieldnames=fields); writer.writeheader()
        for c in analysis['cases']:
            writer.writerow({**{k:c.get(k) for k in ('path','Kn','status')},
                             **{k:c['numerics'].get(k) for k in ('n_cells','particle_weight','planned_steps')},
                             **{k:c.get('assessment',{}).get(k) for k in ('max_edge_over_mfp','min_mean_parcels')},
                             'breakdown_status':c.get('predicted_breakdown',{}).get('status')})
    lines = ['# Radial DSMC sizing study','',analysis['scope'],'',
             '| Kn | Cells | Status | Remaining evidence or correction |','|---:|---:|---|---|']
    for i,c in enumerate(analysis['cases']):
        lines.append(f"| {c['Kn']:g} | {c['numerics']['n_cells']} | {c['status']} | {c.get('reason',c['status'])} |")
        write_case_report(root/f'case-{i:03d}',c,plots=plots)
    lines += ['','## Available comparisons','',json.dumps(analysis['comparisons'],indent=2),
              '', 'A satisfied case requires no further mesh sweep. Missing comparisons do not fail the measured mesh criterion.','']
    (root/'report.md').write_text('\n'.join(lines))
    return sorted(p for p in root.rglob('*') if p.is_file())
