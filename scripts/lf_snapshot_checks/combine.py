"""Combine offline scale comparisons and plot parallel increments at the supplied scales."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def combine(paths, destination):
    tables = []
    snapshots = set()
    columns = []
    keys = set()
    required = {'snapshot','time','cycle','nu_coll','scale','factor','direction',
                'method','family','reference','half_structure_function_rate',
                'relative_rms_difference','resolution_flag'}
    for path in paths:
        source = Path(path)/'scales.csv'
        with source.open(newline='') as stream:
            reader = csv.DictReader(stream)
            if not required.issubset(reader.fieldnames or []):
                raise ValueError(f'Missing scale columns in {source}')
            for name in reader.fieldnames:
                if name not in columns:
                    columns.append(name)
            rows = list(reader)
        if not rows:
            raise ValueError(f'Empty scale table {source}')
        identity = {(r['snapshot'],r['time'],r['cycle']) for r in rows}
        if len(identity) != 1 or identity & snapshots:
            raise ValueError(f'Mixed or duplicate snapshot in {source}')
        snapshots |= identity
        if len({r['nu_coll'] for r in rows}) != 1:
            raise ValueError(f'Inconsistent collision setting in {source}')
        for row in rows:
            key = tuple(row[k] for k in ('snapshot','time','cycle','scale','factor','direction','method','family','reference'))
            if key in keys:
                raise ValueError(f'Duplicate comparison in {source}')
            keys.add(key)
        tables.append(rows)
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=False)
    with (destination/'scales.csv').open('w',newline='') as stream:
        writer = csv.DictWriter(stream,fieldnames=columns)
        writer.writeheader()
        for rows in tables:
            writer.writerows(rows)
    center = [[r for r in rows if r['direction']=='parallel' and np.isclose(float(r['factor']),1)] for rows in tables]
    names = sorted({r['scale'] for rows in center for r in rows})
    if not names:
        raise ValueError('No parallel, unit-factor comparisons')
    def shortlist(rows):
        rows = [r for r in rows if r['method'] != r['reference']]
        receiver = [r for r in rows if r['family']=='receiver']
        fixed = [r for r in receiver if r['method']=='fixed_k']
        traces = [r for r in receiver if r['method']!='fixed_k']
        receiver = fixed + traces[:1] + (traces[-1:] if len(traces)>1 else [])
        uniform = [r for r in rows if r['family']=='uniform kernel']
        if uniform:
            uniform = [max(uniform, key=lambda r: float(r['method'].split('R=')[-1]))]
        optical = [r for r in rows if r['family']=='path scattering']
        return receiver + uniform + optical
    size = (max(7,6*len(names)),max(3.5,3.7*len(tables)))
    fig, axes = plt.subplots(len(tables),len(names),figsize=size,squeeze=False,layout='constrained')
    transfer, taxes = plt.subplots(len(tables),len(names),figsize=size,squeeze=False,layout='constrained')
    for i,rows in enumerate(center):
        for j,name in enumerate(names):
            ax,tax = axes[i,j],taxes[i,j]
            chosen = shortlist([r for r in rows if r['scale']==name])
            if not chosen:
                for axis in (ax,tax):
                    axis.text(.5,.5,'No non-reference comparisons',ha='center',transform=axis.transAxes)
                continue
            labels = [r['method']+' → '+r['reference'] for r in chosen]
            positions = np.arange(len(chosen))
            values = [100*float(r['relative_rms_difference']) if r['relative_rms_difference'] else np.nan for r in chosen]
            colors = ['#b45526' if r['resolution_flag'].lower() not in ('','false','ok','resolved','0') else '#287d88' for r in chosen]
            ax.barh(positions,values,color=colors)
            ax.set_yticks(positions,labels,fontsize=7)
            ax.invert_yaxis()
            ax.set_title(f"{name}, ν={chosen[0]['nu_coll']}, t={chosen[0]['time']}")
            ax.grid(axis='x',alpha=.2)
            # Select the same methods as the central panel, then retain every factor.
            allrows = [r for r in tables[i] if r['scale']==name and r['direction']=='parallel']
            series = [(r['method'],r['family'],r['reference']) for r in chosen]
            refkeys = {(r['reference'],r['family'],r['reference']) for r in chosen}
            for method,family,ref in series + sorted(refkeys):
                group = sorted([r for r in allrows if (r['method'],r['family'],r['reference'])==(method,family,ref)],
                               key=lambda r:float(r['factor']))
                if not group:
                    continue
                factors = [float(r['factor']) for r in group]
                values = [float(r['half_structure_function_rate']) for r in group]
                errors = [float(r['transfer_batch_standard_error']) if r.get('transfer_batch_standard_error') else 0 for r in group]
                label = method.replace('fixed_k','Fixed k').replace('SHD_receiver_','Field-line ').replace('rational_uniform_','Rational ').replace('SHD_uniform_','Uniform reference ')
                line=tax.errorbar(factors,values,yerr=errors,fmt='o-',ms=4,lw=1,capsize=2,label=label)
                flagged=[r['resolution_flag'].lower() not in ('','false','ok','resolved','0') for r in group]
                for factor,value,flag in zip(factors,values,flagged):
                    if flag:
                        tax.plot(factor,value,'x',color='#b45526',ms=9,mew=1.8)
            tax.set_title(f"{name}, ν={chosen[0]['nu_coll']}, t={chosen[0]['time']}")
            tax.set_xticks([.5,1,2])
            tax.set_xlabel('Separation / supplied scale')
            tax.set_ylabel('Signed half-structure-function rate')
            tax.grid(alpha=.2)
            tax.legend(fontsize=6,frameon=False,loc='best')
            tax.axhline(0,color='0.5',lw=.7)
            ax.set_xlabel('Tendency-increment RMS difference (%)')
    caption='Parallel chords at supplied scales; orange carries a resolution flag'
    fig.suptitle(caption,fontsize=10)
    transfer.suptitle(caption+'\nPoints and descriptive batch errors; orange crosses flag geometry or resolution limits',fontsize=10)
    fig.savefig(destination/'scale_comparison.png',dpi=160)
    transfer.savefig(destination/'scale_transfer.png',dpi=160)
    plt.close(fig)
    plt.close(transfer)
    return destination


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results',type=Path,nargs='+')
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    print(combine(args.results,args.out))
