"""Summarize paired, locally oriented increments of conservative LF tendencies.

Caller supplies independently chosen anchors and their physical paired endpoints.
Orientations describe the displacement at the anchor; curved-field endpoints must
be supplied by the tracer. This module never calls Cartesian shells parallel scales.
"""
import numpy as np


def build_pairs(fields, spacing, scales, n, seed=1729):
    """Return xyz integer pairs[S,D,N,2,3] and displacement metadata.

    Scales are a name-to-physical-length mapping. Each is sampled at .5,1,2
    times its supplied length. Endpoints use nearest cell centers. Chords are
    oriented with the field at the first endpoint and are not curved traces.
    """
    fields = np.asarray(fields)
    spacing = np.asarray(spacing, dtype=float)
    if fields.ndim != 4 or fields.shape[-1] != 6 or spacing.shape != (3,):
        raise ValueError('Expected fields[z,y,x,6] and spacing[x,y,z]')
    if not np.isfinite(spacing).all() or np.any(spacing <= 0) or not scales or any(not np.isfinite(v) or v <= 0 for v in scales.values()):
        raise ValueError('Positive grid spacing and physical scales required')
    shape = np.array(fields.shape[:3][::-1])
    rng = np.random.default_rng(seed)
    count = min(int(n), int(np.prod(shape)))
    if count < 2:
        raise ValueError('At least two anchors required')
    ids = rng.choice(int(np.prod(shape)), count, replace=False)
    anchors = np.column_stack(np.unravel_index(ids, fields.shape[:3]))[:, ::-1]
    B = fields[anchors[:,2], anchors[:,1], anchors[:,0], :3]
    norm = np.linalg.norm(B,axis=1)
    if np.any(norm <= 0) or not np.isfinite(B).all():
        raise ValueError('Anchor has invalid magnetic direction')
    parallel = B/norm[:,None]
    basis = np.eye(3)[np.argmin(abs(parallel),axis=1)]
    perpendicular = np.cross(parallel,basis)
    perpendicular /= np.linalg.norm(perpendicular,axis=1)[:,None]
    second = np.cross(parallel, perpendicular)
    azimuth = rng.uniform(0, 2*np.pi, count)
    perpendicular = np.cos(azimuth)[:,None]*perpendicular + np.sin(azimuth)[:,None]*second
    pairs, metadata, lengths = [], [], []
    for name, length in scales.items():
        for factor in (.5,1.,2.):
            target = factor*float(length)
            oriented = []
            record = {'scale':str(name),'factor':factor,'target_separation':target,'directions':[]}
            for direction, vector in zip(('parallel','perpendicular'),(parallel,perpendicular)):
                displacement = np.rint(target*vector/spacing).astype(int)
                actual = displacement*spacing
                distance = np.linalg.norm(actual,axis=1)
                cosine = np.divide(np.sum(actual*vector,axis=1),distance,out=np.zeros(count),where=distance>0)
                endpoints = (anchors+displacement)%shape
                oriented.append(np.stack((anchors,endpoints),axis=1))
                record['directions'].append(dict(direction=direction,
                    actual_length_min=float(distance.min()),actual_length_median=float(np.median(distance)),
                    actual_length_max=float(distance.max()),
                    actual_to_target_quantiles=np.quantile(distance/target,[0,.1,.5,.9,1]).tolist(),
                    zero_displacement_fraction=float(np.mean(distance==0)),
                    median_angle_error_degrees=float(np.median(np.degrees(np.arccos(np.clip(cosine,-1,1))))),
                    below_four_grid_spacings=bool(target<4*max(spacing)),
                    longer_than_half_shortest_box=bool(target>.5*min(shape*spacing))))
            pairs.append(oriented); metadata.append(record); lengths.append(target)
    return np.asarray(pairs,dtype=np.int64), dict(scales=metadata,lengths=lengths,
        directions=['parallel','perpendicular'],anchors=count,seed=int(seed),
        coordinate_order='xyz',array_order='zyx',
        geometry='nearest-cell straight chord oriented to local field at anchor; perpendicular azimuth uniform and held fixed across scales',
        uncertainty='anchors can be correlated; batch standard errors are descriptive')


def summarize(delta, tendency, methods, lengths, directions, reference=-1, batches=8):
    """delta[S,D,N,2], tendency[M,S,D,N,2], lengths[S] in physical box units."""
    delta = np.asarray(delta, dtype=float)
    tendency = np.asarray(tendency, dtype=float)
    if delta.ndim != 4 or delta.shape[-1] != 2 or tendency.ndim != 5 or tendency.shape[1:] != delta.shape:
        raise ValueError('Expected delta[S,D,N,2] and tendency[M,S,D,N,2]')
    if len(methods) != len(tendency) or len(lengths) != delta.shape[0] or len(directions) != delta.shape[1]:
        raise ValueError('Labels do not match array dimensions')
    if not np.isfinite(delta).all() or not np.isfinite(tendency).all():
        raise ValueError('Nonfinite paired values')
    if not len(methods) or not np.isfinite(lengths).all() or min(lengths) <= 0 or delta.shape[2] < 2 or batches < 1:
        raise ValueError('Positive physical separations and at least two pairs required')
    if not -len(methods) <= reference < len(methods):
        raise ValueError('Reference index is outside method list')
    reference = reference % len(methods)
    inc = delta[..., 1] - delta[..., 0]
    rate = tendency[..., 1] - tendency[..., 0]
    output = {'reference': str(methods[reference]), 'rows': [],
              'definition': 'mean(delta_Delta_p * delta_dDelta_p_dt); half the LF contribution to the increment structure-function derivative',
              'uncertainty': 'standard error across consecutive anchor batches; descriptive only because turbulent anchors can be spatially correlated'}
    for s, length in enumerate(lengths):
        for d, direction in enumerate(directions):
            ref = rate[reference, s, d]
            scale = np.sqrt(np.mean(ref**2))
            significant = abs(ref) > max(1e-30, 1e-4*scale)
            for m, method in enumerate(methods):
                q = rate[m, s, d]
                transfer = inc[s, d] * q
                groups = np.array_split(transfer, min(batches, len(q)//2))
                means = np.array([x.mean() for x in groups])
                paired = inc[s,d]*(q-ref)
                paired_means = np.array([x.mean() for x in np.array_split(paired,len(groups))])
                paired_error = float(paired_means.std(ddof=1)/np.sqrt(len(groups))) if len(groups)>1 else None
                paired_mean = float(paired.mean())
                error = np.sqrt(np.mean((q-ref)**2))
                row = dict(separation=float(length), direction=str(direction), method=str(method), pairs=len(q),
                           increment_variance=float(np.mean(inc[s,d]**2)),
                           half_structure_function_rate=float(transfer.mean()),
                           transfer_batch_standard_error=float(means.std(ddof=1)/np.sqrt(len(means))) if len(means)>1 else None,
                           transfer_difference_from_reference=paired_mean,
                           paired_transfer_difference_standard_error=paired_error,
                           paired_difference_compatible_with_zero=bool(abs(paired_mean)<=2*paired_error) if paired_error is not None else None,
                           sampling_warning=('Too few batches to estimate sampling scatter' if len(groups)<4 else
                               'Paired difference is within two descriptive batch standard errors' if paired_error and abs(paired_mean)<=2*paired_error else None),
                           tendency_increment_rms=float(np.sqrt(np.mean(q*q))),
                           difference_rms=float(error), relative_rms_difference=float(error/scale) if scale else None,
                           opposite_sign_fraction=float(np.mean(q[significant]*ref[significant]<0)) if significant.any() else None)
                output['rows'].append(row)
    return output


def selfcheck():
    phase = np.linspace(0, 2*np.pi, 128, endpoint=False)
    delta = np.stack((np.sin(phase), np.sin(phase+.3)), axis=-1)[None,None]
    rate = np.stack((-delta, -2*delta))
    result = summarize(delta, rate, ['one','two'], [.3], ['parallel'])
    a,b = result['rows']
    assert np.isclose(b['half_structure_function_rate'], 2*a['half_structure_function_rate'])
    assert a['half_structure_function_rate'] < 0
    assert np.isclose(a['relative_rms_difference'], .5)
    assert b['relative_rms_difference'] == 0 and b['opposite_sign_fraction'] == 0
    same = summarize(delta,rate[:1],['one'],[.3],['parallel'])['rows'][0]
    assert same['difference_rms']==0 and same['transfer_difference_from_reference']==0
    fields = np.ones((32,32,32,6)); fields[...,:3] = [0,0,1]
    pairs, meta = build_pairs(fields,[1/32]*3,{'test':.25},1024,38)
    pairs2,_ = build_pairs(fields,[1/32]*3,{'test':.25},1024,38)
    assert np.array_equal(pairs,pairs2)
    disp=(pairs[1,1,:,1]-pairs[1,1,:,0]+16)%32-16
    assert np.all(disp[:,2]==0) and abs(np.mean(disp[:,0]**2)-np.mean(disp[:,1]**2))<3
    return True


if __name__ == '__main__':
    selfcheck()
    print('Scale-analysis signed sinusoid check passed.')
