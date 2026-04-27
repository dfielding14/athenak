# Regression test for face-centered div(B) preservation under moving AMR.
#
# The divb_amr pgen initializes B from a discrete vector potential, then forces a moving
# AMR refinement pattern.  This script checks the direct face-centered divergence history
# in 1D, 2D, 3D, and a turbulence-driven 3D stress case.

import logging
import glob
import math
import os
import sys

import scripts.utils.athena as athena

sys.path.insert(0, '../vis/python')
import athena_read  # noqa

athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])

_CASES = [
    ('1D', 'tests/divb_amr_1d.athinput', 'DivBAMR1D', 64, 24, []),
    ('2D', 'tests/divb_amr_2d.athinput', 'DivBAMR2D', 48*48, 20, []),
]
_CASES += [
    ('3D-L{0}'.format(level), 'tests/divb_amr_3d.athinput',
     'DivBAMR3DL{0}'.format(level), 32*32*32, 8,
     ['mesh_refinement/num_levels={0}'.format(level + 1),
      'problem/refine_levels={0}'.format(level),
      'mesh_refinement/max_nmb_per_rank=16384',
      'time/nlim=8'])
    for level in range(1, 6)
]
_CASES += [
    ('3D+turb', 'tests/divb_amr_turb_3d.athinput', 'DivBAMR3DTurb', 32*32*32, 12, []),
]

_MAX_NDIV_TOL = 2.0e-11
_L1_NDIV_TOL = 2.0e-12
_L2_NDIV_TOL = 5.0e-12

_FIELD_DIAGNOSTICS = [
    ('2D', 'DivBAMR2D', [('x3', 'bfield_slice_x3', 'divb_slice_x3', 'j2_slice_x3')]),
    ('3D-L5', 'DivBAMR3DL5', [('x3', 'bfield_slice_x3', 'divb_slice_x3', 'j2_slice_x3'),
                              ('x2', 'bfield_slice_x2', 'divb_slice_x2', 'j2_slice_x2')]),
    ('3D+turb', 'DivBAMR3DTurb', [('x3', 'bfield_slice_x3', 'divb_slice_x3', 'j2_slice_x3'),
                                  ('x2', 'bfield_slice_x2', 'divb_slice_x2',
                                   'j2_slice_x2')]),
]

_GAMMA = 5.0/3.0
_DIVB_BNORM = max(1.0, abs(0.7) + abs(0.2) + abs(-0.15) + 6.0*abs(0.25))


def run(**kwargs):
    logger.debug('Running test ' + __name__)
    for _, input_file, basename, _, _, arguments in _CASES:
        hst_file = os.path.join('build', 'src', basename + '.user.hst')
        if os.path.exists(hst_file):
            os.remove(hst_file)
        for fname in glob.glob(os.path.join('build', 'src', 'bin',
                                            basename + '.*slice*.bin')):
            os.remove(fname)
        athena.run(input_file, [
            'job/basename=' + basename,
            'output1/dcycle=1',
        ] + arguments)


def _latest_binary(basename, file_id):
    pattern = os.path.join('build', 'src', 'bin',
                           '{0}.{1}.*.bin'.format(basename, file_id))
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError('missing diagnostic binary output: ' + pattern)
    return files[-1]


def _slice_array(data, name):
    arr = data[name]
    if arr.shape[0] == 1:
        return (arr[0, :, :], data['x1v'], data['x2v'],
                [data['x1f'][0], data['x1f'][-1], data['x2f'][0], data['x2f'][-1]],
                r'$x_1$', r'$x_2$')
    if arr.shape[1] == 1:
        return (arr[:, 0, :], data['x1v'], data['x3v'],
                [data['x1f'][0], data['x1f'][-1], data['x3f'][0], data['x3f'][-1]],
                r'$x_1$', r'$x_3$')
    if arr.shape[2] == 1:
        return (arr[:, :, 0], data['x2v'], data['x3v'],
                [data['x2f'][0], data['x2f'][-1], data['x3f'][0], data['x3f'][-1]],
                r'$x_2$', r'$x_3$')
    raise ValueError('diagnostic output is not a slice: shape={0}'.format(arr.shape))


def _plot_magnetic_slice(label, basename, plane, field_id, divb_id, j2_id):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np

    sys.path.insert(0, '../vis/python')
    import bin_convert_new  # noqa

    field = bin_convert_new.read_binary_as_athdf(
        _latest_binary(basename, field_id), quantities=['bcc1', 'bcc2', 'bcc3', 'eint'],
        return_levels=True)
    divb = bin_convert_new.read_binary_as_athdf(
        _latest_binary(basename, divb_id), quantities=['divb'])
    j2 = bin_convert_new.read_binary_as_athdf(
        _latest_binary(basename, j2_id), quantities=['j2'])

    bx, x, y, extent, xlabel, ylabel = _slice_array(field, 'bcc1')
    by = _slice_array(field, 'bcc2')[0]
    bz = _slice_array(field, 'bcc3')[0]
    eint = _slice_array(field, 'eint')[0]
    levels = _slice_array(field, 'Levels')[0]
    divb_slice = _slice_array(divb, 'divb')[0]
    jmag = np.sqrt(np.maximum(_slice_array(j2, 'j2')[0], 0.0))

    b2 = bx*bx + by*by + bz*bz
    bmag = np.sqrt(np.maximum(b2, 0.0))
    pmag = 0.5*b2
    beta = np.where(b2 > 0.0, 2.0*(_GAMMA - 1.0)*eint/b2, np.nan)
    dx = min([np.min(np.diff(coord)) for coord in (field['x1f'], field['x2f'], field['x3f'])
              if len(coord) > 2])
    ndivb = np.abs(divb_slice)*dx/_DIVB_BNORM

    panels = [
        (bx, r'$B_x$', 'RdBu_r', True),
        (by, r'$B_y$', 'RdBu_r', True),
        (bz, r'$B_z$', 'RdBu_r', True),
        (bmag, r'$|B|$', 'magma', False),
        (pmag, r'$P_B=B^2/2$', 'viridis', False),
        (np.log10(np.maximum(beta, 1.0e-30)), r'$\log_{10}\beta$', 'plasma', False),
        (np.log10(np.maximum(jmag, 1.0e-30)), r'$\log_{10}|J|$', 'inferno', False),
        (np.log10(np.maximum(ndivb, 1.0e-30)), r'$\log_{10}(|\nabla\cdot B|\Delta x/B)$',
         'cividis', False),
        (levels, 'AMR level', 'tab20', False),
    ]

    fig, axes = plt.subplots(3, 3, figsize=(13.5, 11.0), constrained_layout=True)
    for ax, (arr, title, cmap, symmetric) in zip(axes.flat, panels):
        finite = arr[np.isfinite(arr)]
        if finite.size == 0:
            vmin = vmax = None
        elif symmetric:
            vmax = np.percentile(np.abs(finite), 99.5)
            vmax = vmax if vmax > 0.0 else 1.0
            vmin = -vmax
        else:
            vmin = np.percentile(finite, 0.5)
            vmax = np.percentile(finite, 99.5)
            if vmin == vmax:
                vmin = vmax = None
        image = ax.imshow(arr, origin='lower', extent=extent, cmap=cmap,
                          vmin=vmin, vmax=vmax, interpolation='nearest', aspect='equal')
        if np.nanmax(levels) > np.nanmin(levels):
            contour_levels = np.arange(np.nanmin(levels) + 0.5, np.nanmax(levels) + 0.5)
            ax.contour(x, y, levels, levels=contour_levels, colors='k',
                       linewidths=0.25, alpha=0.35)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        fig.colorbar(image, ax=ax, shrink=0.82)

    fig.suptitle('{0} moving-AMR magnetic diagnostics, {1} slice, t={2:g}, cycle={3}'.
                 format(label, plane, field['Time'], field['NumCycles']))
    outname = os.path.join('build', 'src',
                           '{0}_magnetic_{1}_diagnostics.png'.format(basename, plane))
    fig.savefig(outname, dpi=180)
    plt.close(fig)
    logger.info('wrote magnetic AMR diagnostic slice: ' + outname)


def _make_magnetic_slice_plots():
    try:
        for label, basename, planes in _FIELD_DIAGNOSTICS:
            for plane in planes:
                _plot_magnetic_slice(label, basename, *plane)
    except ImportError as err:
        logger.warning('Skipping magnetic AMR slice plots; plotting dependency missing: {0}'.
                       format(err))
    except FileNotFoundError as err:
        logger.warning(str(err))
        raise


def analyze():
    logger.debug('Analyzing test ' + __name__)
    analyze_status = True

    for label, _, basename, root_ncell, min_rows, _ in _CASES:
        hst_file = os.path.join('build', 'src', basename + '.user.hst')
        data = athena_read.hst(hst_file)

        max_ndiv = max(data['max_ndiv'])
        l1_ndiv = max(s/v for s, v in zip(data['sum_ndiv'], data['vol']))
        l2_ndiv = max(math.sqrt(s/v) for s, v in zip(data['sum_n2'], data['vol']))
        max_ncell = max(data['ncell'])

        if len(data['time']) < min_rows:
            logger.warning('{0} AMR div(B) history is too short: {1} rows, expected {2}'.
                           format(label, len(data['time']), min_rows))
            analyze_status = False
        if max_ncell <= root_ncell:
            logger.warning('{0} AMR div(B) test did not refine: max_ncell={1:g}, '
                           'root_ncell={2:g}'.
                           format(label, max_ncell, root_ncell))
            analyze_status = False
        if max_ndiv > _MAX_NDIV_TOL:
            logger.warning('{0} max normalized div(B) too large: {1:g} threshold {2:g}'.
                           format(label, max_ndiv, _MAX_NDIV_TOL))
            analyze_status = False
        if l1_ndiv > _L1_NDIV_TOL:
            logger.warning('{0} L1 normalized div(B) too large: {1:g} threshold {2:g}'.
                           format(label, l1_ndiv, _L1_NDIV_TOL))
            analyze_status = False
        if l2_ndiv > _L2_NDIV_TOL:
            logger.warning('{0} L2 normalized div(B) too large: {1:g} threshold {2:g}'.
                           format(label, l2_ndiv, _L2_NDIV_TOL))
            analyze_status = False

    try:
        _make_magnetic_slice_plots()
    except FileNotFoundError:
        analyze_status = False

    return analyze_status
