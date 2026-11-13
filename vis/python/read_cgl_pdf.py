"""Read a shared or node-sharded AthenaK CGL ASCII volume histogram."""
from pathlib import Path
import numpy as np


def read_cgl_pdf(path):
    """Return (time, weights), including both open tails along each histogram axis."""
    path = Path(path)
    kind, _, index = path.parent.name.partition('_')
    if kind in ('node', 'rank') and len(index) == 8 and index.isdigit():
        folders = sorted(path.parent.parent.glob(kind + '_' + '[0-9]'*8))
        paths = [folder/path.name for folder in folders]
        if not paths:
            raise FileNotFoundError(path)
    else:
        paths = [path]
    total = None
    for file in paths:
        with file.open() as handle:
            stamp = handle.readline().split('=')
            if len(stamp) != 2 or stamp[0].strip() != '# time':
                raise ValueError(f'Invalid CGL PDF header in {file}')
            time = float(stamp[1])
            values = np.loadtxt(handle)
        if not np.isfinite(time) or not np.isfinite(values).all() or (values < 0).any():
            raise ValueError(f'Invalid histogram in {file}')
        if total is None:
            saved_time, total = time, values
        elif time != saved_time or total.shape != values.shape:
            raise ValueError(f'Inconsistent PDF shard {file}')
        else:
            total += values
    return saved_time, total
