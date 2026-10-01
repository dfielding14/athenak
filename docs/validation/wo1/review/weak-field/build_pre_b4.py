"""Link an isolated pre-B4 flux control using an existing Makefiles CPU build.

Usage: python build_pre_b4.py REPOSITORY OUTPUT_DIRECTORY
"""

from pathlib import Path
import shlex
import subprocess
import sys

root = Path(sys.argv[1]).resolve()
base = Path(sys.argv[2]).resolve() / 'pre-b4'
base.mkdir(parents=True, exist_ok=True)
headers = ['mhd/rsolvers/hlle_cgl.hpp', 'mhd/rsolvers/llf_mhd_singlestate.hpp']
for header in headers:
    path = base / 'include' / header
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(subprocess.check_output(
        ['git', 'show', '7d1f79557^:src/' + header], cwd=root,
    ))

build = root / 'tst/build/src'
flags = {
    line.split(' = ')[0]: shlex.split(line.split(' = ')[1])
    for line in (build / 'CMakeFiles/athena.dir/flags.make').read_text().splitlines()
    if ' = ' in line
}
objects = {}
for name in ('mhd/mhd_fluxes', 'mhd/mhd_fofc', 'pgen/tests/cgl_fofc'):
    output = base / (Path(name).name + '.o')
    command = [
        '/usr/bin/c++', '-I' + str(base / 'include'), *flags['CXX_DEFINES'],
        *flags['CXX_INCLUDES'], *flags['CXX_FLAGS'], '-c',
        str(root / 'src' / f'{name}.cpp'), '-o', str(output),
    ]
    print(shlex.join(command), flush=True)
    subprocess.run(command, check=True)
    objects['CMakeFiles/athena.dir/' + name + '.cpp.o'] = str(output)

command = shlex.split((build / 'CMakeFiles/athena.dir/link.txt').read_text())
command = [objects.get(part, part) for part in command]
command[command.index('-o') + 1] = str(base / 'athena')
print(shlex.join(command), flush=True)
subprocess.run(command, cwd=build, check=True)
