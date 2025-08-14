# Running Simulations

## Basic Execution

```bash
./athena -i input.athinput
```

## Parallel Execution

### MPI
```bash
mpirun -np 16 ./athena -i input.athinput
```

### OpenMP
```bash
export OMP_NUM_THREADS=8
./athena -i input.athinput
```

### Hybrid MPI+OpenMP
```bash
export OMP_NUM_THREADS=4
mpirun -np 4 ./athena -i input.athinput
```

## Command Line Options

| Option | Description |
|--------|-------------|
| `-i` | Input file (required) |
| `-c` | Show configuration |
| `-n` | Parse input without running |
| `-m` | Show mesh structure |
| `-r` | Restart from checkpoint |
| `-d` | Output directory |
| `-t` | Final time (override input) |

## Monitoring Progress

AthenaK outputs diagnostic information during runtime:

```
cycle=0 time=0.0000e+00 dt=1.0000e-03
cycle=100 time=1.0000e-01 dt=1.0234e-03
```

## Output Files

### VTK Format
- Readable by ParaView, VisIt
- Files: `basename.block*.NNNNN.vtk`

### Binary Format
- Native binary format for efficient I/O
- Files: `basename.block*.NNNNN.bin`
- Use `single_file_per_rank=true` for parallel I/O

### Restart Files
- Binary format for exact restart
- Files: `basename.NNNNN.rst`

## Performance Tips

### GPU Execution
- Use multiple MeshBlocks per GPU
- Typical: 8-32 blocks per device
- Set environment: `export CUDA_VISIBLE_DEVICES=0,1`

### CPU Optimization
- Enable vectorization: `-march=native`
- Use OpenMP for shared memory
- Pin MPI processes to cores

## Troubleshooting

### Timestep Too Small
- Reduce `cfl_number`
- Check for unresolved features
- Increase resolution

### Out of Memory
- Reduce MeshBlock size
- Use more MPI ranks
- Enable restart dumps

### NaN in Solution
- Check initial conditions
- Verify EOS parameters
- Try first-order reconstruction