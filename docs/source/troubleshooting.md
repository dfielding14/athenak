# Troubleshooting Guide

## Build Issues

### Kokkos Not Found
```bash
git submodule update --init --recursive
```

### CUDA Not Detected
```bash
export CUDA_ROOT=/usr/local/cuda
export PATH=$CUDA_ROOT/bin:$PATH
```

### MPI Issues
```bash
# For GPU-aware MPI
export MPICH_GPU_SUPPORT_ENABLED=1
```

## Runtime Errors

### Segmentation Fault
- Check array bounds
- Verify ghost zone access
- Ensure proper memory allocation

### Wrong Results
- Verify input parameters
- Check boundary conditions
- Confirm reconstruction method

### Performance Issues
- Increase MeshBlock size for GPUs
- Check load balancing
- Profile with Kokkos tools

## Common Problems

### div(B) Growth in MHD
- Use constrained transport
- Enable FOFC if needed
- Check prolongation settings

### AMR Instabilities
- Reduce refinement frequency
- Adjust refinement criteria
- Enable first-order flux correction

### Particle Issues
- Check particle boundary conditions
- Verify pusher stability
- Monitor particle counts

## Getting Help

- Check existing documentation
- Review example problems
- Post issues on GitHub
- Contact developers

## See Also
- [Migration Guide](migration/from_athena_plus_plus.md)
- [Common Gotchas](migration/common_gotchas.md)
- [Configuration Guide](configuration.md)