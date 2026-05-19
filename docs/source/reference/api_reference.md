# API Reference

## Core Classes

### Mesh
```cpp
class Mesh {
public:
  // Construction
  explicit Mesh(ParameterInput *pin);
  
  // Member variables
  int nmb_total;        // Total number of MeshBlocks across all levels/ranks
  int nmb_thisrank;     // Number of MeshBlocks on this MPI rank
  int nmb_maxperrank;   // Max allowed MBs per device (memory limit for AMR)
  int nmb_rootx1, nmb_rootx2, nmb_rootx3; // # of MeshBlocks at root level
  
  // Mesh properties
  RegionSize mesh_size;    // Physical size of mesh
  RegionIndcs mesh_indcs;  // Cell indices in mesh
  RegionIndcs mb_indcs;    // Cell indices in MeshBlocks
  BoundaryFlag mesh_bcs[6]; // Boundary conditions at 6 faces
  
  // Dimensionality flags
  bool one_d, two_d, three_d;
  bool multilevel;      // true for SMR and AMR
  bool adaptive;        // true only for AMR
};
```

### MeshBlock
```cpp
class MeshBlock {
public:
  // Block properties
  int mb_gid;           // Global ID of this MeshBlock
  int mb_lid;           // Local ID on this rank
  LogicalLoc lloc;      // Logical location in tree
  
  // Mesh refinement
  MeshRefinement *pmr;  // Pointer to mesh refinement
  
  // Block boundaries
  Real mb_size[3];      // Physical size of MeshBlock
};
```

### Driver
```cpp
class Driver {
public:
  // Main evolution
  void Execute();
  
  // Task management
  void BuildTaskList();
  void ExecuteTaskList();
  
  // Timestep
  Real CalculateNewTimestep();
};
```

## Physics Modules

### Hydro
```cpp
class Hydro {
public:
  // Constructor
  Hydro(MeshBlockPack *ppack, ParameterInput *pin);
  
  // Task assembly
  void AssembleHydroTasks(std::map<std::string, std::shared_ptr<TaskList>> tl);
  
  // Core methods (as tasks)
  TaskStatus InitRecv(Driver *pdrive, int stage);
  TaskStatus CopyCons(Driver *pdrive, int stage);
  TaskStatus Fluxes(Driver *pdrive, int stage);
  TaskStatus SendFlux(Driver *pdrive, int stage);
  TaskStatus RecvFlux(Driver *pdrive, int stage);
  TaskStatus ExpRKUpdate(Driver *pdrive, int stage);
  
  // Variables
  DvceArray5D<Real> u0;  // Conserved variables
  DvceArray5D<Real> w0;  // Primitive variables
  DvceArray5D<Real> u1;  // Conserved vars at intermediate step
  DvceFaceFld5D<Real> uflx; // Fluxes on cell faces
};
```

### MHD
```cpp
class MHD {
public:
  // Constructor
  MHD(MeshBlockPack *ppack, ParameterInput *pin);
  
  // Task assembly
  void AssembleMHDTasks(std::map<std::string, std::shared_ptr<TaskList>> tl);
  
  // Core MHD tasks
  TaskStatus InitRecv(Driver *pdrive, int stage);
  TaskStatus Fluxes(Driver *pdrive, int stage);
  TaskStatus CornerE(Driver *pdrive, int stage);
  TaskStatus CT(Driver *pdrive, int stage);  // Constrained Transport
  
  // Variables
  DvceArray5D<Real> u0;      // Conserved variables including B
  DvceArray5D<Real> w0;      // Primitive variables
  DvceFaceFld5D<Real> b0;    // Face-centered B-field
  DvceEdgeFld5D<Real> efld;  // Edge-centered electric field
};
```

## Kokkos Views

### Device Arrays
```cpp
// 5D array for cell-centered data (m,n,k,j,i)
// m = MeshBlock index, n = variable index, k,j,i = spatial indices
using DvceArray5D = Kokkos::View<T*****, LayoutWrapper, DevMemSpace>;

// 5D array for face-centered fields
using DvceFaceFld5D = Kokkos::View<T*****, LayoutWrapper, DevMemSpace>;

// 5D array for edge-centered fields
using DvceEdgeFld5D = Kokkos::View<T*****, LayoutWrapper, DevMemSpace>;
```

### Parallel Patterns
```cpp
// Parallel for loop
par_for("kernel_name", DevExeSpace(),
        kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
KOKKOS_LAMBDA(int k, int j, int i) {
  // Kernel body
});

// Parallel reduction
par_reduce("reduce", DevExeSpace(),
           kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
KOKKOS_LAMBDA(int k, int j, int i, Real &sum) {
  sum += value(k,j,i);
}, sum_total);
```

## Task System

### Task Registration
```cpp
// Tasks are registered to named task lists using string keys
TaskID none;
id.irecv = tl["before_stagen"]->AddTask(&Hydro::InitRecv, this, none);
id.flux = tl["stagen"]->AddTask(&Hydro::Fluxes, this, id.copyu);
id.update = tl["after_stagen"]->AddTask(&Hydro::ExpRKUpdate, this, id.flux);
```

### Task List Names
- `"before_stagen"` - Before stage N of RK integrator
- `"stagen"` - Stage N of RK integrator (n = 1,2,3,4)
- `"after_stagen"` - After stage N of RK integrator

## Problem Generator API

### Required Function
```cpp
void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (restart) return;
  
  // Get MeshBlockPack
  auto &mbp = pmy_mesh_->pmb_pack;
  
  // Set initial conditions using Kokkos parallel_for
  par_for("pgen", DevExeSpace(), 0, mbp->nmb_thispack-1,
          mbp->is, mbp->ie, mbp->js, mbp->je, mbp->ks, mbp->ke,
  KOKKOS_LAMBDA(int m, int i, int j, int k) {
    // Initialize cell (m,i,j,k)
  });
}
```

### Optional Functions
```cpp
// User source terms (defined in problem generator files)
void UserSource(Mesh* pm, const Real bdt);

// User-defined boundary conditions
// Set via user_bcs function pointer in problem generators
```

## Input/Output

### ParameterInput
```cpp
// Read parameters
Real value = pin->GetReal("block", "parameter");
int ivalue = pin->GetInteger("block", "parameter");
std::string svalue = pin->GetString("block", "parameter");
bool bvalue = pin->GetBoolean("block", "parameter");

// With defaults
Real value = pin->GetOrAddReal("block", "parameter", default_value);
```

### Output Classes
```cpp
class BaseOutput {
public:
  virtual void WriteOutputData(Mesh *pm) = 0;
};
```

## Boundary Conditions

### Standard Types (string names in input file)
- `"periodic"` - Periodic boundary
- `"outflow"` - Outflow/open boundary
- `"reflect"` - Reflecting boundary (NOT "reflecting")
- `"inflow"` - Inflow boundary
- `"diode"` - Diode (outflow-only) boundary
- `"vacuum"` - Vacuum boundary
- `"user"` - User-defined boundary
- `"shear_periodic"` - Shearing box periodic

### Setting in Input File
```ini
<mesh>
ix1_bc = outflow
ox1_bc = outflow
ix2_bc = periodic
ox2_bc = periodic
```

## Reconstruction

### Interface Functions
```cpp
template<typename T>
KOKKOS_INLINE_FUNCTION
void PLM(const T &q, T &ql, T &qr,
         const int il, const int iu);
```

## Riemann Solvers

### Solver Interface
```cpp
template<typename T>
KOKKOS_INLINE_FUNCTION
void HLLCRiemannSolver(const T &prim_l,
                       const T &prim_r,
                       T &flux);
```

## Constants and Enums

### Variable Indices
```cpp
// From athena.hpp:51-52
enum VariableIndex {
  IDN = 0,       // Density (both conserved and primitive)
  IM1 = 1,       // Momentum x1 (conserved)
  IVX = 1,       // Velocity x1 (primitive)
  IM2 = 2,       // Momentum x2 (conserved)
  IVY = 2,       // Velocity x2 (primitive)
  IM3 = 3,       // Momentum x3 (conserved)
  IVZ = 3,       // Velocity x3 (primitive)
  IEN = 4,       // Total energy (conserved)
  ITM = 4,       // Temperature (primitive, for some EOS)
  IPR = 4,       // Pressure (primitive)
  IYF = 5        // Electron fraction (for tabulated EOS)
};

// Magnetic field components
enum BFieldIndex {IBX=0, IBY=1, IBZ=2, NMAG=3};
```

## Macros

### Loop Bounds
```cpp
KOKKOS_INLINE_FUNCTION
void kernel(...) {
  // With ghost zones
  for (int k=kb.s; k<=kb.e; ++k) {
    for (int j=jb.s; j<=jb.e; ++j) {
      for (int i=ib.s; i<=ib.e; ++i) {
        // ...
      }
    }
  }
}
```

## Error Handling

### Assertions
```cpp
#if DEBUG
  assert(condition && "Error message");
#endif
```

### Kokkos Error Check
```cpp
KOKKOS_ASSERT(condition);
```

## See Also
- [File Reference](file_reference.md) - Complete file listing
- [Input Parameters](input_parameters.md) - Configuration options
- [Module Documentation](../modules/index.md) - Detailed module guides