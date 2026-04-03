# PDF Output Enhancement Plan

This document describes the implementation plan for extending the PDF output functionality in AthenaK to support:
1. Binning by spatial coordinates (Cartesian, spherical, cylindrical)
2. Arbitrary dimension PDFs (1D, 2D, 3D, 4D)
3. New derived variables for mass and energy fluxes

---

## Implementation Philosophy

- **Efficiency and robustness**: Write the most efficient and robust code possible
- **Ask, don't guess**: When uncertain about any design decision, implementation detail, or physics, ask for clarification rather than making assumptions
- **Test incrementally**: Test each piece as it is added before moving to the next
- **Match codebase style**: Study existing code patterns in AthenaK and adopt the same style and approach; look through similar implementations and copy that approach

---

## Design Decisions

- **Coordinate system**: Simulations are Cartesian only; derived coordinates (spherical r, θ, φ; cylindrical R, φ, z) computed from Cartesian x, y, z
- **Coordinate positions**: Use global domain coordinates via `CellCenterX()` function
- **Maximum dimensions**: Capped at 4D
- **Output format**: Binary files with separate ASCII header
- **Backward compatibility**: `variable` treated as alias for `variable_1`
- **Stride ordering**: C-style row-major (last dimension varies fastest, stride[ndim-1] = 1)
- **Angle ranges**: θ ∈ [0, π], φ ∈ [0, 2π]
- **MPI**: Requires GPU-aware MPI on GPU systems (standard on Frontier)

---

## Part 1: Add Coordinate Variables to Derived Variables

### 1.1 Add to `var_choice[]` in `outputs.hpp`

Add new coordinate variable names to the `var_choice` array:

```cpp
// Coordinate derived variables (add after existing entries, update NOUTPUT_CHOICES)
"coord_x", "coord_y", "coord_z",           // Cartesian
"coord_r", "coord_theta", "coord_phi",     // Spherical (r, θ, φ)
"coord_cyl_R", "coord_cyl_phi", "coord_cyl_z"  // Cylindrical (R, φ, z)
```

Update `NOUTPUT_CHOICES` accordingly (currently 158, add 9 → 167).

### 1.2 Register Coordinate Variables in `basetype_output.cpp`

In `BaseTypeOutput::BaseTypeOutput()`, add handling for each coordinate variable. Pattern follows existing derived variables like `temperature`:

```cpp
// Cartesian x coordinate
if (variable.compare("coord_x") == 0) {
  out_params.contains_derived = true;
  out_params.n_derived += 1;
  outvars.emplace_back("coord_x", out_params.i_derived, &derived_var);
}

// Cartesian y coordinate
if (variable.compare("coord_y") == 0) {
  out_params.contains_derived = true;
  out_params.n_derived += 1;
  outvars.emplace_back("coord_y", out_params.i_derived, &derived_var);
}

// Cartesian z coordinate
if (variable.compare("coord_z") == 0) {
  out_params.contains_derived = true;
  out_params.n_derived += 1;
  outvars.emplace_back("coord_z", out_params.i_derived, &derived_var);
}

// Spherical radius: r = sqrt(x² + y² + z²)
if (variable.compare("coord_r") == 0) {
  out_params.contains_derived = true;
  out_params.n_derived += 1;
  outvars.emplace_back("coord_r", out_params.i_derived, &derived_var);
}

// Spherical polar angle: θ = acos(z/r), range [0, π]
if (variable.compare("coord_theta") == 0) {
  out_params.contains_derived = true;
  out_params.n_derived += 1;
  outvars.emplace_back("coord_theta", out_params.i_derived, &derived_var);
}

// Spherical azimuthal angle: φ = atan2(y, x) shifted to [0, 2π]
if (variable.compare("coord_phi") == 0) {
  out_params.contains_derived = true;
  out_params.n_derived += 1;
  outvars.emplace_back("coord_phi", out_params.i_derived, &derived_var);
}

// Cylindrical radius: R = sqrt(x² + y²)
if (variable.compare("coord_cyl_R") == 0) {
  out_params.contains_derived = true;
  out_params.n_derived += 1;
  outvars.emplace_back("coord_cyl_R", out_params.i_derived, &derived_var);
}

// Cylindrical azimuthal angle: φ = atan2(y, x) shifted to [0, 2π]
if (variable.compare("coord_cyl_phi") == 0) {
  out_params.contains_derived = true;
  out_params.n_derived += 1;
  outvars.emplace_back("coord_cyl_phi", out_params.i_derived, &derived_var);
}

// Cylindrical z (same as Cartesian z, kept for user convenience)
if (variable.compare("coord_cyl_z") == 0) {
  out_params.contains_derived = true;
  out_params.n_derived += 1;
  outvars.emplace_back("coord_cyl_z", out_params.i_derived, &derived_var);
}
```

### 1.3 Compute Coordinate Variables in `derived_variables.cpp`

Add computation kernels in `BaseTypeOutput::ComputeDerivedVariable()`. Key pattern: branch on variable name **outside** the lambda, then use `CellCenterX()` inside.

```cpp
// Include at top of file (already present)
#include "coordinates/cell_locations.hpp"

// In ComputeDerivedVariable():

// Cartesian x coordinate
if (name.compare("coord_x") == 0) {
  if (derived_var.extent(4) <= 1)
    Kokkos::realloc(derived_var, nmb, n_dv, n3, n2, n1);
  auto dv = derived_var;
  int nx1 = indcs.nx1;
  par_for("coord_x", DevExeSpace(), 0, (nmb-1), ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    dv(m, i_dv, k, j, i) = CellCenterX(i-is, nx1, x1min, x1max);
  });
  i_dv += 1;
}

// Cartesian y coordinate
if (name.compare("coord_y") == 0) {
  if (derived_var.extent(4) <= 1)
    Kokkos::realloc(derived_var, nmb, n_dv, n3, n2, n1);
  auto dv = derived_var;
  int nx2 = indcs.nx2;
  par_for("coord_y", DevExeSpace(), 0, (nmb-1), ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    dv(m, i_dv, k, j, i) = CellCenterX(j-js, nx2, x2min, x2max);
  });
  i_dv += 1;
}

// Cartesian z coordinate
if (name.compare("coord_z") == 0) {
  if (derived_var.extent(4) <= 1)
    Kokkos::realloc(derived_var, nmb, n_dv, n3, n2, n1);
  auto dv = derived_var;
  int nx3 = indcs.nx3;
  par_for("coord_z", DevExeSpace(), 0, (nmb-1), ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    dv(m, i_dv, k, j, i) = CellCenterX(k-ks, nx3, x3min, x3max);
  });
  i_dv += 1;
}

// Spherical radius: r = sqrt(x² + y² + z²)
if (name.compare("coord_r") == 0) {
  if (derived_var.extent(4) <= 1)
    Kokkos::realloc(derived_var, nmb, n_dv, n3, n2, n1);
  auto dv = derived_var;
  int nx1 = indcs.nx1;
  int nx2 = indcs.nx2;
  int nx3 = indcs.nx3;
  par_for("coord_r", DevExeSpace(), 0, (nmb-1), ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real x = CellCenterX(i-is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    Real y = CellCenterX(j-js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    Real z = CellCenterX(k-ks, nx3, size.d_view(m).x3min, size.d_view(m).x3max);
    dv(m, i_dv, k, j, i) = sqrt(x*x + y*y + z*z);
  });
  i_dv += 1;
}

// Spherical polar angle: θ = acos(z/r), range [0, π]
if (name.compare("coord_theta") == 0) {
  if (derived_var.extent(4) <= 1)
    Kokkos::realloc(derived_var, nmb, n_dv, n3, n2, n1);
  auto dv = derived_var;
  int nx1 = indcs.nx1;
  int nx2 = indcs.nx2;
  int nx3 = indcs.nx3;
  par_for("coord_theta", DevExeSpace(), 0, (nmb-1), ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real x = CellCenterX(i-is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    Real y = CellCenterX(j-js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    Real z = CellCenterX(k-ks, nx3, size.d_view(m).x3min, size.d_view(m).x3max);
    Real r = sqrt(x*x + y*y + z*z);
    dv(m, i_dv, k, j, i) = (r > 0.0) ? acos(z / r) : 0.0;
  });
  i_dv += 1;
}

// Spherical/Cylindrical azimuthal angle: φ, range [0, 2π]
if (name.compare("coord_phi") == 0 || name.compare("coord_cyl_phi") == 0) {
  if (derived_var.extent(4) <= 1)
    Kokkos::realloc(derived_var, nmb, n_dv, n3, n2, n1);
  auto dv = derived_var;
  int nx1 = indcs.nx1;
  int nx2 = indcs.nx2;
  Real two_pi = 2.0 * M_PI;
  par_for("coord_phi", DevExeSpace(), 0, (nmb-1), ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real x = CellCenterX(i-is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    Real y = CellCenterX(j-js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    Real phi = atan2(y, x);
    if (phi < 0.0) phi += two_pi;  // shift from [-π,π] to [0,2π]
    dv(m, i_dv, k, j, i) = phi;
  });
  i_dv += 1;
}

// Cylindrical radius: R = sqrt(x² + y²)
if (name.compare("coord_cyl_R") == 0) {
  if (derived_var.extent(4) <= 1)
    Kokkos::realloc(derived_var, nmb, n_dv, n3, n2, n1);
  auto dv = derived_var;
  int nx1 = indcs.nx1;
  int nx2 = indcs.nx2;
  par_for("coord_cyl_R", DevExeSpace(), 0, (nmb-1), ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real x = CellCenterX(i-is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    Real y = CellCenterX(j-js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    dv(m, i_dv, k, j, i) = sqrt(x*x + y*y);
  });
  i_dv += 1;
}

// Cylindrical z (alias for Cartesian z, kept for user convenience)
if (name.compare("coord_cyl_z") == 0) {
  if (derived_var.extent(4) <= 1)
    Kokkos::realloc(derived_var, nmb, n_dv, n3, n2, n1);
  auto dv = derived_var;
  int nx3 = indcs.nx3;
  par_for("coord_cyl_z", DevExeSpace(), 0, (nmb-1), ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    dv(m, i_dv, k, j, i) = CellCenterX(k-ks, nx3, x3min, x3max);
  });
  i_dv += 1;
}
```

---

## Part 2: Generalize to N-D PDFs (up to 4D)

### 2.1 Modify `OutputParameters` in `outputs.hpp`

Replace scalar bin parameters with arrays. Keep old names for backward compatibility:

```cpp
struct OutputParameters {
  // ... existing fields ...

  // PDF parameters - generalized for N-D (max 4D)
  static const int MAX_PDF_DIM = 4;
  int pdf_ndim = 0;                              // actual number of dimensions (1-4)
  std::string pdf_variables[MAX_PDF_DIM];        // variable names for each dimension
  Real bin_min_arr[MAX_PDF_DIM];                 // min values for each dimension
  Real bin_max_arr[MAX_PDF_DIM];                 // max values for each dimension
  int nbins_arr[MAX_PDF_DIM] = {0,0,0,0};        // number of bins for each dimension
  bool logscale_arr[MAX_PDF_DIM] = {false,false,false,false}; // log scale flags

  // Keep old scalar fields for backward compatibility (aliases to arr[0], arr[1])
  // These are populated during parsing for legacy input files
  std::string variable;      // alias for pdf_variables[0]
  std::string variable_2;    // alias for pdf_variables[1]
  Real bin_min, bin_max;     // aliases for bin_min_arr[0], bin_max_arr[0]
  Real bin2_min, bin2_max;   // aliases for bin_min_arr[1], bin_max_arr[1]
  int nbin=0, nbin2=0;       // aliases for nbins_arr[0], nbins_arr[1]
  bool logscale=true, logscale2=true;  // aliases for logscale_arr[0], logscale_arr[1]

  bool mass_weighted=false;
  // ... rest of existing fields ...
};
```

### 2.2 Modify `PDFData` in `outputs.hpp`

Replace 2D-specific storage with N-D generalized storage:

```cpp
struct PDFData {
  static const int MAX_PDF_DIM = 4;

  int ndim;                                    // number of dimensions (1-4)
  int nbins[MAX_PDF_DIM];                      // bins per dimension
  int strides[MAX_PDF_DIM];                    // strides for flattened indexing (C-order)
  int total_bins;                              // product of all (nbins+2) with overflow bins

  // Bin edges for each dimension (size nbins[d]+1 each)
  Kokkos::View<Real*> bins[MAX_PDF_DIM];

  // Step sizes for fast bin calculation
  Real step_size[MAX_PDF_DIM];
  Real bin_min[MAX_PDF_DIM];                   // cached for kernel use
  bool logscale[MAX_PDF_DIM];

  bool bins_written;
  bool mass_weighted;

  // Flattened 1D result array (includes overflow bins: +2 per dimension)
  DvceArray1D<Real> result_;
  Kokkos::Experimental::ScatterView<Real*, LayoutWrapper> scatter_result;

  // Constructor
  PDFData(int ndim_, const int* nbins_, const Real* bin_min_, const Real* bin_max_,
          const bool* logscale_);
  ~PDFData() = default;
};
```

### 2.3 Modify `PDFOutput` class in `outputs.hpp`

Add member variable for reusable device array:

```cpp
class PDFOutput : public BaseTypeOutput {
 public:
  PDFOutput(ParameterInput *pin, Mesh *pm, OutputParameters oparams);

  PDFData pdf_data;

  // Reusable device array for output variables (avoids repeated allocation)
  DvceArray5D<Real> outvars_device_;
  bool outvars_device_allocated_ = false;

  void LoadOutputData(Mesh *pm) override;
  void WriteOutputFile(Mesh *pm, ParameterInput *pin) override;
};
```

### 2.4 Update Input Parsing in `outputs.cpp`

Parse generalized N-D input with backward compatibility:

```cpp
// In the PDF parsing section of Outputs::Outputs():

if (opar.file_type.compare("pdf") == 0) {
  // Determine number of dimensions by checking which variable_N parameters exist
  opar.pdf_ndim = 0;

  // Check for new-style input: variable_1, variable_2, variable_3, variable_4
  for (int d = 0; d < OutputParameters::MAX_PDF_DIM; ++d) {
    std::string var_key = "variable_" + std::to_string(d + 1);
    std::string min_key = "bin" + std::to_string(d + 1) + "_min";
    std::string max_key = "bin" + std::to_string(d + 1) + "_max";
    std::string nbin_key = "nbin" + std::to_string(d + 1);
    std::string log_key = "logscale" + std::to_string(d + 1);

    if (pin->DoesParameterExist(opar.block_name, var_key)) {
      opar.pdf_variables[d] = pin->GetString(opar.block_name, var_key);
      opar.bin_min_arr[d] = pin->GetReal(opar.block_name, min_key);
      opar.bin_max_arr[d] = pin->GetReal(opar.block_name, max_key);
      opar.nbins_arr[d] = pin->GetInteger(opar.block_name, nbin_key);
      opar.logscale_arr[d] = pin->GetOrAddBoolean(opar.block_name, log_key, false);
      opar.pdf_ndim = d + 1;
    }
  }

  // Backward compatibility: check for old-style input (variable, variable_2, etc.)
  if (opar.pdf_ndim == 0) {
    // Old style: "variable" is dimension 1
    if (pin->DoesParameterExist(opar.block_name, "variable")) {
      opar.pdf_variables[0] = pin->GetString(opar.block_name, "variable");
      opar.bin_min_arr[0] = pin->GetReal(opar.block_name, "bin_min");
      opar.bin_max_arr[0] = pin->GetReal(opar.block_name, "bin_max");
      opar.nbins_arr[0] = pin->GetInteger(opar.block_name, "nbin");
      opar.logscale_arr[0] = pin->GetOrAddBoolean(opar.block_name, "logscale", false);
      opar.pdf_ndim = 1;

      // Old style: "variable_2" is dimension 2 (only if nbin2 > 0)
      if (pin->DoesParameterExist(opar.block_name, "variable_2")) {
        int nbin2_val = pin->GetOrAddInteger(opar.block_name, "nbin2", 0);
        if (nbin2_val > 0) {
          opar.pdf_variables[1] = pin->GetString(opar.block_name, "variable_2");
          opar.bin_min_arr[1] = pin->GetReal(opar.block_name, "bin2_min");
          opar.bin_max_arr[1] = pin->GetReal(opar.block_name, "bin2_max");
          opar.nbins_arr[1] = nbin2_val;
          opar.logscale_arr[1] = pin->GetOrAddBoolean(opar.block_name, "logscale2", false);
          opar.pdf_ndim = 2;
        }
      }
    }
  }

  // Populate legacy aliases for code that still uses them
  opar.variable = opar.pdf_variables[0];
  opar.variable_2 = (opar.pdf_ndim >= 2) ? opar.pdf_variables[1] : "";
  opar.bin_min = opar.bin_min_arr[0];
  opar.bin_max = opar.bin_max_arr[0];
  opar.nbin = opar.nbins_arr[0];
  opar.logscale = opar.logscale_arr[0];
  if (opar.pdf_ndim >= 2) {
    opar.bin2_min = opar.bin_min_arr[1];
    opar.bin2_max = opar.bin_max_arr[1];
    opar.nbin2 = opar.nbins_arr[1];
    opar.logscale2 = opar.logscale_arr[1];
  }

  opar.mass_weighted = pin->GetOrAddBoolean(opar.block_name, "mass_weighted", false);

  // Validate
  if (opar.pdf_ndim == 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "PDF output requires at least one variable" << std::endl;
    exit(EXIT_FAILURE);
  }
}
```

### 2.5 Update `PDFData` Constructor

C-style row-major ordering: last dimension has stride 1, first dimension has largest stride.

```cpp
PDFData::PDFData(int ndim_, const int* nbins_, const Real* bin_min_,
                 const Real* bin_max_, const bool* logscale_)
    : ndim(ndim_), bins_written(false), mass_weighted(false) {

  // Initialize all dimensions
  for (int d = 0; d < MAX_PDF_DIM; ++d) {
    nbins[d] = 0;
    strides[d] = 0;
    bin_min[d] = 0.0;
    step_size[d] = 0.0;
    logscale[d] = false;
  }

  // Calculate strides in C-order (last dimension varies fastest)
  // For shape [n0, n1, n2, n3]: strides = [n1*n2*n3, n2*n3, n3, 1]
  // With overflow bins: each dimension has size (nbins[d]+2)

  // First pass: store nbins and compute bin edges
  for (int d = 0; d < ndim; ++d) {
    nbins[d] = nbins_[d];
    bin_min[d] = bin_min_[d];
    logscale[d] = logscale_[d];

    // Calculate step size
    if (logscale[d]) {
      step_size[d] = (std::log10(bin_max_[d]) - std::log10(bin_min_[d])) / nbins[d];
    } else {
      step_size[d] = (bin_max_[d] - bin_min[d]) / nbins[d];
    }

    // Allocate and populate bins for this dimension
    bins[d] = Kokkos::View<Real*>("bins_" + std::to_string(d), nbins[d] + 1);
    auto bins_host = Kokkos::create_mirror_view(bins[d]);

    if (logscale[d]) {
      Real logmin = std::log10(bin_min_[d]);
      Real logmax = std::log10(bin_max_[d]);
      for (int i = 0; i <= nbins[d]; ++i) {
        bins_host(i) = std::pow(10.0, logmin + i * (logmax - logmin) / nbins[d]);
      }
    } else {
      Real step = (bin_max_[d] - bin_min_[d]) / nbins[d];
      for (int i = 0; i <= nbins[d]; ++i) {
        bins_host(i) = bin_min_[d] + i * step;
      }
    }
    Kokkos::deep_copy(bins[d], bins_host);
  }

  // Second pass: compute strides in C-order
  // stride[ndim-1] = 1, stride[d] = stride[d+1] * (nbins[d+1]+2)
  strides[ndim - 1] = 1;
  for (int d = ndim - 2; d >= 0; --d) {
    strides[d] = strides[d + 1] * (nbins[d + 1] + 2);
  }

  // Total size
  total_bins = strides[0] * (nbins[0] + 2);

  Kokkos::fence();
}
```

### 2.6 Update `PDFOutput` Constructor in `pdf.cpp`

```cpp
PDFOutput::PDFOutput(ParameterInput *pin, Mesh *pm, OutputParameters op)
    : BaseTypeOutput(pin, pm, op),
      pdf_data(op.pdf_ndim, op.nbins_arr, op.bin_min_arr, op.bin_max_arr, op.logscale_arr) {

  // Create directory for outputs
  std::string dir_name = "pdf_" + op.file_id;
  mkdir(dir_name.c_str(), 0775);

  pdf_data.mass_weighted = op.mass_weighted;

  // Validate logscale with positive bin_min
  for (int d = 0; d < op.pdf_ndim; ++d) {
    if (op.logscale_arr[d] && op.bin_min_arr[d] <= 0.0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "logscale" << (d+1) << " is true but bin"
                << (d+1) << "_min <= 0.0" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  // Allocate flattened result array
  pdf_data.result_ = DvceArray1D<Real>("pdf_result", pdf_data.total_bins);
  pdf_data.scatter_result = Kokkos::Experimental::ScatterView<Real*, LayoutWrapper>(
      pdf_data.result_);
}
```

### 2.7 Update `PDFOutput::LoadOutputData()` in `pdf.cpp`

Inlined bin computation (no nested lambdas for CUDA compatibility):

```cpp
void PDFOutput::LoadOutputData(Mesh *pm) {
  // Compute derived variables for all dimensions
  for (int d = 0; d < out_params.pdf_ndim; ++d) {
    if (out_params.contains_derived) {
      ComputeDerivedVariable(out_params.pdf_variables[d], pm);
    }
  }

  // Get physics pointer for mass weighting
  DvceArray5D<Real> *u0_ptr = nullptr;
  if (pm->pmb_pack->phydro != nullptr) {
    u0_ptr = &(pm->pmb_pack->phydro->u0);
  } else if (pm->pmb_pack->pmhd != nullptr) {
    u0_ptr = &(pm->pmb_pack->pmhd->u0);
  }
  if (u0_ptr == nullptr && pdf_data.mass_weighted) {
    std::cout << "### FATAL ERROR: mass_weighted PDF requires hydro or MHD" << std::endl;
    exit(EXIT_FAILURE);
  }

  auto &size = pm->pmb_pack->pmb->mb_size;
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb = pm->pmb_pack->nmb_thispack;
  int nx1 = indcs.nx1 + 2*indcs.ng;
  int nx2 = indcs.nx2 + 2*indcs.ng;
  int nx3 = indcs.nx3 + 2*indcs.ng;

  // Allocate/resize reusable device array for output variables
  int ndim = out_params.pdf_ndim;
  if (!outvars_device_allocated_ ||
      outvars_device_.extent(0) != static_cast<size_t>(ndim) ||
      outvars_device_.extent(1) != static_cast<size_t>(nmb)) {
    outvars_device_ = DvceArray5D<Real>("outvars_device", ndim, nmb, nx3, nx2, nx1);
    outvars_device_allocated_ = true;
  }

  // Copy output variable data to device
  for (int d = 0; d < ndim; ++d) {
    auto d_slice = Kokkos::subview(*(outvars[d].data_ptr),
        Kokkos::ALL(), outvars[d].data_index, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
    auto d_target_slice = Kokkos::subview(outvars_device_, d,
        Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
    Kokkos::deep_copy(d_target_slice, d_slice);
  }
  Kokkos::fence();

  // Reset scatter view
  pdf_data.scatter_result.reset();
  Kokkos::deep_copy(pdf_data.result_, 0);
  Kokkos::fence();

  // Capture references for kernel
  auto result = pdf_data.result_;
  auto scatter = pdf_data.scatter_result;
  auto outvars_dev = outvars_device_;
  bool mass_weighted = pdf_data.mass_weighted;

  // Capture per-dimension parameters as scalars (required for GPU)
  int ndim_cap = ndim;
  int nbin0 = pdf_data.nbins[0], nbin1 = pdf_data.nbins[1];
  int nbin2 = pdf_data.nbins[2], nbin3 = pdf_data.nbins[3];
  int stride0 = pdf_data.strides[0], stride1 = pdf_data.strides[1];
  int stride2 = pdf_data.strides[2], stride3 = pdf_data.strides[3];
  Real step0 = pdf_data.step_size[0], step1 = pdf_data.step_size[1];
  Real step2 = pdf_data.step_size[2], step3 = pdf_data.step_size[3];
  Real min0 = pdf_data.bin_min[0], min1 = pdf_data.bin_min[1];
  Real min2 = pdf_data.bin_min[2], min3 = pdf_data.bin_min[3];
  bool log0 = pdf_data.logscale[0], log1 = pdf_data.logscale[1];
  bool log2 = pdf_data.logscale[2], log3 = pdf_data.logscale[3];
  auto bins0 = pdf_data.bins[0], bins1 = pdf_data.bins[1];
  auto bins2 = pdf_data.bins[2], bins3 = pdf_data.bins[3];

  // Capture u0 only if needed
  DvceArray5D<Real> u0_local;
  if (u0_ptr != nullptr) {
    u0_local = *u0_ptr;
  }

  par_for("pdf_nd", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    int global_idx = 0;

    // Dimension 0 (always present)
    {
      Real val = outvars_dev(0, m, k, j, i);
      int bin_idx;
      if (val < bins0(0)) {
        bin_idx = 0;  // underflow
      } else if (val >= bins0(nbin0)) {
        bin_idx = nbin0 + 1;  // overflow
      } else if (log0) {
        bin_idx = static_cast<int>(log10(val / min0) / step0) + 1;
      } else {
        bin_idx = static_cast<int>((val - min0) / step0) + 1;
      }
      global_idx += bin_idx * stride0;
    }

    // Dimension 1
    if (ndim_cap >= 2) {
      Real val = outvars_dev(1, m, k, j, i);
      int bin_idx;
      if (val < bins1(0)) {
        bin_idx = 0;
      } else if (val >= bins1(nbin1)) {
        bin_idx = nbin1 + 1;
      } else if (log1) {
        bin_idx = static_cast<int>(log10(val / min1) / step1) + 1;
      } else {
        bin_idx = static_cast<int>((val - min1) / step1) + 1;
      }
      global_idx += bin_idx * stride1;
    }

    // Dimension 2
    if (ndim_cap >= 3) {
      Real val = outvars_dev(2, m, k, j, i);
      int bin_idx;
      if (val < bins2(0)) {
        bin_idx = 0;
      } else if (val >= bins2(nbin2)) {
        bin_idx = nbin2 + 1;
      } else if (log2) {
        bin_idx = static_cast<int>(log10(val / min2) / step2) + 1;
      } else {
        bin_idx = static_cast<int>((val - min2) / step2) + 1;
      }
      global_idx += bin_idx * stride2;
    }

    // Dimension 3
    if (ndim_cap >= 4) {
      Real val = outvars_dev(3, m, k, j, i);
      int bin_idx;
      if (val < bins3(0)) {
        bin_idx = 0;
      } else if (val >= bins3(nbin3)) {
        bin_idx = nbin3 + 1;
      } else if (log3) {
        bin_idx = static_cast<int>(log10(val / min3) / step3) + 1;
      } else {
        bin_idx = static_cast<int>((val - min3) / step3) + 1;
      }
      global_idx += bin_idx * stride3;
    }

    // Compute weight (volume, optionally mass-weighted)
    Real weight = size.d_view(m).dx1 * size.d_view(m).dx2 * size.d_view(m).dx3;
    if (mass_weighted) {
      weight *= u0_local(m, IDN, k, j, i);
    }

    // Atomic add to histogram
    auto res = scatter.access();
    res(global_idx) += weight;
  });

  Kokkos::Experimental::contribute(result, scatter);
  Kokkos::fence();

  // MPI reduction (works with GPU-aware MPI on Frontier)
#if MPI_PARALLEL_ENABLED
  if (global_variable::my_rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, result.data(), pdf_data.total_bins,
               MPI_ATHENA_REAL, MPI_SUM, 0, MPI_COMM_WORLD);
  } else {
    MPI_Reduce(result.data(), result.data(), pdf_data.total_bins,
               MPI_ATHENA_REAL, MPI_SUM, 0, MPI_COMM_WORLD);
  }
#endif
}
```

### 2.8 Update `PDFOutput::WriteOutputFile()` in `pdf.cpp`

Binary output format with ASCII header file:

```cpp
void PDFOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  // Only rank 0 writes
  if (global_variable::my_rank != 0) {
    out_params.file_number++;
    if (out_params.last_time < 0.0) {
      out_params.last_time = pm->time;
    } else {
      out_params.last_time += out_params.dt;
    }
    return;
  }

  std::string dir_name = "pdf_" + out_params.file_id;

  // Write header file (once) - ASCII with all metadata
  if (!pdf_data.bins_written) {
    std::string header_fname = dir_name + "/" + out_params.file_basename + ".header";
    FILE *hfile = std::fopen(header_fname.c_str(), "w");
    if (hfile == nullptr) {
      std::cout << "### FATAL ERROR: Could not open " << header_fname << std::endl;
      exit(EXIT_FAILURE);
    }

    std::fprintf(hfile, "# AthenaK PDF Header\n");
    std::fprintf(hfile, "ndim = %d\n", out_params.pdf_ndim);
    std::fprintf(hfile, "total_bins = %d\n", pdf_data.total_bins);
    std::fprintf(hfile, "mass_weighted = %d\n", pdf_data.mass_weighted ? 1 : 0);
    std::fprintf(hfile, "precision = %s\n", sizeof(Real) == 4 ? "float32" : "float64");
    std::fprintf(hfile, "byte_order = %s\n", IsBigEndian() ? "big" : "little");

    // Shape and strides
    std::fprintf(hfile, "shape =");
    for (int d = 0; d < out_params.pdf_ndim; ++d) {
      std::fprintf(hfile, " %d", pdf_data.nbins[d] + 2);
    }
    std::fprintf(hfile, "\n");

    std::fprintf(hfile, "strides =");
    for (int d = 0; d < out_params.pdf_ndim; ++d) {
      std::fprintf(hfile, " %d", pdf_data.strides[d]);
    }
    std::fprintf(hfile, "\n");

    // Per-dimension info
    for (int d = 0; d < out_params.pdf_ndim; ++d) {
      std::fprintf(hfile, "\n[dim%d]\n", d);
      std::fprintf(hfile, "variable = %s\n", out_params.pdf_variables[d].c_str());
      std::fprintf(hfile, "nbins = %d\n", pdf_data.nbins[d]);
      std::fprintf(hfile, "logscale = %d\n", pdf_data.logscale[d] ? 1 : 0);
      std::fprintf(hfile, "bin_min = %.15e\n", pdf_data.bin_min[d]);
      std::fprintf(hfile, "step_size = %.15e\n", pdf_data.step_size[d]);

      // Write bin edges
      auto bins_host = Kokkos::create_mirror_view(pdf_data.bins[d]);
      Kokkos::deep_copy(bins_host, pdf_data.bins[d]);
      Kokkos::fence();

      std::fprintf(hfile, "bin_edges =");
      for (int n = 0; n <= pdf_data.nbins[d]; ++n) {
        std::fprintf(hfile, " %.15e", bins_host(n));
      }
      std::fprintf(hfile, "\n");
    }

    std::fclose(hfile);
    pdf_data.bins_written = true;
  }

  // Write binary data file
  char number[6];
  std::snprintf(number, sizeof(number), "%05d", out_params.file_number);
  std::string data_fname = dir_name + "/" + out_params.file_basename + "."
                         + std::string(number) + ".bin";

  FILE *pfile = std::fopen(data_fname.c_str(), "wb");
  if (pfile == nullptr) {
    std::cout << "### FATAL ERROR: Could not open " << data_fname << std::endl;
    exit(EXIT_FAILURE);
  }

  // Write time as first value
  Real time_val = pm->time;
  std::fwrite(&time_val, sizeof(Real), 1, pfile);

  // Copy result to host and write
  auto result_host = Kokkos::create_mirror_view(pdf_data.result_);
  Kokkos::deep_copy(result_host, pdf_data.result_);
  Kokkos::fence();

  std::fwrite(result_host.data(), sizeof(Real), pdf_data.total_bins, pfile);
  std::fclose(pfile);

  // Update counters
  out_params.file_number++;
  if (out_params.last_time < 0.0) {
    out_params.last_time = pm->time;
  } else {
    out_params.last_time += out_params.dt;
  }
  pin->SetInteger(out_params.block_name, "file_number", out_params.file_number);
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
}
```

---

## Part 3: Update Variable Registration in `basetype_output.cpp`

Modify the variable loading loop to handle N-D PDFs:

```cpp
// In BaseTypeOutput::BaseTypeOutput(), update the PDF variable vector construction:

std::vector<std::string> variables;
variables.push_back(out_params.variable);

if (out_params.file_type == "pdf") {
  // Clear and rebuild with all PDF dimensions
  variables.clear();
  for (int d = 0; d < out_params.pdf_ndim; ++d) {
    variables.push_back(out_params.pdf_variables[d]);
  }
}

// The existing for-loop over variables handles the rest
for (const auto& variable : variables) {
  // ... existing variable matching code ...
}
```

---

## Example Input File Formats

### New Style (4D PDF example):

```
<output3>
file_type   = pdf
variable_1  = coord_r
bin1_min    = 0.1
bin1_max    = 10.0
nbin1       = 50
logscale1   = true

variable_2  = coord_theta
bin2_min    = 0.0
bin2_max    = 3.14159265359
nbin2       = 30
logscale2   = false

variable_3  = hydro_w_d
bin3_min    = 1e-6
bin3_max    = 1e2
nbin3       = 100
logscale3   = true

variable_4  = temperature
bin4_min    = 1e4
bin4_max    = 1e8
nbin4       = 80
logscale4   = true

mass_weighted = true
dt          = 0.1
```

### Old Style (backward compatible 2D PDF):

```
<output2>
file_type   = pdf
variable    = hydro_w_d
bin_min     = 1e-6
bin_max     = 1e2
nbin        = 100
logscale    = true

variable_2  = temperature
bin2_min    = 1e4
bin2_max    = 1e8
nbin2       = 80
logscale2   = true

mass_weighted = false
dt          = 0.1
```

---

## Summary of Files to Modify

| File | Changes |
|------|---------|
| `src/outputs/outputs.hpp` | Add coordinate vars to `var_choice[]`, update `NOUTPUT_CHOICES`, modify `OutputParameters` and `PDFData` structs, add `outvars_device_` member to `PDFOutput` |
| `src/outputs/outputs.cpp` | Update PDF input parsing for N-D with backward compatibility |
| `src/outputs/basetype_output.cpp` | Add coordinate variable registration, update PDF variable loading loop |
| `src/outputs/derived_variables.cpp` | Add computation kernels for all coordinate variables |
| `src/outputs/pdf.cpp` | Rewrite constructor, `LoadOutputData()`, and `WriteOutputFile()` for N-D with binary output |

---

## Testing Plan

1. **Unit tests for coordinates**: Create a simple test problem, output `coord_x`, `coord_y`, `coord_z`, `coord_r` as derived variables, verify values match expected positions
2. **1D PDF regression**: Run existing 1D PDF test with old-style input, verify identical output
3. **2D PDF regression**: Run existing 2D PDF test with old-style input, verify identical output
4. **2D PDF new style**: Same test with new-style `variable_1`/`variable_2` input
5. **3D PDF**: Create test with 3 variables, verify output shape and values
6. **4D PDF**: Create test with 4 variables (e.g., coord_r, coord_theta, density, temperature)
7. **Coordinate binning**: Test `coord_r` binning produces expected radial profiles
8. **MPI correctness**: Run multi-rank test, verify reduction is correct

---

## Python Reader Example

```python
import numpy as np
import struct

def read_pdf_header(header_file):
    """Read PDF header file and return metadata dict."""
    meta = {'dims': []}
    current_dim = None

    with open(header_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            if line.startswith('[dim'):
                current_dim = int(line[4:-1])
                meta['dims'].append({})
                continue

            if '=' in line:
                key, val = line.split('=', 1)
                key = key.strip()
                val = val.strip()

                if current_dim is not None:
                    # Per-dimension data
                    if key == 'bin_edges':
                        meta['dims'][current_dim][key] = np.array([float(x) for x in val.split()])
                    elif key in ('nbins', 'logscale'):
                        meta['dims'][current_dim][key] = int(val)
                    elif key == 'variable':
                        meta['dims'][current_dim][key] = val
                    else:
                        meta['dims'][current_dim][key] = float(val)
                else:
                    # Global metadata
                    if key in ('ndim', 'total_bins', 'mass_weighted'):
                        meta[key] = int(val)
                    elif key == 'shape':
                        meta[key] = [int(x) for x in val.split()]
                    elif key == 'strides':
                        meta[key] = [int(x) for x in val.split()]
                    elif key in ('precision', 'byte_order'):
                        meta[key] = val

    return meta


def read_pdf_data(data_file, meta):
    """Read binary PDF data file and return time + reshaped array."""
    dtype = np.float32 if meta['precision'] == 'float32' else np.float64

    with open(data_file, 'rb') as f:
        time = np.frombuffer(f.read(dtype().nbytes), dtype=dtype)[0]
        data = np.frombuffer(f.read(), dtype=dtype)

    # Reshape to N-D array
    data = data.reshape(meta['shape'])

    return time, data


def read_pdf(basename, file_number):
    """
    Read PDF output files.

    Parameters
    ----------
    basename : str
        Path to PDF directory + base filename (e.g., "pdf_rho/turb")
    file_number : int
        Output file number

    Returns
    -------
    dict with keys: 'time', 'data', 'meta', 'bin_centers'
    """
    import os
    dir_name = os.path.dirname(basename)
    base = os.path.basename(basename)

    header_file = f"{dir_name}/{base}.header"
    data_file = f"{dir_name}/{base}.{file_number:05d}.bin"

    meta = read_pdf_header(header_file)
    time, data = read_pdf_data(data_file, meta)

    # Compute bin centers for convenience (excluding under/overflow)
    bin_centers = []
    for d in range(meta['ndim']):
        edges = meta['dims'][d]['bin_edges']
        centers = 0.5 * (edges[:-1] + edges[1:])
        bin_centers.append(centers)

    return {
        'time': time,
        'data': data,
        'meta': meta,
        'bin_centers': bin_centers
    }


# Example usage:
# result = read_pdf("pdf_density/turb", 0)
# print(f"Time: {result['time']}")
# print(f"Shape: {result['data'].shape}")
# # Access inner bins (excluding overflow): data[1:-1, 1:-1, ...]
```

---

## Part 3: Mass and Energy Flux Derived Variables

### 3.1 Overview

Add derived variables for computing mass and energy fluxes, designed to be binned by radius (spherical) or height (vertical) using the PDF capabilities. These variables output the flux contribution per cell (flux × volume), which when summed in radial/vertical bins and divided by the bin width gives the physical flux.

**Spherical fluxes** (for radial profiles):
- `mdot_sph` — mass flux: ρ v_r dV
- `mdot_sph_out` — outward mass flux only (v_r > 0)
- `mdot_sph_in` — inward mass flux only (v_r < 0)
- `edot_sph` — energy flux: F_E,r dV
- `edot_sph_out` — outward energy flux only (v_r > 0)
- `edot_sph_in` — inward energy flux only (v_r < 0)

**Vertical fluxes** (for height profiles):
- `mdot_vert` — mass flux: ρ v_z sign(z) dV
- `mdot_vert_out` — outward mass flux only (v_z sign(z) > 0)
- `mdot_vert_in` — inward mass flux only (v_z sign(z) < 0)
- `edot_vert` — energy flux: F_E,z sign(z) dV
- `edot_vert_out` — outward energy flux only (v_z sign(z) > 0)
- `edot_vert_in` — inward energy flux only (v_z sign(z) < 0)

### 3.2 Physics Definitions

**Radial velocity:**
```
v_r = (v_x * x + v_y * y + v_z * z) / r
```
where `r = sqrt(x² + y² + z²)`

**Radial magnetic field (for MHD):**
```
B_r = (B_x * x + B_y * y + B_z * z) / r
```

**Mass flux contribution (spherical):**
```
mdot_sph = ρ * v_r * dV
```

**Energy flux (MHD):**
```
F_E,r = [γP/(γ-1) + ρv²/2 + B²] * v_r - B_r(v·B)
```
where:
- First term: enthalpy flux (thermal + kinetic + magnetic pressure)
- Second term: Poynting flux contribution

**Energy flux (Hydro):**
```
F_E,r = [γP/(γ-1) + ρv²/2] * v_r
```

**Energy flux contribution:**
```
edot_sph = F_E,r * dV
```

**Vertical fluxes:**
Same formulas but using v_z and B_z directly, with a sign(z) factor to make "outward" mean "away from midplane":
```
mdot_vert = ρ * v_z * sign(z) * dV
edot_vert = F_E,z * sign(z) * dV
```

### 3.3 Add to `var_choice[]` in `outputs.hpp`

```cpp
// Mass and energy flux derived variables (add after coordinate variables)
"mdot_sph", "mdot_sph_out", "mdot_sph_in",
"edot_sph", "edot_sph_out", "edot_sph_in",
"mdot_vert", "mdot_vert_out", "mdot_vert_in",
"edot_vert", "edot_vert_out", "edot_vert_in"
```

Update `NOUTPUT_CHOICES`: 167 + 12 → 179.

### 3.4 Register Flux Variables in `basetype_output.cpp`

```cpp
// Spherical mass flux
if (variable.compare("mdot_sph") == 0) {
  out_params.contains_derived = true;
  out_params.n_derived += 1;
  outvars.emplace_back("mdot_sph", out_params.i_derived, &derived_var);
}
if (variable.compare("mdot_sph_out") == 0) {
  out_params.contains_derived = true;
  out_params.n_derived += 1;
  outvars.emplace_back("mdot_sph_out", out_params.i_derived, &derived_var);
}
if (variable.compare("mdot_sph_in") == 0) {
  out_params.contains_derived = true;
  out_params.n_derived += 1;
  outvars.emplace_back("mdot_sph_in", out_params.i_derived, &derived_var);
}

// Spherical energy flux
if (variable.compare("edot_sph") == 0) {
  out_params.contains_derived = true;
  out_params.n_derived += 1;
  outvars.emplace_back("edot_sph", out_params.i_derived, &derived_var);
}
if (variable.compare("edot_sph_out") == 0) {
  out_params.contains_derived = true;
  out_params.n_derived += 1;
  outvars.emplace_back("edot_sph_out", out_params.i_derived, &derived_var);
}
if (variable.compare("edot_sph_in") == 0) {
  out_params.contains_derived = true;
  out_params.n_derived += 1;
  outvars.emplace_back("edot_sph_in", out_params.i_derived, &derived_var);
}

// Vertical mass flux
if (variable.compare("mdot_vert") == 0) {
  out_params.contains_derived = true;
  out_params.n_derived += 1;
  outvars.emplace_back("mdot_vert", out_params.i_derived, &derived_var);
}
if (variable.compare("mdot_vert_out") == 0) {
  out_params.contains_derived = true;
  out_params.n_derived += 1;
  outvars.emplace_back("mdot_vert_out", out_params.i_derived, &derived_var);
}
if (variable.compare("mdot_vert_in") == 0) {
  out_params.contains_derived = true;
  out_params.n_derived += 1;
  outvars.emplace_back("mdot_vert_in", out_params.i_derived, &derived_var);
}

// Vertical energy flux
if (variable.compare("edot_vert") == 0) {
  out_params.contains_derived = true;
  out_params.n_derived += 1;
  outvars.emplace_back("edot_vert", out_params.i_derived, &derived_var);
}
if (variable.compare("edot_vert_out") == 0) {
  out_params.contains_derived = true;
  out_params.n_derived += 1;
  outvars.emplace_back("edot_vert_out", out_params.i_derived, &derived_var);
}
if (variable.compare("edot_vert_in") == 0) {
  out_params.contains_derived = true;
  out_params.n_derived += 1;
  outvars.emplace_back("edot_vert_in", out_params.i_derived, &derived_var);
}
```

### 3.5 Compute Flux Variables in `derived_variables.cpp`

#### Spherical Mass Flux (mdot_sph, mdot_sph_out, mdot_sph_in)

```cpp
// Spherical mass flux: mdot_sph = ρ * v_r * dV
if (name.compare("mdot_sph") == 0 ||
    name.compare("mdot_sph_out") == 0 ||
    name.compare("mdot_sph_in") == 0) {
  if (derived_var.extent(4) <= 1)
    Kokkos::realloc(derived_var, nmb, n_dv, n3, n2, n1);
  auto dv = derived_var;
  int nx1 = indcs.nx1;
  int nx2 = indcs.nx2;
  int nx3 = indcs.nx3;

  // Determine which variant
  bool all_flux = (name.compare("mdot_sph") == 0);
  bool out_only = (name.compare("mdot_sph_out") == 0);
  // in_only is the remaining case

  // Get velocity data (works for both hydro and MHD)
  DvceArray5D<Real> u0;
  if (pmbp->phydro != nullptr) {
    u0 = pmbp->phydro->u0;
  } else if (pmbp->pmhd != nullptr) {
    u0 = pmbp->pmhd->u0;
  }

  par_for("mdot_sph", DevExeSpace(), 0, (nmb-1), ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    // Get cell center coordinates
    Real x = CellCenterX(i-is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    Real y = CellCenterX(j-js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    Real z = CellCenterX(k-ks, nx3, size.d_view(m).x3min, size.d_view(m).x3max);
    Real r = sqrt(x*x + y*y + z*z);

    // Get density and velocity
    Real rho = u0(m, IDN, k, j, i);
    Real vx = u0(m, IM1, k, j, i) / rho;
    Real vy = u0(m, IM2, k, j, i) / rho;
    Real vz = u0(m, IM3, k, j, i) / rho;

    // Radial velocity
    Real v_r = (r > 0.0) ? (vx*x + vy*y + vz*z) / r : 0.0;

    // Cell volume
    Real dV = size.d_view(m).dx1 * size.d_view(m).dx2 * size.d_view(m).dx3;

    // Compute flux contribution based on variant
    Real mdot = rho * v_r * dV;

    if (all_flux) {
      dv(m, i_dv, k, j, i) = mdot;
    } else if (out_only) {
      dv(m, i_dv, k, j, i) = (v_r > 0.0) ? mdot : 0.0;
    } else {  // in_only
      dv(m, i_dv, k, j, i) = (v_r < 0.0) ? mdot : 0.0;
    }
  });
  i_dv += 1;
}
```

#### Spherical Energy Flux (edot_sph, edot_sph_out, edot_sph_in)

```cpp
// Spherical energy flux: edot_sph = F_E,r * dV
if (name.compare("edot_sph") == 0 ||
    name.compare("edot_sph_out") == 0 ||
    name.compare("edot_sph_in") == 0) {
  if (derived_var.extent(4) <= 1)
    Kokkos::realloc(derived_var, nmb, n_dv, n3, n2, n1);
  auto dv = derived_var;
  int nx1 = indcs.nx1;
  int nx2 = indcs.nx2;
  int nx3 = indcs.nx3;

  bool all_flux = (name.compare("edot_sph") == 0);
  bool out_only = (name.compare("edot_sph_out") == 0);

  bool is_mhd = (pmbp->pmhd != nullptr);

  // Get data arrays and EOS
  DvceArray5D<Real> u0;
  DvceFaceFld4D<Real> b0;
  Real gamma;

  if (is_mhd) {
    u0 = pmbp->pmhd->u0;
    b0 = pmbp->pmhd->b0;
    gamma = pmbp->pmhd->peos->eos_data.gamma;
  } else {
    u0 = pmbp->phydro->u0;
    gamma = pmbp->phydro->peos->eos_data.gamma;
  }

  Real gm1 = gamma - 1.0;

  par_for("edot_sph", DevExeSpace(), 0, (nmb-1), ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    // Get cell center coordinates
    Real x = CellCenterX(i-is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    Real y = CellCenterX(j-js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    Real z = CellCenterX(k-ks, nx3, size.d_view(m).x3min, size.d_view(m).x3max);
    Real r = sqrt(x*x + y*y + z*z);

    // Get primitives
    Real rho = u0(m, IDN, k, j, i);
    Real vx = u0(m, IM1, k, j, i) / rho;
    Real vy = u0(m, IM2, k, j, i) / rho;
    Real vz = u0(m, IM3, k, j, i) / rho;
    Real v2 = vx*vx + vy*vy + vz*vz;

    // Radial velocity
    Real v_r = (r > 0.0) ? (vx*x + vy*y + vz*z) / r : 0.0;

    // Cell volume
    Real dV = size.d_view(m).dx1 * size.d_view(m).dx2 * size.d_view(m).dx3;

    Real F_Er;
    if (is_mhd) {
      // Cell-centered B (average of face values)
      Real Bx = 0.5*(b0.x1f(m,k,j,i) + b0.x1f(m,k,j,i+1));
      Real By = 0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i));
      Real Bz = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
      Real B2 = Bx*Bx + By*By + Bz*Bz;

      // Radial B
      Real B_r = (r > 0.0) ? (Bx*x + By*y + Bz*z) / r : 0.0;

      // v · B
      Real vdotB = vx*Bx + vy*By + vz*Bz;

      // Total energy E = P/(γ-1) + ρv²/2 + B²/2
      Real E = u0(m, IEN, k, j, i);

      // Pressure from E: P = (γ-1) * (E - ρv²/2 - B²/2)
      Real P = gm1 * (E - 0.5*rho*v2 - 0.5*B2);

      // Energy flux: F_E,r = [γP/(γ-1) + ρv²/2 + B²] * v_r - B_r(v·B)
      F_Er = (gamma/gm1*P + 0.5*rho*v2 + B2) * v_r - B_r * vdotB;
    } else {
      // Hydro case
      Real E = u0(m, IEN, k, j, i);
      Real P = gm1 * (E - 0.5*rho*v2);

      // Energy flux: F_E,r = [γP/(γ-1) + ρv²/2] * v_r
      F_Er = (gamma/gm1*P + 0.5*rho*v2) * v_r;
    }

    Real edot = F_Er * dV;

    if (all_flux) {
      dv(m, i_dv, k, j, i) = edot;
    } else if (out_only) {
      dv(m, i_dv, k, j, i) = (v_r > 0.0) ? edot : 0.0;
    } else {  // in_only
      dv(m, i_dv, k, j, i) = (v_r < 0.0) ? edot : 0.0;
    }
  });
  i_dv += 1;
}
```

#### Vertical Mass Flux (mdot_vert, mdot_vert_out, mdot_vert_in)

```cpp
// Vertical mass flux: mdot_vert = ρ * v_z * sign(z) * dV
if (name.compare("mdot_vert") == 0 ||
    name.compare("mdot_vert_out") == 0 ||
    name.compare("mdot_vert_in") == 0) {
  if (derived_var.extent(4) <= 1)
    Kokkos::realloc(derived_var, nmb, n_dv, n3, n2, n1);
  auto dv = derived_var;
  int nx3 = indcs.nx3;

  bool all_flux = (name.compare("mdot_vert") == 0);
  bool out_only = (name.compare("mdot_vert_out") == 0);

  DvceArray5D<Real> u0;
  if (pmbp->phydro != nullptr) {
    u0 = pmbp->phydro->u0;
  } else if (pmbp->pmhd != nullptr) {
    u0 = pmbp->pmhd->u0;
  }

  par_for("mdot_vert", DevExeSpace(), 0, (nmb-1), ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real z = CellCenterX(k-ks, nx3, size.d_view(m).x3min, size.d_view(m).x3max);
    Real sign_z = (z >= 0.0) ? 1.0 : -1.0;

    Real rho = u0(m, IDN, k, j, i);
    Real vz = u0(m, IM3, k, j, i) / rho;

    Real dV = size.d_view(m).dx1 * size.d_view(m).dx2 * size.d_view(m).dx3;

    // "Outward" flux = flux away from midplane
    Real v_out = vz * sign_z;
    Real mdot = rho * v_out * dV;

    if (all_flux) {
      dv(m, i_dv, k, j, i) = mdot;
    } else if (out_only) {
      dv(m, i_dv, k, j, i) = (v_out > 0.0) ? mdot : 0.0;
    } else {  // in_only
      dv(m, i_dv, k, j, i) = (v_out < 0.0) ? mdot : 0.0;
    }
  });
  i_dv += 1;
}
```

#### Vertical Energy Flux (edot_vert, edot_vert_out, edot_vert_in)

```cpp
// Vertical energy flux: edot_vert = F_E,z * sign(z) * dV
if (name.compare("edot_vert") == 0 ||
    name.compare("edot_vert_out") == 0 ||
    name.compare("edot_vert_in") == 0) {
  if (derived_var.extent(4) <= 1)
    Kokkos::realloc(derived_var, nmb, n_dv, n3, n2, n1);
  auto dv = derived_var;
  int nx3 = indcs.nx3;

  bool all_flux = (name.compare("edot_vert") == 0);
  bool out_only = (name.compare("edot_vert_out") == 0);

  bool is_mhd = (pmbp->pmhd != nullptr);

  DvceArray5D<Real> u0;
  DvceFaceFld4D<Real> b0;
  Real gamma;

  if (is_mhd) {
    u0 = pmbp->pmhd->u0;
    b0 = pmbp->pmhd->b0;
    gamma = pmbp->pmhd->peos->eos_data.gamma;
  } else {
    u0 = pmbp->phydro->u0;
    gamma = pmbp->phydro->peos->eos_data.gamma;
  }

  Real gm1 = gamma - 1.0;

  par_for("edot_vert", DevExeSpace(), 0, (nmb-1), ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real z = CellCenterX(k-ks, nx3, size.d_view(m).x3min, size.d_view(m).x3max);
    Real sign_z = (z >= 0.0) ? 1.0 : -1.0;

    Real rho = u0(m, IDN, k, j, i);
    Real vx = u0(m, IM1, k, j, i) / rho;
    Real vy = u0(m, IM2, k, j, i) / rho;
    Real vz = u0(m, IM3, k, j, i) / rho;
    Real v2 = vx*vx + vy*vy + vz*vz;

    Real dV = size.d_view(m).dx1 * size.d_view(m).dx2 * size.d_view(m).dx3;

    Real F_Ez;
    if (is_mhd) {
      Real Bx = 0.5*(b0.x1f(m,k,j,i) + b0.x1f(m,k,j,i+1));
      Real By = 0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i));
      Real Bz = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
      Real B2 = Bx*Bx + By*By + Bz*Bz;

      Real vdotB = vx*Bx + vy*By + vz*Bz;

      Real E = u0(m, IEN, k, j, i);
      Real P = gm1 * (E - 0.5*rho*v2 - 0.5*B2);

      // F_E,z = [γP/(γ-1) + ρv²/2 + B²] * v_z - B_z(v·B)
      F_Ez = (gamma/gm1*P + 0.5*rho*v2 + B2) * vz - Bz * vdotB;
    } else {
      Real E = u0(m, IEN, k, j, i);
      Real P = gm1 * (E - 0.5*rho*v2);

      F_Ez = (gamma/gm1*P + 0.5*rho*v2) * vz;
    }

    // "Outward" = away from midplane
    Real v_out = vz * sign_z;
    Real edot = F_Ez * sign_z * dV;

    if (all_flux) {
      dv(m, i_dv, k, j, i) = edot;
    } else if (out_only) {
      dv(m, i_dv, k, j, i) = (v_out > 0.0) ? edot : 0.0;
    } else {  // in_only
      dv(m, i_dv, k, j, i) = (v_out < 0.0) ? edot : 0.0;
    }
  });
  i_dv += 1;
}
```

### 3.6 Example Usage

**Computing radial mass flux profile:**

Input file:
```
<output5>
file_type   = pdf
variable_1  = coord_r
bin1_min    = 0.1
bin1_max    = 100.0
nbin1       = 100
logscale1   = true

variable_2  = mdot_sph
bin2_min    = -1e10
bin2_max    = 1e10
nbin2       = 1
logscale2   = false

mass_weighted = false
dt          = 0.1
```

Python analysis:
```python
result = read_pdf("pdf_mdot/turb", 0)
bin_edges = result['meta']['dims'][0]['bin_edges']
delta_r = bin_edges[1:] - bin_edges[:-1]

# Sum over the mdot dimension (which has only 1 bin + 2 overflow)
# Inner bin is index 1
mdot_sum = result['data'][1:-1, 1]  # exclude r overflow bins, take inner mdot bin

# Divide by radial bin width to get Ṁ(r)
mdot_profile = mdot_sum / delta_r

# Plot
r_centers = result['bin_centers'][0]
plt.loglog(r_centers, mdot_profile)
plt.xlabel('r')
plt.ylabel('Ṁ(r)')
```

### 3.7 Files to Modify (Part 3)

| File | Changes |
|------|---------|
| `src/outputs/outputs.hpp` | Add 12 flux variable names to `var_choice[]`, update `NOUTPUT_CHOICES` to 179 |
| `src/outputs/basetype_output.cpp` | Add registration for all 12 flux variables |
| `src/outputs/derived_variables.cpp` | Add computation kernels for all 12 flux variables |

### 3.8 Testing Plan (Part 3)

1. **Unit test**: Uniform density and radial velocity field, verify `mdot_sph` integrates to correct Ṁ
2. **Symmetry test**: Symmetric inflow/outflow, verify `mdot_sph_in` + `mdot_sph_out` = `mdot_sph`
3. **Hydro vs MHD**: Run same problem with B=0 in MHD, verify identical results to hydro
4. **Energy conservation**: Verify `edot_sph` integrates correctly for known solutions
5. **Vertical flux**: Test with vertical outflow, verify sign conventions work correctly
