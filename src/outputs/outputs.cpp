//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file outputs.cpp
//! \brief implements Outputs class constructor
//!
//! The number and types of outputs are all controlled by the number and values of
//! parameters specified in <output[n]> blocks in the input file.  Each output block must
//! be labelled by a unique integer "n".  Following the convention of the parser
//! implemented in the ParameterInput class, a second output block with the same integer
//! "n" of an earlier block will silently overwrite the values read by the first block.
//! Numbering of output blocks does not need to be consecutive, and blocks may appear
//! in any order in the input file. A new output type will be created for each and every
//! <output[n]> block in the input file.
//!
//! Required parameters that must be specified in an <output[n]> block are:
//!   - variable  = [list of currently implemented strings for specifing output variables
//!                  is defined at start of outputs.hpp file]
//!   - file_type = tab,vtk,hst,bin,rst
//!   - dt        = problem time between outputs
//!
//! EXAMPLE of an <output[n]> block for a TAB dump:
//!   <output3>
//!   file_type   = tab       # Tabular data dump
//!   variable    = prim      # variables to be output
//!   data_format = %12.5e    # Optional data format string
//!   dt          = 0.01      # time increment between outputs
//!   slice_x2    = 0.0       # slice at x2
//!   slice_x3    = 0.0       # slice at x3
//!
//! Each <output[n]> block will result in a new node being created in a linked list of
//! BaseTypeOutput stored in the Outputs class.  During a simulation, outputs are made
//! when the simulation time satisfies the criteria implemented in the Driver class.
//!
//! To implement a new output type, write a new BaseTypeOutput derived class and construct
//! an object of this class in the Outputs constructor at the location indicated by the
//! comment text: 'NEW_OUTPUT_TYPES'.
//========================================================================================

#include <cstdio>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>    // strcmp
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>   // std::string, to_string()

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "outputs.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {

[[noreturn]] void FatalPDFConfiguration(const OutputParameters &op,
                                        const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << "PDF output block '" << op.block_name << "' "
            << message << std::endl;
#if MPI_PARALLEL_ENABLED
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
  std::exit(EXIT_FAILURE);
}

std::string OutputVariable(ParameterInput *pin, const OutputParameters &op) {
  if (op.file_type == "pdf" &&
      pin->DoesParameterExist(op.block_name, "variable_1")) {
    return pin->GetString(op.block_name, "variable_1");
  }
  return pin->GetString(op.block_name, "variable");
}

void ParsePDFParameters(ParameterInput *pin, const InputBlock &block,
                        OutputParameters *op) {
  const auto fail = [&](const std::string &message) {
    FatalPDFConfiguration(*op, message);
  };
  const auto parse_scale = [&](const std::string &key) {
    std::string value = pin->GetString(op->block_name, key);
    if (value == "linear") return PDF_SCALE_LINEAR;
    if (value == "log") return PDF_SCALE_LOG;
    if (value == "symlog") return PDF_SCALE_SYMLOG;
    fail("has invalid " + key + "='" + value +
         "'; expected linear, log, or symlog");
    return PDF_SCALE_LINEAR;
  };
  const auto set_scale = [&](int dimension, const std::string &scale_key,
                             const std::string &log_key,
                             const std::string &linthresh_key,
                             bool legacy_log_default) {
    bool has_scale = pin->DoesParameterExist(op->block_name, scale_key);
    bool has_log = pin->DoesParameterExist(op->block_name, log_key);
    bool has_linthresh = pin->DoesParameterExist(op->block_name, linthresh_key);
    int scale = has_scale ? parse_scale(scale_key) :
        ((has_log ? pin->GetBoolean(op->block_name, log_key) :
                    legacy_log_default) ? PDF_SCALE_LOG : PDF_SCALE_LINEAR);
    if (has_scale && has_log) {
      int legacy_scale = pin->GetBoolean(op->block_name, log_key) ?
          PDF_SCALE_LOG : PDF_SCALE_LINEAR;
      if (scale != legacy_scale) {
        fail("has conflicting " + scale_key + " and " + log_key);
      }
    }
    if (scale == PDF_SCALE_SYMLOG) {
      if (!has_linthresh) {
        fail("requires " + linthresh_key + " when " + scale_key + "=symlog");
      }
      op->pdf_linthresh[dimension] =
          pin->GetReal(op->block_name, linthresh_key);
    } else {
      if (has_linthresh) {
        fail("cannot set " + linthresh_key + " unless " + scale_key + "=symlog");
      }
      op->pdf_linthresh[dimension] = 1.0;
    }
    op->pdf_scale[dimension] = scale;
  };

  bool has_mass_weighted =
      pin->DoesParameterExist(op->block_name, "mass_weighted");
  bool legacy_mass = has_mass_weighted ?
      pin->GetBoolean(op->block_name, "mass_weighted") : false;
  std::string translated_weight = legacy_mass ? "mass" : "volume";
  if (pin->DoesParameterExist(op->block_name, "weight")) {
    op->pdf_weight = pin->GetString(op->block_name, "weight");
    if (has_mass_weighted && op->pdf_weight != translated_weight) {
      fail("has inconsistent mass_weighted and weight settings");
    }
  } else {
    op->pdf_weight = translated_weight;
  }
  if (op->pdf_weight != "volume" && op->pdf_weight != "mass" &&
      op->pdf_weight != "variable") {
    fail("has invalid weight='" + op->pdf_weight +
         "'; expected volume, mass, or variable");
  }
  if (op->pdf_weight == "variable") {
    if (!pin->DoesParameterExist(op->block_name, "weight_variable")) {
      fail("requires weight_variable when weight=variable");
    }
    op->pdf_weight_variable =
        pin->GetString(op->block_name, "weight_variable");
  }
  op->mass_weighted = (op->pdf_weight == "mass");

  bool modern = pin->DoesParameterExist(op->block_name, "variable_1");
  bool requests_modern_storage = modern ||
      pin->DoesParameterExist(op->block_name, "weight") ||
      pin->DoesParameterExist(op->block_name, "scale") ||
      pin->DoesParameterExist(op->block_name, "scale1") ||
      pin->DoesParameterExist(op->block_name, "scale2") ||
      pin->DoesParameterExist(op->block_name, "linthresh") ||
      pin->DoesParameterExist(op->block_name, "linthresh1") ||
      pin->DoesParameterExist(op->block_name, "linthresh2");
  op->pdf_legacy_layout = !requests_modern_storage;

  if (modern) {
    bool gap = false;
    for (int dimension = 0; dimension < OutputParameters::PDF_MAX_DIM; ++dimension) {
      std::string suffix = std::to_string(dimension + 1);
      bool present =
          pin->DoesParameterExist(op->block_name, "variable_" + suffix);
      if (!present) {
        gap = true;
        continue;
      }
      if (gap) fail("has a gap in active variable_N dimensions");
      op->pdf_ndim = dimension + 1;
      op->pdf_variables[dimension] =
          pin->GetString(op->block_name, "variable_" + suffix);
      op->pdf_nbin[dimension] =
          pin->GetInteger(op->block_name, "nbin" + suffix);
      op->pdf_bin_min[dimension] =
          pin->GetReal(op->block_name, "bin" + suffix + "_min");
      op->pdf_bin_max[dimension] =
          pin->GetReal(op->block_name, "bin" + suffix + "_max");
      set_scale(dimension, "scale" + suffix, "logscale" + suffix,
                "linthresh" + suffix, false);
    }
    for (const auto &line : block.line) {
      const std::string prefix = "variable_";
      if (line.param_name.compare(0, prefix.size(), prefix) != 0) continue;
      const std::string suffix = line.param_name.substr(prefix.size());
      bool numeric_suffix = !suffix.empty();
      int dimension = 0;
      for (char ch : suffix) {
        if (ch < '0' || ch > '9') {
          numeric_suffix = false;
          break;
        }
        dimension = 10*dimension + static_cast<int>(ch - '0');
      }
      if (numeric_suffix && dimension > OutputParameters::PDF_MAX_DIM) {
        fail("requests more than four dimensions");
      }
    }
  } else {
    op->pdf_ndim = 1;
    op->pdf_variables[0] = op->variable;
    op->pdf_nbin[0] = pin->GetInteger(op->block_name, "nbin");
    op->pdf_bin_min[0] = pin->GetReal(op->block_name, "bin_min");
    op->pdf_bin_max[0] = pin->GetReal(op->block_name, "bin_max");
    std::string scale_key =
        pin->DoesParameterExist(op->block_name, "scale1") ? "scale1" : "scale";
    std::string log_key = pin->DoesParameterExist(op->block_name, "logscale1") ?
        "logscale1" : "logscale";
    std::string lin_key = pin->DoesParameterExist(op->block_name, "linthresh1") ?
        "linthresh1" : "linthresh";
    set_scale(0, scale_key, log_key, lin_key, true);
    if (pin->DoesParameterExist(op->block_name, "variable_2")) {
      op->pdf_ndim = 2;
      op->pdf_variables[1] = pin->GetString(op->block_name, "variable_2");
      op->pdf_nbin[1] = pin->GetInteger(op->block_name, "nbin2");
      op->pdf_bin_min[1] = pin->GetReal(op->block_name, "bin2_min");
      op->pdf_bin_max[1] = pin->GetReal(op->block_name, "bin2_max");
      set_scale(1, "scale2", "logscale2", "linthresh2", true);
    }
  }

  std::int64_t total_bins = 1;
  for (int dimension = 0; dimension < op->pdf_ndim; ++dimension) {
    const std::string &variable = op->pdf_variables[dimension];
    if (variable == "mhd_w" || variable == "mhd_u" ||
        variable == "hydro_w" || variable == "hydro_u") {
      fail("cannot output variable group '" + variable + "'");
    }
    if (op->pdf_nbin[dimension] <= 0) {
      fail("requires positive nbin for dimension " +
           std::to_string(dimension + 1));
    }
    if (!std::isfinite(op->pdf_bin_min[dimension]) ||
        !std::isfinite(op->pdf_bin_max[dimension])) {
      fail("requires finite bounds for dimension " +
           std::to_string(dimension + 1));
    }
    if (!(op->pdf_bin_min[dimension] < op->pdf_bin_max[dimension])) {
      fail("requires bin_min < bin_max for dimension " +
           std::to_string(dimension + 1));
    }
    if (op->pdf_scale[dimension] == PDF_SCALE_LOG &&
        (op->pdf_bin_min[dimension] <= 0.0 ||
         op->pdf_bin_max[dimension] <= 0.0)) {
      fail("requires positive bounds for logarithmic dimension " +
           std::to_string(dimension + 1));
    }
    if (op->pdf_scale[dimension] == PDF_SCALE_SYMLOG &&
        (!std::isfinite(op->pdf_linthresh[dimension]) ||
         op->pdf_linthresh[dimension] <= 0.0)) {
      fail("requires positive linthresh for symlog dimension " +
           std::to_string(dimension + 1));
    }
    Real transformed_min = PDFTransformValue(
        op->pdf_bin_min[dimension], op->pdf_scale[dimension],
        op->pdf_linthresh[dimension]);
    Real transformed_max = PDFTransformValue(
        op->pdf_bin_max[dimension], op->pdf_scale[dimension],
        op->pdf_linthresh[dimension]);
    Real step_size =
        (transformed_max - transformed_min)/op->pdf_nbin[dimension];
    if (!std::isfinite(transformed_min) || !std::isfinite(transformed_max) ||
        !std::isfinite(step_size) || !(step_size > 0.0)) {
      fail("requires finite transformed bounds and a positive finite bin step for "
           "dimension " + std::to_string(dimension + 1));
    }
    total_bins *= static_cast<std::int64_t>(op->pdf_nbin[dimension]) + 2;
    if (total_bins > std::numeric_limits<int>::max()) {
      fail("has too many total bins for a dense shared histogram");
    }
  }

  op->nbin = op->pdf_nbin[0];
  op->bin_min = op->pdf_bin_min[0];
  op->bin_max = op->pdf_bin_max[0];
  op->logscale = (op->pdf_scale[0] == PDF_SCALE_LOG);
  if (op->pdf_ndim > 1) {
    op->variable_2 = op->pdf_variables[1];
    op->nbin2 = op->pdf_nbin[1];
    op->bin2_min = op->pdf_bin_min[1];
    op->bin2_max = op->pdf_bin_max[1];
    op->logscale2 = (op->pdf_scale[1] == PDF_SCALE_LOG);
  }
}

}  // namespace

//----------------------------------------------------------------------------------------
// Outputs constructor

Outputs::Outputs(ParameterInput *pin, Mesh *pm) {
  // loop over input block names.  Find those that start with "output", read parameters,
  // and add to linked list of BaseTypeOutputs.

  int num_hst=0, num_rst=0, num_log=0; // count # of hst,rst,log outputs
  for (auto it = pin->block.begin(); it != pin->block.end(); ++it) {
    if (it->block_name.compare(0, 6, "output") == 0) {
      OutputParameters opar;  // define temporary OutputParameters struct

      // extract integer number of output block.  Save name and number
      std::string outn = it->block_name.substr(6); // 6 because counting starts at 0!
      opar.block_number = atoi(outn.c_str());
      opar.block_name.assign(it->block_name);

      // set time of last output, time or cycles between outputs
      // when last_time < 0, then outputs will always be made
      opar.last_time = pin->GetOrAddReal(opar.block_name,"last_time", -1.0);
      if (pin->DoesParameterExist(opar.block_name,"dcycle")) {
        opar.dcycle = pin->GetInteger(opar.block_name,"dcycle");
        opar.dt = 0.0;
      } else {
        opar.dt = pin->GetReal(opar.block_name,"dt");
        opar.dcycle = 0;
      }

      if (opar.dcycle == 0 && opar.dt <= 0.0) continue;  // only add output if dt>0

      // set file number, basename, and format
      opar.file_number = pin->GetOrAddInteger(opar.block_name,"file_number",0);
      opar.file_basename = pin->GetString("job","basename");
      opar.file_type = pin->GetString(opar.block_name,"file_type");

      // set output variable and optional file id (default is output variable name)
      // but only for those output types that use them
      if (opar.file_type.compare("hst") != 0 &&
          opar.file_type.compare("rst") != 0 &&
          opar.file_type.compare("log") != 0 &&
          opar.file_type.compare("trk") != 0 &&
          opar.file_type.compare("prtcl_thermo_history") != 0) {
        opar.variable = OutputVariable(pin, opar);
        opar.file_id = pin->GetOrAddString(opar.block_name,"id",opar.variable);
      }

      // read ghost cell option
      opar.include_gzs = pin->GetOrAddBoolean(opar.block_name, "ghost_zones", false);

      // read MeshBlock ID (if specified)
      opar.gid = pin->GetOrAddInteger(opar.block_name, "gid", -1);
      if (opar.gid >= 0 && pm->nmb_total == 1) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << "Cannot specify MeshBlock ID in output block '"
            << opar.block_name << "' when there is only one" << std::endl;
        exit(EXIT_FAILURE);
      }
      if (opar.gid > (pm->nmb_total-1)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << "MeshBlock gid=" << opar.gid << " in output block '"
            << opar.block_name << "' exceeds total number of MeshBlocks" << std::endl;
        exit(EXIT_FAILURE);
      }

      // read slicing options.  Check that slice is within mesh
      if (pin->DoesParameterExist(opar.block_name,"slice_x1")) {
        Real x1 = pin->GetReal(opar.block_name,"slice_x1");
        if (x1 >= pm->mesh_size.x1min && x1 < pm->mesh_size.x1max) {
          opar.slice_x1 = x1;
          opar.slice1 = true;
        } else {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Slice at x1=" << x1 << " in output block '"
              << opar.block_name << "' is out of range of Mesh" << std::endl;
          exit(EXIT_FAILURE);
        }
      } else {
        opar.slice1 = false;
      }

      if (pin->DoesParameterExist(opar.block_name,"slice_x2")) {
        Real x2 = pin->GetReal(opar.block_name,"slice_x2");
        if (x2 >= pm->mesh_size.x2min && x2 < pm->mesh_size.x2max) {
          opar.slice_x2 = x2;
          opar.slice2 = true;
        } else {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Slice at x2=" << x2 << " in output block '"
              << opar.block_name << "' is out of range of Mesh" << std::endl;
          exit(EXIT_FAILURE);
        }
      } else {
        opar.slice2 = false;
      }

      if (pin->DoesParameterExist(opar.block_name,"slice_x3")) {
        Real x3 = pin->GetReal(opar.block_name,"slice_x3");
        if (x3 >= pm->mesh_size.x3min && x3 < pm->mesh_size.x3max) {
          opar.slice_x3 = x3;
          opar.slice3 = true;
        } else {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Slice at x3=" << x3 << " in output block '"
              << opar.block_name << "' is out of range of Mesh" << std::endl;
          exit(EXIT_FAILURE);
        }
      } else {
        opar.slice3 = false;
      }

      // read ghost cell option
      opar.include_gzs = pin->GetOrAddBoolean(opar.block_name, "ghost_zones", false);

      // read MeshBlock ID (if specified)
      opar.gid = pin->GetOrAddInteger(opar.block_name, "gid", -1);
      if (opar.gid >= 0 && pm->nmb_total == 1) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << "Cannot specify MeshBlock ID in output block '"
            << opar.block_name << "' when there is only one" << std::endl;
        exit(EXIT_FAILURE);
      }
      if (opar.gid > (pm->nmb_total-1)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << "MeshBlock gid=" << opar.gid << " in output block '"
            << opar.block_name << "' exceeds total number of MeshBlocks" << std::endl;
        exit(EXIT_FAILURE);
      }

      // set output variable and optional file id (default is output variable name)
      if (opar.file_type.compare("hst") != 0 &&
          opar.file_type.compare("rst") != 0 &&
          opar.file_type.compare("log") != 0 &&
          opar.file_type.compare("prtcl_thermo_history") != 0) {
        opar.variable = OutputVariable(pin, opar);
        opar.file_id = pin->GetOrAddString(opar.block_name,"id",opar.variable);
      } else if (opar.file_type.compare("prtcl_thermo_history") == 0) {
        opar.variable = "prtcl_thermo_history";
        opar.file_id = pin->GetOrAddString(opar.block_name,"id",opar.variable);
      }

      // check that pdf variables are single variables
      // raise error if variable = mhd_w, mhd_u, hydro_w, hydro_u
      if (opar.file_type.compare("pdf") == 0) {
        if (opar.variable.compare("mhd_w") == 0 ||
            opar.variable.compare("mhd_u") == 0 ||
            opar.variable.compare("hydro_w") == 0 ||
            opar.variable.compare("hydro_u") == 0) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "PDF output block '" << opar.block_name
              << "' cannot output variable '" << opar.variable << "'."
              << " The variable must be a single variable not a variable group"
              << std::endl;
          exit(EXIT_FAILURE);
        }
      }

      // set optional boolean to output only user-defined history variables
      if (opar.file_type.compare("hst") == 0) {
        opar.user_hist_only =pin->GetOrAddBoolean(opar.block_name,"user_hist_only",false);
        if (opar.user_hist_only && !(pm->pgen->user_hist)) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "User-history file requested in output block '"
              << opar.block_name << "', but <problem>/user_hist is not true" << std::endl;
          exit(EXIT_FAILURE);
        }
      }

      // set optional data format string used in formatted writes
      opar.data_format = pin->GetOrAddString(opar.block_name, "data_format", "%12.5e");
      opar.data_format.insert(0, " "); // prepend with blank to separate columns

      // Construct new BaseTypeOutput according to file format
      // NEW_OUTPUT_TYPES: Add block to construct new types here
      BaseTypeOutput *pnode;
      if (opar.file_type.compare("tab") == 0) {
        pnode = new FormattedTableOutput(pin,pm,opar);
        pout_list.insert(pout_list.begin(),pnode);
      } else if (opar.file_type.compare("hst") == 0) {
        pnode = new HistoryOutput(pin,pm,opar);
        pout_list.insert(pout_list.begin(),pnode);
        num_hst++;
      } else if (opar.file_type.compare("log") == 0) {
        pnode = new EventLogOutput(pin,pm,opar);
        pout_list.insert(pout_list.begin(),pnode);
        num_log++;
      } else if (opar.file_type.compare("vtk") == 0) {
        pnode = new MeshVTKOutput(pin,pm,opar);
        pout_list.insert(pout_list.begin(),pnode);
      } else if (opar.file_type.compare("pvtk") == 0) {
        pnode = new ParticleVTKOutput(pin,pm,opar);
        pout_list.insert(pout_list.begin(),pnode);
      } else if (opar.file_type.compare("trk") == 0) {
        pnode = new TrackedParticleOutput(pin,pm,opar);
        pout_list.insert(pout_list.begin(),pnode);
      } else if (opar.file_type.compare("prtcl_thermo_history") == 0) {
        pnode = new ParticleThermoHistoryOutput(pin,pm,opar);
        pout_list.insert(pout_list.begin(),pnode);
      } else if (opar.file_type.compare("cbin") == 0) {
        opar.single_file_per_rank = pin->GetOrAddBoolean(opar.block_name,
          "single_file_per_rank", false);
        opar.coarsen_factor = pin->GetInteger(opar.block_name,"coarsen_factor");
        opar.compute_moments = pin->GetOrAddBoolean(opar.block_name,
          "compute_moments", false);
        pnode = new CoarsenedBinaryOutput(pin,pm,opar);
        pout_list.insert(pout_list.begin(),pnode);
      } else if (opar.file_type.compare("pdf") == 0) {
        ParsePDFParameters(pin, *it, &opar);
        pnode = new PDFOutput(pin,pm,opar);
        pout_list.insert(pout_list.begin(),pnode);
      } else if (opar.file_type.compare("bin") == 0) {
        opar.single_file_per_rank = pin->GetOrAddBoolean(opar.block_name,
          "single_file_per_rank", false);
        opar.data_precision = pin->GetOrAddString(opar.block_name,
          "data_precision", "float32");
        if (opar.data_precision.compare("float32") != 0 &&
            opar.data_precision.compare("real") != 0) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Invalid data_precision = '" << opar.data_precision
              << "' in binary output block '" << opar.block_name
              << "'. Expected float32 or real." << std::endl;
          exit(EXIT_FAILURE);
        }
        pnode = new MeshBinaryOutput(pin,pm,opar);
        pout_list.insert(pout_list.begin(),pnode);
      } else if (opar.file_type.compare("cart") == 0) {
        pnode = new CartesianGridOutput(pin,pm,opar);
        pout_list.insert(pout_list.begin(),pnode);
      } else if (opar.file_type.compare("sph") == 0) {
        pnode = new SphericalSurfaceOutput(pin,pm,opar);
        pout_list.insert(pout_list.begin(),pnode);
      } else if (opar.file_type.compare("rst") == 0) {
      // Add restarts to the tail end of BaseTypeOutput list, so file counters for other
      // output types are up-to-date in restart file
        opar.single_file_per_rank = pin->GetOrAddBoolean(opar.block_name,
          "single_file_per_rank", false);
        pnode = new RestartOutput(pin,pm,opar);
        pout_list.push_back(pnode);
        num_rst++;
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << "Unrecognized file format = '" << opar.file_type
            << "' in output block '" << opar.block_name << "'" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  }

  // check there were no more than one history, event log, or restart files requested
  if (num_hst > 1 || num_rst > 1 || num_log > 1) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "More than one history, event log, or restart output block found in "
              << "input file" << std::endl;
    exit(EXIT_FAILURE);
  }
}

//----------------------------------------------------------------------------------------
// destructor

Outputs::~Outputs() {
  // Must manually delete memory assigned to each OutputType object stored in pout_list
  for (BaseTypeOutput* pnode : pout_list) {
    delete pnode;
  }
  pout_list.clear();
}
