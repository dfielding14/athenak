//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file outputs.cpp
//! \brief implements Outputs class constructor
//!
//! Output streams are controlled by <output...> blocks in the input file. Blocks may
//! appear in any order, and their suffixes do not need to be consecutive.
//!
//! Required parameters that must be specified in an <output[n]> block are:
//!   - file_type = tab,hst,log,vtk,pvtk,trk,cbin,pdf,bin,cart,sph,sphslice,rst
//!   - dt or dcycle
//! Most streams also require variable; hst, log, trk, and rst do not.
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
//! Each active block creates one BaseTypeOutput stored in the Outputs vector. During a
//! simulation, outputs are made
//! when the simulation time satisfies the criteria implemented in the Driver class.
//!
//! To implement a new output type, write a new BaseTypeOutput derived class and construct
//! an object of this class in the Outputs constructor at the location indicated by the
//! comment text: 'NEW_OUTPUT_TYPES'.
//========================================================================================

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>    // strcmp
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>   // std::string, to_string()

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "output_file_utils.hpp"
#include "outputs.hpp"

namespace {

[[noreturn]] void FatalOutputsError(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << message << std::endl;
#if MPI_PARALLEL_ENABLED
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
  std::exit(EXIT_FAILURE);
}

FileShardMode ParseShardMode(ParameterInput *pin, const std::string &block_name,
                             bool initialize_node_communicator=true) {
  bool per_rank = pin->GetOrAddBoolean(block_name, "single_file_per_rank", false);
  bool per_node = pin->GetOrAddBoolean(block_name, "single_file_per_node", false);
  if (per_rank && per_node) {
    FatalOutputsError("Output block '" + block_name +
                      "' cannot set both single_file_per_rank=true and "
                      "single_file_per_node=true.");
  }
  if (per_node) {
    if (initialize_node_communicator) {
      global_variable::InitializeNodeCommunicator();
    }
    return FileShardMode::node;
  }
  return per_rank ? FileShardMode::rank : FileShardMode::shared;
}

std::string PartitionTemplate(FileShardMode mode) {
  if (mode == FileShardMode::rank) return "{RANK}/";
  if (mode == FileShardMode::node) return "{NODE}/";
  return "";
}

int ParseFileNumber(ParameterInput *pin, const std::string &block_name) {
  if (!pin->DoesParameterExist(block_name, "file_number")) {
    return pin->GetOrAddInteger(block_name, "file_number", 0);
  }
  std::string text = pin->GetString(block_name, "file_number");
  if (text.empty() ||
      !std::all_of(text.begin(), text.end(), [](unsigned char ch) {
        return std::isdigit(ch) != 0;
      })) {
    FatalOutputsError("Output block '" + block_name +
                      "' requires file_number to be a non-negative integer.");
  }
  std::uint64_t value = 0;
  try {
    value = std::stoull(text);
  } catch (const std::exception &) {
    FatalOutputsError("Output block '" + block_name +
                      "' has an unrepresentable file_number.");
  }
  if (value > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
    FatalOutputsError("Output block '" + block_name +
                      "' file_number is outside the publishable range.");
  }
  return static_cast<int>(value);
}

bool OutputBlockIsActive(ParameterInput *pin, const std::string &block_name) {
  if (pin->DoesParameterExist(block_name, "dcycle")) {
    return pin->GetInteger(block_name, "dcycle") != 0;
  }
  Real dt = pin->GetReal(block_name, "dt");
  if (!std::isfinite(dt)) {
    FatalOutputsError("Output block '" + block_name + "' requires finite dt.");
  }
  return dt > 0.0;
}

std::string OutputVariable(ParameterInput *pin, const std::string &block_name,
                           const std::string &file_type) {
  if (file_type == "pdf" && pin->DoesParameterExist(block_name, "variable_1")) {
    return pin->GetString(block_name, "variable_1");
  }
  return pin->GetString(block_name, "variable");
}

bool UsesVariableAndId(const std::string &file_type) {
  return file_type != "hst" && file_type != "rst" && file_type != "log" &&
      file_type != "trk" && file_type != "prtcl_thermo_history";
}

void ParseVariableAndId(ParameterInput *pin, OutputParameters *op) {
  if (!UsesVariableAndId(op->file_type)) return;
  op->variable = OutputVariable(pin, op->block_name, op->file_type);
  op->file_id = pin->GetOrAddString(op->block_name, "id", op->variable);
}

void ParseMeshBlockSelection(ParameterInput *pin, Mesh *pm, OutputParameters *op) {
  op->gid = pin->GetOrAddInteger(op->block_name, "gid", -1);
  if (op->gid >= 0 && pm->nmb_total == 1) {
    FatalOutputsError("Cannot specify MeshBlock ID in output block '" +
                      op->block_name + "' when there is only one.");
  }
  if (op->gid > (pm->nmb_total - 1)) {
    FatalOutputsError("MeshBlock gid=" + std::to_string(op->gid) +
                      " in output block '" + op->block_name +
                      "' exceeds total number of MeshBlocks.");
  }
}

void ParseSlices(ParameterInput *pin, Mesh *pm, OutputParameters *op) {
  const auto parse_slice = [&](const char *key, Real minimum, Real maximum,
                               Real *coordinate, bool *enabled) {
    if (!pin->DoesParameterExist(op->block_name, key)) return;
    Real value = pin->GetReal(op->block_name, key);
    if (!(value >= minimum && value < maximum)) {
      FatalOutputsError(std::string("Slice at ") + key + "=" +
                        std::to_string(value) + " in output block '" +
                        op->block_name + "' is out of range of Mesh.");
    }
    *coordinate = value;
    *enabled = true;
  };
  parse_slice("slice_x1", pm->mesh_size.x1min, pm->mesh_size.x1max,
              &op->slice_x1, &op->slice1);
  parse_slice("slice_x2", pm->mesh_size.x2min, pm->mesh_size.x2max,
              &op->slice_x2, &op->slice2);
  parse_slice("slice_x3", pm->mesh_size.x3min, pm->mesh_size.x3max,
              &op->slice_x3, &op->slice3);
}

void ParsePdfParameters(ParameterInput *pin, const InputBlock &block,
                        OutputParameters *op) {
  const auto fail_pdf = [&](const std::string &message) {
    FatalOutputsError("PDF output block '" + op->block_name + "' " + message);
  };
  const auto parse_scale = [&](const std::string &key) {
    std::string value = pin->GetString(op->block_name, key);
    if (value == "linear") return PDF_SCALE_LINEAR;
    if (value == "log") return PDF_SCALE_LOG;
    if (value == "symlog") return PDF_SCALE_SYMLOG;
    fail_pdf("has invalid " + key + "='" + value +
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
        fail_pdf("has conflicting " + scale_key + " and " + log_key);
      }
    }
    if (scale == PDF_SCALE_SYMLOG) {
      if (!has_linthresh) {
        fail_pdf("requires " + linthresh_key + " when " + scale_key + "=symlog");
      }
      op->pdf_linthresh[dimension] = pin->GetReal(op->block_name, linthresh_key);
    } else {
      if (has_linthresh) {
        fail_pdf("cannot set " + linthresh_key + " unless " + scale_key + "=symlog");
      }
      op->pdf_linthresh[dimension] = 1.0;
    }
    op->pdf_scale[dimension] = scale;
  };

  if (op->variable == "mhd_w" || op->variable == "mhd_u" ||
      op->variable == "hydro_w" || op->variable == "hydro_u") {
    fail_pdf("cannot output variable '" + op->variable +
             "'. The variable must be a single variable not a variable group");
  }

  bool has_mass_weighted = pin->DoesParameterExist(op->block_name, "mass_weighted");
  bool legacy_mass = has_mass_weighted ?
      pin->GetBoolean(op->block_name, "mass_weighted") : false;
  std::string translated_weight = legacy_mass ? "mass" : "volume";
  if (pin->DoesParameterExist(op->block_name, "weight")) {
    op->pdf_weight = pin->GetString(op->block_name, "weight");
    if (has_mass_weighted && op->pdf_weight != translated_weight) {
      fail_pdf("has inconsistent mass_weighted and weight settings");
    }
  } else {
    op->pdf_weight = translated_weight;
  }
  if (op->pdf_weight != "volume" && op->pdf_weight != "mass" &&
      op->pdf_weight != "variable") {
    fail_pdf("has invalid weight='" + op->pdf_weight +
             "'; expected volume, mass, or variable");
  }
  if (op->pdf_weight == "variable") {
    if (!pin->DoesParameterExist(op->block_name, "weight_variable")) {
      fail_pdf("requires weight_variable when weight=variable");
    }
    op->pdf_weight_variable = pin->GetString(op->block_name, "weight_variable");
  }
  op->mass_weighted = (op->pdf_weight == "mass");

  bool modern = pin->DoesParameterExist(op->block_name, "variable_1");
  op->shard_mode = ParseShardMode(pin, op->block_name);
  bool requests_modern_storage = modern || IsSharded(op->shard_mode) ||
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
    for (int dimension = 0; dimension < op->PDF_MAX_DIM; ++dimension) {
      std::string suffix = std::to_string(dimension + 1);
      bool present =
          pin->DoesParameterExist(op->block_name, "variable_" + suffix);
      if (!present) {
        gap = true;
        continue;
      }
      if (gap) {
        fail_pdf("has a gap in active variable_N dimensions");
      }
      op->pdf_ndim = dimension + 1;
      op->pdf_variables[dimension] =
          pin->GetString(op->block_name, "variable_" + suffix);
      op->pdf_nbin[dimension] = pin->GetInteger(op->block_name, "nbin" + suffix);
      op->pdf_bin_min[dimension] =
          pin->GetReal(op->block_name, "bin" + suffix + "_min");
      op->pdf_bin_max[dimension] =
          pin->GetReal(op->block_name, "bin" + suffix + "_max");
      set_scale(dimension, "scale" + suffix, "logscale" + suffix,
                "linthresh" + suffix, false);
    }
    for (const auto &line : block.line) {
      const std::string prefix = "variable_";
      if (line.param_name.compare(0, prefix.size(), prefix) != 0) {
        continue;
      }
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
        fail_pdf("requests more than four dimensions");
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
    std::string log_key =
        pin->DoesParameterExist(op->block_name, "logscale1") ?
        "logscale1" : "logscale";
    std::string lin_key =
        pin->DoesParameterExist(op->block_name, "linthresh1") ?
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
    if (op->pdf_nbin[dimension] <= 0) {
      fail_pdf("requires positive nbin for dimension " +
               std::to_string(dimension + 1));
    }
    if (!std::isfinite(op->pdf_bin_min[dimension]) ||
        !std::isfinite(op->pdf_bin_max[dimension])) {
      fail_pdf("requires finite bounds for dimension " +
               std::to_string(dimension + 1));
    }
    if (!(op->pdf_bin_min[dimension] < op->pdf_bin_max[dimension])) {
      fail_pdf("requires bin_min < bin_max for dimension " +
               std::to_string(dimension + 1));
    }
    if (op->pdf_scale[dimension] == PDF_SCALE_LOG &&
        (op->pdf_bin_min[dimension] <= 0.0 || op->pdf_bin_max[dimension] <= 0.0)) {
      fail_pdf("requires positive bounds for logarithmic dimension " +
               std::to_string(dimension + 1));
    }
    if (op->pdf_scale[dimension] == PDF_SCALE_SYMLOG &&
        (!std::isfinite(op->pdf_linthresh[dimension]) ||
         op->pdf_linthresh[dimension] <= 0.0)) {
      fail_pdf("requires positive linthresh for symlog dimension " +
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
      fail_pdf("requires finite transformed bounds and a positive finite bin "
               "step for dimension " + std::to_string(dimension + 1));
    }
    total_bins *= static_cast<std::int64_t>(op->pdf_nbin[dimension]) + 2;
    if (total_bins > std::numeric_limits<int>::max()) {
      fail_pdf("has too many total bins for a dense shared histogram");
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

void ReserveOutputNamespaces(ParameterInput *pin) {
  const std::string basename = pin->GetString("job", "basename");
  std::map<std::string, std::string> owner_by_target;
  const auto reserve = [&](const std::string &target, const std::string &owner) {
    const std::string normalized_target =
        output_file_utils::LexicallyNormalTarget(target);
    auto inserted = owner_by_target.emplace(normalized_target, owner);
    if (!inserted.second && inserted.first->second != owner) {
      FatalOutputsError("Output blocks '" + inserted.first->second + "' and '" + owner +
                        "' resolve to the same public target family '" +
                        normalized_target + "'.");
    }
  };

  for (const auto &block : pin->block) {
    const std::string &name = block.block_name;
    if (name.compare(0, 6, "output") != 0 || !OutputBlockIsActive(pin, name)) {
      continue;
    }
    const std::string type = pin->GetString(name, "file_type");
    ParseFileNumber(pin, name);
    FileShardMode shard_mode = FileShardMode::shared;
    if (type == "bin" || type == "cbin" || type == "pdf" ||
        type == "sphslice" || type == "rst") {
      shard_mode = ParseShardMode(pin, name, false);
    }
    const std::string partition = PartitionTemplate(shard_mode);

    if (type == "hst") {
      reserve("history:" + basename, name);
    } else if (type == "log") {
      reserve(basename + ".log", name);
    } else if (type == "trk") {
      reserve("trk/" + basename + ".trk", name);
    } else if (type == "prtcl_thermo_history") {
      const std::string id = pin->DoesParameterExist(name, "id")
          ? pin->GetString(name, "id") : "prtcl_thermo_history";
      output_file_utils::ValidatePathComponent(
          id, "Output block '" + name + "' id", FatalOutputsError);
      reserve("prtcl_thermo_history/" + basename + "." + id + ".thp", name);
    } else if (type == "rst") {
      reserve("restart:", name);
      if (shard_mode == FileShardMode::node) {
        reserve("rst/" + basename + ".{SEQ}.rst", name);
        reserve("rst/{NODE}/" + basename + ".{SEQ}.g{GEN}.payload.rst", name);
      } else {
        reserve("rst/" + partition + basename + ".{SEQ}.rst", name);
      }
    } else {
      const std::string variable = OutputVariable(pin, name, type);
      const std::string id = pin->DoesParameterExist(name, "id")
          ? pin->GetString(name, "id") : variable;
      output_file_utils::ValidatePathComponent(
          id, "Output block '" + name + "' id", FatalOutputsError);
      if (type == "tab") {
        reserve("tab/" + basename + "." + id + ".{SEQ}.tab", name);
      } else if (type == "vtk") {
        const int configured_gid =
            pin->DoesParameterExist(name, "gid") ? pin->GetInteger(name, "gid") : -1;
        std::string gid =
            configured_gid >= 0 ? "." + std::to_string(configured_gid) : "";
        reserve("vtk/" + basename + "." + id + gid + ".{SEQ}.vtk", name);
      } else if (type == "pvtk") {
        const int configured_gid =
            pin->DoesParameterExist(name, "gid") ? pin->GetInteger(name, "gid") : -1;
        std::string gid =
            configured_gid >= 0 ? "." + std::to_string(configured_gid) : "";
        reserve("pvtk/" + basename + "." + id + gid + ".{SEQ}.part.vtk", name);
      } else if (type == "bin") {
        reserve("bin/" + partition + basename + "." + id + ".{SEQ}.bin", name);
      } else if (type == "cbin") {
        std::string factor = std::to_string(pin->GetInteger(name, "coarsen_factor"));
        reserve("cbin_" + id + "_" + factor + "/" + partition + basename + "." + id +
                ".{SEQ}.cbin", name);
      } else if (type == "cart") {
        reserve("cart/" + basename + "." + id + ".{SEQ}.bin", name);
      } else if (type == "sph") {
        std::ostringstream radius;
        radius << std::fixed << std::setprecision(2) << pin->GetReal(name, "radius");
        reserve("sph/" + basename + ".r=" + radius.str() + "." + id + ".{SEQ}.vtk",
                name);
      } else if (type == "sphslice") {
        reserve("bin/" + partition + basename + "." + id + "." +
                output_file_utils::FormatSphericalSliceRadius(
                    pin->GetReal(name, "slice_r")) +
                ".{SEQ}.sph.bin", name);
      } else if (type == "pdf") {
        std::string directory = "pdf_" + id;
        for (int dimension = 2; dimension <= OutputParameters::PDF_MAX_DIM; ++dimension) {
          std::string key = "variable_" + std::to_string(dimension);
          if (!pin->DoesParameterExist(name, key)) break;
          std::string component = pin->GetString(name, key);
          output_file_utils::ValidatePathComponent(
              component, "PDF output block '" + name + "' " + key, FatalOutputsError);
          directory += "_" + component;
        }
        bool modern = pin->DoesParameterExist(name, "variable_1") ||
            IsSharded(shard_mode) || pin->DoesParameterExist(name, "weight") ||
            pin->DoesParameterExist(name, "scale") ||
            pin->DoesParameterExist(name, "scale1") ||
            pin->DoesParameterExist(name, "scale2") ||
            pin->DoesParameterExist(name, "linthresh") ||
            pin->DoesParameterExist(name, "linthresh1") ||
            pin->DoesParameterExist(name, "linthresh2");
        reserve(directory + "/" + partition + basename +
                (modern ? ".header.pdf" : ".bins.pdf"), name);
        reserve(directory + "/" + partition + basename + ".{SEQ}.pdf", name);
      }
    }
  }
}

void ValidateCoarsenFactor(Mesh *pm, const std::string &block_name, int factor) {
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  if (pm->multilevel) {
    FatalOutputsError(
        "Coarsened-binary output supports uniform meshes only; static refinement "
        "and AMR are not supported.");
  }
  if (indcs.nx2 <= 1 || indcs.nx3 <= 1) {
    FatalOutputsError(
        "Coarsened-binary output supports three-dimensional meshes only.");
  }
  int shortest = std::min({indcs.nx1, indcs.nx2, indcs.nx3});
  if (factor < 2 || (factor & (factor - 1)) != 0 || factor > shortest) {
    FatalOutputsError("Coarsened-binary output block '" + block_name +
                      "' requires coarsen_factor to be a power of two between 2 and "
                      "the shortest MeshBlock dimension (" +
                      std::to_string(shortest) + ").");
  }
}

}  // namespace

//----------------------------------------------------------------------------------------
// Outputs constructor

Outputs::Outputs(ParameterInput *pin, Mesh *pm) {
  // loop over input block names.  Find those that start with "output", read parameters,
  // and append the corresponding BaseTypeOutput owner to the output vector.

  ReserveOutputNamespaces(pin);
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
      opar.file_number = ParseFileNumber(pin, opar.block_name);
      opar.file_basename = pin->GetString("job","basename");
      opar.file_type = pin->GetString(opar.block_name,"file_type");

      ParseVariableAndId(pin, &opar);
      if (opar.file_type == "prtcl_thermo_history") {
        opar.variable = "prtcl_thermo_history";
        opar.file_id = pin->GetOrAddString(opar.block_name, "id", opar.variable);
      }
      opar.include_gzs = pin->GetOrAddBoolean(opar.block_name, "ghost_zones", false);
      ParseMeshBlockSelection(pin, pm, &opar);
      ParseSlices(pin, pm, &opar);

      // set optional boolean to output only user-defined history variables
      if (opar.file_type.compare("hst") == 0) {
        opar.user_hist_only =pin->GetOrAddBoolean(opar.block_name,"user_hist_only",false);
        if (opar.user_hist_only && !(pm->pgen->user_hist)) {
          FatalOutputsError("User-history file requested in output block '" +
                            opar.block_name +
                            "', but <problem>/user_hist is not true.");
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
        opar.shard_mode = ParseShardMode(pin, opar.block_name);
        opar.coarsen_factor = pin->GetInteger(opar.block_name,"coarsen_factor");
        ValidateCoarsenFactor(pm, opar.block_name, opar.coarsen_factor);
        opar.compute_moments = pin->GetOrAddBoolean(opar.block_name,
          "compute_moments", false);
        pnode = new CoarsenedBinaryOutput(pin,pm,opar);
        pout_list.insert(pout_list.begin(),pnode);
      } else if (opar.file_type.compare("pdf") == 0) {
        ParsePdfParameters(pin, *it, &opar);
        pnode = new PDFOutput(pin,pm,opar);
        pout_list.insert(pout_list.begin(),pnode);
      } else if (opar.file_type.compare("bin") == 0) {
        opar.shard_mode = ParseShardMode(pin, opar.block_name);
        opar.data_precision = pin->GetOrAddString(opar.block_name,
          "data_precision", "float32");
        if (opar.data_precision.compare("float32") != 0 &&
            opar.data_precision.compare("real") != 0) {
          FatalOutputsError("Invalid data_precision = '" + opar.data_precision +
                            "' in binary output block '" + opar.block_name +
                            "'. Expected float32 or real.");
        }
        pnode = new MeshBinaryOutput(pin,pm,opar);
        pout_list.insert(pout_list.begin(),pnode);
      } else if (opar.file_type.compare("cart") == 0) {
        pnode = new CartesianGridOutput(pin,pm,opar);
        pout_list.insert(pout_list.begin(),pnode);
      } else if (opar.file_type.compare("sph") == 0) {
        pnode = new SphericalSurfaceOutput(pin,pm,opar);
        pout_list.insert(pout_list.begin(),pnode);
      } else if (opar.file_type.compare("sphslice") == 0) {
        opar.shard_mode = ParseShardMode(pin, opar.block_name);
        pnode = new SphericalSliceOutput(pin,pm,opar);
        pout_list.insert(pout_list.begin(),pnode);
      } else if (opar.file_type.compare("rst") == 0) {
      // Add restarts to the tail end of BaseTypeOutput list, so file counters for other
      // output types are up-to-date in restart file
        opar.shard_mode = ParseShardMode(pin, opar.block_name);
        pnode = new RestartOutput(pin,pm,opar);
        pout_list.push_back(pnode);
        num_rst++;
      } else {
        FatalOutputsError("Unrecognized file format = '" + opar.file_type +
                          "' in output block '" + opar.block_name + "'.");
      }
    }
  }

  // check there were no more than one history, event log, or restart files requested
  if (num_hst > 1 || num_rst > 1 || num_log > 1) {
    FatalOutputsError(
        "More than one history, event log, or restart output block found in input file.");
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
