//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cartgrid.cpp
//! \brief writes data on a Cartesian sub-grid in binary format

#include <fstream>
#include <string>
#include <sstream>

#include "athena.hpp"
#include "globals.hpp"
#include "mpi_utils.hpp"
#include "mesh/mesh.hpp"
#include "output_file_utils.hpp"
#include "outputs.hpp"
#include "parameter_input.hpp"
#include "utils/cart_grid.hpp"

namespace {

[[noreturn]] void FatalCartesianGridError(const std::string &message) {
  mpi_utils::AbortWorld(std::string("### FATAL ERROR in ") + __FILE__ +
                        " at line " + std::to_string(__LINE__) + "\n" + message);
}

}  // namespace

CartesianGridOutput::CartesianGridOutput(ParameterInput *pin, Mesh *pm,
                                         OutputParameters op)
    : BaseTypeOutput(pin, pm, op) {
  Real center[3], extent[3];
  int numpoints[3];
  bool is_cheb;

  if (out_params.contains_derived) {
    FatalCartesianGridError(
        "Cartesian-grid derived-field interpolation is not supported until "
        "ghost-zone-safe sampling is implemented.");
  }
  output_file_utils::EnsureDirectory("cart", 0755, "Cartesian grid output",
                                     FatalCartesianGridError);

  center[0] = pin->GetOrAddReal(op.block_name, "center_x", 0.0);
  center[1] = pin->GetOrAddReal(op.block_name, "center_y", 0.0);
  center[2] = pin->GetOrAddReal(op.block_name, "center_z", 0.0);
  extent[0] = pin->GetOrAddReal(op.block_name, "extent_x", 2.0);
  extent[1] = pin->GetOrAddReal(op.block_name, "extent_y", 2.0);
  extent[2] = pin->GetOrAddReal(op.block_name, "extent_z", 2.0);
  numpoints[0] = pin->GetOrAddInteger(op.block_name, "numpoints_x", 32);
  numpoints[1] = pin->GetOrAddInteger(op.block_name, "numpoints_y", 32);
  numpoints[2] = pin->GetOrAddInteger(op.block_name, "numpoints_z", 32);
  is_cheb = pin->GetOrAddBoolean(op.block_name, "chebyshev", false);

  pcart = new CartesianGrid(pm->pmb_pack, center, extent, numpoints, is_cheb);

  for (int d = 0; d < 3; ++d) {
    md.center[d] = center[d];
    md.extent[d] = extent[d];
    md.numpoints[d] = numpoints[d];
  }
  md.is_cheb = is_cheb;
}

CartesianGridOutput::~CartesianGridOutput() { delete pcart; }

void CartesianGridOutput::LoadOutputData(Mesh *pm) {
  // If AMR is enabled we need to reset the CartesianGrid
  if (pm->adaptive) {
    pcart->SetInterpolationIndices();
    pcart->SetInterpolationWeights();
  }

  int nout_vars = outvars.size();
  Kokkos::realloc(outarray, nout_vars, 1, md.numpoints[0], md.numpoints[1],
                  md.numpoints[2]);

  // Calculate derived variables, if required
  if (out_params.contains_derived) {
    ComputeDerivedVariable(out_params.variable, pm);
  }

  for (int n = 0; n < nout_vars; ++n) {
    pcart->InterpolateToGrid(outvars[n].data_index, *(outvars[n].data_ptr));
    auto v_slice =
        Kokkos::subview(outarray, n, 0, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
    Kokkos::deep_copy(v_slice, pcart->interp_vals.h_view);
  }
#if MPI_PARALLEL_ENABLED
  // Note that InterpolateToGrid will set zero values for points not owned by
  // current rank
  int count = nout_vars * md.numpoints[0] * md.numpoints[1] * md.numpoints[2];
  if (0 == global_variable::my_rank) {
    mpi_utils::CheckMpi(
        MPI_Reduce(MPI_IN_PLACE, outarray.data(), count, MPI_ATHENA_REAL, MPI_SUM,
                   0, MPI_COMM_WORLD),
        "MPI_Reduce for Cartesian grid output values");
  } else {
    mpi_utils::CheckMpi(
        MPI_Reduce(outarray.data(), outarray.data(), count, MPI_ATHENA_REAL,
                   MPI_SUM, 0, MPI_COMM_WORLD),
        "MPI_Reduce for Cartesian grid output values");
  }
#endif
}

void CartesianGridOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
#if MPI_PARALLEL_ENABLED
  if (0 == global_variable::my_rank) {
#endif

    // Assemble filename
    std::string sequence = output_file_utils::FormatSequence(
        out_params.file_number, "Cartesian grid output", FatalCartesianGridError);
    std::string fname = "cart/" + out_params.file_basename + "." + out_params.file_id +
                        "." + sequence + ".bin";

    // Open file
    std::ofstream ofile(fname, std::ios::binary);
    if (!ofile.is_open()) {
      FatalCartesianGridError("Cannot open Cartesian grid output '" + fname + "'.");
    }

    // Write metadata
    md.cycle = pm->ncycle;
    md.time = pm->time;
    md.noutvars = outvars.size();
    ofile.write(reinterpret_cast<char *>(&md), sizeof(MetaData));

    // Write list of variables
    {
      std::stringstream msg;
      for (int n = 0; n < md.noutvars - 1; ++n) {
        msg << outvars[n].label << " ";
      }
      msg << outvars[md.noutvars - 1].label;
      std::string smsg = msg.str();
      int len = smsg.size();
      ofile.write(reinterpret_cast<char *>(&len), sizeof(int));
      ofile.write(smsg.c_str(), len);
    }

    // Write actual data
    for (int n = 0; n < md.noutvars; ++n) {
      for (int k = 0; k < md.numpoints[2]; ++k) {
        for (int j = 0; j < md.numpoints[1]; ++j) {
          for (int i = 0; i < md.numpoints[0]; ++i) {
            // Note that we are accessing the array with the convention of
            // CartesianGrid which is opposite from the one used in the rest of
            // the code, but we write the output as k, j, i
            float var = outarray(n, 0, i, j, k);
            ofile.write(reinterpret_cast<char *>(&var), sizeof(float));
          }
        }
      }
    }
    if (!ofile.good()) {
      ofile.close();
      FatalCartesianGridError("Could not write Cartesian grid output '" + fname +
                              "' completely.");
    }
    ofile.close();
    if (ofile.fail()) {
      FatalCartesianGridError("Could not close Cartesian grid output '" + fname + "'.");
    }

#if MPI_PARALLEL_ENABLED
  }
#endif

  // increment counters
  out_params.file_number = output_file_utils::AdvanceFileNumber(
      out_params.file_number, "Cartesian grid output", FatalCartesianGridError);
  if (out_params.last_time < 0.0) {
    out_params.last_time = pm->time;
  } else {
    out_params.last_time += out_params.dt;
  }
  pin->SetInteger(out_params.block_name, "file_number", out_params.file_number);
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
}
