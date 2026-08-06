//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file formatted_table.cpp
//  \brief writes output data as a formatted (ASCI) table.  Since outputing data in this
//  format is very slow and creates large files, it cannot be used for anything other than
//  1D slices.  Code will issue error if this format is selected for 2D or 3D outputs.
//  Output is written to a single file even with multiple MeshBlocks and MPI ranks.

#include <cstdio>      // fwrite(), fclose(), fopen(), fnprintf(), snprintf()
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "mpi_utils.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "output_file_utils.hpp"
#include "outputs.hpp"

namespace {

[[noreturn]] void FatalFormattedTableError(const std::string &message) {
  mpi_utils::AbortWorld(std::string("### FATAL ERROR in ") + __FILE__ +
                        " at line " + std::to_string(__LINE__) + "\n" +
                        message);
}

void CheckedFormattedTablePrint(std::FILE *output, const std::string &filename,
                                const char *text) {
  if (std::fputs(text, output) == EOF) {
    FatalFormattedTableError("Could not write formatted table output '" +
                             filename + "'.");
  }
}

template <typename Arg, typename... Args>
void CheckedFormattedTablePrint(std::FILE *output, const std::string &filename,
                                const char *format, Arg arg, Args... args) {
  if (std::fprintf(output, format, arg, args...) < 0) {
    FatalFormattedTableError("Could not write formatted table output '" +
                             filename + "'.");
  }
}

void CheckedFormattedTableFlush(std::FILE *output, const std::string &filename) {
  if (std::fflush(output) != 0) {
    FatalFormattedTableError("Could not flush formatted table output '" +
                             filename + "'.");
  }
}

void CheckedFormattedTableClose(std::FILE *output, const std::string &filename) {
  if (std::fclose(output) != 0) {
    FatalFormattedTableError("Could not close formatted table output '" +
                             filename + "'.");
  }
}

}  // namespace

//----------------------------------------------------------------------------------------
// ctor: also calls BaseTypeOutput base class constructor

FormattedTableOutput::FormattedTableOutput(ParameterInput *pin, Mesh *pm,
                                           OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  // check that 1D slice specified, otherwise issue warning and quit
  if (pm->multi_d) {
    if (!(out_params.slice1) && !(out_params.slice2)) {
      FatalFormattedTableError("Formatted table outputs can only contain 1D slices. "
                               "Please add additional slice planes.");
    }
  }
  if (pm->three_d) {
    if ((!(out_params.slice2) && !(out_params.slice3)) ||
        (!(out_params.slice1) && !(out_params.slice3))) {
      FatalFormattedTableError("Formatted table outputs can only contain 1D slices. "
                               "Please add additional slice planes.");
    }
  }
  // create directories for outputs. Comments in binary.cpp constructor explain why
  output_file_utils::EnsureDirectory("tab", 0775, "formatted table output",
                                     FatalFormattedTableError);
}

//----------------------------------------------------------------------------------------
//! \fn void FormattedTableOutput:::WriteOutputFile(Mesh *pm)
//  \brief writes output_data_ to file in tabular format using C style std::fprintf

void FormattedTableOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  // create filename: "tab/file_basename" + "." + "file_id" + "." + XXXXX + ".tab"
  // where XXXXX = file_number with a minimum width of 5 digits
  std::string fname;
  std::string number = output_file_utils::FormatSequence(
      out_params.file_number, "formatted table output", FatalFormattedTableError);

  fname.assign("tab/");
  fname.append(out_params.file_basename);
  fname.append(".");
  fname.append(out_params.file_id);
  fname.append(".");
  fname.append(number);
  fname.append(".tab");

  // master rank creates file and writes header (even though it may not have any actual
  // data to write below)
  if (global_variable::my_rank == 0) {
    FILE *pfile;
    if ((pfile = std::fopen(fname.c_str(),"w")) == nullptr) {
      FatalFormattedTableError("Output file '" + fname + "' could not be opened.");
    }

    // print file header
    CheckedFormattedTablePrint(pfile, fname, "# Athena++ data at time=%e", pm->time);
    CheckedFormattedTablePrint(pfile, fname, "  cycle=%d \n", pm->ncycle);

    // write one of x1, x2, x3 column headers
    CheckedFormattedTablePrint(pfile, fname, "# gid  ");
    if (!(out_params.slice1)) {
      CheckedFormattedTablePrint(pfile, fname, " i       x1v     ");
    }
    if (!(out_params.slice2)) {
      CheckedFormattedTablePrint(pfile, fname, " j       x2v     ");
    }
    if (!(out_params.slice3)) {
      CheckedFormattedTablePrint(pfile, fname, " k       x3v     ");
    }

    // write data col headers from outvars vector
    for (auto it : outvars) {
      CheckedFormattedTablePrint(pfile, fname, "    %s     ", it.label.c_str());
    }
    CheckedFormattedTablePrint(pfile, fname, "\n"); // terminate line
    CheckedFormattedTableClose(pfile, fname);   // don't forget to close the output file
  }
#if MPI_PARALLEL_ENABLED
  mpi_utils::CheckMpi(MPI_Barrier(MPI_COMM_WORLD),
                      "MPI_Barrier after formatted table header publication");
#endif

  // now all ranks open file and append data
  FILE *pfile;
  if ((pfile = std::fopen(fname.c_str(),"a")) == nullptr) {
    FatalFormattedTableError("Output file '" + fname + "' could not be opened.");
  }
  for (int r=0; r<global_variable::nranks; ++r) {
    // MPI ranks append data one-at-a-time in order, due to MPI_Barrier at end of loop
    // This could be slow for very large numbers of ranks, however this is not a regime
    // where .tab files are expected to be used very much.
    if (r == global_variable::my_rank) {
      // loop over output MeshBlocks, output all data
      int nout_vars = outvars.size();
      int nout_mbs = (outmbs.size());
      for (int m=0; m<nout_mbs; ++m) {
        auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
        auto &size  = pm->pmb_pack->pmb->mb_size;
        MeshBlock* pmb = pm->pmb_pack->pmb;
        int idx = pm->FindMeshBlockIndex(outmbs[m].mb_gid);
        int &is = indcs.is;
        int &js = indcs.js;
        int &ks = indcs.ks;
        int &ois = outmbs[m].ois;
        int &oie = outmbs[m].oie;
        int &ojs = outmbs[m].ojs;
        int &oje = outmbs[m].oje;
        int &oks = outmbs[m].oks;
        int &oke = outmbs[m].oke;
        Real &x1min = size.h_view(idx).x1min;
        Real &x1max = size.h_view(idx).x1max;
        Real &x2min = size.h_view(idx).x2min;
        Real &x2max = size.h_view(idx).x2max;
        Real &x3min = size.h_view(idx).x3min;
        Real &x3max = size.h_view(idx).x3max;
        int &nx1 = indcs.nx1;
        int &nx2 = indcs.nx2;
        int &nx3 = indcs.nx3;
        for (int k=oks; k<=oke; ++k) {
          for (int j=ojs; j<=oje; ++j) {
            for (int i=ois; i<=oie; ++i) {
              CheckedFormattedTablePrint(pfile, fname, "%05d", pmb->mb_gid.h_view(idx));
              // write x1, x2, x3 indices and coordinates
              if (oie != ois) {
                // note extra space for formatting
                CheckedFormattedTablePrint(pfile, fname, " %04d", i);
                Real x1cc = CellCenterX(i-is,nx1,x1min,x1max);
                CheckedFormattedTablePrint(pfile, fname, out_params.data_format.c_str(),
                                           x1cc);
              }
              if (oje != ojs) {
                // note extra space for formatting
                CheckedFormattedTablePrint(pfile, fname, " %04d", j);
                Real x2cc = CellCenterX(j-js,nx2,x2min,x2max);
                CheckedFormattedTablePrint(pfile, fname, out_params.data_format.c_str(),
                                           x2cc);
              }
              if (oke != oks) {
                // note extra space for formatting
                CheckedFormattedTablePrint(pfile, fname, " %04d", k);
                Real x3cc = CellCenterX(k-ks,nx3,x3min,x3max);
                CheckedFormattedTablePrint(pfile, fname, out_params.data_format.c_str(),
                                           x3cc);
              }

              // write each output variable on same line
              for (int n=0; n<nout_vars; ++n) {
                CheckedFormattedTablePrint(pfile, fname, out_params.data_format.c_str(),
                                           outarray(n,m,k-oks,j-ojs,i-ois));
              }
              CheckedFormattedTablePrint(pfile, fname, "\n"); // terminate line
            }
          }
        }
      }  // end loop over MeshBlocks
    }
    CheckedFormattedTableFlush(pfile, fname);
#if MPI_PARALLEL_ENABLED
    mpi_utils::CheckMpi(MPI_Barrier(MPI_COMM_WORLD),
                        "MPI_Barrier for ordered formatted table append");
#endif
  }

  CheckedFormattedTableClose(pfile, fname);   // don't forget to close the output file

  // increment counters
  out_params.file_number = output_file_utils::AdvanceFileNumber(
      out_params.file_number, "formatted table output", FatalFormattedTableError);
  if (out_params.last_time < 0.0) {
    out_params.last_time = pm->time;
  } else {
    out_params.last_time += out_params.dt;
  }
  // store filenumber and time into ParameterInput for restarts
  pin->SetInteger(out_params.block_name, "file_number", out_params.file_number);
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);

  return;
}
