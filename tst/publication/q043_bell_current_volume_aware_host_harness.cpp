//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================

#define Q043_BELL_CURRENT_VOLUME_AWARE_HOST_CONTRACT 1
#include "../../src/pgen/tests/q043_bell_current_volume_aware.cpp"

#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char **argv) {
  using q043_bell_current_volume_aware::DepositedJOverC;
  using q043_bell_current_volume_aware::EigenmodeAtPhase;
  using q043_bell_current_volume_aware::HasRequiredDepositedJOverC;
  using q043_bell_current_volume_aware::HasRequiredSpeciesChargeOverMass;
  using q043_bell_current_volume_aware::ModeParameters;
  using q043_bell_current_volume_aware::RequiredDepositQScale;
  using q043_bell_current_volume_aware::RequiredDepositedJOverC;
  using q043_bell_current_volume_aware::RootCellVolume;
  using q043_bell_current_volume_aware::SourceMode;
  using q043_bell_current_volume_aware::SourceModeAmplitudeIsValid;
  using q043_bell_current_volume_aware::SourceModeSpeciesMassIsValid;

  if (argc == 6) {
    const std::string mode_name(argv[1]);
    SourceMode mode;
    if (mode_name == "uniform_current_oracle") {
      mode = SourceMode::uniform_current_oracle;
    } else if (mode_name == "corrected_linear_eigenmode") {
      mode = SourceMode::corrected_linear_eigenmode;
    } else {
      std::cout << "invalid\n";
      return 2;
    }
    const double amplitude = std::stod(argv[2]);
    const double species_mass = std::stod(argv[3]);
    const double species_charge = std::stod(argv[4]);
    const double expected_q_over_mc = std::stod(argv[5]);
    const bool valid = SourceModeAmplitudeIsValid(mode, amplitude) &&
        SourceModeSpeciesMassIsValid(mode, species_mass) &&
        HasRequiredSpeciesChargeOverMass(
            species_mass, species_charge, expected_q_over_mc);
    std::cout << (valid ? "valid" : "invalid") << '\n';
    return valid ? 0 : 2;
  }

  const double b_g = 1.0;
  const double k0 = 2.0*std::acos(-1.0);
  const double charge = k0*1.0e-6;
  const double v_cr = 2.5;
  const double light_speed = 2500.0;
  std::cout << std::setprecision(17);
  struct RootGrid {
    int dimension;
    double x1_extent;
    int nx1;
    double x2_extent;
    int nx2;
    double x3_extent;
    int nx3;
  };
  const RootGrid root_grids[] = {
    {1, 1.0, 32, 1.0, 4, 1.0, 1},
    {1, 1.0, 64, 1.0, 4, 1.0, 1},
    {2, std::sqrt(5.0), 64, std::sqrt(1.25), 32, 1.0, 1},
    {2, std::sqrt(5.0), 128, std::sqrt(1.25), 64, 1.0, 1},
    {3, std::sqrt(21.0), 128, std::sqrt(5.25), 64, std::sqrt(1.3125), 32},
    {3, std::sqrt(21.0), 256, std::sqrt(5.25), 128, std::sqrt(1.3125), 64},
  };
  for (const auto &grid : root_grids) {
    const double root_cell_volume = RootCellVolume(
        grid.x1_extent, grid.nx1, grid.x2_extent, grid.nx2,
        grid.x3_extent, grid.nx3);
    const double qscale = RequiredDepositQScale(
        1.0, charge, v_cr, root_cell_volume, b_g, k0);
    std::cout << "root_grid " << grid.dimension << ' ' << grid.nx1 << ' '
              << grid.nx2 << ' ' << grid.nx3 << ' ' << root_cell_volume << ' '
              << qscale << ' '
              << DepositedJOverC(1.0, qscale, charge, v_cr, root_cell_volume) << ' '
              << HasRequiredDepositedJOverC(
                     1.0, qscale, charge, v_cr, root_cell_volume, b_g, k0)
              << '\n';
    std::cout << "light_speed_multiplied " << grid.dimension << ' '
              << root_cell_volume << ' '
              << DepositedJOverC(
                     1.0, qscale*light_speed, charge, v_cr, root_cell_volume) << ' '
              << HasRequiredDepositedJOverC(
                     1.0, qscale*light_speed, charge, v_cr,
                     root_cell_volume, b_g, k0)
              << '\n';
  }
  const double decomposition_root_cell_volume =
      RootCellVolume(std::sqrt(5.0), 64, std::sqrt(1.25), 32, 1.0, 1);
  for (const int meshblock_nx1 : {16, 32, 64}) {
    for (const int meshblock_nx2 : {8, 16, 32}) {
      std::cout << "decomposition " << meshblock_nx1 << ' ' << meshblock_nx2
                << ' ' << decomposition_root_cell_volume << '\n';
    }
  }
  for (const double artificial_light_speed : {25.0, 2500.0, 250000.0}) {
    std::cout << "target " << artificial_light_speed << ' '
              << RequiredDepositedJOverC(b_g, k0) << '\n';
  }

  const ModeParameters parameters = {
    2, 0.4, 1.0e-6, 1.0, 1.0, 1.0, 1.0, k0
  };
  const auto sample = EigenmodeAtPhase(parameters, 0.37);
  std::cout << "shared_carrier " << sample.magnetic.x1 << ' '
            << sample.magnetic.x2 << ' ' << sample.magnetic.x3 << ' '
            << sample.velocity.x1 << ' ' << sample.velocity.x2 << ' '
            << sample.velocity.x3 << '\n';

  const ModeParameters uniform_parameters = {
    2, 0.4, 0.0, 1.0, 1.0, 1.0, 1.0, k0
  };
  const auto uniform = EigenmodeAtPhase(uniform_parameters, 0.37);
  std::cout << "uniform_carrier " << uniform.magnetic.x1 << ' '
            << uniform.magnetic.x2 << ' ' << uniform.magnetic.x3 << ' '
            << uniform.velocity.x1 << ' ' << uniform.velocity.x2 << ' '
            << uniform.velocity.x3 << '\n';
  std::cout << "source_modes "
            << SourceModeAmplitudeIsValid(SourceMode::uniform_current_oracle, 0.0) << ' '
            << SourceModeAmplitudeIsValid(SourceMode::uniform_current_oracle, 1.0e-6)
            << ' '
            << SourceModeAmplitudeIsValid(SourceMode::corrected_linear_eigenmode, 0.0)
            << ' '
            << SourceModeAmplitudeIsValid(
                   SourceMode::corrected_linear_eigenmode, 1.0e-6)
            << '\n';
  std::cout << "source_mode_masses "
            << SourceModeSpeciesMassIsValid(SourceMode::uniform_current_oracle, 2.0)
            << ' '
            << SourceModeSpeciesMassIsValid(SourceMode::uniform_current_oracle, 7.0)
            << ' '
            << SourceModeSpeciesMassIsValid(SourceMode::corrected_linear_eigenmode, 2.0)
            << ' '
            << SourceModeSpeciesMassIsValid(SourceMode::corrected_linear_eigenmode, 1.0)
            << '\n';
}
