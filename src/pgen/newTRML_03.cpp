#include <iostream> // cout

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "pgen.hpp"

Real glob_rho_cold;
Real glob_rho_hot;
Real glob_pres;
Real glob_t_cool_min;
Real glob_T_cutoff_over_T_cold;
Real glob_beta;

void CoolingSrc(Mesh* pm, Real bdt);
void CoolingTimestep(Mesh* pm);

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) 
{
    //for restarting
    if (restart) return;
    //Grab meshblock pointers
    MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
    auto &indcs = pmy_mesh_->mb_indcs;
    EOS_Data &eos = pmbp->phydro->peos->eos_data;

    if(eos.use_e)
    {
        printf("Using E");
    }
    if(eos.use_t)
    {
        printf("Using T");
    }

    //Indices
    int &is = indcs.is; int &ie = indcs.ie;
    int &js = indcs.js; int &je = indcs.je;
    int &ks = indcs.ks; int &ke = indcs.ke;

    glob_rho_cold = pin->GetReal("problem","rho_cold");
    glob_rho_hot = pin->GetReal("problem","rho_hot");
    glob_pres = pin->GetReal("problem","pres");
    glob_t_cool_min = pin->GetReal("problem","t_cool_min");
    glob_T_cutoff_over_T_cold = pin->GetReal("problem","T_cutoff_over_T_cold");
    glob_beta = pin->GetReal("problem","beta");

    Real rho_cold = glob_rho_cold;
    Real rho_hot = glob_rho_hot;
    Real pres = glob_pres;
    Real t_cool_min = glob_t_cool_min;
    Real T_cutoff_over_T_cold = glob_T_cutoff_over_T_cold;
    Real beta = glob_beta;

    Real x1min = pin->GetReal("mesh", "x1min");
    Real x1max = pin->GetReal("mesh", "x1max");
    Real x2min = pin->GetReal("mesh", "x2min");
    Real x2max = pin->GetReal("mesh", "x2max");
    Real x3min = pin->GetReal("mesh", "x3min");
    Real x3max = pin->GetReal("mesh", "x3max");

    user_srcs_func = CoolingSrc;
    user_time_step_func = CoolingTimestep;

    if (pmbp->phydro != nullptr) {
        auto &u0 = pmbp->phydro->u0;
        EOS_Data &eos = pmbp->phydro->peos->eos_data;
        Real gm1 = eos.gamma - 1.0;
        Real dens = rho_cold/10.0;

        // Set initial conditions
        par_for("pgen_turb", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
        KOKKOS_LAMBDA(int m, int k, int j, int i) 
        {
            u0(m,IDN,k,j,i) = dens;
            u0(m,IM1,k,j,i) = 0.0;
            u0(m,IM2,k,j,i) = 0.0;
            u0(m,IM3,k,j,i) = 0.0;
            if (eos.is_ideal) {
                u0(m,IEN,k,j,i) = pres/gm1;
            }
        });
    }

    return;
}

void CoolingSrc(Mesh* pm, Real bdt)
{
    MeshBlockPack *pmbp = pm->pmb_pack;

    auto &indcs = pm->mb_indcs;
    int &is = indcs.is; int &ie = indcs.ie;
    int &js = indcs.js; int &je = indcs.je;
    int &ks = indcs.ks; int &ke = indcs.ke;
    auto &u0 = pmbp->phydro->u0;
    auto &w0 = pmbp->phydro->w0;
    EOS_Data &eos = pmbp->phydro->peos->eos_data;
    Real gm1 = eos.gamma - 1.0;

    Real rho_cold = glob_rho_cold;
    Real rho_hot = glob_rho_hot;
    Real pres = glob_pres;
    Real t_cool_min = glob_t_cool_min;
    Real T_cutoff_over_T_cold = glob_T_cutoff_over_T_cold;
    Real beta = glob_beta;

    Real T_cold = pres/rho_cold;
    Real T_hot = pres/rho_hot;

    par_for("pgen_turb", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) 
    {
        Real dens = w0(m,IDN,k,j,i);
        Real eint = w0(m,IEN,k,j,i)-0.5*(SQR(w0(m,IVX,k,j,i))+SQR(w0(m,IVY,k,j,i))+SQR(w0(m,IVZ,k,j,i)))*dens;
        Real temp = eint/dens*gm1;

        Real temp_at_t_cool_min = (beta>0)?T_cold:(T_cutoff_over_T_cold*T_cold);
        Real Edot = eint*pow(temp/temp_at_t_cool_min,-beta)/t_cool_min;

        if((temp/T_cold)<1 || (temp/T_cold)>T_cutoff_over_T_cold) Edot=0.0;

        u0(m,IEN,k,j,i) -= bdt * Edot;

    });
}

void CoolingTimestep(Mesh* pm)
{
    pm->pgen->dtnew = glob_t_cool_min;
}

