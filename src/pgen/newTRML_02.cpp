#include <iostream> // cout

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "pgen.hpp"
#include "globals.hpp"

Real glob_rho_cold;
Real glob_rho_hot;
Real glob_pres;
Real glob_t_cool_min;
Real glob_T_cutoff_over_T_cold;
Real glob_beta;

Real glob_cooling_rate;

void CoolingSrc(Mesh* pm, Real bdt);
void CoolingSrc2(Mesh* pm, Real bdt);
void CoolingSrc3(Mesh* pm, Real bdt);
void HistoryOutput(HistoryData *pdata, Mesh *pm);

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) 
{
    //for restarting
    if (restart) return;
    //Grab meshblock pointers
    MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
    auto &indcs = pmy_mesh_->mb_indcs;

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

    glob_cooling_rate = 0.0;

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

    user_srcs_func = CoolingSrc3;
    user_hist_func = HistoryOutput;

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
        Real dens = u0(m,IDN,k,j,i);
        Real eint = u0(m,IEN,k,j,i)-0.5*(SQR(u0(m,IM1,k,j,i))+SQR(u0(m,IM2,k,j,i))+SQR(u0(m,IM3,k,j,i)))/dens;
        Real temp = eint/dens*gm1;

        Real temp_at_t_cool_min = (beta>0)?T_cold:(T_cutoff_over_T_cold*T_cold);
        Real Edot = eint*pow(temp/temp_at_t_cool_min,-beta)/t_cool_min;

        if((temp/T_cold)<1 || (temp/T_cold)>T_cutoff_over_T_cold) Edot=0.0;

        u0(m,IEN,k,j,i) -= bdt * Edot;

    });
}

void CoolingSrc2(Mesh* pm, Real bdt)
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
        //This is the exact solution for the cooling
        Real Edot; 
        if(beta!=0.0)
            Edot = eint*pow(1-pow(temp/temp_at_t_cool_min,-beta)*bdt/t_cool_min*beta,1/beta)-eint;
        else
            Edot = eint*exp(-bdt/t_cool_min)-eint;

        if((temp/T_cold)<1 || (temp/T_cold)>T_cutoff_over_T_cold) Edot=0.0;

        u0(m,IEN,k,j,i) += Edot;

    });
}

void CoolingSrc3(Mesh* pm, Real bdt)
{
    MeshBlockPack *pmbp = pm->pmb_pack;

    auto &indcs = pm->mb_indcs;
    auto &size = pmbp->pmb->mb_size;
    int &is = indcs.is; int &ie = indcs.ie;
    int &js = indcs.js; int &je = indcs.je;
    int &ks = indcs.ks; int &ke = indcs.ke;
    int nx1 = indcs.nx1;
    int nx2 = indcs.nx2;
    int nx3 = indcs.nx3;
    const int nmkji = (pmbp->nmb_thispack)*nx3*nx2*nx1;
    const int nkji = nx3*nx2*nx1;
    const int nji  = nx2*nx1;

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

    Real net_cool=0.0;
    Kokkos::parallel_reduce("cooling_src", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
    KOKKOS_LAMBDA(const int &idx, Real &net_cooling) 
    {
        int m = (idx)/nkji;
        int k = (idx - m*nkji)/nji;
        int j = (idx - m*nkji - k*nji)/nx1;
        int i = (idx - m*nkji - k*nji - j*nx1) + is;
        k += ks;
        j += js;

        Real dens = w0(m,IDN,k,j,i);
        Real eint = w0(m,IEN,k,j,i)-0.5*(SQR(w0(m,IVX,k,j,i))+SQR(w0(m,IVY,k,j,i))+SQR(w0(m,IVZ,k,j,i)))*dens;
        Real temp = eint/dens*gm1;

        Real temp_at_t_cool_min = (beta>0)?T_cold:(T_cutoff_over_T_cold*T_cold);
        //This is the exact solution for the cooling
        Real deltaE; 
        if(beta!=0.0)
            deltaE = eint*pow(1-pow(temp/temp_at_t_cool_min,-beta)*bdt/t_cool_min*beta,1/beta)-eint;
        else
            deltaE = eint*exp(-bdt/t_cool_min)-eint;

        if((temp/T_cold)<1 || (temp/T_cold)>T_cutoff_over_T_cold) deltaE=0.0;

        net_cooling += deltaE/bdt*size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;

        u0(m,IEN,k,j,i) += deltaE;
    }, Kokkos::Sum<Real>(net_cool));
    #if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, &net_cool, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    #endif
    glob_cooling_rate = net_cool;
}

void HistoryOutput(HistoryData *pdata, Mesh *pm)
{
    pdata->nhist = 1;
    pdata->label[0] = "cooling_rate";
    if (global_variable::my_rank == 0) 
    {
        pdata->hdata[0] = glob_cooling_rate;
    }
    else
    {
        pdata->hdata[0] = 0.0;
    }
}