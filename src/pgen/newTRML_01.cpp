#include <iostream> // cout

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "pgen.hpp"

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
    auto &size = pmbp->pmb->mb_size;

    Real rho_cold = pin->GetReal("problem","rho_cold");
    Real rho_hot = pin->GetReal("problem","rho_hot");
    Real pres = pin->GetReal("problem","pres");
    Real velocity = pin->GetReal("problem","velocity");
    // Real x1min = pin->GetReal("mesh", "x1min");
    // Real x1max = pin->GetReal("mesh", "x1max");
    // Real x2min = pin->GetReal("mesh", "x2min");
    // Real x2max = pin->GetReal("mesh", "x2max");
    // Real x3min = pin->GetReal("mesh", "x3min");
    // Real x3max = pin->GetReal("mesh", "x3max");

    if (pmbp->phydro != nullptr) {
        auto &u0 = pmbp->phydro->u0;
        EOS_Data &eos = pmbp->phydro->peos->eos_data;
        Real gm1 = eos.gamma - 1.0;

        // Set initial conditions
        par_for("pgen_turb", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
        KOKKOS_LAMBDA(int m, int k, int j, int i) 
        {
            Real &x1min = size.d_view(m).x1min;
            Real &x1max = size.d_view(m).x1max;
            Real &x2min = size.d_view(m).x2min;
            Real &x2max = size.d_view(m).x2max;
            Real &x3min = size.d_view(m).x3min;
            Real &x3max = size.d_view(m).x3max;
            
            Real coordx = CellCenterX(i-is,indcs.nx1,x1min,x1max);
            Real coordy = CellCenterX(j-js,indcs.nx2,x2min,x2max);
            Real coordz = CellCenterX(k-ks,indcs.nx3,x3min,x3max);
            Real dens = (coordz<=0?rho_cold:rho_hot);
            u0(m,IDN,k,j,i) = dens;
            u0(m,IM1,k,j,i) = dens*velocity*(coordz<=0?0.5:(-0.5));
            u0(m,IM2,k,j,i) = 0.0;
            u0(m,IM3,k,j,i) = 0.0;
            if (eos.is_ideal) {
                u0(m,IEN,k,j,i) = pres/gm1 +
                0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i)) +
                SQR(u0(m,IM3,k,j,i)))/u0(m,IDN,k,j,i);
            }

            for(int i0=0;i0<3;i0++)
            {
                for(int j0=0;j0<3;j0++)
                {
                    for(int k0=0;k0<3;k0++)
                    {
                        u0(m,IM3,k,j,i)+=velocity*0.01*sin(coordz*6.28*i0)*sin(coordy*6.28*j0)*sin(coordx*6.28*k0)*u0(m,IDN,k,j,i);
                    }
                }
            }
        });
    }

    return;
}