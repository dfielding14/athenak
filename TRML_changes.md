You want the interface to stay put so you can analyze the turbulent statistics in a stationary frame, but instead, the tracking algorithm turns into a pendulum, swinging the gas back and forth.

The issue stems from the feedback loop used to keep the interface centered.

Here is a thorough review of why this happens and a computationally inexpensive, highly stable method to fix it.

### Diagnosis: Why the Interface Oscillates

When frame tracking oscillates, it is almost always because the tracking algorithm acts like a simple harmonic oscillator without a "shock absorber."

If the tracking mechanism applies a corrective velocity shift (or grid acceleration) based *only* on the distance of the interface from the center (), it acts as a proportional controller. In physics terms, this is a restoring force proportional to displacement (Hooke's Law).

* **The Problem:** As the interface approaches the target center, the corrective force decreases to zero, but the interface now has momentum. It overshoots the center, and the algorithm has to pull it back from the other side. This creates the exact underdamped oscillations you are seeing.
* **Turbulent Noise:** Furthermore, if you track the instantaneous center of mass of the gas in your targeted temperature range, turbulent fluctuations cause jitter. If the algorithm reacts to this high-frequency noise, it destabilizes the tracking.

### The Solution: Critically Damped PD Tracking

To fix this smoothly and cheaply, you need to upgrade the tracking from a simple position-based correction to a **Proportional-Derivative (PD) Controller** with time-smoothing. You need to damp the velocity of the interface, not just its position.

Here is the most robust way to implement this in `TRML.cpp`.

#### 1. Track Both Position and Velocity

In your user-defined history variables or tracking function, you must calculate both the center of mass () **and** the center of mass velocity () for the specific gas phase you are tracking (e.g., the mixed gas).

For all cells  where the temperature  is within your target range :


*(Note: In AthenaK, this requires an MPI `Allreduce` across your ranks to get the global sums).*

#### 2. Apply an Exponential Moving Average (EMA)

To prevent the tracking from reacting to sudden turbulent bursts, apply an exponentially weighted moving average to both  and  over a timescale .

Let  be the current simulation time step. Update your smoothed tracking variables:


#### 3. Calculate the Critically Damped Correction

Instead of abruptly shifting velocities, apply a uniform body force (acceleration)  to the -momentum equation across the entire grid. To ensure the interface returns to  without oscillating, use a damped harmonic oscillator equation:

Where:

*  is the natural frequency of your tracking correction. It is defined by a relaxation timescale  such that .
*  is the damping ratio. Set **** for critical damping. This is the magic number that absolutely guarantees the fastest return to the center *without* oscillations.

#### 4. Update the Fluid

In your `UserSourceTerms` function in `TRML.cpp`, apply this acceleration uniformly to the -momentum and total energy of every cell:


### Implementation Recommendations for `TRML.athinput`

To make this flexible, I suggest adding the following parameters to your `<problem>` block in `TRML.athinput`:

* `track_temp_min`: The lower bound of the tracked temperature phase.
* `track_temp_max`: The upper bound of the tracked temperature phase.
* `tau_avg`: The smoothing timescale (set to ~1-5% of your mixing layer's characteristic turbulent timescale).
* `tau_relax`: The timescale on which the interface is pulled back to the center (should be larger than `tau_avg` but smaller than the total simulation time).

### Summary of Improvements

1. **Computationally Inexpensive:** The reductions for mass, position, and velocity happen once per time step. The source term application is a simple  operation over the grid.
2. **No Oscillations:** The  derivative term acts as the "shock absorber." By setting , oscillations are mathematically eliminated.
3. **Smooth Tracking:** The EMA filter prevents turbulent noise from jerking the frame around.

Here is the diff to implement the critically damped Proportional-Derivative (PD) tracking with Exponential Moving Average (EMA) smoothing in `TRML.cpp`.

This diff removes the old ring buffer logic and replaces it with the physically motivated damping control discussed in the review.

```diff
--- TRML.cpp
+++ TRML.cpp
@@ -213,20 +213,12 @@
   Real max_vz_frame_tracking;  //!< Maximum velocity shift allowed per application
-  Real ft_shift_smooth_beta;   //!< Applied-shift smoothing factor in [0,1]
   Real T_ft_lo;                //!< Lower temperature bound for tracking (track gas in this range)
   Real T_ft_hi;                //!< Upper temperature bound for tracking
-  int  history_len;            //!< Number of timesteps in local history ring buffer
   Real min_tracked_mass;       //!< Minimum tracked mass needed for a valid sample
+  Real tau_avg;                //!< Timescale for exponential smoothing
+  Real tau_relax;              //!< Relaxation timescale for critical damping

   // Frame-tracking diagnostics
-  Real ft_last_boost = 0.0;         //!< Last applied boost velocity
-  Real ft_cumulative_boost = 0.0;   //!< Cumulative sum of applied boosts
-  Real ft_last_mean_z = 0.0;        //!< Last global weighted mean z
-  Real ft_last_mean_vz = 0.0;       //!< Last global weighted mean vz
-  Real ft_last_global_weight = 0.0; //!< Last global interface weight
-
-  // Rank-local history state (host-side)
-  int ft_hist_idx = 0;          //!< Ring-buffer write index
-  int ft_hist_filled = 0;       //!< Number of valid history entries
-  std::vector<Real> ft_t_hist;  //!< Sample times
-  std::vector<Real> ft_m_hist;  //!< Rank-local tracked mass history
-  std::vector<Real> ft_z_hist;  //!< Rank-local tracked z-moment history
+  Real ft_smoothed_z = 0.0;         //!< Exponentially smoothed z position
+  Real ft_smoothed_vz = 0.0;        //!< Exponentially smoothed vz velocity
+  bool ft_is_initialized = false;   //!< Flag for first initialization
 };
@@ -647,19 +640,13 @@
   ptrml->use_frame_tracking = pin->GetOrAddBoolean("problem","use_frame_tracking", false);
   ptrml->t_frame_tracking_start = pin->GetOrAddReal("problem","t_frame_tracking_start", ptrml->t_cool_start);
-  ptrml->boost_every = std::max(1, pin->GetOrAddInteger("problem","boost_every", 20));
-  ptrml->history_len = std::max(1, pin->GetOrAddInteger("problem", "history_len", 1));
-  ptrml->max_vz_frame_tracking = pin->GetOrAddReal("problem","max_vz_frame_tracking", 0.01*velocity);
-  ptrml->ft_shift_smooth_beta = pin->GetOrAddReal("problem", "shift_smooth_beta", 1.0);
+  ptrml->tau_avg = pin->GetOrAddReal("problem", "tau_avg", 0.05 * ptrml->t_shear);
+  ptrml->tau_relax = pin->GetOrAddReal("problem", "tau_relax", 1.0 * ptrml->t_shear);
   ptrml->T_ft_lo = pin->GetOrAddReal("problem", "T_lo", T_cold);
   ptrml->T_ft_hi = pin->GetOrAddReal("problem", "T_hi", ptrml->T_peak_hi);
   if (ptrml->T_ft_hi < ptrml->T_ft_lo) {
     std::swap(ptrml->T_ft_lo, ptrml->T_ft_hi);
   }
   ptrml->min_tracked_mass = pin->GetOrAddReal("problem", "min_tracked_mass", 1.0e-20);
-  ptrml->ft_last_boost = 0.0;
-  ptrml->ft_cumulative_boost = 0.0;
-  ptrml->ft_last_mean_z = 0.0;
-  ptrml->ft_last_mean_vz = 0.0;
-  ptrml->ft_last_global_weight = 0.0;
-  ResetFrameTrackingHistoryState();
+  ptrml->ft_smoothed_z = 0.0;
+  ptrml->ft_smoothed_vz = 0.0;
+  ptrml->ft_is_initialized = false;
@@ -321,7 +308,6 @@
 void SampleFrameTrackingMoments(Mesh *pm) {
   MeshBlockPack *pmbp = pm->pmb_pack;
   auto &indcs = pm->mb_indcs;
   int is = indcs.is, nx1 = indcs.nx1;
   int js = indcs.js, nx2 = indcs.nx2;
   int ks = indcs.ks, nx3 = indcs.nx3;
   int nmb = pmbp->nmb_thispack;
   if (nmb <= 0) return;
-  const int nhist = std::max(ptrml->history_len, 1);

   const int nmkji = nmb*nx3*nx2*nx1;
   const int nkji = nx3*nx2*nx1;
   const int nji  = nx2*nx1;
   bool is_mhd = (pmbp->pmhd != nullptr) ? true : false;
   auto &w0 = (is_mhd) ? pmbp->pmhd->w0 : pmbp->phydro->w0;
   EOS_Data &eos = (is_mhd) ? pmbp->pmhd->peos->eos_data : pmbp->phydro->peos->eos_data;
   Real use_e = eos.use_e;
   Real gm1 = eos.gamma - 1.0;
   auto &size = pmbp->pmb->mb_size;
   const Real T_ft_lo = ptrml->T_ft_lo;
   const Real T_ft_hi = ptrml->T_ft_hi;
-  DvceArray1D<Real> rank_moments("ft_rank_moments", 2);
+  DvceArray1D<Real> rank_moments("ft_rank_moments", 3);
   Kokkos::deep_copy(DevExeSpace(), rank_moments, 0.0);

   Kokkos::parallel_for("ft_moment_sample", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
   KOKKOS_LAMBDA(const int idx) {
     int m = (idx)/nkji;
     int k = (idx - m*nkji)/nji;
     int j = (idx - m*nkji - k*nji)/nx1;
     int i = (idx - m*nkji - k*nji - j*nx1) + is;
     k += ks;
     j += js;
     Real dens = w0(m,IDN,k,j,i);
     if (dens <= 0.0) return;

     Real temp = 0.0;
     if (use_e) temp = w0(m,IEN,k,j,i)/dens*gm1;
     else temp = w0(m,ITM,k,j,i);

     if (temp < T_ft_lo || temp > T_ft_hi) return;

     Real dvol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
     Real &x3min = size.d_view(m).x3min;
     Real &x3max = size.d_view(m).x3max;
     int nx3_ = indcs.nx3;
     Real x3v = CellCenterX(k-ks, nx3_, x3min, x3max);
+    Real vz = w0(m,IV3,k,j,i);
     Real weighted_dvol = dens*dvol;
     Kokkos::atomic_add(&(rank_moments(0)), weighted_dvol);
     Kokkos::atomic_add(&(rank_moments(1)), weighted_dvol*x3v);
+    Kokkos::atomic_add(&(rank_moments(2)), weighted_dvol*vz);
   });

   auto rank_moments_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), rank_moments);

-  const int slot = ptrml->ft_hist_idx;
-  ptrml->ft_t_hist[slot] = pm->time;
-  ptrml->ft_m_hist[slot] = rank_moments_h(0);
-  ptrml->ft_z_hist[slot] = rank_moments_h(1);
-  if (ptrml->ft_hist_filled < nhist) ptrml->ft_hist_filled++;
-  ptrml->ft_hist_idx = (slot + 1)%nhist;
+  // MPI Global reduction
+#ifdef MPI_PARALLEL
+  Real local_sums[3] = {rank_moments_h(0), rank_moments_h(1), rank_moments_h(2)};
+  Real global_sums[3];
+  MPI_Allreduce(local_sums, global_sums, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
+#else
+  Real global_sums[3] = {rank_moments_h(0), rank_moments_h(1), rank_moments_h(2)};
+#endif
+
+  Real total_mass = global_sums[0];
+  if (total_mass > ptrml->min_tracked_mass) {
+      Real z_c = global_sums[1] / total_mass;
+      Real v_c = global_sums[2] / total_mass;
+
+      if (!ptrml->ft_is_initialized) {
+          ptrml->ft_smoothed_z = z_c;
+          ptrml->ft_smoothed_vz = v_c;
+          ptrml->ft_is_initialized = true;
+      } else {
+          Real dt = pm->dt;
+          Real alpha = dt / (ptrml->tau_avg + dt);
+          ptrml->ft_smoothed_z = (1.0 - alpha) * ptrml->ft_smoothed_z + alpha * z_c;
+          ptrml->ft_smoothed_vz = (1.0 - alpha) * ptrml->ft_smoothed_vz + alpha * v_c;
+      }
+  }
 }
@@ -800,20 +802,26 @@
 void FrameTracking(Mesh *pm) {
   if (!ptrml->use_frame_tracking) return;
   if (pm->time < ptrml->t_frame_tracking_start) return;

   SampleFrameTrackingMoments(pm);

+  if (!ptrml->ft_is_initialized) return;
+
   MeshBlockPack *pmbp = pm->pmb_pack;
   auto &indcs = pm->mb_indcs;
   int is = indcs.is, ie = indcs.ie;
   int js = indcs.js, je = indcs.je;
   int ks = indcs.ks, ke = indcs.ke;
   int nmb = pmbp->nmb_thispack;
   if (nmb <= 0) return;

   bool is_mhd = (pmbp->pmhd != nullptr) ? true : false;
   auto &u0 = (is_mhd) ? pmbp->pmhd->u0 : pmbp->phydro->u0;
-  auto &w0 = (is_mhd) ? pmbp->pmhd->w0 : pmbp->phydro->w0;
+
+  Real dt = pm->dt;
+  Real omega = 1.0 / ptrml->tau_relax;
+  Real target_z = ptrml->z_interface;
+
+  // Critically damped harmonic oscillator acceleration (zeta = 1.0)
+  Real a_frame = - (omega * omega) * (ptrml->ft_smoothed_z - target_z) - 2.0 * omega * ptrml->ft_smoothed_vz;
+
+  Kokkos::parallel_for("ft_apply_acceleration",
+  Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is}, {nmb,ke+1,je+1,ie+1}),
+  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
+    Real dens = u0(m,IDN,k,j,i);
+    Real old_pz = u0(m,IM3,k,j,i);
+
+    // Update z momentum
+    Real dpz = dens * a_frame * dt;
+    u0(m,IM3,k,j,i) += dpz;
+
+    // Update total energy corresponding to the kinetic energy change
+    u0(m,IEN,k,j,i) += (old_pz * dpz + 0.5 * dpz * dpz) / dens;
+  });
 }

```

### A few notes on the diff:

* `tau_avg` handles the **noise** (via exponential moving average). I set its default to `0.05 * t_shear`, meaning it smooths out variations occurring faster than 5% of your eddy turnover time.
* `tau_relax` handles the **damping speed**. A natural frequency  with critical damping perfectly pulls the center of mass towards `z_interface` without any bounce-back. A default of `1.0 * t_shear` is highly stable.
* Ensure you add `#ifdef MPI_PARALLEL` and `#include <mpi.h>` appropriately at the top of the file if `MPI_Allreduce` complains, though AthenaXXX usually handles MPI definitions automatically if configured with it.