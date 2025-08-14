#ifndef UTILS_SN_HPP_
#define UTILS_SN_HPP_

KOKKOS_INLINE_FUNCTION
Real GetNthSNTime(Real M, Real t_start, Real unit_time, int n) {
    // Constants from FIRE-3 Paper
    const Real a1 = 0.39;
    const Real a2 = 0.51;
    const Real a3 = 0.18;
    const Real a4 = 0.0083;
    const Real t1 = 3.7;  // Myr
    const Real t2 = 7.0;  // Myr
    const Real t3 = 44.0; // Myr
    const Real s1 = log(a2/a1)/log(t2/t1);
    const Real s2 = log(a3/a2)/log(t3/t2);
    const Real s3 = -1.1;
    const Real s_to_Myr = 3.171e-14;

    // Helper function to calculate cumulative number of SNe at time t
    auto n_sn = [&](Real t) -> Real {
        if (t < t1) {
            return 0.0;
        } else if (t <= t2) {
            return ((a1*t1)/(s1+1))*(pow(t/t1,s1+1)-1)/1000.0;
        } else if (t <= t3) {
            Real n_t2 = ((a1*t1)/(s1+1))*(pow(t2/t1,s1+1)-1)/1000.0;
            return n_t2 + ((a2*t2)/(s2+1))*(pow(t/t2,s2+1)-1)/1000.0;
        } else {
            Real n_t2 = ((a1*t1)/(s1+1))*(pow(t2/t1,s1+1)-1)/1000.0;
            Real n_t3 = n_t2 + ((a2*t2)/(s2+1))*(pow(t3/t2,s2+1)-1)/1000.0;
            return n_t3 + ((a4*t3)/(s3+1))*(pow(t/t3,s3+1)-1)/1000.0;
        }
    };

    // Calculate transition points
    Real n1 = n_sn(t2) * M;  // End of first phase
    Real n2 = n_sn(t3) * M;  // End of second phase

    // Calculate time in Myr based on which phase the nth SN falls in
    Real t_Myr;
    if (n < n1) {
        // First phase (Core-Collapse, early)
        t_Myr = t1 * pow(((1000.0 * n * (s1+1)) / (M * a1 * t1)) + 1, 1.0/(s1+1));
    } else if (n < n2) {
        // Second phase (Core-Collapse, late)
        t_Myr = t2 * pow(((1000.0 * (n - n1) * (s2+1)) / (M * a2 * t2)) + 1, 1.0/(s2+1));
    } else {
        // Third phase (Type Ia)
        t_Myr = t3 * pow(((1000.0 * (n - n2) * (s3+1)) / (M * a4 * t3)) + 1, 1.0/(s3+1));
    }

    // Convert from Myr to code units and add starting time
    return t_start + (t_Myr / (unit_time * s_to_Myr));
}

#endif // UTILS_SN_HPP_
