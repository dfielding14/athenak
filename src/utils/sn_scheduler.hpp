#ifndef UTILS_SN_HPP_
#define UTILS_SN_HPP_

class SNInj {
public:
    SNInj();
    ~SNInj();
    
    std::vector<Real> GetSNTimes(Real M, Real t_start, Real t_end, Real unit_time);

private:
    Real a1, a2, a3, a4;
    Real t1, t2, t3;
    Real s1, s2, s3;
    
    Real n_sn(Real t);
    Real get_sn_t(Real n, Real n1, Real n2, Real M);
};

// Constructor
SNInj::SNInj() {
    // From FIRE-3 Paper
    a1 = 0.39;
    a2 = 0.51;
    a3 = 0.18;
    a4 = 0.0083;

    t1 = 3.7; // Myr
    t2 = 7.0;
    t3 = 44.0;

    s1 = std::log(a2/a1)/log(t2/t1);
    s2 = std::log(a3/a2)/std::log(t3/t2);
    s3 = -1.1;
    
    return;
}

// Destructor
SNInj::~SNInj() {
    return;
}

// Number of SN in a given time
Real SNInj::n_sn(Real t) { 
    // Core-Collapse
    if (t < t1) {
        return 0.0;
    } else if (t <= t2) {
        return ((a1*t1)/(s1+1))*(std::pow(t/t1,s1+1)-1)/1000.0;
    } else if (t <= t3) {
        return n_sn(t2) + ((a2*t2)/(s2+1))*(std::pow(t/t2,s2+1)-1)/1000.0;
    }
    
    // Ia
    return n_sn(t3) + ((a4*t3)/(s3+1))*(std::pow(t/t3,s3+1)-1)/1000;
}

// Get the time of the nth SN              
Real SNInj::get_sn_t(Real n, Real n1, Real n2, Real M) {
    if (n < n1) {
        return t1*std::pow(((1000*n*(s1+1))/(M*a1*t1))+1,1/(s1+1));
    } else if (n < n2) {
        return t2*std::pow(((1000*(n-n1)*(s2+1))/(M*a2*t2))+1,1/(s2+1));
    } else {
        return t3*std::pow(((1000*(n-n2)*(s3+1))/(M*a4*t3))+1,1/(s3+1));
    }
    
    return 1e20; // Should never reach here
}

// Get list of SNe. t_start and t_end are in code units. Returns in code units.
std::vector<Real> SNInj::GetSNTimes(Real M, Real t_start, Real t_end, Real unit_time) {
    const Real s_to_Myr = 3.171e-14;

    Real n1 = n_sn(t2)*M;
    Real n2 = n_sn(t3)*M;
    Real n_end = n_sn((t_end-t_start)*unit_time*s_to_Myr)*M;

    std::vector<Real> times;
    times.reserve(static_cast<size_t>(n_end));
    for (int i=1; i < n_end; i++) {
        times.push_back(t_start+(get_sn_t(i, n1, n2, M)/(unit_time*s_to_Myr)));
    }
    
    times.push_back(1e20); // Dummy SN at t=inf
    
    return times;
}

#endif // UTILS_SN_HPP_
