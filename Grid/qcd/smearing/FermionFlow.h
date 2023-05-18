#pragma once

NAMESPACE_BEGIN(Grid);

template<class Gimpl, class GaugeAction, class SpinorField>
class FermionFlow: public GradientFlowBase<Gimpl,GaugeAction> {
private:
    int Nstep;
    int Nckpoints;
    RealD epsilon;
    mutable std::vector<typename Gimpl::GaugeField> ckpoint_fields;
    mutable std::vector<typename Gimpl::GaugeField> cached_fields;
    //add 2nd level cache, saving only a single block
    void evolve_step(typename Gimpl::GaugeField &U, RealD &tau) const;
    void evolve_step(typename Gimpl::GaugeField &U, SpinorField &chi, RealD &tau) const;
    void fermion_laplacian(typename Gimpl::GaugeField &U, SpinorField &chi_in, SpinorField &chi_out) const;
    void evolve_step_adjoint(typename Gimpl::GaugeField &U, SpinorField &chi, RealD &tau) const;


public:
    INHERIT_GIMPL_TYPES(Gimpl)

    FermionFlow(const RealD Epsilon, const int Nstep, const GaugeField &U, 
    unsigned int meas_interval = 1, int Nckpoints = 10): GradientFlowBase<Gimpl,GaugeAction>(meas_interval), Nstep(Nstep),
    Nckpoints(Nckpoints), epsilon(Epsilon), ckpoint_fields(Nckpoints,U.Grid()), cached_fields(Nstep/Nckpoints,U.Grid()) {}

    void smear(GaugeField& out, SpinorField& chi_out, const GaugeField& in, const SpinorField& chi_in) const;
    void smear(GaugeField& out, const GaugeField& in) const override {};
    void smear_adjoint(SpinorField& chi_out, const SpinorField& chi_in) const;
    
};

template<class Gimpl, class GaugeAction, class SpinorField>
void FermionFlow<Gimpl, GaugeAction, SpinorField>::fermion_laplacian(typename Gimpl::GaugeField &U, SpinorField &chi_in, SpinorField &chi_out) const {
    chi_out = 0;
    for (int mu = 0; mu < Nd; mu++) {
        GaugeLinkField Ulink(U.Grid());
        Ulink = PeekIndex<LorentzIndex>(U,mu);
        chi_out += Ulink*Cshift(chi_in,mu,1)-2.0*chi_in+adj(Gimpl::CshiftLink(Ulink,mu,-1))*Cshift(chi_in,mu,-1);
    }
    return;
}

template <class Gimpl, class GaugeAction, class SpinorField>
void FermionFlow<Gimpl, GaugeAction, SpinorField>::evolve_step(typename Gimpl::GaugeField &U, RealD &tau) const{
    GaugeField Z(U.Grid());
    GaugeField tmp(U.Grid());
    this->SG.deriv(U, Z);
    Z *= 0.25;                                  // Z0 = 1/4 * F(U)
    Gimpl::update_field(Z, U, -2.0*epsilon);    // U = W1 = exp(ep*Z0)*W0

    Z *= -17.0/8.0;
    this->SG.deriv(U, tmp); Z += tmp;                 // -17/32*Z0 +Z1
    Z *= 8.0/9.0;                               // Z = -17/36*Z0 +8/9*Z1
    Gimpl::update_field(Z, U, -2.0*epsilon);    // U_= W2 = exp(ep*Z)*W1

    Z *= -4.0/3.0;
    this->SG.deriv(U, tmp); Z += tmp;                 // 4/3*(17/36*Z0 -8/9*Z1) +Z2
    Z *= 3.0/4.0;                               // Z = 17/36*Z0 -8/9*Z1 +3/4*Z2
    Gimpl::update_field(Z, U, -2.0*epsilon);    // V(t+e) = exp(ep*Z)*W2
    tau += epsilon;
}


template<class Gimpl, class GaugeAction, class SpinorField>
void FermionFlow<Gimpl, GaugeAction, SpinorField>::evolve_step(typename Gimpl::GaugeField &U, SpinorField &chi, RealD &tau) const {
    GaugeField Z(U.Grid());
    GaugeField tmp(U.Grid());
    SpinorField phi0(U.Grid());
    SpinorField phi1(U.Grid());
    SpinorField phi2(U.Grid());
    //W0
    

    this->SG.deriv(U,Z);
    fermion_laplacian(U,chi,phi0);
    Z*=0.25;
    Gimpl::update_field(Z, U, -2.0*epsilon);
    phi1 = chi + 0.25*epsilon*phi0; //phi1 = phi0 + 1/4eps*laplace(W0)phi0
    //U is now W1
    
    Z *= -17.0/8.0;
    this->SG.deriv(U, tmp);
    
    

    fermion_laplacian(U,phi1,phi2); //phi2 is laplace(W1)phi1   
    //phi2:
    phi2 = chi + 8.0/9.0*epsilon*phi2-2.0/9.0*epsilon*phi0; //phi2 = phi0 + 8/9*eps*laplace(W1)phi1-2/9*epsilon*laplace(W0)phi0
    // std::cout << GridLogMessage << "Testing phi0 " << norm2(PeekIndex<SpinorIndex>(phi0,0)) << std::endl;
    Z += tmp;

    Z *= 8.0/9.0;
    Gimpl::update_field(Z, U, -2.0*epsilon);
    // U is now W2

    Z *= -4.0/3.0;
    this->SG.deriv(U, tmp);
    
    Z += tmp;

    Z *= 3.0/4.0;

    fermion_laplacian(U,phi2,phi0);
    
    chi = phi1 + 3.0/4.0*epsilon*phi0; ///phi3 = phi + 3/4*eps*laplace(W2)phi2
    

    Gimpl::update_field(Z, U, -2.0*epsilon);


    tau+=epsilon;

    return;


}

//evolve adjoint equation from t+epsilon to t (needed for e.g. chiral condensate)
//GaugeField passed to this function is situated at time t
template<class Gimpl, class GaugeAction, class SpinorField>
void FermionFlow<Gimpl, GaugeAction, SpinorField>::evolve_step_adjoint(typename Gimpl::GaugeField &U, SpinorField &chi, RealD &tau) const {
    SpinorField lambda3(U.Grid());
    SpinorField lambda2(U.Grid());
    SpinorField lambda1(U.Grid());
    SpinorField laplace_tmp(U.Grid());
  
    GaugeField W1(U.Grid());
    GaugeField W2(U.Grid());
    W2 = U;
    W1 = U;
    //lambda3 = chi(t+epsilon)
    lambda3 = chi;

    //prepare W2
    GaugeField Z(U.Grid());
    GaugeField tmp(U.Grid());
    this->SG.deriv(W2, Z);
    Z *= 0.25;                                  // Z0 = 1/4 * F(U)
    Gimpl::update_field(Z, W2, -2.0*epsilon);    // U = W1 = exp(ep*Z0)*W0

    Z *= -17.0/8.0;
    this->SG.deriv(W2, tmp); Z += tmp;                 // -17/32*Z0 +Z1
    Z *= 8.0/9.0;                               // Z = -17/36*Z0 +8/9*Z1
    Gimpl::update_field(Z, W2, -2.0*epsilon);    // U_= W2 = exp(ep*Z)*W1


    //lambda2 = 3/4*eps*laplace(W2)lambda3
    fermion_laplacian(W2,lambda3,laplace_tmp);
    lambda2 = 3.0/4.0*epsilon*laplace_tmp;
    
    //prepare W1
    this->SG.deriv(W1,tmp);
    tmp *= 0.25;
    Gimpl::update_field(tmp,W1, -2.0*epsilon);

    //lambda1 = lambda3 + 8/9*eps*laplace(W1)lambda2
    fermion_laplacian(W1,lambda2,laplace_tmp);
    lambda1 = lambda3 + 8.0/9.0*epsilon*laplace_tmp;
    
    //lambda0 = lambda1 + lambda2 + 1/4*eps*laplace(W0)(lambda1-8/9lambda2)
    lambda3 = lambda1 - 8.0/9.0*lambda2;
    fermion_laplacian(U,lambda3,laplace_tmp);

    chi = lambda1 + lambda2 + 0.25*epsilon*laplace_tmp;

    tau-=epsilon;
}

template<class Gimpl, class GaugeAction, class SpinorField>
void FermionFlow<Gimpl, GaugeAction, SpinorField>::smear(GaugeField& out, SpinorField& chi_out, const GaugeField& in, const SpinorField& chi_in) const {
    std::cout << GridLogMessage
        << "[FermionFlow] Nstep   : " << Nstep<< std::endl;
    std::cout << GridLogMessage
        << "[FermionFlow] epsilon   : " << epsilon << std::endl;
    std::cout << GridLogMessage
        << "[FermionFlow] full trajectory : " << Nstep * epsilon << std::endl;

    out = in;
    chi_out = chi_in;
    ckpoint_fields[0] = in;
    RealD taus = 0;
    for (unsigned int step = 1; step <= Nstep; step++) { //step indicates the number of smearing steps applied at the time of measurement
    auto start = std::chrono::high_resolution_clock::now();
    evolve_step(out, chi_out, taus);
    std::cout << GridLogMessage << "[FermionFlow]  chi norm after step = " << norm2(PeekIndex<SpinorIndex>(chi_out,0)) << std::endl;
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
#ifdef WF_TIMING
    std::cout << "Time to evolve " << diff.count() << " s\n";
#endif
    //Checkpoint gaugefield for solving adjoint equation in reverse
    if (step % (Nstep/Nckpoints) == 0 && step < Nstep) {
        ckpoint_fields[step/(Nstep/Nckpoints)] = out;
    }

    //Perform measurements
    for(auto const &meas : this->functions)
      if( step % meas.first == 0 ) meas.second(step,taus,out);
  }
}

template<class Gimpl, class GaugeAction, class SpinorField>
void FermionFlow<Gimpl, GaugeAction, SpinorField>::smear_adjoint(SpinorField& chi_out, const SpinorField& chi_in) const {
    std::cout << GridLogMessage
        << "[FermionFlow] Nstep   : " << Nstep << std::endl;
    std::cout << GridLogMessage
        << "[FermionFlow] epsilon   : " << epsilon << std::endl;
    std::cout << GridLogMessage
        << "[FermionFlow] full trajectory : " << Nstep * epsilon << std::endl;
    std::cout << GridLogMessage
        << "[FermionFlow] solving adjoint eq from t=" << Nstep*epsilon << " to 0" << std::endl;

    
    chi_out = chi_in;
    RealD taus = epsilon*Nstep;

//ck idx 0            1           2          17               18              19
//      [0] 1 2 3 4 [5] 6 7 8 9 [10] 11 ... [90] 91 92 93 94 [95] 96 97 98 99 [100]    <-- U(t+i*epsilon)
//                                           |------------>|
 

    for (int step = Nstep-1; step >= 0; step--) {
        
        std::cout << GridLogMessage << "[Adjoint FermionFlow] reverse step " << step << std::endl;
        int reconstruction_index = step%Nckpoints;
        int ckp_index = step/(Nstep/Nckpoints);


        GaugeField tmp_U(ckpoint_fields[ckp_index].Grid());
        //tmp_U = ckpoint_fields[ckp_index];

        if (reconstruction_index == Nstep/Nckpoints - 1) {
            //reconstruct and cache
            for (int i = 0; i < reconstruction_index; i++) {
                RealD inner_tau = 0;
                evolve_step(tmp_U,inner_tau);
                cached_fields[i]=tmp_U;
            }
            
        }
        else {
            tmp_U = cached_fields[reconstruction_index];
        }

        // for (int i = 0; i < reconstruction_index; i++) {
        //     RealD inner_tau = 0;
        //     evolve_step(tmp_U,inner_tau);
        // }
        evolve_step_adjoint(tmp_U, chi_out, taus); //goes from taus to taus-epsilon
        std::cout << GridLogMessage << "[Adjoint FermionFlow]  chi norm after (reverse) step = " << norm2(PeekIndex<SpinorIndex>(chi_out,0)) << std::endl;
    
    }

}


NAMESPACE_END(Grid);