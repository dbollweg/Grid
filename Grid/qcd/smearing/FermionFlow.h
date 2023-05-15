#pragma once

NAMESPACE_BEGIN(Grid);

template<class Gimpl, class GaugeAction, class SpinorField>
class FermionFlow: public GradientFlowBase<Gimpl,GaugeAction> {
private:
    int Nstep;
    RealD epsilon;

    void evolve_step(typename Gimpl::GaugeField &U, SpinorField &chi, RealD &tau) const;
    void fermion_laplacian(typename Gimpl::GaugeField &U, SpinorField &chi_in, SpinorField &chi_out) const;

public:
    INHERIT_GIMPL_TYPES(Gimpl)

    FermionFlow(const RealD Epsilon, const int Nstep, unsigned int meas_interval = 1): GradientFlowBase<Gimpl,GaugeAction>(meas_interval), Nstep(Nstep), epsilon(Epsilon) {}

    void smear(GaugeField& out, SpinorField& chi_out, const GaugeField& in, const SpinorField& chi_in) const;
    void smear(GaugeField& out, const GaugeField& in) const override {};
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

template<class Gimpl, class GaugeAction, class SpinorField>
void FermionFlow<Gimpl, GaugeAction, SpinorField>::smear(GaugeField& out, SpinorField& chi_out, const GaugeField& in, const SpinorField& chi_in) const {
    std::cout << GridLogMessage
        << "[FermionFlow] Nstep   : " << Nstep<< std::endl;
    std::cout << GridLogMessage
        << "[FermionFlow] epsilon   : " << epsilon << std::endl;
    std::cout << GridLogMessage
        << "[FermionFlow] full trajectory : " << Nstep * epsilon << std::endl;

    SpinorField tmp(chi_out.Grid());
    out = in;
    chi_out = chi_in;

    RealD taus = 0;
    for (unsigned int step = 1; step <= Nstep; step++) { //step indicates the number of smearing steps applied at the time of measurement
    auto start = std::chrono::high_resolution_clock::now();
    evolve_step(out, chi_out, taus);
    tmp=chi_out;
    std::cout << GridLogMessage << "[FermionFlow]  chi norm after step = " << norm2(PeekIndex<SpinorIndex>(tmp,0)) << std::endl;
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
#ifdef WF_TIMING
    std::cout << "Time to evolve " << diff.count() << " s\n";
#endif
    //Perform measurements
    for(auto const &meas : this->functions)
      if( step % meas.first == 0 ) meas.second(step,taus,out);
  }
}



NAMESPACE_END(Grid);