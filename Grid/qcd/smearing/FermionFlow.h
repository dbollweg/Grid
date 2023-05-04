#pragma once

NAMESPACE_BEGIN(Grid);

template<class Gimpl, class SpinorField>
class FermionFlow: public ZeuthenFlowBase<Gimpl> {
private:
    int Nstep;
    RealD epsilon;

    void evolve_step(typename Gimpl::GaugeField &U, SpinorField chi, RealD &tau) const;
    void fermion_laplacian(typename Gimpl::GaugeField &U, SpinorField chi_in, SpinorField chi_out) const;

public:
    INHERIT_GIMPL_TYPES(Gimpl)

    FermionFlow(const RealD Epsilon, const int Nstep, unsigned int meas_interval = 1): ZeuthenFlowBase<Gimpl>(meas_interval), Nstep(Nstep), epsilon(Epsilon) {}

    void smear(GaugeField& out, SpinorField& chi_out, const GaugeField& in, const SpinorField& chi_in) const;
    void smear(GaugeField& out, const GaugeField& in) const override {};
};

template<class Gimpl, class SpinorField>
void FermionFlow<Gimpl, SpinorField>::fermion_laplacian(typename Gimpl::GaugeField &U, SpinorField chi_in, SpinorField chi_out) const {
    chi_out = 0;
    for (int mu = 0; mu < Nd; mu++) {
        GaugeLinkField Ulink(U.Grid());
        Ulink = PeekIndex<LorentzIndex>(U,mu);
        chi_out += Ulink*Cshift(chi_in,mu,1)-2.0*chi_in+adj(Gimpl::CshiftLink(Ulink,mu,-1))*Cshift(chi_in,mu,-1);
    }
}

template<class Gimpl, class SpinorField>
void FermionFlow<Gimpl,SpinorField>::evolve_step(typename Gimpl::GaugeField &U, SpinorField chi, RealD &tau) const {
    GaugeField Z(U.Grid());
    GaugeField tmp(U.Grid());
    SpinorField phi0(U.Grid());
    SpinorField phi1(U.Grid());
    SpinorField phi2(U.Grid());
    //W0
    phi0=chi;

    this->SG.deriv(U,Z);
    this->zeuthen_force(U,Z);
    fermion_laplacian(U,chi,phi0);
    Z*=0.25;
    Gimpl::update_field(Z, U, -2.0*epsilon);
    phi1 = chi + 0.25*epsilon*phi0;
    //U is now W1

    Z *= -17.0/8.0;
    this->SG.deriv(U, tmp);
    this->zeuthen_force(U, tmp);
    

    fermion_laplacian(U,phi1,phi2); //phi2 is laplace(W1)phi1
    //phi2:
    phi2 = chi + 8.0/9.0*epsilon*phi2-2.0/9.0*phi0;
    Z += tmp;

    Z *= 8.0/9.0;
    Gimpl::update_field(Z, U, -2.0*epsilon);
    // U is now W2

    Z *= -4.0/3.0;
    this->SG.deriv(U, tmp);
    this->zeuthen_force(U, tmp);
    Z += tmp;

    Z *= 3.0/4.0;

    fermion_laplacian(U,phi2,phi2);
    chi = phi1 + 3.0/4.0*epsilon*phi2;

    Gimpl::update_field(Z, U, -2.0*epsilon);


    tau+=epsilon;




}

template<class Gimpl, class SpinorField>
void FermionFlow<Gimpl,SpinorField>::smear(GaugeField& out, SpinorField& chi_out, const GaugeField& in, const SpinorField& chi_in) const {
    std::cout << GridLogMessage
        << "[FermionFlow + ZeuthenFlow] Nstep   : " << Nstep<< std::endl;
    std::cout << GridLogMessage
        << "[FermionFlow + ZeuthenFLow] epsilon   : " << epsilon << std::endl;
    std::cout << GridLogMessage
        << "[FermionFlow + ZeuthenFlow] full trajectory : " << Nstep * epsilon << std::endl;


    out = in;
    chi_out = chi_in;

    RealD taus = 0;
    for (unsigned int step = 1; step <= Nstep; step++) { //step indicates the number of smearing steps applied at the time of measurement
    auto start = std::chrono::high_resolution_clock::now();
    evolve_step(out, chi_out, taus);
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