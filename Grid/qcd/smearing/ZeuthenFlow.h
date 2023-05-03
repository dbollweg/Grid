#pragma once

NAMESPACE_BEGIN(Grid);

template<class Gimpl>
class ZeuthenFlowBase: public Smear<Gimpl>{
public:
    typedef std::function<void(int, RealD, const typename Gimpl::GaugeField &)> FunctionType;

protected:
    std::vector< std::pair<int, FunctionType> > functions;

    mutable SymanzikGaugeAction<Gimpl> SG;

public:
    INHERIT_GIMPL_TYPES(Gimpl)

    explicit ZeuthenFlowBase(unsigned int meas_interval = 1):
        SG(SymanzikGaugeAction<Gimpl>(3.0)) {

            setDefaultMeasurements(meas_interval);
        }

    void resetActions() { functions.clear(); }

    void addMeasurement(int meas_interval, FunctionType meas) {
        functions.push_back({meas_interval, meas});
    }

    void setDefaultMeasurements(int topq_meas_interval = 1);

    void derivative(GaugeField&, const GaugeField&, const GaugeField &) const override{
        assert(0);
    }

    static RealD energyDensityPlaquette(const RealD t, const GaugeField& U);

    std::vector<RealD> flowMeasureEnergyDensityPlaquette(GaugeField& V, const GaugeField& U, int measure_interval = 1);

    std::vector<RealD> flowMeasureEnergyDensityPlaquette(const GaugeField& U, int measure_interval = 1);

};


template <class Gimpl>
class ZeuthenFlow: public ZeuthenFlowBase<Gimpl>{
private:
    int Nstep;
    RealD epsilon;

    void evolve_step(typename Gimpl::GaugeField &U, RealD &tau) const;

    void zeuthen_force(typename Gimpl::GaugeField &U, typename Gimpl::GaugeField &Z) const;

public:
    INHERIT_GIMPL_TYPES(Gimpl)

    ZeuthenFlow(const RealD epsilon, const int Nstep, unsigned int meas_interval = 1): ZeuthenFlowBase<Gimpl>(meas_interval), Nstep(Nstep), epsilon(epsilon) {}

    void smear(GaugeField& out, const GaugeField& in) const override;
};

template <class Gimpl>
RealD ZeuthenFlowBase<Gimpl>::energyDensityPlaquette(const RealD t, const GaugeField& U){
  static WilsonGaugeAction<Gimpl> SG(3.0);
  return 2.0 * t * t * SG.S(U)/U.Grid()->gSites();
}

template <class Gimpl>
std::vector<RealD> ZeuthenFlowBase<Gimpl>::flowMeasureEnergyDensityPlaquette(GaugeField &V, const GaugeField& U, int measure_interval){
  std::vector<RealD> out;
  resetActions();
  addMeasurement(measure_interval, [&out](int step, RealD t, const typename Gimpl::GaugeField &U){ 
      std::cout << GridLogMessage << "[ZeuthenFlow] Computing plaquette energy density for step " << step << std::endl;
      out.push_back( energyDensityPlaquette(t,U) );
    });      
  smear(V,U);
  return out;
}

template <class Gimpl>
std::vector<RealD> ZeuthenFlowBase<Gimpl>::flowMeasureEnergyDensityPlaquette(const GaugeField& U, int measure_interval){
  GaugeField V(U);
  return flowMeasureEnergyDensityPlaquette(V,U, measure_interval);
}

template <class Gimpl>
void ZeuthenFlowBase<Gimpl>::setDefaultMeasurements(int topq_meas_interval){
  addMeasurement(1, [](int step, RealD t, const typename Gimpl::GaugeField &U){
      std::cout << GridLogMessage << "[ZeuthenFlow] Energy density (plaq) : "  << step << "  " << t << "  " << energyDensityPlaquette(t,U) << std::endl;
    });
  addMeasurement(topq_meas_interval, [](int step, RealD t, const typename Gimpl::GaugeField &U){
      std::cout << GridLogMessage << "[ZeuthenFlow] Top. charge           : "  << step << "  " << WilsonLoops<Gimpl>::TopologicalCharge(U) << std::endl;
    });
}

template<class Gimpl>
void ZeuthenFlow<Gimpl>::zeuthen_force(typename Gimpl::GaugeField &U, typename Gimpl::GaugeField &Z) const{
    GaugeField Zprime(U.Grid());
    for (int mu = 0; mu < Nd; mu++) {
        GaugeLinkField Ulink(U.Grid());
        GaugeLinkField Zp(U.Grid());
        Ulink = PeekIndex<LorentzIndex>(U,mu);
        //U(x)*F(x+mu)*U^{dagger}(x) + U^{dagger}(x-mu)F(x-mu)U(x-mu)
        Zp = Ulink * Gimpl::CshiftLink(PeekIndex<LorentzIndex>(Z,mu), mu, 1) * adj(Ulink) + adj(Gimpl::CshiftLink(Ulink,mu,-1)) * Gimpl::CshiftLink(PeekIndex<LorentzIndex>(Z,mu),mu,-1) * Gimpl::CshiftLink(Ulink,mu,-1);        

        PokeIndex<LorentzIndex>(Zprime, Zp, mu);
    }
    Z = (5.0/6.0)*Z + (1.0/12.0)*Zprime;

}

template<class Gimpl>
void ZeuthenFlow<Gimpl>::evolve_step(typename Gimpl::GaugeField &U, RealD &tau) const{
    GaugeField Z(U.Grid());
    GaugeField tmp(U.Grid());
    this->SG.deriv(U, Z);

    //in contrast to Wilson flow, we now need to apply the cov. derivatives
    //from the zeuthen flow equation!

    zeuthen_force(U, Z);

    Z*= 0.25;
    Gimpl::update_field(Z, U, -2.0*epsilon);

    Z *= -17.0/8.0;
    this->SG.deriv(U, tmp);
    zeuthen_force(U, tmp);
    Z += tmp;

    Z *= 8.0/9.0;
    Gimpl::update_field(Z, U, -2.0*epsilon);

    Z *= -4.0/3.0;
    this->SG.deriv(U, tmp);
    zeuthen_force(U,tmp);
    Z += tmp;

    Z *= 3.0/4.0;

    Gimpl::update_field(Z, U, -2.0*epsilon);

    tau+=epsilon;
}



template <class Gimpl>
void ZeuthenFlow<Gimpl>::smear(GaugeField& out, const GaugeField& in) const{
  std::cout << GridLogMessage
	    << "[ZeuthenFlow] Nstep   : " << Nstep << std::endl;
  std::cout << GridLogMessage
	    << "[ZeuthenFlow] epsilon : " << epsilon << std::endl;
  std::cout << GridLogMessage
	    << "[ZeuthenFlow] full trajectory : " << Nstep * epsilon << std::endl;

  out = in;
  RealD taus = 0.;
  for (unsigned int step = 1; step <= Nstep; step++) { //step indicates the number of smearing steps applied at the time of measurement
    auto start = std::chrono::high_resolution_clock::now();
    evolve_step(out, taus);
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
