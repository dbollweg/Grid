/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/modules/plaquette.h

Copyright (C) 2017

Author: Guido Cossu <guido.cossu@ed.ac.uk>
Author: Christopher Kelly <ckelly@bnl.gov>
Author: Dennis Bollweg <dbollweg@bnl.gov>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
			   /*  END LEGAL */

#pragma once

NAMESPACE_BEGIN(Grid);

template <class Gimpl, class GaugeAction>
class GradientFlowBase: public Smear<Gimpl>{
public: 

    typedef std::function<void(int, RealD, const typename Gimpl::GaugeField &)> FunctionType;

protected:

    std::vector<std::pair<int, FunctionType> > functions;

    mutable GaugeAction SG;

public:

    INHERIT_GIMPL_TYPES(Gimpl)

    explicit GradientFlowBase(unsigned int meas_interval = 1): SG(3.0) {
            setDefaultMeasurements(meas_interval);
        }

    void resetActions() { functions.clear(); }

    void addMeasurement(int meas_interval, FunctionType meas) { functions.push_back({meas_interval, meas}); }

    void setDefaultMeasurements(int topq_meas_interval = 1);

    void derivative(GaugeField&, const GaugeField&, const GaugeField&) const override {
        assert(0);
    }


    //Compute t^2 <E(t)> for time t from the plaquette
    static RealD energyDensityPlaquette(const RealD t, const GaugeField& U);

    //Compute t^2 <E(t)> for time t from the 1x1 cloverleaf form
    //t is the Wilson flow time
    static RealD energyDensityCloverleaf(const RealD t, const GaugeField& U);
    
    //Evolve the gauge field by Nstep steps of epsilon and return the energy density computed every interval steps
    //The smeared field is output as V
    std::vector<RealD> flowMeasureEnergyDensityPlaquette(GaugeField &V, const GaugeField& U, int measure_interval = 1);

    //Version that does not return the smeared field
    std::vector<RealD> flowMeasureEnergyDensityPlaquette(const GaugeField& U, int measure_interval = 1);


    //Evolve the gauge field by Nstep steps of epsilon and return the Cloverleaf energy density computed every interval steps
    //The smeared field is output as V
    std::vector<RealD> flowMeasureEnergyDensityCloverleaf(GaugeField &V, const GaugeField& U, int measure_interval = 1);

    //Version that does not return the smeared field
    std::vector<RealD> flowMeasureEnergyDensityCloverleaf(const GaugeField& U, int measure_interval = 1);
};

template <class Gimpl, class GaugeAction>
class GradientFlow: public GradientFlowBase<Gimpl, GaugeAction> {
private:
    int Nstep;
    RealD epsilon;

    void evolve_step(typename Gimpl::GaugeField &U, RealD &tau) const;

public:
    INHERIT_GIMPL_TYPES(Gimpl)

    GradientFlow(const RealD epsilon, const int Nstep, unsigned int meas_interval = 1): 
        GradientFlowBase<Gimpl, GaugeAction>(meas_interval), Nstep(Nstep), epsilon(epsilon){}

    void smear(GaugeField& out, const GaugeField& in) const override;

};

template <class Gimpl, class GaugeAction>
class GradientFlowAdaptive: public GradientFlowBase<Gimpl, GaugeAction> {
private:
    RealD init_epsilon;
    RealD maxTau;
    RealD tolerance;

    int evolve_step_adaptive(typename Gimpl::GaugeField& U, RealD &tau, RealD &eps) const;
public:
    INHERIT_GIMPL_TYPES(Gimpl);

    GradientFlowAdaptive(const RealD init_epsilon, const RealD maxTau, const RealD tolerance, unsigned int meas_interval = 1):
        GradientFlowBase<Gimpl,GaugeAction>(meas_interval), init_epsilon(init_epsilon), maxTau(maxTau), tolerance(tolerance){}

    void smear(GaugeField& out, const GaugeField& in) const override;    
};

template <class Gimpl, class GaugeAction>
RealD GradientFlowBase<Gimpl, GaugeAction>::energyDensityPlaquette(const RealD t, const GaugeField& U){
  static GaugeAction SG(3.0);
  return 2.0 * t * t * SG.S(U)/U.Grid()->gSites();
}

//Compute t^2 <E(t)> for time from the 1x1 cloverleaf form
template <class Gimpl, class GaugeAction>
RealD GradientFlowBase<Gimpl, GaugeAction>::energyDensityCloverleaf(const RealD t, const GaugeField& U){
  typedef typename Gimpl::GaugeLinkField GaugeMat;
  typedef typename Gimpl::GaugeField GaugeLorentz;

  assert(Nd == 4);
  //E = 1/2 tr( F_munu F_munu )
  //However as  F_numu = -F_munu, only need to sum the trace of the squares of the following 6 field strengths:
  //F_01 F_02 F_03   F_12 F_13  F_23
  GaugeMat F(U.Grid());
  LatticeComplexD R(U.Grid());
  R = Zero();
  
  for(int mu=0;mu<3;mu++){
    for(int nu=mu+1;nu<4;nu++){
      WilsonLoops<Gimpl>::FieldStrength(F, U, mu, nu);
      R = R + trace(F*F);
    }
  }
  ComplexD out = sum(R);
  out = t*t*out / RealD(U.Grid()->gSites());
  return -real(out); //minus sign necessary for +ve energy
}


template <class Gimpl, class GaugeAction>
std::vector<RealD> GradientFlowBase<Gimpl, GaugeAction>::flowMeasureEnergyDensityPlaquette(GaugeField &V, const GaugeField& U, int measure_interval){
  std::vector<RealD> out;
  resetActions();
  addMeasurement(measure_interval, [&out](int step, RealD t, const typename Gimpl::GaugeField &U){ 
      std::cout << GridLogMessage << "[GradientFlow] Computing plaquette energy density for step " << step << std::endl;
      out.push_back( energyDensityPlaquette(t,U) );
    });      
  smear(V,U);
  return out;
}

template <class Gimpl, class GaugeAction>
std::vector<RealD> GradientFlowBase<Gimpl, GaugeAction>::flowMeasureEnergyDensityPlaquette(const GaugeField& U, int measure_interval){
  GaugeField V(U);
  return flowMeasureEnergyDensityPlaquette(V,U, measure_interval);
}

template <class Gimpl, class GaugeAction>
std::vector<RealD> GradientFlowBase<Gimpl, GaugeAction>::flowMeasureEnergyDensityCloverleaf(GaugeField &V, const GaugeField& U, int measure_interval){
  std::vector<RealD> out;
  resetActions();
  addMeasurement(measure_interval, [&out](int step, RealD t, const typename Gimpl::GaugeField &U){ 
      std::cout << GridLogMessage << "[GradientFlow] Computing Cloverleaf energy density for step " << step << std::endl;
      out.push_back( energyDensityCloverleaf(t,U) );
    });      
  smear(V,U);
  return out;
}

template <class Gimpl, class GaugeAction>
std::vector<RealD> GradientFlowBase<Gimpl, GaugeAction>::flowMeasureEnergyDensityCloverleaf(const GaugeField& U, int measure_interval){
  GaugeField V(U);
  return flowMeasureEnergyDensityCloverleaf(V,U, measure_interval);
}

template <class Gimpl, class GaugeAction>
void GradientFlowBase<Gimpl, GaugeAction>::setDefaultMeasurements(int topq_meas_interval){
  addMeasurement(1, [](int step, RealD t, const typename Gimpl::GaugeField &U){
      std::cout << GridLogMessage << "[GradientFlow] Energy density (plaq) : "  << step << "  " << t << "  " << energyDensityPlaquette(t,U) << std::endl;
    });
  addMeasurement(topq_meas_interval, [](int step, RealD t, const typename Gimpl::GaugeField &U){
      std::cout << GridLogMessage << "[GradientFlow] Top. charge           : "  << step << "  " << WilsonLoops<Gimpl>::TopologicalCharge(U) << std::endl;
    });
}




template <class Gimpl, class GaugeAction>
void GradientFlow<Gimpl, GaugeAction>::evolve_step(typename Gimpl::GaugeField &U, RealD &tau) const{
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


template <class Gimpl, class GaugeAction>
void GradientFlow<Gimpl, GaugeAction>::smear(GaugeField& out, const GaugeField& in) const{
  std::cout << GridLogMessage
	    << "[GradientFlow] Nstep   : " << Nstep << std::endl;
  std::cout << GridLogMessage
	    << "[GradientFlow] epsilon : " << epsilon << std::endl;
  std::cout << GridLogMessage
	    << "[GradientFlow] full trajectory : " << Nstep * epsilon << std::endl;

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

template <class Gimpl, class GaugeAction>
int GradientFlowAdaptive<Gimpl, GaugeAction>::evolve_step_adaptive(typename Gimpl::GaugeField &U, RealD &tau, RealD &eps) const{
    if (maxTau - tau < eps){
        eps = maxTau-tau;
    }
    //std::cout << GridLogMessage << "Integration epsilon : " << epsilon << std::endl;
    GaugeField Z(U.Grid());
    GaugeField Zprime(U.Grid());
    GaugeField tmp(U.Grid()), Uprime(U.Grid()), Usave(U.Grid());
    Uprime = U;
    Usave = U;

    this->SG.deriv(U, Z);
    Zprime = -Z;
    Z *= 0.25;                                  // Z0 = 1/4 * F(U)
    Gimpl::update_field(Z, U, -2.0*eps);    // U = W1 = exp(ep*Z0)*W0

    Z *= -17.0/8.0;
    this->SG.deriv(U, tmp); Z += tmp;                 // -17/32*Z0 +Z1
    Zprime += 2.0*tmp;
    Z *= 8.0/9.0;                               // Z = -17/36*Z0 +8/9*Z1
    Gimpl::update_field(Z, U, -2.0*eps);    // U_= W2 = exp(ep*Z)*W1
        

    Z *= -4.0/3.0;
    this->SG.deriv(U, tmp); Z += tmp;                 // 4/3*(17/36*Z0 -8/9*Z1) +Z2
    Z *= 3.0/4.0;                               // Z = 17/36*Z0 -8/9*Z1 +3/4*Z2
    Gimpl::update_field(Z, U, -2.0*eps);    // V(t+e) = exp(ep*Z)*W2

    // Ramos arXiv:1301.4388
    Gimpl::update_field(Zprime, Uprime, -2.0*eps); // V'(t+e) = exp(ep*Z')*W0

    // Compute distance using Ramos' definition
    GaugeField diffU = U - Uprime;
    RealD max_dist = 0;

    for(int mu=0;mu<Nd;mu++){
        typename Gimpl::GaugeLinkField diffU_mu = PeekIndex<LorentzIndex>(diffU, mu);
        RealD dist_mu = sqrt( maxLocalNorm2(diffU_mu) ) /Nc/Nc; //maximize over sites
        max_dist = std::max(max_dist, dist_mu); //maximize over mu
    }
    
    int ret;
    if(max_dist < tolerance) {
        tau += eps;
        ret = 1;
    } else {
        U = Usave;
        ret = 0;
    }
    eps = eps*0.95*std::pow(tolerance/max_dist,1./3.);
    std::cout << GridLogMessage << "Adaptive smearing : Distance: "<< max_dist <<" Step successful: " << ret << " New epsilon: " << eps << std::endl; 

    return ret;
}



template <class Gimpl, class GaugeAction>
void GradientFlowAdaptive<Gimpl,GaugeAction>::smear(GaugeField& out, const GaugeField& in) const{
  std::cout << GridLogMessage
	    << "[GradientFlow] initial epsilon : " << init_epsilon << std::endl;
  std::cout << GridLogMessage
	    << "[GradientFlow] full trajectory : " << maxTau << std::endl;
  std::cout << GridLogMessage
	    << "[GradientFlow] tolerance   : " << tolerance << std::endl;
  out = in;
  RealD taus = 0.;
  RealD eps = init_epsilon;
  unsigned int step = 0;
  do{
    int step_success = evolve_step_adaptive(out, taus, eps); 
    step += step_success; //step will not be incremented if the integration step fails

    //Perform measurements
    if(step_success)
      for(auto const &meas : this->functions)
	if( step % meas.first == 0 ) meas.second(step,taus,out);
  } while (taus < maxTau);
}


//Wrapper for Zeuthen flow. Modifies deriv such that the (1+1/a^2 laplace) operation 
//is applied after dS/dU calculation. 
template<class Gimpl>
class ZeuthenAction : public SymanzikGaugeAction<Gimpl> {
public:
    INHERIT_GIMPL_TYPES(Gimpl);
    ZeuthenAction(RealD beta) : SymanzikGaugeAction<Gimpl>(beta) {};

    void zeuthen_force(const typename Gimpl::GaugeField &U, typename Gimpl::GaugeField &Z) const {
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

    virtual void deriv(const GaugeField& Umu, GaugeField& dSdU) {
        SymanzikGaugeAction<Gimpl>::deriv(Umu,dSdU);
        this->zeuthen_force(Umu,dSdU);
    }    
};

template<class Gimpl>
class WilsonFlow : public GradientFlow<Gimpl,WilsonGaugeAction<Gimpl> > {
    public:
    INHERIT_GIMPL_TYPES(Gimpl);
    WilsonFlow(const RealD epsilon, const int Nstep, unsigned int meas_interval = 1):
        GradientFlow<Gimpl,WilsonGaugeAction<Gimpl> >(epsilon, Nstep, meas_interval) {}
};

template<class Gimpl>
class ZeuthenFlow : public GradientFlow<Gimpl,ZeuthenAction<Gimpl> > {
    public:
    INHERIT_GIMPL_TYPES(Gimpl);
    ZeuthenFlow(const RealD epsilon, const int Nstep, unsigned int meas_interval = 1):
        GradientFlow<Gimpl,ZeuthenAction<Gimpl> >(epsilon, Nstep, meas_interval) {}
};

template<class Gimpl>
class WilsonFlowAdaptive : public GradientFlowAdaptive<Gimpl,WilsonGaugeAction<Gimpl> > {
    public:
    INHERIT_GIMPL_TYPES(Gimpl);
    WilsonFlowAdaptive(const RealD init_epsilon, const RealD maxTau, const RealD tolerance, unsigned int meas_interval = 1):
        GradientFlowAdaptive<Gimpl,WilsonGaugeAction<Gimpl> >(init_epsilon, maxTau, tolerance, meas_interval) {}
};

template<class Gimpl>
class ZeuthenFlowAdaptive : public GradientFlowAdaptive<Gimpl,ZeuthenAction<Gimpl> > {
    public:
    INHERIT_GIMPL_TYPES(Gimpl);
    ZeuthenFlowAdaptive(const RealD init_epsilon, const RealD maxTau, const RealD tolerance, unsigned int meas_interval = 1):
        GradientFlowAdaptive<Gimpl,ZeuthenAction<Gimpl> >(init_epsilon, maxTau, tolerance, meas_interval) {}
};


NAMESPACE_END(Grid);