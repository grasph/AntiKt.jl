// Based on code from Fastjet, Copyright (c) 2005-2020, Matteo Cacciari,
// Gavin P. Salam and Gregory Soyez (hep-ph/0512210 and in the manual,
// arXiv:1111.6097) 
//----------------------------------------------------------------------
//
//  FastJet is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  The algorithms that underlie FastJet have required considerable
//  development. They are described in the original FastJet paper,
//  . If you use
//  FastJet as part of work towards a scientific publication, please
//  quote the version you use and include a citation to the manual and
//  optionally also to hep-ph/0512210.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------


#include "PseudoJet.h"
#include<iostream>
#include<sstream>
#include<cmath>
#include<algorithm>
#include <cstdarg>

using namespace std;


//----------------------------------------------------------------------
// another constructor...
PseudoJet::PseudoJet(const double px_in, const double py_in, const double pz_in,
		     const double E_in) {

  _E  = E_in ;
  _px = px_in;
  _py = py_in;
  _pz = pz_in;

  this->_finish_init();

  // some default values for the history and user indices
  _reset_indices();

}


//----------------------------------------------------------------------
/// do standard end of initialisation
void PseudoJet::_finish_init () {
  _kt2 = this->px()*this->px() + this->py()*this->py();
  _phi = pseudojet_invalid_phi;
  // strictly speaking, _rap does not need initialising, because
  // it's never used as long as _phi == pseudojet_invalid_phi
  // (and gets set when _phi is requested). However ATLAS
  // 2013-03-28 complained that they sometimes have NaN's in
  // _rap and this interferes with some of their internal validation.
  // So we initialise it; penalty is about 0.3ns per PseudoJet out of about
  // 10ns total initialisation time (on a intel Core i7 2.7GHz)
  _rap = pseudojet_invalid_rap;
}

//----------------------------------------------------------------------
void PseudoJet::_set_rap_phi() const {

  if (_kt2 == 0.0) {
    _phi = 0.0; }
  else {
    _phi = atan2(this->py(),this->px());
  }
  if (_phi < 0.0) {_phi += twopi;}
  if (_phi >= twopi) {_phi -= twopi;} // can happen if phi=-|eps<1e-15|?
  if (this->E() == abs(this->pz()) && _kt2 == 0) {
    // Point has infinite rapidity -- convert that into a very large
    // number, but in such a way that different 0-pt momenta will have
    // different rapidities (so as to lift the degeneracy between
    // them) [this can be relevant at parton-level]
    double MaxRapHere = MaxRap + abs(this->pz());
    if (this->pz() >= 0.0) {_rap = MaxRapHere;} else {_rap = -MaxRapHere;}
  } else {
    // get the rapidity in a way that's modestly insensitive to roundoff
    // error when things pz,E are large (actually the best we can do without
    // explicit knowledge of mass)
    double effective_m2 = max(0.0,m2()); // force non tachyonic mass
    //    double E_plus_pz    = _E + abs(_pz); // the safer of p+, p-
    double E_plus_pz    = _E + fabs(_pz); // the safer of p+, p-
    // p+/p- = (p+ p-) / (p-)^2 = (kt^2+m^2)/(p-)^2
    _rap = 0.5*log((_kt2 + effective_m2)/(E_plus_pz*E_plus_pz));
    if (_pz > 0) {_rap = - _rap;}
  }

}


//----------------------------------------------------------------------
// return the pseudorapidity
double PseudoJet::eta() const {
  if (px() == 0.0 && py() ==0.0) return MaxRap;
  if (pz() == 0.0) return 0.0;

  double theta = atan(pt()/pz());
  if (theta < 0) theta += pi;
  return -log(tan(theta/2));
}

#if 0
//----------------------------------------------------------------------
// return "sum" of two pseudojets
PseudoJet operator+ (const PseudoJet & jet1, const PseudoJet & jet2) {
  //return PseudoJet(jet1.four_mom()+jet2.four_mom());
  return PseudoJet(jet1.px()+jet2.px(),
		   jet1.py()+jet2.py(),
		   jet1.pz()+jet2.pz(),
		   jet1.E() +jet2.E()  );
}


//----------------------------------------------------------------------
// return difference of two pseudojets
PseudoJet operator- (const PseudoJet & jet1, const PseudoJet & jet2) {
  //return PseudoJet(jet1.four_mom()-jet2.four_mom());
  return PseudoJet(jet1.px()-jet2.px(),
		   jet1.py()-jet2.py(),
		   jet1.pz()-jet2.pz(),
		   jet1.E() -jet2.E()  );
}


//----------------------------------------------------------------------
// return the product, coeff * jet
PseudoJet operator* (double coeff, const PseudoJet & jet) {
  // see the comment in operator*= about ensuring valid rap phi
  // before a multiplication to handle case of multiplication by
  // zero, while maintaining rapidity and phi
  jet._ensure_valid_rap_phi();
  //return PseudoJet(coeff*jet.four_mom());
  // the following code is hopefully more efficient
  PseudoJet coeff_times_jet(jet);
  coeff_times_jet *= coeff;
  return coeff_times_jet;
}

//----------------------------------------------------------------------
// return the product, coeff * jet
PseudoJet operator* (const PseudoJet & jet, double coeff) {
  return coeff*jet;
}

//----------------------------------------------------------------------
// return the ratio, jet / coeff
PseudoJet operator/ (const PseudoJet & jet, double coeff) {
  return (1.0/coeff)*jet;
}

//----------------------------------------------------------------------
/// multiply the jet's momentum by the coefficient
PseudoJet & PseudoJet::operator*=(double coeff) {
  // operator*= aims to maintain the rapidity and azimuth
  // for the PseudoJet; if they have already been evaluated
  // this is fine, but if they haven't and coeff is sufficiently
  // small as to cause a zero or underflow result, then a subsequent
  // invocation of rap or phi will lead to a non-sensical result.
  // So, here, we preemptively ensure that rapidity and phi
  // are correctly cached
  _ensure_valid_rap_phi();
  _px *= coeff;
  _py *= coeff;
  _pz *= coeff;
  _E  *= coeff;
  _kt2*= coeff*coeff;
  // phi and rap are unchanged
  return *this;
}

//----------------------------------------------------------------------
/// divide the jet's momentum by the coefficient
PseudoJet & PseudoJet::operator/=(double coeff) {
  (*this) *= 1.0/coeff;
  return *this;
}


//----------------------------------------------------------------------
/// add the other jet's momentum to this jet
PseudoJet & PseudoJet::operator+=(const PseudoJet & other_jet) {
  _px += other_jet._px;
  _py += other_jet._py;
  _pz += other_jet._pz;
  _E  += other_jet._E ;
  _finish_init(); // we need to recalculate phi,rap,kt2
  return *this;
}


//----------------------------------------------------------------------
/// subtract the other jet's momentum from this jet
PseudoJet & PseudoJet::operator-=(const PseudoJet & other_jet) {
  _px -= other_jet._px;
  _py -= other_jet._py;
  _pz -= other_jet._pz;
  _E  -= other_jet._E ;
  _finish_init(); // we need to recalculate phi,rap,kt2
  return *this;
}

//----------------------------------------------------------------------
bool operator==(const PseudoJet & a, const PseudoJet & b) {
  if (a.px() != b.px()) return false;
  if (a.py() != b.py()) return false;
  if (a.pz() != b.pz()) return false;
  if (a.E () != b.E ()) return false;

  if (a.user_index()    != b.user_index()) return false;
  if (a.cluster_hist_index() != b.cluster_hist_index()) return false;
  if (a.user_info_ptr() != b.user_info_ptr()) return false;
  if (a.structure_ptr() != b.structure_ptr()) return false;

  return true;
}

//----------------------------------------------------------------------
// check if the jet has zero momentum
bool operator==(const PseudoJet & jet, const double val) {
  if (val != 0)
    throw Error("comparing a PseudoJet with a non-zero constant (double) is not allowed.");
  return (jet.px() == 0 && jet.py() == 0 &&
	  jet.pz() == 0 && jet.E() == 0);
}

#endif


#ifdef EXTENDED_PSEUDOJET
//----------------------------------------------------------------------
void PseudoJet::reset_PtYPhiM(double pt_in, double y_in, double phi_in, double m_in) {
  assert(phi_in < 2*twopi && phi_in > -twopi);
  double ptm = (m_in == 0) ? pt_in : sqrt(pt_in*pt_in+m_in*m_in);
  double exprap = exp(y_in);
  double pminus = ptm/exprap;
  double pplus  = ptm*exprap;
  double px_local = pt_in*cos(phi_in);
  double py_local = pt_in*sin(phi_in);
  reset_momentum(px_local,py_local,0.5*(pplus-pminus),0.5*(pplus+pminus));
  set_cached_rap_phi(y_in,phi_in);
}

//----------------------------------------------------------------------
/// returns other.phi() - this.phi(), i.e. the phi distance to
/// other, constrained to be in range -pi .. pi
double PseudoJet::delta_phi_to(const PseudoJet & other) const {
  double dphi = other.phi() - phi();
  if (dphi >  pi) dphi -= twopi;
  if (dphi < -pi) dphi += twopi;
  return dphi;
}
#endif
