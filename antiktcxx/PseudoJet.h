//FJSTARTHEADER
// $Id: PseudoJet.hh 4442 2020-05-05 07:50:11Z soyez $
//
// Based on PseudoJet from Matteo Cacciari, Gavin P. Salam and Gregory Soyez
// (hep-ph/0512210, arXiv:1111.6097)
//----------------------------------------------------------------------
//
//  FastJet is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  The algorithms that underlie FastJet have required considerable
//  development. They are described in the original FastJet paper,
//  hep-ph/0512210 and in the manual, . If you use
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
//FJENDHEADER


#ifndef __PSEUDOJET_HH__
#define __PSEUDOJET_HH__

#include<vector>
#include<cassert>
#include<cmath>
#include<iostream>

/// Used to protect against parton-level events where pt can be zero
/// for some partons, giving rapidity=infinity. KtJet fails in those cases.
const double MaxRap = 1e5;

/// default value for phi, meaning it (and rapidity) have yet to be calculated)
const double pseudojet_invalid_phi = -100.0;
const double pseudojet_invalid_rap = -1e200;

static const double pi = 4.*atan(1.);
static const double twopi = 2*pi;

/// @ingroup basic_classes
/// \class PseudoJet
/// Class to contain pseudojets, including minimal information of use to
/// jet-clustering routines.
class PseudoJet {

 public:
  //----------------------------------------------------------------------
  /// @name Constructors and destructor
  //\{
  /// default constructor, which as of FJ3.0 provides an object for
  /// which all operations are now valid and which has zero momentum
  ///
  // (cf. this is actually OK from a timing point of view and in some
  // cases better than just having the default constructor for the
  // internal shared pointer: see PJtiming.cc and the notes therein)
  //PseudoJet() : _px(0), _py(0), _pz(0), _E(0) {_finish_init(); _reset_indices();}
  //PseudoJet() : _px(0), _py(0), _pz(0), _E(0), _phi(pseudojet_invalid_phi), _rap(pseudojet_invalid_rap), _kt2(0) {_reset_indices();}
  /// construct a pseudojet from explicit components
  PseudoJet(const double px, const double py, const double pz, const double E);

  /// default (virtual) destructor
  virtual ~PseudoJet(){};
  //\} ---- end of constructors and destructors --------------------------

  //----------------------------------------------------------------------
  /// @name Kinematic access functions
  //\{
  //----------------------------------------------------------------------
  inline double E()   const {return _E;}
  inline double px()  const {return _px;}
  inline double py()  const {return _py;}
  inline double pz()  const {return _pz;}

  /// returns phi (in the range 0..2pi)
  inline double phi() const {return phi_02pi();}
  
  /// returns phi in the range -pi..pi
  inline double phi_std()  const {
    _ensure_valid_rap_phi();
    return _phi > pi ? _phi-twopi : _phi;}

  /// returns phi in the range 0..2pi
  inline double phi_02pi() const {
    _ensure_valid_rap_phi();
    return _phi;
  }

  /// returns the rapidity or some large value when the rapidity
  /// is infinite
  inline double rap() const {
    _ensure_valid_rap_phi();
    return _rap;
  }

  double eta() const;


  /// returns the squared transverse momentum
  inline double pt2() const {return _kt2;}

  /// returns the scalar transverse momentum
  inline double  pt() const {return sqrt(_kt2);}

  /// returns the squared invariant mass // like CLHEP
  inline double  m2() const {return (_E+_pz)*(_E-_pz)-_kt2;}
  /// returns the invariant mass
  /// (If m2() is negative then -sqrt(-m2()) is returned, as in CLHEP)
  inline double  m() const;

#ifdef EXTENDED_PSEUDOJET
  /// returns the squared transverse mass = kt^2+m^2
  inline double mt2() const {return (_E+_pz)*(_E-_pz);}
  /// returns the transverse mass = sqrt(kt^2+m^2)
  inline double mt() const {return sqrt(std::abs(mperp2()));}

  /// return the transverse energy
  inline double Et() const {return (_kt2==0) ? 0.0 : _E/sqrt(1.0+_pz*_pz/_kt2);}

  /// return the transverse energy squared
  inline double Et2() const {return (_kt2==0) ? 0.0 : _E*_E/(1.0+_pz*_pz/_kt2);}

  /// return the cylinder (rap-phi) distance between this jet and another,
  /// \f$\Delta_R = \sqrt{\Delta y^2 + \Delta \phi^2}\f$.
  inline double delta_R(const PseudoJet & other) const {
    return sqrt(squared_distance(other));
  }

  /// returns other.phi() - this.phi(), constrained to be in
  /// range -pi .. pi
  double delta_phi_to(const PseudoJet & other) const;
#endif //EXTENDED_PSEUDOJET

  //\}  ------- end of kinematic access functions


  //----------------------------------------------------------------------
  /// @name Kinematic modification functions
  //\{
  //----------------------------------------------------------------------

 #ifdef EXTENDED_PSEUDOJET
  PseudoJet & operator*=(double);
  PseudoJet & operator/=(double);
  PseudoJet & operator+=(const PseudoJet &);
  PseudoJet & operator-=(const PseudoJet &);

  /// reset the 4-momentum according to the supplied components and
  /// put the user and history indices back to their default values
  inline void reset(double px, double py, double pz, double E);

  /// reset the PseudoJet according to the specified pt, rapidity,
  /// azimuth and mass (also resetting indices, etc.)
  /// (phi should satisfy -2pi<phi<4pi)
  inline void reset_PtYPhiM(double pt_in, double y_in, double phi_in, double m_in=0.0) {
    reset_momentum_PtYPhiM(pt_in, y_in, phi_in, m_in);
    _reset_indices();
  }
#endif

  //\} --- end of kin mod functions ------------------------------------



  //----------------------------------------------------------------------
  /// @name Members mainly intended for internal use
  //----------------------------------------------------------------------
  //\{
  /// return the cluster_hist_index, intended to be used by clustering
  /// routines.
  inline int cluster_hist_index() const {return _cluster_hist_index;}
  /// set the cluster_hist_index, intended to be used by clustering routines.
  inline void set_cluster_hist_index(const int index) {_cluster_hist_index = index;}

  //\} ---- end of internal use functions ---------------------------

 protected:

 private:
  // NB: following order must be kept for things to behave sensibly...
  double _px,_py,_pz,_E;
  mutable double _phi, _rap;
  double _kt2;
  int    _cluster_hist_index, _user_index;

  /// calculate phi, rap, kt2 based on the 4-momentum components
  void _finish_init();

  /// set the indices to default values
  void _reset_indices();

  /// ensure that the internal values for rapidity and phi
  /// correspond to 4-momentum structure
  inline void _ensure_valid_rap_phi() const {
    if (_phi == pseudojet_invalid_phi) _set_rap_phi();
  }

  /// set cached rapidity and phi values
  void _set_rap_phi() const;

};


//----------------------------------------------------------------------
// routines for basic binary operations

#if 0
PseudoJet operator+(const PseudoJet &, const PseudoJet &);
PseudoJet operator-(const PseudoJet &, const PseudoJet &);
PseudoJet operator*(double, const PseudoJet &);
PseudoJet operator*(const PseudoJet &, double);
PseudoJet operator/(const PseudoJet &, double);


/// returns true if the 4 momentum components of the two PseudoJets
/// are identical and all the internal indices (user, cluster_history)
/// + structure and user-info shared pointers are too
bool operator==(const PseudoJet &, const PseudoJet &);

/// inequality test which is exact opposite of operator==
inline bool operator!=(const PseudoJet & a, const PseudoJet & b) {return !(a==b);}

/// Can only be used with val=0 and tests whether all four
/// momentum components are equal to val (=0.0)
bool operator==(const PseudoJet & jet, const double val);
inline bool operator==(const double val, const PseudoJet & jet) {return jet == val;}

/// Can only be used with val=0 and tests whether at least one of the
/// four momentum components is different from val (=0.0)
inline bool operator!=(const PseudoJet & a, const double val)  {return !(a==val);}
inline bool operator!=( const double val, const PseudoJet & a) {return !(a==val);}

/// returns the 4-vector dot product of a and b
inline double dot_product(const PseudoJet & a, const PseudoJet & b) {
  return a.E()*b.E() - a.px()*b.px() - a.py()*b.py() - a.pz()*b.pz();
}

/// returns the cosine of the angle between a and b
inline double cos_theta(const PseudoJet & a, const PseudoJet & b) {
  double dot_3d = a.px()*b.px() + a.py()*b.py() + a.pz()*b.pz();
  return std::min(1.0, std::max(-1.0, dot_3d/sqrt(a.modp2()*b.modp2())));
}

/// returns the angle between a and b
inline double theta(const PseudoJet & a, const PseudoJet & b) {
  return acos(cos_theta(a,b));
}
#endif

//----------------------------------------------------------------------
inline void PseudoJet::_reset_indices() {
  set_cluster_hist_index(-1);
}


// taken literally from CLHEP
inline double PseudoJet::m() const {
  double mm = m2();
  return mm < 0.0 ? -std::sqrt(-mm) : std::sqrt(mm);
}

#ifdef EXTENDED_PSEUDOJET
inline void PseudoJet::reset(double px_in, double py_in, double pz_in, double E_in) {
  _px = px_in;
  _py = py_in;
  _pz = pz_in;
  _E  = E_in;
  _finish_init();
  _reset_indices();
}
#endif

#endif // __PSEUDOJET_HH__
