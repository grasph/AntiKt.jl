// Code extended from fastjet ClusterSequence_Tiled_N2.cc
//
// Copyright (c) 2005-2020, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet.
//
//  FastJet is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  The algorithms that underlie FastJet have required considerable
//  development. They are described in the original FastJet paper,
//  hep-ph/0512210 and in the manual, arXiv:1111.6097. If you use
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


// The tiled N^2 part of the ClusterSequence class -- separated out
// from the rest of the class implementation so as to speed up
// compilation of this particular part while it is under test.

#include<iostream>
#include<vector>
#include<cmath>
#include<algorithm>
#include <tuple>
#include <math.h>
#include <stdexcept>
#include "PseudoJet.h"

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/ReaderAscii.h"

#include <sys/time.h> //for code timing
#include <unistd.h>//for code timing
#include <cstring>
#include <cstdio>

double R = 0.4;
double ptmin = 5.0;

using namespace std;

//static const double pi = 4*atan(1.);
//static const double twopi = 2*pi;


struct history_element{
  int parent1; /// index in _history where first parent of this jet
  /// was created (InexistentParent if this jet is an
  /// original particle)

  int parent2; /// index in _history where second parent of this jet
  /// was created (InexistentParent if this jet is an
  /// original particle); BeamJet if this history entry
  /// just labels the fact that the jet has recombined
  /// with the beam)

  int child;   /// index in _history where the current jet is
  /// recombined with another jet to form its child. It
  /// is Invalid if this jet does not further
  /// recombine.

  int jetp_index; /// index in the _jets vector where we will find the
  /// PseudoJet object corresponding to this jet
  /// (i.e. the jet created at this entry of the
  /// history). NB: if this element of the history
  /// corresponds to a beam recombination, then
  /// jetp_index=Invalid.

  double dij;  /// the distance corresponding to the recombination
  /// at this stage of the clustering.

  double max_dij_so_far; /// the largest recombination distance seen
  /// so far in the clustering history.
};

enum JetType {Invalid=-3, InexistentParent = -2, BeamJet = -1};

/// structure analogous to BriefJet, but with the extra information
/// needed for dealing with tiles
struct TiledJet {
  double     eta, phi, kt2, NN_dist;
  TiledJet * NN, *previous, * next;
  int        _jets_index, tile_index, diJ_posn;
};

//template <class J>
inline double _tj_dist(const TiledJet* const jetA,
		       const TiledJet* const jetB){
  double dphi = pi-std::abs(pi-std::abs(jetA->phi - jetB->phi));
  double deta = (jetA->eta - jetB->eta);
  return dphi*dphi + deta*deta;
}

//----------------------------------------------------------------------
//template <class J>
inline double _tj_diJ(const TiledJet* const jet){
  double kt2 = jet->kt2;
  if (jet->NN != NULL) {if (jet->NN->kt2 < kt2) {kt2 = jet->NN->kt2;}}
  return jet->NN_dist * kt2;
}



struct Tile {
  static const int n_tile_neighbours = 9;

  /// pointers to neighbouring tiles, including self
  Tile* begin_tiles[n_tile_neighbours];

  /// neighbouring tiles, excluding selfo
  //Tile** surrounding_tiles;

  /// half of neighbouring tiles, no self
  Tile** RH_tiles;

  /// just beyond end of tiles
  Tile** end_tiles;

  /// start of list of TiledJets contained in this tile
  TiledJet * head;

  /// sometimes useful to be able to tag a tile
  bool tagged;
};

struct Tiling{
  std::vector<Tile> _tiles;
  double _tiles_eta_min, _tiles_eta_max;
  double _tile_size_eta, _tile_size_phi;
  int    _n_tiles_phi,_tiles_ieta_min,_tiles_ieta_max;
};



struct ClusterSequence {
  /// This contains the physical PseudoJets; for each PseudoJet one
  /// can find the corresponding position in the _history by looking
  /// at _jets[i].cluster_hist_index().
  std::vector<PseudoJet> _jets;

  /// this vector will contain the branching history; for each stage,
  /// _history[i].jetp_index indicates where to look in the _jets
  /// vector to get the physical PseudoJet.
  std::vector<history_element> _history;

  // make a local copy of the jet definition (for future use)
  //_jet_def = jet_def_in; //needed ??

  //  _writeout_combinations = writeout_combinations;


  double _Qtot;

  Tiling tiling;

  // initialised the wrapper to the current CS
  //_structure_shared_ptr.reset(new ClusterSequenceStructure(this));

  //SharedPtr<PseudoJetStructureBase> _structure_shared_ptr;
};


// reasonably robust return of tile index given ieta and iphi, in particular
// it works even if iphi is negative
int _tile_index (const Tiling& tiling, int ieta, int iphi){
  // note that (-1)%n = -1 so that we have to add _n_tiles_phi
  // before performing modulo operation
  return (ieta-tiling._tiles_ieta_min)*tiling._n_tiles_phi
    + (iphi+tiling._n_tiles_phi) % tiling._n_tiles_phi;
}


//----------------------------------------------------------------------
/// return the tile index corresponding to the given eta,phi point
int _tile_index(const Tiling& tiling, const double eta, const double phi) {
  int ieta, iphi;
  if      (eta <= tiling._tiles_eta_min) {ieta = 0;}
  else if (eta >= tiling._tiles_eta_max) {ieta = tiling._tiles_ieta_max-tiling._tiles_ieta_min;}
  else {
    //ieta = int(floor((eta - _tiles_eta_min) / _tile_size_eta));
    ieta = int(((eta - tiling._tiles_eta_min) / tiling._tile_size_eta));
    // following needed in case of rare but nasty rounding errors
    if (ieta > tiling._tiles_ieta_max-tiling._tiles_ieta_min) {
      ieta = tiling._tiles_ieta_max-tiling._tiles_ieta_min;}
  }
  // allow for some extent of being beyond range in calculation of phi
  // as well
  //iphi = (int(floor(phi/_tile_size_phi)) + _n_tiles_phi) % _n_tiles_phi;
  // with just int and no floor, things run faster but beware
  iphi = int((phi+twopi)/tiling._tile_size_phi) % tiling._n_tiles_phi;
  return (iphi + ieta * tiling._n_tiles_phi);
}


std::tuple<double, double>
determine_rapidity_extent(const vector<PseudoJet> & particles) {
  // Have a binning of rapidity that goes from -nrap to nrap
  // in bins of size 1; the left and right-most bins include
  // include overflows from smaller/larger rapidities
  int nrap = 20;
  int nbins = 2*nrap;
  std::vector<double> counts(nbins, 0);

  // get the minimum and maximum rapidities and at the same time bin
  // the multiplicities as a function of rapidity to help decide how
  // far out it's worth going
  double minrap =  numeric_limits<double>::max();
  double maxrap = -numeric_limits<double>::max();
  int ibin;
  for (unsigned i = 0; i < particles.size(); i++) {
    // ignore particles with infinite rapidity
    if (particles[i].E() == abs(particles[i].pz())) continue;
    double rap = particles[i].rap();
    if (rap < minrap) minrap = rap;
    if (rap > maxrap) maxrap = rap;
    // now bin the rapidity to decide how far to go with the tiling.
    // Remember the bins go from ibin=0 (rap=-infinity..-19)
    // to ibin = nbins-1 (rap=19..infinity for nrap=20)
    ibin = int(rap+nrap);
    if (ibin < 0) ibin = 0;
    if (ibin >= nbins) ibin = nbins - 1;
    counts[ibin]++;
  }

  // now figure out the particle count in the busiest bin
  double max_in_bin = 0;
  for (ibin = 0; ibin < nbins; ibin++) {
    if (max_in_bin < counts[ibin]) max_in_bin = counts[ibin];
  }

  // and find minrap, maxrap such that edge bin never contains more
  // than some fraction of busiest, and at least a few particles; first do
  // it from left. NB: the thresholds chosen here are largely
  // guesstimates as to what might work.
  //
  // 2014-07-17: in some tests at high multiplicity (100k) and particles going up to
  //             about 7.3, anti-kt R=0.4, we found that 0.25 gave 20% better run times
  //             than the original value of 0.5.
  const double allowed_max_fraction = 0.25;
  // the edge bins should also contain at least min_multiplicity particles
  const double min_multiplicity = 4;
  // now calculate how much we can accumulate into an edge bin
  double allowed_max_cumul = floor(max(max_in_bin * allowed_max_fraction, min_multiplicity));
  // make sure we don't require more particles in a bin than max_in_bin
  if (allowed_max_cumul > max_in_bin) allowed_max_cumul = max_in_bin;

  // start scan over rapidity bins from the left, to find out minimum rapidity of tiling
  double cumul_lo = 0;
  double cumul2 = 0;
  for (ibin = 0; ibin < nbins; ibin++) {
    cumul_lo += counts[ibin];
    if (cumul_lo >= allowed_max_cumul) {
      double y = ibin-nrap;
      if (y > minrap) minrap = y;
      break;
    }
  }
  assert(ibin != nbins); // internal consistency check that you found a bin
  cumul2 += cumul_lo*cumul_lo;

  // ibin_lo is the index of the leftmost bin that should be considered
  int ibin_lo = ibin;

  // then do it from right, to find out maximum rapidity of tiling
  double cumul_hi = 0;
  for (ibin = nbins-1; ibin >= 0; ibin--) {
    cumul_hi += counts[ibin];
    if (cumul_hi >= allowed_max_cumul) {
      double y = ibin-nrap+1; // +1 here is the rapidity bin width
      if (y < maxrap) maxrap = y;
      break;
    }
  }
  assert(ibin >= 0); // internal consistency check that you found a bin

  // ibin_hi is the index of the rightmost bin that should be considered
  int ibin_hi = ibin;

  // consistency check
  assert(ibin_hi >= ibin_lo);

  // now work out cumul2
  if (ibin_hi == ibin_lo) {
    // if there is a single bin (potentially including overflows
    // from both sides), cumul2 is the square of the total contents
    // of that bin, which we obtain from cumul_lo and cumul_hi minus
    // the double counting of part that is contained in both
    // (putting double
    cumul2 = pow(double(cumul_lo + cumul_hi - counts[ibin_hi]), 2);
  } else {
    // otherwise we have a straightforward sum of squares of bin
    // contents
    cumul2 += cumul_hi*cumul_hi;

    // now get the rest of the squared bin contents
    for (ibin = ibin_lo+1; ibin < ibin_hi; ibin++) {
      cumul2 += counts[ibin]*counts[ibin];
    }
  }
  return std::make_tuple(minrap, maxrap);
}


//----------------------------------------------------------------------
void _tj_remove_from_tiles(Tiling& tiling, TiledJet * const jet) {
  Tile * tile = & tiling._tiles[jet->tile_index];

  if (jet->previous == NULL) {
    // we are at head of the tile, so reset it.
    // If this was the only jet on the tile then tile->head will now be NULL
    tile->head = jet->next;
  } else {
    // adjust link from previous jet in this tile
    jet->previous->next = jet->next;
  }
  if (jet->next != NULL) {
    // adjust backwards-link from next jet in this tile
    jet->next->previous = jet->previous;
  }
}

//----------------------------------------------------------------------
/// Set up the tiles:
///  - decide the range in eta
///  - allocate the tiles
///  - set up the cross-referencing info between tiles
///
/// The neighbourhood of a tile is set up as follows
///
///	      LRR
///           LXR
///           LLR
///
/// such that tiles is an array containing XLLLLRRRR with pointers
///                                         |   \ RH_tiles
///                                         \ surrounding_tiles
///
/// with appropriate precautions when close to the edge of the tiled
/// region.
///
Tiling initialise_tiling(const std::vector<PseudoJet>& particles, double Rparam) {

  Tiling tiling;

  // first decide tile sizes (with a lower bound to avoid huge memory use with
  // very small R)
  double default_size = max(0.1, Rparam);
  tiling._tile_size_eta = default_size;
  // it makes no sense to go below 3 tiles in phi -- 3 tiles is
  // sufficient to make sure all pair-wise combinations up to pi in
  // phi are possible
  tiling._n_tiles_phi   = max(3,int(floor(twopi/default_size)));
  tiling._tile_size_phi = twopi / tiling._n_tiles_phi; // >= Rparam and fits in 2pi

  std::tie(tiling._tiles_eta_min, tiling._tiles_eta_max)  = determine_rapidity_extent(particles);


  // now adjust the values
  tiling._tiles_ieta_min = int(floor(tiling._tiles_eta_min/tiling._tile_size_eta));
  tiling._tiles_ieta_max = int(floor(tiling._tiles_eta_max/tiling._tile_size_eta));
  tiling._tiles_eta_min = tiling._tiles_ieta_min * tiling._tile_size_eta;
  tiling._tiles_eta_max = tiling._tiles_ieta_max * tiling._tile_size_eta;

  // allocate the tiles
  tiling._tiles.resize((tiling._tiles_ieta_max-tiling._tiles_ieta_min+1)*tiling._n_tiles_phi);

  // now set up the cross-referencing between tiles
  for (auto ieta = tiling._tiles_ieta_min; ieta <= tiling._tiles_ieta_max; ieta++) {
    for (auto iphi = 0; iphi < tiling._n_tiles_phi; iphi++) {
      auto tile = &tiling._tiles[_tile_index(tiling, ieta,iphi)];
      // no jets in this tile yet
      tile->head = NULL; // first element of tiles points to itself
      tile->begin_tiles[0] =  tile;
      auto pptile = & (tile->begin_tiles[0]);
      pptile++;
      //
      // set up L's in column to the left of X
      //tile->surrounding_tiles = pptile;
      if (ieta > tiling._tiles_ieta_min) {
	// with the _tile_index subroutine, we can safely run tiles from
	// idphi=-1 to idphi=+1, because it takes care of
	// negative and positive boundaries
	for (auto idphi = -1; idphi <=+1; idphi++) {
	  *pptile = & tiling._tiles[_tile_index(tiling, ieta-1,iphi+idphi)];
	  pptile++;
	}
      }
      // now set up last L (below X)
      *pptile = & tiling._tiles[_tile_index(tiling, ieta,iphi-1)];
      pptile++;
      // set up first R (above X)
      tile->RH_tiles = pptile;
      *pptile = & tiling._tiles[_tile_index(tiling, ieta,iphi+1)];
      pptile++;
      // set up remaining R's, to the right of X
      if (ieta < tiling._tiles_ieta_max) {
	for (auto idphi = -1; idphi <= +1; idphi++) {
	  *pptile = & tiling._tiles[_tile_index(tiling, ieta+1,iphi+idphi)];
	  pptile++;
	}
      }
      // now put semaphore for end tile
      tile->end_tiles = pptile;
      // finally make sure tiles are untagged
      tile->tagged = false;
    }
  }
  return tiling;
}


//----------------------------------------------------------------------
// Set up a TiledJet
inline void _tj_set_jetinfo(TiledJet* const jet,
			    ClusterSequence& cs,
			    const int jets_index, double R2) {


  jet->eta  = cs._jets[jets_index].rap();
  jet->phi  = cs._jets[jets_index].phi_02pi();
  jet->kt2  = cs._jets[jets_index].pt2() > 1.e-300 ? 1. / cs._jets[jets_index].pt2(): 1.e300;
  jet->_jets_index = jets_index;
  // initialise NN info as well
  jet->NN_dist = R2;
  jet->NN      = NULL;

  // Find out which tile it belonds to
  jet->tile_index = _tile_index(cs.tiling, jet->eta, jet->phi);

  // Insert it into the tile's linked list of jets
  Tile * tile = &cs.tiling._tiles[jet->tile_index];
  jet->previous   = NULL;
  jet->next       = tile->head;
  if (jet->next != NULL) {jet->next->previous = jet;}
  tile->head      = jet;
}


//----------------------------------------------------------------------
/// output the contents of the tiles
void _print_tiles(const Tiling& tiling, TiledJet * briefjets ){
  for (vector<Tile>::const_iterator tile = tiling._tiles.begin();
       tile < tiling._tiles.end(); tile++) {
    cout << "Tile " << tile - tiling._tiles.begin()<<" = ";
    vector<int> list;
    for (TiledJet * jetI = tile->head; jetI != NULL; jetI = jetI->next) {
      list.push_back(jetI-briefjets);
      //cout <<" "<<jetI-briefjets;
    }
    sort(list.begin(),list.end());
    for (unsigned int i = 0; i < list.size(); i++) {cout <<" "<<list[i];}
    cout <<"\n";
  }
}


//----------------------------------------------------------------------
/// Like _add_neighbours_to_tile_union, but only adds neighbours if
/// their "tagged" status is false; when a neighbour is added its
/// tagged status is set to true.
///
/// Note that with a high level of warnings (-pedantic -Wextra -ansi,
/// gcc complains about tile_index maybe being used uninitialised for
/// oldB in _minheap_faster_tiled_N2_cluster(). We
/// have explicitly checked that it was harmless so we could disable
/// the gcc warning by hand using the construct below
///
///  #pragma GCC diagnostic push
///  #pragma GCC diagnostic ignored "-Wpragmas"
///  #pragma GCC diagnostic ignored "-Wuninitialized"
///  #pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
///    ...
///  #pragma GCC diagnostic pop
///
/// the @GCC diagnostic push/pop directive was only introduced in
/// gcc-4.6, so for broader usage, we'd need to insert #pragma GCC
/// diagnostic ignored "-Wpragmas" at the top of this file
inline void _add_untagged_neighbours_to_tile_union(Tiling& tiling,
						   const int tile_index,
						   vector<int>& tile_union,
						   int& n_near_tiles)  {
  for (Tile ** near_tile = tiling._tiles[tile_index].begin_tiles;
       near_tile != tiling._tiles[tile_index].end_tiles; near_tile++){
    if (! (*near_tile)->tagged) {
      (*near_tile)->tagged = true;
      // get the tile number
      tile_union[n_near_tiles] = *near_tile - &(tiling._tiles[0]);
      n_near_tiles++;
    }
  }
}


//----------------------------------------------------------------------
// initialise the history in a standard way
void _fill_initial_history (ClusterSequence& cs) {

  // reserve sufficient space for everything
  //_jets.reserve(_jets.size()*2);
  cs._history.reserve(cs._jets.size()*2);

  cs._Qtot = 0;

  for (int i = 0; i < static_cast<int>(cs._jets.size()) ; i++) {
    history_element element;
    element.parent1 = InexistentParent;
    element.parent2 = InexistentParent;
    element.child   = Invalid;
    element.jetp_index = i;
    element.dij     = 0.0;
    element.max_dij_so_far = 0.0;

    cs._history.push_back(element);

    // do any momentum preprocessing needed by the recombination scheme
    //_jet_def.recombiner()->preprocess(_jets[i]); No preprocessinhg for E-scheme.

    // get cross-referencing right from PseudoJets
    cs._jets[i].set_cluster_hist_index(i);
    //cs._set_structure_shared_ptr(cs._jets[i]); //FIXME. Mandatory ?!?

    // determine the total energy in the event
    cs._Qtot += cs._jets[i].E();
  }
  //cs._initial_n = cs._jets.size(); //FIXME: needed?!?
  //cs._deletes_self_when_unused = false; //FIXME this field can be removed ?
}


//----------------------------------------------------------------------
// initialise the history in a standard way
void _add_step_to_history(ClusterSequence& cs,
			  const int parent1,
			  const int parent2, const int jetp_index,
			  const double dij) {

  history_element element;
  element.parent1 = parent1;
  element.parent2 = parent2;
  element.jetp_index = jetp_index;
  element.child = Invalid;
  element.dij   = dij;
  element.max_dij_so_far = max(dij,cs._history[cs._history.size()-1].max_dij_so_far);
  cs._history.push_back(element);

  int local_step = cs._history.size()-1;
  //#ifndef __NO_ASSERTS__
  //assert(local_step == step_number);
  //#endif

  // sanity check: make sure the particles have not already been recombined
  //
  // Note that good practice would make this an assert (since this is
  // a serious internal issue). However, we decided to throw an
  // InternalError so that the end user can decide to catch it and
  // retry the clustering with a different strategy.

  assert(parent1 >= 0);
  if (cs._history[parent1].child != Invalid){
    throw std::runtime_error("Internal error. Trying to recombine an object that has previsously been recombined");
  }
  cs._history[parent1].child = local_step;
  if (parent2 >= 0) {
    if (cs._history[parent2].child != Invalid){
      throw std::runtime_error("Internal error. Trying to recombine an object that has previsously been recombined");
    }
    cs._history[parent2].child = local_step;
  }

  // get cross-referencing right from PseudoJets
  if (jetp_index != Invalid) {
    assert(jetp_index >= 0);
    //cout << _jets.size() <<" "<<jetp_index<<"\n";
    cs._jets[jetp_index].set_cluster_hist_index(local_step);
    //_set_structure_shared_ptr(_jets[jetp_index]); FIXME: Needed ???
  }

  //if (_writeout_combinations) {
  //  cout << local_step << ": "
  //	 << parent1 << " with " << parent2
  //	 << "; y = "<< dij<<endl;
  //}

}


//======================================================================
/// carries out the bookkeeping associated with the step of recombining
/// jet_i and jet_j (assuming a distance dij) and returns the index
/// of the recombined jet, newjet_k.
void _do_ij_recombination_step(ClusterSequence& cs,
			       const int jet_i, const int jet_j,
			       const double dij,
			       int & newjet_k) {

  // Create the new jet by recombining the first two with
  // the E-scheme
  //

  //  cout << "Recombine jets with pts "
  //       << sqrt(cs._jets[jet_i].px()* cs._jets[jet_i].px()
  //	       + cs._jets[jet_i].py()* cs._jets[jet_i].py())
  //       << " and "
  //       << sqrt(cs._jets[jet_j].px()* cs._jets[jet_j].px()
  //	       + cs._jets[jet_j].py()*cs._jets[jet_j].py()) << "\n";

  cs._jets.emplace_back(cs._jets[jet_i].px() + cs._jets[jet_j].px(),
			cs._jets[jet_i].py() + cs._jets[jet_j].py(),
			cs._jets[jet_i].pz() + cs._jets[jet_j].pz(),
			cs._jets[jet_i].E()  + cs._jets[jet_j].E());
  // get its index
  newjet_k = cs._jets.size()-1;

  // get history index
  int newstep_k = cs._history.size();
  // and provide jet with the info
  cs._jets[newjet_k].set_cluster_hist_index(newstep_k);

  // finally sort out the history
  int hist_i = cs._jets[jet_i].cluster_hist_index();
  int hist_j = cs._jets[jet_j].cluster_hist_index();

  _add_step_to_history(cs, min(hist_i, hist_j), max(hist_i,hist_j),
		       newjet_k, dij);
}

//======================================================================
/// carries out the bookkeeping associated with the step of recombining
/// jet_i with the beam
void _do_iB_recombination_step(ClusterSequence& cs,
			       const int jet_i, const double diB) {
  // recombine the jet with the beam
  _add_step_to_history(cs,
		       cs._jets[jet_i].cluster_hist_index(), BeamJet,
		       Invalid, diB);

  // // get history index
  // int newstep_k = _history.size();
  //
  // _add_step_to_history(newstep_k,_jets[jet_i].cluster_hist_index(),BeamJet,
  //		       Invalid, diB);

}

//----------------------------------------------------------------------
// return all inclusive jets of a ClusterSequence with pt > ptmin
std::vector<PseudoJet> inclusive_jets (ClusterSequence& cs,
				       const double ptmin = 0.){
  double dcut = ptmin*ptmin;
  int i = cs._history.size() - 1; // last jet
  std::vector<PseudoJet> jets_local;
  // for inclusive jets with a plugin algorithm, we make no
  // assumptions about anything (relation of dij to momenta,
  // ordering of the dij, etc.)
  while (i >= 0) {
    if (cs._history[i].parent2 == BeamJet) {
      int parent1 = cs._history[i].parent1;
      const PseudoJet & jet = cs._jets[cs._history[parent1].jetp_index];
      if (jet.pt2() >= dcut) {jets_local.push_back(jet);}
    }
    i--;
  }
  return jets_local;
}


//----------------------------------------------------------------------
/// run a tiled clustering
std::vector<PseudoJet>
_faster_tiled_N2_cluster(const std::vector<PseudoJet>& particles,
			 double Rparam, double ptmin = 0.) {

  double R2 = Rparam * Rparam;
  double invR2 = 1./R2;

  ClusterSequence cs;

  // this will ensure that we can point to jets without difficulties
  // arising.
  cs._jets.reserve(particles.size()*2);

  // insert initial jets this way so that any type L that can be
  // converted to a pseudojet will work fine (basically PseudoJet
  // and any type that has [] subscript access to the momentum
  // components, such as CLHEP HepLorentzVector).
  //std::copy(particles.begin(), particles.end(), std::back_inserter(cs._jets.begin()));
  for(const auto p: particles) cs._jets.push_back(p);

  _fill_initial_history(cs);

  cs.tiling = initialise_tiling(particles, Rparam);

  int n = cs._jets.size();
  //FIXEME: use a vector for an automatic memory deallocation ?
  TiledJet * tiledjets = new TiledJet[n];
  TiledJet * jetA = tiledjets, * jetB;
  TiledJet oldB;
  oldB.tile_index=0; // prevents a gcc warning

  // will be used quite deep inside loops, but declare it here so that
  // memory (de)allocation gets done only once
  std::vector<int> tile_union(3*Tile::n_tile_neighbours);

  // initialise the basic jet info
  for (int i = 0; i< n; i++) {
    _tj_set_jetinfo(jetA, cs, i, R2);
    //cout << i<<": "<<jetA->tile_index<<"\n";
    jetA++; // move on to next entry of tiledjets
  }
  TiledJet * head = tiledjets; // a nicer way of naming start

  // set up the initial nearest neighbour information
  for (auto& tile:  cs.tiling._tiles){
    // first do it on this tile
    for (jetA = tile.head; jetA != NULL; jetA = jetA->next) {
      for (jetB = tile.head; jetB != jetA; jetB = jetB->next) {
	double dist = _tj_dist(jetA,jetB);
	if (dist < jetA->NN_dist) {jetA->NN_dist = dist; jetA->NN = jetB;}
	if (dist < jetB->NN_dist) {jetB->NN_dist = dist; jetB->NN = jetA;}
      }
    }
    // then do it for RH tiles
    for (Tile ** RTile = tile.RH_tiles; RTile != tile.end_tiles; RTile++) {
      for (jetA = tile.head; jetA != NULL; jetA = jetA->next) {
	for (jetB = (*RTile)->head; jetB != NULL; jetB = jetB->next) {
	  double dist = _tj_dist(jetA,jetB);
	  if (dist < jetA->NN_dist) {jetA->NN_dist = dist; jetA->NN = jetB;}
	  if (dist < jetB->NN_dist) {jetB->NN_dist = dist; jetB->NN = jetA;}
	}
      }
    }
    // no need to do it for LH tiles, since they are implicitly done
    // when we set NN for both jetA and jetB on the RH tiles.
  }

  // now create the diJ (where J is i's NN) table -- remember that
  // we differ from standard normalisation here by a factor of R2
  // (corrected for at the end).
  struct diJ_plus_link {
    double     diJ; // the distance
    TiledJet * jet; // the jet (i) for which we've found this distance
    // (whose NN will the J).
  };
  //FIXME: use a vector for automatic memory deallocation?
  diJ_plus_link * diJ = new diJ_plus_link[n];
  jetA = head;
  for (int i = 0; i < n; i++) {
    diJ[i].diJ = _tj_diJ(jetA); // kt distance * R^2
    diJ[i].jet = jetA;  // our compact diJ table will not be in
    jetA->diJ_posn = i; // one-to-one corresp. with non-compact jets,
    // so set up bi-directional correspondence here.
    //    cout << "clust_dist of jet " << (jetA-head)
    //	 << " with pt " << cs._jets[jetA->_jets_index].pt() <<": "
    //	 << diJ[i].diJ * invR2 << "\n";
    //    if(jetA->NN) cout << "Pt of the nearest jet: " << cs._jets[jetA->NN->_jets_index].pt()
    //		      << "Geom dist.: " <<  jetA->NN_dist << "\n";
    jetA++; // have jetA follow i
  }

  // now run the recombination loop
  int history_location = n-1;
  while (n > 0) {

    // find the minimum of the diJ on this round
    diJ_plus_link * best, *stop; // pointers a bit faster than indices
    // could use best to keep track of diJ min, but it turns out to be
    // marginally faster to have a separate variable (avoids n
    // dereferences at the expense of n/2 assignments).
    double diJ_min = diJ[0].diJ; // initialise the best one here.
    best = diJ;                  // and here
    stop = diJ+n;
    for (diJ_plus_link * here = diJ+1; here != stop; here++) {
      if (here->diJ < diJ_min) {best = here; diJ_min  = here->diJ;}
    }

    // do the recombination between A and B
    history_location++;
    jetA = best->jet;
    jetB = jetA->NN;
    // put the normalisation back in
    diJ_min *= invR2;

    if (jetB != NULL) {
      // jet-jet recombination
      // If necessary relabel A & B to ensure jetB < jetA, that way if
      // the larger of them == newtail then that ends up being jetA and
      // the new jet that is added as jetB is inserted in a position that
      // has a future!
      if (jetA < jetB) {
	std::swap(jetA,jetB);
      }

      int nn; // new jet index
      _do_ij_recombination_step(cs, jetA->_jets_index, jetB->_jets_index, diJ_min, nn);
      _tj_remove_from_tiles(cs.tiling, jetA);
      oldB = * jetB;  // take a copy because we will need it...
      _tj_remove_from_tiles(cs.tiling, jetB);
      _tj_set_jetinfo(jetB, cs, nn, R2); // cause jetB to become _jets[nn]
      //                               (also registers the jet in the tiling)
    } else {
      // jet-beam recombination
      // get the hist_index
      _do_iB_recombination_step(cs, jetA->_jets_index, diJ_min);
      _tj_remove_from_tiles(cs.tiling, jetA);
    }

    // first establish the set of tiles over which we are going to
    // have to run searches for updated and new nearest-neighbours --
    // basically a combination of vicinity of the tiles of the two old
    // and one new jet.
    int n_near_tiles = 0;
    _add_untagged_neighbours_to_tile_union(cs.tiling,
					   jetA->tile_index,
					   tile_union, n_near_tiles);
    if (jetB != NULL) {
      if (jetB->tile_index != jetA->tile_index) {
	_add_untagged_neighbours_to_tile_union(cs.tiling,
					       jetB->tile_index,
					       tile_union,n_near_tiles);
      }
      if (oldB.tile_index != jetA->tile_index &&
	  oldB.tile_index != jetB->tile_index) {
	_add_untagged_neighbours_to_tile_union(cs.tiling,
					       oldB.tile_index,
					       tile_union,n_near_tiles);
      }
    }

    // now update our nearest neighbour info and diJ table
    // first reduce size of table
    n--;
    // then compactify the diJ by taking the last of the diJ and copying
    // it to the position occupied by the diJ for jetA
    diJ[n].jet->diJ_posn = jetA->diJ_posn;
    diJ[jetA->diJ_posn] = diJ[n];

    // Initialise jetB's NN distance as well as updating it for
    // other particles.
    // Run over all tiles in our union
    for (int itile = 0; itile < n_near_tiles; itile++) {
      Tile * tile_ptr = &cs.tiling._tiles[tile_union[itile]];
      tile_ptr->tagged = false; // reset tag, since we're done with unions
      // run over all jets in the current tile
      for (TiledJet * jetI = tile_ptr->head; jetI != NULL; jetI = jetI->next) {
	// see if jetI had jetA or jetB as a NN -- if so recalculate the NN
	if (jetI->NN == jetA || (jetI->NN == jetB && jetB != NULL)) {
	  jetI->NN_dist = R2;
	  jetI->NN      = NULL;
	  // now go over tiles that are neighbours of I (include own tile)
	  for (Tile ** near_tile  = tile_ptr->begin_tiles;
	       near_tile != tile_ptr->end_tiles; near_tile++) {
	    // and then over the contents of that tile
	    for (TiledJet * jetJ  = (*near_tile)->head;
		 jetJ != NULL; jetJ = jetJ->next) {
	      double dist = _tj_dist(jetI,jetJ);
	      if (dist < jetI->NN_dist && jetJ != jetI) {
		jetI->NN_dist = dist; jetI->NN = jetJ;
	      }
	    }
	  }
	  diJ[jetI->diJ_posn].diJ = _tj_diJ(jetI); // update diJ kt-dist
	}
	// check whether new jetB is closer than jetI's current NN and
	// if jetI is closer than jetB's current (evolving) nearest
	// neighbour. Where relevant update things
	if (jetB!= NULL) {
	  double dist = _tj_dist(jetI,jetB);
	  if (dist < jetI->NN_dist) {
	    if (jetI != jetB) {
	      jetI->NN_dist = dist;
	      jetI->NN = jetB;
	      diJ[jetI->diJ_posn].diJ = _tj_diJ(jetI); // update diJ...
	    }
	  }
	  if (dist < jetB->NN_dist) {
	    if (jetI != jetB) {
	      jetB->NN_dist = dist;
	      jetB->NN      = jetI;}
	  }
	}
      }
    } //next itile

      // finally, register the updated kt distance for B
    if (jetB != NULL) {diJ[jetB->diJ_posn].diJ = _tj_diJ(jetB);}

  }

  // final cleaning up;
  delete[] diJ;
  delete[] tiledjets;

  return inclusive_jets(cs, ptmin);
}

void in_mem_process(const char* fname, long long maxevents = -1,
		    long long ntests = 1){

  HepMC3::ReaderAscii input_file (fname);

  int events_parsed = 0;

  std::vector<std::vector<PseudoJet>> events;

  struct timeval t0, t1;
  gettimeofday(&t0, 0);
  while(!input_file.failed()) {

    if(maxevents >= 0 && events_parsed >= maxevents) break;

    std::vector<PseudoJet> input_particles;

    HepMC3::GenEvent evt(HepMC3::Units::GEV, HepMC3::Units::MM);

    // Read event from input file
    input_file.read_event(evt);

    // If reading failed - exit loop
    if( input_file.failed() ) break;

    ++events_parsed;
    input_particles.clear();
    input_particles.reserve(evt.particles().size());
    for(auto p: evt.particles()){
      if(p->status() == 1){
	input_particles.emplace_back(p->momentum().px(),
				     p->momentum().py(),
				     p->momentum().pz(),
				     p->momentum().e());
      }
    }

    events.push_back(input_particles);
  }
  gettimeofday(&t1, 0);
  std::cout << "Time to read " << events.size() << " event(s): " << 1.e-3*(1.e6*(t1.tv_sec-t0.tv_sec)
									   + (t1.tv_usec-t0.tv_usec))
	    << "ms\n\n";

  if(events.size() == 0) return;


  //warm up
  std::vector<PseudoJet> inclusive_jets =
    _faster_tiled_N2_cluster(events[0], R, ptmin);
  // label the columns
  printf("Jets in the first processed event:\n");
  printf("%5s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt");

  // print out the details for each jet
  for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
    printf("%5u %15.8f %15.8f %15.8f\n",
	   i, inclusive_jets[i].rap(), inclusive_jets[i].phi(),
	   inclusive_jets[i].pt());
  }
  printf("\n");

  int njet_acc = 0;

  gettimeofday(&t0, 0);
  for(int itest = 0; itest < ntests; ++itest){
    for(const auto& evt: events){
      auto jets = _faster_tiled_N2_cluster(evt, R, ptmin);
      njet_acc += jets.size();
    }
  }
  gettimeofday(&t1, 0);
  std::cout << "Duration per event: " << (1.e6*(t1.tv_sec-t0.tv_sec)
					  + (t1.tv_usec-t0.tv_usec)) / ntests / events.size()
	    << " us\n";

  std::cout << "Number of processed events: " << events.size() << "\n";
  std::cout << "Sum of jet multiplicity over the events and processing loops: " <<  njet_acc << "\n";

  //  std::sort(inclusive_jets.begin(), inclusive_jets.end(),
  //	    [](auto a, auto b){ return a.pt() > b.pt(); });

}

void print_help(const char* progname){
    std::cout << "Usage: " << progname << " <HepMC3_input_file> [max_events [ntests]]" << std::endl;
}

int main(int argc, char* argv[]){

  if( argc < 2 ) {
    print_help(argv[0]);
    return 1;
  }

  if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0){
    print_help(argv[0]);
    return 0;
  }

  long long maxevents = -1;
  long long ntests = 1;

  if(argc > 2) maxevents = strtoul(argv[2], 0, 0);
  if(argc > 3) ntests = strtoul(argv[3], 0, 0);

  FILE* f = fopen(argv[1], "r");
  if(f == nullptr){
    std::cerr << "Cannot open file '" << argv[1] << "': " << strerror(errno) << ".\n";
    return 1;
  } else{
    fclose(f);
  }
      
  
  in_mem_process(argv[1], maxevents, ntests);
}
