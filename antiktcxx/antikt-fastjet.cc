//----------------------------------------------------------------------
/// \file
/// \page Example01 01 - basic usage example
///
/// fastjet basic example program:
///   simplest illustration of the usage of the basic classes:
///   fastjet::PseudoJet, fastjet::JetDefinition and 
///   fastjet::ClusterSequence
///
/// run it with    : ./01-basic < data/single-event.dat
///
/// Source code: 01-basic.cc
//----------------------------------------------------------------------

//STARTHEADER
// $Id: 01-basic.cc 4354 2018-04-22 07:12:37Z salam $
//
// Copyright (c) 2005-2018, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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
//  development and are described in hep-ph/0512210. If you use
//  FastJet as part of work towards a scientific publication, please
//  include a citation to the FastJet paper.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//ENDHEADER

#include "fastjet/ClusterSequence.hh"
#include <iostream> // needed for io
#include <cstdio>   // needed for io
#include <sys/time.h>

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/ReaderAscii.h"


using namespace std;
using namespace fastjet;

double R = 0.4;
double ptmin = 5.0;
int n_tests = 1000;


/// an example program showing how to use fastjet
//int main(){
//  
//  // read in input particles
//  //----------------------------------------------------------
//  vector<fastjet::PseudoJet> input_particles;
//  
//  double px, py , pz, E;
//  while (cin >> px >> py >> pz >> E) {
//    // create a fastjet::PseudoJet with these components and put it onto
//    // back of the input_particles vector
//    input_particles.push_back(fastjet::PseudoJet(px,py,pz,E)); 
//  }
//  
//
//  // create a jet definition: 
//  // a jet algorithm with a given radius parameter
//  //----------------------------------------------------------
//  double R = 0.4;
//  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
//
//
//  // run the jet clustering with the above jet definition
//  //----------------------------------------------------------
//  fastjet::ClusterSequence clust_seq(input_particles, jet_def);
//
//
//  // get the resulting jets ordered in pt
//  //----------------------------------------------------------
//  double ptmin = 5.0;
//  vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));
//
//
//  // tell the user what was done
//  //  - the description of the algorithm used
//  //  - extract the inclusive jets with pt > 5 GeV
//  //    show the output as 
//  //      {index, rap, phi, pt}
//  //----------------------------------------------------------
//  cout << "Ran " << jet_def.description() << endl;
//
//  // label the columns
//  printf("%5s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt");
// 
//  // print out the details for each jet
//  for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
//    printf("%5u %15.8f %15.8f %15.8f\n",
//	   i, inclusive_jets[i].rap(), inclusive_jets[i].phi(),
//	   inclusive_jets[i].perp());
//  }
//
//  return 0;
//}




void in_mem_process(const char* fname, long long maxevents = -1){
  
  HepMC3::ReaderAscii input_file (fname);

  int events_parsed = 0;

  std::vector<std::vector<PseudoJet>> events;

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
  
  if(events.size() == 0) return;

  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R, E_scheme, N2Tiled);
  
  //warm up
  fastjet::ClusterSequence clust_seq(events[0], jet_def);
  auto inclusive_jets = clust_seq.inclusive_jets(ptmin);
  
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
  struct timeval t0, t1;
  gettimeofday(&t0, 0);
  for(int itest = 0; itest < n_tests; ++itest){
    for(const auto& evt: events){
      fastjet::ClusterSequence clust_seq(evt, jet_def);
      auto jets = clust_seq.inclusive_jets(ptmin);
      njet_acc += jets.size();
    }
  }
  gettimeofday(&t1, 0);
  std::cout << "Duration: " << (1.e6*(t1.tv_sec-t0.tv_sec)
				+ (t1.tv_usec-t0.tv_usec)) / n_tests / events.size()
	    << "us\n";
  
  std::cout << "Number of processed events: " << events.size() << "\n";
  std::cout << "Sum of jet multiplicity over the events and processing loops: " <<  njet_acc << "\n";

  std::sort(inclusive_jets.begin(), inclusive_jets.end(),
	    [](auto a, auto b){ return a.pt() > b.pt(); });
  
}



int main(int argc, char* argv[]){

  if( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " <HepMC3_input_file> [max_events]" << std::endl;
    exit(-1);
  }

  long long maxevents = -1;

  if(argc > 2) maxevents = strtoul(argv[2], 0, 0);

  //  process_file(argv[1], maxevents);
  in_mem_process(argv[1], maxevents);
}
