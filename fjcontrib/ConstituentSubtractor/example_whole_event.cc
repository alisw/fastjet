// $Id: example_whole_event.cc 1141 2018-06-28 12:29:14Z berta $
//
//----------------------------------------------------------------------
// Example on how to do pileup correction on the whole event
//
// run it with
//  ./example_whole_event < ../data/Pythia-Zp2jets-lhc-pileup-1ev.dat
//----------------------------------------------------------------------
//
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------


#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "ConstituentSubtractor.hh" // In external code, this should be fastjet/contrib/ConstituentSubtractor.hh 

#include "functions.hh"


using namespace std;
using namespace fastjet;

//----------------------------------------------------------------------
int main(){
  // set up before event loop:
  double max_eta=4; // specify the maximal pseudorapidity for the input particles. It is used for the subtraction. Particles with eta>|max_eta| are removed and not used during the subtraction (they are not returned). The same parameter should be used for the GridMedianBackgroundEstimator as it is demonstrated in this example. If JetMedianBackgroundEstimator is used, then lower parameter should be used  (to avoid including particles outside this range). 
  double max_eta_jet=3; // the maximal pseudorapidity for selected jets. Not used for the subtraction - just for the final output jets.
  contrib::ConstituentSubtractor subtractor;
  subtractor.set_distance_type(contrib::ConstituentSubtractor::deltaR); // free parameter for the type of distance between particle i and ghost k. There are two options: "deltaR" or "angle" which are defined as deltaR=sqrt((y_i-y_k)^2+(phi_i-phi_k)^2) or Euclidean angle between the momenta
  subtractor.set_max_distance(0.3); // free parameter for the maximal allowed distance between particle i and ghost k
  subtractor.set_alpha(1);  // free parameter for the distance measure (the exponent of particle pt). The larger the parameter alpha, the more are favoured the lower pt particles in the subtraction process
  subtractor.set_ghost_area(0.01); // free parameter for the density of ghosts. The smaller, the better - but also the computation is slower.
  subtractor.set_max_eta(max_eta); // parameter for the maximal eta cut

  // example selector for ConstituentSubtractor:
  Selector sel_CS_correction = SelectorPhiRange(0,3) * SelectorEtaMin(-1.5) * SelectorEtaMax(0);
  // subtractor.set_selector(&sel_CS_correction);             // uncomment in case you want to use selector. Only ghosts fulfilling this selector will be constructed. No selection on input particles is done. The selection on input particles is the responsibility of the user. However, the maximal eta cut (specified with "set_max_eta" function) is done always for both, ghosts and input particles.

  subtractor.initialize();



  // in event loop:
  // read in input particles
  vector<PseudoJet> hard_event, full_event;
  read_event(hard_event, full_event);

  hard_event = SelectorAbsEtaMax(max_eta)(hard_event);
  full_event = SelectorAbsEtaMax(max_eta)(full_event);

  cout << "# read an event with " << hard_event.size() << " signal particles and " << full_event.size() - hard_event.size() << " background particles with pseudo-rapidity |eta|<4" << endl;

 
  // jet definition and selector for jets
  JetDefinition jet_def(antikt_algorithm, 0.7);
  Selector sel_jets = SelectorNHardest(3) * SelectorAbsEtaMax(max_eta_jet);

  // clustering of the hard-scatter event and uncorrected event
  ClusterSequence clust_seq_hard(hard_event, jet_def);
  ClusterSequence clust_seq_full(full_event, jet_def);
  vector<PseudoJet> hard_jets = sel_jets(clust_seq_hard.inclusive_jets());
  vector<PseudoJet> full_jets = sel_jets(clust_seq_full.inclusive_jets());

  // background estimator
  GridMedianBackgroundEstimator bge_rho(max_eta,0.5);
  bge_rho.set_particles(full_event);
  subtractor.set_background_estimator(&bge_rho);

  // this sets the same background estimator to be used for deltaMass density, rho_m, as for pt density, rho:
  subtractor.set_common_bge_for_rho_and_rhom(true); // it does not make any difference for massless input particles (rho_m is always zero)

  // print info (optional)
  cout << subtractor.description() << endl;

  // the correction of the whole event with ConstituentSubtractor
  vector<PseudoJet> corrected_event=subtractor.subtract_event(full_event);

  // clustering of the corrected event
  ClusterSequence clust_seq_corr(corrected_event, jet_def);
  vector<PseudoJet> corrected_jets = sel_jets(clust_seq_corr.inclusive_jets());


  ios::fmtflags f( cout.flags() );
  cout << setprecision(4) << fixed;
  cout << endl << "Corrected particles in the whole event:" << endl;
  for (unsigned int i=0; i<corrected_event.size(); i++){
    const PseudoJet &particle = corrected_event[i];
    cout << "pt = " << particle.pt()
         << ", phi = " << particle.phi()
         << ", rap = " << particle.rap()
         << ", |mass| = " << fabs(particle.m()) << endl;
  }
  cout << endl;

  // shape variables:
  //----------------------------------------------------------
  JetWidth width;

  // subtract and print the result
  //----------------------------------------------------------
  cout.flags( f );
  cout << setprecision(4);
  cout << "# original hard jets" << endl;
  for (unsigned int i=0; i<hard_jets.size(); i++){
    const PseudoJet &jet = hard_jets[i];
    cout << "pt = " << jet.pt()
	 << ", rap = " << jet.rap()
	 << ", mass = " << jet.m()
	 << ", width = " << width(jet) << endl;
  }
  cout << endl;

  cout << "# unsubtracted full jets" << endl;
  for (unsigned int i=0; i<full_jets.size(); i++){
    const PseudoJet &jet = full_jets[i];
    cout << "pt = " << jet.pt()
	 << ", rap = " << jet.rap()
	 << ", mass = " << jet.m()
	 << ", width = " << width(jet) << endl;
  }
  cout << endl;

  cout << "# subtracted full jets" << endl;
  for (unsigned int i=0; i<corrected_jets.size(); i++){
    const PseudoJet &jet = corrected_jets[i];

    cout << "pt = " << jet.pt()
	 << ", rap = " << jet.rap()
	 << ", mass = " << jet.m()
	 << ", width = " << width(jet) << endl;
	 }
  cout << endl;

  return 0;
}



