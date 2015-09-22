//STARTHEADER
// $Id: MassDropTagger.cc 3249 2013-10-28 15:06:25Z salam $
//
// Copyright (c) 2005-2011, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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

#include <fastjet/tools/MassDropTagger.hh>
#include <fastjet/ClusterSequence.hh>
#include <sstream>

FASTJET_BEGIN_NAMESPACE

LimitedWarning MassDropTagger::_warnings_nonca;

using namespace std;

//----------------------------------------------------------------------
// MassDropTagger class implementation
//----------------------------------------------------------------------

//------------------------------------------------------------------------
// description of the tagger
string MassDropTagger::description() const{ 
  ostringstream oss;
  oss << "MassDropTagger with mu=" << _mu << " and ycut=" << _ycut;
  return oss.str();
}

//------------------------------------------------------------------------
// returns the tagged PseudoJet if successful, 0 otherwise
//  - jet   the PseudoJet to tag
PseudoJet MassDropTagger::result(const PseudoJet & jet) const{
  PseudoJet j = jet;

  // issue a warning if the jet is not obtained through a C/A
  // clustering
  if ((! j.has_associated_cluster_sequence()) ||
      (j.validated_cs()->jet_def().jet_algorithm() != cambridge_algorithm))
    _warnings_nonca.warn("MassDropTagger should only be applied on jets from a Cambridge/Aachen clustering; use it with other algorithms at your own risk.");


  PseudoJet j1, j2;
  bool had_parents;

  // we just ask that we can "walk" in the cluster sequence.
  // appropriate errors will be thrown automatically if this is not
  // the case
  while ((had_parents = j.has_parents(j1,j2))) {
    // make parent1 the more massive jet
    if (j1.m2() < j2.m2()) std::swap(j1,j2);

    // if we pass the conditions on the mass drop and its degree of
    // asymmetry (kt_dist/m^2 > rtycut [where kt_dist/m^2 \sim
    // z/(1-z)), then we've found something interesting, so exit the
    // loop
    if ( (j1.m2() < _mu*_mu*j.m2()) && (j1.kt_distance(j2) > _ycut*j.m2()) )
      break;
    else
      j = j1;
  }
    
  if (!had_parents)
    // no Higgs found, return an empty PseudoJet
    return PseudoJet();

  // create the result and its structure
  PseudoJet result_local = j;
  MassDropTaggerStructure * s = new MassDropTaggerStructure(result_local);
//  s->_original_jet = jet;
  s->_mu = (j.m2()!=0.0) ? sqrt(j1.m2()/j.m2()) : 0.0;
  s->_y  = (j.m2()!=0.0) ? j1.kt_distance(j2)/j.m2() : 0.0;

  result_local.set_structure_shared_ptr(SharedPtr<PseudoJetStructureBase>(s));

  return result_local;
}

FASTJET_END_NAMESPACE

