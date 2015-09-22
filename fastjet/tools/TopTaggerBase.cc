//STARTHEADER
// $Id: TopTaggerBase.cc 2689 2011-11-14 14:51:06Z soyez $
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

#include <fastjet/tools/TopTaggerBase.hh>

FASTJET_BEGIN_NAMESPACE

using namespace std;

// compute the W helicity angle
//
// The helicity angle is a standard observable in top decays, used to
// determine the Lorentz structure of the top- W coupling [13]. It is
// defined as the angle, measured in the rest frame of the
// reconstructed W, between the reconstructed top's flight direction
// and one of the W decay products. Normally, it is studied in
// semi-leptonic top decays, where the charge of the lepton uniquely
// identifies these decay products. In hadronic top decays there is an
// ambiguity which we resolve by choosing the lower pT subjet, as
// measured in the lab frame.
//
// The jet passed to this function is expected to already have
// the structure of a top, including a functional "W()" call;
// the W must be made of two pieces.
double TopTaggerBase::_cos_theta_W(const PseudoJet & res) const{
  // the two jets of interest: top and lower-pt prong of W
  const PseudoJet & W  = res.structure_of<TopTaggerBase>().W();
  vector<PseudoJet> W_pieces = W.pieces();
  assert(W_pieces.size() == 2);
  //assert(W_pieces[0].perp2() >= W_pieces[1].perp2());
  //PseudoJet W2  = W_pieces[1];
  // extract the softer of the two W pieces.
  PseudoJet W2 =  (W_pieces[0].perp2() < W_pieces[1].perp2())
                    ? W_pieces[0] 
                    : W_pieces[1];
  PseudoJet top = res;
  
  // transform these jets into jets in the rest frame of the W
  W2.unboost(W);
  top.unboost(W);

  return (W2.px()*top.px() + W2.py()*top.py() + W2.pz()*top.pz())/
    sqrt(W2.modp2() * top.modp2());
}


FASTJET_END_NAMESPACE
