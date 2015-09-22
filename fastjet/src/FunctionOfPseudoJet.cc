//STARTHEADER
// $Id: FunctionOfPseudoJet.cc 2577 2011-09-13 15:11:38Z salam $
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

#include <fastjet/FunctionOfPseudoJet.hh>
#include <string>
#include <sstream>

using namespace std;

FASTJET_BEGIN_NAMESPACE

// //----------------------------------------------------------------------
// /// helper for the selector "FunctionOfPJ(jet) > ref"
// class SW_FofPJ_larger : public SelectorWorker {
// public:
//   /// ctor with specification of the number of objects to keep
//   SW_FofPJ_larger(const FunctionOfPseudoJet<double> & fn, const double &ref)
//     : _fn(fn), _ref(ref){}
// 
//   /// return true if the jet is a pure-ghost jet
//   virtual bool pass(const PseudoJet & jet) const { return _fn(jet) > _ref; }
//   
//   /// returns a description of the worker
//   virtual string description() const { 
//     ostringstream oss;
//     oss << "(" << _fn.description() << " > " << _ref << ")";
//     return oss.str();
//   }
// 
// protected:
//   const FunctionOfPseudoJet<double> & _fn;
//   double _ref;
// };
// 
// 
// // select objects that are (or are only made of) ghosts
// Selector operator >(const FunctionOfPseudoJet<double> & fn, const double & ref){
//   return Selector(new SW_FofPJ_larger(fn, ref));
// }
// 
// 
// 
// //----------------------------------------------------------------------
// /// helper for the selector "FunctionOfPJ(jet) < ref"
// class SW_FofPJ_smaller : public SelectorWorker {
// public:
//   /// ctor with specification of the number of objects to keep
//   SW_FofPJ_smaller(const FunctionOfPseudoJet<double> & fn, const double &ref)
//     : _fn(fn), _ref(ref){}
// 
//   /// return true if the jet is a pure-ghost jet
//   virtual bool pass(const PseudoJet & jet) const { return _fn(jet) < _ref; }
//   
//   /// returns a description of the worker
//   virtual string description() const { 
//     ostringstream oss;
//     oss << "(" << _fn.description() << " < " << _ref << ")";
//     return oss.str();
//   }
// 
// protected:
//   const FunctionOfPseudoJet<double> & _fn;
//   double _ref;
// };
// 
// 
// // select objects that are (or are only made of) ghosts
// Selector operator <(const FunctionOfPseudoJet<double> & fn, const double & ref){
//   return Selector(new SW_FofPJ_smaller(fn, ref));
// }
// 
// 
// //----------------------------------------------------------------------
// /// helper for the selector "FunctionOfPJ(jet) >= ref"
// class SW_FofPJ_larger_equal : public SelectorWorker {
// public:
//   /// ctor with specification of the number of objects to keep
//   SW_FofPJ_larger_equal(const FunctionOfPseudoJet<double> & fn, const double &ref)
//     : _fn(fn), _ref(ref){}
// 
//   /// return true if the jet is a pure-ghost jet
//   virtual bool pass(const PseudoJet & jet) const { return _fn(jet) >= _ref; }
//   
//   /// returns a description of the worker
//   virtual string description() const { 
//     ostringstream oss;
//     oss << "(" << _fn.description() << " >= " << _ref << ")";
//     return oss.str();
//   }
// 
// protected:
//   const FunctionOfPseudoJet<double> & _fn;
//   double _ref;
// };
// 
// 
// // select objects that are (or are only made of) ghosts
// Selector operator >=(const FunctionOfPseudoJet<double> & fn, const double & ref){
//   return Selector(new SW_FofPJ_larger_equal(fn, ref));
// }
// 
// 
// 
// //----------------------------------------------------------------------
// /// helper for the selector "FunctionOfPJ(jet) <= ref"
// class SW_FofPJ_smaller_equal : public SelectorWorker {
// public:
//   /// ctor with specification of the number of objects to keep
//   SW_FofPJ_smaller_equal(const FunctionOfPseudoJet<double> & fn, const double &ref)
//     : _fn(fn), _ref(ref){}
// 
//   /// return true if the jet is a pure-ghost jet
//   virtual bool pass(const PseudoJet & jet) const { return _fn(jet) <= _ref; }
//   
//   /// returns a description of the worker
//   virtual string description() const { 
//     ostringstream oss;
//     oss << "(" << _fn.description() << " <= " << _ref << ")";
//     return oss.str();
//   }
// 
// protected:
//   const FunctionOfPseudoJet<double> & _fn;
//   double _ref;
// };
// 
// 
// // select objects that are (or are only made of) ghosts
// Selector operator <=(const FunctionOfPseudoJet<double> & fn, const double & ref){
//   return Selector(new SW_FofPJ_smaller_equal(fn, ref));
// }
// 

FASTJET_END_NAMESPACE
