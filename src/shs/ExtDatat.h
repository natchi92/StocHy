/*
 * ExtDatat.h
 *
 *  Created on: 11 Jan 2018
 *      Author: nathalie
 */

#include <armadillo>

#ifndef EXTDATAT_H_
#define EXTDATAT_H_
// To setup store for reading external
// inputs and initial
// values for both continous variables
// and discrete locations for the
// simulator

using 3DMatrix = std::vector<arma::mat>;
using Matrix = arma::mat;
using ZerosMatrix = arma::zeros<arma::mat>
using emptyMatrix = arma::zeros<Matrix>(1,1);

class SimulationData
{
  SimulationData() {};
  SimulationData(3DMatrix _cTVar, Matrix _aTVar, Matrix _dTVar, 3DMatrix _initMode, 3DMatrix _initDis) :
    continuousVariables{_cTVar}, actionVariables{_aTVar}, disturbanceVariables{_dTVar}, initialMode{_initMode},
    x_k{_cTVar.size() > 0 ? _cTVar[0]: emptyMatrix}, u_k{_aTVar.size() > 0 ? _cTVar[0]: emptyMatrix}, d_k{_dTVar.size() > 0 ? _cTVar[0]: emptyMatrix},
    initialModeDistribution{_initDis} {};

private:
  3DMatrix continuousVariables; // Continuous variables: dim - x_dim,T
  Matrix actionVariables;              // Exogenous inputs: dim - u_dim,T
  Matrix disturbanceVariables;              // Disturbance signals: dim - d_dim,T
  3DMatrix initialMode; // initial  value discrete modes
  Matrix x_k;                 // continuous variable evolution
  Matrix u_k;                // input variable evolution
  Matrix d_k;                // disturbance variable evolution
  3DMatrix initialModeDistribution; // initial disturbance

  virtual ~SimulationData(){};
};

#endif /* EXTDATAT_H_ */
