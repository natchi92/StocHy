/*
 * TaskSpec.h
 *
 *  Created on: 17 Feb 2018
 *      Author: nathalie
 */

#ifndef TASKSPEC_H_
#define TASKSPEC_H_
#include <armadillo>
#include "utility.h"

class TaskSpecification
{
	TaskSpecification(): {};
	TaskSpecification(int _timeHorizon, int _runs) :  task(SIMULATION), timeHorizon(_timeHorizon), runs(_runs) {};
	TaskSpecification(int _timeHorizon, int _propertySpec, Matrix _safeSet, Matrix _actionSet, Matrix _targetSet, double _eps) : task(MDP_ABSTRACTION), timeHorizon(_timeHorizon), propertySpec(_propertySpec),
		safeSet(_safeSet), targetSet(_targetSet), eps(_eps) {};
	TaskSpecification(int _timeHorizon, int _propertySpec, Matrix _boundary, Matrix _gridSize, Matrix _refTol) : task(IMDP_ABSTRACTION), timeHorizon(_timeHorizon), propertySpec(_propertySpec),
		boundary(_boundary), gridsize(_gridSize), refTol(_refTol) {};
	void setGridType(const GridType _gridType) { gridType = _gridType; }
	void setKernelAssumption(const KernelType _kernelAssumptions) {kernelAssumptions = _kernelAssumptions; }
	void setControlType(const int _controlType) { controlType = _controlType; }
	bool isControlled { return controlType == 1; }

	private:
		LibraryType task;
		int timeHorizon;									// Time horizon to perform tool
		int runs;											// Number of runs
		// Model checking related
		PropertyType propertySpec;

		// synthesis FAUST related inputs
		arma::mat safeSet;
		arma::mat targetSet;
		arma::mat inputSet;
		double eps;
		GridType gridType;									// 1- Uniform Grid, 2- Local Adaptive, 3 - MCMC approx
		KernelType kernelAssumptions;						// 1- Lipschitz via Integral, 2- Lipschitz via sample,
		// 3- Max-min
		int controlType;
		// BMDP related
		arma::mat boundary;
		arma::mat gridSize;
		arma::mat refTol;

	public:
		virtual ~TaskSpecification();
};
#endif														/* TASKSPEC_H_ */
