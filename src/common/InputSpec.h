/*a
 * InputSpec.h
 *
 *  Created on: 17 Feb 2018
 *      Author: nathalie
 */

#ifndef INPUTSPEC_H_
#define INPUTSPEC_H_

#include "SHS.h"
#include "TaskSpec.h"

template <class T, class T2> class inputSpec_t
{
	inputSpec_t(): stoch(true), sigm(false), par(false) {};
	inputSpec_t(shs_t<T, T2> _model, taskSpec_t _task): stoch(true), sigm(false), par(false), myModel(_model), myTask(_task) {};
	private:
		bool stoch;											// identify whether system is deterministic or not
		bool sigm;											// To indicate to use the sigmoid function as the discrete
		// transition kernel
		bool par;											// To indicate multiple initial states
		shs_t<T, T2> myModel;								// obtain state space model
		taskSpec_t myTask;									// store task to be performed

	public:
		virtual ~inputSpec_t() {}
};
#endif
