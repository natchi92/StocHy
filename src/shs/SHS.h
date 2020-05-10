/*
 * shs_t.h
 *
 *  Created on: 14 Nov 2017
 *      Author: nathalie
 */

#include "SSModels.h"
#include <ginac/ginac.h>
#include <random>

#ifndef STOCHY_COMMON_SHS_H
#define STOCHY_COMMON_SHS_H

// General class template to cater for all the different SHS model definitions
template <class T, class T2> class shs_t
{
	shs_t() noexcept: shs{1},
		continuousSpaceDimension
	{
		1
	}
	, discreteModes
	{
	}
	, discreteKernel
	{
	}
	, Pk
	{
		Matrix(discreteModes, continuousSpaceDimension)
	},
		initialDiscreteMode{Matrix(discreteModes, 1)},
		stateSpaceModelsPerMode
	{
	}
	, continuousStateData
	{
	}
	{
	}

	shs_t(const int _continuousSpaceDimension, const T2 _discreteModes, const _discreteKernel, const SimulationData  &_continuousStateData,
		const std::vector<ssmodels_t> &_stateSpaceModelsPerMode) noexcept:
	shs														// TODO(ncauchi) remove this variable very clumsy
	{
		1
	},
		continuousSpaceDimension
	{
		_continuousSpaceDimension
	}
	, discreteModes
	{
		_discreteModes
	},
		discreteKernel{_discreteKernel},
		Pk													// TODO(ncauchi) make better
	{
		_continuousStateData.initialModeDistribution[0]
	},
		initialDiscreteMode{_continuousStateData.initialMode[0]},
		stateSpaceModelsPerMode
	{
		_stateSpaceModelsPerMode
	}
	, continuousStateData
	{
		_continuousStateData
	}
	{
	}

	shs_t(const int _continuousSpaceDimension, const T2 _discreteModes, const _discreteKernel,
		const std::vector<ssmodels_t> &_stateSpaceModelsPerMode) noexcept:
	shs														// TODO(ncauchi) remove this variable very clumsy
	{
		1
	},
		continuousSpaceDimension
	{
		_continuousSpaceDimension
	}
	, discreteModes
	{
		_discreteModes
	},
		discreteKernel{_discreteKernel},
		Pk{},
		initialDiscreteMode
	{
	},
		stateSpaceModelsPerMode
	{
		_stateSpaceModelsPerMode
	}
	, continuousStateData
	{
	}
	{
	}

	shs_t(const _discreteKernel, const SimulationData  &_continuousStateData) noexcept :
	shs														// TODO(ncauchi) remove this variable very clumsy
	{
		1
	},
		continuousSpaceDimension
	{
	}
	, discreteModes
	{
	},
		discreteKernel{_discreteKernel},
		Pk													// TODO(ncauchi) make better
	{
		_continuousStateData.initialModeDistribution[0]
	},
		initialDiscreteMode{_continuousStateData.initialMode[0]},
		stateSpaceModelsPerMode
	{
	}
	, continuousStateData
	{
		_continuousStateData
	}
	{
	}

	public:
		bool shs;											// 1- SHS, 0- HS TODO(ncauchi) remove the need for this
		int continuousSpaceDimension						// Dimension of continuous state space in given mode
			T2 discreteModes;								// Discrete modes
		T discreteKernel;									// Discrete kernel- can represent either guards or probabilities
		Matrix Pk;
		Matrix initialDiscreteMode;							// Modes
		std::vector<ssmodels_t> stateSpaceModelsPerMode;	// container of models
		SimulationData continuousStateData;
		virtual ~shs_t() {}
}


// TODO(ncauchi) remove references to deleted function checkData

/*************************************************************************************************/

// class template specialisation
// for the case when have a hybrid model and the transitions are governed by
// logical guards not probabilities
template<> class shs_t<std::vector<std::string>, int>
{
	/*shs_t() noexcept: shs{1}, continuousSpaceDimension{1}, discreteModes{}, discreteKernel{}, Pk{Matrix(discreteModes, continuousSpaceDimension)},
			  initialDiscreteMode{Matrix(discreteModes, 1)}, stateSpaceModelsPerMode{}, continuousStateData{}  {}

	shs_t(const int _continuousSpaceDimension, const T2 _discreteModes, const _discreteKernel, const SimulationData  &_continuousStateData,
		  const std::vector<ssmodels_t> &_stateSpaceModelsPerMode) noexcept:
		  shs{1}, // TODO(ncauchi) remove this variable very clumsy
	continuousSpaceDimension{_continuousSpaceDimension}, discreteModes{_discreteModes},
	discreteKernel{_discreteKernel},
	Pk{_continuousStateData.initialModeDistribution[0]}, // TODO(ncauchi) make better
	initialDiscreteMode{_continuousStateData.initialMode[0]},
	stateSpaceModelsPerMode{_stateSpaceModelsPerMode}, continuousStateData{_continuousStateData} {}

	shs_t(const int _continuousSpaceDimension, const T2 _discreteModes, const _discreteKernel,
	const std::vector<ssmodels_t> &_stateSpaceModelsPerMode) noexcept:
	shs{1}, // TODO(ncauchi) remove this variable very clumsy
	continuousSpaceDimension{_continuousSpaceDimension}, discreteModes{_discreteModes},
	discreteKernel{_discreteKernel},
	Pk{}, initialDiscreteMode{},
	stateSpaceModelsPerMode{_stateSpaceModelsPerMode}, continuousStateData{} {}

	shs_t(const _discreteKernel, const SimulationData  &_continuousStateData) noexcept :
	shs{1}, // TODO(ncauchi) remove this variable very clumsy
	continuousSpaceDimension{}, discreteModes{},
	discreteKernel{_discreteKernel},
	Pk{_continuousStateData.initialModeDistribution[0]}, // TODO(ncauchi) make better
	initialDiscreteMode{_continuousStateData.initialMode[0]},
	stateSpaceModelsPerMode{}, continuousStateData{_continuousStateData} {}*/
															//TODO(ncauchi) check if these are needed

	shs_t(const char *fn, SimulationData &_continuousStateData)
	{
		// Obtain number of discrete modes and obtain discreteKernel
		if (! obtainTqfromMat(fn, *this))
		{
			// Get transition probabilities associated with the guards
			// and store in pk
			Pk = {Matrix(discreteModes, discreteModes)};
			continuousStateData = _continuousStateData;
															// Initial mode            //  Pk = data.initialModeDistribution;
			initialDiscreteMode = continuousStateData.initialMode[0];

			// create array containing the ssmodels
			std::vector<ssmodels_t> models;
			for (int i = 1; i <= discreteModes; i++)
			{
				int x_dimension = _continuousStateData.x_k.n_rows;
				int u_dimension = _continuousStateData.u_k.n_rows;
				int d_dimension = _continuousStateData.d_k.n_rows;

				ssmodels_t stateSpaceModelsForCurrentMode(x_dimension, u_dimension, d_dimension);

				// Reading models
				std::cout << "INFO: Initialising model of continuous variables in discrete mode "
					<< i << std::endl;

				stateSpaceModelsForCurrentMode.obtainSSfromMat(fn, mod, i);
				stateSpaceModelsForCurrentMode.checkModel();
				stateSpaceModelsPerMode.push_back(stateSpaceModelsForCurrentMode);
			}

			n = stateSpaceModelsPerMode[0].x_dim;			// TODO(ncauchi) make flexible to have  different continuous variable per mode

			// Check dimensions of data
			if (checkData == 3)
			{
				std::cout << "ERROR: Incorrect size of disturbance vector" << std::endl;
			}
		}
		public:
			bool shs;										// 1- SHS, 0- HS TODO(ncauchi) remove the need for this
			int continuousSpaceDimension					// Dimension of continuous state space in given mode
				int discreteModes;							// Discrete modes
			std::vector<std::string> discreteKernel;		// Discrete kernel- can represent either guards or probabilities
			Matrix Pk;
			Matrix initialDiscreteMode;						// Modes
			std::vector<ssmodels_t> stateSpaceModelsPerMode;// container of models
			SimulationData continuousStateData;
			virtual ~shs_t() {}
	}

	//TODO(ncauchi) check if these are needed
	shs_t(const char *fn)
	{
		// Obtain number of discrete modes and obtain discreteKernel
		if (!obtainTqfromMat(fn, *this))
		{
			// Get transition probabilities associated with the guards
			// and store in pk
			Pk = {Matrix(discreteModes, discreteModes)};
			// create array containing the ssmodels
			std::vector<ssmodels_t> models;
			for (int i = 1; i <= discreteModes; i++)
			{
				int x_dimension = stateSpaceModelsPerMode[0].A.n_rows;
				int u_dimension = stateSpaceModelsPerMode[0].B.n_rows
					int d_dimension = stateSpaceModelsPerMode[0].F.n_rows

					ssmodels_t stateSpaceModelsForCurrentMode(x_dimension, u_dimension, d_dimension);

				// Reading models
				std::cout << "INFO: Initialising model of continuous variables in discrete mode "
					<< i << std::endl;

				stateSpaceModelsForCurrentMode.obtainSSfromMat(fn, mod, i);
				stateSpaceModelsForCurrentMode.checkModel();
				stateSpaceModelsPerMode.push_back(stateSpaceModelsForCurrentMode);
			}

			n = stateSpaceModelsPerMode[0].x_dim;			// TODO(ncauchi) make flexible to have  different continuous variable per mode

			// Check dimensions of data
			if (checkData == 3)
			{
				std::cout << "ERROR: Incorrect size of disturbance vector" << std::endl;
			}
		}

		public:
			bool shs;										// 1- SHS, 0- HS TODO(ncauchi) remove the need for this
			int continuousSpaceDimension					// Dimension of continuous state space in given mode
				int discreteModes;							// Discrete modes
			std::vector<std::string> discreteKernel;		// Discrete kernel- can represent either guards or probabilities
			Matrix Pk;
			Matrix initialDiscreteMode;						// Modes
			std::vector<ssmodels_t> stateSpaceModelsPerMode;// container of models
			SimulationData continuousStateData;
			virtual ~shs_t() {}
	}

	// TODO(ncauchi) continue from here

	void step_hyb(shs_t &old, int n, int cond)
	{
		// TODO: Update to cater for multiple discrete modes
		int xdim = old.continuousSpaceDimension
		// Matrix x_old=Matrix(old.stateSpaceModelsPerMode[(int) old.initialDiscreteMode(n,0)].x_dim,1);
		// x_old = old.continuousStateData.continuousVariables[n];
			Matrix x_old(old.continuousStateData.continuousVariables[n]);
		Matrix x_new = Matrix(xdim, 1);

		double q_old = old.initialDiscreteMode(n, 0);
		double q_new = q_old;
		// std::vector<double> vec(x_old.size());
		// Matrix x_old(vec);
		// Eigen::Map<Matrix>(vec.data(), x_old.n_rows, x_old.n_cols) =x_old;
		// //TODO fix

		Matrix u_old = Matrix(old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode(n, 0)].u_dim, 1);
		u_old = old.continuousStateData.u_k;
		Matrix d_old = Matrix(old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode(n, 0)].d_dim, 1);
		d_old = old.continuousStateData.d_k;

		double curMode = q_old;
		std::cout << "Current mode: " << curMode << std::endl;
		std::cout << "Current temp: " << x_old << std::endl;

		if (old.discreteModes > 1)
		{
			// Get discrete transition
			// Number of symbolic variables needed
			int kmax = old.NumSymbols(old);
			std::vector<std::string> x = old.discreteKernel;
			// Generate list of symbols depending on
			// kmax, check type of guards whether simply a function of
			// x or is also of u and d
			GiNaC::lst syms = this->generateListofSymbols(x[kmax]);

			// For each transition from current mode, generate guard
			// Check if take guard or stay
			double index = curMode;
			double updated = 0;
			while (index < pow(old.discreteModes, 2) && old.discreteModes > 1)
			{
				std::cout << "discreteKernel" << x[index] << std::endl;
				bool guard = this->getCurrentGuard(x[index], x_old, syms);
				if (guard && !updated)
				{
					// q_new = std::floor(index/old.discreteModes);
					std::cout << "prob:" << old.Pk[0] << std::endl;
					std::cout << "q-old:" << q_old << std::endl;

					std::vector<int> q = t_q(old, q_old, 0);
					std::cout << "q" << q[0] << " " << q[1] << std::endl;
					q_new = q[1];
					std::cout << "q_new: " << q_new
						<< std::endl;						//<< "index+1: "<< index+1
					updatepk(old, index);
					updated = 1;
				}
				index += old.discreteModes;
			}
		}
		if (old.stateSpaceModelsPerMode[(int)q_old].d_dim == 0)
		{
			x_new =
				old.stateSpaceModelsPerMode[(int)q_new].updateLTI(old.stateSpaceModelsPerMode[(int)q_old], x_old, u_old);
			// TODO Generalise for other models
			if (((q_old == 0 && q_new == 2) || (q_old == 1 && q_new == 0) ||
				(q_old == 2 && q_new == 0)) &&
				cond)
			{
				x_new(1, 0) = 0;
			}
		}
		else
		{
			x_new = old.stateSpaceModelsPerMode[(int)q_new].updateLTIad(old.stateSpaceModelsPerMode[(int)q_old], x_old,
				u_old, d_old);
		}
		old.continuousStateData.x_k = x_new;
		// Append result to continuousVariables in old object (columns)
		old.continuousStateData.continuousVariables.push_back(x_new);
		// Append result to discreteModes in old object (columns)
		old.initialDiscreteMode(n + 1, 0) = q_new;

	}

	void step_hyb_ad(shs_t &old, int n, int cond, int steps)
	{
		Matrix x_old(old.continuousStateData.continuousVariables[n]);
		Matrix x_new = Matrix(old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode(n, 0)].x_dim, steps);
		double q_old = old.initialDiscreteMode(n, 0);
		double q_new = q_old;
		Matrix u_old = Matrix(old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode(n, 0)].u_dim, 1);
		u_old = old.continuousStateData.u_k;
		Matrix d_old = Matrix(old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode(n, 0)].d_dim, 1);
		d_old = old.continuousStateData.d_k;

		if (old.discreteModes == 1)
		{
			x_new = old.stateSpaceModelsPerMode[0].updateLTIst(old.stateSpaceModelsPerMode[0], x_old, u_old, d_old);
		}
		else
		{
			for (int j = 0; j < steps; j++)
			{
				// Get discrete transition
				// Number of symbolic variables needed
				int kmax = old.NumSymbols(old);
				std::vector<std::string> x = old.discreteKernel;
				// Generate list of symbols depending on
				// kmax, check type of guards whether simply a function of
				// x or is also of u and d
				GiNaC::lst syms = this->generateListofSymbols(x[kmax]);

				// For each transition from current mode, generate guard
				// Check if take guard or stay
				double index = q_old;

				double updated = 0;
				while (index < pow(old.discreteModes, 2) && old.discreteModes > 1)
				{
					//  std::cout << "discreteKernel" << x[index] <<std::endl;
					bool guard = this->getCurrentGuard(x[index], x_old, syms);
					if (guard && !updated)
					{
						std::vector<int> q = t_q(old, q_old, j);
						q_new = q[1];
						updatepk(old, index, j);
						updated = 1;
					}
					index += old.discreteModes;
				}
				// Append result to discreteModes in old object (columns)
				old.initialDiscreteMode(n + 1, j) = q_new;
				x_new.col(j) = old.stateSpaceModelsPerMode[(int)q_new].updateLTIst(
					old.stateSpaceModelsPerMode[(int)q_old], x_old.col(j), u_old, d_old);
			}
		}
		// Append result to continuousVariables in old object (columns)
		old.continuousStateData.continuousVariables.push_back(x_new);

	}

	void run(shs_t &old, int N, int cond, int monte)
	{
		// Start simulation timers
		clock_t begin, end;
		begin = clock();
		double time = 0;
		// For deterministic case perform 1 run
		// For stochastic version perform Monte Carlo + compute mean
		int i = 0;
		while (i < N)
		{
			if (i == 0)
			{
				old.obtainpk(*this);
			}
			double a = old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode(i, 0)].sigma.n_rows;
			double b = old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode(i, 0)].sigma.n_cols;
			if ((a + b) > 3)
			{
				old.step_hyb_ad(old, i, cond, monte);
			}
			else
			{
				old.step_hyb(old, i, cond);
			}
			if ((unsigned)old.continuousStateData.actionVariables.n_cols == (unsigned)N)
				old.continuousStateData.u_k = old.continuousStateData.actionVariables.row(i);
			if ((unsigned)old.continuousStateData.disturbanceVariables.n_cols == (unsigned)N)
				old.continuousStateData.d_k = old.continuousStateData.disturbanceVariables.row(i);
			i++;
		}
		Matrix y = old.continuousStateData.continuousVariables[0].t();

		end = clock();
		time = (double)(end - begin) / CLOCKS_PER_SEC;
		std::ostringstream oss;
		auto t = std::time(nullptr);
		auto tm = *std::localtime(&t);

		oss << std::put_time(&tm, "%d-%m-%Y-%H-%M-%S");
		auto str = oss.str();

		std::cout << std::endl;
		std::cout << "--------------------------------------" <<std::endl;
		std::cout << " Simulation time                      " << std::endl;
		std::cout << "--------------------------------------" <<std::endl;
		std::cout << " " << time << std::endl;
		std::cout << "--------------------------------------" <<std::endl;
		std::cout << std::endl;

		// Option to export to file results
		std::ofstream myfile;
		std::string exportOpt;
		std::string str0("y");
		std::string str1("yes");
		std::cout << "Would you like to store simulation results [y- yes, n - no] " << std::endl;
		std::cin >> exportOpt;
		if ((exportOpt.compare(str0) == 0) || (exportOpt.compare(str1) == 0))
		{
			// check if results folder exists:
			if(checkFolderExists("../results") == -1)
			{
				if(mkdir("../results", 0777) == -1)
				{
					std::cerr << "Error cannot create results directory: " <<std::strerror(errno) <<std::endl;
					exit(0);
				}
			}
			// Store results in file
			std::string f_name = "../results/Simulationtime_" + str + ".txt";
			myfile.open(f_name);
			myfile << time <<std::endl;
			myfile.close();
			std::string y_name = "../results/y_" + str + ".txt";
			y.save(y_name, arma::raw_ascii);
			Matrix q = old.initialDiscreteMode;
			std::string q_name = "../results/modes_" + str + ".txt";
			q.save(q_name, arma::raw_ascii);
			// TODO plotting
			//      old.createPySimPlots(y, modes, x_dim);
		}

	}
	bool obtainTqfromMat(const char *fn, shs_t &init)
	{
		// Reading model file input in .Matrix format
		// and storing into ssmodel class
		mat_t *matf;
		matvar_t *matvar, *contents;
		bool error = 0;
		// Read .Matrix file
		matf = Mat_Open(fn, MAT_ACC_RDONLY);
		if (matf)											// if successful in reading file
		{
			// read each variable within file and populate
			// state space model based on variable name
			contents = Mat_VarRead(matf, "discreteKernel");
			if (contents == NULL)
			{
				std::cout << "Variable discreteKernel not found in file" << std::endl;
				std::cout << "Number of modes set to 1" << std::endl;
				init.discreteModes = 1;
			}
			else
			{
				init.populateTq(*contents);
				// Mat_VarFree(matvar);
				matvar = NULL;
				contents = NULL;
			}
			Mat_Close(matf);

		} else												// unsuccessfull in opening file
		{
			return 1;										// throw "Error opening mat file";
		}
	}
	void obtainpk(shs_t &init)
	{
		// Traverse discreteKernel and get equialent probability matrix
		double index = 0;
		std::vector<std::string> str = init.discreteKernel;

		Matrix p = Matrix(init.discreteModes, init.discreteModes);
		for (int j = 0; j < init.discreteModes; j++)
		{
			for (int k = 0; k < init.discreteModes; k++)
			{
				std::vector<std::string> spl = splitStr(str[index], ':');
				if (spl.size() == 1)
				{
					p(k, j) = 0;

				}
				else
				{
					p(k, j) = stod(spl[0]);
					init.discreteKernel[index] = spl[1];
				}
				index += 1;
			}
		}
		p = checkStochasticity(p);
		std::cout << p << std::endl;
		for (unsigned int i = 0; i < init.Pk.size(); i++)
		{
			init.Pk[i] = p;
		}
	}
	void updatepk(shs_t &init, int currentGuardindex)
	{
		// Traverse discreteKernel and get equialent probability matrix
		std::vector<std::string> str = init.discreteKernel;
		Matrix p = init.Pk[0];
		int row = currentGuardindex % init.discreteModes;
		std::cout << "row : " << row << std::endl;
		int col = (int)currentGuardindex / init.discreteModes;
		std::cout << "col: " << col << std::endl;
		int found1 = -1;
		for (unsigned j = 0; j < p.n_cols; j++)
		{
			if (p(row, j) == 1)
			{
				found1 = j;
			}
		}
		if (found1 != -1)
		{
			p(row, found1) = 0;
			p(row, col) = 1;
		}

		std::cout << "p: " << p << std::endl;
		p = checkStochasticity(p);
		std::cout << p << std::endl;

		init.Pk[0] = p;
	}
	void updatepk(shs_t &init, int currentGuardindex, int step)
	{
		// Traverse discreteKernel and get equivalent probability matrix
		std::vector<std::string> str = init.discreteKernel;
		Matrix p = init.Pk[step];
		int row = currentGuardindex % init.discreteModes;
		std::cout << "row : " << row << std::endl;
		int col = (int)(currentGuardindex) / init.discreteModes;
		std::cout << "col: " << col << std::endl;
		std::cout << "p: " << p << std::endl;

		int found1 = -1;
		for (unsigned j = 0; j < p.n_cols; j++)
		{
			if (p(row, j) == 1)
			{
				found1 = j;
			}
		}
		if (found1 != -1)
		{
			p(row, found1) = 0;
			p(row, col) = 1;
		}

		std::cout << "p: " << p << std::endl;
		p = checkStochasticity(p);
		std::cout << p << std::endl;

		init.Pk[step] = p;
	}

	private:
		GiNaC::symbol &get_symbol(const std::string &s)
		{
			static std::map<std::string, GiNaC::symbol> directory;
			std::map<std::string, GiNaC::symbol>::iterator i = directory.find(s);
			if (i != directory.end())
			{
				return i->second;
			}
			else
			{
				return directory.insert(make_pair(s, GiNaC::symbol(s))).first->second;
			}
		}
		GiNaC::lst generateListofSymbols(std::string str)
		{
			GiNaC::lst symbols = {};

			// Check if there are 'x','u' or 'd' characters
			// If there are compute locations
			std::vector<int> x_symb = findLocation(str, 'x');
			std::vector<int> u_symb = findLocation(str, 'u');
			std::vector<int> d_symb = findLocation(str, 'd');
			// Convert strings to symbols as needed
			if (x_symb.size() > 0)
			{
				for (unsigned int i = 0; i < x_symb.size(); i++)
				{
					std::string str2 = str;
					if (str2.find('&') | str2.find('|'))
					{
						std::vector<std::string> x = splitStr(str2, '&');
						std::vector<std::string> y = splitStr(str2, '|');
						if (x.size() != str2.size())
						{
							for (unsigned int j = 0; j < x.size(); j++)
							{
								//	   ;
								// 	   std::cout << x[j].substr(0,2)  << std::endl;
								GiNaC::symbol t = get_symbol(x[j].substr(0, 2));
								symbols.append(t);
							}
						}
						if (y.size() >= str.size() - 1)
						{
							for (unsigned int m = 0; m < y.size(); m++)
							{

								GiNaC::symbol t = get_symbol(y[m].substr(0, 2));
								symbols.append(t);
							}
						}
					}
					else
					{
						str2.substr(x_symb[i], 2);
						GiNaC::symbol t = get_symbol(str2);
						symbols.append(t);
					}
				}
			}
			if (u_symb.size() > 0)
			{
				int k = 0;
				for (unsigned int i = x_symb.size() - 1;
					i < x_symb.size() + u_symb.size() - 1; i++)
				{
					std::string str2 = str;					//.substr(u_symb[k],str.size());
					if (str2.find('&') | str2.find('|'))
					{
						std::vector<std::string> x = splitStr(str2, '&');
						std::vector<std::string> y = splitStr(str2, '|');
						for (unsigned int j = 0; j < x.size(); j++)
						{

							GiNaC::symbol t = get_symbol(x[j].substr(0, 2));
							symbols.append(t);
						}
						for (unsigned int k = 0; k < y.size(); k++)
						{

							GiNaC::symbol t = get_symbol(y[k].substr(0, 2));
							symbols.append(t);
						}
					}
					else
					{
						str2.substr(u_symb[k], 2);
						GiNaC::symbol t = get_symbol(str2);
						symbols.append(t);
					}
					k++;
				}
			}

			if (d_symb.size() > 0)
			{
				int l = 0;
				for (unsigned int i = x_symb.size() + u_symb.size() - 1;
					i < x_symb.size() + u_symb.size() + d_symb.size() - 1; i++)
				{
					std::string str2 = str;					//.substr(d_symb[l],str.size());
					if (str2.find('&') | str2.find('|'))
					{
						std::vector<std::string> x = splitStr(str2, '&');
						std::vector<std::string> y = splitStr(str2, '|');
						for (unsigned int j = 0; j < x.size(); j++)
						{

							GiNaC::symbol t = get_symbol(x[j].substr(0, 2));
							symbols.append(t);
						}
						for (unsigned int k = 0; k < y.size(); k++)
						{

							GiNaC::symbol t = get_symbol(y[k].substr(0, 2));
							symbols.append(t);
						}
					}
					else
					{
						str2.substr(d_symb[l], 2);
						GiNaC::symbol t = get_symbol(str2);
						symbols.append(t);
					}
					l++;
				}
			}
			return symbols;
		}
		// Function to determine the type of guard
		// and output the lhs and rhs conditions
		//
		bool getCurrentGuard(std::string x, Matrix x_old, GiNaC::lst syms)
		{
			std::vector<bool> guard = {};
			std::vector<std::string> y;
			std::vector<std::string> y0;					/// = splitStr(x,'&');
			std::vector<std::string> y01;					// = splitStr(x,'|');
			bool logAND = false;
			bool logOR = false;
			// Determine type of guard whether it contains:
			// ' ', '0', '>','<','<=','<=' or a combination
			if (x.size() == 1)								// No action
			{
				guard.push_back(0);							// stay same

			}
			else
			{
				std::vector<std::string> y = splitStr(x, ':');
				std::vector<std::string> y0 = splitStr(x, '&');
				std::vector<std::string> y01 = splitStr(x, '|');
				std::vector<std::string> xnew;
				int steps = 2;								// TODO make a function of number of symbols
				if (y0.size() > 1 && y01.size() == 1)
				{
					steps = y0.size();
					xnew = y0;
					logAND = true;
				}
				else if (y0.size() == 1 && y01.size() > 1)
				{
					steps = y01.size();
					xnew = y01;
					logAND = true;
				}
				else
				{
					steps = 1;								// y0.size()+y01.size();
					xnew = y0;
					//  xnew.push_back(y01[0]);//TODO FIX
					logAND = false;
					logOR = false;
				}
				for (int j = 0; j < steps; j++)
				{
					std::cout << "xold" << x_old(j, 0) << std::endl;

					std::vector<std::string> y1 = splitStr(xnew[j], '=');
					std::vector<std::string> y2 = splitStr(xnew[j], '>');
					std::vector<std::string> y3 = splitStr(xnew[j], '<');
					int totalSize = 0;
					if (y3[0].length() == x.length() && y2[0].length() != x.length() &&
						y1[0].length() != x.length())
						totalSize = y1.size() + y2.size();
					else if (y3[0].length() != x.length() && y2[0].length() == x.length() &&
						y1[0].length() != x.length())
						totalSize = y1.size() + y3.size();
					else if (y3[0].length() != x.length() && y2[0].length() != x.length() &&
						y1[0].length() == x.length())
						totalSize = y2.size() + y3.size();
					else
						totalSize = y1.size() + y2.size() + y3.size();
					// std::cout << "Total size: "<< totalSize <<std::endl;
					// Filter strings from logical operators
					for (unsigned int i = 0; i < y1.size(); i++)
					{
						// std::cout<< "Original expression: " << y1[i] << std::endl;
						int a = y1[i].find("x"), b = y1[i].find(">"), c = y1[i].find("<"),
							l = y1[i].length();

						if (a > -1 && (b > -1 || c > -1) && l > 3)
						{
							if (y1[i].find("x") == 0)
							{
								y1[i].erase(0, 3);
							}
							else
							{
								if (b > -1 && c > -1)
								{
									y1[i].erase(y1[i].find('>'), y1[i].length() - y1[i].find('>'));
									y1[i].erase(y1[i].find('<'), y1[i].length() - y1[i].find('<'));
								}
								else if (b == -1 && c > -1)
								{
									y1[i].erase(y1[i].find('<'), y1[i].length() - y1[i].find('<'));
								}
								else if (b > -1 && c == -1)
								{
									y1[i].erase(y1[i].find('>'), y1[i].length() - y1[i].find('>'));
								}
							}
						}
						else
						{
							y1[i].erase(std::remove(y1[i].begin(), y1[i].end(), '>'),
								y1[i].end());
							y1[i].erase(std::remove(y1[i].begin(), y1[i].end(), '<'),
								y1[i].end());
						}
						// std::cout<< "Clean expression: " << y1[i] << std::endl;
					}
					for (unsigned int i = 0; i < y2.size(); i++)
					{
						// std::cout<< "Original expression: " << y2[i] << std::endl;
						int index = y2[i].find("<=");
						int index2 = y2[i].find("x");
						if (index > -1)
						{
							if (index2 > index)
							{
								y2[i].erase(index, y2[i].length() - index);
							}
							else
							{
								y2[i].erase(0, index + 2);
							}
						}
						y2[i].erase(std::remove(y2[i].begin(), y2[i].end(), '='),
							y2[i].end());
						// std::cout<< "Clean expression: " << y2[i] << std::endl;
					}
					for (unsigned int i = 0; i < y3.size(); i++)
					{
						// std::cout<< "Original expression: " << y3[i] << std::endl;
						int index = y3[i].find(">=");
						int index2 = y3[i].find("x");
						if (index > -1)
						{
							if (index2 > index)
							{
								y3[i].erase(index, y3[i].length() - index);
							}
							else
							{
								y3[i].erase(0, index + 2);
							}
						}
						y3[i].erase(std::remove(y3[i].begin(), y3[i].end(), '='),
							y3[i].end());

						// std::cout<< "Clean expression: " << y3[i] << std::endl;
					}
					// To compute guard
					std::cout << syms[0] << std::endl;
					bool x1exist = y1[0].find("x1") != std::string::npos;
					bool x2exist = y1[0].find("x2") != std::string::npos;
					bool x1y2exist = y2[0].find("x1") != std::string::npos;
					bool x2y2exist = y2[0].find("x2") != std::string::npos;
					bool x1y3exist = y3[0].find("x1") != std::string::npos;
					bool x2y3exist = y3[0].find("x2") != std::string::npos;
					std::cout << x1exist << x2exist << x1y2exist << x2y2exist << x1y3exist
						<< x2y3exist << std::endl;
					switch (totalSize)
					{
						case 2:								// = or > or <
						{
							if (y1.size() == 2)
							{
								GiNaC::ex ex1(y1[0], syms);
								GiNaC::ex ex2(y1[1], syms);
								GiNaC::ex lhs;
								GiNaC::ex rhs;
								if (x1exist)
								{
									lhs = ex1.subs(syms[0] == x_old(0, 0));
									rhs = ex2.subs(syms[0] == x_old(0, 0));
								}
								else
								{
									lhs = ex1.subs(syms[1] == x_old(1, 0));
									rhs = ex2.subs(syms[1] == x_old(1, 0));
								}
								guard.push_back((lhs == rhs));
								break;
							}
							if (y2.size() == 2)
							{
								GiNaC::ex ex1(y2[0], syms);
								GiNaC::ex ex2(y2[1], syms);
								GiNaC::ex lhs;
								GiNaC::ex rhs;
								if (x1y2exist)
								{
									lhs = ex1.subs(syms[0] == x_old(0, 0));
									rhs = ex2.subs(syms[0] == x_old(0, 0));
								}
								else
								{
									lhs = ex1.subs(syms[1] == x_old(1, 0));
									rhs = ex2.subs(syms[1] == x_old(1, 0));
								}
								guard.push_back(lhs > rhs);
								break;
							}
							if (y3.size() == 2)
							{
								GiNaC::ex ex1(y3[0], syms);
								GiNaC::ex ex2(y3[1], syms);
								GiNaC::ex lhs;
								GiNaC::ex rhs;
								if (x1y3exist)
								{
									lhs = ex1.subs(syms[0] == x_old(0, 0));
									rhs = ex2.subs(syms[0] == x_old(0, 0));
								}
								else
								{
									lhs = ex1.subs(syms[1] == x_old(1, 0));
									rhs = ex2.subs(syms[1] == x_old(1, 0));
								}
								guard.push_back(lhs < rhs);
								break;
							}
						} break;
						case 3:								// a<x<b, a>x>b
							if (y2.size() == 3)
							{
								GiNaC::ex ex1(y2[0], syms);
								GiNaC::ex ex2(y2[1], syms);
								GiNaC::ex ex3(y2[2], syms);
								GiNaC::ex lhs;
								GiNaC::ex rhs;
								GiNaC::ex b;
								if (x1y2exist)
								{
									lhs = ex1.subs(syms[0] == x_old(0, 0));
									rhs = ex2.subs(syms[0] == x_old(0, 0));
									b = ex3.subs(syms[0] == x_old(0, 0));
								}
								else
								{
									lhs = ex1.subs(syms[1] == x_old(1, 0));
									rhs = ex2.subs(syms[1] == x_old(1, 0));
									b = ex3.subs(syms[1] == x_old(1, 0));
								}

								guard.push_back((lhs > rhs) || (rhs > b));
								break;
							}
							if (y3.size() == 3)
							{
								GiNaC::ex ex1(y3[0], syms);
								GiNaC::ex ex2(y3[1], syms);
								GiNaC::ex ex3(y3[2], syms);
								GiNaC::ex lhs;
								GiNaC::ex rhs;
								GiNaC::ex b;
								if (x1y3exist)
								{
									lhs = ex1.subs(syms[0] == x_old(0, 0));
									rhs = ex2.subs(syms[0] == x_old(0, 0));
									b = ex3.subs(syms[0] == x_old(0, 0));
								}
								else
								{
									lhs = ex1.subs(syms[1] == x_old(1, 0));
									rhs = ex2.subs(syms[1] == x_old(1, 0));
									b = ex3.subs(syms[1] == x_old(1, 0));
								}
								guard.push_back((lhs < rhs) || (rhs < b));
								break;
							}
							break;
						case 4:								// a>=x, a<=x,a>x<b, a<x>b
							if (y1.size() == 2 && y2.size() == 2)
							{
								GiNaC::ex ex1(y1[0], syms);
								GiNaC::ex ex2(y1[1], syms);
								GiNaC::ex lhs;
								GiNaC::ex rhs;
								if (x1exist)
								{
									lhs = ex1.subs(syms[0] == x_old(0, 0));
									rhs = ex2.subs(syms[0] == x_old(0, 0));
								}
								else
								{
									lhs = ex1.subs(syms[1] == x_old(1, 0));
									rhs = ex2.subs(syms[1] == x_old(1, 0));
								}
								int where = x.find('x');
								int mid = x.length() / 2;

								std::cout << lhs << " " << rhs << std::endl;
								if (mid - where > 0)
									guard.push_back((lhs == rhs) || (lhs > rhs));
								else
									guard.push_back((rhs == rhs) || (rhs > lhs));
								break;
							}
							if (y1.size() == 2 && y3.size() == 2)
							{
								GiNaC::ex ex1(y1[0], syms);
								GiNaC::ex ex2(y1[1], syms);
								GiNaC::ex lhs;
								GiNaC::ex rhs;
								if (x1exist)
								{
									lhs = ex1.subs(syms[0] == x_old(0, 0));
									rhs = ex2.subs(syms[0] == x_old(0, 0));
								}
								else
								{
									std::cout << x_old(1, 0) << std::endl;
									lhs = ex1.subs(syms[1] == x_old(1, 0));
									rhs = ex2.subs(syms[1] == x_old(1, 0));
								}
								int where = x.find('x');
								int mid = x.length() / 2;
								bool one = lhs == rhs;
								bool two = (lhs < rhs);
								std::cout << one << " " << two << std::endl;
								if (mid - where > 0)
									guard.push_back((lhs == rhs) || (lhs < rhs));
								else
									guard.push_back((rhs == rhs) || (rhs < lhs));
								break;
							}
							if (y2.size() == 2 && y3.size() == 2)
							{
								GiNaC::ex ex1(y2[0], syms);
								GiNaC::ex ex2(y2[1], syms);
								GiNaC::ex ex3(y3[0], syms);
								GiNaC::ex ex4(y3[1], syms);

								GiNaC::ex lhs, b, c;
								GiNaC::ex rhs;
								if (x1y2exist)
								{
									lhs = ex1.subs(syms[0] == x_old(0, 0));
									rhs = ex2.subs(syms[0] == x_old(0, 0));
								}
								else
								{
									lhs = ex1.subs(syms[1] == x_old(1, 0));
									rhs = ex2.subs(syms[1] == x_old(1, 0));
								}
								if (x1y3exist)
								{
									b = ex3.subs(syms[0] == x_old(0, 0));
									c = ex4.subs(syms[0] == x_old(0, 0));
								}
								else
								{
									b = ex3.subs(syms[1] == x_old(1, 0));
									c = ex4.subs(syms[1] == x_old(1, 0));
								}
								guard.push_back((lhs > rhs) || (b < c));
								break;
							}
							break;
						case 5:								// a>=x>b, a<=x<b,
							if (y1.size() == 2 && y2.size() == 1 && y3.size() == 2)
							{
								GiNaC::ex ex1(y1[0], syms);
								GiNaC::ex ex2(y1[1], syms);
								GiNaC::ex ex3(y3[0], syms);
								GiNaC::ex ex4(y3[1], syms);
								GiNaC::ex lhs, b, c;
								GiNaC::ex rhs;
								if (x1exist)
								{
									lhs = ex1.subs(syms[0] == x_old(0, 0));
									rhs = ex2.subs(syms[0] == x_old(0, 0));
								}
								else
								{
									lhs = ex1.subs(syms[1] == x_old(1, 0));
									rhs = ex2.subs(syms[1] == x_old(1, 0));
								}
								if (x1y3exist)
								{
									b = ex3.subs(syms[0] == x_old(0, 0));
									c = ex4.subs(syms[0] == x_old(0, 0));
								}
								else
								{
									b = ex3.subs(syms[1] == x_old(1, 0));
									c = ex4.subs(syms[1] == x_old(1, 0));
								}
								bool a = (lhs == rhs), d = (b < c);
								std::cout << "==" << a << " <" << d << std::endl;
								guard.push_back((lhs == rhs) || (b < c));
							}
							if (y1.size() == 2 && y2.size() == 2 && y3.size() == 1)
							{
								GiNaC::ex ex1(y1[0], syms);
								GiNaC::ex ex2(y1[1], syms);
								GiNaC::ex ex3(y2[0], syms);
								GiNaC::ex ex4(y2[1], syms);
								GiNaC::ex lhs, b, c;
								GiNaC::ex rhs;
								if (x1exist)
								{
									lhs = ex1.subs(syms[0] == x_old(0, 0));
									rhs = ex2.subs(syms[0] == x_old(0, 0));
								}
								else
								{
									lhs = ex1.subs(syms[1] == x_old(1, 0));
									rhs = ex2.subs(syms[1] == x_old(1, 0));
								}
								if (x1y2exist)
								{
									b = ex3.subs(syms[0] == x_old(0, 0));
									c = ex4.subs(syms[0] == x_old(0, 0));
								}
								else
								{
									b = ex3.subs(syms[1] == x_old(1, 0));
									c = ex4.subs(syms[1] == x_old(1, 0));
								}
								guard.push_back((lhs == rhs) || (b > c));
							}
							if (y1.size() == 2 && y2.size() == 3)
							{
								GiNaC::ex ex1(y1[0], syms);
								GiNaC::ex ex2(y1[1], syms);
								GiNaC::ex ex3(y2[0], syms);
								GiNaC::ex ex4(y2[1], syms);
								GiNaC::ex ex5(y2[2], syms);
								GiNaC::ex lhs, b, c;
								GiNaC::ex rhs, d;
								if (x1exist)
								{
									lhs = ex1.subs(syms[0] == x_old(0, 0));
									rhs = ex2.subs(syms[0] == x_old(0, 0));
								}
								else
								{
									lhs = ex1.subs(syms[1] == x_old(1, 0));
									rhs = ex2.subs(syms[1] == x_old(1, 0));
								}
								if (x1y2exist)
								{
									b = ex3.subs(syms[0] == x_old(0, 0));
									c = ex5.subs(syms[0] == x_old(0, 0));
									d = ex4.subs(syms[0] == x_old(0, 0));
								}
								else
								{
									b = ex3.subs(syms[1] == x_old(1, 0));
									c = ex5.subs(syms[1] == x_old(1, 0));
									d = ex4.subs(syms[1] == x_old(1, 0));
								}
								int where = x.find('=');
								int mid = x.length() / 2;
								if (mid - where > 0)
									guard.push_back((lhs == c) || (b > c) || (c > d));
								else
									guard.push_back((rhs == c) || (b > c) || (c > d));
								break;
							}
							if (y1.size() == 2 && y3.size() == 3)
							{
								GiNaC::ex ex1(y1[0], syms);
								GiNaC::ex ex2(y1[1], syms);
								GiNaC::ex ex3(y3[0], syms);
								GiNaC::ex ex4(y3[1], syms);
								GiNaC::ex ex5(y3[2], syms);
								GiNaC::ex lhs, b, c;
								GiNaC::ex rhs, d;
								if (x1exist)
								{
									lhs = ex1.subs(syms[0] == x_old(0, 0));
									rhs = ex2.subs(syms[0] == x_old(0, 0));
								}
								else
								{
									lhs = ex1.subs(syms[1] == x_old(1, 0));
									rhs = ex2.subs(syms[1] == x_old(1, 0));
								}
								if (x1y3exist)
								{
									b = ex3.subs(syms[0] == x_old(0, 0));
									c = ex5.subs(syms[0] == x_old(0, 0));
									d = ex4.subs(syms[0] == x_old(0, 0));
								}
								else
								{
									b = ex3.subs(syms[1] == x_old(1, 0));
									c = ex5.subs(syms[1] == x_old(1, 0));
									d = ex4.subs(syms[1] == x_old(1, 0));
								}
								int where = x.find('=');
								int mid = x.length() / 2;
								if (mid - where > 0)
									guard.push_back((lhs == c) || (b < c) || (c < d));
								else
									guard.push_back((rhs == c) || (b < c) || (c < d));
								break;
							}
							break;
						case 6:								// a<=x>b, a>x<=b,a>=x<b,a<x>=b
						{
							GiNaC::ex ex1(y1[0], syms);
							GiNaC::ex ex2(y1[1], syms);
							GiNaC::ex ex3(y2[0], syms);
							GiNaC::ex ex4(y2[1], syms);
							GiNaC::ex ex5(y3[0], syms);
							GiNaC::ex ex6(y3[1], syms);
							GiNaC::ex lhs, b, c;
							GiNaC::ex rhs, d, e;
							if (x1exist)
							{
								lhs = ex1.subs(syms[0] == x_old(0, 0));
								rhs = ex2.subs(syms[0] == x_old(0, 0));
							}
							else
							{
								lhs = ex1.subs(syms[1] == x_old(1, 0));
								rhs = ex2.subs(syms[1] == x_old(1, 0));
							}
							if (x1y2exist)
							{
								b = ex3.subs(syms[0] == x_old(0, 0));
								c = ex4.subs(syms[0] == x_old(0, 0));
							}
							else
							{
								b = ex3.subs(syms[1] == x_old(1, 0));
								c = ex4.subs(syms[1] == x_old(1, 0));
							}
							if (x1y3exist)
							{
								d = ex5.subs(syms[0] == x_old(0, 0));
								e = ex6.subs(syms[0] == x_old(0, 0));
							}
							else
							{
								d = ex5.subs(syms[1] == x_old(1, 0));
								e = ex6.subs(syms[1] == x_old(1, 0));
							}
							int where = x.find('=');
							int mid = x.length() / 2;
							if (mid - where > 0)
								guard.push_back((lhs == c) || (b > c) || (c < d));
							else
								guard.push_back((rhs == c) || (b > c) || (c < d));
							break;
						}

						case 7:								// a>=x<=b, a<=x>=b
						{
							GiNaC::ex ex1(y1[0], syms);
							GiNaC::ex ex2(y1[1], syms);
							GiNaC::ex ex2a(y1[2], syms);
							GiNaC::ex ex3(y2[0], syms);
							GiNaC::ex ex4(y2[1], syms);
							GiNaC::ex ex5(y3[0], syms);
							GiNaC::ex ex6(y3[1], syms);
							GiNaC::ex lhs, a, b, c;
							GiNaC::ex rhs, d, e;
							if (x1exist)
							{
								lhs = ex1.subs(syms[0] == x_old(0, 0));
								rhs = ex2.subs(syms[0] == x_old(0, 0));
								a = ex2a.subs(syms[0] == x_old(0, 0));
							}
							else
							{
								lhs = ex1.subs(syms[1] == x_old(1, 0));
								rhs = ex2.subs(syms[1] == x_old(1, 0));
								a = ex2a.subs(syms[1] == x_old(0, 0));
							}
							if (x1y2exist)
							{
								b = ex3.subs(syms[0] == x_old(0, 0));
								c = ex4.subs(syms[0] == x_old(0, 0));
							}
							else
							{
								b = ex3.subs(syms[1] == x_old(1, 0));
								c = ex4.subs(syms[1] == x_old(1, 0));
							}
							if (x1y3exist)
							{
								d = ex5.subs(syms[0] == x_old(0, 0));
								e = ex6.subs(syms[0] == x_old(0, 0));
							}
							else
							{
								d = ex5.subs(syms[1] == x_old(1, 0));
								e = ex6.subs(syms[1] == x_old(1, 0));
							}
							guard.push_back((lhs == rhs) || (rhs == a) || (b > c) || (c < d));
							break;
						}
						default:
						{
							guard.push_back(0);
							break;
						}
					}
				}
			}
			bool answer;
			std::cout << guard[0] << ", " << guard[1] << std::endl;
			if (logAND && logOR)							// TODO FIX AND GENERALISE
			{
				answer = guard[0] && guard[1];
			}
			else if (logAND && !logOR)
			{
				answer = guard[0] && guard[1];
			}
			// else if(!logAND && logOR)
			// {
			else
			{
				answer = guard[0] || guard[1];
			}
			//}
			// else
			//{
			// answer = guard[0];
			// }

			return answer;
		}

		int checkData()
		{
			int error = 0;
			if ((unsigned)this->stateSpaceModelsPerMode[0].x_dim != (unsigned)this->continuousStateData.continuousVariables[0].n_rows)
			{
				error = 1;
			}
			if ((unsigned)this->stateSpaceModelsPerMode[0].u_dim != (unsigned)this->continuousStateData.actionVariables.n_cols &&
				this->stateSpaceModelsPerMode[0].u_dim > 0)
			{
				error = 2;
			}
			if ((unsigned)this->stateSpaceModelsPerMode[0].d_dim != (unsigned)this->continuousStateData.disturbanceVariables.n_cols &&
				this->stateSpaceModelsPerMode[0].d_dim > 0)
			{
				error = 3;
			}
			return error;
		}
		int NumSymbols(shs_t &old)
		{
			// To find how many symbolic continuous variables are needed
			// find string with maximum length
			int kmax = 0;
			for (int k = 0; k < pow(old.discreteModes, 2); k++)
			{
				if (sizeof(old.discreteKernel[k]) > sizeof(old.discreteKernel[kmax]))
				{
					kmax = k;
				}
			}
			// std::cout << "kmax: " << kmax << std::endl;
			return kmax;
		}
		void populateTq(matvar_t &content)
		{
			ssmodels_t container;

			// discreteKernel can be of two version
			// Cells containing strings for guards or
			// Numeric containing probabilities
			if (content.data != NULL)
			{
				std::cout << "discreteKernel: " << content.data_type << std::endl;
				// Reading from cells
				std::string str = container.readCells(content);
				std::cout << "discreteKernel: " << str << std::endl;
															// TODO: CHANGEseparater
				std::vector<std::string> x = splitStr(str, ' ');
				int numEl = x.size();
				std::cout << numEl << std::endl;
				int q = sqrt(numEl);
				std::cout << "modes: " << q << std::endl;
				this->discreteModes = q;
				this->discreteKernel = x;					// Recall being stored in column format
			}
			else
			{
				std::cout << "discreteKernel field not input in Matrix file" << std::endl;
			}
		}
		std::vector<int> t_q(shs_t &old, int q_old, int monte)
		{
			int steps = 1;
			std::vector<int> modes(steps + 1);

			modes[0] = q_old;
			// if(old.stateSpaceModelsPerMode[0].sigma.isZero(0) )
			//{
			std::random_device rd;
			std::mt19937 gen(rd());
			int count;
			double sum, actionVariables;
			for (int i = 0; i < steps; i++)
			{
				count = 0;
				sum = 0;
				actionVariables = std::generate_canonical<double, 10>(gen);
				while (sum < actionVariables)
				{
					sum += old.Pk[monte](modes[i], count);
					if (sum > actionVariables)
					{
						modes[i + 1] = count;
					}
					count++;
				}
			}

			return modes;
		}
};

/*************************************************************************************************/
// Single initial state SHS with fixed or sigmoidal probabilities governing the
// transitioning between modes`
template <> class shs_t<Matrix, int>
{
	public:
		int discreteModes;									// Discrete modes
		Matrix discreteKernel;								// Discrete kernel- can represent either guards or probabilities
		int continuousSpaceDimension						// Dimension of continuous state space in given mode
			Matrix Pk;
		Matrix initialDiscreteMode;							// Modes
		std::vector<ssmodels_t> stateSpaceModelsPerMode;	// container of models
		SimulationData continuousStateData;

	private:
		bool sigmoid;

	public:
		shs_t()
		{
			n = 1;											// Dimension of continuous state space in given mode
			Pk = Matrix(discreteModes, n);
			initialDiscreteMode = Matrix(discreteModes, 1);	// Initial mode
			std::vector<ssmodels_t> stateSpaceModelsPerMode;
			discreteModes = 1;								// Discrete modes
			SimulationData continuousStateData;
			discreteKernel = Matrix(n, n);
			sigmoid = false;
		}
		shs_t(int disc, int num, Matrix tq, std::vector<ssmodels_t> &old,
			SimulationData &data)
		{
			discreteModes = disc;							// Discrete modes
			n = num;										// Dimension of continuous state space in given mode
			Pk = data.initialModeDistribution[0];
			initialDiscreteMode = data.initialMode[0];		// Initial mode
			stateSpaceModelsPerMode = old;
			continuousStateData = data;
			discreteKernel = tq;
			sigmoid = false;
			try
			{
				// Check dimensions of data
				int err = this->checkData();
				switch (err)
				{
					case 3:
					{
						throw "Incorrect size of disturbance vector";
					}
					default: { std::cout << "Correct input vectors size" << std::endl; }
				}
			}
			catch (const char *msg)
			{
				std::cerr << msg << std::endl;
				exit(0);
			}
		}
		shs_t(int num, std::vector<ssmodels_t> &old, SimulationData &data)
		{
			discreteModes = data.actionVariables.n_cols;	// Discrete modes
			n = num;										// Dimension of continuous state space in given mode
			Pk = data.initialModeDistribution[0];
			initialDiscreteMode = data.initialMode[0];		// Initial mode
			stateSpaceModelsPerMode = old;
			continuousStateData = data;
			discreteKernel = ZeroMatrix(1, 1);
			discreteKernel.reset();
			sigmoid = false;
			// update actionVariables dimension
			for (size_t i = 0; i < discreteModes; i++)
			{
				old[i].u_dim = data.actionVariables.n_cols;
			}
			try
			{
				// Check dimensions of data
				int err = this->checkData();
				switch (err)
				{
					case 3:
					{
						throw "Incorrect size of disturbance vector";
					}
					default: { std::cout << "Correct input vectors size" << std::endl; }
				}
			}
			catch (const char *msg)
			{
				std::cerr << msg << std::endl;
				exit(0);
			}
		}
		shs_t(int num, std::vector<ssmodels_t> &old)
		{
			discreteModes = old.size();						// Discrete modes
			n = num;										// Dimension of continuous state space in given mode
			Pk = ZeroMatrix(1, 1);
			initialDiscreteMode = ZeroMatrix(1, 1);			// Initial mode
			stateSpaceModelsPerMode = old;
			SimulationData data;
			continuousStateData = data;
			discreteKernel = ZeroMatrix(1, 1);
			discreteKernel.reset();
			sigmoid = false;
		}
		shs_t(std::vector<ssmodels_t> &old)
		{
			discreteModes = old.size();						// Discrete modes
			n = old[0].A.n_cols;							// Dimension of continuous state space in given mode
			Pk = ZeroMatrix(1, 1);
			initialDiscreteMode = ZeroMatrix(1, 1);			// Initial mode
			stateSpaceModelsPerMode = old;
			SimulationData data;
			continuousStateData = data;
			discreteKernel = ZeroMatrix(1, 1);
			discreteKernel.reset();
			sigmoid = false;
		}
		shs_t(Matrix tq, std::vector<ssmodels_t> &old, SimulationData &data)
		{
			discreteModes = tq.n_rows;						// Discrete modes
			n = old[0].x_dim;								// Dimension of continuous state space in given mode
			Pk = data.initialModeDistribution[0];
			initialDiscreteMode = data.initialMode[0];		// Initial mode
			stateSpaceModelsPerMode = old;
			continuousStateData = data;
			discreteKernel = tq;
			discreteKernel = checkStochasticity(discreteKernel);
			sigmoid = false;
			for (size_t i = 0; i < discreteModes; i++)
			{
				old[i].u_dim = data.actionVariables.n_cols;
			}
			try
			{
				// Check dimensions of data
				int err = this->checkData();
				switch (err)
				{
					case 3:
					{
						throw "Incorrect size of disturbance vector";
					}
					default: { std::cout << "Correct input vectors size" << std::endl; }
				}
			}
			catch (const char *msg)
			{
				std::cerr << msg << std::endl;
				exit(0);
			}
		}
		shs_t(Matrix tq, std::vector<ssmodels_t> &old)
		{
			discreteModes = tq.n_rows;						// Discrete modes
			n = old[0].x_dim;								// Dimension of continuous state space in given mode
			Pk = ZeroMatrix(1, 1);
			initialDiscreteMode = ZeroMatrix(1, 1);			// Initial mode
			stateSpaceModelsPerMode = old;
			SimulationData data;
			continuousStateData = data;
			discreteKernel = tq;
			discreteKernel = checkStochasticity(discreteKernel);
			sigmoid = false;
			try
			{
				// Check dimensions of data
				int err = this->checkData();
				switch (err)
				{
					case 3:
					{
						throw "Incorrect size of disturbance vector";
					}
					default: { std::cout << "Correct input vectors size" << std::endl; }
				}
			}
			catch (const char *msg)
			{
				std::cerr << msg << std::endl;
				exit(0);
			}
		}
		shs_t(const char *fn, SimulationData &data)
		{
			// Obtain number of discrete modes and obtain discreteKernel
			bool error = this->obtainTqfromMat(fn, *this);

			if (!error)
			{
				discreteModes = this->discreteModes;
				discreteKernel = this->discreteKernel;
				sigmoid = false;
				continuousStateData = data;
				Pk = data.initialModeDistribution[0];
				initialDiscreteMode = data.initialMode[0];	// Initial mode
				// create array containing the ssmodels
				std::vector<ssmodels_t> models;
				Matrix dummy = -1 * arma::ones<Matrix>(1, 1);
				for (int i = 1; i <= discreteModes; i++)
				{
					int x_dimension = data.x_k.n_rows;
					int u_dimension = data.u_k.n_rows;
					int d_dimension = data.d_k.n_rows;
					if (u_dimension == 1 && d_dimension == 1)
					{
						double lhs = arma::accu(data.u_k == dummy);
						double rhs = arma::accu(data.d_k == dummy);
						if ((lhs > 0) && (rhs > 0))
						{
							u_dimension = 0;
							d_dimension = 0;
						}
						else if ((lhs > 0) && (rhs == 0))
						{
							u_dimension = 0;
						}
					}
					ssmodels_t mod(x_dimension, u_dimension, d_dimension);
					if ((mod.C.n_rows == 1) && (mod.C.n_cols == 1) && x_dimension > 1)
					{
						mod.C.resize(x_dimension, x_dimension);
						mod.C = arma::eye<Matrix>(x_dimension, x_dimension);
					}
															// TODO: update to handle different dimensions
					this->n = x_dimensiocontinuousSpaceDimension
					// for each state (Store)
					// Reading models
						std::cout
						<< "Initialising model of continuous variables in discrete mode "
						<< i << std::endl;
					mod.obtainSSfromMat(fn, mod, i);
					if ((unsigned)x_dimension != (unsigned)mod.x_dim ||
						(unsigned)x_dimension != (unsigned)mod.A.n_rows)
					{
						mod.obtainSSfromMat(fn, mod);
					}
					mod.checkModel(mod);
					models.push_back(mod);
				}
				stateSpaceModelsPerMode = models;
				n = models[0].x_dim;

				// Check dimensions of data
				int err = this->checkData();
				switch (err)
				{
					case 3:
					{
						std::cout << "Incorrect size of disturbance vector" << std::endl;
						exit(0);
					}
					default: { std::cout << "Correct input vectors size" << std::endl; }
				}
			}
			else
			{
				std::cout << "File " << fn << " not found." << std::endl;
				exit(0);
			}
		}
		// generate shs from file
		shs_t(const char *fn)
		{
			// Obtain number of discrete modes and obtain discreteKernel
			bool error = this->obtainTqfromMat(fn, *this);
			if (!error)
			{
				discreteModes = this->discreteModes;
				discreteKernel = this->discreteKernel;
				sigmoid = false;
				// create array containing the ssmodels
				std::vector<ssmodels_t> models;
				Matrix dummy = -1 * arma::ones<Matrix>(1, 1);
				for (int i = 1; i <= discreteModes; i++)
				{

					ssmodels_t mod;
					// Reading models
					std::cout
						<< "Initialising model of continuous variables in discrete mode "
						<< i << std::endl;
					mod.obtainSSfromMat(fn, mod, i);
					this->n = mod.A.n_rows;					// TODO: update to handle different dimensions
					// for each state (Store)

					mod.checkModel(mod);
					models.push_back(mod);
				}
				stateSpaceModelsPerMode = models;
				n = models[0].x_dim;
			}
			else
			{
				std::cout << "File " << fn << " not found." << std::endl;
				exit(0);
			}
		}
		shs_t(const char *fn, int bmdp, int modes)
		{
			if (bmdp == 1)
			{
				std::vector<ssmodels_t> models;
				for (int i = 1; i <= modes; i++)
				{
					ssmodels_t mod(1, 1, 1);
					mod.obtainBMDPfromMat(fn, mod, i);
					mod.x_dim = mod.A.n_rows;
					mod.u_dim = 0;
					mod.d_dim = 0;
					models.push_back(mod);
				}
				this->stateSpaceModelsPerMode = models;
				this->n = models[0].x_dim;
			}
			else
			{
				// Obtain number of discrete modes and obtain discreteKernel
				bool error = this->obtainTqfromMat(fn, *this);
				if (!error)
				{
					discreteModes = this->discreteModes;
					discreteKernel = this->discreteKernel;
					sigmoid = false;
					// create array containing the ssmodels
					std::vector<ssmodels_t> models;
					Matrix dummy = -1 * arma::ones<Matrix>(1, 1);
					for (int i = 1; i <= discreteModes; i++)
					{
						int x_dimension = this->stateSpaceModelsPerMode[i].A.n_rows;
						int u_dimension = this->stateSpaceModelsPerMode[i].B.n_rows;
						int d_dimension = this->stateSpaceModelsPerMode[i].F.n_rows;
						ssmodels_t mod(x_dimension, u_dimension, d_dimension);
															// TODO: update to handle different dimensions
						this->n = x_dimensiocontinuousSpaceDimension
						// for each state (Store)
						// Reading models
							mod.obtainSSfromMat(fn, mod, i);
						if ((unsigned)x_dimension != (unsigned)mod.x_dim ||
							(unsigned)x_dimension != (unsigned)mod.A.n_rows)
						{
							mod.obtainSSfromMat(fn, mod);
						}
						mod.checkModel(mod);
						models.push_back(mod);
					}
					stateSpaceModelsPerMode = models;
					n = models[0].x_dim;
				}
				else
				{
					std::cout << "File " << fn << " not found." << std::endl;
					exit(0);
				}
			}
		}
		shs_t(const char *fn, SimulationData &data, int NumModes)
		{
			// Obtain number of discrete modes and obtain discreteKernel
			discreteModes = NumModes;
			sigmoid = true;
			discreteKernel = arma::eye<Matrix>(NumModes, NumModes);
			continuousStateData = data;
			Pk = data.initialModeDistribution[0];
			initialDiscreteMode = data.initialMode[0];		// Initial mode
			// create array containing the ssmodels
			std::vector<ssmodels_t> models;
			Matrix dummy = -1 * arma::ones<Matrix>(1, 1);
			for (int i = 1; i <= discreteModes; i++)
			{
				int x_dimension = data.x_k.n_rows;
				int u_dimension = data.u_k.n_rows;
				int d_dimension = data.d_k.n_rows;
				double lhs = arma::accu(data.u_k == dummy);
				double rhs = arma::accu(data.d_k == dummy);
				if ((lhs > 0) && (rhs > 0))
				{
					u_dimension = 0;
					d_dimension = 0;
				}
				else if ((lhs > 0) && (rhs == 0))
				{
					u_dimension = 0;
				}
				ssmodels_t mod(x_dimension, u_dimension, d_dimension);
				this->n = x_dimensiocontinuousSpaceDimension// TODO: update to handle different dimensions for
				// each state (Store)
				// Reading models
					mod.obtainSSfromMat(fn, mod, i);
				if ((unsigned)x_dimension != (unsigned)mod.x_dim ||
					(unsigned)x_dimension != (unsigned)mod.A.n_rows)
				{
					mod.obtainSSfromMat(fn, mod);
				}
				mod.checkModel(mod);
				models.push_back(mod);
			}
			stateSpaceModelsPerMode = models;

			try
			{
				// Check dimensions of data
				int err = this->checkData();
				switch (err)
				{
					case 3:
					{
						throw "Incorrect size of disturbance vector";
					}
					default: { throw "Correct input vectors size"; }
				}
			}
			catch (const char *msg)
			{
				std::cerr << msg << std::endl;
				exit(0);
			}
		}
		virtual ~shs_t() {}
		std::vector<int> t_q(shs_t &old, int q_old)
		{
			int steps = 1;
			std::vector<int> modes(steps + 1);

			modes[0] = q_old;
			std::random_device rd;
			std::mt19937 gen(rd());
			int count;
			double sum, actionVariables;
			for (int i = 0; i < steps; i++)
			{
				count = 0;
				sum = 0;
				actionVariables = std::generate_canonical<double, 10>(gen);
				while (sum < actionVariables)
				{
					sum += old.Pk(modes[i], count);
					if (sum > actionVariables)
					{
						modes[i + 1] = count;
					}
					count++;
				}
			}

			return modes;
		}
		// Dynamic model of continuous state
		void step(shs_t &old, int n, int steps)
		{
			Matrix x_old(old.continuousStateData.continuousVariables[n]);
			Matrix x_new = Matrix(old.stateSpaceModelsPerMode[0].x_dim, steps);
			Matrix u_old = old.continuousStateData.u_k;
			Matrix d_old = Matrix(old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode(n, 0)].d_dim, steps);
			d_old = old.continuousStateData.d_k;
			if (old.discreteModes == 1)
			{
				if (old.stateSpaceModelsPerMode[0].N.is_empty())
				{
					if(!u_old.is_empty())
					{
						if(u_old.n_cols == old.stateSpaceModelsPerMode[0].B.n_rows)
						{
							u_old = u_old.t();
						}
						u_old = arma::repmat(old.continuousStateData.u_k.col(n),1,steps);
					}
					if(!d_old.is_empty())
					{
						if(d_old.n_cols == old.stateSpaceModelsPerMode[0].F.n_rows)
						{
							d_old = d_old.t();
						}
						d_old = arma::repmat(d_old.col(n),1,steps);
					}
					x_new = old.stateSpaceModelsPerMode[0].updateLTIst(old.stateSpaceModelsPerMode[0], x_old, u_old, d_old);
				}
				else
				{
					x_new =
						old.stateSpaceModelsPerMode[0].updateBi(old.stateSpaceModelsPerMode[0].A, old.stateSpaceModelsPerMode[0].B,
						old.stateSpaceModelsPerMode[0].N, old.stateSpaceModelsPerMode[0].discreteModes, x_old, u_old);
				}
			}
			else
			{
				double lhs =
					arma::accu(old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode(n, 0)].sigma ==
					ZeroMatrix(
					(double)old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode(n, 0)].sigma.n_rows,
					(double)old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode(n, 0)].sigma.n_cols));
				if (lhs > 0)
				{
					for (int j = 0; j < steps; j++)
					{
						if (old.sigmoid)
						{
							// Generalise: assumes for now that only have pairs of discrete
							// modes
							old.Pk(0, 0) = sigmoidCompute(x_old(0, 0), 100, 19.5);
							old.Pk(0, 1) = 1 - old.Pk(0, 0);
							old.Pk(1, 1) = 1 - old.Pk(0, 0);
							old.Pk(1, 0) = old.Pk(1, 1);
							std::cout << "Trans: " << old.Pk << std::endl;
						}
						else
						{
							if (n > 0 && !old.discreteKernel.is_empty())
							{
								old.Pk = old.Pk * old.discreteKernel;
							}
						}
						if (!old.discreteKernel.is_empty())
						{
							int initialDiscreteMode = (int)old.initialDiscreteMode(n, j);
							std::vector<int> q = t_q(old, initialDiscreteMode);
							int q_new = q[1];

							if (!old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode(n, j)].N.is_empty())
							{
								x_new = old.stateSpaceModelsPerMode[q_new].updateBi(old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode(n, j)].A,
									old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode(n, j)].B,
									old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode(n, j)].N,
									old.stateSpaceModelsPerMode[old.initialDiscreteMode(n, j)].discreteModes,
									x_old, u_old);
							}
							else if (old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode(n, j)].d_dim == 0)
							{
								x_new.col(j) = old.stateSpaceModelsPerMode[q_new].updateLTI(
									old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode(n, j)], x_old.col(j), u_old.col(j));
							}
							else
							{
								x_new.col(j) = old.stateSpaceModelsPerMode[q_new].updateLTIad(
									old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode(n, j)], x_old.col(j), u_old.col(j), d_old.col(j));
							}
							// Append result to discreteModes in old object (columns)
							old.initialDiscreteMode(n + 1, j) = q_new;
						}
						else
						{
							int q_new = old.continuousStateData.actionVariables(n + 1, 0);
							Matrix emptyMat = ZeroMatrix(1, 1);
							emptyMat.reset();
							x_new = old.stateSpaceModelsPerMode[q_new].updateLTIst(
								old.stateSpaceModelsPerMode[(int)u_old(0, 0)], x_old.col(j), emptyMat, emptyMat);
							// Append result to discreteModes in old object (columns)
							old.initialDiscreteMode(n + 1, 0) = q_new;
						}
					}

				}
				else
				{
					// Get next state by updating Pk
					// Assumes continuous variables are affected by white noise
					// described using gaussian distribution with mean x[k] = x_old,
					// variance sigma
					// Compute conditional distribution for each mode
					for (int j = 0; j < steps; j++)
					{
						int initialDiscreteMode = (int)old.initialDiscreteMode(n, j);
						if (old.sigmoid)
						{
							old.Pk(0, 0) = sigmoidCompute(x_old(0, 0), 1, 19.5) *
								sigmoidCompute(x_old(1, 0), 1, 21.25);
							old.Pk(0, 1) = 1 - old.Pk(0, 0);
							old.Pk(1, 1) = old.Pk(0, 0);
							old.Pk(1, 0) = 1 - old.Pk(0, 0);
						}
						else
						{
							computeConditional(old, x_old, j, initialDiscreteMode);
						}
						// Sample from conditional distribution to get q_new
						std::vector<int> q = t_q(old, initialDiscreteMode);
						int q_new = q[1];
						// Append result to discreteModes in old object (columns)
						old.initialDiscreteMode(n + 1, j) = q_new;
						// Assuming no reset kernels
						Matrix x_up = Matrix(old.stateSpaceModelsPerMode[q_new].x_dim, 1);
						if (!old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode(n, j)].N.is_empty())
						{
							x_up = old.stateSpaceModelsPerMode[q_new].updateBi(old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode(n, j)].A,
								old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode(n, j)].B,
								old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode(n, j)].N,
								old.stateSpaceModelsPerMode[old.initialDiscreteMode(n, j)].discreteModes, x_old,
								u_old);
							for (int k = 0; k < old.stateSpaceModelsPerMode[q_new].x_dim; k++)
							{
								x_new(k, j) =
									getSampleNormal(x_up(k, 0), old.stateSpaceModelsPerMode[q_new].sigma(k, 0));
							}

						}
						else if (old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode(n, j)].d_dim == 0)
						{

							x_up = old.stateSpaceModelsPerMode[q_new].updateLTI(old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode(n, j)],
								x_old.col(j), u_old.col(j));
							for (int k = 0; k < old.stateSpaceModelsPerMode[q_new].x_dim; k++)
							{
								x_new(k, j) =
									getSampleNormal(x_up(k, 0), old.stateSpaceModelsPerMode[q_new].sigma(k, 0));
							}
						}
						else
						{
							x_up = old.stateSpaceModelsPerMode[q_new].updateLTIad(old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode(n, j)],
								x_old.col(j), u_old.col(j), d_old.col(j));
							for (int k = 0; k < old.stateSpaceModelsPerMode[q_new].x_dim; k++)
							{
								x_new(k, j) =
									getSampleNormal(x_up(k, 0), old.stateSpaceModelsPerMode[q_new].sigma(k, 0));
							}
						}
					}
				}
			}
			// Append result to continuousVariables in old object (3rd dimension)
			old.continuousStateData.continuousVariables.push_back(x_new);
		}
		void computeConditional(shs_t &old, Matrix x_old, int index, int q_old)
		{
			double min = -1, max = 1;
			double val;										//,err, xmin[1] ={}, xmax[1]={};
			std::vector<double> v = {};
			for (int i = 0; i < 2 * old.stateSpaceModelsPerMode[q_old].x_dim; i++)
			{
				if (i < old.stateSpaceModelsPerMode[q_old].x_dim)
				{
					v.push_back(x_old(i, index));			//
				}
				else
				{
					// TODO: Case when sigma is correlated
					v.push_back(old.stateSpaceModelsPerMode[q_old].sigma(i - old.stateSpaceModelsPerMode[q_old].x_dim, 0));
				}
			}
			val = 1;
			// Update Pk
			if (n == 0)
			{
				old.Pk = old.discreteKernel;
			}
			else
			{
				old.Pk = old.Pk * val * old.discreteKernel;
			}
			//    std::cout << "New kernel: " << old.Pk << std::endl;
		}
		void createPySimPlots(Matrix y, Matrix modes, int x_dim)
		{
			int T = y.n_rows / x_dim;
			std::cout << "Creating python simulation plot file" << std::endl;
			std::ofstream myfile;
			// If folder does not already exist then
			// create it else store in already existant
			// folder
			if(checkFolderExists("../results") == -1)
			{
				if(mkdir("../results", 0777) == -1)
				{
					std::cerr << "Error cannot create results directory: " <<std::strerror(errno) <<std::endl;
					exit(0);
				}
				else
				{
					std::cout <<" Results direcrtory created at ../results" <<std::endl;
				}
			}

			myfile.open("../results/simPlots.py");

			myfile << "from mpl_toolkits.mplot3d import Axes3D" << std::endl;
			myfile << "import numpy as np" << std::endl;
			myfile << "import matplotlib.pyplot as plt" << std::endl;
			myfile << std::endl;

			// Get separate cont var evolution
			3DMatrix y_all;
			Matrix y1 = ZeroMatrix(T, y.n_cols);

			if (x_dim == 1)
			{
				y1 = y;
			}
			else
			{
				for(size_t j = 0; j < x_dim; j++)
				{
					int count =0;
					for (unsigned i = j; i < y.n_rows; i = i + x_dim)
					{
						y1.row(count) = y.row(i);
						count++;
					}
					y_all.push_back(y1);
					y1 = ZeroMatrix(T, y.n_cols);
				}
			}

			for (unsigned i = 0; i < x_dim; i++)
			{
				myfile << "x" << i << "= [";				// Creating variable for 1 trace
				for (unsigned j = i; j < y.n_rows; j = j + x_dim)
				{
					if (j == (y.n_rows - x_dim + i))
					{
						myfile << y(j, 0) << "]" << std::endl;
					}
					else
					{
						myfile << y(j, 0) << ",";
					}
				}

				myfile << "plt.subplot("<<x_dim<< ",3," << i + 1 << ")" << std::endl;
				myfile << "plt.plot(x" << i << ")" << std::endl;
				myfile << "plt.title('Sample trace of continuous variable $x_" << i + 1
					<< "$',fontsize=18)" << std::endl;
				myfile << "plt.xlabel('Time steps') " << std::endl;
				myfile << "plt.ylabel('Continuous variable $x_" << i + 1 << "$')"
					<< std::endl;
				myfile << std::endl;
			}

			// Now we need to print the discrete modes
			myfile << "modes = [";							// Creating variable for mean values
			for (unsigned j = 0; j < T; j++)
			{
				if (j == (T - 1))
				{
					myfile << modes(j, 0) << "]" << std::endl;
				}
				else
				{
					myfile << modes(j, 0) << ",";
				}
			}
			myfile << "plt.subplot("<<x_dim<< ",3," << x_dim + 1 << ")" << std::endl;
			myfile << "plt.plot(modes,marker='o',drawstyle='steps')" << std::endl;
			myfile << "plt.yticks(np.arange(min(modes),max(modes)+1,1))" << std::endl;
			myfile << "plt.title('Sample trace of discrete modes',fontsize=18)"
				<< std::endl;
			myfile << "plt.xlabel('Time steps') " << std::endl;
			myfile << "plt.ylabel('Discrete modes')" << std::endl;
			myfile << std::endl;

			// Plotting of histograms
			int n_add = x_dim + 1;							// Number of additional subplots

			// Defining time horizon
			for (unsigned i = 0; i < n_add; i++)
			{
				myfile << "ax = plt.subplot("<<x_dim<< ",3," << i + x_dim + 2
					<< ", projection = '3d')" << std::endl;
				// Compute number of bins
				int binNo = 1;
				if (i < n_add -1)
				{
					myfile << "data_2d= [ [";				// Creating data set for historgram from y1
					for (unsigned p = 0; p < y_all[i].n_cols; p++)
					{
						for (unsigned j = 0; j < y_all[i].n_rows; j++)
						{
							if (j == (y_all[i].n_rows - 1))
							{
								myfile << y_all[i](j, p) << "]," << std::endl;
							}
							else
							{
								myfile << y_all[i](j, p) << ",";
							}
						}
						if (p < (y_all[i].n_cols - 1))
						{
							myfile << "[";
						}
					}
					int y1min = (int)arma::min(arma::min(y_all[i]));
					int y1max = (int)arma::max(arma::max(y_all[i]));
					double delta1 = (int)(y1max - y1min) / 5;
					if (delta1 == 0)
					{
						delta1 = 1;
					}

					myfile << "]" << std::endl;
					myfile << "data_array = np.array(data_2d)" << std::endl;
					myfile << std::endl;
					myfile << "# This is  the colormap I'd like to use." << std::endl;
					myfile << "cm = plt.cm.get_cmap('viridis')" << std::endl;
					myfile << std::endl;
					myfile << "#Create 3D histogram" << std::endl;
					myfile << "for z in range(" << T << "):" << std::endl;
					myfile << "     j = np.transpose(data_2d)" << std::endl;
					myfile << "     y_data = j[z][1:]" << std::endl;
					myfile << "     hist, bins = np.histogram(y_data,bins=75)"
						<< std::endl;
					myfile << "     center = (bins[:-1] + bins[1:])/2" << std::endl;
					myfile << "     x_span = bins.max() - bins.min()" << std::endl;
					myfile << "     C = [cm((x-bins.min())/x_span) for x in bins]"
						<< std::endl;
					myfile << "     ax.bar(center, hist, zs =z, zdir= 'y',color =C, fill "
						"= 'true')"
						<< std::endl;
					myfile << "     ax.set_ylabel('Time steps') " << std::endl;
					myfile << "     ax.set_xlabel('Continuous variable $x_" << i + 1 << "$')"
						<< std::endl;
					myfile << "     ax.set_zlabel('Count')" << std::endl;
					myfile << "plt.ylim(" << 1 << "," << T << ")" << std::endl;
					myfile << "plt.xticks(np.arange(" << y1min << " ," << y1max
						<< ",step=" << delta1 << "))" << std::endl;
					myfile << std::endl;

				}
				else if (i == (n_add-1))
				{
					myfile << "data_2d2= [ [";				// Creating data set for historgram from y1
					for (unsigned p = 0; p < modes.n_cols; p++)
					{
						for (unsigned j = 0; j < modes.n_rows; j++)
						{
							if (j == (modes.n_rows - 1))
							{
								myfile << modes(j, p) << "]," << std::endl;
							}
							else
							{
								myfile << modes(j, p) << ",";
							}
						}
						if (p < (y1.n_cols - 1))
						{
							myfile << "[";
						}
					}

					myfile << "]" << std::endl;
					int mmax = (int)arma::max(arma::max(modes));
					myfile << "data_array2 = np.array(data_2d2)" << std::endl;
					myfile << std::endl;
					myfile << "#Create 3D histogram" << std::endl;
					myfile << "for z2 in range(" << T << "):" << std::endl;
					myfile << "     j2 = np.transpose(data_2d2)" << std::endl;
					myfile << "     y_data2 = j2[z2][1:]" << std::endl;
					myfile << "     hist2, bins2 = np.histogram(y_data2,bins=" << mmax + 1
						<< ")" << std::endl;
					myfile << "     center2 = (bins2[:-1] + bins2[1:])/2" << std::endl;
					myfile << "     x_span2 = bins2.max() - bins2.min()" << std::endl;
					myfile << "     C = [cm((x2-bins2.min())/x_span2) for x2 in bins2]"
						<< std::endl;
					myfile << "     ax.bar(center2, hist2, zs =z2, zdir= 'y', fill = "
						"'true', color=C)"
						<< std::endl;
					myfile << "     ax.set_ylabel('Time steps') " << std::endl;
					myfile << "     ax.set_xlabel('Discrete variable $q$')" << std::endl;
					myfile << "     ax.set_zlabel('Count')" << std::endl;
					myfile << "plt.ylim(" << 1 << "," << T << ")" << std::endl;
					myfile << "plt.xticks(np.arange(0," << mmax + 1 << ",step=1))"
						<< std::endl;

				}
			}
			myfile << "plt.savefig(\"SimulationRun.svg\",dpi=150)\n";
			myfile << "plt.show()" << std::endl;
			myfile.close();

		}
		void run(shs_t &old, int N, int steps)
		{
			clock_t begin, end;
			begin = clock();
			double time = 0;

			// For deterministic case perform 1 run
			// For stochastic version perform Monte Carlo + compute mean
			int i = 0, x_dim = old.stateSpaceModelsPerMode[0].A.n_rows;

			Matrix y = ZeroMatrix(N * x_dim, steps);
			Matrix modes = ZeroMatrix(N, steps);
			old.Pk.set_size(size(old.discreteKernel));
			old.Pk = old.discreteKernel;
			int count = 0;
			while (i < N)
			{
				old.step(old, i, steps);
				if (old.stateSpaceModelsPerMode[0].u_dim > 0)
				{
					if ((unsigned)old.continuousStateData.actionVariables.n_rows == (unsigned)N)
					{
						old.continuousStateData.u_k = old.continuousStateData.actionVariables.row(i);
					}
				}
				if (old.stateSpaceModelsPerMode[0].d_dim > 0)
				{
					if ((unsigned)old.continuousStateData.disturbanceVariables.n_rows == (unsigned)N)
					{
						old.continuousStateData.d_k = old.continuousStateData.disturbanceVariables.row(i);
					}
				}
				Matrix tempX = old.continuousStateData.continuousVariables[i];
				if(x_dim == 1)
				{
					y.row(count) = tempX;
				}
				else
				{
					y.rows(count, count + x_dim-1) = tempX;
				}
				modes.row(i) = old.initialDiscreteMode.row(i);
				count += x_dim;
				i++;
			}
			end = clock();
			time = (double)(end - begin) / CLOCKS_PER_SEC;
			std::ostringstream oss;
			auto t = std::time(nullptr);
			auto tm = *std::localtime(&t);

			oss << std::put_time(&tm, "%d-%m-%Y-%H-%M-%S");
			auto str = oss.str();

			std::cout << std::endl;
			std::cout << "--------------------------------------" <<std::endl;
			std::cout << " Simulation time                      " << std::endl;
			std::cout << "--------------------------------------" <<std::endl;
			std::cout << " " << time << std::endl;
			std::cout << "--------------------------------------" <<std::endl;
			std::cout << std::endl;

			// Option to export to file results
			std::ofstream myfile;
			std::string exportOpt;
			std::string str0("y");
			std::string str1("yes");
			std::cout << "Would you like to store simulation results [y- yes, n - no] " << std::endl;
			std::cin >> exportOpt;
			if ((exportOpt.compare(str0) == 0) || (exportOpt.compare(str1) == 0))
			{
				// check if results folder exists:
				if(checkFolderExists("../results") == -1)
				{
					if(mkdir("../results", 0777) == -1)
					{
						std::cerr << "Error cannot create results directory: " <<std::strerror(errno) <<std::endl;
						exit(0);
					}
				}
				// Store results in file
				std::string f_name = "../results/Simulationtime_" + str + ".txt";
				myfile.open(f_name);
				myfile << time <<std::endl;
				myfile.close();
				std::string y_name = "../results/y_" + str + ".txt";
				y.save(y_name, arma::raw_ascii);
				Matrix q = old.initialDiscreteMode;
				std::string q_name = "../results/modes_" + str + ".txt";
				q.save(q_name, arma::raw_ascii);
				old.createPySimPlots(y, modes, x_dim);
			}

		}
		bool obtainTqfromMat(const char *fn, shs_t &init)
		{
			// Reading model file input in .Matrix format
			// and storing into ssmodel class
			mat_t *matf;
			matvar_t *matvar, *contents;
			// Read .mat file
			bool error = 0;

			matf = Mat_Open(fn, MAT_ACC_RDONLY);
			if (matf)										// if successful in reading file
			{
				// read each variable within file and populate
				// state space model based on variable name
				contents = Mat_VarRead(matf, "discreteKernel");
				if (contents == NULL)
				{
					init.discreteModes = 1;
				}
				else
				{
					init.populateTq(*contents);
					contents = NULL;
				}
				Mat_Close(matf);

			} else											// unsuccessfull in opening file
			{
				return 1;
			}
		}

	private:
		int checkData()
		{
			int error = 0;
			if ((unsigned)this->stateSpaceModelsPerMode[0].x_dim != (unsigned)this->continuousStateData.continuousVariables[0].n_rows)
			{
				error = 1;
			}
			if ((unsigned)this->stateSpaceModelsPerMode[0].u_dim != (unsigned)this->continuousStateData.actionVariables.n_cols &&
				(unsigned)this->stateSpaceModelsPerMode[0].u_dim > 0)
			{
				error = 2;
			}
			if ((unsigned)this->stateSpaceModelsPerMode[0].d_dim != (unsigned)this->continuousStateData.disturbanceVariables.n_cols &&
				(unsigned)this->stateSpaceModelsPerMode[0].d_dim > 0)
			{
				error = 3;
			}
			return error;
		}
		void populateTq(matvar_t &content)
		{
			ssmodels_t container;

			// discreteKernel can be of two version
			// Cells containing strings for guards or
			// Numeric containing probabilities
			if (content.data != NULL)
			{
				if (content.data_type == MAT_T_DOUBLE)
				{
					std::string str;
					size_t stride = Mat_SizeOf(content.data_type);
					char *data = (char *)content.data;
					unsigned i, j = 0;
					for (i = 0; i < content.dims[0]; i++)
					{
						for (j = 0; j < content.dims[1]; j++)
						{
							size_t idx = content.dims[0] * j + i;
							void *t = data + idx * stride;
							char substr[100];
							sprintf(substr, "%g",
								*(double *)t);				// Assumes values are of type double
							str.append(substr);
							str.append(" ");
						}
						str.append(";");
					}
					std::vector<std::string> x = splitStr(str, ';');
					int numEl = x.size();
					this->discreteModes = numEl;

					// Check stochasticity of kernel
					Matrix tq = strtodMatrix(x);
					tq = checkStochasticity(tq);
					this->discreteKernel = tq;
				}
				else
				{
					std::cout << "Incorrect discreteKernel format" << std::endl;
				}
			}
			else
			{
				std::cout << "discreteKernel field not input in Matrix file" << std::endl;
			}
		}
};

/*************************************************************************************************/
// SHS with multiple initial modes
template <> class shs_t<Matrix, std::vector<int>>
{
	public:
		std::vector<int> discreteModes;						// Discrete modes
		3DMatrix
			discreteKernel;									// Discrete kernel- can represent either guards or probabilities
		int continuousSpaceDimension						// Dimension of continuous state space in given mode
			3DMatrix Pk;
		3DMatrix initialDiscreteMode;						// Modes
		std::vector<ssmodels_t> stateSpaceModelsPerMode;	// container of models
		SimulationData continuousStateData;

	private:
		bool sigmoid;

	public:
		shs_t()
		{
			discreteModes =									// Discrete modes
			{
				1
			};
			n = 1;											// Dimension of continuous state space in given mode
			Pk = {Matrix(discreteModes[0], n)};
			initialDiscreteMode =							// Initial mode
			{
				Matrix(discreteModes[0], 1)
			};
			std::vector<ssmodels_t> stateSpaceModelsPerMode;
			SimulationData continuousStateData;
			discreteKernel = {Matrix(n, n)};
			sigmoid = false;
		}
		shs_t(std::vector<int> disc, int num, 3DMatrix tq,
			std::vector<ssmodels_t> &old, SimulationData &data)
		{
			discreteModes = disc;							// Discrete modes
			n = num;										// Dimension of continuous state space in given mode
			Pk = data.initialModeDistribution;
			initialDiscreteMode = data.initialMode;			// Initial mode
			stateSpaceModelsPerMode = old;
			continuousStateData = data;
			discreteKernel = tq;
			sigmoid = false;
			try
			{
				// Check dimensions of data
				int err = this->checkData();
				switch (err)
				{
					case 3:
					{
						throw "Incorrect size of disturbance vector";
					}
					default: { throw "Correct input vectors size"; }
				}
			}
			catch (const char *msg)
			{
				std::cerr << msg << std::endl;
				exit(0);
			}
		}
		shs_t(3DMatrix tq, std::vector<ssmodels_t> &old,
			SimulationData &data)
		{
			discreteModes =									// Discrete modes
			{
				1
			};
			n = old[0].x_dim;								// Dimension of continuous state space in given mode
			Pk = data.initialModeDistribution;
			initialDiscreteMode = data.initialMode;			// Initial mode
			stateSpaceModelsPerMode = old;
			continuousStateData = data;
			discreteKernel = tq;
			sigmoid = false;
			try
			{
				// Check dimensions of data
				int err = this->checkData();
				switch (err)
				{
					case 3:
					{
						throw "Incorrect size of disturbance vector";
					}
					default: { throw "Correct input vectors size"; }
				}
			}
			catch (const char *msg)
			{
				std::cerr << msg << std::endl;
				exit(0);
			}
		}
		shs_t(const char *fn, SimulationData &data, std::vector<int> NumModes)
		{
			// Obtain number of discrete modes and obtain discreteKernel
			std::vector<int> num = {2, 2};

			discreteModes = num;

			sigmoid = true;
			for (unsigned int m = 0; m < num.size(); m++)
			{
				discreteKernel.push_back(arma::eye<Matrix>(num[m], num[m]));
				Pk.push_back(data.initialModeDistribution[m]);
															// Initial mode
				initialDiscreteMode.push_back(data.initialMode[m]);
			}
			continuousStateData = data;

			// create array containing the ssmodels
			std::vector<ssmodels_t> models;
			Matrix dummy = -1 * arma::ones<Matrix>(1, 1);
			for (int i = 1; i <= discreteModes[0]; i++)
			{
				int x_dimension = data.x_k.n_rows;
				int u_dimension = data.u_k.n_rows;
				int d_dimension = data.d_k.n_rows;
				double lhs = arma::accu(data.u_k == dummy);
				double rhs = arma::accu(data.d_k == dummy);
				if ((lhs > 0) && (rhs > 0))
				{
					u_dimension = 0;
					d_dimension = 0;
				}
				else if ((lhs > 0) && (rhs == 0))
				{
					u_dimension = 0;
				}
				ssmodels_t mod(x_dimension, u_dimension, d_dimension);
				this->n = x_dimensiocontinuousSpaceDimension// TODO: update to handle different dimensions for
				// each state (Store)
				// Reading models
					mod.obtainSSfromMat(fn, mod, i);
				if ((unsigned)x_dimension != (unsigned)mod.x_dim ||
					(unsigned)x_dimension != (unsigned)mod.A.n_rows)
				{
					mod.obtainSSfromMat(fn, mod);
				}
				mod.checkModel(mod);
				models.push_back(mod);
			}
			stateSpaceModelsPerMode = models;
			try
			{
				// Check dimensions of data
				int err = this->checkData();
				switch (err)
				{
					case 3:
					{
						throw "Incorrect size of disturbance vector";
					}
					default: { std::cout << "Correct input vectors size" << std::endl; }
				}
			}
			catch (const char *msg)
			{
				std::cerr << msg << std::endl;
				exit(0);
			}
		}

		shs_t(const char *fn, std::vector<int> NumModes)
		{
			// Obtain number of discrete modes and obtain discreteKernel
			std::vector<int> num = {2, 2};

			discreteModes = num;

			sigmoid = true;
			for (unsigned int m = 0; m < num.size(); m++)
			{
				discreteKernel.push_back(arma::eye<Matrix>(num[m], num[m]));
			}
			// create array containing the ssmodels
			std::vector<ssmodels_t> models;
			Matrix dummy = -1 * arma::ones<Matrix>(1, 1);
			for (int i = 1; i <= discreteModes[0]; i++)
			{
				int x_dimension = this->stateSpaceModelsPerMode[i].A.n_rows;
				int u_dimension = this->stateSpaceModelsPerMode[i].B.n_rows;
				int d_dimension = this->stateSpaceModelsPerMode[i].F.n_rows;
				ssmodels_t mod(x_dimension, u_dimension, d_dimension);
				this->n = x_dimensiocontinuousSpaceDimension// TODO: update to handle different dimensions for
				// each state (Store)
				// Reading models
					std::cout
					<< "Initialising model of continuous variables in discrete mode " << i
					<< std::endl;
				mod.obtainSSfromMat(fn, mod, i);
				if ((unsigned)x_dimension != (unsigned)mod.x_dim ||
					(unsigned)x_dimension != (unsigned)mod.A.n_rows)
				{
					mod.obtainSSfromMat(fn, mod);
				}
				mod.checkModel(mod);
				models.push_back(mod);
			}
			stateSpaceModelsPerMode = models;
		}

		virtual ~shs_t() {}
		std::vector<int> t_q(shs_t &old, int q_old, int index)
		{
			int steps = 1;
			std::vector<int> modes(steps + 1);

			modes[0] = q_old;
			// if(old.stateSpaceModelsPerMode[0].sigma.isZero(0) )
			//{
			std::random_device rd;
			std::mt19937 gen(rd());
			int count;
			double sum, actionVariables;
			for (int i = 0; i < steps; i++)
			{
				count = 0;
				sum = 0;
				actionVariables = std::generate_canonical<double, 10>(gen);
				while (sum < actionVariables)
				{
					sum += old.Pk[index](modes[i], count);
					if (sum > actionVariables)
					{
						modes[i + 1] = count;
					}
					count++;
				}
			}

			return modes;
		}
		// Dynamic model of continuous state
		void step(shs_t &old, int n, int steps)
		{
			Matrix x_old(old.continuousStateData.continuousVariables[n]);
			// Matrix x_old = Matrix(old.stateSpaceModelsPerMode[(int) old.initialDiscreteMode[0](n,0)].x_dim,1);
			// x_old = old.continuousStateData.continuousVariables[n];
			Matrix x_new = Matrix(old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode[0](n, 0)].x_dim, steps);

			// std::vector<double> vec(x_old.size());
			// Eigen::Map<Matrix>(vec.data(), x_old.n_rows, x_old.n_cols) =x_old; //

			Matrix u_old = Matrix(old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode[0](n, 0)].u_dim, 1);
			u_old = old.continuousStateData.u_k;
			Matrix d_old = Matrix(old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode[0](n, 0)].d_dim, 1);
			d_old = old.continuousStateData.d_k;

			if (old.discreteModes[0] == 1)
			{
				x_new = old.stateSpaceModelsPerMode[0].updateLTIst(old.stateSpaceModelsPerMode[0], x_old, u_old, d_old);
			}
			else
			{

				double equal =
					arma::accu(old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode[0](n, 0)].sigma ==
					ZeroMatrix(
					old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode[0](n, 0)].sigma.n_rows,
					old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode[0](n, 0)].sigma.n_cols));
				if (equal > 0)
				{
					for (int j = 0; j < steps; j++)
					{
						int initialDiscreteMode = (int)old.initialDiscreteMode[0](n, j);
						if (old.sigmoid)
						{
							// Generalise: assumes for now that only have pairs of discrete
							// modes
							old.Pk[0](0, 0) = sigmoidCompute(x_old(0, 0), 100, 19.5);
							old.Pk[0](0, 1) = 1 - old.Pk[0](0, 0);
							old.Pk[0](1, 1) = 1 - old.Pk[0](0, 0);
							old.Pk[0](1, 0) = old.Pk[0](1, 1);
							std::cout << "Trans: " << old.Pk[0] << std::endl;
						}
						else
						{
							if (n > 0)
							{
								old.Pk[0] = old.Pk[0] * old.discreteKernel[0];
							}
						}
						std::vector<int> q = t_q(old, initialDiscreteMode, 0);
						int q_new = q[1];
						if (old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode[0](n, j)].d_dim == 0)
						{
							x_new.col(j) = old.stateSpaceModelsPerMode[q_new].updateLTI(
								old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode[0](n, j)], x_old.col(j), u_old);
						}
						else
						{
							x_new.col(j) = old.stateSpaceModelsPerMode[q_new].updateLTIad(
								old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode[0](n, j)], x_old.col(j), u_old, d_old);
						}
						// Append result to discreteModes in old object (columns)
						old.initialDiscreteMode[0](n + 1, j) = q_new;
					}
				}
				else
				{
					// Get next state by updating Pk
					// Assumes continuous variables are affected by white noise
					// described using gaussian distribution with mean x[k] = x_old,
					// variance sigma
					// Compute conditional distribution for each mode
					for (int j = 0; j < steps; j++)
					{
						std::vector<int> initialDiscreteMode;
						initialDiscreteMode.push_back(old.initialDiscreteMode[0](n, j));
						initialDiscreteMode.push_back(old.initialDiscreteMode[1](n, j));
						std::cout << "initialDiscreteMode: " << initialDiscreteMode[0] << ", x_old: " << x_old << " ,j: " << j
							<< std::endl;
						std::cout << "initialDiscreteMode: " << initialDiscreteMode[1] << ", x_old: " << x_old << " ,j: " << j
							<< std::endl;

						if (old.sigmoid)
						{

							// TODO:Generalise
							old.Pk[0](0, 0) = sigmoidCompute(
								x_old(0, 0), 10, 19.5);		//*sigmoidCompute(x_old(1,0),1,21.25);
							old.Pk[0](0, 1) = 1 - old.Pk[0](0, 0);
							old.Pk[0](1, 1) = 1 - old.Pk[0](0, 0);
							old.Pk[0](1, 0) = old.Pk[0](0, 0);
							old.Pk[1](0, 0) = sigmoidCompute(
								x_old(1, 0), 10, 21.25);	//*sigmoidCompute(x_old(0,0),1,19.5);
							old.Pk[1](0, 1) = 1 - old.Pk[1](0, 0);
							old.Pk[1](1, 1) = 1 - old.Pk[1](0, 0);
							old.Pk[1](1, 0) = old.Pk[1](0, 0);
							std::cout << "Trans: " << old.Pk[0] << std::endl;
							std::cout << "Trans: " << old.Pk[1] << std::endl;
						}
						else
						{
							computeConditional(old, x_old, j, initialDiscreteMode[0]);
						}

						// Sample from conditional distribution to get q_new
						std::vector<int> initialDiscreteMode = t_q(old, initialDiscreteMode[0], 0);
						std::vector<int> q1 = t_q(old, initialDiscreteMode[1], 1);
						std::vector<int> q_new = {initialDiscreteMode[1], q1[1]};
						// Append result to discreteModes in old object (columns)
						old.initialDiscreteMode[0](n + 1, j) = q_new[0];
						old.initialDiscreteMode[1](n + 1, j) = q_new[1];

						// Assuming no reset kernels
						Matrix x_up = Matrix(1, 1);

						if (old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode[0](n, j)].d_dim == 0)
						{
							for (int k = 0; k < old.stateSpaceModelsPerMode[q_new[0]].x_dim; k++)
							{

								x_up = old.stateSpaceModelsPerMode[q_new[k]].updateLTI(
									old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode[k](n, j)].A.row(k),
									old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode[k](n, j)].discreteModes.row(k), x_old.col(j));
								x_new(k, j) =
									getSampleNormal(x_up(0, 0), old.stateSpaceModelsPerMode[q_new[k]].sigma(k, 0));
							}
						}
						else
						{
							for (int k = 0; k < old.stateSpaceModelsPerMode[q_new[0]].x_dim; k++)
							{
								x_up = old.stateSpaceModelsPerMode[q_new[k]].updateLTIad(
									old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode[k](n, j)].A.row(k),
									old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode[k](n, j)].B.row(k),
									old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode[k](n, j)].F.row(k),
									old.stateSpaceModelsPerMode[(int)old.initialDiscreteMode[k](n, j)].discreteModes.row(k), x_old.col(j),
									u_old.row(k), d_old.row(k));
								x_new(k, j) =
									getSampleNormal(x_up(0, 0), old.stateSpaceModelsPerMode[q_new[k]].sigma(k, 0));
							}
						}
					}
				}
				// Append result to continuousVariables in old object (3rd dimension)
				old.continuousStateData.continuousVariables.push_back(x_new);
			}
		}
		void computeConditional(shs_t &old, Matrix x_old, int index, int q_old)
		{
			double min = -1, max = 1;
			double val, err, xmin[1] =
			{
			}
			, xmax[1] =
			{
			};
			std::vector<double> v = {};
			for (int i = 0; i < 2 * old.stateSpaceModelsPerMode[q_old].x_dim; i++)
			{
				if (i < old.stateSpaceModelsPerMode[q_old].x_dim)
				{
					xmin[i] = micontinuousSpaceDimension
						xmax[i] = max;
					v.push_back(x_old(i, index));			//
					std::cout << "v " << v[i] << std::endl;
				}
				else
				{
					// TODO: Case when sigma is correlated
					std::cout << "i-old.stateSpaceModelsPerMode[q_old].x_dim: " << old.stateSpaceModelsPerMode[q_old].x_dim
						<< std::endl;
					std::cout << "q_old " << q_old << std::endl;

					v.push_back(old.stateSpaceModelsPerMode[q_old].sigma(i - old.stateSpaceModelsPerMode[q_old].x_dim, 0));
					std::cout << "v " << v[i] << std::endl;
				}
			}
			val = 1;										// hcubature(1, f_Gauss, &v,old.stateSpaceModelsPerMode[n].x_dim, xmin,xmax,
			// 0,0,1e-12, ERROR_INDIVIDUAL, &val,&err);
			// std::cout << "Integral: " << val << std::endl;
			// Update Pk
			if (n == 0)
			{
				old.Pk[0] = old.discreteKernel[0];
			}
			else
			{
				old.Pk[0] = old.Pk[0] * val * old.discreteKernel[0];
			}
			// std::cout << "New kernel: " << old.Pk[0] << std::endl;
		}
		void run(shs_t &old, int N, int steps)
		{
			// Start simulation timers
			clock_t begin, end;
			begin = clock();
			double time = 0;
			end = clock();

			// For deterministic case perform 1 run
			// For stochastic version perform Monte Carlo + compute mean
			int i = 0;
			int x_dim =old.stateSpaceModelsPerMode[0].x_dim ;
			Matrix y = ZeroMatrix(N * x_dim, steps);
			Matrix modes = ZeroMatrix(N, steps);
			old.Pk = old.discreteKernel;
			int count = 0;
			while (i < N)
			{
				old.step(old, i, steps);
				if (old.stateSpaceModelsPerMode[0].u_dim > 0)
				{
					if ((unsigned)old.continuousStateData.actionVariables.n_rows == (unsigned)N)
					{
						old.continuousStateData.u_k = old.continuousStateData.actionVariables.row(i);
					}
				}
				if (old.stateSpaceModelsPerMode[0].d_dim > 0)
				{
					if ((unsigned)old.continuousStateData.disturbanceVariables.n_rows == (unsigned)N)
					{
						old.continuousStateData.d_k = old.continuousStateData.disturbanceVariables.row(i);
					}
				}
				Matrix tempX = old.continuousStateData.continuousVariables[i];
				if(old.stateSpaceModelsPerMode[0].x_dim == 1)
				{
					y.row(count) = tempX;
				}
				else
				{
					y.rows(count, count + 1) = tempX;
				}
				modes.row(i) = old.initialDiscreteMode[0].row(i);
				count += old.stateSpaceModelsPerMode[0].x_dim ;
				i++;
			}

			time = (double)(end - begin) / CLOCKS_PER_SEC;
			std::ostringstream oss;
			auto t = std::time(nullptr);
			auto tm = *std::localtime(&t);

			oss << std::put_time(&tm, "%d-%m-%Y-%H-%M-%S");
			auto str = oss.str();

			std::cout << std::endl;
			std::cout << "--------------------------------------" <<std::endl;
			std::cout << " Simulation time                      " << std::endl;
			std::cout << "--------------------------------------" <<std::endl;
			std::cout << " " << time << std::endl;
			std::cout << "--------------------------------------" <<std::endl;
			std::cout << std::endl;

			// Option to export to file results
			std::ofstream myfile;
			std::string exportOpt;
			std::string str0("y");
			std::string str1("yes");
			std::cout << "Would you like to store simulation results [y- yes, n - no] " << std::endl;
			std::cin >> exportOpt;
			if ((exportOpt.compare(str0) == 0) || (exportOpt.compare(str1) == 0))
			{
				// check if results folder exists:
				if(checkFolderExists("../results") == -1)
				{
					if(mkdir("../results", 0777) == -1)
					{
						std::cerr << "Error cannot create results directory: " <<std::strerror(errno) <<std::endl;
						exit(0);
					}
				}
				// Store results in file
				std::string f_name = "../results/Simulationtime_" + str + ".txt";
				myfile.open(f_name);
				myfile << time <<std::endl;
				myfile.close();
				std::string y_name = "../results/y_" + str + ".txt";
				y.save(y_name, arma::raw_ascii);
				Matrix q = old.initialDiscreteMode[0];
				std::string q_name = "../results/modes_" + str + ".txt";
				q.save(q_name, arma::raw_ascii);
			}
		}
		void obtainTqfromMat(const char *fn, shs_t &init)
		{
			// Reading model file input in .Matrix format
			// and storing into ssmodel class
			mat_t *matf;
			matvar_t *matvar, *contents;
			// Read .Matrix file
			try
			{
				matf = Mat_Open(fn, MAT_ACC_RDONLY);
				if (matf)									// if successful in reading file
				{
					// read each variable within file and populate
					// state space model based on variable name
					contents = Mat_VarRead(matf, "discreteKernel");
					if (contents == NULL)
					{
						std::cout << "Variable discreteKernel not found in file" << std::endl;
						std::cout << "Number of modes set to 1" << std::endl;
						init.discreteModes[0] = 1;
					}
					else
					{
						init.populateTq(*contents);
						// Mat_VarFree(matvar);
						matvar = NULL;
						contents = NULL;
					}
				} else										// unsuccessfull in opening file
				{
					throw "Error opening mat file";
				}
				Mat_Close(matf);
			}
			catch (const char *msg)
			{
				std::cout << msg << std::endl;
				exit(0);
			}
		}

	private:
		int checkData()
		{
			int error = 0;
			if ((unsigned)this->stateSpaceModelsPerMode[0].x_dim != (unsigned)this->continuousStateData.continuousVariables[0].n_rows)
			{
				error = 1;
			}
			if ((unsigned)this->stateSpaceModelsPerMode[0].u_dim != (unsigned)this->continuousStateData.U.n_rows &&
				(unsigned)this->stateSpaceModelsPerMode[0].u_dim > 0)
			{
				error = 2;
			}
			if ((unsigned)this->stateSpaceModelsPerMode[0].d_dim != (unsigned)this->continuousStateData.disturbanceVariables.n_rows &&
				(unsigned)this->stateSpaceModelsPerMode[0].d_dim > 0)
			{
				error = 3;
			}
			return error;
		}
		void populateTq(matvar_t &content)
		{
			ssmodels_t container;

			// discreteKernel can be of two version
			// Cells containing strings for guards or
			// Numeric containing probabilities
			if (content.data != NULL)
			{
				if (content.data_type == MAT_T_DOUBLE)
				{
					std::string str;
					size_t stride = Mat_SizeOf(content.data_type);
					char *data = (char *)content.data;
					unsigned i, j = 0;
					for (i = 0; i < content.dims[0]; i++)
					{
						for (j = 0; j < content.dims[1]; j++)
						{
							size_t idx = content.dims[0] * j + i;
							void *t = data + idx * stride;
							char substr[100];
							sprintf(substr, "%g",
								*(double *)t);				// Assumes values are of type double
							str.append(substr);
							str.append(" ");
						}
						str.append(";");
					}
					std::vector<std::string> x = splitStr(str, ';');
					int numEl = x.size();
					std::cout << numEl << std::endl;
					std::cout << "modes: " << numEl << std::endl;
					this->discreteModes[0] = numEl;

					// Check stochasticity of kernel
					Matrix tq = strtodMatrix(x);
					tq = checkStochasticity(tq.t());
					this->discreteKernel[0] = tq;
				}
				else
				{
					std::cout << "Incorrect discreteKernel format" << std::endl;
				}
			}
			else
			{
				std::cout << "discreteKernel field not input in Matrix file" << std::endl;
			}
		}
};
#endif
