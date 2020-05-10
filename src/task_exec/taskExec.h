/*
 * taskExec.h
 *
 *  Created on: 19 Feb 2018
 *      Author: nathalie
 */

#ifndef TASKEXEC_H_
#define TASKEXEC_H_

#include "Bmdp.h"
#include "MDP.h"
#include "InputSpec.h"

using Matrix = Matrix;
using Identity = arma::eye<arma::mat>;
using ZeroMatrix = arma::zeros<arma::mat>;

static void performTask(inputSpec_t<Matrix, int> input)
{
	switch (input.myTask.task)
	{
		case SIMULATIOtimeHorizon:							// Perform simulation depending on model type
		{
			if (input.myTask.timeHorizon <= 0)
			{
				std::cout << "ERROR: Incorrect length of time horizon, value needs to be "
					"positive and greater then 0" << std::endl;
				exit(0);
			}
			if (input.myTask.runs <= 0)
			{
				std::cout << "ERROR: Incorrect number of monte carlo simulations, value needs "
					"to be positive and greater then 0"
					<< std::endl;
				std::cout << "INFO:: Will proceed with default number of simulations i.e. 5000 "
					"simulation runs" << std::endl;
				input.myTask.runs = 5000;
			}

			// Check model
			int num_cont = input.myModel.continuousSpaceDimension;
			for (int i = 0; i < input.myModel.discreteModes; ++i)
			{
				if (input.myModel.stateSpaceModelsPerMode[i].A.n_rows != num_cont)
				{
					std::cout << "ERROR: Different number of continuous variables for each mode, "
						"currently not supported. Work in progress"
						<< std::endl;
					exit(0);
				}
				if (input.myModel.stateSpaceModelsPerMode[i].Sigma.n_rows != num_cont)
				{
					std::cout << "ERROR: Different number of continuous variables for each mode, "
						"currently not supported. Work in progress"
						<< std::endl;
					exit(0);
				}
			}
			if (input.myModel.discreteKernel.n_rows != input.myModel.stateSpaceModelsPerMode.size())
			{
				std::cout << "ERROR: Incorrect transition matrix. timeHorizoneed to have number of rows "
					"== number of columns == number of models representing the "
					"evolution of the continuous variables."
					<< std::endl;
				exit(0);
			}
			input.myModel.run(input.myModel, input.myTask.timeHorizon, input.myTask.runs);
			break;
		}
		case MDP_ABSTRACTION:								// Perform verification task depending on model and tool to be interfaced with
		{
			faust_t myF;
			shs_t<Matrix, int> myModel = input.myModel;
			myF.model = myModel;
			try
			{
				// Check if dynamics are ill-scaled
				// Obtaining problem definition
				clock_t begin, end;
				begin = clock();
				double time = 0;

				PropertyType Problem = input.myTask.propertySpec;
				GridType Gridding = input.myTask.gridType;
				KernelType Distribution = input.myTask.assumptionsKernel;
				bool controlled = input.myTask.iscontrolledled();

				// controlled variables
				double epsilon = input.myTask.eps;
				int timeHorizon = input.myTask.timeHorizon;
				shs_t<Matrix, int> model = input.myModel;

				// Get Safe, Input and Target set
				Matrix SafeSet = input.myTask.safeSet;
				Matrix InputSet = input.myTask.inputSet;
				Matrix TargetSet = input.myTask.targetSet;

				faust_t taskMDP = myF;
				// Check correctness of SafeSet input
				if (SafeSet.n_cols != 2)
				{
					std::cout << "ERROR: There is no correctly defined Safe Set" << std::endl;
					exit(0);
				}
				double temp = arma::accu((SafeSet.col(1) - SafeSet.col(0)) < 0);
				if (temp > 0)
				{
					std::cout << "ERROR: The edges of the Safe Set must have positive length. "
						"Make sure that the first column is the lower bound and "
						"the second column is the upper bound" << std::endl;
					exit(0);
				}

				// Check if Safe Set needs rescaling
				// Check if need to rescale boundary
				Matrix max_Ss = arma::max(SafeSet);
				Matrix min_Ss = arma::min(SafeSet);
				double diff_Ss = max_Ss(1) - min_Ss(0);
				Matrix j = Identity(SafeSet.n_rows, SafeSet.n_rows);
				if(diff_Ss > 50)
				{
					// Identify row with smallest values to rescale to that
					Matrix min_row = arma::min(SafeSet,0);
					Matrix OrSS = SafeSet ;					// The original Safeset definition
					Matrix to_inv = OrSS;
					for(unsigned i =0; i < OrSS.n_rows; i++)
					{
						if(min_row(0,1) != SafeSet(i,1))
						{
							SafeSet.row(i) = min_row;
						}
					}
					Matrix y_inv = arma::diagmat(min_row);
					Matrix j_bar = OrSS/y_inv;
					j = j_bar.replace(arma::datum::inf, 0);
					for(unsigned k = 0; k < model.stateSpaceModelsPerMode.size(); k++)
					{
						taskMDP.model.stateSpaceModelsPerMode[k].Sigma = arma::inv(j)*model.stateSpaceModelsPerMode[k].Sigma;
					}
				}
				if (controlled)
				{
					if (model.stateSpaceModelsPerMode[0].B.n_cols != InputSet.n_rows)
					{
						std::cout << "ERROR: There is no correctly defined Input Set" << std::endl;
						exit(0);
					}
				}

				// Check if dimensions are correct
				if ((unsigned)model.stateSpaceModelsPerMode[0].A.n_cols != SafeSet.n_rows)
				{
					std::cout << "ERROR: The dimension of the Kernel does not match the "
						"dimensions of the Safe Set" << std::endl;
					exit(0);
				}

				// Check correctness of target set if problem of reach and avoid
				if (Problem == VERIFY_REACH_AVOID)
				{
					if (TargetSet.n_cols != 2)
					{
						std::cout << "ERROR: There is no correctly defined Target Set";
						exit(0);
					}
					double temp = arma::accu((TargetSet.col(1) - TargetSet.col(0)) < 0);
					if (temp > 0)
					{
						std::cout << "ERROR: The edges of the Target Set must have positive length. "
							"Make sure that the first column is the lower bound and "
							"the second column is the upper bound";
						exit(0);
					}
					// Check if Target Set is inside Safe Set
					arma::umat uv = ((SafeSet - TargetSet) > 0);
					arma::umat temp1 = uv;
					arma::umat temp2 = arma::zeros<arma::umat>(SafeSet.n_rows, 1);
					temp2 =
						arma::join_horiz(temp2, arma::ones<arma::umat>(SafeSet.n_rows, 1));
					if (!arma::approx_equal(temp1, temp2, "both", 2, 0.1))
					{
						arma::umat temp3 = (((SafeSet - TargetSet) >= 0));
						if (!arma::approx_equal(temp3, temp2, "both", 2, 0.1))
						{
							std::cout << "ERROR: The Target Set cannot be outside the Safe Set" << std::endl;
							exit(0);
						}
					}
				}

				// Solving the problem
				// Because of taking the center point, the error will be twice as small.
				// This allows to make epsilon twice as large.
				epsilon = 2 * epsilon;
				if (epsilon <= 0)
				{
					std::cout << "ERROR: Incorrect maximum abstraction error, needs to be > 0"
						<< std::endl;
					exit(0);
				}
															// To avoid using list of if-condition
				int task2Solve = 100 * static_cast<int>(Problem) + 10* static_cast<int>(Gridding) + static_cast<int>(Distribution);
				// statements since it is slower
				// For multiple modes store the discreteKernel for each mode in a vector
				std::vector<std::vector<Matrix>> Tp_all;
				for(size_t i = 0; i < input.myModel.discreteModes; i++)
				{
					// Get current model to perform abstraction and task on
					input.myModel.stateSpaceModelsPerMode[0] = input.myModel.stateSpaceModelsPerMode[i];
					taskMDP.myKernel(input.myModel);
					if (controlled)
					{
						switch (task2Solve)
						{
							case 111:
							{
								taskMDP.Uniform_grid(epsilon, timeHorizon, SafeSet);
								break;
							}
							case 121:
							{
								taskMDP.Uniform_grid(10 * epsilon, timeHorizon, SafeSet);
								taskMDP.Adaptive_grid_multicell(epsilon, timeHorizon);
								break;
							}
							case 131:
							{
								taskMDP.Uniform_grid(10 * epsilon, timeHorizon, SafeSet);
								taskMDP.Adaptive_grid_multicell_semilocal(epsilon, timeHorizon, SafeSet);
								break;
							}
							case 211:
							{
								taskMDP.Uniform_grid_ReachAvoid(epsilon, timeHorizon, SafeSet, TargetSet);
								break;
							}
							case 221:
							{
								taskMDP.Uniform_grid_ReachAvoid(10 * epsilon, timeHorizon, SafeSet,
									TargetSet);
								taskMDP.Adaptive_grid_ReachAvoid(epsilon, timeHorizon, SafeSet, TargetSet);
								break;
							}
							case 231:
							{
								taskMDP.Uniform_grid_ReachAvoid(10 * epsilon, timeHorizon, SafeSet,
									TargetSet);
								taskMDP.Adaptive_grid_ReachAvoid_semilocal(epsilon, timeHorizon, SafeSet,
									TargetSet);
								break;
							}
							case 112:
							{
								taskMDP.Uniform_grid_MCapprox(epsilon, timeHorizon, SafeSet);
								break;
							}
							case 212:
							{
								taskMDP.Uniform_grid_ReachAvoid_MCapprox(epsilon, timeHorizon, SafeSet,
									TargetSet);
								break;
							}
							case 122:
							{
								std::cout << "This options is not available yet. Work in progress.";
								taskMDP.Adaptive_grid_MCapprox(epsilon, timeHorizon, SafeSet);
								break;
							}
							case 222:
							{
								std::cout << "This options is not available yet. Work in progress.";
								taskMDP.Adaptive_grid_ReachAvoidMCapprox(epsilon, timeHorizon, SafeSet,
									TargetSet);
								break;
							}
							default:
							{
								std::cout << "This options is not available yet. Work in progress.";
								exit(0);
								break;
							}
						}
						std::cout << "The abstraction consists of " << taskMDP.X.n_rows
							<< " representative points." << std::endl;
						if (taskMDP.X.n_rows > 1000000)
						{
							std::cout << "ERROR: Abstraction is too large, need more memory"
								<< std::endl;
							exit(0);
						}
						// Because of taking the center points as representative points the
						// resulting error is half of the outcome error.
						taskMDP.E = 0.5 * taskMDP.E;

						// Creation of Markov Chain
						if (Distribution == 2)
						{
							taskMDP.MCapprox(epsilon);
						}
						else
						{
							taskMDP.MCcreator(epsilon);
						}
						// Calculation of the resulting problem
						if(input.myModel.discreteModes == 1)
						{
							switch (Problem)
							{
								case VERIFY_SAFETY:
								{
									taskMDP.StandardProbSafety(timeHorizon);
									end = clock();
									time = (double)(end - begin) / CLOCKS_PER_SEC;
									break;
								}
								case VERIFY_REACH_AVOID:
								{
									taskMDP.StandardReachAvoid(TargetSet, timeHorizon);
									end = clock();
									time = (double)(end - begin) / CLOCKS_PER_SEC;
									break;
								}
								default:
								{
									std::cout << "ERROR: This options is not available yet. Work in progress."<< std::endl;
									exit(0);
								}
								break;
							}
						}
						else
						{
							Tp_all.push_back(taskMDP.Tp);
						}
					}
					else
					{
						switch (task2Solve)
						{
							case 111:
							{
								taskMDP.Uniform_grid_Contr(epsilon, timeHorizon, SafeSet, InputSet);
								break;
							}
							case 121:
							{
								taskMDP.Uniform_grid_Contr(10 * epsilon, timeHorizon, SafeSet, InputSet);
								taskMDP.Adaptive_grid_multicell_Contr(epsilon, timeHorizon, SafeSet,
									InputSet);
								break;
							}
							case 131:
							{
								taskMDP.Uniform_grid_Contr(10 * epsilon, timeHorizon, SafeSet, InputSet);
								taskMDP.Adaptive_grid_semilocal_Contr(epsilon, timeHorizon, SafeSet,
									InputSet);
								break;
							}
							case 211:
							{
								taskMDP.Uniform_grid_ReachAvoid_Contr(epsilon, timeHorizon, SafeSet, InputSet,
									TargetSet);
								break;
							}
							case 221:
							{
								taskMDP.Uniform_grid_ReachAvoid_Contr(10 * epsilon, timeHorizon, SafeSet,
									InputSet, TargetSet);
								taskMDP.Adaptive_grid_ReachAvoid(epsilon, timeHorizon, SafeSet, TargetSet);
								break;
							}
							case 231:
							{
								std::cout << "This options is not available yet. Work in progress.";
								exit(0);
								taskMDP.Uniform_grid_ReachAvoid_Contr(10 * epsilon, timeHorizon, SafeSet,
									InputSet, TargetSet);
								//    taskMDP.Adaptive_grid_ReachAvoid_semilocal(epsilon,timeHorizon,SafeSet,TargetSet);
								break;
							}
							case 112:
							{
								taskMDP.Uniform_grid_MCapprox_Contr(epsilon, timeHorizon, SafeSet, InputSet);
								break;
							}
							case 212:
							{
								taskMDP.Uniform_grid_ReachAvoidMCapprox_Contr(epsilon, timeHorizon, SafeSet,
									InputSet, TargetSet);
								break;
							}
							case 222:
							{
								taskMDP.Adaptive_grid_ReachAvoidMCapprox(epsilon, timeHorizon, SafeSet,
									TargetSet);
								break;
							}
							default:
							{
								std::cout << "This options is not available yet. Work in progress.";
								exit(0);
								break;
							}
						}

						std::cout << "INFO: The abstraction consists of " << taskMDP.X.n_rows
							<< " state representative points." << std::endl;
						std::cout << "INFO: The abstraction consists of " << taskMDP.U.n_rows
							<< " input representative points." << std::endl;

						if (taskMDP.X.n_rows * taskMDP.U.n_rows > 100000)
						{
							std::cout << "ERROR: Abstraction is too large, need more memory"
								<< std::endl;
							exit(0);
						}

						// Because of taking the center points as representative points the
						// resulting error is half of the outcome error.
						taskMDP.E = 0.5 * taskMDP.E;
						// Creation of Markov Chain
						if (Distribution == SAMPLE)
						{
							taskMDP.MCapprox_Contr(epsilon, model);
						}
						else
						{
							taskMDP.MCcreator_Contr(epsilon);
						}
						// Calculation of the resulting problem
						if(input.myModel.discreteModes == 1)
						{
							switch (Problem)
							{
								case VERIFY_SAFETY:
								{
									taskMDP.StandardProbSafety_Contr(timeHorizon);
									end = clock();
									time = (double)(end - begin) / CLOCKS_PER_SEC;
									break;
								}
								case VERIFY_REACH_AVOID:
								{
									taskMDP.StandardReachAvoid_Contr(TargetSet, timeHorizon);
									end = clock();
									time = (double)(end - begin) / CLOCKS_PER_SEC;
									break;
								}
								default:
								{
									std::cout << "This options is not available yet. Work in progress.";
									exit(0);
								}
								break;
							}
						}
						else
						{
							Tp_all.push_back(taskMDP.Tp);
						}
					}
				}

				// Get current time to time stamp outputs
				auto t = std::time(nullptr);
				auto tm = *std::localtime(&t);
				std::ostringstream oss;
				oss << std::put_time(&tm, "%d-%m-%Y-%H-%M-%S");
				auto str = oss.str();

				// Rescale grid axis to original axis
				arma::vec inter_d  = arma::ones<arma::vec>(taskMDP.X.n_rows);
				for(unsigned d = 0; d < taskMDP.X.n_cols/2; d++)
				{
					inter_d = j(d,d)*inter_d;
					taskMDP.X.col(d) = inter_d%taskMDP.X.col(d);
					inter_d = arma::ones<arma::vec>(taskMDP.X.n_rows);
				}
				if(input.myModel.discreteModes == 1)
				{
					taskMDP.formatOutput(time, str, Problem, timeHorizon);
				}
			}
			// For multiple modes only now that I have obtained all the
			//  Transition probabilities I need to perform parallel composition
			// and compute overall MDP_ABSTRACTIOtimeHorizon
			// Tp_all contains all the Tp matrices
			// For Safety ---------------------------------------
			if(input.myModel.discreteModes > 1)
			{
				Matrix discreteKernel = input.myModel.discreteKernel;
				// for each mode compute V
				// first get Tp for first control
				std::vector<Matrix> Tpall;
				for(size_t qd = 0; qd < input.myModel.discreteModes; qd++)
				{
					Tpall.push_back(Tp_all[qd][0]);
				}
				// Compute V for first one
				taskMDP.StandardProbSafety(input.myModel.discreteModes, discreteKernel, Tpall, timeHorizon);
				// Where we are storing all the Value functions for each control
				int u_dim = 1;
				if(!taskMDP.U.is_empty())
				{
					u_dim= taskMDP.U.n_rows;
				}
				Matrix Vall(taskMDP.V.n_rows, u_dim);
				// store first computed V for first control in first column
				Vall.col(0) =taskMDP.V;

				// Iteratively compute for all remaining control
				Tpall.clear();
				for(size_t ud = 1; ud < u_dim; ud++)
				{
					for(size_t qd = 0; qd < input.myModel.discreteModes; qd++)
					{
						Tpall.push_back(Tp_all[qd][ud]);
					}
					// Compute the next value function
					taskMDP.StandardProbSafety(input.myModel.discreteModes, discreteKernel, Tpall, timeHorizon);
					// Store result
					Vall.col(ud) = taskMDP.V;
				}
				// timeHorizonow to get final value  and optimal policy need to find
				// V with max and corresponding column index
				Matrix newV = arma::max(Vall, 1);
				taskMDP.V = newV;

				arma::ucolvec index = arma::index_max(Vall, 1);
				taskMDP.OptimalPol = index;

				//-------------------------------------------------------
				auto t = std::time(nullptr);
				auto tm = *std::localtime(&t);
				std::ostringstream oss;
				oss << std::put_time(&tm, "%d-%m-%Y-%H-%M-%S");
				auto str = oss.str();
				taskMDP.formatOutput(time, str, Problem, timeHorizon);
			}
		}
		catch (const std::bad_alloc &e)
		{
			std::cerr << e.what() << std::endl;
			exit(0);
		}
		break;
	}
	case IMDP_ABSTRACTION									// model checking using BMDP
	{
		// Initialise timers
		clock_t begin, end;
		double time;
		begin = clock();

		// Get model definitions
		input.myModel.stateSpaceModelsPerMode[0].F = arma::eye(input.myModel.continuousSpaceDimension, input.myModel.continuousSpaceDimension);
		bmdp_t taskBMDP(input.myTask, input.myModel);
		// Check input model
		int num_cont = input.myModel.continuousSpaceDimension;
		for (int i = 0; i < input.myModel.discreteModes; ++i)
		{
			if (input.myModel.stateSpaceModelsPerMode[i].A.n_rows != num_cont)
			{
				std::cout << "ERROR: Different number of continuous variables for each mode, "
					"currently not supported. Work in progress"
					<< std::endl;
				exit(0);
			}
			if (input.myModel.stateSpaceModelsPerMode[i].Sigma.n_rows != num_cont)
			{
				std::cout << "ERROR: Different number of continuous variables for each mode, "
					"currently not supported. Work in progress"
					<< std::endl;
				exit(0);
			}
		}
		// Set number of actions
		taskBMDP.N = input.myModel.discreteModes;			// TODO(ncauchi) check
		for (int i = 0; i < input.myModel.discreteModes; ++i)
		{
			if (i == 0)
			{
				Matrix j = ZeroMatrix(1, input.myModel.continuousSpaceDimension);
				taskBMDP.desc.mean[0] = j;
				taskBMDP.desc.cov[0] = taskBMDP.desc.dyn.dynamics[i].Sigma;
			}
			else
			{
				Matrix j = ZeroMatrix(1, input.myModel.continuousSpaceDimension);
				taskBMDP.desc.mean.push_back(j);
				taskBMDP.desc.cov.push_back(taskBMDP.desc.dyn.dynamics[i].Sigma);
			}
		}
		taskBMDP.desc.boundary = input.myTask.boundary;
		taskBMDP.desc.gridsize = input.myTask.gridsize;
		if (taskBMDP.desc.boundary.n_rows != taskBMDP.desc.gridsize.n_cols)
		{
			std::cout << "ERROR: Incorrect boundary or grid size" << std::endl;
			exit(0);
		}
		taskBMDP.desc.reftol = input.myTask.reftol;
		if (taskBMDP.desc.reftol.n_cols != taskBMDP.desc.gridsize.n_cols)
		{
			std::cout << "ERROR: Incorrect relative tolerance parameter input" << std::endl;
			exit(0);
		}
		int RA = 0;
		if (input.myTask.propertySpec == VERIFY_REACH_AVOID ||input.myTask.propertySpec == SYNTHESIS_REACH_AVOID )
		{
			RA = 1;
		}
		// Check if need to rescale boundary
		Matrix proj = Identity(taskBMDP.desc.boundary.n_rows, taskBMDP.desc.boundary.n_rows);

		taskBMDP.eps = input.myTask.eps;
		// Constructing the abstraction
		std::cout << "INFO: Constructing the IMDP abstraction " << std::endl;

		// Identify if performing synthesis or verification
		// and whether safety or reach avoid
		taskBMDP.bmdpAbstraction(input.myTask.timeHorizon, RA);
		std::cout << "INFO: Done with abstraction construction" << std::endl;
		if (input.myTask.propertySpec == VERIFY_SAFETY || input.myTask.propertySpec == SYNTHESIS_SAFETY)
		{
			arma::uvec phi1 = arma::ones<arma::uvec>(taskBMDP.Stepsmax.n_cols);
			phi1(phi1.n_rows - 1, 0) = 0;
			arma::uvec labels = arma::zeros<arma::uvec>(taskBMDP.Stepsmax.n_cols);
			labels(labels.n_rows - 1, 0) = 1;
			std::cout << "Performing  model checking " << std::endl;
			taskBMDP.createSynthFile(phi1, labels);
			if (input.myTask.propertySpec == VERIFY_SAFETY)
			{
				taskBMDP.runSafety(1e-4, input.myTask.timeHorizon);
			}
			else
			{
				taskBMDP.runSynthesis(1e-4, input.myTask.timeHorizon);
			}
		}
		else
		{
			// Get coordinates of labels for phi1 and phi 2
			arma::uvec phi1 =
				taskBMDP.getLabels("phi1.txt", input.myModel.continuousSpaceDimension, true);

			// Since not phi1 need to negate phi1
			for (unsigned i = 0; i < phi1.n_rows; i++)
			{
				if (phi1(i) == 1)
				{
					phi1(i) = 0;
				}
				else
				{
					phi1(i) = 1;
				}
			}

			arma::uvec phi2 =
				taskBMDP.getLabels("phi2.txt", input.myModel.continuousSpaceDimension, false);

			std::cout << "INFO: Performing  model checking " << std::endl;
			taskBMDP.createSynthFile(phi1, phi2);

			if (input.myTask.propertySpec == 2)
			{
				taskBMDP.runSafety(1e-4, input.myTask.timeHorizon);
			}
			else
			{
				taskBMDP.runSynthesis(1e-4, input.myTask.timeHorizon);
			}
		}

		end = clock();
		time = (double)(end - begin) / CLOCKS_PER_SEC;
		auto t = std::time(nullptr);
		auto tm = *std::localtime(&t);
		std::ofstream myfile;
		std::ostringstream oss;
		oss << std::put_time(&tm, "%d-%m-%Y-%H-%M-%S");
		auto str = oss.str();

		taskBMDP.formatOutput(time, str);

		// Simulation of BMDP
		if (input.myTask.propertySpec == 3 && input.myModel.continuousSpaceDimension)
		{
			std::string exportOpt;
			std::string str("y");
			std::string str1("yes");
			std::string str2("YES");
			std::cout << "Would you like to simulate evolution under synthesised "
				"policy [y- yes, n - no]"
				<< std::endl;
			std::cin >> exportOpt;
			if ((exportOpt.compare(str) == 0) || (exportOpt.compare(str1) == 0)|| (exportOpt.compare(str2) == 0))
			{
				std::cout << "Simulation of model using optimal policy" << std::endl;
				Matrix init(2, 1);
				init(0, 0) = -0.5;
				init(1, 0) = -1;
				int mode = 1;
				Matrix policy = taskBMDP.Policy;

				// Find IMDP state you are in; we start from mode 1
				int init_state = 0;
				std::random_device rand_dev;
				std::mt19937 generator(rand_dev());
				std::normal_distribution<double> d2{0, 1};
				for (unsigned i = 0; i < taskBMDP.mode[0].vertices.size(); i++)
				{
					Matrix crnt_vertex = taskBMDP.mode[0].vertices[i];
					int inpoly = pnpoly(1, crnt_vertex.row(0), crnt_vertex.row(1),
						init(0, 0), init(1, 0));
					if (inpoly)
					{
						init_state = i;
						break;
					}
				}

				Matrix prev(2, 10);
				prev.col(0) = init;
				Matrix W(2, 1);
				for (unsigned t = 0; t < 9; t++)
				{
					for (unsigned j = 0; j < 2; j++)
					{
						W(j, 0) = d2(generator);
					}
					Matrix x_new =
						input.myModel.stateSpaceModelsPerMode[policy(init_state, 0)].A * prev.col(t) +
						input.myModel.stateSpaceModelsPerMode[policy(init_state, 0)].Sigma * W;
					if (x_new(0, 0) > arma::max(taskBMDP.desc.boundary.row(0)))
					{
						x_new(0, 0) = arma::max(taskBMDP.desc.boundary.row(0));
					}
					if (x_new(1, 0) > arma::max(taskBMDP.desc.boundary.row(1)))
					{
						x_new(1, 0) = arma::max(taskBMDP.desc.boundary.row(1));
					}
					if (x_new(1, 0) < arma::min(taskBMDP.desc.boundary.row(1)))
					{
						x_new(1, 0) = arma::min(taskBMDP.desc.boundary.row(1));
					}
					if (x_new(0, 0) < arma::min(taskBMDP.desc.boundary.row(0)))
					{
						x_new(0, 0) = arma::min(taskBMDP.desc.boundary.row(0));
					}
					for (unsigned i = 0;
						i < taskBMDP.mode[policy(init_state, 0)].vertices.size(); i++)
					{
						Matrix crnt_vertex =
							taskBMDP.mode[policy(init_state, 0)].vertices[i];
						int inpoly = pnpoly(1, crnt_vertex.row(0), crnt_vertex.row(1),
							x_new(0, 0), x_new(1, 0));
						if (inpoly)
						{
							init_state = i;
							break;
						}
					}
					prev.col(t + 1) = x_new;
				}
				std::string sim_name = "../results/Simulation_" + str + ".txt";
				prev = prev.t();
				prev.save(sim_name, arma::raw_ascii);
			}
			else
			{
				std::cout << "Done!" << std::endl;
			}
		}
		break;
	}
	default:
	{
		std::cout << "timeHorizono task specification given" << std::endl;
		exit(0);
	}
}


}
#endif
