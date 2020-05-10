/*
 * taskExec.h
 *
 *  Created on: 19 Feb 2018
 *      Author: nathalie
 */

#ifndef TASKEXEC_H_
#define TASKEXEC_H_

#include "Bmdp.h"
#include "FAUST.h"
#include "InputSpec.h"

static void performTask(inputSpec_t<arma::mat, int> input)
{
  switch (input.myTask.task)
  {
    case SIMULATION: // Perform simulation depending on model type
    {
      if (input.myTask.T <= 0)
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
        if (input.myModel.stateSpaceModelsPerMode[i].sigma.n_rows != num_cont) {
          std::cout << "ERROR: Different number of continuous variables for each mode, "
                       "currently not supported. Work in progress"
                    << std::endl;
          exit(0);
        }
      }
      if (input.myModel.discreteKernel.n_rows != input.myModel.stateSpaceModelsPerMode.size()) {
        std::cout << "ERROR: Incorrect transition matrix. Need to have number of rows "
                     "== number of columns == number of models representing the "
                     "evolution of the continuous variables."
                  << std::endl;
        exit(0);
      }
      input.myModel.run(input.myModel, input.myTask.T, input.myTask.runs);


      break;
    }
    case MDP_ABSTRACTION: // Perform verification task depending on model and tool to be
            // interfaced with
    {
      faust_t myF;
      shs_t<arma::mat, int> myModel = input.myModel;
      myF.model = myModel;
      try {
        // Check if dynamics are ill-scaled
        // Obtaining problem definition
        clock_t begin, end;
        begin = clock();
        double time = 0;

        int Problem = 100 * input.myTask.propertySpec;
        int Gridding = 10 * input.myTask.typeGrid;
        int Distribution = input.myTask.assumptionsKernel;
        int Control = input.myTask.Controlled;


        // Problem variables
        double epsilon = input.myTask.eps;
        int N = input.myTask.T;
        shs_t<arma::mat, int> model = input.myModel;

        // Get Safe, Input and Target set
        arma::mat SafeSet = input.myTask.safeSet;
        arma::mat InputSet = input.myTask.inputSet;
        arma::mat TargetSet = input.myTask.targetSet;

        faust_t taskFAUST = myF;
        // Check correctness of SafeSet input
        if (SafeSet.n_cols != 2) {
          std::cout << "There is no correctly defined Safe Set";
          exit(0);
        }
        double temp = arma::accu((SafeSet.col(1) - SafeSet.col(0)) < 0);
        if (temp > 0) {
          std::cout << "The edges of the Safe Set must have positive length. "
                       "Make sure that the first column is the lower bound and "
                       "the second column is the upper bound";
          exit(0);
        }
        // Check if Safe Set needs rescaling
        // Check if need to rescale boundary
        arma::mat max_Ss = arma::max(SafeSet);
        arma::mat min_Ss = arma::min(SafeSet);
        double diff_Ss = max_Ss(1) - min_Ss(0);
        arma::mat j = arma::eye<arma::mat>(SafeSet.n_rows, SafeSet.n_rows);
        if(diff_Ss > 50) {
          // Identify row with smallest values to rescale to that
          arma::mat min_row = arma::min(SafeSet,0);
          arma::mat OrSS = SafeSet ; // The original Safeset definition
          arma::mat to_inv = OrSS;
          for(unsigned i =0; i < OrSS.n_rows; i++) {
            if(min_row(0,1) != SafeSet(i,1)) {
              SafeSet.row(i) = min_row;
            }
          }
          arma::mat y_inv = arma::diagmat(min_row);
          arma::mat j_bar = OrSS/y_inv;
          j = j_bar.replace(arma::datum::inf, 0);
          for(unsigned k = 0; k < model.stateSpaceModelsPerMode.size(); k++) {
            taskFAUST.model.stateSpaceModelsPerMode[k].sigma = arma::inv(j)*model.stateSpaceModelsPerMode[k].sigma;
          }
        }
        if (Control == 0) {
          InputSet.clear();
        } else {
          if (model.stateSpaceModelsPerMode[0].B.n_cols != InputSet.n_rows) {
            std::cout << "There is no correctly defined Input Set";
            exit(0);
          }
        }

        // Check if dimensions are correct
        if ((unsigned)model.stateSpaceModelsPerMode[0].A.n_cols != SafeSet.n_rows) {
          std::cout << "The dimension of the Kernel does not match the "
                       "dimensions of the Safe Set";
          exit(0);
        }

        // Check correctness of target set if problem of reach and avoid
        if (Problem == 2) {
          if (TargetSet.n_cols != 2) {
            std::cout << "There is no correctly defined Target Set";
            exit(0);
          }
          double temp = arma::accu((TargetSet.col(1) - TargetSet.col(0)) < 0);
          if (temp > 0) {
            std::cout << "The edges of the Target Set must have positive length. "
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
          if (!arma::approx_equal(temp1, temp2, "both", 2, 0.1)) {
            arma::umat temp3 = (((SafeSet - TargetSet) >= 0));
            if (!arma::approx_equal(temp3, temp2, "both", 2, 0.1)) {
              std::cout << "The Target Set cannot be outside the Safe Set";
              exit(0);
            }
          }
        }

        // Solving the problem
        // Because of taking the center point, the error will be twice as small.
        // This allows to make epsilon twice as large.
        epsilon = 2 * epsilon;
        if (epsilon <= 0) {
          std::cout << "Incorrect maximum abstraction error, needs to be > 0"
                    << std::endl;
          exit(0);
        }
        int task2Solve = Problem + Gridding +
                         Distribution; // To avoid using list of if-condition
                                       // statements since it is slower
        // For multiple modes store the discreteKernel for each mode in a vector
        std::vector<std::vector<arma::mat>> Tp_all;
        for(size_t i = 0; i < input.myModel.discreteModes; i++) {
          // Get current model to perform abstraction and task on
          input.myModel.stateSpaceModelsPerMode[0] = input.myModel.stateSpaceModelsPerMode[i];
          taskFAUST.myKernel(input.myModel);
          switch (Control) {
          case 0: {
            switch (task2Solve) {
            case 111: {
              taskFAUST.Uniform_grid(epsilon, N, SafeSet);
              break;
            }
            case 121: {
              taskFAUST.Uniform_grid(10 * epsilon, N, SafeSet);
              taskFAUST.Adaptive_grid_multicell(epsilon, N);
              break;
            }
            case 131: {
              taskFAUST.Uniform_grid(10 * epsilon, N, SafeSet);
              taskFAUST.Adaptive_grid_multicell_semilocal(epsilon, N, SafeSet);
              break;
            }
            case 211: {
              taskFAUST.Uniform_grid_ReachAvoid(epsilon, N, SafeSet, TargetSet);
              break;
            }
            case 221: {
              taskFAUST.Uniform_grid_ReachAvoid(10 * epsilon, N, SafeSet,
                                                TargetSet);
              taskFAUST.Adaptive_grid_ReachAvoid(epsilon, N, SafeSet, TargetSet);
              break;
            }
            case 231: {
              taskFAUST.Uniform_grid_ReachAvoid(10 * epsilon, N, SafeSet,
                                                TargetSet);
              taskFAUST.Adaptive_grid_ReachAvoid_semilocal(epsilon, N, SafeSet,
                                                           TargetSet);
              break;
            }
            case 112: {
              taskFAUST.Uniform_grid_MCapprox(epsilon, N, SafeSet);
              break;
            }
            case 212: {
              taskFAUST.Uniform_grid_ReachAvoid_MCapprox(epsilon, N, SafeSet,
                                                         TargetSet);
              break;
            }
            case 122: {
              std::cout << "This options is not available yet. Work in progress.";
              taskFAUST.Adaptive_grid_MCapprox(epsilon, N, SafeSet);
              break;
            }
            case 222: {
              std::cout << "This options is not available yet. Work in progress.";
              taskFAUST.Adaptive_grid_ReachAvoidMCapprox(epsilon, N, SafeSet,
                                                         TargetSet);
              break;
            }
            default: {
              std::cout << "This options is not available yet. Work in progress.";
              exit(0);
              break;
            }
            }
            std::cout << "The abstraction consists of " << taskFAUST.X.n_rows
                      << " representative points." << std::endl;
            if (taskFAUST.X.n_rows > 1000000) {
              std::cout << "Abstraction is too large, need more memory"
                        << std::endl;
              exit(0);
            }
            // Because of taking the center points as representative points the
            // resulting error is half of the outcome error.
            taskFAUST.E = 0.5 * taskFAUST.E;

    	       // Creation of Markov Chain
            if (Distribution == 2) {
              taskFAUST.MCapprox(epsilon);
            } else {
              taskFAUST.MCcreator(epsilon);
            }
            // Calculation of the resulting problem
            if(input.myModel.discreteModes == 1) {
              switch (Problem) {
                case 100: {
                  taskFAUST.StandardProbSafety(N);
                  end = clock();
                  time = (double)(end - begin) / CLOCKS_PER_SEC;
                  break;
              }
              case 200: {
                taskFAUST.StandardReachAvoid(TargetSet, N);
                end = clock();
                time = (double)(end - begin) / CLOCKS_PER_SEC;
                break;
               }
               default: {
                 std::cout << "This options is not available yet. Work in progress.";
                 exit(0);
               }
               break;
             }
            }
            else {
              Tp_all.push_back(taskFAUST.Tp);
            }
          //  }
            break;
          }
          case 1: {
            switch (task2Solve) {
            case 111: {
              taskFAUST.Uniform_grid_Contr(epsilon, N, SafeSet, InputSet);
              break;
            }
            case 121: {
              taskFAUST.Uniform_grid_Contr(10 * epsilon, N, SafeSet, InputSet);
              taskFAUST.Adaptive_grid_multicell_Contr(epsilon, N, SafeSet,
                                                      InputSet);
              break;
            }
            case 131: {
              taskFAUST.Uniform_grid_Contr(10 * epsilon, N, SafeSet, InputSet);
              taskFAUST.Adaptive_grid_semilocal_Contr(epsilon, N, SafeSet,
                                                      InputSet);
              break;
            }
            case 211: {
              taskFAUST.Uniform_grid_ReachAvoid_Contr(epsilon, N, SafeSet, InputSet,
                                                      TargetSet);
              break;
            }
            case 221: {
              taskFAUST.Uniform_grid_ReachAvoid_Contr(10 * epsilon, N, SafeSet,
                                                      InputSet, TargetSet);
              taskFAUST.Adaptive_grid_ReachAvoid(epsilon, N, SafeSet, TargetSet);
              break;
            }
            case 231: {
              std::cout << "This options is not available yet. Work in progress.";
              exit(0);
              taskFAUST.Uniform_grid_ReachAvoid_Contr(10 * epsilon, N, SafeSet,
                                                      InputSet, TargetSet);
              //    taskFAUST.Adaptive_grid_ReachAvoid_semilocal(epsilon,N,SafeSet,TargetSet);
              break;
            }
            case 112: {
              taskFAUST.Uniform_grid_MCapprox_Contr(epsilon, N, SafeSet, InputSet);
              break;
            }
            case 212: {
              taskFAUST.Uniform_grid_ReachAvoidMCapprox_Contr(epsilon, N, SafeSet,
                                                              InputSet, TargetSet);
              break;
            }
            case 222: {
              taskFAUST.Adaptive_grid_ReachAvoidMCapprox(epsilon, N, SafeSet,
                                                         TargetSet);
              break;
            }
            default: {
              std::cout << "This options is not available yet. Work in progress.";
              exit(0);
              break;
            }
            }
            std::cout << "The abstraction consists of " << taskFAUST.X.n_rows
                      << " state representative points." << std::endl;
            std::cout << "The abstraction consists of " << taskFAUST.U.n_rows
                      << " input representative points." << std::endl;
            if (taskFAUST.X.n_rows * taskFAUST.U.n_rows > 100000) {
              std::cout << "Abstraction is too large, need more memory"
                        << std::endl;
              exit(0);
            }
            // Because of taking the center points as representative points the
            // resulting error is half of the outcome error.
            taskFAUST.E = 0.5 * taskFAUST.E;
            // Creation of Markov Chain
            if (Distribution == 2) {
              taskFAUST.MCapprox_Contr(epsilon, model);
            } else {
              taskFAUST.MCcreator_Contr(epsilon);
            }
            // Calculation of the resulting problem
            if(input.myModel.discreteModes == 1) {
              switch (Problem) {
                case 100: {
                  taskFAUST.StandardProbSafety_Contr(N);
                  end = clock();
                  time = (double)(end - begin) / CLOCKS_PER_SEC;
                break;
              }
              case 200: {
                taskFAUST.StandardReachAvoid_Contr(TargetSet, N);
                end = clock();
    	          time = (double)(end - begin) / CLOCKS_PER_SEC;
    	           break;
               }
               default: {
                 std::cout << "This options is not available yet. Work in progress.";
                 exit(0);
               }
               break;
             }
            }
            else {
              Tp_all.push_back(taskFAUST.Tp);
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
         arma::vec inter_d  = arma::ones<arma::vec>(taskFAUST.X.n_rows);
         for(unsigned d = 0; d < taskFAUST.X.n_cols/2; d++) {
           inter_d = j(d,d)*inter_d;
           taskFAUST.X.col(d) = inter_d%taskFAUST.X.col(d);
           inter_d = arma::ones<arma::vec>(taskFAUST.X.n_rows);
         }
         if(input.myModel.discreteModes == 1) {
           taskFAUST.formatOutput(time, str, Problem, N);
         }
        }
        // For multiple modes only now that I have obtained all the
        //  Transition probabilities I need to perform parallel composition
        // and compute overall MDP_ABSTRACTION
        // Tp_all contains all the Tp matrices
        // For Safety ---------------------------------------
        // TODO Check > 2D for safety + REACH AVOID
        if(input.myModel.discreteModes > 1) {
          arma::mat discreteKernel = input.myModel.discreteKernel;
          // for each mode compute V
          // first get Tp for first control
          std::vector<arma::mat> Tpall;
          for(size_t qd = 0; qd < input.myModel.discreteModes; qd++) {
            Tpall.push_back(Tp_all[qd][0]);
          }
          // Compute V for first one
          taskFAUST.StandardProbSafety(input.myModel.discreteModes, discreteKernel, Tpall, N);
          // Where we are storing all the Value functions for each control
          int u_dim = 1;
          if(!taskFAUST.U.is_empty()) {
            u_dim= taskFAUST.U.n_rows;
          }
          arma::mat Vall(taskFAUST.V.n_rows, u_dim);
          // store first computed V for first control in first column
          Vall.col(0) =taskFAUST.V;

          // Iteratively compute for all remaining control
          Tpall.clear();
          for(size_t ud = 1; ud < u_dim; ud++) {
            for(size_t qd = 0; qd < input.myModel.discreteModes; qd++) {
              Tpall.push_back(Tp_all[qd][ud]);
            }
            // Compute the next value function
            taskFAUST.StandardProbSafety(input.myModel.discreteModes, discreteKernel, Tpall, N);
            // Store result
            Vall.col(ud) = taskFAUST.V;
          }
          // Now to get final value  and optimal policy need to find
          // V with max and corresponding column index
          arma::mat newV = arma::max(Vall, 1);
          taskFAUST.V = newV;

          arma::ucolvec index = arma::index_max(Vall, 1);
          taskFAUST.OptimalPol = index;

          //-------------------------------------------------------
          auto t = std::time(nullptr);
          auto tm = *std::localtime(&t);
          std::ostringstream oss;
          oss << std::put_time(&tm, "%d-%m-%Y-%H-%M-%S");
          auto str = oss.str();
          taskFAUST.formatOutput(time, str, Problem, N);
        }
      } catch (const std::bad_alloc &e) {
        std::cerr << e.what() << std::endl;
        exit(0);
      }
      break;
    }
    case IMDP_ABSTRACTION: // model checking using BMDP
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
      for (int i = 0; i < input.myModel.discreteModes; ++i) {
        if (input.myModel.stateSpaceModelsPerMode[i].A.n_rows != num_cont) {
          std::cout << "Different number of continuous variables for each mode, "
                       "currently not supported. Work in progress"
                    << std::endl;
          exit(0);
        }
        if (input.myModel.stateSpaceModelsPerMode[i].sigma.n_rows != num_cont) {
          std::cout << "Different number of continuous variables for each mode, "
                       "currently not supported. Work in progress"
                    << std::endl;
          exit(0);
        }
      }
      // Set number of actions
      taskBMDP.actNum = input.myModel.discreteModes;
      for (int i = 0; i < input.myModel.discreteModes; ++i) {
        if (i == 0) {
          arma::mat j = arma::zeros<arma::mat>(1, input.myModel.continuousSpaceDimension);
          taskBMDP.desc.mean[0] = j;
          taskBMDP.desc.cov[0] = taskBMDP.desc.dyn.dynamics[i].sigma;
        } else {
          arma::mat j = arma::zeros<arma::mat>(1, input.myModel.continuousSpaceDimension);
          taskBMDP.desc.mean.push_back(j);
          taskBMDP.desc.cov.push_back(taskBMDP.desc.dyn.dynamics[i].sigma);
        }
      }
      taskBMDP.desc.boundary = input.myTask.boundary;
      taskBMDP.desc.gridsize = input.myTask.gridsize;
      if (taskBMDP.desc.boundary.n_rows != taskBMDP.desc.gridsize.n_cols) {
        std::cout << "Incorrect boundary or grid size" << std::endl;
        exit(0);
      }
      taskBMDP.desc.reftol = input.myTask.reftol;
      if (taskBMDP.desc.reftol.n_cols != taskBMDP.desc.gridsize.n_cols) {
        std::cout << "Incorrect relative tolerance parameter input" << std::endl;
        exit(0);
      }
      int RA = 0;
      if (input.myTask.propertySpec == 2 ||input.myTask.propertySpec == 4 ) {
        RA = 1;
      }
      // Check if need to rescale boundary
      arma::mat proj = arma::eye<arma::mat>(taskBMDP.desc.boundary.n_rows, taskBMDP.desc.boundary.n_rows);

      taskBMDP.eps = input.myTask.eps;
      // Constructing the abstraction
      std::cout << "Constructing the IMDP abstraction " << std::endl;

      // Identify if performing synthesis or verification
      // and whether safety or reach avoid
      taskBMDP.bmdpAbstraction(input.myTask.T, RA);
      std::cout << "Done with abstraction construction" << std::endl;
      if (input.myTask.propertySpec == 1 || input.myTask.propertySpec == 3) {
        arma::uvec phi1 = arma::ones<arma::uvec>(taskBMDP.Stepsmax.n_cols);
        phi1(phi1.n_rows - 1, 0) = 0;
        arma::uvec labels = arma::zeros<arma::uvec>(taskBMDP.Stepsmax.n_cols);
        labels(labels.n_rows - 1, 0) = 1;
        std::cout << "Performing  model checking " << std::endl;
        taskBMDP.createSynthFile(phi1, labels);
        if (input.myTask.propertySpec == 1) {
            taskBMDP.runSafety(1e-4, input.myTask.T);
        } else {
          taskBMDP.runSynthesis(1e-4, input.myTask.T);
        }
      } else {
        // Get coordinates of labels for phi1 and phi 2
        arma::uvec phi1 =
            taskBMDP.getLabels("phi1.txt", input.myModel.continuousSpaceDimension, true);

        // Since not phi1 need to negate phi1
        for (unsigned i = 0; i < phi1.n_rows; i++) {
          if (phi1(i) == 1) {
            phi1(i) = 0;
          } else {
            phi1(i) = 1;
          }
        }

        arma::uvec phi2 =
            taskBMDP.getLabels("phi2.txt", input.myModel.continuousSpaceDimension, false);

        std::cout << "Performing  model checking " << std::endl;
        taskBMDP.createSynthFile(phi1, phi2);

        if (input.myTask.propertySpec == 2) {
          taskBMDP.runSafety(1e-4, input.myTask.T);
        } else {
          taskBMDP.runSynthesis(1e-4, input.myTask.T);
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
      if (input.myTask.propertySpec == 3 && input.myModel.continuousSpaceDimension) {
        std::string exportOpt;
        std::string str("y");
        std::string str1("yes");
        std::string str2("YES");
        std::cout << "Would you like to simulate evolution under synthesised "
                     "policy [y- yes, n - no]"
                  << std::endl;
        std::cin >> exportOpt;
        if ((exportOpt.compare(str) == 0) || (exportOpt.compare(str1) == 0)|| (exportOpt.compare(str2) == 0)) {
          std::cout << "Simulation of model using optimal policy" << std::endl;
          arma::mat init(2, 1);
          init(0, 0) = -0.5;
          init(1, 0) = -1;
          int mode = 1;
          arma::mat policy = taskBMDP.Policy;

          // Find IMDP state you are in; we start from mode 1
          int init_state = 0;
          std::random_device rand_dev;
          std::mt19937 generator(rand_dev());
          std::normal_distribution<double> d2{0, 1};
          for (unsigned i = 0; i < taskBMDP.mode[0].vertices.size(); i++) {
            arma::mat crnt_vertex = taskBMDP.mode[0].vertices[i];
            int inpoly = pnpoly(1, crnt_vertex.row(0), crnt_vertex.row(1),
                                init(0, 0), init(1, 0));
            if (inpoly) {
              init_state = i;
              break;
            }
          }

          arma::mat prev(2, 10);
          prev.col(0) = init;
          arma::mat W(2, 1);
          for (unsigned t = 0; t < 9; t++) {
            for (unsigned j = 0; j < 2; j++) {
              W(j, 0) = d2(generator);
            }
            arma::mat x_new =
                input.myModel.stateSpaceModelsPerMode[policy(init_state, 0)].A * prev.col(t) +
                input.myModel.stateSpaceModelsPerMode[policy(init_state, 0)].sigma * W;
            if (x_new(0, 0) > arma::max(taskBMDP.desc.boundary.row(0))) {
              x_new(0, 0) = arma::max(taskBMDP.desc.boundary.row(0));
            }
            if (x_new(1, 0) > arma::max(taskBMDP.desc.boundary.row(1))) {
              x_new(1, 0) = arma::max(taskBMDP.desc.boundary.row(1));
            }
            if (x_new(1, 0) < arma::min(taskBMDP.desc.boundary.row(1))) {
              x_new(1, 0) = arma::min(taskBMDP.desc.boundary.row(1));
            }
            if (x_new(0, 0) < arma::min(taskBMDP.desc.boundary.row(0))) {
              x_new(0, 0) = arma::min(taskBMDP.desc.boundary.row(0));
            }
            for (unsigned i = 0;
                 i < taskBMDP.mode[policy(init_state, 0)].vertices.size(); i++) {
              arma::mat crnt_vertex =
                  taskBMDP.mode[policy(init_state, 0)].vertices[i];
              int inpoly = pnpoly(1, crnt_vertex.row(0), crnt_vertex.row(1),
                                  x_new(0, 0), x_new(1, 0));
              if (inpoly) {
                init_state = i;
                break;
              }
            }
            prev.col(t + 1) = x_new;
          }
          std::string sim_name = "../results/Simulation_" + str + ".txt";
          prev = prev.t();
          prev.save(sim_name, arma::raw_ascii);
        } else {
          std::cout << "Done!" << std::endl;
        }
      }
      break;
    }
    default: {
      std::cout << "No task specification given" << std::endl;
      exit(0);
    }
    }
}

#endif
