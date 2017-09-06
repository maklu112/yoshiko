/*
 * CoreAlgorithm.cpp
 *
 *  Created on: Aug 11, 2017
 *      Author: philipp
 */

#include "CoreAlgorithm.h"

using namespace std;
using namespace lemon;

namespace ysk {

	ClusterEditingSolutions* CoreAlgorithm::run()
	{

		//TODO: Add more 'breakpoints' where the task checks if it is already cancelled

		if (verbosity >= 5)
			cout << _instance << endl;
	
		// warn user about permanent and forbidden edges
		const WorkingCopyGraph g = _instance->getWorkingCopyInstance().getGraph();
		if (verbosity > 1) {
			cout << endl
					<< "Interpreting the following edges as forbidden/permanent: "
					<< endl;
			for (WorkingCopyGraph::EdgeIt e(g); e != INVALID; ++e) {
				if (_instance->getWorkingCopyInstance().isForbidden(e))
					cout << _instance->getEdgeName(e) << "\tforbidden" << endl;
	
				if (_instance->getWorkingCopyInstance().isPermanent(e))
					cout << _instance->getEdgeName(e) << "\tpermanent" << endl;
			}
		}
	
		bool conserveMultipleSolutions = false;
	
		if (_parameter.nrOptimalSolutions > 1) {
			conserveMultipleSolutions = true;
		}
	
		if (verbosity > 1) {
			cout << endl << "applying FPT reduction rules..." << endl;
			cout << "number of nodes:\t"
					<< countNodes(_instance->getWorkingCopyInstance().getGraph())
					<< endl;
		}
	
		//Initialize the flags for the solution
		SolutionFlags flags;


		//Apply reduction rules
		bitset<NUMBER_OF_REDUCTION_RULES> rules(_parameter.rulesBitMask);
	
		ClusterEditingReduction cer(
				rules,
				_parameter.multiplicativeFactor,
				conserveMultipleSolutions
				);
		cer.perform(*_instance);
	
		vector<ClusterReductionInstance*>& reduced = cer.getInstances();
	
		if (verbosity > 1) {
			cout << "=========================" << endl;
			cout << "FPT reduction rules applied exhaustively." << endl;
			cout << "time:\t" << clk << endl;
		}
		flags.totalCost += cer.getTotalCost();
	
		if (verbosity > 1) {
			cout << "total cost:\t" << cer.getTotalCost() << endl;
			cout << "number of instances:\t" << reduced.size() << endl;
	
			cout << endl << "==================================" << endl
					<< "==================================" << endl;
			cout << endl << "solving (reduced) instances..." << endl;
			cout << "-------------------------" << endl;
		}
	

		int j = 0;
	
		long totalNumberOfSolutions = 1;
	
		vector<vector<vector<vector<int> > > > instances;
	
		for (vector<ClusterReductionInstance*>::iterator it = reduced.begin();
				it != reduced.end(); it++, j++) {
			if (verbosity > 1) {
				cout << "solving instance 'no " << j << "'..." << endl;
			}
	
			ClusterEditingInstance& i = *(*it)->getInstance();
	
			ClusterEditingSolutions s;

			long numberOfSolutions = 1;
	
			//Solve the remaining and reduced instance either with the heuristic or the ILP
			if (!_parameter.useHeuristic) {

				//ILP
				_solver = ILPSolver(
						_parameter.separateTriangles,
						_parameter.separatePartitionCuts,
						_parameter.nrOptimalSolutions
				);

				if (_useInformer){
					_solver.registerInformer(_informer);
				}

				_solverActive = true;
				try {
					numberOfSolutions = _solver.solve(i, s, flags);
				} catch (IloException &e) {
					cout << "CPLEX error: " << e.getMessage() << endl;
				}
				_solver.terminate();
				_solverActive = false;
	
			} else {
				//HEURISTIC
				ParameterizedInstance parameterizedInstance(
						i.getWorkingCopyInstance(),
						numeric_limits<double>::signaling_NaN());
				parameterizedInstance.init();
				InducedCostsHeuristic h(parameterizedInstance);
				h.start();
				flags.totalCost += h.getSolution(s);
			}
	
			if (verbosity > 1)
				cout << "time:\t" << clk << endl;
	
			totalNumberOfSolutions = totalNumberOfSolutions * numberOfSolutions;

			//expand solutions: replace merged nodes by cluster
			vector<vector<vector<int> > > partitions;
			expandSolutions(i, s, partitions);

			instances.push_back(partitions);
	
			if (verbosity > 1)
				cout << endl << "==================================" << endl
						<< endl;
	
			delete (*it)->getInstance();
		}
		if (verbosity > 1) {
			cout << "all instances solved." << endl;
			cout << "total number of solutions:\t" << totalNumberOfSolutions
					<< endl;
			cout << "total cost (data reduction + ILP or heuristic):\t";
		}
		if (verbosity > 0)
			cout << flags.totalCost << endl;
	

		_result->resize(totalNumberOfSolutions);
	
		size_t k = 0;
		vector<vector<int> > p;
		mergeSolutions(0, k, p, *_result, instances);

		//Restore timeout flag
		_result->setFlags(flags);
		return _result;
	}
	
	void CoreAlgorithm::registerCplexInformer(yskLib::CplexInformer* informer){
		if (verbosity > 2) cout << "Registering CplexInformer @ CoreAlgorithm" << endl;
		_informer = informer;
		_useInformer = true;
	}

	void CoreAlgorithm::cancel(){
		//There is really not a whole lot to do but if we have the ILP running we need to cancel it and retrieve incomplete solutions
		if (_solverActive){
			_solver.terminate();
		}
	}

	void expandSolutions(ClusterEditingInstance& cei,
			ClusterEditingSolutions& ces,
			vector<vector<vector<int> > >& partitions) {
		for (size_t k = 0; k < ces.getNumberOfSolutions(); k++) {
			vector<vector<int> > partition;
			for (size_t l = 0; l < ces.getNumberOfClusters(k); ++l) {
				vector<int> cluster;
				for (vector<int>::const_iterator cit = ces.getCluster(k, l).begin();
						cit != ces.getCluster(k, l).end(); ++cit) {
					cluster.insert(cluster.end(),
							cei.getCluster(cei.getOrig().nodeFromId(*cit))->begin(),
							cei.getCluster(cei.getOrig().nodeFromId(*cit))->end());
				}
				partition.push_back(cluster);
			}
			partitions.push_back(partition);
		}
	}

	void mergeSolutions(size_t i, size_t& k, vector<vector<int> >& partition,
			ClusterEditingSolutions& solutions,
			vector<vector<vector<vector<int> > > >& instances) {
		if (i == instances.size()) {
			solutions.getSolution(k).insert(solutions.getSolution(k).end(),
					partition.begin(), partition.end());
			k++;
			return;
		}
		for (size_t j = 0; j < instances[i].size(); j++) {
			vector<vector<int> > p;
			p.insert(p.end(), partition.begin(), partition.end());
			p.insert(p.end(), instances[i][j].begin(), instances[i][j].end());

			mergeSolutions(i + 1, k, p, solutions, instances);
		}
	}

}
