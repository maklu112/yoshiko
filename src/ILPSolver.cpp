#include "ILPSolver.h"

using namespace std;
using namespace yskLib;
using namespace lemon;

namespace ysk {

/*************************************************************************
 * This callback procedure is called by Cplex during the computation
 * of the ILP.
 *************************************************************************/
ILOMIPINFOCALLBACK1(gap_callback,yskLib::CplexInformer*, informer) {
	double gap = getMIPRelativeGap();
	if (verbosity >= 3) cout << "Informing listener about new gap value: " << gap << endl;
    informer->updateStatus(SOLVING_ILP,gap);
    return;
}

/*************************************************************************
 * This callback procedure is called by Cplex during the computation
 * of the ILP.
 *************************************************************************/
ILOLAZYCONSTRAINTCALLBACK4(triangle_callback, const ClusterEditingInstance&, inst,  const IloBoolVarArray&, x, const IloBoolVarArray&, y,  const int, clusterCount) {

    IloEnv env = getEnv();
    //const double epsilon = IloCplex::EpInt;  gunnar: EpInt ist falsch, aber vielleicht was anderes
    //cout << "in callback "  << endl;
    try {

        int no_added_triplet_cuts = 0; // no cuts added yet in this iteration
        //bool fractional = false; //TODO: @Gunnar, is this still needed?

        //We fetch the variables for further reference
        IloNumArray x_vals(env);
        IloNumArray y_vals(env);
        //cout << "Fetching x and y values" << endl;
        getValues(x_vals, x);
        //cout << "[X DONE]" << endl;
        getValues(y_vals, y);
        //cout << "[Y DONE]" << endl;

        const FullGraph g = inst.getOrig();
        const int n = g.nodeNum();

        // round LP solution a bit:

        for (FullGraph::EdgeIt e(g); e != INVALID; ++e) {
            if (x_vals[g.id(e)] < eps) x_vals[g.id(e)] = 0.0;
            if (x_vals[g.id(e)] > 1.0 - eps) x_vals[g.id(e)] = 1.0;
            //if (x_vals[g.id(e)] > eps && x_vals[g.id(e)] < 1.0 - eps) fractional = true;
        }


        //Round k-cluster stuff as well ???
        for (int i = 0; i < clusterCount; i++){
            //TODO: Maybe?
        }

        unsigned long long no_checks = (n+clusterCount) * (n+clusterCount-1) * (n+clusterCount-2) / 6, no_checked = 0; //TODO: Why long long?

        //	   cout << " versus " << no_checks << endl;

        for (FullGraph::NodeIt i(g); i != INVALID; ++i) {
            if (time_limit != -1) {
                if (clk.systemTime() > time_limit) //TOOD: Coherent use of time_limit in ALL of yoshiko
                    throw("Time limit exceeded in yoshiko triangle callback");

                cout << clk << endl;
            }

            FullGraph::NodeIt j(g); j = i;
            for (++j; j != INVALID; ++j) {
                if (verbosity > 1) // && ((no_checked++ % 10000) == 0))
                    cout << "\rchecking triangle inequalities... (" << (no_checked)/double(no_checks) * 100
                    << " %) added " << no_added_triplet_cuts << "              " << flush;

                FullGraph::NodeIt k(g); k = j;
                //      for (++k; k != INVALID && no_added < 100000; ++k) { // if (k != i && k != j) {
                for (++k; k != INVALID; ++k) { // if (k != i && k != j) {
                    ++no_checked;
                    if (x_vals[g.id(g.edge(i, j))] + x_vals[g.id(g.edge(j, k))] - x_vals[g.id(g.edge(i, k))] > 1 + 3 * eps) {
                        add(x[g.id(g.edge(i, j))] + x[g.id(g.edge(j, k))] - x[g.id(g.edge(i, k))] <= 1);
                        ++no_added_triplet_cuts;
                    }

                    if (x_vals[g.id(g.edge(i, j))] - x_vals[g.id(g.edge(j, k))] + x_vals[g.id(g.edge(i, k))] > 1 + 3 * eps) {
                        add(x[g.id(g.edge(i, j))] - x[g.id(g.edge(j, k))] + x[g.id(g.edge(i, k))] <= 1);
                        ++no_added_triplet_cuts;
                    }

                    if (-x_vals[g.id(g.edge(i, j))] + x_vals[g.id(g.edge(j, k))] + x_vals[g.id(g.edge(i, k))] > 1 + 3 * eps) {
                        add(-x[g.id(g.edge(i, j))] + x[g.id(g.edge(j, k))] + x[g.id(g.edge(i, k))] <= 1);
                        ++no_added_triplet_cuts;
                    }
                }

                //Additional conditions for virtual nodes / "cluster centers"
                if (clusterCount != 0){ //This condition is probably redundant and only for easier readability
                    for (int k = 0; k <  clusterCount; ++k) {
                        //We increment the number of checked inequalities (to be precise we check 3 on each iteration)
                        ++no_checked;

                        if (y_vals[k*clusterCount + g.id(j)] + x_vals[g.id(g.edge(j, i))] - y_vals[k*clusterCount + g.id(i)] > 1 + 3 * eps) {
                            add( y[k*clusterCount + g.id(j)] + x[g.id(g.edge(j, i))] - y[k*clusterCount + g.id(i)] <= 1);
                            ++no_added_triplet_cuts;
                        }

                        if (y_vals[k*clusterCount + g.id(j)] - x_vals[g.id(g.edge(j, i))] + y_vals[k*clusterCount + g.id(i)] > 1 + 3 * eps) {
                            add( y[k*clusterCount + g.id(j)] - x[g.id(g.edge(j, i))] + y[k*clusterCount + g.id(i)] <= 1);
                            ++no_added_triplet_cuts;
                        }

                        if (y_vals[k*clusterCount + g.id(j)] + x_vals[g.id(g.edge(j, i))] + y_vals[k*clusterCount + g.id(i)] > 1 + 3 * eps) {
                            add(-y[k*clusterCount + g.id(j)] + x[g.id(g.edge(j, i))] + y[k*clusterCount + g.id(i)] <= 1);
                            ++no_added_triplet_cuts;
                        }
                    }
                }

            }
        }

        if (verbosity > 1)
            cout << "done." << endl;

        return;
    }
    catch (IloException& e) {
        env.out() << e.getMessage() << endl;
        throw -1;
    }
    catch (Exception& e) {
        env.out() << e.what() << endl;
        throw -1;
    }
    catch (...) {
        env.out() << "something strange in separation callback" << endl;
        throw -1;
    }
}

/*************************************************************************
 * This callback procedure is called by Cplex during the computation
 * of the ILP.
 *************************************************************************/
ILOUSERCUTCALLBACK2(partition_cut_callback, const ClusterEditingInstance&, inst,  const IloBoolVarArray&, x) {

    IloEnv env = getEnv();
    //const double epsilon = IloCplex::EpInt;  gunnar: EpInt ist falsch, aber vielleicht was anderes
    cout << "in callback "  << endl;
    try {

        int no_added = 0; // no cuts added yet in this iteration
        //bool fractional = false;

        IloNumArray x_vals(env);
        getValues(x_vals, x);

        const FullGraph g = inst.getOrig();

        // round LP solution a bit:

        for (FullGraph::EdgeIt e(g); e != INVALID; ++e) {
            if (x_vals[g.id(e)] < eps) x_vals[g.id(e)] = 0.0;
            if (x_vals[g.id(e)] > 1.0 - eps) x_vals[g.id(e)] = 1.0;
            //if (x_vals[g.id(e)] > eps && x_vals[g.id(e)] < 1.0 - eps) fractional = true;
        }

        // we do have fractional values:
        //cout << "fractional solution in callback..." << endl;

        // heuristic 2
        int no_added_partition_cuts = 0; int cnt_i = 0;
        for (FullGraph::NodeIt i(g); i != INVALID; ++i) {
            list<FullGraph::Node> W;
            for (FullGraph::NodeIt j(g); j != INVALID; ++j) if (i < j)
                if (x_vals[g.id(g.edge(i, j))] > 0.0)
                    W.push_back(j);

            // permute??? warum???

            for (list<FullGraph::Node>::const_iterator cit_j = W.begin(); cit_j != W.end(); ++cit_j) {
                FullGraph::Node j = *cit_j;
                IloExpr cut(env);
                double cut_value = 0.0;

                list<FullGraph::Node> T;
                T.push_back(j);
                int idx_ij = g.id(g.edge(i, j));
                cut += x[idx_ij];
                cut_value += x_vals[idx_ij];

                for (list<FullGraph::Node>::const_iterator cit_k = W.begin(); cit_k != W.end(); ++cit_k) if (*cit_k != j) {
                    FullGraph::Node k = *cit_k;
                    double take_k = 0.0;

                    for (list<FullGraph::Node>::const_iterator cit_l = T.begin(); cit_l != T.end(); ++cit_l)
                        take_k += x_vals[g.id(g.edge(k, *cit_l))];

                    //obacht: check .05
                    if (x_vals[g.id(g.edge(i, k))] - take_k > .05) {
                        T.push_back(k);
                        cut += x[g.id(g.edge(i, k))];
                        cut_value += x_vals[g.id(g.edge(i, k))];
                    }
                }

                for (list<FullGraph::Node>::const_iterator cit_t1 = T.begin(); cit_t1 != T.end(); ++cit_t1)
                    for (list<FullGraph::Node>::const_iterator cit_t2 = T.begin(); cit_t2 != T.end(); ++cit_t2) if (*cit_t1 < *cit_t2) {
                        cut -= x[g.id(g.edge(*cit_t1, *cit_t2))];
                        cut_value -= x_vals[g.id(g.edge(*cit_t1, *cit_t2))];
                    }

                // check whether ineq is violated:
                if (cut_value > 1.0 + eps && cut_value/((double) T.size()) > .3) {
                    //cout <<  cut_value/((double) T.size()) << endl;
                    //cout << cut << " = " << cut_value << endl;
                    //cout << "+";
                    no_added++; ++no_added_partition_cuts;
                    add(cut <= 1);
                    cut.end();
                }
            }
            if (verbosity > 1) {
                cout << "\rchecked partition inequalities... (" << ++cnt_i << "/" << g.nodeNum() << ")   cuts added: " << no_added_partition_cuts << flush;
            }
        }
        if (verbosity > 1)
            cout << endl;

        return;
    }
    catch (IloException& e) {
        env.out() << e << endl;
        throw -1;
    }
    catch (Exception& e) {
        env.out() << e.what() << endl;
        throw -1;
    }
    catch (...) {
        env.out() << "something strange in separation callback" << endl;
        throw -1;
    }
}

void ILPSolver::registerInformer(yskLib::CplexInformer* informer){
	if (verbosity > 2) cout << "Registered CplexInformer @ ILPSolver" << endl;
	_informer = informer;
}


/**
 * Gracefully terminates Cplex by calling the aborter
 */
void ILPSolver::terminate(){
	if (verbosity > 2){
		cout << "Terminating CPLEX Environment" << endl;
	}
	if (_cplexInitialized){
		_aborter.abort();
	}
}

/**
 * Attempts to solve the given ClusterEditingInstance using Cplex and parses the solution into the provided CES
 * @param inst The problem instance that is to be solved
 * @param s The solution in which the results are stored (Should usually be a new/empty instance)
 * @param flags A flags object, which is modified to reflect events occuring during the solving process
 */
long ILPSolver::solveCPLEX(const ClusterEditingInstance& inst, ClusterEditingSolutions& s, SolutionFlags& flags) {

    if(inst.isDirty()) {
	//TODO: Better error handling, when would this happen?
        cerr << "Fatal error: ClusterEditingInstance is dirty."<<endl;
        exit(-1);
    }

    //Fetch the full graph
    const FullGraph g = inst.getOrig();

    //skip if size is 1 because cplex throws an error in this case.
    //@gunnar, vielleicht kann man das auch schöner umgehen... //TODO: Still relevant?
    if(g.nodeNum() <= 1) {
        if (verbosity > 1) {
            cout << "size of instance is 1. No further computation is necessary."<<endl;
        }
        s.resize(1);
        s.setSolution(0, inst);
    	flags.solvedInstances +=1;
        return 1;
    }

    // first build ILOG model...
    IloEnv cplexEnv; // get ILOG environment
    IloModel M(cplexEnv); // get model
    IloCplex cplex(cplexEnv); // get cplex

    _aborter = IloCplex::Aborter(cplexEnv); //generate Aborter object and link to cplex
    cplex.use(_aborter);
    _cplexInitialized = true; //flag for external objects to check if a cplex instance is initialized (and can thus be stopped)

    // shut up cplex (or don't)
    if (verbosity < 3) {
    	cplex.setOut(cplexEnv.getNullStream());
    	cplex.setWarning(cplexEnv.getNullStream());
    	cplex.setError(cplexEnv.getNullStream());
    }
    else{
    	cplex.setParam(IloCplex::MIPDisplay,5);
    	cplex.setOut(cout);
    }

    //If applicable, set thread number
    if (no_threads != -1)
    	cplex.setParam(IloCplex::Threads, no_threads);
    if (verbosity > 3)
        cout << "Running CPLEX with " << no_threads << " threads ... " << endl;

    // set CPU time limit
    cplex.setParam(IloCplex::ClockType, 1);
    if (time_limit != -1)
    	cplex.setParam(IloCplex::TiLim, time_limit);


    //Limit Memory usage as this runs in a GUI and we should still be able to use the GUI
    //cplex.setParam(IloCplex::NodeFileInd,2);


    //TODO: Potential for research, observe relations between cuts / running time / instance properties

    //Set all generic cuts off
//     cplex.setParam(IloCplex::Cliques, -1);
//     cplex.setParam(IloCplex::Covers, -1);
//     cplex.setParam(IloCplex::DisjCuts, -1);
//     cplex.setParam(IloCplex::FlowCovers, -1);
//     cplex.setParam(IloCplex::FracCuts, -1);
//     cplex.setParam(IloCplex::GUBCovers, -1);
//     cplex.setParam(IloCplex::Cliques, -1);
//     cplex.setParam(IloCplex::ImplBd, -1);
//     cplex.setParam(IloCplex::MIRCuts, -1);
//     cplex.setParam(IloCplex::FlowPaths, -1);

    int n = g.nodeNum();


    //Some output because why not
    if (verbosity > 1) {
        cout << n << " nodes" << endl;
    }

        if (_useKCluster) cout << "k fixed" << endl;

    //Generate variables:

    IloBoolVarArray x(cplexEnv, g.edgeNum());
    for (FullGraph::EdgeIt e(g); e != INVALID; ++e) {
        std::stringstream var_name;
        var_name << "x_" << g.id(g.u(e)) << "_" << g.id(g.v(e));
        x[g.id(e)] = IloBoolVar(cplexEnv, var_name.str().c_str());
        M.add(x[g.id(e)]);
        //cout << g.id(g.source(a)) << " --> " << g.id(g.target(a)) << "\t" << g.id(a) << " " << g.edge_no(g.id(g.source(a)), g.id(g.target(a))) << endl;
    }

    IloBoolVarArray y(cplexEnv, _clusterCount * g.nodeNum());
    if (_useKCluster && _sep_triangles) M.add(y); //Explicit addition to model (as not used previously)
    if (_useKCluster) { // add k auxiliary vertices and edges to all the others
        for (int i = 0; i <  _clusterCount; ++i) {
            for (FullGraph::NodeIt v(g); v != INVALID; ++v) {
                std::stringstream var_name;
                var_name << "y_" << i  << "_" << g.id(v);
                y[i*g.nodeNum() + g.id(v)] = IloBoolVar(cplexEnv, var_name.str().c_str());
                M.add(y[i*g.nodeNum() + g.id(v)]);
            }
        }
    }

    //Build objective function
    IloExpr obj_expr(cplexEnv);
    for (FullGraph::EdgeIt e(g); e != INVALID; ++e){
        if (!inst.isForbidden(e) && !inst.isPermanent(e)) {
            obj_expr -= inst.getWeight(e) * x[g.id(e)]; //-= is overload for adding a linear expression
            if (inst.getWeight(e) > 0)
                obj_expr += inst.getWeight(e);
        }
    }

    M.add(IloObjective(cplexEnv, obj_expr, IloObjective::Minimize));

    //Add inequalities for forbidden/permanent edges:
    for (FullGraph::EdgeIt e(g); e != INVALID; ++e) {
        if (inst.isForbidden(e))
            M.add(x[g.id(e)] == 0);

        if (inst.isPermanent(e))
            M.add(x[g.id(e)] == 1);
    }

    //Set edges in between centers to forbidden (Only in k-cluster mode)
    //if (_useKCluster){
        //for (int i = 0; i <  _clusterCount; ++i) {
            //for (int j = 0; j <  _clusterCount; ++j) {
                //TODO: As those edges don't even exist as variables, is this needed at all?
            //}
        //}
    //}

    //CALLBACKS

    if (_sep_triangles) cplex.use(triangle_callback(cplexEnv, inst, x,y,_useKCluster ? _clusterCount : 0)); // use triangle callback --> lazy constraints!
    if (_sep_partition_cuts) cplex.use(partition_cut_callback(cplexEnv, inst, x)); // use partition callback --> user cuts!
    //if we're using an informer to report to any outside source we will register it here
    if (_informer != nullptr){
    	cplex.use(gap_callback(cplexEnv,_informer));
    	//TODO: Maybe use even more callbacks to provide more info
    }

    if (!_sep_triangles) {
        if (verbosity > 1) {
            cout << "adding triangle inequalities... " << flush;
        }

        for (FullGraph::NodeIt i(g); i != INVALID; ++i) {
            FullGraph::NodeIt j(g); j = i;
            for (++j; j != INVALID; ++j) {
                FullGraph::NodeIt k(g); k = j;
                for (++k; k != INVALID; ++k) {
                    M.add( x[g.id(g.edge(i, j))] + x[g.id(g.edge(j, k))] - x[g.id(g.edge(i, k))] <= 1);
                    M.add( x[g.id(g.edge(i, j))] - x[g.id(g.edge(j, k))] + x[g.id(g.edge(i, k))] <= 1);
                    M.add(-x[g.id(g.edge(i, j))] + x[g.id(g.edge(j, k))] + x[g.id(g.edge(i, k))] <= 1);
                }
            }
        }

        //Additional Triangle Equalities for K-Cluster Problem

        if (_useKCluster){
            for (int i = 0; i <  _clusterCount; ++i) {
                for (FullGraph::NodeIt j(g); j != INVALID; ++j) {
                    FullGraph::NodeIt k(g); k = j;
                    for (++k; k != INVALID; ++k) {
                        M.add( y[i*g.nodeNum() + g.id(j)] + x[g.id(g.edge(j, k))] - y[i*g.nodeNum() + g.id(k)] <= 1);
                        M.add( y[i*g.nodeNum() + g.id(j)] - x[g.id(g.edge(j, k))] + y[i*g.nodeNum() + g.id(k)] <= 1);
                        M.add(-y[i*g.nodeNum() + g.id(j)] + x[g.id(g.edge(j, k))] + y[i*g.nodeNum() + g.id(k)] <= 1);
                    }
                }
            }
        }
    }

    //Additional Constraints for K-Cluster problem
    if(_useKCluster){
        //First, we want to assure, that every node is connected to EXACTLY one cluster center,
        //In other words we want to make sure that the number of edges from a given node to cluster centers is = 1
        for (FullGraph::NodeIt i(g); i != INVALID; ++i) {
            IloExpr node_association(cplexEnv);
            for (int j = 0; j <  _clusterCount; ++j) {
                node_association += y[j*g.nodeNum() + g.id(i)];
            }
            M.add(node_association == 1);
        }

        //Then we want to also guarantee that no cluster centers are "empty" meaning they should have at least one node associated
        for (int j = 0; j <  _clusterCount; ++j) {
            IloExpr cluster_association(cplexEnv);
            for (FullGraph::NodeIt i(g); i != INVALID; ++i) {
                cluster_association += y[j*g.nodeNum()+ g.id(i)];
            }
            M.add(cluster_association >= 1);
        }

    }

    if (verbosity > 1)
        cout << "done." << endl;

    //cout << "extracting (and exporting) ILP... " << flush; (TODO: Sinn? is anything happening here that I am missing?)

    if (verbosity > 1) {
        cout << "extracting ILP... " << flush;
    }

    cplex.extract(M);

    //if (verbose) cplex.exportModel("yoshiko.lp");

    if (verbosity > 1) {
        cout << "done." << endl;
    }

    if (verbosity > 1) {
        cout << "solving ILP (computing one optimal solution)... " << flush;
    }

    //Inform informer that the ILP solving process starts
    if (_informer != nullptr){
    	_informer->updateStatus(YoshikoState::SOLVING_ILP);
    }

    //Start Solver, returns true if the solution is feasible (not optimal!)
    flags.ilpGenerated = true; //Mark solution as being generated by ILP
    bool feasible = cplex.solve();

    //Time-out handling
    if (cplex.getCplexStatus() == IloCplex::AbortTimeLim){
    	bool resume = false;
    	//If we use the informer we have the possibility of asking the listener to resume the solver
    	if (_informer != nullptr){
    		if (verbosity > 1) cout << "Requesting further instructions ..." << endl;
    		if (_informer->continueOnTimeout()){
    			if (verbosity > 2) cout << "Continuing at users request!" << endl;
    			resume = true;
    			cplex.setParam(IloCplex::TiLim, IloCplex::IntParam_MAX);
    			feasible = cplex.solve();
    		}
    	}
    	if (!resume){
        	cout << "[CPLEX TIME OUT]" << endl;
        	flags.timedOut =  true;
    	}

    }

    //For some reason the solver terminated without providing an optimal solution
    if (!feasible) {
        cout << endl << endl << endl;
        cout << cplex.getStatus() << endl;
        cout << endl << endl << endl;
        cerr << "yoshiko: Optimization problems. CPLEX status code " << cplex.getStatus() << endl;
    }
    else if (cplex.getCplexStatus() != IloCplex::Optimal){
    	flags.optimal = false;
    }
    else{
    	flags.solvedInstances +=1;
    }

    //Assign cost to the solution
    double z = cplex.getObjValue();
    flags.totalCost = z;
    if(verbosity >1){
        cout << "Total Cost: " << z << endl;
    }

    //Add Gap for instance
    flags.lastGap = cplex.getMIPRelativeGap();

    //Some output
    if (verbosity > 1) {
        cout << "CPLEX status code " << cplex.getStatus() << endl;
        cout << "upper bound: " << z << endl;
        cout << "lower bound: " << cplex.getBestObjValue() << endl;
        cout << "gap (%):     " << cplex.getMIPRelativeGap()*100 << endl;
    }


    //Fetch additional solution
    long numsol = 1;
    if (_num_opt_sol > 1) {
    	cplex.setParam(IloCplex::SolnPoolGap, 0.0); //
    	cplex.setParam(IloCplex::SolnPoolIntensity, 4); // "all" (sub)optimal solutions
    	cplex.setParam(IloCplex::PopulateLim, _num_opt_sol);


        if (verbosity > 1)
            cout << "solving ILP (finding more optimal solutions)..." << flush;
        cplex.populate();

        if (verbosity > 1)
            cout << "done." << endl;

        numsol = cplex.getSolnPoolNsolns();
        if (verbosity > 1)
            cout  << numsol << " optimal solutions." << endl;
    }

    //Add solutions to return instance
    s.resize(numsol);
    for (int k = 0; k < numsol; ++k) {
        IloNumArray x_vals(cplexEnv);
        cplex.getValues(x_vals, x, k);
        s.setSolution(k, x_vals, inst);
    }


    //CLEANUP! Important to keep memory consumption reasonable
    cplexEnv.end();

    return numsol;
}

long ILPSolver::solveCOIN(const ClusterEditingInstance& inst, ClusterEditingSolutions& s, SolutionFlags& flags) {

//TODO: Error handling and verbosity stuff

	const FullGraph g = inst.getOrig();

// Ich kopiert hier einfach mal die Beiden Fehlerabfragen aus solveCPLEX
	if(inst.isDirty()) {
//TODO: Better error handling, when would this happen?
			cerr << "Fatal error: ClusterEditingInstance is dirty."<<endl;
			exit(-1);
	}

	if(g.nodeNum() <= 1) {
			if (verbosity > 1) {
					cout << "size of instance is 1. No further computation is necessary."<<endl;
			}
			s.resize(1);
			s.setSolution(0, inst);
			flags.solvedInstances +=1;
			return 1;
	}

	int n = g.nodeNum();
	int indizes[n][n],indizesY[_clusterCount][n];
  int m,x,y,z=0; //Zählervariablen

  OsiSymSolverInterface si;
  CoinBuild bo = CoinBuild(0); // build object to build the problem model

	cout << "Initialisiere " << n << " Spalten" << endl;
  // set bounds so edges are between 0 and 1, later we set them as integers and they become 0 or 1.
  // Also initialize the Columns
  for(m=0;m<(n*(n-1))/2;m++){
    si.addCol(0,{},{},0,1,1);
  }

	for(x=0;x<n;x++){
		for(y=x+1;y<n;y++){
			indizes[x][y] = z;
			z++;
		}
	}

	if(_useKCluster){
		for(m=0;m<_clusterCount;m++){
			si.addCol(0,{},{},0,1,0);
		}
		for(m=0;m<_clusterCount;m++){
			for(y=0;y<n;y++){
				indizesY[m][y] = z;
				z++;
			}
		}
	}

	std::list<int> weights;
	double ele[][3] = {{1,1,-1},{1,-1,1},{-1,1,1}};

	cout << "Erstelle Zeilen" << endl;
	x = 0;
	for (FullGraph::NodeIt i(g);i != INVALID ;++i) {
			FullGraph::NodeIt j(g); j = i;
			for (++j,y=x+1; j != INVALID, y<n; ++j, y++) {
					FullGraph::NodeIt k(g); k = j;
					if(inst.isForbidden(g.edge(i,j))){
						si.setColLower(indizes[x][y],0);
						si.setColUpper(indizes[x][y],0);
					}
					if(inst.isPermanent(g.edge(i,j))){
						si.setColLower(indizes[x][y],1);
						si.setColUpper(indizes[x][y],1);
					}
					weights.push_back(inst.getWeight(g.edge(i,j)));
					for (++k, z=y+1; k != INVALID, z<n; ++k,z++) {
						int col[] = {indizes[x][y],indizes[y][z],indizes[x][z]};
						if(!inst.isForbidden(g.edge(i,j))&&!inst.isPermanent(g.edge(i,j))&&!inst.isForbidden(g.edge(j,k))&&!inst.isPermanent(g.edge(j,k))&&!inst.isForbidden(g.edge(i,k))&&!inst.isPermanent(g.edge(i,k)))
							bo.addRow(3,col,ele[0],-1,1);
							bo.addRow(3,col,ele[1],-1,1);
							bo.addRow(3,col,ele[2],-1,1);
					}
					// Inequalities for K-Clusterin
					if(_useKCluster){
						for(m=0;m<_clusterCount;m++){
							int colY[] = {indizes[x][y],indizesY[m][x],indizesY[m][y]};
							bo.addRow(3,colY,ele[0],-1,1);
							bo.addRow(3,colY,ele[1],-1,1);
							bo.addRow(3,colY,ele[2],-1,1);
						}
					}
			}
			// Jeder Knoten hat genau einen Clusterknoten
			if(_useKCluster){
				int noEmptyNode[_clusterCount];
				double eleY[_useKCluster];
				for(m=0;m<_clusterCount;m++){
					noEmptyNode[m] = indizesY[m][x];
					eleY[m] = 1;
				}
				bo.addRow(_clusterCount,noEmptyNode,eleY,1,1);
			}

			// Clustersize, each node needs a min or max number of edges;
			if(_minCluster > 1 || _maxCluster < INT_MAX){
				int neighbors[n-1];
				double newEle[n-1];
				for(m=0;m<x;m++){
					neighbors[m] = indizes[m][x];
					newEle[m] = 1;
				}
				for(m=x+1;m<n;m++){
					neighbors[m-1] = indizes[x][m];
					newEle[m-1] = 1;
				}
				bo.addRow(n-1,neighbors,newEle,_minCluster,_maxCluster);
			}
			x++;
	}

	// Every Cluster needs at least one Node
	if(_useKCluster){
		int noEmptyCluster[n];
		double eleY[n];
		for(m=0;m<_clusterCount;m++){
			for(x=0;x<n;x++){
				noEmptyCluster[x] = indizesY[m][x];
				eleY[x] = 1;
			}
			bo.addRow(n,noEmptyCluster,eleY,1,n);
		}
	}

  printf("Variablen: %d  Inequalities: %d \n",si.getNumCols(),bo.numberRows());

	try {
		si.addRows(bo);
		} catch (CoinError e) {
		cout << "Fehler bei Modellerstellung" << endl;
		cout << "COIN-Error: " << e.message() << endl;
	}
  double sumWeight = 0;
	std::list<int>::iterator ptr;
  for(m=0,ptr = weights.begin();m<(n*n-1)/2,ptr != weights.end();m++,ptr++){
    si.setInteger(m);
		cout << "set " << m << " as integer" << endl;
    si.setObjCoeff(m,-1*(*ptr));
    if(*ptr>0){
      //printf("%f %f \n",&ptr, sumWeight);
      sumWeight += *ptr;
    }
  }

	if(_useKCluster){
		for(m=(n*(n-1))/2;m<((n*(n-1))/2)+_clusterCount*n;m++){
			si.setInteger(m);
		}
	}

  si.addCol(0,{},{},1,1,sumWeight);

  si.findIntegers(false);
	//si.writeMps("tmpoutput");
	try {
		si.initialSolve();
		} catch (CoinError e) {
		cout << "Fehler beim Solving" << endl;
		cout << "COIN-Error: " << e.message() << endl;
	}
	const double *results = si.getColSolution();
	flags.totalCost = si.getObjValue();

	if(si.isProvenOptimal()){
		return 1;
	}
	else{
		return 0;
	}

}

} // namespace ysk
