#include "ClusterEditingSolutions.h"


using namespace std;
using namespace lemon;

namespace ysk {

	size_t ClusterEditingSolutions::getNumberOfClusters(size_t i) const {
	  return _solutions[i].size();
	}

	vector<int>& ClusterEditingSolutions::getCluster(size_t i, size_t k) {
	  return _solutions[i][k];
	}

	vector<vector<int> >& ClusterEditingSolutions::getSolution(size_t i) {
		//Out of bounds
		if(i >= _solutions.size()){
			//throw new Exception("Attempted to access solution with index: "+i);
			cerr << "Attempted to access solution with index: " <<i << endl;
			exit(-1);
		}
		return _solutions[i];
	}

	void ClusterEditingSolutions::resize(long numberOfSolutions) {
	  _solutions.resize(numberOfSolutions);
	}

	size_t ClusterEditingSolutions::getNumberOfSolutions() const {
	  return _solutions.size();
	}

	void ClusterEditingSolutions::setSolution(int k, const ClusterEditingInstance &inst) {
	  _solutions[k].resize(countNodes(inst.getOrig()));
	  int i = 0;
	  for(FullGraph::NodeIt u(inst.getOrig()); u != INVALID; ++u, i++) {
		_solutions[k][i].push_back(inst.getOrig().id(u));
	  }

	}

	void ClusterEditingSolutions::setSolution(int k, const WorkingCopyInstance& inst) {
	  _solutions[k].resize(countNodes(inst.getGraph()));

	  int i = 0;
	  for(WorkingCopyGraph::NodeIt u(inst.getGraph()); u != INVALID; ++u, i++) {
		for(vector<int>::iterator it = inst.getClusters()[u]->begin(); it != inst.getClusters()[u]->end(); ++it) {
		  _solutions[k][i].push_back(*it);
		}
	  }
	}

	void ClusterEditingSolutions::setSolution(int k, const IloNumArray &x_vals, const ClusterEditingInstance &i) {
	  // compute clusters
	  if (verbosity > 1)
		cout << "computing clusters for solution " << k << "... " << flush;

	  const FullGraph g = i.getOrig();

	  ListGraph c; // graph that will contain clusters as (fully) connected components

	  FullGraph::NodeMap<ListGraph::Node> A(g);

	  for (FullGraph::NodeIt v(g); v != INVALID; ++v)
		A[v] = c.addNode();

	  for (FullGraph::EdgeIt e(g); e != INVALID; ++e)
		if (x_vals[g.id(e)] > 1 - eps)
		  c.addEdge(A[g.u(e)], A[g.v(e)]);

	  ListGraph::NodeMap<int> comp_num(c);
	  _solutions[k].resize(connectedComponents(c, comp_num));

	  for (FullGraph::NodeIt v(g); v != INVALID; ++v)
		_solutions[k][comp_num[A[v]]].push_back(g.id(v));

	  if (verbosity > 1)
		cout << "done." << endl;
	}

	// Same as above for Symphony
	void ClusterEditingSolutions::setSolution(int* indizes,int n,const double results[], const ClusterEditingInstance &i) {
	  // compute clusters
	  if (verbosity > 1)
		cout << "computing clusters for solution " << 1 << "... " << flush;

	  const FullGraph g = i.getOrig();
		int x,y;

	  ListGraph c; // graph that will contain clusters as (fully) connected components

	  FullGraph::NodeMap<ListGraph::Node> A(g);

	  for (FullGraph::NodeIt v(g); v != INVALID; ++v)
		A[v] = c.addNode();

		x = 0;
	  for (FullGraph::NodeIt i(g); i != INVALID; ++i){
			FullGraph::NodeIt j(g); j = i;
			for (++j,y=x+1; j != INVALID, y<n; ++j, y++){
				cout << x*n+y << endl;
				if (results[indizes[x*n+y]] > 1 - eps)
		  		c.addEdge(A[i], A[j]);
			}
			x++;
		}

		cout << "echt jetzt?" << endl;

	  ListGraph::NodeMap<int> comp_num(c);
	  _solutions[0].resize(connectedComponents(c, comp_num));

	  for (FullGraph::NodeIt v(g); v != INVALID; ++v)
		_solutions[0][comp_num[A[v]]].push_back(g.id(v));

	  if (verbosity > 1)
		cout << "done." << endl;
	}

	SolutionFlags ClusterEditingSolutions::getFlags(){
		return _flags;
	}

	void ClusterEditingSolutions::setFlags(SolutionFlags f){
		_flags = f;
	}

	void ClusterEditingSolutions::printSolution(size_t index){
		int idx = 0;
		for(auto &entry : _solutions[index]){
			cout << "Cluster " << idx << ": ";
			for (auto i = entry.begin(); i != entry.end(); ++i)
			    cout << *i<<" ";
			cout <<endl;
			idx ++;
		}
		cout << endl;
	}


} // namespace ysk
