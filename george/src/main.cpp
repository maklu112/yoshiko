#include <iostream>

#include <lemon/arg_parser.h>
#include <string>

#include "GraphGenerator.h"

using namespace lemon;
using namespace std;





int main(int argc, char * const argv[]){

  ArgParser ap(argc, argv);

  int n;
  int r = 0;
  int k = 3;
  string filename = "graph";

  ap.refOption("r", "randon seed",r,false);
  ap.refOption("n", "Number of Nodes", n, true);
  ap.refOption("k", "relative Number", k, false);
  ap.refOption("f", "output filename",filename,false);

  ap.parse();

  GraphGenerator g = GraphGenerator(n,k);

  if(n>=k)
    g.newGraph(filename,r);

  return 1;
}
