#include <iostream>

#include <lemon/arg_parser.h>
#include <string>

#include "GraphGenerator.h"

using namespace lemon;
using namespace std;





int main(int argc, char * const argv[]){

  ArgParser ap(argc, argv);

  int n;
  int s = 0;
  int k = 3;
  bool random = false;
  string filename = "graph";

  ap.refOption("r", "completly random",random,false);
  ap.refOption("s", "randon seed",s,false);
  ap.refOption("n", "Number of Nodes", n, true);
  ap.refOption("k", "relative Number", k, false);
  ap.refOption("f", "output filename",filename,false);

  ap.parse();

  GraphGenerator g = GraphGenerator(n,k);

  if(random)
    g.newRandomGraph(filename,s);
  else{
    if(n>=k)
      g.newGraph(filename,s);
  }

  return 1;
}
