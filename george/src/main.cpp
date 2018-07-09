#include <iostream>

#include <lemon/arg_parser.h>

#include "GraphGenerator.h"

using namespace lemon;





int main(int argc, char * const argv[]){

  ArgParser ap(argc, argv);

  int n;
  int k = 3;

  ap.refOption("n", "Number of Nodes", n, true);
  ap.refOption("k", "relative Number", k, false);

  ap.parse();

  GraphGenerator g = GraphGenerator(n,k);

  g.newGraph();

  return 1;
}
