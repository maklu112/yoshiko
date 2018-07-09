#include <lemon/full_graph.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>

class GraphGenerator{
  public:

    GraphGenerator(int n, int k)
      :_size(n),
       _var(k)
    {};

    void newGraph();

  private:
    int _size;
    int _var;
};
