#include <lemon/full_graph.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <string>

class GraphGenerator{
  public:

    GraphGenerator(int n, int k)
      :_size(n),
       _var(k)
    {};

    void newGraph(std::string &filename,int s);
    void newRandomGraph(std::string &filename, int s);

  private:
    int _size;
    int _var;
};
