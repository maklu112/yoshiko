#include "GraphGenerator.h"

using namespace std;
using namespace lemon;

void GraphGenerator::newGraph(string &filename,int r){
  // Graph, random seed und Zählervariablen

  srand(time(NULL)+r*_size);
  int x,y;

  // Fülle leeren Graph mit schwachen kanten oder starken nicht-Kanten 1/3
  int weights[_size][_size];
  for(x=0;x<_size;x++){
    for(y=x+1;y<_size;y++){
      if(rand()% 3 == 1) weights[x][y] = 2;
      else weights[x][y] = -10;
    }
  }

  // Erstes Cluster zwischen n/k und n.
  int startpoint = (rand() % (_size - (_size / _var))) + (_size / _var);

  // Cluster sollen entweder starke Kanten oder schwache nicht-Kanten haben 1/3
  for(x=startpoint;x<_size;x++){
    for(y=x+1;y<_size;y++){
      if(rand() % 3 == 1){
        weights[x][y] = -2;
      }
      else weights[x][y] = 10;
    }
  }

  // Nach dem n/k -n sollen die restlichen cluster zufällig entstehen 0-i;
  int i = rand() % startpoint;

  while(i > 0){
    for(x=i;x<startpoint;x++){
      for(y=x+1;y<startpoint;y++){
        if(rand() % 3 == 1){
          weights[x][y] = -2;
        }
        else weights[x][y] = 10;
      }
    }
    startpoint = i;
    i = rand() % startpoint;
  }

  // wenn i=0 ist muss trotzdem noch 0-startpoint gefüllt werden
  for(x=i;x<startpoint;x++){
    for(y=x+1;y<startpoint;y++){
      if(rand() % 3 == 1){
        weights[x][y] = -2;
      }
      else weights[x][y] = 10;
    }
  }
/*
  for(x=0;x<_size;x++){
    for(y=x+1;y<_size;y++){
      printf("%d ",weights[x][y]);
    }
  } */
  //cout << "wtf" << endl;

  ofstream output;
  filename += ".jena";
  const char *chr = filename.c_str();

  output.open(chr);
  output << _size << "\n";
  for(x=0;x<_size;x++){
    output << x << "\n";
  }
  for(x=0;x<_size;x++){
    for(y=x+1;y<_size;y++){
      output << weights[x][y] << " ";
    }
    output << "\n";
  }

  output.close();
}
