# yoshiko

(Weighted) Cluster Editing


## Compilation instructions

### Dependencies

* LEMON 1.3
* ILOG CPLEX (>= 12.0)

### Compiling

Get yoshiko from github:

    git clone <HTTPS clone URL (see on the right side of this page)>


First, LEMON 1.3 needs to be installed:

Note: It is required to have LEMON as a shared library. If you obtained lemon prior to this installation and get errors during the compilation process consider recompiling it yourself as a shared library.

    wget http://lemon.cs.elte.hu/pub/sources/lemon-1.3.tar.gz
    tar xvzf lemon-1.3.tar.gz
    cd lemon-1.3
    cmake -DCMAKE_INSTALL_PREFIX=~/lemon -DBUILD_SHARED_LIBS=TRUE
    make install

You can remove the LEMON sources now, i.e., `rm -rf lemon-1.3`.
Make sure you move the created library file with ending (for instance .dll on windows, .dylib on mac, .so on linux) to the shared library folder of your system.

CPLEX is a commercial product owned by IBM. For academic purposes it can be obtained at no charge via IBM's Academic Initiative programme:

  http://www-03.ibm.com/ibm/university/academic/pub/page/membership

Next, yoshiko can be compiled:
Note: If you want to compile it as a stand-alone executable, use: cmake -DBUILD_AS_EXEC=ON

    mkdir build
    cd build
    cmake ..
    make

In case auto-detection of LEMON or CPLEX fails, do

    cmake \
    -DLIBLEMON_ROOT=~/lemon \
    -DCPLEX_INC_DIR=~/ILOG/cplex/include/ \
    -DCPLEX_LIB_DIR=~/ILOG/cplex/lib/x86-64_osx/static_pic \
    -DCONCERT_LIB_DIR=~/ILOG/concert/lib/x86-64_osx/static_pic \
    -DCONCERT_INC_DIR=~/ILOG/concert/include/ ..

## Running yoshiko

To run yoshiko on the test instance:

    ./yoshiko -f ../data/test/fusarium.txt -v 4

Get a list of options:

    ./yoshiko -h


## Documentation

### Data formats

The tool supports different line-oriented file formats. Each line is terminated by a line feed.

* The **Jena** file format is as follows (see also the example in `data/fusarium.txt`) and can be used to input weighted or unweighted instances of the cluster editing problem.

  First line: number of nodes, say n.

  Next n lines: user-defined node labels.

  Next n-1 lines: upper triangular matrix specifying the node weights. That is, the `j`th entry of line `i` specifies the weight `w_i,i+j`. As conventional, we interpret positive values as deletion costs of existing edges and negative values as insertion costs of non-existing edges.

* The **SIF** file format is the file format used by Cytoscape and can be used to input an unweighted version of cluster editing, which is just given as an unweighted graph. For documentation, see http://wiki.cytoscape.org/Cytoscape_User_Manual/Network_Formats. Here is a small example:

        ./yoshiko -f data/test/simple_sif_example.sif -F 2 -v 4
