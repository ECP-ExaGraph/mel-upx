#include <sys/resource.h>
#include <sys/time.h>
#include <pthread.h>
#include <unistd.h>

#include <cassert>
#include <cstdlib>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <mpi.h>
#include <upcxx/upcxx.hpp>

#if defined(USE_UPX_RPC)
#include "maxematch_upx_rpc.hpp"
#else
#define USE_UPX_RMA 1
#include "maxematch_upx.hpp"
#endif

static std::string inputFileName;
static int me, nprocs;
static int ranksPerNode = 1;
static GraphElem nvRGG = 0;
static int generateGraph = 0;
static int randomEdgePercent = 0;
static bool randomNumberLCG = false;
static bool readBalanced = false;
static double threshold = 1.0E-6;
static bool isUnitEdgeWeight = true;

// parse command line parameters
static void parseCommandLine(const int argc, char * const argv[]);



int main(int argc, char *argv[])
{
    double t0, t1, td, td0, td1;

#if defined(USE_UPX_RPC) || defined(USE_UPX_RMA)
    upcxx::init(); // init UPC++
#endif

    // init MPI, if necessary
    int mpi_already_init;
    MPI_Initialized(&mpi_already_init);
    if (!mpi_already_init) MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    parseCommandLine(argc, argv);
 
    Graph* g = nullptr;

    if (me == 0) {
        std::cout << "Using " << nprocs << " ranks" << std::endl;
    }
    
    td0 = MPI_Wtime();

    if (generateGraph) 
    { 
        GenerateRGG gr(nvRGG);
        g = gr.generate(randomNumberLCG, isUnitEdgeWeight, randomEdgePercent);

        if (me == 0) 
        {
            std::cout << "Generated Random Geometric Graph with d: " << gr.get_d() << std::endl;
            const GraphElem nv = g->get_nv();
            const GraphElem ne = g->get_ne();
            std::cout << "Number of vertices: " << nv << std::endl;
            std::cout << "Number of edges: " << ne << std::endl;
            //std::cout << "Sparsity: "<< (double)((double)nv / (double)(nvRGG*nvRGG))*100.0 <<"%"<< std::endl;
            std::cout << "Average degree: " << (ne / nv) << std::endl;
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
    else 
    {
        BinaryEdgeList rm;
        if (readBalanced == true)
            g = rm.read_balanced(me, nprocs, ranksPerNode, inputFileName);
        else
            g = rm.read(me, nprocs, ranksPerNode, inputFileName);
    }
        
    g->print_dist_stats();
    assert(g != nullptr);

    MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG_PRINTF  
    assert(g);
#endif
    td1 = MPI_Wtime();
    td = td1 - td0;

    double tdt = 0.0;
    MPI_Reduce(&td, &tdt, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (me == 0)  
    {
        if (!generateGraph)
            std::cout << "Time to read input file and create distributed graph (in s): " 
                << tdt << std::endl;
        else
            std::cout << "Time to generate distributed graph of " 
                << nvRGG << " vertices (in s): " << tdt << std::endl;
    }

    // start graph matching

#if defined(USE_UPX_RPC)
    if (me == 0)
        std::cout << "MPI and UPCXX (RPC): ";
    MaxEdgeMatchUPXRPC mt(g);
#else
    if (me == 0)
        std::cout << "MPI and UPCXX (RMA): ";
    MaxEdgeMatchUPX mt(g);
    // create dist_object<...> per process
    upcxx::global_ptr<GraphElem> arr = upcxx::new_array<GraphElem>(
            mt.get_nelems());
    assert(arr);
    upcxx::dist_object<upcxx::global_ptr<GraphElem>> dobj(arr);
    mt.create_upx_ptr(dobj);
#endif

    upcxx::barrier();
    MPI_Barrier(MPI_COMM_WORLD);

    // invoke matching
    t0 = MPI_Wtime();
    std::vector<EdgeTuple> const& M = mt();
    
#if defined(USE_UPX_RPC) || defined(USE_UPX_RMA)
    upcxx::barrier();
#endif
    MPI_Barrier(MPI_COMM_WORLD);
    
    t1 = MPI_Wtime();
    double p_tot = t1 - t0, t_tot = 0.0;
    
    MPI_Reduce(&p_tot, &t_tot, 1, MPI_DOUBLE, 
            MPI_SUM, 0, MPI_COMM_WORLD);
    if (me == 0)
        std::cout << "Execution time (in s) for maximal edge matching: " 
            << (double)(t_tot/(double)nprocs) << std::endl;

#if defined(CHECK_RESULTS)    
    mt.check_results();
#endif
#if defined(PRINT_RESULTS)    
    mt.print_M();
#endif
    
#if defined(USE_UPX_RPC) || defined(USE_UPX_RMA)
    upcxx::barrier();
#endif
    MPI_Barrier(MPI_COMM_WORLD);

#if defined(USE_UPX_RMA)
    mt.destroy_upx_ptr();
    upcxx::barrier();
#endif
    mt.clear();

#if defined(USE_UPX_RPC) || defined(USE_UPX_RMA)
    upcxx::barrier();
#endif
    MPI_Barrier(MPI_COMM_WORLD);


    // Finalize MPI only if we initted it
    if (!mpi_already_init) MPI_Finalize();

#if defined(USE_UPX_RPC) || defined(USE_UPX_RMA)
    upcxx::finalize(); // finalize UPC++
#endif

    return 0;
}

void parseCommandLine(const int argc, char * const argv[])
{
  int ret;

  while ((ret = getopt(argc, argv, "f:br:t:n:wlp:")) != -1) {
    switch (ret) {
    case 'f':
      inputFileName.assign(optarg);
      break;
    case 'b':
      readBalanced = true;
      break;
    case 'r':
      ranksPerNode = atoi(optarg);
      break;
    case 't':
      threshold = atof(optarg);
      break;
    case 'n':
      nvRGG = atol(optarg);
      if (nvRGG > 0)
          generateGraph = true; 
      break;
    case 'w':
      isUnitEdgeWeight = false;
      break;
    case 'l':
      randomNumberLCG = true;
      break;
    case 'p':
      randomEdgePercent = atoi(optarg);
      break;
    default:
      assert(0 && "Should not reach here!!");
      break;
    }
  }

  if (me == 0 && (argc == 1)) {
      std::cerr << "Must specify some options." << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -99);
  }
  
  if (me == 0 && !generateGraph && inputFileName.empty()) {
      std::cerr << "Must specify a binary file name with -f or provide parameters for generating a graph." << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -99);
  }
   
  if (me == 0 && !generateGraph && randomNumberLCG) {
      std::cerr << "Must specify -g for graph generation using LCG." << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -99);
  } 
   
  if (me == 0 && !generateGraph && randomEdgePercent) {
      std::cerr << "Must specify -g for graph generation first to add random edges to it." << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -99);
  } 
  
  if (me == 0 && !generateGraph && !isUnitEdgeWeight) {
      std::cerr << "Must specify -g for graph generation first before setting edge weights." << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -99);
  }
  
  if (me == 0 && generateGraph && ((randomEdgePercent < 0) || (randomEdgePercent >= 100))) {
      std::cerr << "Invalid random edge percentage for generated graph!" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -99);
  }
  
  if (me == 0 && readBalanced && generateGraph) {
      std::cout << "Balanced graph distribution is only applicable to real-world graphs, and not applicable to synthetic graphs." << std::endl;
  }
} // parseCommandLine
