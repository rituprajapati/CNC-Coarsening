#include <cnc/cnc.h>
#include <utility>
#include <iostream>

#include <cmath>
#include <stdlib.h>

#include "Vector.h"
#include "Matrix.h"
#include "Node.h"
#include "Twoscalecoeffs.h"
#include "Quadrature.h"
#include "Test_functions.h"
#include <cstdio>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <stack>
// http://stackoverflow.com/questions/22387586/measuring-execution-time-of-a-function-in-c
#include <chrono>
using namespace std;
// http://stackoverflow.com/questions/22387586/measuring-execution-time-of-a-function-in-c
using namespace std::chrono;

unordered_set<string> st;
unordered_map<int, int> height_map;
int thread = 0;

/**********************************************************************/
/*  Forward declaration
/**********************************************************************/

struct CnCContext;
Vector sub(const Vector &v1, const Vector &v2);
double sub_scale_factor(double n);
void init_twoscale(int );
int global_co_height;


/**********************************************************************/
/* This struct is used by different steps to put the 
/* output in item collection
/**********************************************************************/

template<typename K, typename V>
struct OutputTerminal {

  CnC::item_collection<K, V> *item_collection;
  std::vector<CnC::tag_collection<K> *> next_op_tags;

  /*----------------------------------------------------------------*/
  /* Output terminal constructor
  /*----------------------------------------------------------------*/
  OutputTerminal(CnC::item_collection<K, V> *item_collection,
                 std::vector<CnC::tag_collection<K> *> next_op_tags)
                : 
                item_collection(item_collection), next_op_tags(next_op_tags) {}


  /*----------------------------------------------------------------*/
  /* This method is invoked to put item in the output_item_collection 
  /* as well as putting the tag to the next tag_collections
  /*----------------------------------------------------------------*/
  void put(K key, V value) const {
    item_collection->put(key, value);

    if( key.first == 0 )
      for (CnC::tag_collection<K> *next_op_tag: next_op_tags) {
        next_op_tag->put(key);
      }
  }

  /*----------------------------------------------------------------*/
  /* This method is invoked to put tag to the next tag_collections
  /* --> mainly used for recursion purposes
  /*----------------------------------------------------------------*/
  void put(K key) const {
    for (CnC::tag_collection<K> *next_op_tag: next_op_tags) {
      next_op_tag->put(key);
    }
  }
};



/**********************************************************************/
/* Struct  - step_Base
/* Base struct to be used by different operators   
/**********************************************************************/
struct step_Base {

public:
  std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *> input_terminals;
  std::vector<OutputTerminal<std::pair< int, pair<int, int>>, Node>> output_terminals;

  step_Base( 
            std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *> input_terminals,
            std::vector<OutputTerminal<std::pair<int, pair<int, int>>, Node>> output_terminals)
            :
            input_terminals(input_terminals),
            output_terminals(output_terminals) {}
};



/**********************************************************************/
/* Struct  - Project
/* Used By - projectA_step, projectB_step   
/**********************************************************************/
struct Project : step_Base {

  double (*func)(double); 

  int execute(const std::pair<int, pair<int, int>> &node, CnCContext &context) const;

  Vector sValue(int n, int l, CnCContext &context) const;

  /*----------------------------------------------------------------*/
  /* Project Constructor
  /*----------------------------------------------------------------*/
  Project(
          double (*func)(double), 
          std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *> input_terminals,
          std::vector<OutputTerminal<std::pair<int,pair<int, int>>, Node>> output_terminals)
          :
          func(func),
          step_Base(input_terminals, output_terminals)
          {}
};


/**********************************************************************/
/* Struct  - Compress_Prolog
/* Used By - compress_prolog_FuncA_step, compress_prolog_FuncB_step
/**********************************************************************/
struct Compress_Prolog :step_Base{

  int execute(const std::pair< int, pair<int, int>> &node, CnCContext &context) const;

  /*----------------------------------------------------------------*/
  /* Compress_Prolog Constructor
  /*----------------------------------------------------------------*/
  Compress_Prolog( 
                  std::vector<CnC::item_collection<std::pair< int, pair<int, int>>, Node> *> input_terminals,
                  std::vector<OutputTerminal<std::pair< int, pair<int, int>>, Node>> output_terminals)
                  : 
                  step_Base(input_terminals, output_terminals) 
                  {}
};


/**********************************************************************/
/* Struct  - Compress_doIt
/* Used By - compress_doIt_funcA_step, compress_doIt_funcB_step
/**********************************************************************/
struct Compress_doIt :step_Base{

  int execute(const std::pair< int, pair<int, int>> &node, CnCContext &context) const;

  /*----------------------------------------------------------------*/
  /* Compress_doIt Constructor
  /*----------------------------------------------------------------*/
  Compress_doIt( 
                std::vector<CnC::item_collection<std::pair< int, pair<int, int>>, Node> *> input_terminals,
                std::vector<OutputTerminal<std::pair< int, pair<int, int>>, Node>> output_terminals)
                : 
                step_Base(input_terminals, output_terminals) 
                {}
};



/**********************************************************************/
/* Struct  - GaxpyOp
/* Used By - gaxpyOp_step
/**********************************************************************/
struct GaxpyOp :step_Base{

  double alpha;
  double beta;

  int execute(const std::pair< int, pair<int, int>> &node, CnCContext &context) const;
  /*----------------------------------------------------------------*/
  /* GaxpyOp Constructor
  /*----------------------------------------------------------------*/
  GaxpyOp( 
          double alpha,
          double beta,
          std::vector<CnC::item_collection<std::pair< int, pair<int, int>>, Node> *> input_terminals,
          std::vector<OutputTerminal<std::pair<int, pair<int, int>>, Node>> output_terminals)
          : 
          alpha( alpha ),
          beta( beta),
          step_Base(input_terminals, output_terminals) 
          {}
};



/**********************************************************************/
/* Struct  - Reconstruct_Prolog
/* Used By - reconstruct_prolog_step
/**********************************************************************/
struct Reconstruct_Prolog : step_Base{

  int execute(const std::pair<int, pair<int, int>> &node, CnCContext &context) const;

  /*----------------------------------------------------------------*/
  /* Reconstruct_Prolog Constructor
  /*----------------------------------------------------------------*/
  Reconstruct_Prolog( 
                    std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *> input_terminals,
                    std::vector<OutputTerminal<std::pair< int, pair<int, int>>, Node>> output_terminals)
                    : 
                    step_Base(input_terminals, output_terminals) 
                    {}
};


/**********************************************************************/
/* Struct  - Reconstruct_doIt
/* Used By - reconstruct_doIt_step
/**********************************************************************/
struct Reconstruct_doIt : step_Base{

  int execute(const std::pair< int,pair<int, int>> &node, CnCContext &context) const;

  /*----------------------------------------------------------------*/
  /* Reconstruct_doIt Constructor
  /*----------------------------------------------------------------*/
  Reconstruct_doIt( 
                  std::vector<CnC::item_collection<std::pair< int, pair<int, int>>, Node> *> input_terminals,
                  std::vector<OutputTerminal<std::pair< int, pair<int, int>>, Node>> output_terminals)
                  : 
                  step_Base(input_terminals, output_terminals) 
                  {}
};

/**********************************************************************/
/* Struct  - Norm2
/* Used By - norm2_step
/**********************************************************************/
struct Norm2 : step_Base{

  int execute(const std::pair<int, pair<int, int>> &node, CnCContext &context) const;

  /*----------------------------------------------------------------*/
  /* Norm2 Constructor
  /*----------------------------------------------------------------*/
  Norm2( 
        std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *> input_terminals,
        std::vector<OutputTerminal<std::pair<int, pair<int, int>>, Node>> output_terminals)
        : 
        step_Base(input_terminals, output_terminals) 
        {}
};

/**********************************************************************/
/* Struct  - Diff_Prolog
/* Used By - diff_prolog_step
/**********************************************************************/
struct Diff_Prolog :step_Base{

  int execute(const std::pair< int,pair<int, int>> &node, CnCContext &context) const;

  /*----------------------------------------------------------------*/
  /* Diff_Prolog Constructor
  /*----------------------------------------------------------------*/
  Diff_Prolog( 
              std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *> input_terminals,
              std::vector<OutputTerminal<std::pair<int, pair<int, int>>, Node>> output_terminals)
              : 
              step_Base(input_terminals, output_terminals) 
              {}
};

/**********************************************************************/
/* Struct  - Diff_doIt
/* Used By - Diff_doIt_step
/**********************************************************************/
struct Diff_doIt : step_Base {

  int execute(const std::pair<int, pair<int, int>> &node, CnCContext &context) const;

  Vector unfilter(const Vector &inputVector, int k, const Matrix * hg) const;

  /*----------------------------------------------------------------*/
  /* Diff_doIt Constructor
  /*----------------------------------------------------------------*/
  Diff_doIt(  
            std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *> input_terminals,
            std::vector<OutputTerminal<std::pair<int, pair<int, int>>, Node>> output_terminals)
            :
            step_Base(input_terminals, output_terminals)
            {}
};

/**********************************************************************/
/* Struct  - BinaryOp
/* Used By - subtract_1_step
/**********************************************************************/
struct BinaryOp : step_Base{

  using funcT = Vector (*)(const Vector &, const Vector&);
  funcT func;   

  double (*scale_factor)(double);

  Vector unfilter(Vector inputVector, int k, Matrix * hg) const;
  int execute(const std::pair< int, pair<int, int>> &node, CnCContext &context) const;

  /*----------------------------------------------------------------*/
  /* BinaryOp Constructor
  /*----------------------------------------------------------------*/
  BinaryOp(  
          const funcT &func, 
          double (*scale_factor)(double), 
          std::vector<CnC::item_collection<std::pair< int, pair<int, int>>, Node> *> input_terminals,
          std::vector<OutputTerminal<std::pair< int, pair<int, int>>, Node>> output_terminals)
          : 
          func(func), 
          scale_factor(scale_factor), 
          step_Base( input_terminals, output_terminals) 
          {}
};

/*********************************************************************
/* Struct  - Printer
/* Used By - printer_step
/**********************************************************************/
struct Printer : step_Base{

  int execute(const std::pair<int, pair<int, int>> &node, CnCContext &context) const;

  /*----------------------------------------------------------------*/
  /* Printer Constructor
  /*----------------------------------------------------------------*/
  Printer( 
          std::vector<CnC::item_collection<std::pair< int, pair<int, int>>, Node> *> input_terminals,
          std::vector<OutputTerminal<std::pair< int, pair<int, int>>, Node>> output_terminals)
          : 
          step_Base(input_terminals, output_terminals) 
          {}
};


/**********************************************************************/
/* Struct  - CnCContext
/**********************************************************************/
struct CnCContext : public CnC::context<CnCContext> {

public:
  int max_height;

  int k, quad_npt, max_level;

  double thresh;

  double (*a_function)(double);
  double (*b_function)(double);

  Matrix *hg, *hg0, *hg1, *hgT;
  Matrix *rm, *r0, *rp;

  Vector *quad_w, *quad_x;
  Matrix *quad_phi, *quad_phiT, *quad_phiw;

  tbb::concurrent_vector<double> inner_results;
  tbb::concurrent_vector<double> norm2_result;

  /*----------------------------------------------------------------*/
  /* CnCContext constructor
  /*----------------------------------------------------------------*/
  CnCContext(int k, double thresh, int max_level, int max_height)
  : 
  CnC::context<CnCContext>(), 
  max_height(max_height),
  k(k), 
  thresh(thresh),
  norm2_result(1, 0.0),
  max_level(max_level)
  {

    for( int i = 0, j = max_height; i < max_height; i++){
      height_map[i] = j;

      if( j == 1)
        j = max_height;
      else
        j--;
    }

    global_co_height = max_height;
    init_twoscale(k);
    init_quadrature(k);
    make_dc_periodic();
  }

  /*----------------------------------------------------------------*/
  /* init_twoscale
  /*----------------------------------------------------------------*/
  void init_twoscale(int k) {

    double  (*hgInput)[22] = twoscalecoeffs(k);

    hg = new Matrix(2*k, 2*k);
    hg0 = new Matrix(2*k, k);
    hg1 = new Matrix(2*k, k);
    hgT = new Matrix(2*k, 2*k);

    for (int i = 0; i < 2 * k; ++i) {
      for (int j = 0; j < 2 * k; ++j) {
        hg->set_item(i, j, hgInput[i][j]);
        hgT->set_item(i, j, hgInput[j][i]);
      }
    }

    for (int i = 0; i < 2 * k; ++i) {
      for (int j = 0; j < k; ++j) {
        hg0->set_item(i, j, hgInput[i][j]);
        hg1->set_item(i, j, hgInput[i][j+k]);
      }
    }
  }


  /*----------------------------------------------------------------*/
  /* init_quadrature
  /*----------------------------------------------------------------*/
  void init_quadrature(int order) {
    double *x = gauss_legendre_point(order);
    double *w = gauss_legendre_weight(order);

    quad_w = new Vector(w, 0, order);
    quad_x = new Vector(x, 0, order);

    int npt = order;
    quad_npt = npt;

    quad_phi = new Matrix(npt, k);
    quad_phiT = new Matrix(k, npt);
    quad_phiw = new Matrix(npt, k);

    for (int i = 0; i < npt; ++i) {
      double * p = phi((*quad_x)[i], k);
      for (int m = 0; m < k; ++m) {
        quad_phi->set_item(i, m, p[m]);
        quad_phiT->set_item(m, i, p[m]);
        quad_phiw->set_item(i, m, w[i] * p[m]);
      }
    }
  }


  /*----------------------------------------------------------------*/
  /* make_dc_periodic
  /*----------------------------------------------------------------*/
  void make_dc_periodic() {

    rm = new Matrix(k, k);
    r0 = new Matrix(k, k);
    rp = new Matrix(k, k);

    double iphase = 1.0;

    for (int i = 0; i < k; ++i) {
      double jphase = 1.0;

      for (int j = 0; j < k; ++j) {
        double gammaij = sqrt(( 2 * i + 1) * ( 2 * j + 1));
        double Kij;
        if ((( i -  j ) > 0) && (((i - j ) % 2) == 1 )) {
          Kij = 2.0;
        } 
        else {
          Kij = 0.0;
        }

        r0->set_item(i, j, (0.5 * (1.0 - iphase * jphase - 2.0 * Kij) * gammaij));
        rm->set_item(i, j, (0.5 * jphase * gammaij));
        rp->set_item(i, j, (-0.5 * iphase * gammaij));

        jphase = -1 * jphase;
      }
      iphase = -1 * iphase;
    }
  }

};


/**********************************************************************/
/* This method is used to instantiate the general struct BinaryOp to 
/* be specifically subtraction of two mathematical functions
/**********************************************************************/
Vector sub(const Vector &v1, const Vector &v2) {

  Vector result(v1);
  for (unsigned int i = 0; i < v2.length(); ++i) {
    result.data[i] -= v2.data[i];
  }
  return result;
}


/**********************************************************************/
/* This method is used to instantiate the general struct BinaryOp to 
/* be specifically multiplication of two mathematical functions
/**********************************************************************/
Vector mul(const Vector &v1, const Vector &v2) {

  Vector result(v1);
  for (unsigned int i = 0; i < v2.length(); ++i) {
    result.data[i] *= v2.data[i];
  }
  return result;
}

/**********************************************************************/
/* used in the method execute of class BinaryOp for the 
/* subtraction instance 
/**********************************************************************/
double sub_scale_factor(double n) {
   return 1.0;
}

/**********************************************************************/
/* used in the method execute of class BinaryOp for the 
/* multiplication instance 
/**********************************************************************/
double mul_scale_factor(double n) {
   return ( sqrt( pow ( 2.0, n) ) );
}

/**********************************************************************/
/* Mathematical test functions test1 and test2 would use 
/* the following guassian method
/**********************************************************************/
double gaussian(double x, double a, double coeff) {
    return coeff*exp(-a*x*x);
}



/**********************************************************************/
/* Project::sValue
/**********************************************************************/
Vector Project::sValue(int n, int l, CnCContext &context) const {

  Vector s(context.k);
  Vector &quad_x_ref = *(context.quad_x);
  Matrix &quad_phiw_ref = *(context.quad_phiw);

  double h = pow(0.5, n);
  double scale = sqrt(h);

  for (int mu = 0; mu < context.quad_npt; ++mu) {
    double x = (l + quad_x_ref[mu]) * h;
    double fValue = func(x);

    for (int i = 0; i < (context.k); ++i) {
      s[i] = s[i] + (scale * fValue * (quad_phiw_ref.get_item(mu, i)));
    }
  }
  return s;
}


/**********************************************************************/
/* Project::execute
/**********************************************************************/
int Project::execute(const std::pair< int, std::pair<int, int>> &root_node, CnCContext &context) const {

  queue<std::pair<int, std::pair<int, int>>> subtree_queue;
  subtree_queue.push(root_node);
  int k = context.k;

  //Current subtree will be executed inside this thread only
  while( !subtree_queue.empty() ){
    std::pair< int, std::pair<int, int>> curr_node = subtree_queue.front();
    subtree_queue.pop();

    int n, l, h, new_h;
    h = curr_node.first;
    n = curr_node.second.first;
    l = curr_node.second.second; 
    new_h = (h + 1) % context.max_height;

    //Get s0 and s1 for current node
    Vector s0 = sValue(n + 1, 2 * l, context);
    Vector s1 = sValue(n + 1, 2 * l + 1, context);

    // concatenation of s0 and s1
    Vector s(s0 | s1); 
    Vector d(s * (*(context.hgT)));

    //If we have reached threshold or maximum levels, we need to put the leaves
    if (d.normf(k, 2 * k) < context.thresh || n >= context.max_level - 1) {
      output_terminals[1].put(curr_node, Node(n, l, k, Vector(), Vector(), true));
      output_terminals[1].put(std::make_pair(new_h, std::make_pair(n + 1, 2 * l)), Node(n + 1, 2 * l, k, s0, Vector(), false));
      output_terminals[1].put(std::make_pair(new_h, std::make_pair(n + 1, 2 * l +1)), Node(n + 1, 2 * l + 1, k, s1, Vector(), false));
    }

    //If not maximum level or threshold continue making the tree
    else {
      output_terminals[1].put(curr_node, Node(n, l, k, Vector(), Vector(), true));

      if( new_h == 0 ){
        output_terminals[0].put(std::make_pair(new_h, std::make_pair(n + 1, 2 * l)));
        output_terminals[0].put(std::make_pair(new_h, std::make_pair(n + 1, 2 * l+1)));  
      }
      else{
        subtree_queue.push( std::make_pair(new_h, std::make_pair(n + 1, 2 * l)) );
        subtree_queue.push( std::make_pair(new_h, std::make_pair(n + 1, 2 * l+1)) );
      }
    }
  }

  return CnC::CNC_Success;
}


/**********************************************************************/
/* BinaryOp::unfilter
/**********************************************************************/
Vector BinaryOp::unfilter(Vector inputVector, int k, Matrix * hg) const {

  Vector inputVector_copy(inputVector);
  Vector vector_d(2 * k);
  vector_d.set_slice_from_another_vector(0, k, inputVector);
  Vector vector_s2 = (vector_d * (*hg));
  return vector_s2;
}

/**********************************************************************/
/* BinaryOp::execute
/**********************************************************************/
int BinaryOp::execute(const std::pair< int, pair<int, int>> &root_node, CnCContext &context) const {

  queue<std::pair<int, std::pair<int, int>>> subtree_queue;
  subtree_queue.push(root_node);
  int k = context.k;

  while( !subtree_queue.empty() ){

    std::pair< int, std::pair<int, int>> node = subtree_queue.front();
    subtree_queue.pop();

    Node left;
    Node right;
    int n, l, h, new_h;
    h = node.first;
    n = node.second.first;
    l = node.second.second; 
    new_h = (h + 1) % context.max_height;

    input_terminals[0]->get(node, left);
    input_terminals[1]->get(node, right);

    // If both of them are at the leaf level
    if (left.s.length() != 0 && right.s.length() != 0) { 
      double scale_fact = scale_factor(n);
      Vector f_vector(left.s * (*(context.quad_phiT)));
      Vector g_vector(right.s * (*(context.quad_phiT)));

      Vector temp = func(f_vector, g_vector);
      Vector resultVector((temp * (*context.quad_phiw)).scale(scale_fact));

      output_terminals[1].put(node, Node(n, l, k, resultVector, Vector(), false));
    }

    //Else unify them
    else {
      if (left.s.length() != 0) {
        Vector left_unfiltered = unfilter(left.s, k, context.hg);
        output_terminals[0].put(std::make_pair( new_h, make_pair(n + 1, 2 * l)), Node(n + 1, 2 * l, k, left_unfiltered.get_slice(0, k), Vector(), false));
        output_terminals[0].put(std::make_pair( new_h, make_pair(n + 1, 2 * l + 1)), Node(n + 1, 2 * l + 1, k, left_unfiltered.get_slice(k, 2 * k), Vector(), false));

      }
      else if (right.s.length() != 0) {
        Vector right_unfiltered = unfilter(right.s, k, context.hg);
        output_terminals[2].put(std::make_pair( new_h, make_pair(n + 1, 2 * l)), Node(n + 1, 2 * l, k, right_unfiltered.get_slice(0, k), Vector(), false));
        output_terminals[2].put(std::make_pair( new_h, make_pair(n + 1, 2 * l + 1)), Node(n + 1, 2 * l + 1, k, right_unfiltered.get_slice(k, 2 * k), Vector(), false));
      }

      output_terminals[1].put(node, Node(n, l, k, Vector(), Vector(), true));

      //add next level nodes in queue if it will not be called using thread
      if( new_h != 0 ){
        subtree_queue.push( std::make_pair(new_h, std::make_pair(n+1, 2*l)) );
        subtree_queue.push( std::make_pair(new_h, std::make_pair(n+1, 2*l+1)) );
      }
    }
  }

  return CnC::CNC_Success;
}

struct compress_output {
  pair<int, pair<int, int>> node;
  Node nodeInfo;
};

string key(int h, int n, int l){
  return to_string(h)+"+"+to_string(n)+"+"+to_string(l);
}
/**********************************************************************/
/* Compress_doIt::execute
/**********************************************************************/
int Compress_doIt::execute( const std::pair< int, pair<int, int>> &root_node, CnCContext &context ) const {

  // cout << "Tile started " << root_node.first << "  " << root_node.second.first << "  " << root_node.second.second << "\n";
  std::pair< int, pair<int, int>> tile_parent = make_pair(context.max_height - 1, make_pair(root_node.second.first-1, root_node.second.second/2));
  stack<std::pair< int, pair<int, int>>> subtree_stack;
  queue<std::pair< int, pair<int, int>>> levels;
  unordered_map<string, Node> mp_left, mp_right, mp_input;
  vector<struct compress_output> output;

  levels.push(root_node);
  int k = context.k;

  //put the inverted subtree in stack 
  while( !levels.empty() ){
    Node nodeInfo;
    std::pair< int, pair<int, int>> curr_node = levels.front();
    levels.pop();
    subtree_stack.push(curr_node);

    int n, l, h, new_h;
    h = curr_node.first;
    n = curr_node.second.first;
    l = curr_node.second.second;
    new_h = (h+1) % context.max_height;

    input_terminals[1]->get(curr_node, nodeInfo);
    string curr_key = key(h, n, l);
    mp_input[curr_key] = nodeInfo;

    //If curr node's children are supposed to be executed by this thread put them in queue
    if( new_h != 0 &&  nodeInfo.has_children ){
      levels.push( std::make_pair(new_h, make_pair(n+1, 2*l)) );
      levels.push( std::make_pair(new_h, make_pair(n+1, 2*l+1)) );
    }
  }

  //Process each node in stack 
  while( !subtree_stack.empty() ){

    Node nodeInfo, left, right;
    std::pair< int, pair<int, int>> node = subtree_stack.top();
    subtree_stack.pop();

    int n, l, h, new_h;
    h = node.first;
    n = node.second.first;
    l = node.second.second;
    new_h = h == 0 ? context.max_height - 1 : h-1;

    string node_key = key(h, n, l);
    string p_key = key(new_h, n-1, l/2);

    struct compress_output new_output;
    new_output.node = node;
    nodeInfo = mp_input[node_key];

    //current node is a leaf node
    if( !nodeInfo.has_children){

        new_output.nodeInfo = Node( n, l, k, Vector(), Vector(k), false );

        if( l & 0x1uL){
          mp_right[p_key] = Node( n-1, l/2, k, nodeInfo.s, Vector(), false );
        }
        else{
          mp_left[p_key] = Node( n-1, l/2, k, nodeInfo.s, Vector(), false );
        }
    }
    else {
      if( mp_left.find(node_key) != mp_left.end() ){
        left =  mp_left[node_key];
      } else {
        input_terminals[0]->get(node, left);
      }

      if( mp_right.find(node_key) != mp_right.end() ){
        right =  mp_right[node_key];
      } else {
        input_terminals[0]->get(node, right);
      }
      Vector s( left.s | right.s );
      Vector d(s * (*(context.hgT)));

      Vector sValue(d.data, 0, k);
      Vector dValue(d.data, k, 2 * k);

      //If at root node place both s and d
      if( n == 0){
        new_output.nodeInfo = Node( n, l, k, sValue, dValue, true);
      }
      else {
        new_output.nodeInfo = Node( n, l, k, Vector(), dValue, true);

        if( l & 0x1uL){
          mp_right[p_key] = Node( n-1, l/2, k, sValue, Vector(), true );
        }
        else{
          mp_left[p_key] = Node( n-1, l/2, k, sValue, Vector(), true );
        }   
      }
    }
    output.push_back(new_output);
  }
  
  //Put the results
  string tile_p_key = key(tile_parent.first, tile_parent.second.first, tile_parent.second.second);

  if( root_node.second.second & 0x1uL){
    output_terminals[2].put(tile_parent, mp_right[tile_p_key]);
  }
  else{
    output_terminals[0].put(tile_parent, mp_left[tile_p_key]);
  }   

  // cout << output.size() << "\n";
  for( int i = 0 ; i < output.size() ; i++){
    output_terminals[1].put(output[i].node, output[i].nodeInfo);
  }

  return CnC::CNC_Success;
}


/**********************************************************************/
/* GaxpyOp::execute
/**********************************************************************/
int GaxpyOp::execute( const std::pair< int, pair<int, int>> &root_node, CnCContext &context ) const {

  queue<std::pair<int, std::pair<int, int>>> subtree_queue;
  subtree_queue.push(root_node);
  int k = context.k;

  while( !subtree_queue.empty() ){
    std::pair< int, std::pair<int, int>> node = subtree_queue.front();
    subtree_queue.pop();

    Node left;
    Node right;

    int n, l, h, new_h;
    h = node.first;
    n = node.second.first;
    l = node.second.second;
    new_h = (h + 1) % context.max_height;

    input_terminals[0]->get(node, left);
    input_terminals[1]->get(node, right);

    Vector tempD(left.d);
    tempD.gaxpy( alpha, right.d, beta);
    Vector tempS;

    if( n == 0 && l == 0){
      tempS = left.s;
      tempS.gaxpy( alpha, right.s, beta);
    }

    output_terminals[1].put( node, Node( n, l, k, tempS, tempD, left.has_children || right.has_children ) );

    if( left.has_children && !right.has_children ){
      output_terminals[2].put( make_pair( new_h, make_pair(n + 1, l* 2)), Node(n+1, l * 2, k, Vector(), Vector(k), false));
      output_terminals[2].put( make_pair( new_h, make_pair(n + 1, l * 2 + 1)), Node(n+1, l*2+1, k, Vector(), Vector(k), false));
    }

    if( !left.has_children && right.has_children ){
      output_terminals[0].put( make_pair( new_h, make_pair(n + 1, l * 2)), Node(n+1, l*2, k, Vector(), Vector(k), false));
      output_terminals[0].put( make_pair( new_h, make_pair(n + 1, l * 2 + 1)), Node(n+1, l*2+1, k, Vector(), Vector(k), false));
    }

    if( new_h != 0 && (left.has_children || right.has_children)){
      subtree_queue.push( std::make_pair(new_h, std::make_pair(n+1, 2*l)) );
      subtree_queue.push( std::make_pair(new_h, std::make_pair(n+1, 2*l+1)) );
    }
  }

  return CnC::CNC_Success;
}


/**********************************************************************/
/* Reconstruct_Prolog::execute
/**********************************************************************/
int Reconstruct_Prolog::execute( const std::pair< int, pair<int, int>> &node, CnCContext &context ) const {

  Node input;
  if( node.second.first == 0) {
    input_terminals[0]->get( node, input);
    output_terminals[0].put( node, input);
  }
  return CnC::CNC_Success;
}

/**********************************************************************/
/* Reconstruct_doIt::execute
/**********************************************************************/
int Reconstruct_doIt::execute( const std::pair< int, pair<int, int>> &root_node, CnCContext &context ) const {

  queue<std::pair<int, std::pair<int, int>>> subtree_queue;
  subtree_queue.push(root_node);
  int k = context.k;

  while( !subtree_queue.empty() ){
    std::pair< int, std::pair<int, int>> node = subtree_queue.front();
    subtree_queue.pop();

    Node s_coeff;
    Node node_information;

    int n, l, h, new_h, out_height;
    h = node.first;
    n = node.second.first;
    l = node.second.second;
    new_h = (h+1) % context.max_height;

    input_terminals[0]->get( node, s_coeff);
    input_terminals[1]->get( node, node_information);

    Vector s = s_coeff.s;

    if( node_information.has_children ){

      Vector v1( s| node_information.d );
      Vector v2( v1 * (*context.hg) );

      Vector leftChildS(v2.data, 0, k);
      Vector rightChildS(v2.data, k, 2 * k);

      output_terminals[0].put(make_pair( new_h, make_pair(n + 1, l * 2)), Node( n + 1, l * 2, k, leftChildS, Vector(), false));
      output_terminals[0].put(make_pair( new_h, make_pair(n + 1, l * 2 + 1)), Node( n + 1, l * 2 + 1, k, rightChildS, Vector(), false));

      output_terminals[1].put( node, Node( n, l, k, Vector(), Vector(), true));

      if( new_h != 0 ){
        subtree_queue.push( std::make_pair(new_h, std::make_pair(n+1, 2*l)) );
        subtree_queue.push( std::make_pair(new_h, std::make_pair(n+1, 2*l+1)) );
      }
    }
    else{
      output_terminals[1].put( node, Node( n, l, k, s, Vector(), false));
    }
  }

  return CnC::CNC_Success;
}


/**********************************************************************/
/* Norm2::execute
/**********************************************************************/
int Norm2::execute( const std::pair<int, pair<int, int>> &root_node, CnCContext &context ) const {

  queue<std::pair<int, std::pair<int, int>>> subtree_queue;
  subtree_queue.push(root_node);

  while( !subtree_queue.empty() ){
    std::pair< int, std::pair<int, int>> node = subtree_queue.front();
    subtree_queue.pop();

    int n, l, h, new_h;
    h = node.first;
    n = node.second.first;
    l = node.second.second;
    new_h = (h+1) % context.max_height;

    Node nodeInfo;
    input_terminals[0]->get(node, nodeInfo);

    context.norm2_result[0] +=  pow( nodeInfo.s.normf(), 2);

    if( nodeInfo.has_children ){

      if( new_h == 0 ){
        output_terminals[0].put( std::make_pair( new_h, std::make_pair(n + 1, l * 2)));
        output_terminals[0].put( std::make_pair( new_h, std::make_pair(n + 1, l * 2 + 1)));
      }
      else{
        subtree_queue.push( std::make_pair(new_h, std::make_pair(n+1, 2*l )) );
        subtree_queue.push( std::make_pair(new_h, std::make_pair(n+1, 2*l+1)) );

      }
    }
  }
  return CnC::CNC_Success;
}

/**********************************************************************/
/* Diff_Prolog::execute
/**********************************************************************/
int Diff_Prolog::execute( const std::pair<int, pair<int, int>> &root_node, CnCContext &context ) const {

  queue<std::pair<int, std::pair<int, int>>> subtree_queue;
  subtree_queue.push(root_node);
  int k = context.k;

  while( !subtree_queue.empty() ){
    std::pair< int, std::pair<int, int>> node = subtree_queue.front();
    subtree_queue.pop();

    int n, l, h, new_h;
    h = node.first;
    n = node.second.first;
    l = node.second.second;
    new_h = (h+1) % context.max_height;

    Node nodeInfo;
    input_terminals[0]->get(node, nodeInfo);

    output_terminals[2].put( make_pair( h, make_pair( n, l==0ul ? (1ul<<n)-1 : l-1)), nodeInfo);
    output_terminals[1].put( node, nodeInfo);
    output_terminals[0].put( make_pair( h, make_pair( n, l==((1ul<<n)-1) ? 0 : l+1)), nodeInfo);

    if( new_h != 0 && nodeInfo.has_children){
      subtree_queue.push( std::make_pair(new_h, std::make_pair(n+1, 2*l)) );
      subtree_queue.push( std::make_pair(new_h, std::make_pair(n+1, 2*l+1)) );
    }
  }

  return CnC::CNC_Success;
}

/**********************************************************************/
/* Diff_doIt::unfilter
/**********************************************************************/
Vector Diff_doIt::unfilter(const Vector &inputVector, int k, const Matrix * hg) const {
   
  Vector vector_d(2 * k);
  vector_d.set_slice_from_another_vector(0, k, inputVector);
  Vector vector_s2 = (vector_d * (*hg));
  return vector_s2;
}


/**********************************************************************/
/* Diff_doIt::execute
/**********************************************************************/
int Diff_doIt::execute( const std::pair<int, pair<int, int>> &root_node, CnCContext &context ) const {

  queue<std::pair<int, std::pair<int, int>>> subtree_queue;
  subtree_queue.push(root_node);
  int k = context.k;


  while( !subtree_queue.empty() ){
    std::pair< int, std::pair<int, int>> node = subtree_queue.front();
    subtree_queue.pop();

    int n, l, h, new_h;
    h = node.first;
    n = node.second.first;
    l = node.second.second;
    new_h = (h+1) % context.max_height;

    Node left, center, right;
    input_terminals[0]->get(node, left);
    input_terminals[1]->get(node, center);
    input_terminals[2]->get(node, right);

    if ((left.s.length() != 0 && left.has_children) || (left.s.length() == 0 && !left.has_children)) {
      std::cout << "ERROR at left" << std::endl;
    }

    if ((center.s.length() != 0 && center.has_children)|| (center.s.length() == 0 && !center.has_children)) {
      std::cout << "ERROR at center" << std::endl;
    }

    if ((right.s.length() != 0 && right.has_children) ||(right.s.length() == 0 && !right.has_children)) {
      std::cout << "ERROR at right" << std::endl;
    }


    if (left.s.length() != 0 && center.s.length() != 0 && right.s.length() != 0) {
      Vector r = ((*context.rp) * left.s) + ((*context.r0) * center.s) + ((*context.rm) * right.s);
      output_terminals[3].put(node, Node( n, l, k, r.scale( pow(2.0, n)), Vector(), false));
    }
    else {
      output_terminals[3].put(node, Node( n, l, k, Vector(), Vector(), true));

       if (left.s.length() != 0) {   
          Vector unfiltered = unfilter(left.s, k, context.hg);
          output_terminals[0].put(make_pair( new_h, make_pair( n+1, l*2)), Node( n +1, (l-1)*2 +1, k, unfiltered.get_slice(k, 2 * k), Vector(), false));
       }
         
       if (center.s.length() != 0) {

         Vector unfiltered = unfilter(center.s, k, context.hg);
         output_terminals[2].put( make_pair( new_h, make_pair( n+1, l*2)), Node( n +1, l*2 +1, k, unfiltered.get_slice(k, 2 * k), Vector(), false));
         output_terminals[0].put( make_pair( new_h,  make_pair( n+1, l*2+1)), Node( n +1, l*2, k, unfiltered.get_slice(0,k), Vector(), false));
         output_terminals[1].put( make_pair( new_h, make_pair( n+1, l*2)), Node( n +1, l*2, k, unfiltered.get_slice(0,k), Vector(), false));
         output_terminals[1].put( make_pair( new_h, make_pair( n+1, l*2+1)), Node( n +1, l*2+1, k, unfiltered.get_slice(k,2*k), Vector(), false));
       }

       if (right.s.length() != 0) {
          Vector unfiltered = unfilter(right.s, k, context.hg);
          output_terminals[2].put( make_pair( new_h, make_pair( n+1, l*2+1)), Node( n +1, (l+1)*2 , k, unfiltered.get_slice(0,k), Vector(), false));
       }

        if( new_h != 0 ){
        subtree_queue.push( std::make_pair(new_h, std::make_pair(n+1, 2*l)) );
        subtree_queue.push( std::make_pair(new_h, std::make_pair(n+1, 2*l+1)) );
      }       
    }
  }
  
  return CnC::CNC_Success;
}


/**********************************************************************/
/* Printer::execute
/**********************************************************************/
int Printer::execute(const std::pair<int, pair<int, int>> &node, CnCContext &context) const {

  queue<std::pair<int, std::pair<int, int>>> subtree_queue;
  subtree_queue.push(node);
  int k = context.k;

  while( !subtree_queue.empty() ){
    std::pair< int, std::pair<int, int>> curr_node = subtree_queue.front();
    subtree_queue.pop();

    int n, l, h, new_h;
    h = curr_node.first;
    n = curr_node.second.first;
    l = curr_node.second.second; 
    new_h = (h + 1) % context.max_height;
    
    Node nodeInfo;
    input_terminals[0]->get(curr_node, nodeInfo);

    // if( !nodeInfo.has_children )  
    std::cout << "Printer:: Node with info: (Key:"
              << h << " (" 
              << n << ", " 
              << l << "), " 
              << nodeInfo.toString() 
              << ")" 
              << std::endl;

    if( nodeInfo.has_children && new_h != 0 ){ 

      subtree_queue.push( std::make_pair(new_h, std::make_pair(n + 1, 2 * l)) );
      subtree_queue.push( std::make_pair(new_h, std::make_pair(n + 1, 2 * l+1)) );
    }
  }
  
  return CnC::CNC_Success;
}



/**********************************************************************/
/* Differentiation test
/**********************************************************************/
struct diff_test: CnCContext{

  /*----------------------------------------------------------------*/
  /* Item Collections
  /*----------------------------------------------------------------*/
  CnC::item_collection<std::pair<int, pair<int, int>>, Node> project_item;
  CnC::item_collection<std::pair<int, pair<int, int>>, Node> diff_left_item;
  CnC::item_collection<std::pair<int, pair<int, int>>, Node> diff_right_item;
  CnC::item_collection<std::pair<int, pair<int, int>>, Node> diff_centre_item;
  CnC::item_collection<std::pair<int, pair<int, int>>, Node> diff_result_item;

  CnC::tag_collection<std::pair<int, pair<int, int>>> project_tag;
  CnC::tag_collection<std::pair<int, pair<int, int>>> diff_prolog_tag;
  CnC::tag_collection<std::pair<int, pair<int, int>>> diff_doIt_tag;
  CnC::tag_collection<std::pair<int, pair<int, int>>> printer_tag;
  CnC::tag_collection<std::pair<int, pair<int, int>>> norm2_tag;

  using OutputTerminalType = OutputTerminal<std::pair<int, pair<int, int>>, Node>;

  CnC::step_collection<Project> project_step;
  CnC::step_collection<Diff_Prolog> diff_prolog_step;
  CnC::step_collection<Diff_doIt> diff_doIt_step;
  CnC::step_collection<Printer> printer_step;
  CnC::step_collection<Norm2> norm2_step;

  diff_test( int k, double thresh, int max_level,  int height, double (*funcA)(double))
  :
  CnCContext( k, thresh, max_level, height),
  project_item(*this),
  diff_left_item(*this),
  diff_right_item(*this),
  diff_centre_item(*this),
  diff_result_item(*this),
  project_tag(*this),
  diff_prolog_tag(*this),
  diff_doIt_tag(*this),
  printer_tag(*this),
  norm2_tag(*this),

  /*----------------------------------------------------------------*/
  /* Declare project_step
  /*----------------------------------------------------------------*/

  project_step(
                *this, 
                "project_step", 
                Project(
                        funcA, 
                        std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *>{},
                        std::vector<OutputTerminalType> {
                            OutputTerminalType(nullptr, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&project_tag}),
                            OutputTerminalType(&project_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&diff_prolog_tag})}
                        )
                ),

  /*----------------------------------------------------------------*/
  /* Declare diff_prolog_step
  /*----------------------------------------------------------------*/

  diff_prolog_step(
                *this, 
                "diff_prolog_step", 
                Diff_Prolog(
                        std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *>{ &project_item },
                        std::vector<OutputTerminalType> {
                            OutputTerminalType(&diff_left_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {}),
                            OutputTerminalType(&diff_centre_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&diff_doIt_tag}),
                            OutputTerminalType(&diff_right_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {})}
                        )
                ),

  /*----------------------------------------------------------------*/
  /* Declare diff_doIt_step
  /*----------------------------------------------------------------*/

  diff_doIt_step(
                *this, 
                "diff_doIt_step", 
                Diff_doIt(
                        std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *>{ &diff_left_item, &diff_centre_item, &diff_right_item },
                        std::vector<OutputTerminalType> {
                            OutputTerminalType(&diff_left_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {}),
                            OutputTerminalType(&diff_centre_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&diff_doIt_tag}),
                            OutputTerminalType(&diff_right_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {}),
                            OutputTerminalType(&diff_result_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&printer_tag})}
                        )
                ),
 
  /*----------------------------------------------------------------*/
  /* Declare norm2_f1_step
  /*----------------------------------------------------------------*/

  norm2_step( 
                *this, 
                "norm2_step", 
                Norm2( 
                    std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *> {&project_item}, 
                    std::vector<OutputTerminalType>{
                      OutputTerminalType(nullptr, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&norm2_tag})})
               ),

  /*----------------------------------------------------------------*/
  /* Declare printer_step
  /*----------------------------------------------------------------*/
  printer_step(
              *this, 
              "printer_step", 
              Printer( 
                std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *>{&diff_result_item}, 
                std::vector<OutputTerminalType>{})
              )
  {

    /*----------------------------------------------------------------*/
    /* Tag Prescription
    /*----------------------------------------------------------------*/
    project_tag.prescribes(project_step, *this);
    diff_prolog_tag.prescribes(diff_prolog_step, *this);
    diff_doIt_tag.prescribes(diff_doIt_step, *this);
    printer_tag.prescribes(printer_step, *this);
    norm2_tag.prescribes(norm2_step, *this);

    /*----------------------------------------------------------------*/
    /* Steps Produce Consume
    /*----------------------------------------------------------------*/

    project_step.produces(project_item);

    diff_prolog_step.consumes(project_item);
    diff_prolog_step.produces(diff_left_item);
    diff_prolog_step.produces(diff_centre_item);
    diff_prolog_step.produces(diff_right_item);

    diff_doIt_step.consumes(diff_left_item);
    diff_doIt_step.consumes(diff_centre_item);
    diff_doIt_step.consumes(diff_right_item);

    diff_doIt_step.produces(diff_left_item);
    diff_doIt_step.produces(diff_centre_item);
    diff_doIt_step.produces(diff_right_item);
    diff_doIt_step.produces(diff_result_item);

    printer_step.consumes(diff_result_item);

    norm2_step.consumes(project_item);
  }

};


/**********************************************************************/
/* Addition test 
/**********************************************************************/
struct addition_test : CnCContext {
 
  /*----------------------------------------------------------------*/
  /* Item Collections
  /*----------------------------------------------------------------*/
   CnC::item_collection<std::pair<int, pair<int, int>>, Node> projectA_item;
   CnC::item_collection<std::pair<int, pair<int, int>>, Node> projectB_item; 
  

   CnC::item_collection<std::pair<int, pair<int, int>>, Node> compressA_left_item;
   CnC::item_collection<std::pair<int, pair<int, int>>, Node> compressA_right_item;
   CnC::item_collection<std::pair<int, pair<int, int>>, Node> funcA_coeff_compressed_item;
   CnC::item_collection<std::pair<int, pair<int, int>>, Node> compressB_left_item;
   CnC::item_collection<std::pair<int, pair<int, int>>, Node> compressB_right_item;
   CnC::item_collection<std::pair<int, pair<int, int>>, Node> funcB_coeff_compressed_item;
   CnC::item_collection<std::pair<int, pair<int, int>>, Node> gaxpy_result_item;
   CnC::item_collection<std::pair<int, pair<int, int>>, Node> reconstruct_result_item;
   CnC::item_collection<std::pair<int, pair<int, int>>, Node> s_coeff_item;


  /*----------------------------------------------------------------*/
  /* Tag Collections
  /*----------------------------------------------------------------*/
   CnC::tag_collection<std::pair<int, pair<int, int>>> projectA_tag;
   CnC::tag_collection<std::pair<int, pair<int, int>>> projectB_tag;
   CnC::tag_collection<std::pair<int, pair<int, int>>> printer_tag;


   CnC::tag_collection<std::pair<int, pair<int, int>>> compress_doIt_funcA_tag;
   CnC::tag_collection<std::pair<int, pair<int, int>>> compress_doIt_funcB_tag;
   CnC::tag_collection<std::pair<int, pair<int, int>>> gaxpyOP_tag;
   CnC::tag_collection<std::pair<int, pair<int, int>>> reconstruct_prolog_tag;
   CnC::tag_collection<std::pair<int, pair<int, int>>> reconstruct_doIt_tag;

   CnC::tag_collection<std::pair<int, pair<int, int>>> norm2_f1_tag;
   CnC::tag_collection<std::pair<int, pair<int, int>>> norm2_f2_tag;
   CnC::tag_collection<std::pair<int, pair<int, int>>> norm2_add_tag;


  /*----------------------------------------------------------------*/
  /* Step Collections
  /*----------------------------------------------------------------*/
   using OutputTerminalType = OutputTerminal<std::pair<int, pair<int, int>>, Node>;

   CnC::step_collection<Project> projectA_step;
   CnC::step_collection<Project> projectB_step;
   CnC::step_collection<Printer> printer_step;

   CnC::step_collection<Compress_doIt> compress_doIt_funcA_step;
   CnC::step_collection<Compress_doIt> compress_doIt_funcB_step;
   CnC::step_collection<GaxpyOp> gaxpyOp_step;
   CnC::step_collection<Reconstruct_Prolog> reconstruct_prolog_step;
   CnC::step_collection<Reconstruct_doIt>  reconstruct_doIt_step;


   CnC::step_collection<Norm2> norm2_f1_step;
   CnC::step_collection<Norm2> norm2_f2_step;
   CnC::step_collection<Norm2> norm2_add_step;

  /*----------------------------------------------------------------*/
  /* addition_test Constructor
  /*----------------------------------------------------------------*/
   addition_test(int k, double thresh, int max_level, int height, double (*funcA)(double), double (*funcB)(double))
   : 
     CnCContext( k, thresh, max_level, height),

     projectA_item(*this),
     projectB_item(*this), 
     projectA_tag(*this), 
     projectB_tag(*this), 
     printer_tag(*this), 
     
     compressA_left_item(*this),
     compressA_right_item(*this),
     funcA_coeff_compressed_item(*this),
     compressB_left_item(*this),
     compressB_right_item(*this),
     funcB_coeff_compressed_item(*this),
     gaxpy_result_item(*this),
     reconstruct_result_item(*this),
     s_coeff_item(*this),

     compress_doIt_funcA_tag(*this),
     compress_doIt_funcB_tag(*this),
     gaxpyOP_tag(*this),
     reconstruct_prolog_tag(*this),
     reconstruct_doIt_tag(*this),

     norm2_f1_tag(*this),
     norm2_f2_tag(*this),
     norm2_add_tag(*this),
     /*----------------------------------------------------------------*/
     /* Declare projectA_step
     /*----------------------------------------------------------------*/
     projectA_step(
                  *this, 
                  "projectA_step", 
                  Project(
                          funcA, 
                          std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *>{},
                          std::vector<OutputTerminalType> {
                              OutputTerminalType(nullptr, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&projectA_tag}),
                              OutputTerminalType(&projectA_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&compress_doIt_funcA_tag})}
                          )
                  ),
 
     /*----------------------------------------------------------------*/
     /* Declare projectB_step
     /*----------------------------------------------------------------*/
     projectB_step(
                  *this, 
                  "projectB_step", 
                  Project(
                          funcB, 
                          std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *>{},
                          std::vector<OutputTerminalType> {
                              OutputTerminalType(nullptr, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&projectB_tag}),
                              OutputTerminalType(&projectB_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&compress_doIt_funcB_tag})}
                          )
                  ),
    
  
    /*----------------------------------------------------------------*/
    /* Declare printer_step
    /*----------------------------------------------------------------*/
    printer_step(
                *this, 
                "printer_step", 
                Printer( 
                  std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *>{&reconstruct_result_item}, 
                  std::vector<OutputTerminalType>{})
                ),

    /*----------------------------------------------------------------*/
    /* Declare compress_doIt_funcA_step
    /*----------------------------------------------------------------*/
    compress_doIt_funcA_step ( 
                                *this, 
                                "compress_doIt_funcA_step", 
                                Compress_doIt( 
                                    std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *> {&compressA_left_item, &projectA_item, &compressA_right_item}, 
                                    std::vector<OutputTerminalType>{
                                    OutputTerminalType(&compressA_left_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {}),
                                    OutputTerminalType(&funcA_coeff_compressed_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&gaxpyOP_tag}),
                                    OutputTerminalType(&compressA_right_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {})})
                              ),


    /*----------------------------------------------------------------*/
    /* Declare compress_doIt_funcB_step
    /*----------------------------------------------------------------*/
    compress_doIt_funcB_step ( 
                                *this, 
                                "compress_doIt_funcB_step", 
                                Compress_doIt( 
                                    std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *> {&compressB_left_item, &projectB_item, &compressB_right_item}, 
                                    std::vector<OutputTerminalType>{
                                    OutputTerminalType(&compressB_left_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {}),
                                    OutputTerminalType(&funcB_coeff_compressed_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {}),
                                    OutputTerminalType(&compressB_right_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {})})
                              ),


    /*----------------------------------------------------------------*/
    /* Declare gaxpyOp_step
    /*----------------------------------------------------------------*/
    gaxpyOp_step ( 
                  *this, 
                  "gaxpyOp_step", 
                  GaxpyOp( 
                      1.0,
                      1.0,
                      std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *> {&funcA_coeff_compressed_item, &funcB_coeff_compressed_item}, 
                      std::vector<OutputTerminalType>{
                      OutputTerminalType(&funcA_coeff_compressed_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&gaxpyOP_tag}),
                      OutputTerminalType(&gaxpy_result_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&reconstruct_prolog_tag}),
                      OutputTerminalType(&funcB_coeff_compressed_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {})})
                ),

    /*----------------------------------------------------------------*/
    /* Declare reconstruct_prolog_step
    /*----------------------------------------------------------------*/
    reconstruct_prolog_step ( 
                            *this, 
                            "reconstruct_prolog_step", 
                            Reconstruct_Prolog( 
                                std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *> {&gaxpy_result_item}, 
                                std::vector<OutputTerminalType>{
                                OutputTerminalType(&s_coeff_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&reconstruct_doIt_tag})})
                          ),

  
    /*----------------------------------------------------------------*/
    /* Declare reconstruct_doIt_step
    /*----------------------------------------------------------------*/
    reconstruct_doIt_step( 
                          *this, 
                          "reconstruct_doIt_step", 
                          Reconstruct_doIt( 
                              std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *> {&s_coeff_item, &gaxpy_result_item}, 
                              std::vector<OutputTerminalType>{
                                OutputTerminalType(&s_coeff_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&reconstruct_doIt_tag}),
                                OutputTerminalType(&reconstruct_result_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&printer_tag})})
                         ),

    /*----------------------------------------------------------------*/
    /* Declare norm2_f1_step
    /*----------------------------------------------------------------*/

    norm2_f1_step( 
                  *this, 
                  "norm2_f1_step", 
                  Norm2( 
                      std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *> {&projectA_item}, 
                      std::vector<OutputTerminalType>{
                        OutputTerminalType(nullptr, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&norm2_f1_tag})})
                 ),

    /*----------------------------------------------------------------*/
    /* Declare norm2_f2_step
    /*----------------------------------------------------------------*/

    norm2_f2_step( 
                  *this, 
                  "norm2_f2_step", 
                  Norm2( 
                      std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *> {&projectB_item}, 
                      std::vector<OutputTerminalType>{
                        OutputTerminalType(nullptr, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&norm2_f2_tag})})
                 ),


    /*----------------------------------------------------------------*/
    /* Declare norm2_add_step
    /*----------------------------------------------------------------*/

    norm2_add_step( 
                  *this, 
                  "norm2_add_step", 
                  Norm2( 
                      std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *> {&reconstruct_result_item}, 
                      std::vector<OutputTerminalType>{
                        OutputTerminalType(nullptr, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&norm2_add_tag})})
                 )


    {
      
      /*----------------------------------------------------------------*/
      /* Tag Prescription
      /*----------------------------------------------------------------*/
      projectA_tag.prescribes(projectA_step, *this);
      projectB_tag.prescribes(projectB_step, *this);
      printer_tag.prescribes(printer_step, *this);
      compress_doIt_funcA_tag.prescribes( compress_doIt_funcA_step, *this);
      compress_doIt_funcB_tag.prescribes( compress_doIt_funcB_step, *this);
      gaxpyOP_tag.prescribes( gaxpyOp_step, *this);
      reconstruct_prolog_tag.prescribes( reconstruct_prolog_step, *this);
      reconstruct_doIt_tag.prescribes( reconstruct_doIt_step, *this);
      norm2_f1_tag.prescribes( norm2_f1_step, *this);
      norm2_f2_tag.prescribes( norm2_f2_step, *this);
      norm2_add_tag.prescribes( norm2_add_step, *this);


      /*----------------------------------------------------------------*/
      /* Steps Produce Consume
      /*----------------------------------------------------------------*/
      projectA_step.produces(projectA_item);

      projectB_step.produces(projectB_item);

      compress_doIt_funcA_step.consumes( compressA_left_item );
      compress_doIt_funcA_step.consumes( compressA_right_item );
      compress_doIt_funcA_step.produces( funcA_coeff_compressed_item);
      compress_doIt_funcA_step.produces( compressA_left_item);
      compress_doIt_funcA_step.produces( compressA_right_item);

      compress_doIt_funcB_step.consumes( compressB_left_item );
      compress_doIt_funcB_step.consumes( compressB_right_item );
      compress_doIt_funcB_step.produces( funcB_coeff_compressed_item);
      compress_doIt_funcB_step.produces( compressB_left_item );
      compress_doIt_funcB_step.produces( compressB_right_item );

      gaxpyOp_step.consumes( funcA_coeff_compressed_item);
      gaxpyOp_step.consumes( funcB_coeff_compressed_item);
      gaxpyOp_step.produces( funcA_coeff_compressed_item);
      gaxpyOp_step.produces( funcB_coeff_compressed_item);
      gaxpyOp_step.produces( gaxpy_result_item);

      reconstruct_prolog_step.consumes( gaxpy_result_item );
      reconstruct_prolog_step.produces( s_coeff_item);

      reconstruct_doIt_step.consumes( s_coeff_item);
      reconstruct_doIt_step.consumes( gaxpy_result_item);
      reconstruct_doIt_step.produces( s_coeff_item);
      reconstruct_doIt_step.produces( reconstruct_result_item);

      printer_step.consumes(reconstruct_result_item);

      norm2_f1_step.consumes( projectA_item);
      norm2_f2_step.consumes( projectB_item);
      norm2_add_step.consumes( reconstruct_result_item);
   }
};


/**********************************************************************/
/* multiplication test
/**********************************************************************/
struct multiplication_test: CnCContext{
   
  /*----------------------------------------------------------------*/
  /* Item Collections
  /*----------------------------------------------------------------*/
  CnC::item_collection<std::pair<int, pair<int, int>>, Node> projectA_item;
  CnC::item_collection<std::pair<int, pair<int, int>>, Node> projectB_item; 
  CnC::item_collection<std::pair<int, pair<int, int>>, Node> mul_item;
  
  /*----------------------------------------------------------------*/
  /* Tag Collections
  /*----------------------------------------------------------------*/
  CnC::tag_collection<std::pair<int, pair<int, int>>> projectA_tag;
  CnC::tag_collection<std::pair<int, pair<int, int>>> projectB_tag;
  CnC::tag_collection<std::pair<int, pair<int, int>>> mul_tag;
  CnC::tag_collection<std::pair<int, pair<int, int>>> printer_tag;
  CnC::tag_collection<std::pair<int, pair<int, int>>> norm2_f1_tag;
  CnC::tag_collection<std::pair<int, pair<int, int>>> norm2_f2_tag;
  CnC::tag_collection<std::pair<int, pair<int, int>>> norm2_mul_tag;

  /*----------------------------------------------------------------*/
  /* Step Collections
  /*----------------------------------------------------------------*/
  using OutputTerminalType = OutputTerminal<std::pair<int, pair<int, int>>, Node>;

  CnC::step_collection<Project> projectA_step;
  CnC::step_collection<Project> projectB_step;
  CnC::step_collection<BinaryOp> mul_step;
  CnC::step_collection<Printer> printer_step;
  CnC::step_collection<Norm2> norm2_f1_step;
  CnC::step_collection<Norm2> norm2_f2_step;
  CnC::step_collection<Norm2> norm2_mul_step;

  /*----------------------------------------------------------------*/
  /* CnCContext Constructor
  /*----------------------------------------------------------------*/
  multiplication_test(int k, double thresh, int max_level, int height, double (*funcA)(double), double (*funcB)(double))
  : 
    CnCContext( k, thresh, max_level, height),
    projectA_item(*this),
    projectB_item(*this), 
    mul_item(*this), 
    projectA_tag(*this), 
    projectB_tag(*this), 
    mul_tag(*this),
    printer_tag(*this), 
    norm2_f1_tag(*this),
    norm2_f2_tag(*this),
    norm2_mul_tag(*this),

    /*----------------------------------------------------------------*/
    /* Declare projectA_step
    /*----------------------------------------------------------------*/
    projectA_step(
                  *this, 
                  "projectA_step", 
                  Project(
                          funcA, 
                          std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *>{},
                          std::vector<OutputTerminalType> {
                              OutputTerminalType(nullptr, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&projectA_tag}),
                              OutputTerminalType(&projectA_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&mul_tag})}
                          )
                 ),

    /*----------------------------------------------------------------*/
    /* Declare projectB_step
    /*----------------------------------------------------------------*/
    projectB_step(
                *this, 
                "projectB_step", 
                Project(
                        funcB, 
                        std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *>{},
                        std::vector<OutputTerminalType> {
                            OutputTerminalType(nullptr, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&projectB_tag}),
                            OutputTerminalType(&projectB_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {})}
                        )
                ),
    
    /*----------------------------------------------------------------*/
    /* Declare mul_step
    /*----------------------------------------------------------------*/
    mul_step(
            *this, 
            "mul_step", 
            BinaryOp(
                    &mul, 
                    &mul_scale_factor, 
                    std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *> {&projectA_item, &projectB_item},
                    std::vector<OutputTerminalType>{
                        OutputTerminalType(&projectA_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&mul_tag}),
                        OutputTerminalType(&mul_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&printer_tag}),
                        OutputTerminalType(&projectB_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {})}
                    )
            ),

    /*----------------------------------------------------------------*/
    /* Declare printer_step
    /*----------------------------------------------------------------*/
    printer_step(
                *this, 
                "printer_step", 
                Printer( 
                  std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *>{&mul_item}, 
                  std::vector<OutputTerminalType>{})
                ),


    /*----------------------------------------------------------------*/
    /* Declare norm2_f1_step
    /*----------------------------------------------------------------*/

    norm2_f1_step( 
                  *this, 
                  "norm2_f1_step", 
                  Norm2( 
                      std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *> {&projectA_item}, 
                      std::vector<OutputTerminalType>{
                        OutputTerminalType(nullptr, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&norm2_f1_tag})})
                 ),

    /*----------------------------------------------------------------*/
    /* Declare norm2_f2_step
    /*----------------------------------------------------------------*/

    norm2_f2_step( 
                  *this, 
                  "norm2_f2_step", 
                  Norm2( 
                      std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *> {&projectB_item}, 
                      std::vector<OutputTerminalType>{
                        OutputTerminalType(nullptr, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&norm2_f2_tag})})
                 ),


    /*----------------------------------------------------------------*/
    /* Declare norm2_add_step
    /*----------------------------------------------------------------*/

    norm2_mul_step( 
                  *this, 
                  "norm2_add_step", 
                  Norm2( 
                      std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *> {&mul_item}, 
                      std::vector<OutputTerminalType>{
                        OutputTerminalType(nullptr, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&norm2_mul_tag})})
                 )


    
    {
      
      /*----------------------------------------------------------------*/
      /* Tag Prescription
      /*----------------------------------------------------------------*/
      projectA_tag.prescribes(projectA_step, *this);
      projectB_tag.prescribes(projectB_step, *this);
      mul_tag.prescribes(mul_step, *this);
      printer_tag.prescribes(printer_step, *this);
      norm2_f1_tag.prescribes(norm2_f1_step, *this);
      norm2_f2_tag.prescribes(norm2_f2_step, *this);
      norm2_mul_tag.prescribes(norm2_mul_step, *this);

      /*----------------------------------------------------------------*/
      /* Steps Produce Consume
      /*----------------------------------------------------------------*/
      projectA_step.produces(projectA_item);

      projectB_step.produces(projectB_item);

      mul_step.consumes(projectA_item);
      mul_step.consumes(projectB_item);
      mul_step.produces(projectA_item);
      mul_step.produces(projectB_item);
      mul_step.produces(mul_item);

      printer_step.consumes(mul_item);

      norm2_f1_step.consumes( projectA_item);
      norm2_f2_step.consumes( projectB_item);
      norm2_mul_step.consumes( mul_item);

    }
};


/**********************************************************************/
/* gaxpy_sub_test test
/**********************************************************************/
struct gaxpy_sub_test: CnCContext{

  /*----------------------------------------------------------------*/
  /* Item Collections
  /*----------------------------------------------------------------*/
   CnC::item_collection<std::pair< int, pair<int, int>>, Node> projectA_item;
   CnC::item_collection<std::pair< int, pair<int, int>>, Node> projectB_item; 
   CnC::item_collection<std::pair< int, pair<int, int>>, Node> subtract_1_item; 
   CnC::item_collection<std::pair< int, pair<int, int>>, Node> subtract_2_item; 

   CnC::item_collection<std::pair< int, pair<int, int>>, Node> compressA_left_item;
   CnC::item_collection<std::pair< int, pair<int, int>>, Node> compressA_right_item;
   CnC::item_collection<std::pair< int, pair<int, int>>, Node> funcA_coeff_compressed_item;
   CnC::item_collection<std::pair< int, pair<int, int>>, Node> compressB_left_item;
   CnC::item_collection<std::pair< int, pair<int, int>>, Node> compressB_right_item;
   CnC::item_collection<std::pair< int, pair<int, int>>, Node> funcB_coeff_compressed_item;
   CnC::item_collection<std::pair< int, pair<int, int>>, Node> gaxpy_result_item;
   CnC::item_collection<std::pair< int, pair<int, int>>, Node> reconstruct_result_item;
   CnC::item_collection<std::pair< int, pair<int, int>>, Node> s_coeff_item;

  using OutputTerminalType = OutputTerminal<std::pair<int, pair<int, int>>, Node>;

  /*----------------------------------------------------------------*/
  /* Tag Collections
  /*----------------------------------------------------------------*/
   CnC::tag_collection<std::pair< int, pair<int, int>>> projectA_tag;
   CnC::tag_collection<std::pair< int, pair<int, int>>> printer_tag;

   CnC::tag_collection<std::pair< int, pair<int, int>>> projectB_tag;
   CnC::tag_collection<std::pair< int, pair<int, int>>> subtract_1_tag;
   CnC::tag_collection<std::pair< int, pair<int, int>>> subtract_2_tag;
   CnC::tag_collection<std::pair< int, pair<int, int>>> compress_doIt_funcA_tag;
   CnC::tag_collection<std::pair< int, pair<int, int>>> compress_doIt_funcB_tag;
   CnC::tag_collection<std::pair< int, pair<int, int>>> gaxpyOP_tag;
   CnC::tag_collection<std::pair< int, pair<int, int>>> reconstruct_prolog_tag;
   CnC::tag_collection<std::pair< int, pair<int, int>>> reconstruct_doIt_tag;
   CnC::tag_collection<std::pair<int, pair<int, int>>> norm2_tag;

   /*----------------------------------------------------------------*/
  /* Step Collections
  /*----------------------------------------------------------------*/

   CnC::step_collection<Project> projectA_step;
   CnC::step_collection<Printer> printer_step;

   CnC::step_collection<Project> projectB_step;
   CnC::step_collection<BinaryOp> subtract_1_step;
   CnC::step_collection<BinaryOp> subtract_2_step;
   CnC::step_collection<Compress_doIt> compress_doIt_funcA_step;
   CnC::step_collection<Compress_doIt> compress_doIt_funcB_step;
   CnC::step_collection<GaxpyOp> gaxpyOp_step;
   CnC::step_collection<Reconstruct_Prolog> reconstruct_prolog_step;
   CnC::step_collection<Reconstruct_doIt>  reconstruct_doIt_step;
   CnC::step_collection<Norm2> norm2_step;

  gaxpy_sub_test( int k, double thresh, int max_level,  int height, double (*funcA)(double), double (*funcB)(double))
  :
  CnCContext( k, thresh, max_level, height),
  projectA_item(*this),
  projectA_tag(*this), 
  printer_tag(*this), 

  projectB_item(*this), 
  projectB_tag(*this), 
  subtract_1_tag(*this), 
  subtract_2_tag(*this), 
  subtract_1_item(*this), 
  subtract_2_item(*this), 
  compressA_left_item(*this),
  compressA_right_item(*this),
  funcA_coeff_compressed_item(*this),
  compressB_left_item(*this),
  compressB_right_item(*this),
  funcB_coeff_compressed_item(*this),
  gaxpy_result_item(*this),
  reconstruct_result_item(*this),
  s_coeff_item(*this),
  compress_doIt_funcA_tag(*this),
  compress_doIt_funcB_tag(*this),
  gaxpyOP_tag(*this),
  reconstruct_prolog_tag(*this),
  reconstruct_doIt_tag(*this),
  norm2_tag(*this),

  /*----------------------------------------------------------------*/
  /* Declare projectA_step
  /*----------------------------------------------------------------*/
  projectA_step(
              *this, 
              "projectA_step", 
              Project(
                      funcA, 
                      std::vector<CnC::item_collection<std::pair< int, pair<int, int>>, Node> *>{},
                      std::vector<OutputTerminalType> {
                          OutputTerminalType(nullptr, std::vector<CnC::tag_collection<std::pair< int, pair<int, int>>> *> {&projectA_tag}),
                          OutputTerminalType(&projectA_item, std::vector<CnC::tag_collection<std::pair< int, pair<int, int>>> *> {&compress_doIt_funcA_tag})}
                          // OutputTerminalType(&projectA_item, std::vector<CnC::tag_collection<std::pair< int, pair<int, int>>> *> {&compress_doIt_funcA_tag, &subtract_1_tag})}
                      )
              ),

  /*----------------------------------------------------------------*/
  /* Declare projectB_step
  /*----------------------------------------------------------------*/
  projectB_step(
              *this, 
              "projectB_step", 
              Project(
                      funcB, 
                      std::vector<CnC::item_collection<std::pair< int, pair<int, int>>, Node> *>{},
                      std::vector<OutputTerminalType> {
                          OutputTerminalType(nullptr, std::vector<CnC::tag_collection<std::pair< int, pair<int, int>>> *> {&projectB_tag}),
                          // OutputTerminalType(&projectB_item, std::vector<CnC::tag_collection<std::pair< int, pair<int, int>>> *>{})}
                          OutputTerminalType(&projectB_item, std::vector<CnC::tag_collection<std::pair< int, pair<int, int>>> *>{&compress_doIt_funcB_tag})}
                      )
              ),


  /*----------------------------------------------------------------*/
  /* Declare subtract_1_step
  /*----------------------------------------------------------------*/
  subtract_1_step(
                *this, 
                "subtract_1_step", 
                BinaryOp(
                        &sub, 
                        &sub_scale_factor, 
                        std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *> {&projectA_item, &projectB_item},
                        std::vector<OutputTerminalType>{
                            OutputTerminalType(&projectA_item, std::vector<CnC::tag_collection<std::pair< int, pair<int, int>>> *> {&subtract_1_tag}),
                            OutputTerminalType(&subtract_1_item, std::vector<CnC::tag_collection<std::pair< int, pair<int, int>>> *> {}),
                            OutputTerminalType(&projectB_item, std::vector<CnC::tag_collection<std::pair< int, pair<int, int>>> *> {})}
                        )
                ),


  /*----------------------------------------------------------------*/
  /* Declare subtract_2_step
  /*----------------------------------------------------------------*/
  subtract_2_step(
                *this, 
                "subtract_2_step", 
                BinaryOp(
                        &sub, 
                        &sub_scale_factor, 
                        std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *> {&reconstruct_result_item, &subtract_1_item},
                        std::vector<OutputTerminalType>{
                            OutputTerminalType(&reconstruct_result_item, std::vector<CnC::tag_collection<std::pair< int, pair<int, int>>> *> {}),
                            OutputTerminalType(&subtract_2_item, std::vector<CnC::tag_collection<std::pair< int, pair<int, int>>> *> {&printer_tag}),
                            OutputTerminalType(&subtract_1_item, std::vector<CnC::tag_collection<std::pair< int, pair<int, int>>> *> {&subtract_2_tag})}
                        )
                ),

  /*----------------------------------------------------------------*/
  /* Declare compress_doIt_funcA_step
  /*----------------------------------------------------------------*/
  compress_doIt_funcA_step ( 
                            *this, 
                            "compress_doIt_funcA_step", 
                            Compress_doIt( 
                                std::vector<CnC::item_collection<std::pair< int, pair<int, int>>, Node> *> {&compressA_left_item, &projectA_item, &compressA_right_item}, 
                                std::vector<OutputTerminalType>{
                                OutputTerminalType(&compressA_left_item, std::vector<CnC::tag_collection<std::pair< int, pair<int, int>>> *> {}),
                                OutputTerminalType(&funcA_coeff_compressed_item, std::vector<CnC::tag_collection<std::pair< int, pair<int, int>>> *> {&printer_tag}),
                                OutputTerminalType(&compressA_right_item, std::vector<CnC::tag_collection<std::pair< int, pair<int, int>>> *> {})})
                          ),


  /*----------------------------------------------------------------*/
  /* Declare compress_doIt_funcB_step
  /*----------------------------------------------------------------*/
  compress_doIt_funcB_step ( 
                            *this, 
                            "compress_doIt_funcB_step", 
                            Compress_doIt( 
                                std::vector<CnC::item_collection<std::pair< int, pair<int, int>>, Node> *> {&compressB_left_item, &projectB_item, &compressB_right_item}, 
                                std::vector<OutputTerminalType>{
                                OutputTerminalType(&compressB_left_item, std::vector<CnC::tag_collection<std::pair< int, pair<int, int>>> *> {}),
                                OutputTerminalType(&funcB_coeff_compressed_item, std::vector<CnC::tag_collection<std::pair< int, pair<int, int>>> *> {}),
                                OutputTerminalType(&compressB_right_item, std::vector<CnC::tag_collection<std::pair< int, pair<int, int>>> *> {})})
                          ),


  /*----------------------------------------------------------------*/
  /* Declare gaxpyOp_step
  /*----------------------------------------------------------------*/
  gaxpyOp_step ( 
              *this, 
              "gaxpyOp_step", 
              GaxpyOp( 
                  1.0,
                  -1.0,
                  std::vector<CnC::item_collection<std::pair< int, pair<int, int>>, Node> *> {&funcA_coeff_compressed_item, &funcB_coeff_compressed_item}, 
                  std::vector<OutputTerminalType>{
                  OutputTerminalType(&funcA_coeff_compressed_item, std::vector<CnC::tag_collection<std::pair< int, pair<int, int>>> *> {&gaxpyOP_tag}),
                  OutputTerminalType(&gaxpy_result_item, std::vector<CnC::tag_collection<std::pair< int, pair<int, int>>> *> {&reconstruct_prolog_tag}),
                  OutputTerminalType(&funcB_coeff_compressed_item, std::vector<CnC::tag_collection<std::pair< int, pair<int, int>>> *> {})})
            ),

  /*----------------------------------------------------------------*/
  /* Declare reconstruct_prolog_step
  /*----------------------------------------------------------------*/
  reconstruct_prolog_step ( 
                        *this, 
                        "reconstruct_prolog_step", 
                        Reconstruct_Prolog( 
                            std::vector<CnC::item_collection<std::pair< int, pair<int, int>>, Node> *> {&gaxpy_result_item}, 
                            std::vector<OutputTerminalType>{
                            OutputTerminalType(&s_coeff_item, std::vector<CnC::tag_collection<std::pair< int, pair<int, int>>> *> {&reconstruct_doIt_tag})})
                      ),


  /*----------------------------------------------------------------*/
  /* Declare reconstruct_doIt_step
  /*----------------------------------------------------------------*/
  reconstruct_doIt_step( 
                      *this, 
                      "reconstruct_doIt_step", 
                      Reconstruct_doIt( 
                          std::vector<CnC::item_collection<std::pair< int, pair<int, int>>, Node> *> {&s_coeff_item, &gaxpy_result_item}, 
                          std::vector<OutputTerminalType>{
                            OutputTerminalType(&s_coeff_item, std::vector<CnC::tag_collection<std::pair< int, pair<int, int>>> *> {&reconstruct_doIt_tag}),
                            OutputTerminalType(&reconstruct_result_item, std::vector<CnC::tag_collection<std::pair< int, pair<int, int>>> *> {&subtract_2_tag})})
                     ),


  /*----------------------------------------------------------------*/
  /* Declare printer_step
  /*----------------------------------------------------------------*/
  printer_step(
            *this, 
            "printer_step", 
            Printer( 
              std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *>{&funcA_coeff_compressed_item}, 
              std::vector<OutputTerminalType>{})
            ),

  /*----------------------------------------------------------------*/
  /* Declare subtract_2_step
  /*----------------------------------------------------------------*/
  norm2_step( 
            *this, 
            "norm2_step", 
            Norm2( 
                std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *> {&subtract_2_item}, 
                std::vector<OutputTerminalType>{
                  OutputTerminalType(nullptr, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&norm2_tag})})
           )

    {

      /*----------------------------------------------------------------*/
      /* Tag Prescription
      /*----------------------------------------------------------------*/
      projectA_tag.prescribes(projectA_step, *this);
      projectB_tag.prescribes(projectB_step, *this);
      printer_tag.prescribes(printer_step, *this);
      subtract_1_tag.prescribes(subtract_1_step, *this);
      subtract_2_tag.prescribes(subtract_2_step, *this);
      compress_doIt_funcA_tag.prescribes( compress_doIt_funcA_step, *this);
      compress_doIt_funcB_tag.prescribes( compress_doIt_funcB_step, *this);
      gaxpyOP_tag.prescribes( gaxpyOp_step, *this);
      reconstruct_prolog_tag.prescribes( reconstruct_prolog_step, *this);
      reconstruct_doIt_tag.prescribes( reconstruct_doIt_step, *this);
      norm2_tag.prescribes( norm2_step, *this);

      /*----------------------------------------------------------------*/
      /* Steps Produce Consume
      /*----------------------------------------------------------------*/
      projectA_step.produces(projectA_item);
      projectB_step.produces(projectB_item);

      compress_doIt_funcA_step.consumes( projectA_item );
      compress_doIt_funcA_step.consumes( compressA_left_item );
      compress_doIt_funcA_step.consumes( compressA_right_item );
      compress_doIt_funcA_step.produces( funcA_coeff_compressed_item);
      compress_doIt_funcA_step.produces( compressA_left_item);
      compress_doIt_funcA_step.produces( compressA_right_item);

      compress_doIt_funcB_step.consumes( projectB_item );
      compress_doIt_funcB_step.consumes( compressB_left_item );
      compress_doIt_funcB_step.consumes( compressB_right_item );
      compress_doIt_funcB_step.produces( funcB_coeff_compressed_item);
      compress_doIt_funcB_step.produces( compressB_left_item );
      compress_doIt_funcB_step.produces( compressB_right_item );

      gaxpyOp_step.consumes( funcA_coeff_compressed_item);
      gaxpyOp_step.consumes( funcB_coeff_compressed_item);
      gaxpyOp_step.produces( funcA_coeff_compressed_item);
      gaxpyOp_step.produces( funcB_coeff_compressed_item);
      gaxpyOp_step.produces( gaxpy_result_item);

      reconstruct_prolog_step.consumes( gaxpy_result_item );
      reconstruct_prolog_step.produces( s_coeff_item);

      reconstruct_doIt_step.consumes( s_coeff_item);
      reconstruct_doIt_step.consumes( gaxpy_result_item);
      reconstruct_doIt_step.produces( s_coeff_item);
      reconstruct_doIt_step.produces( reconstruct_result_item);

      subtract_1_step.consumes(projectA_item);
      subtract_1_step.consumes(projectB_item);
      subtract_1_step.produces(projectA_item);
      subtract_1_step.produces(projectB_item);
      subtract_1_step.produces(subtract_1_item);

      subtract_2_step.consumes(subtract_1_item);
      subtract_2_step.consumes(reconstruct_result_item);
      subtract_2_step.produces(subtract_1_item);
      subtract_2_step.produces(reconstruct_result_item);
      subtract_2_step.produces(subtract_2_item);

      printer_step.consumes(funcA_coeff_compressed_item);
      norm2_step.consumes( subtract_2_item);
  }

  // void check_result( CnC::item_collection<std::pair< int, pair<int, int>>, Node> result_h, int n, int l, int h){
  //   Node node_1, node_h;
  //   int new_h;

  //   if( h == 1)
  //     new_h = global_co_height;
  //   else
  //     h = h-1;

  //   result_1.get( make_pair(1, make_pair(n,l)), node_1);
  //   result_h.get( make_pair(h, make_pair(n,l)), node_h);

  //   cout << "Without coarsening: " << node_1.toString() << "  With coarsening " << node_h.toString() << "\n";

  //   if( node_1.has_children && node_h.has_children ){
  //       check_result(result_1, result_h, n+1, 2*l, new_h);
  //       check_result(result_1, result_h, n+1, 2*l+1, new_h);
  //   }
  //   else if( !node_1.has_children && !node_h.has_children ){
  //     return;
  //   }
  //   else{
  //     cout << "Error!\n";
  //   }
  // }

};
