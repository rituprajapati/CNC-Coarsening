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
// http://stackoverflow.com/questions/22387586/measuring-execution-time-of-a-function-in-c
#include <chrono>
using namespace std;
// http://stackoverflow.com/questions/22387586/measuring-execution-time-of-a-function-in-c
using namespace std::chrono;


/**********************************************************************/
/*  Forward declaration
/**********************************************************************/

struct CnCContext;
Vector sub(const Vector &v1, const Vector &v2);
double sub_scale_factor(double n);
void init_twoscale(int );



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

   void put(int h, K key, V value) const {
      item_collection->put(key, value);
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
int Project::execute(const std::pair< int, std::pair<int, int>> &node, CnCContext &context) const {

  int n, l, h;
  h = node.first;
  n = node.second.first;
  l = node.second.second; 

  Vector s0 = sValue(n + 1, 2 * l, context);
  Vector s1 = sValue(n + 1, 2 * l + 1, context);

  int k = context.k;
  Vector s(s0 | s1); // concatenation of s0 and s1
  Vector d(s * (*(context.hgT)));

  //Case 2
   if (d.normf(k, 2 * k) < context.thresh || n >= context.max_level - 1) {

    if( h == context.max_height )
      output_terminals[1].put(node, Node(n, l, k, Vector(), Vector(), true));
    else
      output_terminals[1].put(h, node, Node(n, l, k, Vector(), Vector(), true));
  
    //If height is 1, then the children needs to trigger new threads
    if( h == 1 ){
      output_terminals[1].put(std::make_pair(context.max_height, std::make_pair(n + 1, 2 * l)), Node(n + 1, 2 * l, k, s0, Vector(), false));
      output_terminals[1].put(std::make_pair(context.max_height, std::make_pair(n + 1, 2 * l +1)), Node(n + 1, 2 * l + 1, k, s1, Vector(), false));
    }
    //Else only put the items
    else{
      output_terminals[1].put(h-1, std::make_pair(h-1, std::make_pair(n + 1, 2 * l)), Node(n + 1, 2 * l, k, s0, Vector(), false));
      output_terminals[1].put(h-1, std::make_pair(h-1, std::make_pair(n + 1, 2 * l +1)), Node(n + 1, 2 * l + 1, k, s1, Vector(), false));
    }
  }

  //Case 3
   else {
    
     if( h == context.max_height )
      output_terminals[1].put(node, Node(n, l, k, Vector(), Vector(), true));
    else
      output_terminals[1].put(h, node, Node(n, l, k, Vector(), Vector(), true));


    if( h == 1){
      output_terminals[0].put(std::make_pair(context.max_height, std::make_pair(n + 1, 2 * l)));
      output_terminals[0].put(std::make_pair(context.max_height, std::make_pair(n + 1, 2 * l+1)));  
    }
    else{
      execute( std::make_pair(h-1, std::make_pair(n + 1, 2 * l)), context);
      execute( std::make_pair(h-1, std::make_pair(n + 1, 2 * l+1)), context);
    }
  }
  return CnC::CNC_Success;
}


/**********************************************************************/
/* Printer::execute
/**********************************************************************/
int Printer::execute(const std::pair<int, pair<int, int>> &node, CnCContext &context) const {

  int n, l, h;
  h = node.first;
  n = node.second.first;
  l = node.second.second; 

  Node nodeInfo;
  input_terminals[0]->get(node, nodeInfo);

  std::cout << "Printer:: Node with info: (Key:"
            << h << " (" 
            << n << ", " 
            << l << "), " 
            << nodeInfo.toString() 
            << ")" 
            << std::endl;


  if( nodeInfo.has_children && h > 1){          
    execute( std::make_pair(h-1, std::make_pair(n + 1, 2 * l)), context);
    execute( std::make_pair(h-1, std::make_pair(n + 1, 2 * l+1)), context);
  }
  return CnC::CNC_Success;
}



/**********************************************************************/
/* Differentiation test
/**********************************************************************/
struct Project_test: CnCContext{

  /*----------------------------------------------------------------*/
  /* Item Collections
  /*----------------------------------------------------------------*/
  CnC::item_collection<std::pair<int, pair<int, int>>, Node> project_item;
  CnC::tag_collection<std::pair<int, pair<int, int>>> project_tag;
  CnC::tag_collection<std::pair<int, pair<int, int>>> printer_tag;

  using OutputTerminalType = OutputTerminal<std::pair<int, pair<int, int>>, Node>;

  CnC::step_collection<Project> project_step;
  CnC::step_collection<Printer> printer_step;

  Project_test( int k, double thresh, int max_level,  int height, double (*funcA)(double))
  :
  CnCContext( k, thresh, max_level, height),
  project_item(*this),
  project_tag(*this),
  printer_tag(*this),

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
                            OutputTerminalType(&project_item, std::vector<CnC::tag_collection<std::pair<int, pair<int, int>>> *> {&printer_tag})}
                        )
                ),

  /*----------------------------------------------------------------*/
  /* Declare printer_step
  /*----------------------------------------------------------------*/
  printer_step(
              *this, 
              "printer_step", 
              Printer( 
                std::vector<CnC::item_collection<std::pair<int, pair<int, int>>, Node> *>{&project_item}, 
                std::vector<OutputTerminalType>{})
              )
  {

    /*----------------------------------------------------------------*/
    /* Tag Prescription
    /*----------------------------------------------------------------*/
    project_tag.prescribes(project_step, *this);
    printer_tag.prescribes(printer_step, *this);

    /*----------------------------------------------------------------*/
    /* Steps Produce Consume
    /*----------------------------------------------------------------*/

    project_step.produces(project_item);

    printer_step.consumes(project_item);

  }

};


int main(int argc, char *argv[]) {
   int k = 5;
   int max_level;
   double thresh;
   int height;

   std::istringstream iss1( argv[1] );
   if( ! (iss1 >> max_level)){
    cout << "Unable to assign max_levels\n";
    return -1;
   }

   std::istringstream iss2( argv[2] );
   if( ! (iss2 >> thresh)){
    cout << "Unable to assign thresh\n";
    return -1;
   }

   std::istringstream iss3( argv[3] );
   if( ! (iss3 >> height)){
    cout << "Unable to assign height\n";
    return -1;
   }

  

   high_resolution_clock::time_point t1 = high_resolution_clock::now();

   Project_test obj(k, thresh, max_level, height, test[0] );
   obj.project_tag.put(std::make_pair(height, make_pair(0, 0)));
   obj.wait();

   high_resolution_clock::time_point t2 = high_resolution_clock::now();
   auto duration = duration_cast<microseconds>( t2 - t1 ).count();
   std::cout << "Time taken during computation: " << duration/1000000.0 << std::endl;
   return 0;
}

