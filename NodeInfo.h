#ifndef _NODE_INFO
#define _NODE_INFO

#include "Vector.h"
#include "Matrix.h"
#include <math.h>

/*#include <iostream>*/



class NodeInfo {
public:
   int n;
   int l;

   int which;	

   NodeInfo(int nInput, int lInput, int whichInput): n(nInput), l(lInput), which(whichInput) {}
   NodeInfo(const NodeInfo &node): n(node.n), l(node.l), which(node.which) {}
 
};

#endif
