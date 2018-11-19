#include "oned_queue.h"


int main(int argc, char *argv[]) {
   int k = 5;
   int max_level;
   double thresh;
   int height; //Rename to tile height

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
   // Project_test obj(k, thresh, max_level, height, test[0], test[1] );
   // obj.projectA_tag.put(std::make_pair(0, make_pair(0, 0)));
   // obj.projectB_tag.put(std::make_pair(0, make_pair(0, 0)));
   // obj.wait();

   diff_test diff_test_obj( k, thresh, max_level, height, test[0]);
   diff_test_obj.project_tag.put( std::make_pair(0, make_pair(0, 0)) );
   diff_test_obj.wait();

   high_resolution_clock::time_point t2 = high_resolution_clock::now();
   auto duration = duration_cast<microseconds>( t2 - t1 ).count();
   std::cout << "Time taken during computation: " << duration/1000000.0 << std::endl;
   return 0;
}

