
#include "oned_queue.h"


void Multiplication( int k, int max_level, double thresh, int height ){

  int npt = 20;

  vector< double > range;
  for (double i = 0.0; i < (double)npt+1.0; i++)
    range.push_back(i);

  /*****************************************************************/
  /* Multiplication test 
  /*****************************************************************/
  for( int j = 0; j < 3; j++ ){
    
      cout << "\n";

      multiplication_test mul_test_obj( k, thresh, max_level, height, test[0], test[j]);
      mul_test_obj.projectA_tag.put(make_pair( 0, std::make_pair(0, 0) ) );
      mul_test_obj.projectB_tag.put(make_pair( 0, std::make_pair(0, 0) )); 
      mul_test_obj.wait();

      mul_test_obj.norm2_f1_tag.put(make_pair( 0, std::make_pair(0, 0) ));
      mul_test_obj.wait();
      cout << "norm of f1 is  " << sqrt( mul_test_obj.norm2_result[0]) << "\n";
      mul_test_obj.norm2_result[0] = 0.0;


      mul_test_obj.norm2_f2_tag.put(make_pair( 0, std::make_pair(0, 0) ));
      mul_test_obj.wait();
      cout << "norm of f2 is  " << sqrt(mul_test_obj.norm2_result[0]) << "\n";
      mul_test_obj.norm2_result[0] = 0.0;


      mul_test_obj.norm2_mul_tag.put(make_pair( 0, std::make_pair(0, 0) ));
      mul_test_obj.wait();
      cout << "norm of Multiplication is  " << sqrt(mul_test_obj.norm2_result[0]) << "\n";
      mul_test_obj.norm2_result[0] = 0.0;

  }

}


void Addition( int k, int max_level, double thresh, int height ){

   int npt = 20;

    vector< double > range;
    for (double i = 0.0; i < (double)npt+1.0; i++)
      range.push_back(i);

  /*****************************************************************/
  /* Addition test which in turn tests gaxpy             */
  /*****************************************************************/
  for( int j = 0; j < 3; j++ ){
    
      cout << "\n";

      addition_test add_obj(k, thresh, max_level, height, test[0], test[j]);
      add_obj.projectA_tag.put(make_pair( 0, std::make_pair(0, 0) ));
      add_obj.projectB_tag.put(make_pair( 0, std::make_pair(0, 0) ));   
      add_obj.wait();

      add_obj.norm2_f1_tag.put(make_pair( 0, std::make_pair(0, 0) ));
      add_obj.wait();
      cout << "norm of f1 is  " << sqrt(add_obj.norm2_result[0]) << "\n";
      add_obj.norm2_result[0] = 0.0;


      add_obj.norm2_f2_tag.put(make_pair( 0, std::make_pair(0, 0) ));
      add_obj.wait();
      cout << "norm of f2 is  " << sqrt(add_obj.norm2_result[0]) << "\n";
      add_obj.norm2_result[0] = 0.0;


      add_obj.norm2_add_tag.put(make_pair( 0, std::make_pair(0, 0) ));
      add_obj.wait();
      cout << "norm of addition is  " << sqrt(add_obj.norm2_result[0]) << "\n";
      add_obj.norm2_result[0] = 0.0;
  }

}



void Differentiation( int k, int max_level, double thresh, int height ){

  int npt = 20;
  vector< double > range;
  for (double i = 0.0; i < (double)npt+1.0; i++)
    range.push_back(i);


  for( int j = 0; j < 3; j++ ){

    diff_test diff_test_obj( k, thresh, max_level, height, test[j]);

    //Initialize the 0th and 1st level with 
    diff_test_obj.project_tag.put( make_pair( 0, std::make_pair(0, 0) ) );
    diff_test_obj.wait();

    diff_test_obj.norm2_tag.put(make_pair( 0, std::make_pair(0, 0) ));
    diff_test_obj.wait();

    cout << "\nNorm of function is " << sqrt(diff_test_obj.norm2_result[0]) << "\n";
  }

}



void Gaxpy_Sub( int k, int max_level, double thresh, int height ){

  int npt = 20;
  vector< double > range;
  for (double i = 0.0; i < (double)npt+1.0; i++)
    range.push_back(i);

  /*****************************************************************/
  /* Multiplication test 
  /*****************************************************************/
  for( int j = 0; j < 3; j++ ){
    
      gaxpy_sub_test test_obj( k, thresh, max_level, height, test[0], test[j]);
      test_obj.projectA_tag.put(make_pair( 0, std::make_pair(0, 0)));
      test_obj.projectB_tag.put(make_pair( 0, std::make_pair(0, 0))); 
      test_obj.wait();

      test_obj.norm2_tag.put(make_pair( 0, std::make_pair(0, 0)));
      test_obj.wait();
      cout << "norm of (test[0]-test[" << j <<  "])- (test[0]-test[" <<j << "]) is  " << sqrt( test_obj.norm2_result[0]) << "\n";
  }

}



int main(int argc, char *argv[]) {
   int k = 5;
   int max_level;
   double thresh;
   int num_threads, height;

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
    cout << "Unable to assig height\n";
    return -1;
   }

   high_resolution_clock::time_point t1 = high_resolution_clock::now();

   /*--------------------------------*/
   /*    DIFFERENTIATION TEST        */
   /*--------------------------------*/
   Differentiation( k, max_level, thresh, height );

   /*--------------------------------*/
   /*    ADDITION TEST               */
   /*--------------------------------*/
   Addition(k, max_level, thresh, height);


   /*--------------------------------*/
   /*    MULTIPLICATION TEST         */
   /*--------------------------------*/
   Multiplication( k, max_level, thresh, height );

   /*--------------------------------*/
   /*    MULTIPLICATION TEST         */
   /*--------------------------------*/
   Gaxpy_Sub(k, max_level, thresh, height);

   high_resolution_clock::time_point t2 = high_resolution_clock::now();
   auto duration = duration_cast<microseconds>( t2 - t1 ).count();
   std::cout << "Time taken during computation: " << duration/1000000.0 << std::endl;
   return 0;
}
