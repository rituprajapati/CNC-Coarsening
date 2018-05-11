const double pi = 3.1415926535897;


/*****************************************************************/
/* gaussian with square normalized to 1              */
/*****************************************************************/

double test1 ( double x ){
  
  double a = 500.0;
  double tmp = 2 * a / pi;
  double tmp2 = -a * pow(( x-0.5 ), 2);
  return pow( tmp , 0.25 ) * exp( tmp2 );
  // return x;

}

/*****************************************************************/
/* superposition of multiple gaussians               */
/*****************************************************************/
double test2 (  double x ) { 

  return ( test1( x-0.3) + test1 ( x ) + test1 ( x+0.3 ) );
}

/*****************************************************************/
/* a more interesting (singular and oscillating) function      */
/* ... note it is never computed exactly at the singularity    */
/* except by the test code below.                  */
/*****************************************************************/

double test3 ( double x ){
  double a = 100.0 * pi;

  if ( x == 0.5 ){
    return 0.0;
  }
  else{
    return cos ( a * ( x-0.5) ) / sqrt ( abs ( x-0.5 ));
  }
}

/*****************************************************************/
/* derivative of test1                       */
/*****************************************************************/

double dtest1 ( double x ) {

  double a = 500.0;
  double tmp1 =  -2.0 * a * (x-0.5);
  double tmp2 =  pow( 2 * a /pi, 0.25 );
  double tmp3 =  exp( -a * pow( (x-0.5), 2 ));
  return ( tmp1 * tmp2 * tmp3);

  // return 1;
}


/*****************************************************************/
/* derivative of test2                       */
/*****************************************************************/

double dtest2 ( double x ){

  return  ( dtest1( x-0.3 ) + dtest1( x )  + dtest1( x+0.3) );
}


/*****************************************************************/
/* derivative of test3                       */
/*****************************************************************/

double dtest3 ( double x ){

  double a = 100.0 * pi;

  if ( x == 0.5 ){
    return 0.0;
  }
  else{

    double s = 1.0;
    if ( x < 0.5 ){
      s = -1.0;
    }
    double tmp1 = -a * sin( a * (x-0.5) );
    double tmp2 = sqrt( abs(x-0.5) );
    double tmp3 = s * 0.5 * cos( a * (x-0.5) );
    double tmp4 = pow( abs( x-0.5 ), 1.5 );

    return ( tmp1 / tmp2 - tmp3 / tmp4 );
  }
} 

double ( *test[3] ) ( double x ) = { test1, test2, test3 };
double ( *dtest[3] ) ( double x ) = { dtest1, dtest2, dtest3 };