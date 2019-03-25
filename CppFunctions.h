#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;


// Function to print vector in cpp to check correctness 
void PrinVect(NumericVector x);

// Function to replicate number equivalent to rep in Rcpp exist for NumericVector 
NumericVector rep_N(long double x, int n);

// Function max of two long double 
long double Sup(long double a,long double b);

// Function min of two long double 
long double Inf(long double a,long double b);

  
// Function returning the infection rate of cattle
long double b_cattle( NumericVector x);

// Function returning the infection rate for people
long double b_people(NumericVector x);

//Function returning the proportion infection of mosquitoes
long double b_mos(NumericVector x);

//Function to perform the calculus of the rates for people, cattle and mosquitoes
NumericVector Rates_Updates(NumericVector A, NumericVector subparam);

//Function combining vectors in a list to return one vector after concatenation of the others
NumericVector combine(const List& list);

// transhumance flood shearing 
NumericVector transhumance(NumericVector t,NumericVector p, double val);
  
// return sinusoidal function of t for dry years
NumericVector SinFun(NumericVector t, NumericVector p);

//Function with the ODEs to be solved 
List ODE(NumericVector t, NumericVector state, NumericVector param);