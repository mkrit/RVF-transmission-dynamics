#include <Rcpp.h>
#include <math.h>
#include <iostream>
#include <limits>

using namespace std;
using namespace Rcpp;


/* Function to print vector in cpp to check correctness */
void PrinVect(NumericVector x){
  int i=0;
  for (i = 0; i < x.size(); i++ ) {
    Rcpp::Rcout << x(i) << std::endl; 
  }
}

/* Function to replicate number equivalent to rep in Rcpp exist for NumericVector */
NumericVector rep_N(long double x, int n){
  NumericVector v(rep(NumericVector::create(x),n));
  return v;
}

//[[Rcpp::export]]
/* Function max of two long double */
long double Sup(long double a,long double b)
{
  return double (a<=b ? b:a);
}

/* Function min of two long double */
long double Inf(long double a,long double b)
{
  return (a<=b ? a:b);
}


//[[Rcpp::export]]
// Function returning the infection rate of cattle
long double b_cattle( NumericVector x){
  long double a = x[0], b = x[1], c = x[2], d = x[3], e = x[4], f = x[5], g = x[6], h = x[7],  i = x[8], j = x[9], k = x[10], l = x[11], m = x[12], n = x[13];
  long double x1=-std::log (std::pow((1.0-b),(c*d/e)) * std::pow((1.0-f),(g*h/e))* std::pow((1.0-i),(j*k/e))* std::pow((1.0-l),(m*n/e)));
  x1=(std::isinf(x1) ? a : x1);
  b=(e==0.0 ? 0.0 : Inf(a,x1));
  
  //std::cout<<"b cattle "<<b<<" "<<x1<<" "<<a<<std::endl;
  return b;
}

//[[Rcpp::export]]
// Function returning the infection rate for people
long double b_people(NumericVector x) {
  long double a=x[0], b=x[1],  c=x[2] , d=x[3] , e=x[4], f=x[5], g=x[6], h=x[7], i=x[8], j=x[9], k=x[10], l=x[11];
  long double m=x[12],n=x[13], o=x[14], p=x[15], q=x[16], r=x[17];
  long double aa=0.0, bb=0.0;
  
  aa=(e==0.0 ? 0.0 : 1.0 - std::pow((1.0 - b),(c*d/e)) * std::pow((1.0-f),(g*h/e)) * std::pow((1.0-i),(j*k/e)) * std::pow((1.0-l),(m*n/e)));
  bb=(q+r==0.0 ? 0.0 : 1.0 - std::pow((1.0 - o),(p*q/(q+r))));
  long double x1 = -std::log (1.0 - (aa + bb - aa*bb));
  x1=(std::isinf(x1) ? a : x1);
  b = Inf(a,x1);
  return b;
}


//[[Rcpp::export]]
//Function returning the proportion infection of mosquitoes
long double b_mos(NumericVector x) {

  long double a = x[0], b = x[1], c = x[2], d = x[3], e = x[4], f = x[5], g = x[6];
  long double x1=(c==0.0 ? 0.0 : a*b/c)+(e+f==0.0 ? 0.0 : d*e/(e+f) )+ g;//ici b_wl is outside ifelse here
  return (x1<=1.0 ? x1 : 1.0);
  
}

// [[Rcpp::export]]
//Function to perform the calculus of the rates for people, cattle and mosquitoes
NumericVector Rates_Updates(NumericVector A, NumericVector subparam1, NumericVector subparam2){
  long double AS=A[0], AI=A[1], BS=A[2], BI=A[3], CS=A[4], CI=A[5], DS=A[6], DI=A[7], MI=A[8], NM=A[9], HI=A[10],NH=A[11];
  long double AA = AS + AI, BA = BS + BI, CA = CS + CI, DA = DS + DI;
  long double O_ = subparam1[0]*(AA)*subparam1[1] + subparam1[2]*(BA)*subparam1[3]+ subparam1[4]*(CA)*subparam1[5]+ subparam1[6]*(DA)*subparam1[7], chi_ = subparam1[8]*NH, temp=chi_/Sup(O_,0.000001), s_ = Inf(1.0,temp);
  long double o_ah = subparam1[0] * subparam1[1] * s_, o_bh = subparam1[2] * subparam1[3] * s_, o_ch = subparam1[4] * subparam1[5] * s_, o_dh =subparam1[6]*subparam1[7]*s_;
  long double O_m= subparam1[9]*(AA)*subparam1[1] + subparam1[10]*(BA)*subparam1[3]+subparam1[11]*(CA)*subparam1[5] + subparam1[12]*(DA)*subparam1[7], chi_m = subparam1[14]*subparam1[13]*(NM+MI);
  temp=chi_m/Sup(O_m,0.000001);
  long double  s_m = Inf(1.0,temp), o_am = subparam1[9] * subparam1[1] * s_m;
  long double o_bm = subparam1[10] * subparam1[3] * s_m, o_cm=subparam1[11] * subparam1[5] * s_m, o_dm=subparam1[12] * subparam1[7] * s_m;
  // o_a1 o_aa2 o_aa3
  long double o_a = o_ah + o_am, o_b = o_bh +o_bm, o_c= o_ch + o_cm, o_d= o_dh +o_dm; 
  // o_a1_2 o_ab1_2 o_c1_2 o_d1_2
  long double  o_aa = subparam1[15]/Sup(0.000001,AA),o_ba = subparam1[15]/Sup(0.000001,BA), o_ca=subparam1[15]/Sup(0.000001,CA), o_da = subparam1[15]/Sup(0.000001,DA); 

  NumericVector p=NumericVector::create(subparam1[16], HI, NH, subparam1[17], MI, NM, subparam2[15]);
  long double b_a = b_mos(p);
  p=NumericVector::create(subparam1[18], HI, NH, subparam1[19], MI, NM, subparam2[15]);
  long double b_b = b_mos(p);
  // b_a1 b_c1 b_d1 b_h1 b_m1
  p=NumericVector::create(subparam2[0], HI, NH, subparam2[1], MI, NM, subparam2[15]);
  long double b_c =b_mos(p);
  p=NumericVector::create(subparam2[2], HI, NH, subparam2[3], MI, NM, subparam2[15]);
  long double b_d =b_mos(p);
  //b_pe1 b_pe2 b_pe3
  
  p=NumericVector::create(subparam2[4], subparam2[5], o_ah, AI, NH, subparam2[6], o_bh, BI, subparam2[7], o_ch, CI,subparam2[8],o_dh,DI,subparam2[9],subparam2[10],MI,NM);
  long double b_h =(b_people(p));
  // b_ca1 b_ca2 b_ca3
  p=NumericVector::create(subparam2[4], subparam2[11], o_am, AI, NM, subparam2[12], o_bm, BI, subparam2[13], o_cm,  CI,subparam2[14], o_dm, DI);
  long double b_m =( b_cattle(p));

  NumericVector L(14);
  L=NumericVector::create(o_a, o_b, o_c, o_d, o_aa, o_ba, o_ca, o_da, b_a, b_b, b_c, b_d, b_h, b_m);
  return(L);
}

// [[Rcpp::export]]
//Function combining vectors to return one vector after concatenation of the others
NumericVector combine(const List& list)
{
  std::size_t n = list.size();
  
  // Figure out the length of the output vector
  std::size_t total_length = 0;
  for (std::size_t i = 0; i < n; ++i)
    total_length += Rf_length(list[i]);
  
  // Allocate the vector
  NumericVector output = no_init(total_length);
  
  // Loop and fill
  std::size_t index = 0;
  for (std::size_t i = 0; i < n; ++i)
  {
    NumericVector el = list[i];
    std::copy(el.begin(), el.end(), output.begin() + index);
    
    // Update the index
    index += el.size();
  }
  
  return output;
  
}


// [[Rcpp::plugins(cpp11)]]
//[[Rcpp::export]]
// transhumance flood shearing 
NumericVector transhumance(NumericVector t,NumericVector p, double val, bool trans){
  NumericVector l_m(t.size());
  for(int it = 0; it != l_m.size(); ++it) {
    l_m[it]=((trans && fmod(t[it],360.0)>=(p[0]-0.4) && fmod(t[it],360.0)<=(p[0]+0.5))?p[1] : val);}
  return l_m;
  }

//[[Rcpp::export]]
// return sinusoidal function of t
NumericVector SinFun(NumericVector t, NumericVector p){
  NumericVector x(t.size());
  for(int i=0;i<t.size();i++) x[i]= (0.5 * cos(4*atan(1)*(t[i] - p[0])/p[1]) + 0.5);
  return(x);
}

//[[Rcpp::export]]
//Function with the ODEs to be solved 
List ODE(NumericVector t, NumericVector state, NumericVector param){
  NumericVector ts, t_a1; 
  long double b_wl = param[2];
  long double max_rate = param[0], flood_prop = param[1], O_alt = param[3], d1 = param[4], d2 = param[5], d3 = param[6], d4 = param[7], d5 = param[8], d6 = param[9], year = param[10], c0 = param[11], mmm = param[12], ds = param[13], nPeak = param[14];
  long double seasonHatch = param[15],g_h= param[16], m_h= param[17], x_h= param[18], a_h= param[19], d_h= param[20], p_mh00= param[21], f_mh1= param[22], f_mh2= param[23], f_mh3= param[24], h_h1= param[25], h_h2= param[26], h_h3= param[27];
  long double p_ha= param[28], p_hb= param[29], p_hc= param[30], p_hd= param[31], r_h= param[32], l_h12= param[33], l_h13= param[34], l_h21= param[35], l_h23= param[36], l_h31= param[37], l_h32= param[38];
  long double g_m_u= param[39], g_m_i= param[40], m_m= param[41], x_m= param[42], a_m= param[43], d_m= param[44], h_m= param[45], p_ma= param[46], p_mb= param[47], p_mc= param[48], p_md= param[49], r_m= param[50], k_m1= param[51], k_m2= param[52], k_m3= param[53];
  long double l_m13= param[54], l_m23= param[55], l_m31= param[56], l_m32= param[57],g_a= param[58], z_a= param[59], m_a= param[60], v_a= param[61], e_ah= param[62], e_am= param[63], p_ah= param[64], p_am= param[65], k_a1= param[66], k_a2= param[67], k_a3= param[68]; 
  long double m_aq1= param[69], m_aq2= param[70], m_aq3= param[71], m_ap1= param[72], m_ap2= param[73], m_ap3= param[74], l_a12= param[75], l_a13= param[76], l_a21= param[77], l_a23= param[78], l_a31= param[79], l_a32= param[80];
  long double g_b= param[81], z_b= param[82], m_b= param[83], t_b= param[84], v_b= param[85], e_bh= param[86], e_bm= param[87], p_bh= param[88], p_bm= param[89], k_b1= param[90], k_b2= param[91], k_b3= param[92];
  long double m_bq1= param[93], m_bq2= param[94], m_bq3= param[95], m_bp1= param[96], m_bp2= param[97], m_bp3= param[98], l_b12= param[99], l_b13= param[100], l_b21= param[101], l_b23= param[102], l_b31= param[103], l_b32= param[104];
  long double g_c= param[105], m_c= param[106], t_c= param[107], v_c= param[108], e_ch= param[109], e_cm= param[110], p_ch= param[111], p_cm= param[112], k_c1= param[113], k_c2= param[114], k_c3= param[115], m_cp1= param[116], m_cp2= param[117], m_cp3= param[118];
  long double l_c12= param[119], l_c13= param[120], l_c21= param[121], l_c23= param[122], l_c31= param[123], l_c32= param[124];
  long double g_d= param[125], m_d= param[126], t_d= param[127], v_d= param[128], e_dh= param[129], e_dm= param[130], p_dh= param[131], p_dm= param[132], k_d1= param[133], k_d2= param[134], k_d3= param[135], m_dp1= param[136], m_dp2= param[137], m_dp3= param[138];
  long double l_d12= param[139], l_d13= param[140], l_d21= param[141], l_d23= param[142], l_d31= param[143], l_d32= param[144];
  long double shearBeg = param[146], shearEnd = param[147], shearUp = param[148];
  bool shearing= param[145], wetDry= param[149], flood= param[150], elNino= param[151], transHum = param[152];
  long double t_a=param[153], l_m21Base=param[154],l_m12Base= param[155];
  NumericVector dry2=param[Range(156,(156+year))];
  long double HS1 = state[0], HE1 = state[1], HI1 = state[2], HR1 = state[3], HS2 = state[4], HE2 = state[5], HI2 = state[6], HR2 = state[7], HS3 = state[8], HE3 = state[9], HI3 = state[10], HR3 = state[11];
  long double MS1 = state[12], ME1 = state[13], MI1 = state[14], MR1 = state[15], MS2 =state[16], ME2 = state[17], MI2 = state[18], MR2 = state[19], MS3 = state[20], ME3 = state[21], MI3 = state[22], MR3 = state[23];
  long double AQ1 = state[24], AP1 = state[25], AS1 = state[26], AI1 = state[27], AQ2 = state[28], AP2 =state[29], AS2 =state[30], AI2 = state[31], AQ3 = state[32], AP3 = state[33], AS3 = state[34], AI3 = state[35];
  long double BQ1 = state[36], BP1 = state[37], BS1 = state[38], BI1 = state[39], BQ2 = state[40], BP2 = state[41], BS2 = state[42], BI2 = state[43], BQ3 = state[44], BP3 = state[45], BS3 = state[46], BI3 = state[47];
  long double CP1 = state[48], CS1 = state[49], CI1 = state[50], CP2 = state[51], CS2 = state[52], CI2 = state[53], CP3 = state[54], CS3 = state[55], CI3 = state[56];
  long double DP1 = state[57], DS1 = state[58], DI1 = state[59], DP2 = state[60], DS2 = state[61], DI2 = state[62], DP3 = state[63], DS3 =state[64], DI3 = state[65];
    
  // population sizes per compartment and per species
  //----------------  human
  long double NH1 = HS1 + HE1 + HI1 + HR1;
  long double NH2 = HS2 + HE2 + HI2 + HR2;
  long double NH3 = HS3 + HE3 + HI3 + HR3;
  //----------------  animal host
  long double NM1 = MS1 + ME1 + MR1;
  long double NM2 = MS2 + ME2 + MR2;
  long double NM3 = MS3 + ME3 + MR3;
  //----------------  vector A
  long double NA1 = AQ1 + AP1 + AS1 + AI1;
  long double NA2 = AQ2 + AP2 + AS2 + AI2;
  long double NA3 = AQ3 + AP3 + AS3 + AI3;
  //----------------  vector B
  long double  NB1 = BQ1 + BP1 + BS1 + BI1;
  long double NB2 = BQ2 + BP2 + BS2 + BI2;
  long double NB3 = BQ3 + BP3 + BS3 + BI3;
  //----------------  vector C
  long double NC1 = CP1 + CS1 + CI1;
  long double NC2 = CP2 + CS2 + CI2;
  long double NC3 = CP3 + CS3 + CI3;
  //----------------  vector D
  long double ND1 = DP1 + DS1 + DI1;
  long double ND2 = DP2 + DS2 + DI2;
  long double ND3 = DP3 + DS3 + DI3;
   
  // increased susceptibility of sheep to mosquito bites because of shearing
  long double shear = 1.0; 
  long double sB=(shearBeg-0.4), sE=(shearEnd+0.5);
  bool a=(shearing && fmod(t[0],360.0)>=sB && fmod(t[0],360.0)<=sE);
  shear = (a? shearUp :  1.0);

  // transhumance
  NumericVector x=NumericVector::create(d1,max_rate);
  NumericVector l_m21 = transhumance(t,x,l_m21Base,transHum);
  x[0]= d2;
  NumericVector l_m12 = transhumance(t,x,l_m12Base,transHum);
 // wet and dry years
  int year_i= t[0]/360;
  long double dry =((wetDry && dry2[year_i] < c0)? mmm: 1.0);
  x[0] = ds, x[1]=nPeak;
  NumericVector test0=SinFun(t,x);
  ts =  dry *(seasonHatch==1? test0: rep_N(1.0,t.size()));
  // seasonal and El Nino hatching of Ae. mcintoshi eggs 
  // annual flooding in march, el nino in december
  NumericVector flooded(t.size());
  NumericVector elNinoed(t.size());
  for(int i=0;i<t.size();i++){
    flooded[i]=(flood && fmod(t[i],360.0)>=(d3-0.4) && fmod(t[i],360.0)<=(d4+0.5)? 1.0 : 0.0);
    elNinoed[i]=(elNino && fmod(t[i],3600.0)>=(d5-0.4) && fmod(t[i],3600.0)<=(d6+0.5)? 1.0 : 0.0);
  }
  t_a1 = flood_prop*max_rate*flooded+ max_rate*elNinoed;
  NumericVector p1=NumericVector::create(e_ah,v_a,e_bh,v_b,e_ch,v_c,e_dh, v_d, h_h1,e_am,e_bm,e_cm, e_dm, shear,h_m, O_alt, p_ha,p_ma, p_hb, p_mb); //19 indice 20 elements
  NumericVector p2=NumericVector::create(p_hc, p_mc, p_hd, p_md, max_rate, p_ah, p_bh, p_ch, p_dh, p_mh00, f_mh1, p_am, p_bm, p_cm, p_dm, 0.0);
  NumericVector X = NumericVector::create(AS1, AI1, BS1,  BI1, CS1, CI1, DS1,  DI1, MI1, NM1, HI1, NH1);
  // biting rates, mortality rates, infection rates floodplain
  NumericVector X1 = Rates_Updates(X, p1, p2);
  long double o_a1=X1[0],o_b1=X1[1],o_c1=X1[2], o_d1=X1[3], o_a1_2=X1[4], o_b1_2=X1[5], o_c1_2=X1[6], o_d1_2=X1[7], b_a1=X1[8],b_b1=X1[9], b_c1=X1[10], b_d1=X1[11], b_h1=X1[12], b_m1=X1[13];
  // Updates h_h1<- h_h2  and f_mh1 <- f_mh2
  p1[8]=h_h2; p2[10] = f_mh2;
  X = NumericVector::create(AS2, AI2, BS2,  BI2, CS2,  CI2, DS2,  DI2, MI2, NM2, HI2,NH2);
  // biting rates, mortality rates, infection rates floodplain
  NumericVector X2 = Rates_Updates(X, p1, p2);
  long double o_a2=X2[0],o_b2=X2[1],o_c2=X2[2], o_d2=X2[3], o_a2_2=X2[4], o_b2_2=X2[5], o_c2_2=X2[6], o_d2_2=X2[7], b_a2=X2[8],b_b2=X2[9], b_c2=X2[10], b_d2=X2[11], b_h2=X2[12], b_m2=X2[13];
  // Updates h_h2<- h_h3  and f_mh2 <- f_mh3
  p1[8]=h_h3; p2[10] = f_mh3, p2[15]=b_wl;
  X = NumericVector::create(AS3, AI3, BS3,  BI3, CS3,  CI3, DS3,  DI3, MI3, NM3, HI3,NH3);
  // biting rates, mortality rates, infection rates floodplain
  NumericVector X3 = Rates_Updates(X, p1, p2);
  long double o_a3=X3[0],o_b3=X3[1],o_c3=X3[2], o_d3=X3[3], o_a3_2=X3[4], o_b3_2=X3[5], o_c3_2=X3[6], o_d3_2=X3[7], b_a3=X3[8],b_b3=X3[9], b_c3=X3[10], b_d3=X3[11], b_h3=X3[12], b_m3=X3[13];
   // differential equations
   //---------------- human
   long double dHS1 = g_h*NH1 +  l_h21*HS2 + l_h31*HS3 + r_h*HR1 - (m_h + b_h1 +      l_h12 + l_h13)*HS1;
   long double dHE1 = b_h1*HS1 + l_h21*HE2 + l_h31*HE3 -           (m_h + x_h +       l_h12 + l_h13)*HE1;
   long double dHI1 = x_h*HE1 +  l_h21*HI2 + l_h31*HI3 -           (m_h + d_h + a_h + l_h12 + l_h13)*HI1;
   long double dHR1 = a_h*HI1 +  l_h21*HR2 + l_h31*HR3 -           (m_h + r_h +       l_h12 + l_h13)*HR1;
   long double dHS2 = g_h*NH2 +  l_h12*HS1 + l_h32*HS3 + r_h*HR2 - (m_h + b_h2 +      l_h21 + l_h23)*HS2;
   long double dHE2 = b_h2*HS2 + l_h12*HE1 + l_h32*HE3 -           (m_h + x_h +       l_h21 + l_h23)*HE2;
   long double dHI2 = x_h*HE2 +  l_h12*HI1 + l_h32*HI3 -           (m_h + d_h + a_h + l_h21 + l_h23)*HI2;
   long double dHR2 = a_h*HI2 +  l_h12*HR1 + l_h32*HR3 -           (m_h + r_h +       l_h21 + l_h23)*HR2;
   long double dHS3 = g_h*NH3 +  l_h13*HS1 + l_h23*HS2 + r_h*HR3 - (m_h + b_h3 +      l_h31 + l_h32)*HS3;
   long double dHE3 = b_h3*HS3 + l_h13*HE1 + l_h23*HE2 -           (m_h + x_h +       l_h31 + l_h32)*HE3;
   long double dHI3 = x_h*HE3 +  l_h13*HI1 + l_h23*HI2 -           (m_h + d_h + a_h + l_h31 + l_h32)*HI3;
   long double dHR3 = a_h*HI3 +  l_h13*HR1 + l_h23*HR2 -           (m_h + r_h +       l_h31 + l_h32)*HR3;
   NumericVector dH=NumericVector::create(dHS1,dHE1,dHI1,dHR1,dHS2,dHE2,dHI2,dHR2,dHS3,dHE3,dHI3,dHR3);
   //---------------- animal host
   NumericVector dMS1 = (g_m_u*NM1 + g_m_i*MI1)*Sup(0.0,1.0-NM1/k_m1) + l_m21*MS2 + l_m31*MS3 + r_m*MR1 - (m_m + b_m1      + l_m12 + l_m13)*MS1;
   NumericVector dME1 = b_m1*MS1                                      + l_m21*ME2 + l_m31*ME3           - (m_m + x_m       + l_m12 + l_m13)*ME1;
   NumericVector dMI1 = x_m*ME1                                       + l_m21*MI2 + l_m31*MI3           - (m_m + d_m + a_m + l_m12 + l_m13)*MI1;
   NumericVector dMR1 = a_m*MI1                                       + l_m21*MR2 + l_m31*MR3           - (m_m + r_m       + l_m12 + l_m13)*MR1;
   NumericVector dMS2 = (g_m_u*NM2 + g_m_i*MI2)*Sup(0.0,1.0-NM2/k_m2) + l_m12*MS1 + l_m32*MS3 + r_m*MR2 - (m_m + b_m2      + l_m21 + l_m23)*MS2;
   NumericVector dME2 = b_m2*MS2 +                                  l_m12*ME1 + l_m32*ME3 -           (m_m + x_m +       l_m21 + l_m23)*ME2;
   NumericVector dMI2 = x_m*ME2 +                                   l_m12*MI1 + l_m32*MI3 -           (m_m + d_m + a_m + l_m21 + l_m23)*MI2;
   NumericVector dMR2 = a_m*MI2 +                                   l_m12*MR1 + l_m32*MR3 -           (m_m + r_m +       l_m21 + l_m23)*MR2;
   long double   dMS3 = (g_m_u*NM3 + g_m_i*MI3)*Sup(0.0,1.0-NM3/k_m3) + l_m13*MS1 +     l_m23*MS2 + r_m*MR3 - (m_m + b_m3 +      l_m31 +     l_m32)*MS3;
   long double   dME3 = b_m3*MS3 +                                  l_m13*ME1 +     l_m23*ME2 -           (m_m + x_m +       l_m31 +     l_m32)*ME3;
   long double   dMI3 = x_m*ME3 +                                   l_m13*MI1 +     l_m23*MI2 -           (m_m + d_m + a_m + l_m31 +     l_m32)*MI3;
   long double   dMR3 = a_m*MI3 +                                   l_m13*MR1 +     l_m23*MR2 -           (m_m + r_m +       l_m31 +     l_m32)*MR3;
   NumericVector dM=combine(List::create(dMS1,dME1,dMI1,dMR1,dMS2,dME2,dMI2,dMR2,dMS3,dME3,dMI3,dMR3));
   //---------------- vector A --- zone 1 dormant
   long double temp=1.0-NA1/k_a1;
   NumericVector dAQ1 = o_a1*g_a*Sup(0.0,temp)*z_a*AI1 -                          (m_aq1 + ts*t_a1)        *AQ1;
   NumericVector dAP1 = g_a*Sup(0.0,temp)*(o_a1*(1.0-z_a)*AI1+(o_a1+o_a1_2)*AS1) -  (m_ap1 + ts*t_a1)        *AP1;
   NumericVector dAS1 = ts*t_a1*AP1 + l_a21*AS2 + l_a31*AS3 -                 (m_a + o_a1*b_a1 + l_a12 + l_a13)*AS1;
   NumericVector dAI1 = ts*t_a1*AQ1 + o_a1*b_a1*AS1 + l_a21*AI2 + l_a31*AI3 - (m_a             + l_a12 + l_a13)*AI1;
   temp = 1.0-NA2/k_a2;
   NumericVector dAQ2 = o_a2*g_a*Sup(0.0,temp)*z_a*AI2 -                          (m_aq2 + ts*t_a)             *AQ2;
   NumericVector dAP2 = g_a*Sup(0.0,temp)*(o_a2*(1.0-z_a)*AI2+(o_a2+o_a2_2)*AS2) -  (m_ap2 + ts*t_a)             *AP2;
   NumericVector dAS2 = ts*t_a*AP2 + l_a12*AS1 + l_a32*AS3 -                      (m_a + o_a2*b_a2 + l_a21 + l_a23)*AS2;
   NumericVector dAI2 = ts*t_a*AQ2 + o_a2*b_a2*AS2 + l_a12*AI1 + l_a32*AI3 -      (m_a             + l_a21 + l_a23)*AI2;
   temp = 1.0-NA3/k_a3;
   NumericVector dAQ3 = o_a3*g_a*Sup(0.0,temp)*z_a*AI3 -                          (m_aq3 + ts*t_a)             *AQ3;
   NumericVector dAP3 = g_a*Sup(0.0,temp)*(o_a3*(1-z_a)*AI3+(o_a3+o_a3_2)*AS3) -  (m_ap3 + ts*t_a)             *AP3;
   NumericVector dAS3 = ts*t_a*AP3 + l_a13*AS1 + l_a23*AS2 -                      (m_a + o_a3*b_a3 + l_a31 + l_a32)*AS3;
   NumericVector dAI3 = ts*t_a*AQ3 + o_a3*b_a3*AS3 + l_a13*AI1 + l_a23*AI2 -      (m_a             + l_a31 + l_a32)*AI3;
   NumericVector xtemp = ts*t_a1*AP1;
   NumericVector dA = combine(List::create(dAQ1,dAP1,dAS1,dAI1,dAQ2,dAP2,dAS2,dAI2,dAQ3,dAP3,dAS3,dAI3));
   //---------------- vector B
   NumericVector dBQ1 = o_b1*g_b*Sup(0.0,1.0-NB1/k_b1)*z_b*BI1 -                         (m_bq1 + ts*t_b)             *BQ1;
   NumericVector dBP1 = g_b*Sup(0.0,1.0-NB1/k_b1)*(o_b1*(1.0-z_b)*BI1+(o_b1+o_b1_2)*BS1) - (m_bp1 + ts*t_b)             *BP1;
   NumericVector dBS1 = ts*t_b*BP1 + l_b21*BS2 + l_b31*BS3 -                     (m_b + o_b1*b_b1 + l_b12 + l_b13)*BS1;
   NumericVector dBI1 = ts*t_b*BQ1 + o_b1*b_b1*BS1 + l_b21*BI2 + l_b31*BI3 -     (m_b             + l_b12 + l_b13)*BI1;
   NumericVector dBQ2 = o_b2*g_b*Sup(0.0,1.0-NB2/k_b2)*z_b*BI2 -                         (m_bq2 + ts*t_b)             *BQ2;
   NumericVector dBP2 = g_b*Sup(0.0,1.0-NB2/k_b2)*(o_b2*(1.0-z_b)*BI2+(o_b2+o_b2_2)*BS2) - (m_bp2 + ts*t_b)             *BP2;
   NumericVector dBS2 = ts*t_b*BP2 + l_b12*BS1 + l_b32*BS3 -                     (m_b + o_b2*b_b2 + l_b21 + l_b23)*BS2;
   NumericVector dBI2 = ts*t_b*BQ2 + o_b2*b_b2*BS2 + l_b12*BI1 + l_b32*BI3 -     (m_b             + l_b21 + l_b23)*BI2;
   NumericVector dBQ3 = o_b3*g_b*Sup(0.0,1.0-NB3/k_b3)*z_b*BI3 -                         (m_bq3 + ts*t_b)             *BQ3;
   NumericVector dBP3 = g_b*Sup(0.0,1.0-NB3/k_b3)*(o_b3*(1.0-z_b)*BI3+(o_b3+o_b3_2)*BS3) - (m_bp3 + ts*t_b)             *BP3;
   NumericVector dBS3 = ts*t_b*BP3 + l_b13*BS1 + l_b23*BS2 -                     (m_b + o_b3*b_b3 + l_b31 + l_b32)*BS3;
   NumericVector dBI3 = ts*t_b*BQ3 + o_b3*b_b3*BS3 + l_b13*BI1 + l_b23*BI2 -     (m_b             + l_b31 + l_b32)*BI3;
   NumericVector dB = combine(List::create(dBQ1,dBP1,dBS1,dBI1,dBQ2,dBP2,dBS2,dBI2,dBQ3,dBP3,dBS3,dBI3));
   //---------------- vector C
   NumericVector dCP1 = g_c*Sup(0.0,1.0-NC1/k_c1)*(o_c1*CI1+(o_c1+o_c1_2)*CS1) - (m_cp1 + ts*t_c)             *CP1;
   NumericVector dCS1 = ts*t_c*CP1 + l_c21*CS2 + l_c31*CS3 -             (m_c + o_c1*b_c1 + l_c12 + l_c13)*CS1;
   long double  dCI1 = o_c1*b_c1*CS1 + l_c21*CI2 + l_c31*CI3 -              (m_c             + l_c12 + l_c13)*CI1;
   NumericVector dCP2 = g_c*Sup(0.0,1.0-NC2/k_c2)*(o_c2*CI2+(o_c2+o_c2_2)*CS2) - (m_cp2 + ts*t_c)             *CP2;
   NumericVector dCS2 = ts*t_c*CP2 + l_c12*CS1 + l_c32*CS3 -             (m_c + o_c2*b_c2 + l_c21 + l_c23)*CS2;
   long double dCI2 = o_c2*b_c2*CS2 + l_c12*CI1 + l_c32*CI3 -              (m_c             + l_c21 + l_c23)*CI2;
   NumericVector dCP3 = g_c*Sup(0.0,1.0-NC3/k_c3)*(o_c3*CI3+(o_c3+o_c3_2)*CS3) - (m_cp3 + ts*t_c)             *CP3;
   NumericVector dCS3 = ts*t_c*CP3 + l_c13*CS1 + l_c23*CS2 -             (m_c + o_c3*b_c3 + l_c31 + l_c32)*CS3;
   long double dCI3 = o_c3*b_c3*CS3 + l_c13*CI1 + l_c23*CI2 -              (m_c             + l_c31 + l_c32)*CI3;
   NumericVector dC = combine(List::create(dCP1,dCS1,dCI1, dCP2,dCS2, dCI2, dCP3,dCS3,dCI3));
   //---------------- vector D
   NumericVector dDP1 = g_d*Sup(0.0,1.0-ND1/k_d1)*(o_d1*DI1+(o_d1+o_d1_2)*DS1) - (m_dp1 + ts*t_d)             *DP1;
   NumericVector dDS1 = ts*t_d*DP1 + l_d21*DS2 + l_d31*DS3 -             (m_d + o_d1*b_d1 + l_d12 + l_d13)*DS1;
   long double dDI1 = o_d1*b_d1*DS1 + l_d21*DI2 + l_d31*DI3 -              (m_d             + l_d12 + l_d13)*DI1;
   NumericVector dDP2 = g_d*Sup(0.0,1.0-ND2/k_d2)*(o_d2*DI2+(o_d2+o_d2_2)*DS2) - (m_dp2 + ts*t_d)             *DP2;
   NumericVector dDS2 = ts*t_d*DP2 + l_d12*DS1 + l_d32*DS3 -             (m_d + o_d2*b_d2 + l_d21 + l_d23)*DS2;
   long double dDI2 = o_d2*b_d2*DS2 + l_d12*DI1 + l_d32*DI3 -              (m_d             + l_d21 + l_d23)*DI2;
   NumericVector dDP3 = g_d*Sup(0.0,1.0-ND3/k_d3)*(o_d3*DI3+(o_d3+o_d3_2)*DS3) - (m_dp3 + ts*t_d)             *DP3;
   NumericVector dDS3 = ts*t_d*DP3 + l_d13*DS1 + l_d23*DS2 -             (m_d + o_d3*b_d3 + l_d31 + l_d32)*DS3;
   long double dDI3 = o_d3*b_d3*DS3 + l_d13*DI1 + l_d23*DI2 -              (m_d             + l_d31 + l_d32)*DI3;
   NumericVector dD = combine(List::create(dDP1,dDS1,dDI1,dDP2,dDS2,dDI2,dDP3,dDS3,dDI3));
   List ret=List::create(dH,dM,dA, dB, dC,dD);
   NumericVector xw = combine(ret);
   return  List::create(xw);
  }
