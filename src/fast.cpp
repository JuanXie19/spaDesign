#include <Rcpp.h>
using namespace Rcpp;
double lnorm(double t,int d) {
 if(t==0) return R::lgammafn(d/2.0)-log(2.0)-d/2.0*log(M_PI);
 return (d/2.0-1)*log(t)-log(R::bessel_i(t,d/2.0-1,2))-d/2.0*log(2*M_PI);
}
// [[Rcpp::export]]
List spa_estep_cpp(NumericMatrix x, NumericVector pi, NumericMatrix c,
 NumericVector r, NumericMatrix phi, NumericVector tau,double sigma) {
 int n=x.nrow(),d=x.ncol(),m=pi.size(); NumericMatrix W(n,m); List ys(m);
 double ll=0; NumericVector terms(m);
 if(sigma<=0) stop("sigma_sq must be positive");
 for(int k=0;k<m;k++) { if(tau[k]<=0) stop("tau must be positive"); ys[k]=NumericMatrix(n,d); }
 for(int i=0;i<n;i++) {
  if(i%256==0) checkUserInterrupt();
  double ws=0, mx=R_NegInf;
  for(int k=0;k<m;k++) {
   NumericMatrix y=ys[k]; double zz=0,nn=0; NumericVector nu(d);
   for(int j=0;j<d;j++) { double z=x(i,j)-c(k,j);zz+=z*z;nu[j]=tau[k]*phi(k,j)+r[k]*z/sigma;nn+=nu[j]*nu[j]; }
   double kap=std::max(sqrt(nn),1e-10);
   double val=exp(lnorm(tau[k],d)-d*log(2*M_PI*sigma)/2-lnorm(kap,d)
       -(std::max(zz,1e-20)+r[k]*r[k])/(2*sigma)+kap-tau[k]);
   W(i,k)=pi[k]*val;ws+=W(i,k);
   if(val<=0 || !R_finite(val)) val=1e-100;
   terms[k]=log(pi[k])+log(val);mx=std::max(mx,terms[k]);
   double ad=exp(log(R::bessel_i(kap,d/2.0,2))-log(R::bessel_i(kap,d/2.0-1,2)));
   ad=std::min(std::max(ad,1e-10),1-1e-10);
   for(int j=0;j<d;j++) y(i,j)=ad*nu[j]/kap;
  }
  double sm=0;for(int k=0;k<m;k++){W(i,k)/=ws;sm+=exp(terms[k]-mx);}ll+=mx+log(sm);
 }
 return List::create(_["W"]=W,_["y_expect_list"]=ys,_["log_likelihood"]=ll);
}
// [[Rcpp::export]]
List spa_mstep_cpp(NumericMatrix x,NumericMatrix W,List ys,NumericMatrix oldc,
 NumericVector oldr,NumericMatrix oldphi,NumericVector oldtau) {
 int n=x.nrow(),d=x.ncol(),m=W.ncol();
 NumericMatrix c=clone(oldc),phi=clone(oldphi);NumericVector r=clone(oldr),tau=clone(oldtau),pi(m);
 double sig=0;
 for(int k=0;k<m;k++){
  checkUserInterrupt();NumericMatrix y=ys[k];double sw=0;NumericVector sy(d);
  for(int i=0;i<n;i++)sw+=W(i,k);pi[k]=sw/n;
  for(int j=0;j<d;j++){double a=0;for(int i=0;i<n;i++){a+=W(i,k)*(x(i,j)-oldr[k]*y(i,j));sy[j]+=W(i,k)*y(i,j);}c(k,j)=a/sw;}
  double a=0;for(int i=0;i<n;i++)for(int j=0;j<d;j++)a+=W(i,k)*(x(i,j)-c(k,j))*y(i,j);
  r[k]=a/sw;double norm=0;for(int j=0;j<d;j++)norm+=sy[j]*sy[j];norm=std::max(sqrt(norm),1e-10);
  for(int j=0;j<d;j++)phi(k,j)=sy[j]/norm;
  double rb=norm/sw;tau[k]=std::max((rb*d-rb*rb*rb)/(1-rb*rb),1e-6);
  for(int i=0;i<n;i++){double zz=0,zy=0;for(int j=0;j<d;j++){double z=x(i,j)-c(k,j);zz+=z*z;zy+=z*y(i,j);}sig+=W(i,k)*(zz-2*r[k]*zy+r[k]*r[k]);}
 }
 return List::create(_["pi_hat"]=pi,_["c_hat"]=c,_["r_hat"]=r,_["phi_hat"]=phi,_["tau_hat"]=tau,_["sigma_sq"]=sig/(n*d));
}
// [[Rcpp::export]]
NumericMatrix spa_gene_cpp(NumericMatrix x,IntegerVector inside) {
 int g=x.nrow(),n=x.ncol();std::vector<bool> in(n,false);
 for(int j:inside)in[j-1]=true;NumericMatrix out(g,5);
 for(int i=0;i<g;i++){
  if(i%256==0)checkUserInterrupt();
  double si=0,so=0;int ni=0;std::vector<double> v;
  for(int j=0;j<n;j++){if(in[j]){si+=x(i,j);ni++;}else{so+=x(i,j);v.push_back(x(i,j));}}
  double mi=ni?si/ni:R_NaN,mo=v.size()?so/v.size():R_NaN,ml=R_NaN;
  if(v.size()){std::sort(v.begin(),v.end());int z=v.size();double med=z%2?v[z/2]:(v[z/2-1]+v[z/2])/2;double sl=0;int nl=0;for(double a:v)if(a<=med){sl+=a;nl++;}ml=sl/nl;}
  out(i,0)=mi-mo;out(i,1)=mi-ml;out(i,2)=mi;out(i,3)=mo;out(i,4)=ml;
 }
 return out;
}
