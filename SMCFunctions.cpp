#include <Rcpp.h>
#include <deque>
#include <vector>
#include <numeric>
using namespace Rcpp;
using namespace std;

//[[Rcpp::export]]
double logFactorial_fun(int k){
  double temp;
  IntegerVector temp_vec;
  if(k==0){
    temp=0;
  }else{
    temp_vec=seq(1,k);
    temp=sum(log(temp_vec));
  }
  return(temp);
}



//[[Rcpp::export]]
List init_fun(double Tend, double hist_end){
  
  List temp;
  for(int i=0;i<Tend;i++){
    temp.push_back(NA_REAL);
  }
  NumericVector tempv={};
  int i=0;
  while(i<=(Tend*7)){
    tempv.push_back(i);
    i=i+7;
  }
  tempv=tempv+(hist_end*7);
  temp.push_back(tempv);
  return (temp);
}
//[[Rcpp::export]]
List initAge_fun(double Tend){

  List temp;
  for(int i=0;i<Tend;i++){
    temp.push_back(NA_REAL);
  }

  return (temp);
}

//[[Rcpp::export]]
List initXSa_fun(int Num_Agroups, double Tend, double hist_end){
  
  List temp;
  for(int j=0;j<Num_Agroups;j++){
    List temp1;
    for(int i=0;i<Tend;i++){
      temp1.push_back(NA_REAL);
    }
    NumericVector tempv={};
    int i=0;
    while(i<=(Tend*7)){
      tempv.push_back(i);
      i=i+7;
    }
    tempv=tempv+(hist_end*7);
    temp1.push_back(tempv);
    temp.push_back(temp1);
  }
  return (temp);
}

// //[[Rcpp::export]]
// NumericVector which_fun(NumericVector source, double key,NumericVector ind){
//   int temp=ind.length();
//   for(int i=temp;i<source.length();i++){
//     if(source[i]<key){
//       ind.push_back(i);
//     }
//   }
//   return(ind);
// }
// 
//[[Rcpp::export]]
int find_interval(double t,double Tend,double hist_end){
  int inda;
  if((t/7-hist_end)>=(Tend-1)){
    inda=Tend-1;
  }else{
    inda=floor(t/7-hist_end);
  }
  return(inda);
}


//[[Rcpp::export]]
double rgammat_fun(int n, NumericVector range, double shape, double rate){
  double F_a=R::pgamma(range[0],shape,1/rate,true,false);
  double F_b=R::pgamma(range[1],shape,1/rate,true,false);
  double u=R::runif(F_a,F_b);
  double temp=R::qgamma(u,shape,1/rate,true,false);

  return (temp);


}


//[[Rcpp::export]]
List ThinningFunCl_prop(int IntervalInd,List T_N,List AgeT_N, double a_g,double b_g,NumericMatrix alpha,List seeds, List age_seeds, NumericMatrix PM,NumericVector age_vec,NumericVector Susc_vec, NumericVector Pop_vec){
  NumericVector IntervLim=T_N[(T_N.length()-1)],range,TN,lambda_vec,range_par,AS,TN1,TN2;
  double Tend=IntervLim[(IntervalInd+1)], Tstart=IntervLim[IntervalInd],T_end=IntervLim[(IntervLim.length()-1)],I_j;
  int Num_off, Num_Agroups, IndAgInf;
  double lambda_par,tstar,pc, lambda_a;
  NumericVector SuscVec=clone(Susc_vec);
  std::deque<double> Q;
  std::deque<double> Qa;
  std::vector <double> TN_prop, AgeTN_prop;
  int Ind_InfInterval= IntervalInd-seeds.length();
  if(Ind_InfInterval<0){
    for(int j=IntervalInd;j<seeds.length();j++){
      TN=seeds[j];
      if(R_IsNA(TN[0])){
        TN={};
      }
      if(TN.length()>0){
        AS=age_seeds[j];
        Q.insert(Q.end(),TN.begin(),TN.end());
        Qa.insert(Qa.end(),AS.begin(),AS.end());
      }
    }
  }
  NumericVector IndInfInterval_vec={0,Ind_InfInterval};
  for(int i=max(IndInfInterval_vec);i<IntervalInd;i++){
    TN=T_N[i];
    if(R_IsNA(TN[0])){
      TN={};
    }
    if(TN.length()>0){
      AS=AgeT_N[i];
      Q.insert(Q.end(),TN.begin(),TN.end());
      Qa.insert(Qa.end(),AS.begin(),AS.end());
    }
  }
  while(!Q.empty()){
    range_par={Tstart,Q.front()};
    lambda_vec={0,Tstart-Q.front()};
    I_j= R::pgamma(Tend-Q.front(),a_g,1/b_g,true,false) - R::pgamma(max(lambda_vec),a_g,1/b_g,true,false);
    Num_Agroups=PM.nrow();
    IndAgInf=Qa.front()-1;
    lambda_a=I_j*sum( alpha(_,IntervalInd) * ( (SuscVec/Pop_vec)*PM(_,IndAgInf) ) );
    NumericVector NumOff_vec={R::rpois(lambda_a),sum(SuscVec)};
    Num_off= min(NumOff_vec);
    range={max(lambda_vec),Tend-Q.front()};
    NumericVector p_sampl=( alpha(_,IntervalInd)*( (SuscVec/Pop_vec)*PM(_,IndAgInf) ) )/sum( alpha(_,IntervalInd)*( (SuscVec/Pop_vec)*PM(_,IndAgInf) ) );
    for(int i=0;i<Num_off;++i){
      tstar=Q.front()+rgammat_fun(1,range,a_g,b_g);
      if(tstar>=max(range_par) & tstar<Tend ){
          //Rprintf("(Q:%f,tstar=%f,L=%f, U=%f, rgamma:%f), ",Q.front(),tstar,range[0],range[1],tg);
        //NumericVector p_sampl=PM(_,IndAgInf)/sum(PM(_,IndAgInf));
        NumericVector flag_g=sample(age_vec,1,true,p_sampl);
        if(Susc_vec[flag_g[0]]>0){
          Q.push_back(tstar);
          Qa.push_back((flag_g[0]+1));
          TN_prop.push_back(tstar);
          AgeTN_prop.push_back((flag_g[0]+1));
          if(flag_g[0]==0){
            TN1.push_back(tstar);
            SuscVec[0]=SuscVec[0]-1;
          }
          if(flag_g[0]==1){
            TN2.push_back(tstar);
            SuscVec[1]=SuscVec[1]-1;
          }
          // if(flag_g[0]==2){
          //   TN3.push_back(tstar);
          //   SuscVec[2]=SuscVec[2]-1;
          // }
          // if(flag_g[0]==3){
          //   TN4.push_back(tstar);
          //   SuscVec[3]=SuscVec[3]-1;
          // }
        }
      }
    }
    Q.pop_front();//erase(Q.begin());
    Qa.pop_front();
  }
  
  //std::sort(TN_prop.begin(),TN_prop.end());
  //auto last = std::unique(TN_prop.begin(), TN_prop.end());
  //TN_prop.erase(last, TN_prop.end());
  NumericVector TN_prop_Nvec(TN_prop.begin(),TN_prop.end());
  NumericVector AgeTN_prop_Nvec(AgeTN_prop.begin(),AgeTN_prop.end());
  
  List results;
  results.push_back(TN_prop_Nvec);
  results.push_back(AgeTN_prop_Nvec);
  results.push_back(TN1);
  results.push_back(TN2);
  //results.push_back(TN3);
  //results.push_back(TN4);
  results.push_back(SuscVec);
  return (results);
}



// 
//[[Rcpp::export]]
List Initialization_FixedParameters(double min_d, double max_d,double min_v,double max_v,int num_particles){
  List results;
  NumericVector theta_d,theta_v;
  double mean_d,sd_d,mean_v,sd_v;
  mean_d=( log(max_d) + log(min_d) )/2;
  sd_d=( log(max_d) - log(min_d) )/8;
  mean_v=( log(max_v) + log(min_v) )/2;
  sd_v=( log(max_v) - log(min_v) )/8 ;
  theta_d=exp( Rcpp::rnorm(num_particles,mean_d, sd_d ) );
  theta_v=exp( Rcpp::rnorm(num_particles,mean_v, sd_v ) );
  //theta_d=Rcpp::runif(num_particles,min_d, max_d );
  //theta_v=Rcpp::runif(num_particles,min_v, max_v );
  results.push_back(theta_d);
  results.push_back(theta_v);
  return(results);
}
//
//
// // //
//[[Rcpp::export]]
List Initialization_step(double a, double b,double a_g,double b_g,int num_particles, int num_intervals,List seeds, List age_seeds, double hist_end, NumericMatrix PM,NumericVector age_vec,NumericVector Susc_vec,NumericVector Pop_vec){
  List results, X_S, X_SAge,X_Sa,X_Susc,X_a;
  List temp,tempa;
  for(int k=0;k<num_particles;k++){
    NumericMatrix Xa(PM.nrow(),num_intervals+1);
    Xa(_,0)=Rcpp::runif(PM.nrow(), a, b);
    List tempa;
  // Rprintf("k=%d\n",k);
    List TN=init_fun(num_intervals+1,2*hist_end);
    List AgeT_N=initAge_fun(num_intervals+1);
    tempa=initXSa_fun(PM.nrow(),num_intervals+1,2*hist_end);
// //List ThinningFunCl_prop(int IntervalInd,List T_N,List AgeT_N, double a_g,double b_g,NumericVector alpha,List seeds, List age_seeds, NumericMatrix PM,NumericVector age_vec,NumericVector Susc_vec, NumericVector Pop_vec){
    temp= ThinningFunCl_prop(0,TN,AgeT_N,a_g,b_g,Xa,seeds,age_seeds, PM,age_vec,Susc_vec,Pop_vec);
    TN[0]=temp[0];
    AgeT_N[0]=temp[1];
    Rcpp::as<Rcpp::List>(tempa[0])[0]=temp[2];
    Rcpp::as<Rcpp::List>(tempa[1])[0]=temp[3];
    //Rcpp::as<Rcpp::List>(tempa[2])[0]=temp[4];
    //Rcpp::as<Rcpp::List>(tempa[3])[0]=temp[5];
    X_S.push_back(TN);
    X_SAge.push_back(AgeT_N);
    X_Sa.push_back(tempa);
    X_Susc.push_back(temp[4]);
    X_a.push_back(Xa);
  }
  results.push_back(X_S);
  results.push_back(X_a);
  results.push_back(X_SAge);
  results.push_back(X_Sa);
  results.push_back(X_Susc);
  return(results);
}

//[[Rcpp::export]]
double findMu_fun(int IntervalInd, double a_g, double b_g, double beta,List T_N, double v_nb,List seeds){
  NumericVector v,vstart,vend,I_j;
  double mu=0,log_LR,lambda_par;

  NumericVector IntervLim= T_N[(T_N.length()-1)];
  double Tend=IntervLim[(IntervalInd+1)], Tstart=IntervLim[IntervalInd];

  int Ind_InfInterval= IntervalInd-seeds.length();
  NumericVector IndInfInterval_vec={0,Ind_InfInterval};

  if(Ind_InfInterval<0){
    for(int j=IntervalInd;j<seeds.length();j++){
      v=seeds[j];
      if(R_IsNA(v[0])){
        v={};
      }
      if(v.length()>0){
        vend=Tend-v;
        vstart=Tstart-v;
        I_j= Rcpp::pgamma(vend,a_g,1/b_g,true,false) - Rcpp::pgamma(vstart,a_g,1/b_g,true,false);
        mu=mu + sum(I_j);
      }
    }
  }

  for(int j=max(IndInfInterval_vec);j<IntervalInd;j++){
    v=T_N[j];
    if(R_IsNA(v[0])){
      v={};
    }
    if(v.length()>0){
      vend=Tend-v;
      vstart=Tstart-v;
      I_j= Rcpp::pgamma(vend,a_g,1/b_g,true,false) - Rcpp::pgamma(vstart,a_g,1/b_g,true,false);
      mu=mu + sum(I_j);
    }
  }
  v=T_N[IntervalInd];
  if(v.length()>0){
    vend=Tend-v;
    I_j= Rcpp::pgamma(vend,a_g,1/b_g,true,false);
    mu=mu + sum(I_j);
  }
  mu= beta*mu;
  return (mu);
}

//[[Rcpp::export]]
double log_likelihood_fun(int IntervalInd, double a_g, double b_g, double beta, double Y,List T_N,double v_nb,List seeds){
  NumericVector v,vstart,vend,I_j;
  double mu=0,log_LR,lambda_par;

  NumericVector IntervLim= T_N[(T_N.length()-1)];
  double Tend=IntervLim[(IntervalInd+1)], Tstart=IntervLim[IntervalInd];

  int Ind_InfInterval= IntervalInd-seeds.length();
  NumericVector IndInfInterval_vec={0,Ind_InfInterval};

  if(Ind_InfInterval<0){
    for(int j=IntervalInd;j<seeds.length();j++){
      v=seeds[j];
      if(R_IsNA(v[0])){
        v={};
      }
      if(v.length()>0){
        vend=Tend-v;
        vstart=Tstart-v;
        I_j= Rcpp::pgamma(vend,a_g,1/b_g,true,false) - Rcpp::pgamma(vstart,a_g,1/b_g,true,false);
        mu=mu + sum(I_j);
      }
    }
  }

  for(int j=max(IndInfInterval_vec);j<IntervalInd;j++){
    v=T_N[j];
    if(R_IsNA(v[0])){
      v={};
    }
    if(v.length()>0){
      vend=Tend-v;
      vstart=Tstart-v;
      I_j= Rcpp::pgamma(vend,a_g,1/b_g,true,false) - Rcpp::pgamma(vstart,a_g,1/b_g,true,false);
      mu=mu + sum(I_j);
    }
  }
  v=T_N[IntervalInd];
  if(v.length()>0){
    vend=Tend-v;
    I_j= Rcpp::pgamma(vend,a_g,1/b_g,true,false);
    mu=mu + sum(I_j);
  }

  mu= beta*mu;

  log_LR = lgamma( Y + (1/v_nb)) - logFactorial_fun(Y) -lgamma(1/v_nb) - ( (1/v_nb)*log(1+(mu*v_nb)) ) + (Y*log(v_nb*mu)) - (Y*log((mu*v_nb)+1)) ;
  if(mu==0){
    log_LR=-740;
  }
  if(mu==0 & Y==0){
    log_LR=0;
  }
  return (log_LR);
}

//[[Rcpp::export]]
NumericVector find_Inf(double t,int indt, List X_Sa, List XSa_seeds,List seeds, NumericMatrix PM){
  
  NumericVector TN,XInf;
  List TNa, XSa;
  int Ind_InfInterval=indt-seeds.length();
  NumericVector IndInfInterval_vec={0,Ind_InfInterval};
  
  for(int i=0; i<PM.nrow();i++){
    TNa=X_Sa[i];
    XSa=XSa_seeds[i];
    double temp=0;
    
    if(Ind_InfInterval<0){
      for(int j=indt;j<seeds.length();j++){
        TN=XSa[j];
        if(R_IsNA(TN[0])){
          TN={};
        }
        if(TN.length()>0){
          temp=temp+TN.length();
        }
      }  
    }
    
    for(int j=max(IndInfInterval_vec);j<indt;j++){
      TN=TNa[j];
      if(R_IsNA(TN[0])){
        TN={};
      }
      if(TN.length()>0){
        temp=temp+TN.length();
      }
    }
    
    TN=TNa[indt];
    for(int j=0;j<TN.length();j++){
      if(TN[j]<=t){
        temp=temp+1;  
      }  
    }
    
    XInf.push_back(temp);
  }
    
  return XInf;
}

// //
//[[Rcpp::export]]
NumericVector weights_fun(NumericVector Y,List X_S,int IntervalInd,double beta, double a_g,double b_g,int num_particles,List seeds,NumericVector theta_v,int Num_Agroups){
  NumericVector logw={}, tempv;
  double temp;
  for(int j=0;j<num_particles;j++){
    temp=0;
    List XSa=X_S[j];
    for(int ia=0;ia<Num_Agroups;ia++){
      List tempa=XSa[ia];
      List seedsa=seeds[ia];
      temp=temp+log_likelihood_fun(IntervalInd, a_g,b_g, beta,Y[ia], tempa,theta_v[j],seedsa);
    }
    logw.push_back(temp);
  }
  return (logw);
}

// // //List ThinningFunCl_prop(int IntervalInd,List T_N, List AgeT_N, double a_g,double b_g,NumericVector alpha,List seeds, List age_seeds, NumericMatrix PM, NumericVector S_Ind, NumericVector Population){
// //
// 
// // //List ThinningFunCl_prop(int IntervalInd,List T_N, List AgeT_N, double a_g,double b_g,NumericVector alpha,List seeds, List age_seeds, NumericMatrix PM, NumericVector S_Ind, NumericVector Population){
// //
//[[Rcpp::export]]
double AuxMultP_fun(NumericVector Y,List T_N,List TN_Age,List TN_a,int IntervalInd,double beta, double a_g,double b_g, double c_g, double d_g, NumericMatrix alpha_mat,double v_nb, List seeds,List age_seeds,List XS_seeds,NumericMatrix PM,double c,double d, NumericVector age_vec, NumericVector Susc_vec, NumericVector Pop_vec){
  NumericVector tempv;
  List TN=clone(T_N);
  List TNAge=clone(TN_Age);
  List TNa=clone(TN_a);
  List res;
  List temp_TN=clone(T_N);
  List temp_TNAge=clone(TN_Age);
  List temp_TNa=clone(TN_a);
  NumericVector SuscVec=clone(Susc_vec);
//for(int j=0;j<1;j++){
  NumericMatrix alpha=clone(alpha_mat);
  for(int j=0;j<PM.nrow();j++){
    alpha(j,IntervalInd)=R::rgamma(c, alpha(j,IntervalInd-1)/d );
  }
  //List ThinningFunCl_prop(int IntervalInd,List T_N,List AgeT_N, double a_g,double b_g,NumericVector alpha,List seeds, List age_seeds, NumericMatrix PM,NumericVector age_vec,NumericVector Susc_vec, NumericVector Pop_vec){
  res=ThinningFunCl_prop(IntervalInd, TN,TNAge,a_g,b_g,alpha,seeds,age_seeds,PM,age_vec,SuscVec,Pop_vec);
  temp_TN[IntervalInd]=res[0];
  temp_TNAge[IntervalInd]=res[1];
  Rcpp::as<Rcpp::List>(temp_TNa[0])[IntervalInd]=res[2];
  Rcpp::as<Rcpp::List>(temp_TNa[1])[IntervalInd]=res[3];
  //Rcpp::as<Rcpp::List>(temp_TNa[2])[IntervalInd]=res[4];
  //Rcpp::as<Rcpp::List>(temp_TNa[3])[IntervalInd]=res[5];
  //
  double temp=0;
  for(int ia=0;ia<PM.nrow();ia++){
    //double log_likelihood_fun(int IntervalInd, double a_g, double b_g, double beta, double Y,List T_N,double alpha, double v_nb,List seeds){
    temp=temp + log_likelihood_fun(IntervalInd, c_g,d_g, beta,Y[ia],temp_TNa[ia],v_nb,XS_seeds[ia]);
  }
  return (temp);
}

//
//[[Rcpp::export]]
NumericVector AuxMult_fun(NumericVector Y,List X_S,List X_SAge,List X_Sa, int IntervalInd,double beta, double a_g,double b_g,double c_g, double d_g,int num_particles,List X_a, List seeds,List age_seeds,List XSa_seeds,NumericMatrix PM, NumericVector mv, NumericVector theta_d,NumericVector age_vec,List X_Susc,NumericVector Pop_vec ){
  NumericVector AuxVec;
 for(int j=0;j<num_particles;j++){
    AuxVec.push_back(AuxMultP_fun(Y,X_S[j],X_SAge[j],X_Sa[j],IntervalInd,beta,a_g,b_g,c_g,d_g,X_a[j],mv[j],seeds,age_seeds,XSa_seeds, PM,theta_d[j],theta_d[j],age_vec,X_Susc[j],Pop_vec));
  }
  return (AuxVec);
}
//[[Rcpp::export]]
List lambda_fun(List T_N,List T_NAge,List T_Na,NumericMatrix alpha, double a_g,double b_g, double t, List seeds,List age_seeds,List XSa_seeds,int Tend,NumericMatrix PM, int a_PM,NumericVector Susc_vec,NumericVector Pop_vec){
  List res;
  double lambda=0,tempp;
  NumericVector a_gv={a_g};
  int inda;
  NumericVector TN,TNAge,ind,temp,IntervLim=T_N[(T_N.length()-1)];
  double T_end=IntervLim[(IntervLim.length()-1)];
  inda=find_interval(t,Tend,2*seeds.length());

  int S_ta=Susc_vec[a_PM];

  for(int j=0;j<inda;j++){
    TN=T_Na[j];
    if(R_IsNA(TN[0])){
      TN={};
    }
    S_ta=S_ta-TN.length();
  }
  TN=T_Na[inda];
  if(R_IsNA(TN[0])){
    TN={};
  }
  for(int j=0;j<TN.length();j++){
    if(TN[j]<t){
      S_ta=S_ta-1;
    }
  }


  int Ind_InfInterval= inda-seeds.length();
  NumericVector IndInfInterval_vec={0,Ind_InfInterval};
  int ia;
  if(Ind_InfInterval<0){
    for(int j=inda;j<seeds.length();j++){
      TN=seeds[j];
      TNAge=age_seeds[j];
      for(int i=0;i<TN.length();i++){
        ia=TNAge[i]-1;
        lambda=lambda + ( alpha(a_PM,inda)* ( pow(b_g,a_g)/gamma(a_gv)[0] ) * pow(t-TN[i],a_g-1)  * exp(-1*b_g*(t-TN[i]) ) * PM(a_PM,ia) ) ;
      }
    }
  }

  for(int j=max(IndInfInterval_vec);j<(T_N.length()-1);j++){
    TN=T_N[j];
    TNAge=T_NAge[j];
    if(R_IsNA(TN[0])){
      TN={};
    }
    if(t>=IntervLim[j+1]){
      if(TN.length()>0 ){
        for(int i=0;i<TN.length();i++){
          ia=TNAge[i]-1;
          lambda=lambda + ( alpha(a_PM,inda)* ( pow(b_g,a_g)/gamma(a_gv)[0] ) *  pow(t-TN[i],a_g-1)  * exp(-1*b_g*(t-TN[i]) ) * PM(a_PM,ia) );
        }
      }
    }else{
      for(int i=0;i<TN.length();i++){
        if(TN[i]<t ){
            ia=TNAge[i]-1;
            lambda=lambda + ( alpha(a_PM,inda)* ( pow(b_g,a_g)/gamma(a_gv)[0] ) *   pow(t-TN[i],a_g-1)  * exp(-1*b_g*(t-TN[i]) ) * PM(a_PM,ia) );
        }
      }
      break;
    }
  }
  lambda=lambda*(S_ta/Pop_vec[a_PM]);
  res.push_back(lambda);
  res.push_back(S_ta);
  return res;
}

 
//[[Rcpp::export]]
List resampling(List X_S,List X_SAge,List X_Sa,List X_a, NumericVector W, int num_particles,NumericVector indi,List X_Susc){
  List res,Xs,XSAge,XSa,XSusc,Xa;
  NumericVector indv=sample(indi, num_particles,true,W);
  Xs=X_S[indv];
  XSAge=X_SAge[indv];
  XSa=X_Sa[indv];
  XSusc=X_Susc[indv];
  Xa=X_a[indv];
  
  res.push_back(Xs);
  res.push_back(Xa);
  res.push_back(XSAge);
  res.push_back(XSa);
  res.push_back(XSusc);
  res.push_back(indv);
  return(res);
}
// // // //
//[[Rcpp::export]]
List RegeneratingPar_fun(NumericVector mv, NumericVector md,NumericVector Ltheta_v,NumericVector Ltheta_d, NumericVector indv, double h_m, NumericVector w){
  double V1=sum(w);
  double V2=sum( pow(w,2) );

  double temp=V1/(pow(V1,2)-V2);

  double mu_d=sum(w*Ltheta_d);
  double Var_d=temp*sum( w* pow( Ltheta_d-mu_d ,2 )  );

  double mu_v=sum(w*Ltheta_v);
  double Var_v=temp*sum( w* pow(Ltheta_v-mu_v,2 ) );
  List results;

  NumericVector th_d, th_v;
  for(int j=0; j<mv.length(); j++){
    int i=indv[j];
    th_d.push_back(  R::rnorm(md[i],  h_m*pow(Var_d,0.5) ) ) ;
    th_v.push_back(  R::rnorm(mv[i],  h_m*pow(Var_v,0.5) )  );
  }

  results.push_back(th_d);
  results.push_back(th_v);
  return (results);

}

//[[Rcpp::export]]
List mu_fun(NumericVector Ltheta_d, NumericVector Ltheta_v,NumericVector theta_d,NumericVector theta_v,NumericVector w,double a_m){

  List results;

  double Lmu_d=sum(w*Ltheta_d);
  NumericVector Lmd= ( a_m*Ltheta_d ) + ( (1-a_m)*Lmu_d ) ;
  double Lmu_v=sum(w*Ltheta_v);
  NumericVector Lmv= ( a_m*Ltheta_v ) + ( (1-a_m)*Lmu_v ) ;

  double mu_d=sum(w*theta_d);
  NumericVector md= ( a_m*theta_d ) + ( (1-a_m)*mu_d ) ;
  double mu_v=sum(w*theta_v);
  NumericVector mv= ( a_m*theta_v ) + ( (1-a_m)*mu_v ) ;

  results.push_back(md);
  results.push_back(mv);
  results.push_back(Lmd);
  results.push_back(Lmv);

  return (results);
}
