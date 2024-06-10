propagation<-function(X_a,X_S,X_SAge,X_Sa,PM,num_particles,IntervalInd,a_g,b_g,seeds,age_seeds,theta_d, age_vec,X_Susc,Pop_vec){
  for(j in seq(1,num_particles,1)){
    c=theta_d[j]
    d=c
    Xa=X_a[[j]]
    for(i in seq(1,dim(PM)[1],1)){
      preva=Xa[i,IntervalInd-1]
      Xa[i,IntervalInd]=rgamma(1,shape=c,rate= d/preva )
    }
    TN=X_S[[j]]
    AgeT_N=X_SAge[[j]]
    XSa=X_Sa[[j]]
    XSusc=X_Susc[[j]]
    #List ThinningFunCl_prop(int IntervalInd,List T_N,List AgeT_N, double a_g,double b_g,NumericVector alpha,List seeds, List age_seeds, NumericMatrix PM,NumericVector age_vec,NumericVector Susc_vec, NumericVector Pop_vec){
    temp= ThinningFunCl_prop(IntervalInd-1,TN,AgeT_N,a_g,b_g,Xa,seeds,age_seeds,PM,age_vec,XSusc,Pop_vec)
    TN[[IntervalInd]]=temp[[1]]
    AgeT_N[[IntervalInd]]=temp[[2]]
    XSa[[1]][[IntervalInd]]=temp[[3]]
    XSa[[2]][[IntervalInd]]=temp[[4]]
    #XSa[[3]][[IntervalInd]]=temp[[5]]
    #XSa[[4]][[IntervalInd]]=temp[[6]]
    XSusc=temp[[5]]
    
    X_S[[j]]=TN
    X_SAge[[j]]=AgeT_N
    X_Sa[[j]]=XSa
    X_a[[j]]=Xa
    X_Susc[[j]]=XSusc
    # 
  }
  res<-list()
  res[[1]]=X_S
  res[[2]]=X_a
  res[[3]]=X_SAge
  res[[4]]=X_Sa
  res[[5]]=X_Susc
  return(res)
  
}
