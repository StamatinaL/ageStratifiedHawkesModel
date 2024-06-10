rm(list=ls())
require(Rcpp)
set.seed(100)
#library(qbld)
#library(R.utils)
sourceCpp("SMCFunctions.cpp")
source("MH_functions.R")

load("Data_Scenario6.RData")
Y<-Saved_Data$Y
alpha<-Saved_Data$alpha
dim(alpha)
T_N<-Saved_Data$T_N
AgeT_N<-Saved_Data$Age_TN
T_Na<-Saved_Data$T_Na
seeds<-Saved_Data$seedsE
age_seeds<-Saved_Data$age_seedsE
XSa_seeds<-Saved_Data$XSa_seedsE
Susc_vec<-Saved_Data$Susc_vecE
Pop_vec<-Saved_Data$Pop_vec
c=Saved_Data$c
d=c
v_nb<-Saved_Data$v_nb
PM<-Saved_Data$PM
d
v_nb
PM
round(PM,2)
sum(rowSums(Y))
colSums(Y)
#Saved_Data$seeds
a=0
b=0.5

length(XSa_seeds[[1]][[1]])+length(XSa_seeds[[1]][[2]])+length(XSa_seeds[[1]][[3]]) 
length(XSa_seeds[[2]][[1]])+length(XSa_seeds[[2]][[2]])+length(XSa_seeds[[2]][[3]]) 


Susc_vec

age_groups<-seq(1,2,1)

num_particles=30000 # number of particles
Tend=17
num_intervals=Tend
# Latent Cases
bg=6.7/(1.8^2)
ag=6.7*bg 
ag/bg
sqrt(ag/(bg^2))

# Obesrved Cases
dg=8.8/(4.1^2)
cg=8.8*dg
cg/dg
sqrt(cg/(dg^2))

# denObs<-density(rgamma(1000000,shape=cg,rate=dg))
# ind<-which(denObs$y==max(denObs$y))
# denObs$x[ind]
# plot(density(rgamma(10000,shape=cg,rate=dg)))


age_vec=seq(0,1,1)
age_vec
num_intervals=Tend

# library(ggplot2)
# df=data.frame(x=1:17,y=Y[,1])
# ggplot(df,aes(x = x, y = y))+theme_bw()+xlab("")+ylab("") + 
#   theme(text = element_text(size=30)) +
#   geom_line(col='red') +
#   geom_line(aes(x=(1:(num_intervals)),y=Y[,2]),color='blue')



hist_end=length(seeds)
hist_end

min_d=10
max_d=20
min_v=0.0001
max_v=0.5
res<-list()
#List Initialization_FixedParameters(double min_d, double max_d,double min_v,double max_v,int num_particles){
res<-Initialization_FixedParameters(min_d,max_d,min_v,max_v,num_particles)
theta_d<-res[[1]]
theta_v<-res[[2]]

Ltheta_d<-log(theta_d)
Ltheta_v<-log(theta_v)
Particles<-list()
res<-list()
#List Initialization_step(double a, double b,double a_g,double b_g,int num_particles, int num_intervals,List seeds, List age_seeds, double hist_end, NumericMatrix PM,NumericVector age_vec,NumericVector Susc_vec,NumericVector Pop_vec){
res<-Initialization_step(a,b,ag,bg,num_particles,num_intervals,seeds,age_seeds,hist_end, PM,age_vec,Susc_vec,Pop_vec)
Particles[[1]]<-res
X_S<-res[[1]]
X_a<-res[[2]]
X_SAge<-res[[3]]
X_Sa<-res[[4]]
X_Susc<-res[[5]]
Susc_vec
beta=0.5
ESS<-c()
ESSw<-c()
logl<-weights_fun(Y[1,],X_Sa,0,beta,cg,dg,num_particles,XSa_seeds,theta_v,length(age_groups))
logwt=logl
maxwt=max(logwt)
temp=logwt-maxwt
wt=exp(temp)
w=wt/sum(wt)
logw<-log(w)
ESSw<-c(ESSw,1/sum(w^2))
ESSw

D_m=0.99
h_m2=1 - ( ( ((3*D_m)-1)/(2*D_m)) ^2 )
h_m=sqrt(h_m2)
a_m2=1-h_m2
a_m=sqrt(a_m2)
a_m
h_m

indi<-seq(0,num_particles-1,1)
indv=indi+1
log_g<-rep(0,num_particles)
#Weights<-c()
X_d<-c()
X_v<-c()
X_dw<-c()
X_vw<-c()
X_dw<-rbind(X_dw,theta_d)
X_vw<-rbind(X_vw,theta_v)
 
for(IntervalInd in 1:(Tend-2)){
  
  print(IntervalInd)
  
  indi<-seq(0,num_particles-1,1)
  indv=indi+1
  
  m_res<-list()
  m_res=mu_fun(Ltheta_d,Ltheta_v,theta_d,theta_v,w,a_m)
  md=m_res[[1]]
  mv=m_res[[2]]
  Lmd=m_res[[3]]
  Lmv=m_res[[4]]
  
  #NumericVector AuxMult_fun(NumericVector Y,List X_S,List X_SAge,List X_Sa, int IntervalInd,double beta, double a_g,double b_g,double c_g, double d_g,int num_particles,NumericMatrix X_a, List seeds,List age_seeds,List XSa_seeds,NumericMatrix PM, NumericVector mv, NumericVector theta_d,NumericVector age_vec,List X_Susc,NumericVector Pop_vec ){
  log_aux<-AuxMult_fun(Y[IntervalInd+1,],X_S,X_SAge,X_Sa,IntervalInd,beta,ag,bg,cg,dg,num_particles,X_a,seeds,age_seeds,XSa_seeds,PM,mv,md,age_vec,X_Susc,Pop_vec)
  log_gt<-log_g+logw+log_aux
  maxgt=max(log_gt)
  temp=log_gt-maxgt
  gt=exp(temp)
  g=gt/sum(gt)
  log_g<-log(g)
  
  ESS<-c(ESS,1/sum(g^2)) 
  
  #List resampling(List X_S,List X_SAge,List X_Sa,List X_a, NumericVector W, int num_particles,NumericVector indi,List X_Susc){
  flagr=0
  #if(ESS[IntervalInd]<(0.8*num_particles)){
    res=resampling(X_S,X_SAge,X_Sa,X_a,g,num_particles,indi,X_Susc)
    X_S=res[[1]]
    X_a=res[[2]]
    X_SAge=res[[3]]
    X_Sa=res[[4]]
    X_Susc=res[[5]]
    indv=res[[6]]+1
    log_g<-rep(0,num_particles)
    flagr=1
  #}

  #List RegeneratingPar_fun(NumericVector mv, NumericVector md,NumericVector theta_v,NumericVector theta_d, NumericVector indv, double h_m, NumericVector w){
  X_d<-rbind(X_d,theta_d)
  X_v<-rbind(X_v,theta_v)
  par_res=list()
  par_res=RegeneratingPar_fun(Lmv,Lmd,Ltheta_v,Ltheta_d,indv,h_m,w)
  Ltheta_d=par_res[[1]]
  Ltheta_v=par_res[[2]]
  
  theta_d=exp(Ltheta_d);
  theta_v=exp(Ltheta_v);
  
  print("propagation\n")
  #propagation<-function(X_a,X_S,X_SAge,X_Sa,PM,num_particles,IntervalInd,a_g,b_g,seeds,age_seeds,theta_d, age_vec,X_Susc,Pop_vec){
  res=propagation(X_a,X_S,X_SAge,X_Sa,PM,num_particles,IntervalInd+1,ag,bg,seeds,age_seeds,theta_d,age_vec,X_Susc,Pop_vec)
  
  X_S=res[[1]]
  X_a=res[[2]]
  X_SAge=res[[3]]
  X_Sa=res[[4]]
  X_Susc=res[[5]]
  #printf("weights, %f\n",max(X_a[IntervalInd,]))
  
  #NumericVector weights_fun(int Y,List X_S,int IntervalInd,double beta, double a_g,double b_g,int num_particles,NumericMatrix X_a,double v_nb)
  #Weights<-rbind(Weights,w)
  if(flagr==1){
    log_aux<-log_aux[indv]
  }
  #NumericVector weights_fun(NumericVector Y,List X_S,int IntervalInd,double beta, double a_g,double b_g,int num_particles,NumericMatrix X_a, List seeds,NumericVector theta_v,int Num_Agroups){
  logl<-weights_fun(Y[IntervalInd+1,],X_Sa,IntervalInd,beta,cg,dg,num_particles,XSa_seeds,theta_v,length(age_groups))
  logwt<-logl-log_aux
  maxwt=max(logwt)
  temp=logwt-maxwt
  wt=exp(temp)
  w=wt/sum(wt)
  logw<-log(w)
  ESSw<-c(ESSw,1/sum(w^2))
  print(ESSw[IntervalInd+1])
  # if(ESSw[(IntervalInd+1)]<(0.8*num_particles)){
  res=resampling(X_S,X_SAge,X_Sa,X_a,w,num_particles,indi,X_Susc)
  ind_sampl=res[[6]]+1
  X_dw<-rbind(X_dw,theta_d[ind_sampl])
  X_vw<-rbind(X_vw,theta_v[ind_sampl])
  Particles[[(IntervalInd+1)]]<-res  
}

X_d<-rbind(X_d,theta_d)
X_v<-rbind(X_v,theta_v)  

 
plot(Y[,1],type="o",col="red",ylim=c(min(Y),max(Y)),lwd=3)
lines(Y[,2],type="o",col="blue",lwd=3)
#lines(Y[,3],type="o",col="pink",lwd=3)
#lines(Y[,4],type="o",col="brown",lwd=3)

sum(rowSums(Y))
colSums(Y)

alpha_est1<-matrix(0,num_intervals-1,num_particles)
alpha_est2<-matrix(0,num_intervals-1,num_particles)
#alpha_est3<-matrix(0,num_intervals-1,num_particles)
#alpha_est4<-matrix(0,num_intervals-1,num_particles)

for(i in 1:(num_intervals-1)){
  if((i+4)<=(Tend-1)){
    temp=Particles[[(i+4)]]
    X_S=temp[[1]]
    X_a=temp[[2]]
  }else{
    temp=Particles[[Tend-1]]
    X_S=temp[[1]]
    X_a=temp[[2]]
  }
  for(j in 1:num_particles){
    Xa=X_a[[j]]
    alpha_est1[i,j]=Xa[1,i]
    alpha_est2[i,j]=Xa[2,i]
 #   alpha_est3[i,j]=Xa[3,i]
  #  alpha_est4[i,j]=Xa[4,i]
  }
}


Saved_RdwH<-list(ESS=ESS,ESSw=ESSw,X_d=X_d,X_v=X_v,X_dw=X_dw,X_vw=X_vw,alpha_est1=alpha_est1,alpha_est2=alpha_est2)
save(Saved_RdwH, file = "DataSc4_RG4H.RData")
######################################################
EstY1<-matrix(0,num_intervals-1,num_particles)
EstY2<-matrix(0,num_intervals-1,num_particles)
#EstY3<-matrix(0,num_intervals-1,num_particles)
#EstY4<-matrix(0,num_intervals-1,num_particles)
for(i in 1:(num_intervals-1)){
  if((i+4)<=(Tend-1)){
    temp=Particles[[(i+4)]]
    X_Sa=temp[[4]]
  }else{
    temp=Particles[[Tend-1]]
    X_Sa=temp[[4]]
  }
  for(j in 1:num_particles){
    XSa=X_Sa[[j]]
    EstY1[i,j]=length(XSa[[1]][[i]])
    EstY2[i,j]=length(XSa[[2]][[i]])
 #   EstY3[i,j]=length(XSa[[3]][[i]])
  #  EstY4[i,j]=length(XSa[[4]][[i]])
  }
}



EstWeeklyHCases<-matrix(0,num_intervals-1,num_particles)
#EstDailyHCases<-matrix(0,112,num_particles)
for(i in 1:(num_intervals-1)){
  if((i+4)<=(Tend-1)){
    temp=Particles[[(i+4)]]
    X_S=temp[[1]]
    X_a=temp[[2]]
  }else{
    temp=Particles[[Tend-1]]
    X_S=temp[[1]]
    X_a=temp[[2]]
  }
  for(j in 1:num_particles){
    EstWeeklyHCases[i,j]=length(X_S[[j]][[i]])
  }
  
}

Saved_InfH<-list(EstY1=EstY1,EstY2=EstY2,EstWeeklyHCases=EstWeeklyHCases)
save(Saved_InfH, file = "DataSc4_InfG4H.RData")

time.points<-seq(42,161,0.5)
time.points<-time.points[-1]
time.points<-time.points[-length(time.points)]
#lambda_est<-lambdaEst_fun(time.points,num_particles,seeds,Tend, a_g,b_g,X_S,X_a)
lambda_est1<-matrix(0,length(time.points),num_particles)
lambda_est2<-matrix(0,length(time.points),num_particles)
#lambda_est3<-matrix(0,length(time.points),num_particles)
#lambda_est4<-matrix(0,length(time.points),num_particles)
#lambda_est<-matrix(0,length(time.points),num_particles)
St1<-matrix(0,length(time.points),num_particles)
St2<-matrix(0,length(time.points),num_particles)
#St3<-matrix(0,length(time.points),num_particles)
#St4<-matrix(0,length(time.points),num_particles)
count_line=0
for(t in time.points){
  print(t)
  count_line=count_line+1
  indt<-find_interval(t,Tend-1,2*hist_end)
  if((indt+4+1)<=(Tend-1)){
    temp=Particles[[(indt+4+1)]]
    X_S=temp[[1]]
    X_SAge=temp[[3]]
    X_a=temp[[2]]
  }else{
    temp=Particles[[Tend-1]]
    X_S=temp[[1]]
    X_SAge=temp[[3]]
    X_a=temp[[2]]
  }
  for(j in 1:num_particles){
    #double lambda_fun(List T_N,List T_NAge,NumericVector alpha, double a_g,double b_g, double t, List seeds,List age_seeds,double Tend,NumericMatrix PM, int a_PM){
    res=lambda_fun(X_S[[j]],X_SAge[[j]],X_Sa[[j]][[1]],X_a[[j]],ag,bg,t,seeds,age_seeds,XSa_seeds,Tend-1,PM,0,Susc_vec,Pop_vec)
    lambda_est1[count_line,j]=res[[1]]
    St1[count_line,j]=res[[2]]
    
    res=lambda_fun(X_S[[j]],X_SAge[[j]],X_Sa[[j]][[2]],X_a[[j]],ag,bg,t,seeds,age_seeds,XSa_seeds,Tend-1,PM,1,Susc_vec,Pop_vec)
    lambda_est2[count_line,j]=res[[1]]
    St2[count_line,j]=res[[2]]
    
    # res=lambda_fun(X_S[[j]],X_SAge[[j]],X_Sa[[j]][[3]],X_a[[j]],ag,bg,t,seeds,age_seeds,XSa_seeds,Tend-1,PM,2,Susc_vec,Pop_vec)
    # lambda_est3[count_line,j]=res[[1]]
    # St3[count_line,j]=res[[2]]
    # 
    # res=lambda_fun(X_S[[j]],X_SAge[[j]],X_Sa[[j]][[4]],X_a[[j]],ag,bg,t,seeds,age_seeds,XSa_seeds,Tend-1,PM,3,Susc_vec,Pop_vec)
    # lambda_est4[count_line,j]=res[[1]]
    # St4[count_line,j]=res[[2]]    
  }
}
Saved_IntH<-list(lambda_est1=lambda_est1,lambda_est2=lambda_est2,St1=St1,St2=St2)
save(Saved_IntH, file = "DataSc4_IntG4H.RData")


################ finding the reproduction numbers #################################
R1<-matrix(0,length(time.points),num_particles)
R2<-matrix(0,length(time.points),num_particles)
count_line=0
for(t in time.points){
  print(t)
  count_line=count_line+1
  indt<-find_interval(t,Tend-1,2*hist_end)
  if((indt+4+1)<=(Tend-1)){
    temp=Particles[[(indt+4+1)]]
    X_S=temp[[1]]
    X_SAge=temp[[3]]
    X_a=temp[[2]]
  }else{
    temp=Particles[[Tend-1]]
    X_S=temp[[1]]
    X_SAge=temp[[3]]
    X_a=temp[[2]]
  }
  for(j in 1:num_particles){
    St_ag<-c( St1[count_line,j], St2[count_line,j])
    tempR<-c()
    for(i in 1:length(age_groups)){
      Xa=X_a[[j]]
      tempR<-c(tempR, sum( ( (St_ag/Pop_vec)*Xa[,indt+1])*PM[,i] ) )
    }
    R1[count_line,j]=tempR[1]
    R2[count_line,j]=tempR[2]
  }
}

Saved_RAgH<-list(R1=R1,R2=R2)
save(Saved_RAgH, file = "DataSc4_RAgH.RData")



