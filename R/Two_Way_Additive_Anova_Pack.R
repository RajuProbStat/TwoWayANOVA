
myfunc1<-function(a,b,level,n,mean,var)
{ 
  B<-10000
  mean_alpha<-rep(NA,a)           #level 1 means
  mean_beta<-rep(NA,b)            #level 2 means
  mean_boot<-array(NA,dim=c(a,b)) # means for bootstrap samples
  var_boot<-array(NA,dim=c(a,b))  # means for bootstrap samples
  alpha<-rep(NA,a)
  mean_alpha_boot<-rep(NA,a)      # level 1 bootstrap means
  mean_beta_boot<-rep(NA,b)       # level 2 bootstrap means
  test1_boot<-rep(NA,B)
  test2_boot<-rep(NA,B)
  for(i in 1:a)
  {
    alpha[i]<-mean(mean[i,])
  }
  for(i in 1:a)
  {
    mean_alpha[i]<-mean(mean[i,])-mean(alpha)
  }
  for(j in 1:b)
  {
    mean_beta[j]<-mean(mean[,j])-mean(alpha)
  }
  var1<-rep(0,a-1)
  var2<-rep(0,b-1)
  T1<-rep(NA,a-1)
  T2<-rep(NA,b-1)
  for(i in 1:a-1)
  {
    var1_0<-0
    for(j in 1:b)
    {
      var1[i]<-var1_0+(var[i+1,j]/n[i+1,j])+(var[i,j]/n[i,j])
      var1_0<-var1[i]
    }
    T1[i]<-(mean_alpha[i+1]-mean_alpha[i])/sqrt(var1[i]/b^2)
  }
  for(j in 1:b-1)
  {
    var2_0<-0
    for(i in 1:a)
    {
      var2[j]<-var2_0+(var[i,j+1]/n[i,j+1])+(var[i,j]/n[i,j])
      var2_0<-var2[j]
    }
    T2[j]<-(mean_beta[j+1]-mean_beta[j])/sqrt(var2[j]/a^2)
  }
  test1_ob<-min(c(min(c(T1)),min(c(T2))))
  test2_ob<-min(c(max(c(T1)),max(c(T2))))
  #for bootstrap samples
  for(l in 1:B)
  {
    for(i in 1:a)
    {
      for(j in 1:b)
      {
        data_boot<-rnorm(n[i,j],0,sqrt(var[i,j]))
        mean_boot[i,j]<-mean(data_boot)
        var_boot[i,j]<-var(data_boot)
      }
      mean_alpha_boot[i]<-mean(mean_boot[i,])
    }
    for(j in 1:b)
    {
      mean_beta_boot[j]<-mean(mean_boot[,j])
    }
    var1_boot<-rep(0,a-1)
    var2_boot<-rep(0,b-1)
    T1_boot<-rep(NA,a-1)
    T2_boot<-rep(NA,b-1)
    for(i in 1:a-1)
    {
      var1_0<-0
      for(j in 1:b)
      {
        var1_boot[i]<-var1_0+(var_boot[i+1,j]/n[i+1,j])+(var_boot[i,j]/n[i,j])
        var1_0<-var1_boot[i]
      }
      T1_boot[i]<-(mean_alpha_boot[i+1]-mean_alpha_boot[i])/sqrt(var1_boot[i]/b^2)
    }
    for(j in 1:b-1)
    {
      var2_0<-0
      for(i in 1:a)
      {
        var2_boot[j]<-var2_0+(var_boot[i,j+1]/n[i,j+1])+(var_boot[i,j]/n[i,j])
        var2_0<-var2_boot[j]
      }
      T2_boot[j]<-(mean_beta_boot[j+1]-mean_beta_boot[j])/sqrt(var2_boot[j]/a^2)
    }
    test1_boot[l]<-min(c(min(c(T1_boot)),min(c(T2_boot))))
    test2_boot[l]<-min(c(max(c(T1_boot)),max(c(T2_boot))))
  }
  critical_test1<-quantile(test1_boot,1-level)
  critical_test2<-quantile(test2_boot,1-level)
  cat("critical point of T-min-min-min, observed value of T-min-min-min, critical point of T-min-max-max and observed value of T-min-max-max are\n")
  return(c(critical_test1,test1_ob,critical_test2,test2_ob))
}


