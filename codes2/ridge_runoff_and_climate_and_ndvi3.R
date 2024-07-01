#
rm(list = ls())
library(MASS)

varset=c("P","NDVI","Es","q")
num_vars=length(varset)
num_stations=7

startyear=1950
endyear=2022
num_years=endyear-startyear+1

xset1=matrix(data=NA,nrow = num_years,ncol = num_vars)
rset=matrix(data = NA,nrow=num_stations,ncol=2)

sensitivity_set=matrix(data=NA,nrow = num_vars-1,ncol = num_stations)
# read data
for (s in 1:num_stations) {
  s1<-as.character(s)
  xset=read.table(paste("G:/variation_and_attribution_for_runoff_and_sediment_flux/codes2/input_data/annual_x_change_relative_to_long_term_mean_",s1,".txt", sep=""))
  xset=xset[,c(1,3,4,5)]
  colnames(xset)=c("P","NDVI","Es","q")
  xset0=as.data.frame(xset)
# build model  
  b<-lm.ridge(q ~ P+NDVI+Es, xset0, lambda = seq(0,1,0.00001))
  index1=which.min(b$GCV)
  opt_coef<-coef(b)[index1,]
  opt_coef1<-as.matrix(opt_coef)
  xset3<-as.matrix(xset)
# calculate fitted values 
  r1=dim(xset3)[1]
  q_esti=matrix(data=NA,nrow = r1,ncol=num_stations)
  q_esti[,s]=xset3[,1:num_vars-1] %*% opt_coef1[2:num_vars]+opt_coef1[1]
  
 # calculate R2 
  corr=cor.test(q_esti[,s],xset3[,num_vars])
  r1=corr[["estimate"]]
  rset[s,1]=r1^2
  rset[s,2]=corr[["p.value"]]
  
  # sensitivity
  sensitivity_set[,s]=as.matrix(opt_coef1[2:num_vars,1])
  
  
}


# output the R2 between simulated and observed runoff
write.table(rset,paste("G:/variation_and_attribution_for_runoff_and_sediment_flux/codes2/output/ridge_regression_r_square",".txt",sep=""),
            col.names = FALSE,row.names = FALSE,quote = FALSE,sep = "\t")
# output the sensitivities
write.table(sensitivity_set,paste("G:/variation_and_attribution_for_runoff_and_sediment_flux/codes2/output/ridge_regression_sensitivity",".txt",sep=""),
            col.names = FALSE,row.names = FALSE,quote = FALSE,sep = "\t")
