rm(list = ls())
library(plspm)
library(ggplot2)
library(export)
library(rgl)
library(rvg)

varset=c("pre","ta","growing_season_noaa_cdr_ndvi","gleam_es","q")
label_set=c("a","b","c","d","e","f","g")

num_vars=length(varset)
num_stations=7

cindex2=c(11,8,9,10,7,12,2)  # # for Q  Adjust the order of sites displayed sequentially in the subplots.
cindex3=c(7,4,5,6,1,3,2)#Adjust the order of sites displayed sequentially in the subplots.


startyear=1950
endyear=2022
num_years=endyear-startyear+1

xset1=matrix(data=NA,nrow = num_years,ncol = num_vars)
rset=matrix(data=NA,nrow=num_stations,ncol=1)

#plotting parameters
plot.layout=c(3,3)
windows(width=8,height=7,xpos = 180,ypos = 25)

par(mfrow=plot.layout)

par(mar=c(1,1,1,1),xpd=TRUE)
par(omi=c(0,0,0,0))

arrpos=matrix(data=0.5,nrow = num_vars,ncol = num_vars)

arrpos[3,5]=0.45
arrpos[2,5]=0.5
arrpos[1,2]=0.5
arrpos[1,3]=0.45

windowsFonts(TNM = windowsFont("Arial"))
par(family="TNM")
par(font=2)
#read data
for (s in 1:num_stations) {
  s1=cindex3[s]
  s2<-as.character(s1)
  xset=read.table(paste("G:/variation_and_attribution_for_runoff_and_sediment_flux/codes2/input_data/annual_x_change_relative_to_long_term_mean_",s2,".txt", sep=""))
  xset=xset[,c(1,2,3,4,5)]
  colnames(xset)=c("P","Ta","NDVI","Es","q")
  xset2=xset

  # set path  
P <-c(0,0,0,0,0)
Ta <- c(0,0,0,0,0)
NDVI <- c(1,1,0,0,0)
Es <- c(1,1,1,0,0)
Q <- c(1,1,1,1,0) 


set_path<-rbind(P,Ta,NDVI,Es,Q)

set_blocks <- list(1,2,3,4,5)
set_modes <- rep("A", 5) 

options(digits = 4)

q_sem=plspm(xset2,set_path,set_blocks,modes=set_modes,
            scaled = TRUE,boot.val = TRUE)
summary(q_sem)
q_sem$path_coefs=round(q_sem$path_coefs,2)
if (s==4){
}
rset[s,1]=q_sem$inner_summary$R2[num_vars]
abs_path_coefs=abs(q_sem$path_coefs)*7
index_coef=which(abs(q_sem$path_coefs)<0.2 & abs(q_sem$path_coefs)>0)
abs_path_coefs[index_coef]=1

##resize arrows
abs_path_coefs1=matrix(data=NA,nrow = num_vars,ncol = num_vars)
for (i1 in 1:num_vars) {
  for (j1 in 1:num_vars) {
    if (abs_path_coefs[i1,j1]==1){
      abs_path_coefs1[num_vars-i1+1,num_vars-j1+1]=0.1
    }else{
    abs_path_coefs1[num_vars-i1+1,num_vars-j1+1]=0.1+abs_path_coefs[i1,j1]*0.04
    }
    
  }
  
}
theme(text = element_text(family = "A",face = "bold"))

# plotting
p1<-plot(q_sem, what = "inner",colpos = "#0571b0",colneg = "#ca0020",
       txt.col="black",lcol = "NA",box.col=c("#f4a582","#92c5de","#92c5de","#92c5de","#92c5de"),
     arr.lwd=abs_path_coefs,cex.txt=1.2,box.cex = 1.4,box.size = 0.1,
     box.lwd=0.2,arr.pos = arrpos,arr.width=abs_path_coefs1,arr.length=abs_path_coefs1,
     arr.tcol="black",txt.font=6)

row=round(s/3,0)+1
col=s%%3

}
par(family="serif")

#      Output the R-squared value between simulated and observed values for each path
write.table(rset,paste("G:/variation_and_attribution_for_runoff_and_sediment_flux/codes2/output/pls-sem_r_square",".txt",sep=""),
            col.names = FALSE,row.names = FALSE,quote = FALSE,sep = "\t")

#     export image
graph2jpg(file='G:/variation_and_attribution_for_runoff_and_sediment_flux/codes2/output/pls_sem_results.jpeg',
          dpi=300,bg="white",font="Times New Roman")


