clarke_compute_pA <- function(expressiondataA, expressiondataAB, pvalA, pvalAB) {
#	R code to implement method for estimation of proportion of
#	one sample type (A) in a two sample mixture (A and B), by method 
#	of clarke et al. 
#
#	Corrected by Gerald Quon, 2011/11/11
#
#	assumes input file of pre-processed gene expression data is tab 
#	delimited text file with header row (sample labels; not input) and 
#	transcript labels in first column; a separate file of the same format
#	contains detection p-values
#	file must contain expression from mixture sample(s) (A and B) and 
#	expression from samples of only one type (A); this method estimates
#	the proportion of sample A in the mixture samples of A and B
#	code designed for use with example data files 
#		E-GEOD-5130-processed-data-titration-2-exp.txt
#		E-GEOD-5130-processed-data-titration-2-pval.txt
#exprs <- as.matrix(read.table("E-GEOD-5130-processed-data-titration-2-exp.txt", header = TRUE, sep = "\t", row.names = 1, as.is = TRUE))
	
	
exprs = cbind(expressiondataA,expressiondataAB);


## p=0.8
# select only expressed genes
idt = intersect(which(pvalsA < 0.05), which(pvalsAB < 0.05));
exprs = exprs[idt,];
	
## can we determine alpha from mean(Ri) or median(Ri)?
# note: if you have more than one sample at each value of p, you can use
# the average across samples in the estimation process 

# guess an initial range for alpha, the parameter for the data transformation
KK=2000; #number of alpha to try
alpha<-seq(from=0.001,to=10,length.out=KK)

#exprs-- 2 column dataset, first column is pure sample A, second column is mixture AB, trying to compute pA, fraction of AB that is A
alpha_data = matrix(NA, KK, 6);

#column names of alpha will be alpha, mean transformed ratio (meantRi), median transformed ratio (mediantRi), minimum transformed ratio (mintRi), arc length s (s), and a statistic we need to compute radius of curvature, 
# || f''(...)||   (radiusstat)
colnames(alpha_data) = c('alpha','meantRi','mediantRi','mintRi','s','radiusstat');
alpha_data = as.data.frame(alpha_data);
alpha_data$alpha = alpha;

for(ii in 1:KK){
	tE<-log(1+ alpha[ii]*exprs)
	tE_A<-tE[,1]
	tE_AB<-tE[,2]
	tRi<-tE_AB/tE_A
	alpha_data$mintRi[ii] = min(tRi)

	alpha_data$meantRi[ii] = mean(tRi)
	alpha_data$mediantRi[ii] = median(tRi)
}

#we need to re-scale mean(tRi) and median(tRi)
alpha_data$meantRi = ((max(alpha_data$alpha)-min(alpha_data$alpha)) * (alpha_data$meantRi-min(alpha_data$meantRi))) / (max(alpha_data$meantRi) - min(alpha_data$meantRi));
alpha_data$mediantRi = ((max(alpha_data$alpha)-min(alpha_data$alpha)) * (alpha_data$mediantRi-min(alpha_data$mediantRi))) / (max(alpha_data$mediantRi) - min(alpha_data$mediantRi));

#########
# we will assume we are going to use meantRi, not mediantRi, to compute pA
#compute s (arc length) values.  needs to be done after rescaling.  
#arc length for 1st alpha is 0
alpha_data$s[1] = 0; 

#we need to compute derivative(mean(tR(\alpha_i))) as in paper
approx_derivatives = (alpha_data$meantRi[2:KK] - alpha_data$meantRi[1:(KK-1)])/(alpha_data$alpha[2:KK] - alpha_data$alpha[1:(KK-1)]);
for(ii in 2:(KK-1)){
	alpha_data$s[ii] = sum(sqrt(1 + (approx_derivatives[2:ii]^2)) * (1/(alpha_data$alpha[2:ii] - alpha_data$alpha[1:(ii-1)])));
}

#compute summary statistic for radius of curvature, || f''(...) ||
#has to wait till we are done computing arc lengths s, because of certain approximations
#also has to start at ii=2 because the center difference approximation needs points before and after the current point
#set first radiusstat to negative infinity so it is never chosen as the best value
alpha_data$radiusstat[1]=-Inf;
for(ii in 2:KK){
	alpha_data$radiusstat[ii] = abs(    (alpha_data$meantRi[ii+1]- 2*alpha_data$meantRi[ii] + alpha_data$meantRi[ii-1])/((alpha_data$s[ii]-alpha_data$s[ii-1])^2)              ); 
}

alpha_best_ix = which.max(alpha_data$radiusstat);

pA = alpha_data$mintRi[alpha_best_ix];
pA

##diagnostics
## make plot of mean Ri across values of alpha
#plot(1:length(alpha_data$alpha),alpha_data$meantRi)
##points(1:length(alpha_data$alpha),alpha_data$mediantRi,col="red")
##legend("bottomright",c("mean","median"),pch=1,col=c("black","red"))
## add vertical line to plot to indicate minimum radius of curvature
#lines(c(alpha_best_ix,alpha_best_ix),c(0,10))
	
	

}

