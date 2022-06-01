##' est_power2
##' 
##' A function to estitamete the power for differential expression analysis of RNA-seq data.
##' With 3 df precision.
##'
##' 
##' @param n Numer of samples.
##' @param alpha alpha level.
##' @return Estimate power
##' @inheritParams sample_size
##' @export
##' @examples n<-63;rho<-2;lambda0<-5;phi0<-0.5;f<-0.01
##' est_power(n=n, rho=rho, lambda0=lambda0, phi0=phi0,f=f)
est_power2<-function(n, w=1, k=1,rho=2, lambda0=5, phi0=1,alpha=0.05,f,m=20000,m1=200){
	if (!missing(f)) {#FDR power
		power_fdr_1000<-est_power_root_fdr2(power=1,n=n, w=w,k=k, rho=rho, lambda0=lambda0, phi0=phi0,m=m,m1=m1,fdr=f)
		if (power_fdr_1000<0.001) {
			return(0)
		} else {
			power_fdr<-uniroot.integer(f=est_power_root_fdr2,interval=c(1,1000),n=n, w=w,k=k, rho=rho, lambda0=lambda0, phi0=phi0,m=m,m1=m1,fdr=f)
			return((power_fdr$root-1)/1000)
		}
	} else {#alpha power
		power<-est_power_root(n=n, w=w,k=k, rho=rho, lambda0=lambda0, phi0=phi0, alpha=alpha)
		#beacuse est_power_root is returning (power-(1-beta)), which can be used in uniroot.integer for fdr based power estimation. 
		#If using it directly for alpha based power, need +0.8 to get correct power as default value for beta is 0.2. 
		return(power+0.8) 
	}
}

est_power_root_fdr2<-function(power,n, w,k, rho, lambda0, phi0,fdr,m,m1,...) {
	alpha_star<-m1 * power/1000*fdr/((m-m1)*(1-fdr))
#	cat(paste0(alpha_star,"\n"))
	beta<-1-power/1000
	est_power_root(n=n, w=w,k=k, rho=rho, lambda0=lambda0, phi0=phi0, alpha=alpha_star,beta=beta,...)
}
