##############################################################################################################################
##R CODE FOR A PERIOD LIFE TABLE, REPRODUCED FROM ONLINE MATERIALS FOR A 2006 FORMAL DEMOGRAPHY WORKSHOP AT STANFORD UNIVERSITY
##THE ORIGINAL POSTING OF THIS CODE IS AVAILABLE AT: 	http://www.stanford.edu/group/heeh/cgi-bin/web/node/75 
##				https://web.archive.org/web/20170209225451/http://web.stanford.edu/group/heeh/cgi-bin/web/node/75
##
##NOTE: I (EDDIE HUNSINGER) REPRODUCED THE CODE HERE TO PROVIDE AN IMMEDIATE LINK TO INPUT DATA FOR QUICK REVIEW BY POTENTIAL USERS, BUT I DID NOT AUTHOR THIS CODE
##FEBRUARY 2011 (LAST UPDATED OCTOBER 2022)
##edyhsgr@gmail.com
##
##A PYTHON (USING JUPYTER NOTEBOOK) VERSION OF THE CODE, I ADDED IN OCTOBER 2022, IS AVAILABLE AT:  
##
##IF YOU WOULD LIKE TO USE, SHARE OR REPRODUCE THIS CODE, BE SURE TO CITE THE SOURCE
##
##EXAMPLE DATA IS LINKED, SO YOU SHOULD BE ABLE TO SIMPLY COPY ALL AND PASTE INTO R
##
##THERE IS NO WARRANTY FOR THIS CODE
##############################################################################################################################

##############################################################################################################################
#STEP 1: Read in and review the population and death data
##############################################################################################################################

females<-read.table(file="https://github.com/AppliedDemogToolbox/StanfordCourseLifeTable/raw/master/StanfordCourseMortalityData.csv",header=TRUE,sep=",")
females

##############################################################################################################################
#STEP 2: Read in or create the fundamental pieces of the life table (age groupings, deaths by age, population by age ->death rates by age
##############################################################################################################################

x <- c(0,1,5,10,15,20,25,35,45,55,65,75,85)
#note that R collapses a single column to a vector when it pulls out the result out of a data.frame
nDx <- females$Death.Count   #other syntax which produces the same result: females[[3]], females[,3], 
nKx <- females$Population
nMx <- nDx / nKx

##############################################################################################################################
#STEP 3: Read in the period life table function
##############################################################################################################################

life.table <- function( x, nMx){
## simple lifetable using Keyfitz and Flieger separation factors and 
## exponential tail of death distribution (to close out life table)
b0 <- 0.07;   b1<- 1.7;      
nmax <- length(x)
#nMx = nDx/nKx   
n <- c(diff(x),999)          		#width of the intervals
nax <- n / 2;		            # default to .5 of interval
nax[1] <- b0 + b1 *nMx[1]    		# from Keyfitz & Flieger(1968)
nax[2] <- ifelse(n[2]==4,1.5,nax[2])  ;           #EddieH note: modified to ifelse to include complete (single-year) tables - October 2022  
nax[nmax] <- 1/nMx[nmax] 	  	# e_x at open age interval
nqx <- (n*nMx) / (1 + (n-nax)*nMx)
nqx<-ifelse( nqx > 1, 1, nqx);	# necessary for high nMx
nqx[nmax] <- 1.0
lx <- c(1,cumprod(1-nqx)) ;  		# survivorship lx
lx <- lx[1:length(nMx)]
ndx <- lx * nqx ;
#Edited by EddieH (below line is my edit of this line) - February 2018 #nLx <- n*lx - nax*ndx;       		# equivalent to n*l(x+n) + (n-nax)*ndx
nLx <- n*lx - (n-nax)*ndx;       		# equivalent to n*l(x+n) + nax*ndx
nLx[nmax] <- lx[nmax]*nax[nmax]
Tx <- rev(cumsum(rev(nLx)))
ex <- ifelse( lx[1:nmax] > 0, Tx/lx[1:nmax] , NA);
lt <- data.frame(x=x,nax=nax,nmx=nMx,nqx=nqx,lx=lx,ndx=ndx,nLx=nLx,Tx=Tx,ex=ex)
return(lt)
}

##############################################################################################################################
#STEP 4: Apply the function to the data, and review the created life table
##############################################################################################################################

females.life.table<-life.table(x,nMx)
females.life.table

#write.table(###, file="G:/###/###.csv", sep=",")

