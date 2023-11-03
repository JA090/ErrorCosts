#' Contains functions to produce Tables and Figures in the paper "An alternative to chasing an elusive alpha that minimizes ill-defined error costs:
#'  simultaneously test a hypothesis at multiple alpha levels" by J Aisbett


require(gt)
require(formattable)
require(gtExtras)
require(tibble)
chromote::default_chromote_object()

#'formatTable() is a function to format tables as in the paper and write out into current directory
#' @param Table is matrix input
#' @param name is the name used for the output word document i.e. name.docx
#' @param repeater defines groups of entries in final table
#' @param minentry is a vector containing column no. with the minimum cost of the single level tests
#' @param n -number of test alpha levels in table
formatTable <- function(Table, name, repeater, alpha, minentry,n) {
  # add column names
  colname = " "
  for (i in 1:n)
    colname = c(colname, paste0("\U03B1", "=", format(alpha[i], scientific =
                                                        F)))
  colname = c(colname, "multi-level", "optimal")
  colnames(Table) = colname
  #convert to gt object
  xx1 = gt(data = as.data.frame(Table))
  # shade columns
  xx1 <-
    tab_style(xx1, cell_fill(color = "snow2"), locations = cells_column_labels())
  xx1 <-
    tab_style(xx1, cell_fill(color = "snow2"), locations = cells_body(columns =
                                                                        1))
  cols_width(xx1, everything() ~ px(100))
  #make lowest cost alpha entry bold
  for (j in 2:length(Table[, 1])) {
    if (is.na(minentry[j - 1]) == F)
      xx1 <-
        tab_style(xx1,
                  cell_text(weight = "bold"),
                  locations = cells_body(columns = minentry[j - 1] + 1, rows = j))
  }
  #put vertical borders on multilevel and optimal columns
  xx1 = gt_add_divider(xx1,columns = c(n + 1, n + 2),color = "black", sides = "right",style = "dashed"
       )
  # fill rows which separate groups of results with white
  k = 1
  while (k < length(Table[, 1])) {
    xx1 <-
      tab_style(xx1, cell_fill(color = "white"), locations = cells_body(rows =k))
    k = k + length(repeater) + 1
  }
  #write out
  gtsave(
    xx1,
    paste0(name, ".docx"),
    vwidth = 600 + 200 * n,
    vheight = 1000,
    expand = 10
  )
}

############################################################
#' Table_1_2 is a function to compute expected total error costs of risk difference tests given two equal groups,
#' at different alpha levels and at simultaneously multiple levels.
#' Output contains Table 1 and Table 2 of the paper. Details of the computations are in that paper

#' Prevalence of true findings is modelled using both zero-mean normal distributions and dichotomous distributions of effect sizes in the research environment
#' Word documents titled Table1 and Table2 are output into current directory.

#' @param alpha # alphas in decreasing order
#' @param R1 risk in standard care/placebo
#' @param sample group sample sizes, assuming equal groups
#' @param R22 vector of true risk in treatment group in dichotomous scenarios
#' @param PP vector of prevalence in dichotomous scenarios
#' @param MM vector of mean  of effect distribution in continuous scenarios
#' @param vv vector of SD of effect distribution in continuous scenarios
#' #costs
#' @param c cost of drug per person
#' @param C av. cost of hospital stay
#'
#' n= number of alphas
#' onetwo = ratio of Type I to Type II costs
#' true effect distribution variance in current computation
#' costs = vector of expected costs for sinle level tests, multilevel and optimal
#' R2 - true effect in treatment group in current calculation with dichotomous scenario
#' P - prevalence in current computation in dichotomous scenario
#' EM -  mean of effect distribution in current computation cts scenario
#' v -  sd of effect distribution in current computation cts scenario
#' CDiff = difference between costs at adjacent test levels in multi-level tests 
#'        (assuming costs are propertional to surprisal values)
#' minentry - vector of positions of lowest cost alphas for each scenario setting
#' row - row in Table

Table_1_2 <- function(
      alpha = c(.25, .05, .001), # alphas in decreasing order
      R1 = .092, #risk, standard care/placebo
      R22 = .092 - c(.025, .05), #risk  in- dichotomous scenarios
      sample = 1000, # group sample sizes, assuming equal groups
      PP = c(.5, .1),#prevalence, dichotomous scenarios
      MM = c(0, -.02), # mean of effect distribution in continuous scenarios
      vv = c(.015, .025), #  SD of effect distribution in continuous scenarios
      cd = 707, #cost of drug per person
      C = 40000 # av. cost of hospital stay
    ) {
  # set up parameters for tests and scenarios
  n <<- length(alpha) # number of alphas
  onetwo <<- cd / C # ratio of treatment to event costs
  #set costs proportional to surprise value of a finding at each alpha level
  c = C * log2(alpha) / log2(alpha[n])
  CDiff = c[1]
  for (j in 2:n)
    CDiff = c(CDiff, c[j] - c[j - 1])  # differences between benefits for different alphas
  ## set up internal computational stuff
  Table2 = matrix(nrow = 0, ncol = n + 3)
  Table1 = matrix(nrow = 0, ncol = n + 3)
  minentry = 100 * 0
  row = 0
  costs = rep(0, n + 2) # will contain expected costs
  II <<- rep(0, n) # will contain total error costs  in cts scenarios
  
  
  #Functions to perform integration over effect sizes in a continuous research scenario
  IntegrateCts <- function(z) {
    for (i in 1:length(z)) {
      zalp <<- z[i]
      I2 <<-
        integrate(ttt2, .00001, R1 - onetwo, stop.on.error = F)$value #only get false negatives if below cut off
      II[i] <<-
        I2 + integrate(ttt1, R1 - onetwo,1, stop.on.error = F)$value # only get false positives afterward
    }
    return(II[1:length(z)])
  }
  
  #functions to compute Type I and Type II errors in continuous research scenarios
  #eh is the risk of hospitalization in treatment group
  
  ttt1 <- function(eh) {
    p2 = eh * (1 - eh) #sd of sample distribution at eh
    sd=sqrt((p1+p2)/sample)
    dnorm(eh - R1, EM, v) * C*onetwo*pnorm((zalp*sd2 -eh + R1 -onetwo)/ sd,0,1)
  }
  
  ttt2 <- function(eh) {
    p2 = eh * (1 - eh) #sd of sample distribution at eh
    factor =  R1 - eh - onetwo # cost weight
    sd=sqrt((p1+p2)/sample)
    dnorm(eh - R1, EM, v) *  C*factor * (1-pnorm((zalp*sd2 -eh + R1-onetwo)/ sd,0,1))
  }
  
  
  #######
  p1 = (1 - R1) * R1 #sample distribution variance of untreated risk
  p22=(1-R1+onetwo)*(R1-onetwo) # sample distribution of treated risk if this is R1-onetwo
  sd2=sqrt((p1+p22)/sample)
  
  # start computations for dichotomous scenario
  row = 0 #start rows for dichotomous table
  # loop for vector of true effects
  for (R2 in R22) {
    # insert scenario risk difference  into table
    Table1 = rbind(Table1, matrix(data = " ", nrow = 1,ncol = n + 3
    ))
    Table1[length(Table1[, 1]), 1] = paste0("RD =", R2 - R1)
    row = row + 1
    #prevalence loop
    for (P in PP) {
      p2 = (1 - R2) * R2 #sample distribution variance if  treated risk = R2
      factor = R1 - R2 - onetwo
      sd = sqrt((p1+p2)/sample)
      zalp= qnorm(alpha, 0, 1)
      # sum type I and type II error costs; the latter
      costs[1:n] = (1-P) * onetwo * C * alpha + P * C * factor * (1 -
                                    pnorm((zalp*sd2+factor)/sd, 0,1)) # single level tests
      costs[n + 1] = sum(costs[1:n] * CDiff) / C # multilevel
      #now find optimal by searching alpha space for minimum of single level tests
      opt2 = rep(0, 3000)
      a = .0001
      jj = 1:10000
      zalp = qnorm(jj * a, 0, 1)
      opt2 = (1-P) * onetwo * C * jj * a + P * C * factor * (1 - pnorm(zalp*sd2,-factor,sd))
      costs[n + 2] = min(opt2)
      optim = round(a * which.min(opt2), 2)
      if (optim == 0)
        optim = formattable(a * which.min(opt2),
                            format = "f",
                            digits = 4)
      # add entry into tables
      rowname = paste0("P=", P) # paste prevalence into table
      Table1 = rbind(Table1, c(rowname, paste0(
        formattable(costs, format = "f", digits = 1)
      )))
      row = row + 1
      Table1[row, n + 3] = paste0(round(costs[length(costs)], 1), " (", optim, ")")
      minentry[row - 1] = which.min(costs[1:n]) # store where minimum cost as this will be bolded in final table
    } # end prevalence loop
  } #end true effect loop
  #output results
  formatTable(Table1, "Table1", R22, alpha, minentry,n)
  
  
  #  computations for cts scenario
  row = 0
  minentry = 100 * 0
  #loop for mean of true effect distribution
  for (EM in MM) {
    row = row + 1
    # insert scenario effect size distribution mean into table
    Table2 = rbind(Table2, matrix(
      data = " ",
      nrow = 1,
      ncol = n + 3
    ))
    Table2[length(Table2[, 1]), 1] = paste0("m = ", EM)
    ## loop over variance of true effect distribution
    for (v in vv) {
      # compute  integrals of costs over effect size distribution in research scenario
      # costs for single alpha tests
      costs[1:n] = IntegrateCts(qnorm(alpha, 0, 1)) # costs
      # costs for multilevel test
      costs[n + 1] = sum(costs[1:n] * CDiff) / C
      #find optimal
      opt = rep(0, 2000)
      j = 1:2000
      a = .0005
      opt = IntegrateCts(qnorm(a * j, 0, 1))
      costs[n + 2] = opt[which.min(opt)]
      optim = round(a * which.min(opt), 2)
      if (optim == 0)
        optim = signif(which.min(opt) * a)
      # write results into table
      rowname1 = paste0("sd = ", round(v, 2))# add sd of scenario effect distribution into table
      Table2 = rbind(Table2, c(rowname1, paste0(
        formattable(costs, format = "f", digits = 1)
      )))
      Table2[row + 1, n + 3] = paste0(round(costs[length(costs)], 1), " (", optim, ")")
      minentry[row] = which.min(costs[1:n]) # store minimum costs from single level tests
      row = row + 1
    }  # end of variance loop
  } # end of loop over research distribution mean
  
  #output and output Table
  formatTable(Table2, "Table2", vv, alpha, minentry,n)
} # end of function Table_1_2

  ############################################################
  #' Table_3 is a function to compute averaged expected total error costs of superiority t-tests given two equal groups,
  #' at different alpha levels and at simultaneously multiple levels.
  #' The averages are over the distribution of anticipated standardized effect size of different research teams
  #' in the research scenario, which is assumed to be normal.
  #' Output contains Table 3 of the paper. Details of the computations are in that paper
  #' Prevalence of true findings is modelled with normal distributions of standardized effect sizes in the research environment
  #' Word document Table3 is output into current directory.
  
  #' @param alpha alphas in decreasing order
  #' @param aa vector of means of true effect size distributions 
  #' @param sd1 standard deviation of true effect size distributions 
  #' @param D mean of distrution of anticipated effect sizes
  #' @param sd2 standard deviation of anticipated effect size distributions 
  #' @param cR vector of cost ratios between Type I and Type II errors.   
  #'
  #' n= group sample size (assuming equal groups)
  #' a = current true distribution mean
  #' cRatio = current Type I to Type II error cost ratio
  #' CDiff = difference between costs at adjacent test levels in multi-level tests 
  #'        (assuming costs are propertional to surprisal values)
  #' minentry - vector of positions of lowest cost alphas for each scenario setting
  #' row - row in Table
  

 Table_3<- function(alpha=c(.25,.025),aa=c(-.1,0,.4),sig1=.2,D=.4,sig2=.1,cR=c(10,4,1))
  {
    
    Table = matrix(nrow = 0, ncol = length(alpha) + 3)
    #########################
    # Type1 (resp Type2) is function to compute Type I (resp. Type II) errors 
    # given effect size e, group size n and true distribution mean a
    
    Type1<-function(e,n){ 
      #prob assumed d*resulting Type I error
      dnorm(e,a,sig1)*(1-pt(-e*sqrt(n/2)+qt(1-alp,2*n-2),2*n-2))
    }
    Type2<-function(e,n){
      dnorm(e,a,sig1)*(pt(-e*sqrt(n/2)+qt(1-alp,2*n-2),2*n-2))
    }
    ########################
    # Type1a (resp Type2a) is a function to compute average Type1 (resp. Type 2) error
    # rates over the distribution of anticipated effect sizes which affects group siZe n
   Type1a<-function(x,n){ 
      n<-as.integer(2*(zj/x)^2)
      dnorm(x,D,sig2)*Integrate(Type1,n,a-4*sig1,0)
    }
    
    Type2a<-function(x,n){
      #prob assumed x*resulting Type II error
      n<-as.integer(2*(zj/x)^2)
      dnorm(x,D,sig2)*Integrate(Type2,n,0,a+4*sig1)
    }
    
    Integrate <- function(f, n,lower, upper)
    {ff<-function(e){f(e,n)}
      integrate(ff,lower,upper,stop.on.error=F)$value
    }
    ##############################
    # start computations
    row=0
    minentry=100*0
    # compute cost weights for multi-level case
    CDiff=rep(log2(alpha[1])/log2(alpha[length(alpha)]),length(alpha))
    for (j in 2:length(alpha)) CDiff[j]=(log2(alpha[j])-log2(alpha[j-1]))/log2(alpha[length(alpha)])
    
    TC1=rep(0, length(alpha))
    TC2=TC1
    # all research teams are assumed to set sample sizes to achieve 80% power at most stringent alpha
    zj<-qnorm(1-alpha[length(alpha)],0,1)+qnorm(.8,0,1)
    
    # start loop for different cost ratios
    for (cRatio in cR)
      {Table = rbind(Table, matrix(data = " ",nrow = 1, ncol = length(alpha) + 3))
      Table[length(Table[, 1]), 1] = paste0("cost ratio = ", cRatio)
      row = row + 1
     # start loop for different means of the true effect size distribution 
      for (a in aa)
        { rowname= "true mean = M" # rowname is entry in column 1
        if (a >0) rowname=paste0("true mean = M + ",a)
            else if (a < 0) rowname=paste0("true mean = M - ",abs(a))
      # compute error costs for single level tests
        for ( j in 1:length(alpha))
            {alp<-alpha[j]
              TC1[j]=cRatio*Integrate(Type1a,0,max(0.03,D-8*sig2),D+8*sig2)
              TC2[j]=Integrate(Type2a,0,max(0.03,D-8*sig2),D+8*sig2)
              }
      #compute multi-level cost
        multi=round(sum((TC1+TC2)*CDiff),2)
       #now find approximate optimal by searching alpha space for minimum of single level tests
        optdel=400
        opt=rep(0,optdel)
        y=1/optdel
        for (j in 1:optdel) {
          alp<-y*j
          opt[j]=cRatio*Integrate(Type1a,0,max(0.03,D-8*sig2),D+8*sig2) +Integrate(Type2a,0,max(0.03,D-8*sig2),D+8*sig2)
        }
      # enter into table
        Table = rbind(Table,c(rowname,round(TC1+TC2,2),multi,paste0(round(min(opt),2)," (",round(y*which.min(opt),2),")")))
        # record lowest cost of single level tests
        minentry[row] = which.min(TC1+TC2)
        row=row+1
      } # end of true mean loop
    }# end of cost ratio loop
    
    formatTable(Table, "Table3", aa, alpha, minentry,length(alpha))
  } # end of function Table_3
  

  #############################################################
#' Tables_S3 is a function to compute expected average total error costs of one sided t-tests given two equal groups,
#' at different alpha levels and at simultaneously multiple levels.
#' Averages are over multiple runs in which costs are randomly generated. However, the expected ratio of Type I to Type II
#' costs is a parameter.
#' #' Prevalence of true findings is modelled using both zero-mean normal distributions and dichotomous distributions of effect sizes in the research environment
#' 
#' Output contains Table S3.1 and Table S3.2 of the paper. (Default alphas produce (b) tables; change for (a) tables.)
#' Details of the computations are in the paper


#' @param alpha # alphas in decreasing order
#' @param noit - number of iterations 
#' @param SS vector of group sample sizes, assuming equal groups
#' @param PP vector of prevalence in dichotomous scenarios
#' @param del vector of difference between "true" and "false" effect sizes in dichotomous scenario
#' @param MM vector of minimum meaningful effect in continuous scenarios
#' @param SD vector of SD of effect distribution in continuous scenarios
#' @param onetwo - on average, the ratio of costs of Type I to Type II errors 
#'
#' n= number of alphas
#' sample = current group sample size
#' M = current minimum meaningful effect
#' alp = current alpha
#' costs = vector of expected costs for sinle level tests, multilevel and optimal
#' P = prevalence in current computation in dichotomous scenario
#' EM =mean of true effect distribution in cts scenario (always set to 0 in current code)
#' minentry = vector of positions of lowest cost alphas for each scenario setting
#' row = row in Table
#' 
Tables_S3 <- function(alpha=c(.025,.0025,.0005), noit=2000, SS=c(12,96), PP=c(.5,.1), del=c(.15,.5),MM=c(0,.64), SD=c(.5,1), onetwo=4) 
  {
  # set up parameters for tests and scenarios
  n <<- length(alpha) # number of alphas
  # initialise
  C =rep(100,n)
  Table3.1 = matrix(nrow = 0, ncol = n + 3)
  Table3.2 = matrix(nrow = 0, ncol = n + 3)
  x=matrix(nrow=4,ncol=n)
  minentry = 100 * 0 # will contain entry with minimum cost per scenario
  row = 0
  costs = rep(0, n + 2) # will contain expected costs
  EM=0 # mean of true effect distribution in cts scenario
  
  
  #Functions to perform integration over effect sizes in a continuous research scenario
  # EM and v are mean and sd of distribution of true effects
  #cutoff is qt(1-alp,df) 

  IntegrateSim<-function(){
  
    I1=integrate(Type1,-10,M)$value #only get false positives if below cut off
    I2=integrate(Type2,M,10)$value # only get false negatives afterwards
    return(c(I1,I2))}
  
  Type1 <- function(ee) {
    dnorm(ee,EM,v)*(1-pt(cutoff+(M-ee)*sqrt(sample/2),df))
  }
  Type2 <- function(ee) {  
    dnorm(ee,EM,v)*pt(cutoff-(ee-M)*sqrt(sample/2),df)
  }
  
  # function to generate costs
  #' @param n number of single level tests = length of alpha
  #' @param onetwo - on average, the ratio of costs of Type I to Type II errors 
  #' Outputs matrix x with rows containing vector of differences in Type I costs between levels,
  #' differences in Type II costs, abs. Type I costs and abs. Type II costs for multi-level tests
  
  randomCosts<- function(n, onetwo){
    cD0=C*runif(n) # difference in costs between alphas
    cD1=C*runif(n)/onetwo
    c0=cD0 # total costs at each alpha
    c1=cD1
    for (i in 2:n){c0[i]=c0[i-1]+cD0[i]
                  c1[i]=c1[i-1]+cD1[i]}
    x[1,]=cD0;x[2,]=cD1;x[3,]=c0;x[4,]=c1
    return(x)
   }
  ########
  # start computations for dichotomous scenario
#
  row = 0 #start rows for dichotomous table
  for (d in del) {
    #vector of true effects
    # insert scenario risk difference  into table
    Table3.1 = rbind(Table3.1, matrix(
      data = " ",
      nrow = 1,
      ncol = n + 3
    ))
    Table3.1[length(Table3.1[, 1]), 1] = paste0("d =", d)
    row = row + 1
    #prevalence loop
    for (P in PP) {
      #loop for group sample size
      for (sample in SS){ 
        df<-2*sample-2
        dd = d*sqrt(sample/2)
        cutoff = qt(1-alpha, df)
        beta=pt(-dd+cutoff,df)
      
      #initialise for iterations
        sum=sumsq=rep(0,n+2)
        for (randomcount in 1:noit)
            { #set random costs
            rC=randomCosts(n,onetwo)
            
            #compute expected costs
            costs[n+1]=sum(rC[1,]*(1-P)*alpha +P*rC[2,]*beta) #multilevel tests
            costs[1:n]=(1-P)*rC[3,n]*alpha+P*rC[4,n]*beta # single level tests
            
             #now find optimal by searching alpha space for minimum of single level tests
            opt = rep(0, 3000)
            a = .0001
            jj = 1:10000
            opt=(1-P)*rC[3,n]*jj*a + P*rC[4,n]*pt(-dd+qt(1-jj*a,df),df)
            costs[n+2]=min(opt)
            
            # add to running sums
            #costs[1:(n+1)]=(costs[1:(n+1)]-costs[n+2])
            sum=sum+costs
            sumsq=sumsq+costs*costs
          } # end of iterations

        sum=sum/noit
        sumsq=round(sqrt(sumsq/noit-sum*sum),1)
      # add entry into tables
        rowname = paste0("P=", round(P, 2),", n=", sample*2) # paste prevalence and sample size into table
        Table3.1 = rbind(Table3.1, c(rowname, paste0(
            formattable(sum, format = "f", digits = 1)," (",sumsq,")")))
      
        minentry[row] = which.min(sum[1:n]) # store where minimum cost as this will be bolded in final table
         row = row + 1
        } # end sample size loop
      }# end prevalence loop
    } #end true effect loop
  
  #output results
  formatTable(Table3.1, "Table3.1", c(del,SS), alpha, minentry,n)
  
  
  #  now do cts scenario
  row = 0
  minentry = 100 * 0
  #loop varying mean of true effect distribution
  for (M in MM) {
    row = row + 1
    # insert scenario effect size distribution mean into table
    Table3.2 = rbind(Table3.2, matrix(
      data = " ",
      nrow = 1,
      ncol = n + 3
    ))
    Table3.2[length(Table3.2[, 1]), 1] = paste0("M = ", M)
    ## loop over variance of true effect distribution
    for (v in SD) {
      # loop over group sample size
      for (sample in SS){
        df<-2*sample-2
        #initialise for iterations
        sum=sumsq=rep(0,n+2)
        for (randomcount in 1:noit)
        { #set random costs
          rC=randomCosts(n,onetwo)
          #compute expected costs
          # compute  integrals of costs over effect size distribution in research scenario
          costs[n+1]=0
          #costs for tests
          for (j in 1:n) {
            alp=alpha[j]
            cutoff=qt(1-alp,df)
            I=IntegrateSim()
            costs[j]=rC[3,n]*I[1]+rC[4,n]*I[2] #single level tests
            costs[n+1]=costs[n+1] + rC[1,j]*I[1]+rC[2,j]*I[2]# multilevel test error costs
            }
           #now find optimal by searching alpha space for minimum of single level tests
          optdel=200
          opt=rep(0,optdel)
          a=1/optdel
          for (j in 1:optdel) {
            alp=a*j
            cutoff=qt(1-alp,df)
            I=IntegrateSim()
            opt[j]=rC[3,n]*I[1]+rC[4,n]*I[2]
          }
          costs[n+2]=opt[which.min(opt)]
          #costs[1:(n+1)]=(costs[1:(n+1)]-costs[n+2])
          sum=sum+costs
          sumsq=sumsq+costs*costs
        } # end of loop for iterations with random costs
        
        sum=sum/noit
        sumsq=round(sqrt(sumsq/noit-sum*sum),1)

      # write results into table
      rowname = paste0("s = ", round(v, 2),", n=", sample*2)# add sd of scenario effect distribution into table
      Table3.2 = rbind(Table3.2, c(rowname, paste0(
        formattable(sum, format = "f", digits = 1)," (",sumsq,")")))
      
      minentry[row] = which.min(costs[1:n]) # store minimum costs from single level tests
      row = row + 1
      } # end of loop for sample size
    }  # end of variance loop
  } # end of loop over  mean
  
  #output Table
  formatTable(Table3.2, "Table3.2", c(MM,SS), alpha, minentry,length(alpha))
}

####end of function Tables_S3


###################################################################
 
  # Code to produce Fig 1
#' @param alp = test level
#' @param a, sig1 = mean and sd of true effect size distribution
#' @param D, sig2 = mean and sd of anticipated effect size distribution
#' @param ytop = y axis maximum
  Fig1a<- function(alp=.025,  a=0, sig1=.2, D=.4,sig2=.1, ytop=4){
  ex=.001*(-1000:1000) # vector of effect sizes
   # First do Fig (a)
  #plot true ES distribution
  xlab=expression(paste("standardized effect size - ", italic("M")))
  plot(ex,dnorm(ex,a,sig1), type="l",xlab=xlab,xlim=c(-1,1),ylim=c(0,ytop),ylab=" ",bty="n")
  # plot distribution of anticipated ES
  lines(ex,dnorm(ex,D,sig2),lty="dashed")
  #plot an example sample distribution assuming anticipated effect size D
  nq =2*((qnorm(1-alp,0,1)+qnorm(.8,0,1))/D)^2 # group size if anticipated effect at mean
  lines(ex,sqrt(nq/2)*dt(ex*sqrt(nq/2),2*nq-2),lty="dotted")
  #put in line at critical value if testing at alp
  tt=qt(1-alp,2*nq-2)/sqrt(nq/2)
  lines(c(tt,tt),c(0,ytop))
  #put in axis
  lines(c(-1,1),c(0,0),lwd=.5)
  }
  #  Fig 1(b)
Fig1b<- function(alp=.025,  a=0, sig1=.2, D=.4,sig2=.1, ytop=4){
  D=.4
  d=(1:1000)/1000 # anticipated effect size 
  
  # compute corresponding total sample sizes
  nn=4*((qnorm(1-alp,0,1)+qnorm(.8,0,1))/d[1])^2 
  for (j in 2:length(d))nn=c(nn, 4*((qnorm(1-alp,0,1)+qnorm(.8,0,1))/d[j])^2)
  nn=as.integer(nn+1) 
  # plot sample size as function of anticipated effect size
  xlab=expression(paste("standardized effect size - ", italic("M")))
  plot(d,nn, type="l",xlab=xlab,xlim=c(0,1),ylim=c(0,4000),ylab="total sample size",bty="n")
  lines(d*1000,dnorm(d,D,sig2))
}

# Fig 1c
Fig1c<- function(alp=.025,  a=0, sig1=.2, D=.4,sig2=.1, ytop=4){
  xlab=expression(paste("anticipated standardized effect size - ",italic(M)))
  plot(nn,dnorm(d,D,sig2), type="l",xlab="total sample size",xlim=c(0,2000),ylim=c(0,4),ylab=" ",bty="n")
  lines(c(0,100000),c(0,0),lwd=.5)
 }
  ####################################

#functions to output supplementary material Figure S4.1a and S4.1b 
#' @param R1 = risk in untreated group
#' @param sample = group sample size
#' @param xtop, ytop respective maximum of axis ranges
#' @param EM, v = vectors of mean and of sd of true distributions
FigS4a <- function(R1 = 0.092,sample =1000, ytop=40,xtop=.1)
{
  M=-707/40000 #- cost ratio is boundary of meaningful effects
  ee = -.02 # effect size at which will draw sampling distribution
  sd=sqrt(((1-R1)*R1 +(1-ee-R1)*(ee+R1)) / sample)
  px = (-100:100) * .01*xtop # difference in cases
  dn = dnorm(px, ee, sd) #sampling distribution at -.02
  plot(px,dn,type = "l",xlab = paste0("risk difference"),xlim = c(-xtop, xtop),
    ylim = c(0, ytop),ylab = " ",main = "",bty = "n",lwd = 2,)
  #critical value at alpha=.62
  cutoff = qnorm(.23, M, sd) 
  lines(c(cutoff, cutoff), c(0, ytop))
  #critical value at 0.05
  cutoff = qnorm(.05, M, sd) 
  lines(c(cutoff, cutoff), c(0, ytop), lty = "dashed",lwd=2)
  #axis
  lines(c(-xtop, xtop), c(0, 0), lwd = 1)
}
################
FigS4b <- function(R1 = 0.092,EM=c(0,-.02),vv=c(0.015,0.025),ytop = 40)
{
  onetwo = 707/40000
  px = (-100:100) * .001 # difference in cases
  py=px
  plot(px, dnorm(px, EM[1], vv[1]),type = "l",xlim = c(-.1, .1),
    ylab = "",ylim = c(0, ytop),xlab = "risk difference",main = " ",
    bty = "n",lwd = 2,lty = "dashed") 
  #research scenario 2
  for (i in 1:length(px))py[i]=max(-R1,dnorm(px,EM[2], vv[2])[i])
  lines(px, py, lwd = 2) 
  #axes
  lines(c(-.1, .1), c(0, 0), lwd = 1)
  lines(c(-onetwo, -onetwo), c(0, ytop))
}
#################################################################
#functions to output Appendix Figures 3.1a and 3.1b

#' @param sample = group size (default has to be changed to get Fig b)
#' @param ytop = upper limit of y axis displayed
#' @param M = meaningful effect boundary
#' 
FigS3 <- function(sample= 6,ytop=1,M=.64)
{
  px=.01*(-200:200) #standardised effect sizes
  yy=sqrt(sample/2)*dt(px*sqrt(sample/2),2*sample-2) #sampling distribution
  par(xpd=F)
  # plot sampling distribution around M
  plot(px+M,yy,type="l",xlab="standardized effect size",xlim=c(-2,2+M),ylim=c(0,ytop),ylab=" ",bty="n")
  
  # plot distribution of true effects in research scenario
  lines(px,dnorm(px,0,.5),lwd=2)
  
  # plot error rates at alphas 0.025 and 0.0025
  ll=length(px+1)/2 # at zero
  cutoff1=qt(1-.025,2*sample-2) # critical value for test at 0.25
  cutoff2=qt(1-.0025,2*sample-2) # critical value for test at 0.25
  #type 1 errors occur when true effect size less than M
  lines(M+px[1:ll], 1-pt(-px[1:ll]*sqrt(sample/2)+cutoff1,2*sample-2),lty="dotted")
   #type 2 errors occur for effect size greater than M
  lines(M-px[1:ll], pt(px[1:ll]*sqrt(sample/2)+cutoff1,2*sample-2),lty="dashed")
  lines(M-px[1:ll], pt(px[1:ll]*sqrt(sample/2)+cutoff2,2*sample-2),lty="dashed",lwd=.5)
  
  # plot line at smallest meaningful
  lines(c(M,M),c(0,ytop),lwd=.5)
  #plot axes
  lines(c(-3,3),c(0,0),lwd=.5)
  lines(c(-2,2),c(1,1),lwd=.5) 
  #put legend outside boundary
  par(mar=c(5,5,5,5),xpd=T)
  legend(x="topright",inset=c(-.3,0),legend=c("true effect dist.","sampling dist.","Type I .025", "Type II .025", "Type II .0025"),
         lty=c(1,1,3,2,2),lwd=c(1.5,1,1,1,.5),cex=.65,bty="n")
  }


###########################################
#sample size calculation in section 3
R1=.092
R2=.044
onetwo=707/40000
factor=R1-R2-onetwo
np1=(1-R1)*R1
np2=(1-R2)*R2
np22=(1-R1+onetwo)*(R1-onetwo)
nsd=sqrt(np1+np2)
nsd2=sqrt(np1+np22)
za=qnorm(.05,0,1);zb=qnorm(.2,0,1)
((za*nsd2+zb*nsd)/factor)^2

pnorm((factor*sqrt(1000)+za*nsd2)/nsd)

    