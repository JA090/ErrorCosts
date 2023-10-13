#' Contains functions to produce Tables and Figures in the paper "Error costs may justify simultaneous testing of one hypothesis at multiple alpha levels"
#' by J Aisbett


require(gt)
require(formattable)
require(gtExtras)
require(tibble)
chromote::default_chromote_object()

#'formatTable() is a function to format and write out tables as .docx
#'Table is matrix input
#'name is the name used for the output word document i.e. name.docx
#' repeater defines groups of entries in final table
#'minentry is a vector containing column no. with the minimum cost of the single level tests
formatTable <- function(Table, name, repeater, alpha, minentry) {
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
  xx1 = gt_add_divider(
    xx1,
    columns = c(n + 1, n + 2),
    color = "black",
    sides = "right",
    style = "dashed"
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
  ccc = c[1]
  for (j in 2:n)
    ccc = c(ccc, c[j] - c[j - 1])  # differences between benefits for different alphas
  ## set up internal computational stuff
  Table2 = matrix(nrow = 0, ncol = n + 3)
  Table1 = matrix(nrow = 0, ncol = n + 3)
  minentry = 100 * 0
  row = 0
  costs = rep(0, n + 2) # will contain expected costs
  II <<- rep(0, n) # will contain total error costs  in cts scenarios
  
  
  #Functions to perform integration over effect sizes in a continuous research scenario
  Integrate <- function(z) {
    for (i in 1:length(z)) {
      zalp <<- z[i]
      I2 <<-
        integrate(ttt2, .00001, R1 - onetwo, stop.on.error = F)$value #only get false negatives if below cut off
      II[i] <<-
        (I2 + integrate(ttt1, R1 - onetwo, 1, stop.on.error = F)$value) # only get false positives afterward
    if (i < 4) print(c(z[i],II[i]-I2,I2))}
    print(c(sd2,p1))
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
  
  
  ########################
  p1 = (1 - R1) * R1 #sample distribution variance of untreated risk
  p22=(1-R1+onetwo)*(R1-onetwo) # sample distribution of treated risk if this is R1-onetwo
  sd2=sqrt((p1+p22)/sample)
  # start computations for dichotomous scenario
  row = 0 #start rows for dichotomous table

  for (R2 in R22) {
    #vector of true effects
    # insert scenario risk difference  into table
    Table1 = rbind(Table1, matrix(
      data = " ",
      nrow = 1,
      ncol = n + 3
    ))
    Table1[length(Table1[, 1]), 1] = paste0("RD =", R2 - R1)
    row = row + 1
    for (P in PP) {
      #prevalence loop
      p2 = (1 - R2) * R2 #sample distribution variance if  treated risk = R2
      factor = R1 - R2 - onetwo
      sd = sqrt((p1+p2)/sample)
      zalp= qnorm(alpha, 0, 1)
      # sum type I and type II error costs; the latter
      costs[1:n] = (1-P) * onetwo * C * alpha + P * C * factor * (1 -
                                                                      pnorm((zalp*sd2+factor)/sd, 0,1)) # single level tests
      costs[n + 1] = sum(costs[1:n] * ccc) / C # multilevel
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
  formatTable(Table1, "Table1", R22, alpha, minentry)
  
  
  #  cts scenario
  row = 0
  minentry = 100 * 0
  for (EM in MM) {
    #mean of true effect distribution
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
      costs[1:n] = Integrate(qnorm(alpha, 0, 1)) # costs
      # costs for multilevel test
      costs[n + 1] = sum(costs[1:n] * ccc) / C
      #find optimal
      opt = rep(0, 2000)
      j = 1:2000
      a = .0005
      opt = Integrate(qnorm(a * j, 0, 1))
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
  #output Table
  formatTable(Table2, "Table2", vv, alpha, minentry)
}


####################################################
#' Table_3_4 is a function to compute expected average total error costs of t-tests given two equal groups,
#' at different alpha levels and at simultaneously multiple levels.
#' Averages are over multiple runs in which costs are randomly generated. However, the expected ratio of Type I to Type II
#' costs is a parameter.
#' Output contains Table 3 and Table 4 of the paper. Details of the computations are in that paper

#' Prevalence of true findings is modelled using both zero-mean normal distributions and dichotomous distributions of effect sizes in the research environment


#' @param alpha # alphas in decreasing order
#' @param noits - number of iterations 
#' @param sample group sample sizes, assuming equal groups
#' @param PP vector of prevalence in dichotomous scenarios
#' @param del vector of difference between "true" and "false" effect sizes in dichotomous scenario
#' @param MM vector of minimum meaningful effect in continuous scenarios
#' @param SD vector of SD of effect distribution in continuous scenarios
#' @param onetwo - on average, the ratio of costs of Type I to Type II errors 
#'
#' n= number of alphas
#' costs = vector of expected costs for sinle level tests, multilevel and optimal
#' P - prevalence in current computation in dichotomous scenario
#' minentry - vector of positions of lowest cost alphas for each scenario setting
#' row - row in Table
Table_3_4 <- function(alpha=c(.025,.0025), noit=2000, SS=c(12,96), PP=c(.5,.1), del=c(.15,.5),MM=c(0,.64), SD=c(.5,1), onetwo=4) {
  # set up parameters for tests and scenarios
  n <<- length(alpha) # number of alphas
  # initialise
  C =rep(100,n)
  Table3 = matrix(nrow = 0, ncol = n + 3)
  Table4 = matrix(nrow = 0, ncol = n + 3)
  x=matrix(nrow=4,ncol=n)
  minentry = 100 * 0 # will contain entry with minimum cost per scenario
  row = 0
  costs = rep(0, n + 2) # will contain expected costs
  EM=0 # mean of true effect distribution in cts scenario
  
  
  #Functions to perform integration over effect sizes in a continuous research scenario

  IntegrateSim<-function(){
  
    I1=integrate(Type1,-10,M)$value #only get false positives if below cut off
    I2=integrate(Type2,M,10)$value # only get false negatives afterwards
    return(c(I1,I2))}
  
  Type1 <- function(ee) {
    dnorm(ee,EM,v)*(1-pt(cutoff+(M-ee)*forSE,df))
  }
  Type2 <- function(ee) {  
    dnorm(ee,EM,v)*pt(cutoff-(ee-M)*forSE,df)
  }
  
  # function to generate costs
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
  ########################
  # start computations for dichotomous scenario
  row = 0 #start rows for dichotomous table
  for (d in del) {
    #vector of true effects
    # insert scenario risk difference  into table
    Table3 = rbind(Table3, matrix(
      data = " ",
      nrow = 1,
      ncol = n + 3
    ))
    Table3[length(Table3[, 1]), 1] = paste0("d =", d)
    row = row + 1
    for (P in PP) {
      #prevalence loop
      for (sample in SS){ #loop for group sample size
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
        sum=sum+costs
        sumsq=sumsq+costs*costs
        } # end of iterations

      sum=sum/noit
      sumsq=round(sqrt(sumsq/noit-sum*sum),1)
      # add entry into tables
      rowname = paste0("P=", round(P, 2),", n=", sample*2) # paste prevalence and sample size into table
      Table3 = rbind(Table3, c(rowname, paste0(
        formattable(sum, format = "f", digits = 1)," (",sumsq,")")))
      
      minentry[row] = which.min(sum[1:n]) # store where minimum cost as this will be bolded in final table
      row = row + 1
      } # end sample size loop
    }# end prevalence loop
  } #end true effect loop
  #output results
  formatTable(Table3, "Table3", c(del,SS), alpha, minentry)
  
  
  #  cts scenario
  row = 0
  minentry = 100 * 0
  for (M in MM) {
    #mean of true effect distribution
    row = row + 1
    # insert scenario effect size distribution mean into table
    Table4 = rbind(Table4, matrix(
      data = " ",
      nrow = 1,
      ncol = n + 3
    ))
    Table4[length(Table4[, 1]), 1] = paste0("M = ", M)
    ## loop over variance of true effect distribution
    for (v in SD) {
      for (sample in SS){
        df<-2*sample-2
        forSE = sqrt(sample/2)
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
          sum=sum+costs
          sumsq=sumsq+costs*costs
        } # end of iterations
        
        sum=sum/noit
        sumsq=round(sqrt(sumsq/noit-sum*sum),1)

      # write results into table
      rowname1 = paste0("s = ", round(v, 2),", n=", sample*2)# add sd of scenario effect distribution into table
      Table4 = rbind(Table4, c(rowname, paste0(
        formattable(sum, format = "f", digits = 1)," (",sumsq,")")))
      
      minentry[row] = which.min(costs[1:n]) # store minimum costs from single level tests
      row = row + 1
      } # end of loop for sample size
    }  # end of variance loop
  } # end of loop over  mean
  #output Table
  formatTable(Table4, "Table4", c(MM,SS), alpha, minentry)
}

####end of function Table_3_4


###################################################################

#functions to output Appendix Figure 1a and b 
#' @param R1 = risk in untreated group
Fig1a <- function(R1 = 0.092)
{
  type = "l"
  M=-707/40000
  ee = M
  sample = 1000
  sd=sqrt(((1-R1)*R1 +(1-ee-R1)*(ee+R1)) / sample)
  ytop = 40
  xtop=.1
  px = (-100:100) * .01*xtop # difference in cases
  dn = dnorm(px, ee, sd) 
  plot(
    px,
    dn,
    type = type,
    xlab = paste0("risk difference"),
    xlim = c(-xtop, xtop),
    ylim = c(0, ytop),
    ylab = " ",
    main = "",
    bty = "n",
    lwd = 2,
    lty = "dashed"
  )
  ee = -.03
  sd=sqrt(((1-R1)*R1 +(1-ee-R1)*(ee+R1)) / sample)
  cutoff = qnorm(.23, M, sd) #cutoff at alpha=.62
  lines(c(cutoff, cutoff), c(0, ytop))
  cutoff = qnorm(.05, M, sd) #cutoff at 0.05
  lines(c(cutoff, cutoff), c(0, ytop), lty = "dotted")
  lines(px, dnorm(px, ee, sd), lwd = 2) #sampling distribution at -0.03
  lines(c(-xtop, xtop), c(0, 0), lwd = 1)
}
########################################################
#' @param R1 = risk in untreated group
Fig1b <- function(R1 = 0.092,EM=c(0,-.02),vv=c(0.015,0.025))
{
  type = "l"
  onetwo = 707/40000
  ytop = 40
  px = (-100:100) * .001 # difference in cases
  py=px
  plot(
    px,
    dnorm(px, EM[1], vv[1]),
    type = type,
    xlim = c(-.1, .1),
    ylab = "",
    ylim = c(0, ytop),
    xlab = "risk difference",
    main = " ",
    bty = "n",
    lwd = 2,
    lty = "dashed"
  ) #research scenario 2
  for (i in 1:length(px))py[i]=max(-.092,dnorm(px,EM[2], vv[2])[i])
  
  lines(px, py, lwd = 2) #research scenario 2
  lines(c(-.1, .1), c(0, 0), lwd = 1)
  lines(c(-onetwo, -onetwo), c(0, ytop))
}
#################################################################
#functions to output Appendix Figures 2a and b
#' @param sample = group size
#' @param ytop = upper limit of y axis displayed
#' 
Fig2 <- function(sample= 6,ytop=1)
{
  type="l"
  M=.64
  px=.01*(-200:200)
  ll=length(px+1)/2
  par(xpd=F)
  # plot sampling distribution around M
  plot(px+M,sqrt(sample/2)*dt(px*sqrt(sample/2),2*sample-2),type=type,xlab="standardised effect size",xlim=c(-2,2+M),ylim=c(0,ytop),ylab=" ",bty="n")
  
  # plot distribution of true effects in research scenario
  lines(px,dnorm(px,0,.5),lwd=2)
  
  # error rates at alphas 0.025 and 0.0025
  # .8 is for step size
  #type 1
  lines(M+px[1:ll], 1-pt(-.8*px[1:ll]*sqrt(sample/2)+qt(1-.025,2*sample-2),2*sample-2),lty="dotted")
   #type 2
  lines(M-px[1:ll], pt(.8*px[1:ll]*sqrt(sample/2)+qt(1-.025,2*sample-2),2*sample-2),lty="dashed")
  lines(M-px[1:ll], pt(.8*px[1:ll]*sqrt(sample/2)+qt(1-.0025,2*sample-2),2*sample-2),lty="dashed",lwd=.5)
  
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

#sample size calculation
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
((za*nsd2+zb*nsd)/(factor))^2

 
    