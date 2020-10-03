mapping <- function(x0=0.4, delta=5, gama=1, beta=0.5){
  x1 <- exp(-delta*(x0*(1-x0**gama))**2)+beta
  return(x1)
}

returnMap <- function(delta=5, gama=1, beta=0.5, 
                      x0min=-0.5, x0max=1.5, ndiv=1000){
  x0 <- seq(x0min, x0max, length=ndiv)
  x1 <- mapping(x0, delta, gama, beta)
  
  #return values
  return(list(x0=x0, x1=x1))
}

cobweb <- function(x0=0.2, delta=5, gama=1, beta=0.5, ntimes=10){
  x <- vector()
  y <- vector()
  npoints <- 0
  for (i in 1:ntimes) {
    x1 <- mapping(x0, delta, gama, beta)
    
    npoints <- npoints+1 
    x[npoints] <- x0
    y[npoints] <- x0
    
    npoints <- npoints+1 
    x[npoints] <- x0
    y[npoints] <- x1
    
    x0 <- x1 
  }
  
  #return values
  return(list(x=x, y=y))
}

bifurcDiagram <- function(x0=0.2, delta=5, gama=seq(1,2,length=1000), beta=0.5, 
                          transient=100, maxiter=100){
  len_delta <- length(delta)
  len_gama  <- length(gama)
  len_beta  <- length(beta)
  
  x0save <- x0
  
  if(len_delta>1 && len_gama==1 && len_beta==1){
    select <- "delta"
    gama   <- seq(gama, gama, length = len_delta)
    beta   <- seq(beta, beta, length = len_delta)
  }else if(len_delta==1 && len_gama>1 && len_beta==1){
    select <- "gama"
    delta  <- seq(delta, delta, length = len_gama)
    beta   <- seq(beta, beta, length = len_gama)
  }else if(len_delta==1 && len_gama==1 && len_beta>1){
    select <- "beta"
    delta  <- seq(delta, delta, length = len_beta)
    gama   <- seq(gama, gama, length = len_beta)
  }else if(len_delta>1 && len_gama>1 && len_beta==1){
    if(len_delta!=len_gama){
      stop("Length of delta and gama must be equal.")
    }
    print("Select if you want to plot 'delta' or 'gama' as the horizontal axis")
    select <- readLines(n=1)
    beta   <- seq(beta, beta, length = len_delta)
  }else if(len_delta>1 && len_gama==1 && len_beta>1){
    if(len_delta!=len_beta){
      stop("Length of delta and beta must be equal.")
    }
    print("Select if you want to plot 'delta' or 'beta' as the horizontal axis")
    select <- readLines(n=1)
    gama   <- seq(gama, gama, length = len_delta)
  }else if(len_delta==1 && len_gama>1 && len_beta>1){
    if(len_delta!=len_beta){
      stop("Length of gama and beta must be equal.")
    }
    print("Select if you want to plot 'gama' or 'beta' as the horizontal axis")
    select <- readLines(n=1)
    delta  <- seq(delta, delta, length = len_gama)
  }else if(len_delta>1 && len_gama>1 && len_beta>1){
    if(len_delta!=len_beta || len_gama!=len_beta){
      stop("Length of gama, delta and beta must be equal.")
    }
  }else{
    return(print("One of the control parameters need to be a vector with length 
                 greater than one."))
  }
  print(c("Which control parameter is changing? ",select))
  
  px <- vector()
  py <- vector()
  npoints <- 0
  for(ix in 1:length(delta)){
    x0 <- x0save
    # removing transient time
    for (iter in 1:transient) {
      x1 <- mapping(x0, delta[ix], gama[ix], beta[ix])
      x0 <- x1
    }
    # save the next points in vectors px and py
    for (iter in 1:maxiter) {
      x1 <- mapping(x0, delta[ix], gama[ix], beta[ix])
      x0 <- x1
      
      npoints  <- npoints+1
      if(select == "delta"){
        px[npoints] <- delta[ix]
      }else if(select == "gama"){
        px[npoints] <- gama[ix]
      }else if(select == "beta"){
        px[npoints] <- beta[ix]
      }
      py[npoints] <- x1
    } 
  }
  #return values
  return(list(parameter=px, x=py))
}

parameterSpace <- function(delta=seq(1e-10,38,length=100), 
                           gama=1, 
                           beta=seq(-1,1,length=100),
                           x0=0.2, transient=50, maxiter=50){
  # Count the number of unique values. One of them must be equal to 1.
  ndelta <- length(unique(delta))
  ngama  <- length(unique(gama))
  nbeta  <- length(unique(beta))
  
  if(ndelta==1 && ngama>1 && ngama==nbeta){
    select <- "beta_vs_gama"
  }else if(ngama==1 && ndelta>1 && ndelta==nbeta){
    select <- "beta_vs_delta"
  }else if(nbeta==1 && ndelta>1 && ndelta==ngama){0
    select <- "gama_vs_delta"
  }else{
    return("ERROR: One control parameter needs to be a constant value, while the other two must change.")
  }
  print(paste("We are going to plot",select))
  
  pars    <- expand.grid(delta=delta, gama=gama, beta=beta)
  npmax   <- length(pars$delta)
  diagram <- vector(length=npmax)
  
  for (np in 1:npmax){
    y0     <- x0[1]
    delta2 <- pars$delta[np]
    gama2  <- pars$gama[np]
    beta2  <- pars$beta[np]
    
    # skipping transient time
    for(iterada in 1:transient[1]){
      y0 <- mapping(y0, delta2, gama2, beta2)
    }
    # calculate Lyapunov exponent for the next iterations
    sum <- 0
    for(iterada in 1:maxiter[1]){
      y0 <- mapping(y0, delta2, gama2, beta2)
      dx <- -2*delta2*y0*(1-y0**gama2)*(1-(gama2+1)*y0**gama2)*exp(-delta2*(y0*(1-y0**gama2))**2)
      sum <- sum+log(abs(dx))
    }
    lyap <- sum/maxiter[1]
    if(abs(lyap)>100 || is.nan(lyap)){
      lyap <- 100
    }
    diagram[np] <- lyap
  }
  if(select == "beta_vs_gama"){
    df  <- data.frame(gama=pars$gama, beta=pars$beta, lyap=diagram)
  }else if(select == "beta_vs_delta"){
    df <- data.frame(delta=pars$delta, beta=pars$beta, lyap=diagram)
  }else if(select == "gama_vs_delta"){
    df <- data.frame(delta=pars$delta, gama=pars$gama, lyap=diagram)
  }
  return(df)
}

plotParSpace <- function(A,zmin=-7, zmax=+7){
  library(ggplot2)
  
  negative <- colorRampPalette(c("red","yellow"))(1000)                      
  positive <- colorRampPalette(c("green","blue"))(1000)
  
  name1 <- colnames(A)[1]
  name2 <- colnames(A)[2]
  
  if(name1=="gama" && name2=="beta"){
    plt <- ggplot(A,mapping=aes(gama,beta))
  }else if(name1=="delta" && name2=="beta"){
    plt <- ggplot(A,mapping=aes(delta,beta))
  }else if(name1=="delta" && name2=="gama"){
    plt <- ggplot(A,mapping=aes(delta,gama))
  }
  plt <- plt+
    geom_raster(aes(fill=lyap))+
    scale_fill_gradientn(guide="colorbar",limits=c(zmin,zmax),colours=c(negative,"white",positive),na.value='red')+
    scale_x_continuous(expand = c(0, 0))+
    scale_y_continuous(expand = c(0, 0))
  return(plt)
}

highlightExtreme <- function(delta=seq(1e-10,38,length=200), 
                             gama=1, 
                             beta=seq(-1,1,length=200),
                             from="xr", to="xm", maxiter=2){
  ndelta <- length(unique(delta))
  ngama  <- length(unique(gama))
  nbeta  <- length(unique(beta))
  
  if(ndelta==1 && ngama>1 && ngama==nbeta){
    select <- "beta_vs_gama"
  }else if(ngama==1 && ndelta>1 && ndelta==nbeta){
    select <- "beta_vs_delta"
  }else if(nbeta==1 && ndelta>1 && ndelta==ngama){0
    select <- "gama_vs_delta"
  }else{
    return("ERROR: One control parameter needs to be a constant value, while the other two must change.")
  }
  print(paste("We are going to plot",select))
  
  
  pars    <- expand.grid(delta=delta, gama=gama, beta=beta)
  npmax   <- length(pars$delta)
  diagram <- vector(length=npmax)
  
  for (np in 1:npmax) {
    if(from=='xl'){
      x0save <- 0
    }else if(from=="xm"){
      x0save <- (1+pars$gama[np])**(-1/pars$gama[np])
    }else if(from=="xr"){
      x0save <- 1
    }
    if(to=='xl'){
      xfinal <- 0
    }else if(to=="xm"){
      xfinal <- (1+pars$gama[np])**(-1/pars$gama[np])
    }else if(to=="xr"){
      xfinal <- 1
    }
    x0 <- x0save
    for(iterada in 1:maxiter){
      x1 <- mapping(x0,pars$delta[np],pars$gama[np],pars$beta[np])
      x0 <- x1
    }
    diagram[np] <- x1-xfinal
  }
  # Plot results using ggplot2 and plotly
  library(ggplot2)
  if(select == "beta_vs_gama"){
    df  <- data.frame(gama=pars$gama, beta=pars$beta, lyapunov=diagram)
    plt <- ggplot(df,mapping=aes(gama,beta))
  }else if(select == "beta_vs_delta"){
    df <- data.frame(delta=pars$delta, beta=pars$beta, lyapunov=diagram)
    plt <- ggplot(df,mapping=aes(delta,beta))
  }else if(select == "gama_vs_delta"){
    df <- data.frame(delta=pars$delta, gama=pars$gama, lyapunov=diagram)
    plt <- ggplot(df,mapping=aes(delta,gama))
  }
  negative <- colorRampPalette(c("red","yellow"))(1000)                      
  positive <- colorRampPalette(c("green","blue"))(1000)
  
  plt <- plt +  
    geom_raster(aes(fill=lyapunov))+
    scale_fill_gradientn(guide="colorbar",limits=c(-5,5),colours=c(negative,"white",positive),na.value='red')+
    scale_x_continuous(expand = c(0, 0))+
    scale_y_continuous(expand = c(0, 0))
  return(plt)
}

extreme <- function(delta=c(0,38), gama=c(1,1), beta=c(-1,1), 
                    p1=c(0,1,-1), p2=c(0,1,0), theta=0,
                    from="xr", to="xm", maxiter=2, step=0.01){
  f <- function(delta, gama, beta, maxiter, from, to){
    if(from=='xl'){
      x0save <- 0
    }else if(from=="xm"){
      x0save <- (1+gama)**(-1/gama)
    }else if(from=="xr"){
      x0save <- 1
    }
    if(to=='xl'){
      xfinal <- 0
    }else if(to=="xm"){
      xfinal <- (1+gama)**(-1/gama)
    }else if(to=="xr"){
      xfinal <- 1
    }
    x0 <- x0save
    for(iterada in 1:maxiter){
      x1 <- mapping(x0,delta,gama,beta)
      x0 <- x1
    }
    return(x1-xfinal)
  }
  # Finding initial point
  fa <- f(p1[1], p1[2], p1[3], maxiter, from, to) #signal of f for p1
  fb <- f(p2[1], p2[2], p2[3], maxiter, from, to) #signal of f for p2
  if(fa*fb>0){
    return(print("Please, check the ranges of p1 and p2"))
  }
  deltaa <- p1[1]
  deltab <- p2[1]
  gamaa  <- p1[2]
  gamab  <- p2[2]
  betaa  <- p1[3]
  betab  <- p2[3]
  for(ntimes in 1:50){
    # Bisection method to find the solution for the initial point
    deltam <- (deltaa+deltab)/2
    gamam  <- (gamaa+gamab)/2
    betam  <- (betaa+betab)/2
    
    fm <- f(deltam, gamam, betam, maxiter, from, to)
    if(abs(fm) < 1e-10){
      delta0 <- deltam
      gama0  <- gamam
      beta0  <- betam
      break
    }
    if(fm*fa > 0){
      deltaa <- deltam
      gamaa  <- gamam
      betaa  <- betam
    }else{
      deltab <- deltam
      gamab  <- gamam
      betab  <- betam
    }
  }
  if(ntimes>=50){
    return(print("ERROR: Check p1 and p2"))
  }
  # Finding second point
  # Another bisection method is considered but now changing the angle theta
  if(is.numeric(theta)){
    tetaa <- theta-pi/2
  }else{
    print("ERROR, unknown direction for theta.")
  }
  tetab <- theta+pi/2
  
  delta1  <- vector()
  gama1   <- vector()
  beta1   <- vector()
  npoints <- 0
  for(nsteps in 1:1000000){
    if(delta[1] == delta[2]){
      deltam <- delta0
      gamam  <- gama0+(gama[2]-gama[1])*step*cos(tetaa)
      betam  <- beta0+(beta[2]-beta[1])*step*sin(tetaa)
    }else if(gama[1] == gama[2]){
      deltam <- delta0+(delta[2]-delta[1])*step*cos(tetaa)
      gamam  <- gama0
      betam  <- beta0+(beta[2]-beta[1])*step*sin(tetaa)
    }else if(beta[1] == beta[2]){
      deltam <- delta0+(delta[2]-delta[1])*step*cos(tetaa)
      gamam  <- gama0+(gama[2]-gama[1])*step*sin(tetaa)
      betam  <- beta0
    }
    fa    <- f(deltam, gamam, betam, maxiter, from, to)
    
    if(delta[1] == delta[2]){
      deltam <- delta0
      gamam  <- gama0+(gama[2]-gama[1])*step*cos(tetab)
      betam  <- beta0+(beta[2]-beta[1])*step*sin(tetab)
    }else if(gama[1] == gama[2]){
      deltam <- delta0+(delta[2]-delta[1])*step*cos(tetab)
      gamam  <- gama0
      betam  <- beta0+(beta[2]-beta[1])*step*sin(tetab)
    }else if(beta[1] == beta[2]){
      deltam <- delta0+(delta[2]-delta[1])*step*cos(tetab)
      gamam  <- gama0+(gama[2]-gama[1])*step*sin(tetab)
      betam  <- beta0
    }
    fb    <- f(deltam, gamam, betam, maxiter, from, to)
    if(fa*fb>0){
      return(print(c("!Something is wrong with tetaa and tetab",nsteps)))
    }
    for(itimes in 1:50){
      tetam  <- (tetaa+tetab)/2
      
      if(delta[1] == delta[2]){
        deltam <- delta0
        gamam  <- gama0+(gama[2]-gama[1])*step*cos(tetam)
        betam  <- beta0+(beta[2]-beta[1])*step*sin(tetam)
      }else if(gama[1] == gama[2]){
        deltam <- delta0+(delta[2]-delta[1])*step*cos(tetam)
        gamam  <- gama0
        betam  <- beta0+(beta[2]-beta[1])*step*sin(tetam)
      }else if(beta[1] == beta[2]){
        deltam <- delta0+(delta[2]-delta[1])*step*cos(tetam)
        gamam  <- gama0+(gama[2]-gama[1])*step*sin(tetam)
        betam  <- beta0
      }
      fm <- f(deltam, gamam, betam, maxiter, from, to)
      
      if(abs(fm)<1e-10){
        npoints <- npoints+1
        delta1[npoints] <- deltam
        gama1[npoints]  <- gamam
        beta1[npoints]  <- betam
        
        if(deltam<delta[1] || deltam>delta[2] || gamam<gama[1] || gamam>gama[2] || betam<beta[1] || betam>beta[2]){
          if(delta[1] == delta[2]){
            df <- data.frame(gama=gama1, beta=beta1)
          }else if(gama[1] == gama[2]){
            df <- data.frame(delta=delta1, beta=beta1)
          }else if(beta[1] == beta[2]){
            df <- data.frame(delta=delta1,gama=gama1)
          }
          return(df)
        }
        delta0 <- deltam
        gama0  <- gamam
        beta0  <- betam
        
        tetaa  <- tetam-pi/2
        tetab  <- tetam+pi/2
        break
      }
      if(fm*fa>0){
        tetaa <- tetam
      }else{
        tetab <- tetam
      }
    }
    if(itimes>=50){
      return(print("ERROR: we do not found a value for (beta1,delta1)"))
    }
  }
  stop(paste("Is something wrong? The simulations is not stopping!"), call. = FALSE, domain = NULL)
  geterrmessage()
}

