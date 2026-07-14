library(MASS)
library(ggplot2)

set.seed(123)


###########################################################
# Moment function
###########################################################

moments <- function(theta,Y,X,Z){
  
  residual <- Y - theta*X - 0.5*X^2
  
  Z * residual
  
}


###########################################################
# Nonlinear GMM estimator
###########################################################

gmm_est <- function(Y,X,Z){
  
  objective <- function(theta){
    
    g <- colMeans(
      moments(theta,Y,X,Z)
    )
    
    sum(g^2)
    
  }
  
  optimize(
    objective,
    interval=c(0,4),
    tol=1e-7
  )$minimum
  
}



###########################################################
# Coverage as q increases
###########################################################

simulate_q <- function(q,
                       R=100,
                       B=200,
                       n=500,
                       theta0=2){
  
  
  coverage <- matrix(
    0,R,3
  )
  
  lengths <- matrix(
    0,R,3
  )
  
  boot_time <- numeric(R)
  fid_time <- numeric(R)
  
  
  for(r in 1:R){
    
    
    #######################################################
    # Generate data
    #######################################################
    
    X <- rnorm(n)
    
    Z <- matrix(
      rnorm(n*q),
      n,q
    )
    
    
    u <- rt(n,df=5)
    
    Y <- theta0*X +
      .5*X^2 +
      u
    
    
    #######################################################
    # Estimate theta
    #######################################################
    
    theta_hat <- gmm_est(
      Y,X,Z
    )
    
    
    #######################################################
    # Estimate moment covariance
    #######################################################
    
    eps <- 1e-5
    
    
    G <- (
      colMeans(
        moments(theta_hat+eps,Y,X,Z)
      )
      -
        colMeans(
          moments(theta_hat-eps,Y,X,Z)
        )
    )/(2*eps)
    
    
    Omega <- cov(
      moments(theta_hat,Y,X,Z)
    )
    
    
    # shrinkage to avoid singularity
    lambda <- 0.01
    
    Omega <- 
      (1-lambda)*Omega +
      lambda*diag(q)
    
    
    Omega_inv <- solve(Omega)
    
    
    V <- as.numeric(
      solve(
        t(G)%*%
          Omega_inv%*%
          G
      )
    )/n
    
    
    se <- sqrt(V)
    
    
    
    #######################################################
    # Wald
    #######################################################
    
    ci <- theta_hat+c(-1,1)*1.96*se
    
    
    coverage[r,1] <-
      theta0>=ci[1] &
      theta0<=ci[2]
    
    
    lengths[r,1] <- diff(ci)
    
    
    
    #######################################################
    # Bootstrap
    #######################################################
    
    start <- Sys.time()
    
    boot <- numeric(B)
    
    
    for(b in 1:B){
      
      id <- sample(
        1:n,
        n,
        replace=TRUE
      )
      
      boot[b] <- gmm_est(
        Y[id],
        X[id],
        Z[id,]
      )
      
    }
    
    
    end <- Sys.time()
    
    boot_time[r] <-
      as.numeric(
        end-start,
        units="secs"
      )
    
    
    ci <- quantile(
      boot,
      c(.025,.975)
    )
    
    
    coverage[r,2] <-
      theta0>=ci[1] &
      theta0<=ci[2]
    
    
    lengths[r,2] <- diff(ci)
    
    
    
    #######################################################
    # SCD
    #######################################################
    
    start <- Sys.time()
    
    fid <- numeric(B)
    
    
    ##################################################
    # Gaussian moment shock distribution
    ##################################################
    
    g_mat <- moments(
      theta_hat,
      Y,
      X,
      Z
    )
    
    Omega_fid <- cov(g_mat)
    
    # shrinkage to avoid singularity when q is large
    lambda <- 0.01
    
    Omega_fid <- 
      (1-lambda)*Omega_fid +
      lambda*diag(q)
    
    
    ##################################################
    # Generate SCD draws
    ##################################################
    
    for(b in 1:B){
      
      ##################################################
      # Gaussian draw of moment uncertainty
      ##################################################
      
      shock <- mvrnorm(
        1,
        mu=rep(0,q),
        Sigma=Omega_fid
      )
      
      target <- shock/sqrt(n)
      
      
      ##################################################
      # invert moment equation
      ##################################################
      
      objective <- function(theta){
        
        g <- colMeans(
          moments(
            theta,
            Y,
            X,
            Z
          )
        )
        
        sum(
          (g-target)^2
        )
        
      }
      
      fid[b] <- optimize(
        objective,
        interval=c(0,4),
        tol=1e-7
      )$minimum
      
    }
    
    
    end <- Sys.time()
    
    fid_time[r] <-
      as.numeric(
        end-start,
        units="secs"
      )
    
    
    ci <- quantile(
      fid,
      c(.025,.975)
    )
    
    
    coverage[r,3] <-
      theta0>=ci[1] &
      theta0<=ci[2]
    
    
    lengths[r,3] <- diff(ci)
    
    
    if(r %% 25 == 0){
      cat("Completed", r, "of", R, "\n")
    }
    
  }
  
  
  data.frame(
    
    q=q,
    
    Method=c(
      "Wald",
      "Bootstrap",
      "SCD"
    ),
    
    Coverage=c(
      mean(coverage[,1]),
      mean(coverage[,2]),
      mean(coverage[,3])
    ),
    
    Avg_Length=c(
      mean(lengths[,1]),
      mean(lengths[,2]),
      mean(lengths[,3])
    ),
    
    Avg_Time=c(
      NA,
      mean(boot_time),
      mean(fid_time)
    )
    
  )
  
}


###########################################################
# Run experiment
###########################################################

q_results <- do.call(
  rbind,
  lapply(
    c(10,50,100,250),
    simulate_q
  )
)


print(q_results)

##########

###########################################################
# Coverage versus n
# Holding q/n = 0.2
###########################################################

n_values <- c(100,250,500,1000)


n_results <- do.call(
  rbind,
  lapply(
    n_values,
    function(n){
      
      simulate_q(
        q = floor(0.2*n),
        n = n,
        R = 100,
        B = 200
      )
      
    }
  )
)


print(n_results)



###########################################################
# Prepare coverage plot
###########################################################

coverage_plot <- data.frame(
  
  n = rep(n_values, each=3),
  
  q = rep(floor(0.2*n_values), each=3),
  
  Method = rep(
    c("Wald",
      "Bootstrap",
      "SCD"),
    times=length(n_values)
  ),
  
  Coverage = n_results$Coverage
  
)



###########################################################
# Plot coverage versus n
###########################################################

ggplot(
  coverage_plot,
  aes(
    x=n,
    y=Coverage,
    group=Method,
    color=Method
  )
)+
  geom_line(size=1)+
  geom_point(size=2)+
  geom_hline(
    yintercept=.95,
    linetype="dashed"
  )+
  ylim(.5,1)+
  theme_minimal()+
  labs(
    y="Coverage Probability",
    x="Sample Size",
    subtitle="Moment dimension q/n held fixed at 0.2"
  )



#######################
###########################################################
# Runtime scaling
# Bootstrap vs SCD 
###########################################################

runtime_q <- function(q,
                      B=500,
                      n=500){
  
  
  X <- rnorm(n)
  
  Z <- matrix(
    rnorm(n*q),
    n,q
  )
  
  
  Y <- 2*X + .5*X^2 + rt(n,5)
  
  
  theta_hat <- gmm_est(
    Y,X,Z
  )
  
  
  ###########################################################
  # Estimate covariance matrix of moments
  ###########################################################
  
  g_mat <- moments(
    theta_hat,
    Y,
    X,
    Z
  )
  
  
  Omega <- cov(g_mat)
  
  
  # covariance regularization
  lambda <- 0.01
  
  Omega <- 
    (1-lambda)*Omega +
    lambda*diag(q)
  
  
  
  ###########################################################
  # Estimate Jacobian
  ###########################################################
  
  eps <- 1e-5
  
  G <- (
    colMeans(
      moments(theta_hat+eps,Y,X,Z)
    )
    -
      colMeans(
        moments(theta_hat-eps,Y,X,Z)
      )
  )/(2*eps)
  
  
  G <- matrix(
    G,
    ncol=1
  )
  
  
  
  ###########################################################
  # Bootstrap
  # Repeated nonlinear GMM estimation
  ###########################################################
  
  boot_time <- system.time({
    
    for(b in 1:B){
      
      id <- sample(
        1:n,
        n,
        replace=TRUE
      )
      
      gmm_est(
        Y[id],
        X[id],
        Z[id,]
      )
      
    }
    
  })[["elapsed"]]
  
  
  
  
  ###########################################################
  # SCD
  # Draws from N(0, Omega_hat)
  ###########################################################
  
  fid_time <- system.time({
    
    
    # Cholesky decomposition of covariance matrix
    L <- chol(Omega)
    
    
    # Linear fiducial transformation
    Omega_inv <- solve(Omega)
    
    A <- solve(
      t(G)%*%
        Omega_inv%*%
        G
    )%*%
      t(G)%*%
      Omega_inv
    
    
    fid_draws <- numeric(B)
    
    
    for(b in 1:B){
      
      # Gaussian moment shock
      Z_star <- t(L)%*%
        rnorm(q)
      
      
      fid_draws[b] <-
        theta_hat +
        as.numeric(
          A %*% Z_star
        )/sqrt(n)
      
    }
    
    
  })[["elapsed"]]
  
  
  
  data.frame(
    
    q=q,
    
    Bootstrap=boot_time,
    
    SCD=fid_time,
    
    Speedup=boot_time/fid_time
    
  )
}



###########################################################
# Run runtime experiment
###########################################################

runtime_results <- do.call(
  rbind,
  lapply(
    c(10,50,100,250,500),
    runtime_q
  )
)


print(runtime_results)



###########################################################
# Reshape for plotting
###########################################################

runtime_plot <- reshape(
  runtime_results,
  varying=c(
    "Bootstrap",
    "SCD"
  ),
  v.names="Time",
  timevar="Method",
  times=c(
    "Bootstrap",
    "SCD"
  ),
  direction="long"
)



###########################################################
# Runtime scaling plot
###########################################################

ggplot(
  runtime_plot,
  aes(
    x=q,
    y=Time,
    color=Method
  )
)+
  geom_line(size=1)+
  geom_point(size=2)+
  scale_y_log10()+
  theme_minimal()+
  labs(
    y="Runtime (seconds, log scale)",
    x="Number of Moments"
  )