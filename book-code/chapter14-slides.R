library(rethinking)


# N households
N <- 25
dyads <- t(combn(N,2))
N_dyads <- nrow(dyads)
# simulate "friendships" in which ties are shared
f <- rbern(N_dyads,0.1) # 10% of dyads are friends
# now simulate directed ties for all individuals
# there can be ties that are not reciprocal
alpha <- (-3) # base rate of ties; -3 ~= 0.05
y <- matrix(NA,N,N) # matrix of ties
p_tie <- NA
for ( i in 1:N ) for ( j in 1:N ) {
  if ( i != j ) {
    # directed tie from i to j
    ids <- sort( c(i,j) )
    the_dyad <- which( dyads[,1]==ids[1] & dyads[,2]==ids[2] )
    p_tie <- f[the_dyad] + (1-f[the_dyad])*inv_logit( alpha )
    y[i,j] <- rbern( 1 , p_tie )
  }
}#ij
# N households


# now simulate gifts
giftsAB <- rep(NA,N_dyads)
giftsBA <- rep(NA,N_dyads)
lambda <- log(c(0.5,2)) # rates of giving for y=0,y=1
for ( i in 1:N_dyads ) {
  A <- dyads[i,1]
  B <- dyads[i,2]
  giftsAB[i] <- rpois( 1 , exp( lambda[1+y[A,B]] ) )
  giftsBA[i] <- rpois( 1 , exp( lambda[1+y[B,A]] ) )
}
# draw net


sim_data <- list(
  GAB = giftsAB,
  GBA = giftsBA,
  N_dyads = N_dyads,
  D = 1:N_dyads
)
# dyad model
f_dyad <- alist(
  GAB ~ poisson( lambdaAB ),
  GBA ~ poisson( lambdaBA ),
  log(lambdaAB) <-  T[D,1] ,
  log(lambdaBA) <-  T[D,2] ,
  a ~ normal(0,1),
  ## dyad effects
  transpars> matrix[N_dyads,2]:T <-
    compose_noncentered( rep_vector(sigma_T,2) , L_Rho_T , Z ),
  matrix[2,N_dyads]:Z ~ normal( 0 , 1 ),
  cholesky_factor_corr[2]:L_Rho_T ~ lkj_corr_cholesky( 2 ),
  sigma_T ~ exponential(1),
  ## compute correlation matrix for dyads
  gq> matrix[2,2]:Rho_T <<- Chol_to_Corr( L_Rho_T )
)
mGD <- ulam( f_dyad , data=sim_data , chains=4 , cores=4 , iter=2000 )

precis(mGD, depth = 3, pars = c("a", "sigma_T", "T"))

post <- extract.samples(mGD)
with(post, exp(a + T))
