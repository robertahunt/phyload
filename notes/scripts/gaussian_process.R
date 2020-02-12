library(rstan)
options(mc.cores=4)

gp.model.code <- "
data{
	int N;
  real n_iid[N];
  real n_epi[N];
  vector[N] mu;
  vector[N] y;
}

parameters{
  real<lower=-1,upper=1> relative_worth;
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
}

transformed parameters{
  // vector[N] mu = rep_vector(0, N);
  real n_eff[N];
  for (i in 1:N) {
    n_eff[i] = n_iid[i] + relative_worth * n_epi[i];
  }
}

model{
  matrix[N, N] L_K;
  matrix[N, N] K = cov_exp_quad(n_eff, alpha, rho);
  real sq_sigma = square(sigma);

  // diagonal elements
  for (n in 1:N)
    K[n, n] = K[n, n] + sq_sigma;

  L_K = cholesky_decompose(K);

  relative_worth ~ uniform(-1,1);
  
  rho ~ inv_gamma(5, 5);
  alpha ~ std_normal();
  sigma ~ std_normal();

  y ~ multi_normal_cholesky(mu, L_K);
}
"

gp.model <- stan_model(model_code=gp.model.code)

file.name <- "~/Downloads/2020-01-06_csv/agg_rf.mean.csv"

sample.indices <- sample.int(675,100)

data <- read.csv(file.name,stringsAsFactors=FALSE)
data <- data[data$X.d. == 0,]
data <- data[sample.indices,]
data <- list(N=dim(data)[1],
             n_iid=data$X.n_i.,
             n_epi=data$X.n_e.,
             y=log(data$value),
             mu=rep(0,dim(data)[1]))

fit.0 <- sampling(gp.model,data=data,iter=1000,chains=4)


data <- read.csv(file.name,stringsAsFactors=FALSE)
data <- data[data$X.d. == 0.5,]
data <- data[sample.indices,]
data <- list(N=dim(data)[1],
             n_iid=data$X.n_i.,
             n_epi=data$X.n_e.,
             y=log(data$value),
             mu=rep(0,dim(data)[1]))

fit.0.5 <- sampling(gp.model,data=data,iter=1000,chains=4)

data <- read.csv(file.name,stringsAsFactors=FALSE)
data <- data[data$X.d. == 2,]
data <- data[sample.indices,]
data <- list(N=dim(data)[1],
             n_iid=data$X.n_i.,
             n_epi=data$X.n_e.,
             y=log(data$value),
             mu=rep(0,dim(data)[1]))

fit.2 <- sampling(gp.model,data=data,iter=1000,chains=4)

data <- read.csv(file.name,stringsAsFactors=FALSE)
data <- data[data$X.d. == 8,]
data <- data[sample.indices,]
data <- list(N=dim(data)[1],
             n_iid=data$X.n_i.,
             n_epi=data$X.n_e.,
             y=log(data$value),
             mu=rep(0,dim(data)[1]))

fit.8 <- sampling(gp.model,data=data,iter=1000,chains=4)

data <- read.csv(file.name,stringsAsFactors=FALSE)
data <- data[data$X.d. == 1000,]
data <- data[sample.indices,]
data <- list(N=dim(data)[1],
             n_iid=data$X.n_i.,
             n_epi=data$X.n_e.,
             y=log(data$value),
             mu=rep(0,dim(data)[1]))

fit.1000 <- sampling(gp.model,data=data,iter=1000,chains=4)

data <- read.csv(file.name,stringsAsFactors=FALSE)
data <- data[data$X.d. == 1000,]
data <- list(N=dim(data)[1],
             n_iid=data$X.n_i.,
             n_epi=data$X.n_e.,
             y=log(data$value),
             mu=rep(0,dim(data)[1]))
vb.1000 <- vb(gp.model,data=data,algorithm="meanfield")


rw.0    <- extract(fit.0,"relative_worth")[[1]]
rw.0.5  <- extract(fit.0.5,"relative_worth")[[1]]
rw.2    <- extract(fit.2,"relative_worth")[[1]]
rw.8    <- extract(fit.8,"relative_worth")[[1]]
rw.1000 <- extract(fit.1000,"relative_worth")[[1]]

par(xpd=TRUE)
historidge(list(rw.0,rw.0.5,rw.2,rw.8,rw.1000),colors=viridis_pal()(5),xlab="relative worth",ylab="")
legend("bottomleft",legend=rev(paste0("d = ",c(0,0.5,2,8,1000))),fill=rev(viridis_pal()(5)),bty="n",border=NA)
