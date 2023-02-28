library("entropy")
library("ggplot2")

set.seed(123)
p <- 0.5
K <- 1000
N <- 12
length_low <- 10
length_up <- 1000
step <- 10

bernoulli_entropy <- -p * log2(p) - (1 - p) * log2(1 - p)
geom_entropy <- -log2(p) - (1 - p) * log2(1 - p) / p

estimate_shannon_entropy <- function(dna_sequence) {
  bases <- table(strsplit(dna_sequence, "")[[1]])
  dist <- bases / sum(bases)
  return(entropy(dist, base=2))
}

inbuild_bern <- c()
inbuild_geom <- c()
entr_deviation_bern  <- c()
entr_deviation_geom  <- c()

for (length in seq(length_low, length_up, step)) {
  for (i in 1:K) {
    bern <- rbinom(length, 1, p)
    geom <- rgeom(length, p)
    
    entr_deviation_bern <- c(entr_deviation_bern, bernoulli_entropy - estimate_shannon_entropy(paste0(bern)))
    entr_deviation_geom <- c(entr_deviation_geom, geom_entropy - estimate_shannon_entropy(paste0(geom)))
  }
  
  inbuild_bern <- c(inbuild_bern, var(entr_deviation_bern))
  inbuild_geom <- c(inbuild_geom, var(entr_deviation_geom))
}

df_bern <- data.frame(length=seq(length_low, length_up, step), entropy=inbuild_bern, distribution="Bernoulli")
df_geom <- data.frame(length=seq(length_low, length_up, step), entropy=inbuild_geom, distribution="Geometric")
df <- rbind(df_bern, df_geom)

ggplot(df, aes(x=length, y=entropy, color=distribution)) + geom_line() + facet_wrap(~distribution, ncol=2) + ggtitle("Variance of Shannon Entropy")

generate_bernoulli <- function(p, n) {
  rand_arr <- runif(n)
  return(ifelse(rand_arr >= 1 - p, 1, 0))
}

generate_geometric <- function(p, n) {
  geom_arr <- c()
  
  while (length(geom_arr) < n) {
    k <- 1
    while (TRUE) {
      rand_arr <- ifelse(runif(1) >= 1 - p, 1, 0)
      if (rand_arr == 1) {
        geom_arr <- c(geom_arr, k)
        break
      } else {
        k <- k + 1
      }
    }
  }
  
  return(geom_arr)
}

brv_bern <- c()
brv_geom <- c()
entr_deviation_bern  <- c()
entr_deviation_geom  <- c()

for (length in seq(length_low, length_up, step)) {
  entr_deviation_bern  <- c()
  entr_deviation_geom  <- c()
  
  for (i in 1:K) {
    bern <- generate_bernoulli(p, length)
    geom <- generate_geometric(p, length)

    entr_deviation_bern <- c(entr_deviation_bern, bernoulli_entropy - estimate_shannon_entropy(paste0(bern)))
    entr_deviation_geom <- c(entr_deviation_geom, geom_entropy - estimate_shannon_entropy(paste0(geom)))
  }

