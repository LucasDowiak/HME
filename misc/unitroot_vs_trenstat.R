
#########################
# Random Walk vs Trend-Stationary Works pretty well at discreminating between the
# two choices

random_walk <- function(n, init=0)
{
  y <- c(init, rnorm(n-1))
  cumsum(y + 0.2)
}

trend_stationary <- function(n, z, init=z[length(z)])
{
  m <- diff(z[c(1, length(z))]) / diff(c(1, length(z)))
  y <- c(init, vector("numeric", n-1))
  trend <- seq_len(n) + n
  for (jj in (seq_len(n-1) + 1)) {
    y[jj] <- 0.1 + m * trend[jj-1] + 0.2 * y[jj-1] + rnorm(1, sd=1.3)
  }
  return(y)
}

a <- random_walk(200)
aa <- trend_stationary(100, a)
yy <- c(a, aa[-1])

plot(a,  type="l")
tmp5 <- multi_call_search("y ~ .", TRACE=1,
                          DATA=data.frame(y=a), INIT_METHOD="model",
                          DEPTH=1, NEXPERTS=2, TOL=1e-5, FIRST_SPLIT=2, MAXEPOCHS = 300,
                          AR=list(list(1,TRUE,TRUE,FALSE),
                                  list(0,TRUE,FALSE,TRUE)),
                          startmaxit=10, maxepochs=150, maxspeed=25)
plot(tmp5)
tmp5 <- tsHME("y ~ .", TRACE=1,
              DATA=data.frame(y=yy), INIT_METHOD = "model",
              DEPTH=1, NEXPERTS=2, TOL=1e-5, FIRST_SPLIT=2, MAXEPOCHS = 300,
              AR=list(list(1,TRUE,TRUE,FALSE),
                      list(0,TRUE,FALSE,TRUE)), MAXIT=10)
#startmaxit=10, maxepochs=150, maxspeed=25)


