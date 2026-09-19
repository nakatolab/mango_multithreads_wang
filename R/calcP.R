# Define a function that calculates P-values of interactions
# added `lower.tail=FALSE` option
calcP <- function(v)
{
  P = pbinom(q=v[1]-1,size=v[2],prob=v[3],lower.tail=FALSE)
  return(P)
}
