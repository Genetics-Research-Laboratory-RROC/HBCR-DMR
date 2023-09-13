t.test.mod <- function(x,y){
    x <- x[!(is.na(x))]
    y <- y[!(is.na(y))]
    if (length(x)==1) x <- rep(x,2)
    if (length(y)==1) y <- rep(y,2)
    if (length(unique(c(x,y))) == 1)  return(1)
    if (length(x)==0 | length(y)==0)  return(1)
    if (length(unique(x)) == 1 & length(unique(y)) == 1){
        x[1] <- x[1] + 1e-10
    }
    return(t.test(x,y)$p.value)
}