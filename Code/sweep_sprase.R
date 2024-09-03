sweep_sparse <- function(x, margin, stats, fun = "/") {
        f <- match.fun(fun)
        if (margin == 1) {
                idx <- x@i + 1
        } else {
                idx <- x@j + 1
        }
        x@x <- f(x@x, stats[idx])
        return(x)
}
