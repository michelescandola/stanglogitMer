#' Internal function to find maximum and minimum of the dependent variable divided by cell.
#'
#' \code{findMaxMin} returns a data frame with the same number of rows
#' of the original data frame with two columns: the minimum and
#' the maximum value obtained by the dependent variable
#' grouped per cell.
#'
#' @param formula the formula to split the dependent variable into cells.
#' @param data the original data frame.

.findMaxMin <- function(formula, data, dependent)  {

  tmp.mat <- model.frame(formula,data)
  tmp.int = interaction(tmp.mat)
  mm      <- list()

  for(i in 1:nrow(tmp.mat)){
    sel <- which(tmp.int[i]==tmp.int)
    mm[[i]] = data.frame(min = min(dependent[sel]),
                         max = max(dependent[sel]))
  }
  mm = do.call("rbind",mm)

  return(mm)
}

#' Internal function to write formulas from a vector of individual strings.
#'
#' @param x is the vector of individuals strings. The first string is the dependent variable.
#' @examples
#' paste.formulas(c("y","Condition","Group"))

.paste.formulas = function(x){
  out=NULL
  for(i in 1:(length(x)-1)) out = paste0(out,x[i],"*")
  out = paste0("~",out,x[length(x)])
  return(out)
}
