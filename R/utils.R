#' lme4-style formulas
#' 
#' Transform the separated fixed and random formulas
#' into the single formula with lme4 syntaxis
#' @export
lme4_fml <- function(fix, rnd, rm_int = TRUE) {
  rnd.terms <- attr(terms(rnd), 'term.labels')
  rnd.terms.lme4 <- paste('(1|', rnd.terms, ')', sep ='')
  int <- ifelse(rm_int, '-1', '')
  rnd.upd <- paste('~ .', int, paste('+', rnd.terms.lme4, collapse = ' '))
  fml.lme4 <- update(fix, rnd.upd)
  return(fml.lme4)
}
