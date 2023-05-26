#' Mediator Plot
#' @param x fit table
#' @export
plot.mrMediator = function(x,...) {
  mrDiagram(fit=x,...)
}

#' Mediator Plot
#' @param fit fit table
#' @param method mediator method either difference or product
#' @importFrom stringr str_split
#' @importFrom dplyr left_join
#' @export
#' @examples
#' \dontrun{
#' fit = mrMediator(data=DTC,
#' exposure = list(name='TSH', type = 'continuous'),
#' mediator = list(name ='SHBG', type = 'continuous'),
#' outcome = list(name='DTC', type ='binary'))
#' mrDiagram(fit =fit, method = 'prod')
#' }
mrDiagram <- function(fit = NULL, method = c('diff', 'prod'),
                      radx = 0.1, rady = 0.04, xmargin = 0.01,
                      box.col = "white", xlim = c(0,1), ylim = NULL) {

  if (is.null(fit)) {
    stop('The estimated table is not found')
  }

  nodes <- nodes
  arrows <- arrows

  method = match.arg(method)

  if (method == 'diff') {
    matchedTable = data.frame(name= c('a','b','c'),
                              labels = c(fit$diff$exposure, fit$diff$mediator, fit$diff$outcome),
                              effect = c(stringr::str_split(fit$diff$unimr.effect, " ", simplify = T)[,1],
                                         stringr::str_split(fit$diff$diff.indirect.effect," ", simplify = T)[,1],
                                         stringr::str_split(fit$diff$mvmr.direct.effect, " ", simplify = T)[,1]),
                              pval = c(fit$diff$pval, fit$diff$indirect.pval, fit$diff$direct.pval))

  } else{
    matchedTable = data.frame(name= c('a','b','c'),
                              labels = c(fit$prod$exposure, fit$prod$mediator, fit$prod$outcome),
                              effect = c(stringr::str_split(fit$diff$unimr.effect, " ", simplify = T)[,1],
                                         stringr::str_split(fit$prod$prod.indirect.effect, " ", simplify = T)[,1],
                                         stringr::str_split(fit$prod$prod.direct.effect, " ", simplify = T)[,1]),
                              pval = c(fit$diff$pval, fit$prod$indirect.pval, fit$prod$direct.pval))
  }

  arrows2 <- dplyr::left_join(arrows, matchedTable, by=c("name"))

  if (is.null(ylim)) {
    ylim <- c(min(nodes$ypos - rady - 0.01), max(nodes$ypos + rady + 0.01))
  }
  if (ylim[1] > 0.2) {
    ylim[1] <- 0.2
  }
  if (ylim[2] < 0.8){
    ylim[2] <- 0.8
  }

  drawMRDiagram(arrows = arrows2, nodes = nodes,
                xmargin = xmargin, radx = radx, rady = rady,
                box.col = box.col, xlim = xlim, ylim = ylim)
}




#' @importFrom diagram openplotmat textrect
#' @export
drawMRDiagram <- function(arrows, nodes,
                          xmargin,radx, rady,
                          box.col = "white",xlim = c(0,1), ylim = c(0,1)) {

  diagram::openplotmat(xlim = xlim, ylim = ylim)
  drawArrows(arrows, nodes, xmargin=xmargin, radx=radx, rady=rady)
  for (i in 1:nrow(nodes)) {
    xpos <- nodes$xpos[i]
    xpos <- adjustxpos(xpos, xmargin, radx)
    mid <- c(xpos, nodes$ypos[i])
    diagram::textrect(mid, radx=radx, rady=rady, lab=arrows$labels[i], box.col=box.col)
  }
}




#' @export
drawArrows <- function(arrows, nodes, xmargin = 0.01, radx = 0.01, rady = 0.04){
  for (i in 1:nrow(arrows)) {
    myarrow2(nodes=nodes, from=arrows$start[i], to=arrows$end[i],
             label=paste0("effect = ",arrows$effect[i]),
             label2 = paste0("p = ", pvalr(arrows$pval[i])) ,xmargin=xmargin,
             radx=radx, rady=rady, label.pos=arrows$labelpos[i],
             arr.pos=arrows$arrpos[i])
  }
}


#' @importFrom diagram straightarrow textplain
#' @export
myarrow <- function(from, to, lwd = 1, adjust = 1, label = "", label2 = '', label.pos = 0.5,
                    radx = 0.10, rady = 0.06, xadj = NULL, yadj = NULL, arr.pos=NULL) {


  mid <- from + label.pos*(to - from)
  xadj1 <- 1
  yadj1 <- -0.3

  if (length(to) > 1) {
    if (from[2] >= to[2]) {
      xadj1 <- 1
    }
    if (from[1] == to[1]) {
      xadj1 <- 1.5
      yadj1 <- 0
    }
  }

  if (is.null(xadj)) xadj <- xadj1
  if (is.null(yadj)) yadj <- yadj1

  diagram::straightarrow(from=from, to=to, lwd=lwd, arr.pos=arr.pos, arr.type="triangle")
  diagram::textplain(mid=mid,lab=label,adj=c(xadj,yadj-1))
  diagram::textplain(mid=mid,lab=label2,adj=c(xadj,yadj+2))
}




#' @export
myarrow2 <- function(nodes, from, to, label="", label2 = '', radx=0.12, rady=0.04, xmargin=0.01,
                     label.pos=0.5, arr.pos=NULL) {

  #from="X";to="Y";label="a"; arr.pos=0.5
  xpos <- nodes$xpos[nodes$name==from]
  xpos <- adjustxpos(xpos,xmargin,radx)
  ypos <- nodes$ypos[nodes$name==from]
  start <- c(xpos, ypos) # from

  xpos <- nodes$xpos[nodes$name==to]
  xpos <- adjustxpos(xpos,xmargin,radx)
  ypos <- nodes$ypos[nodes$name==to]
  end <- c(xpos, ypos) # to

  myarrow(from=start, to=end, label=label, label2 = label2, label.pos=label.pos,
          radx=radx, rady=rady, arr.pos = arr.pos)

}


adjustxpos <- function(xpos, xmargin=0.01, radx=0.12) {
  xspace <- xmargin + 2 * radx
  result <- ifelse(xpos==0.5, 0.5,
                   ifelse(xpos > 0.5, 1-xmargin-radx-(1.0-xpos)*10*xspace,
                          xmargin + radx + (xpos %/% 0.1)*(xmargin+2*radx) +
                            (xpos%%0.1)*10*xspace))
  return(result)
}


pvalr <- function(pvals, sig.limit = .001, digits = 3, html = FALSE) {
  roundr <- function(x, digits = 1) {
    res <- sprintf(paste0('%.', digits, 'f'), x)
    zzz <- paste0('0.', paste(rep('0', digits), collapse = ''))
    res[res == paste0('-', zzz)] <- zzz
    res
  }
  sapply(pvals, function(x, sig.limit) {
    if (x < sig.limit)
      if (html)
        return(sprintf('&lt; %s', format(sig.limit))) else
          return(sprintf('< %s', format(sig.limit)))
    if (x > .1)
      return(roundr(x, digits = 2)) else
        return(roundr(x, digits = digits))
  }, sig.limit = sig.limit)
}




