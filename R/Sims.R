#' @title Simulating Recombination
#'
#'
#'

#TODO make this better and then export
#
# require(tidyverse)
# ##---------------------------------------------------
# ##           Sim Tool
# ##---------------------------------------------------
# recombo_blocker  <- function(n = 1e2, # number of loci
#                              CHROMlen = 1e7,
#                              rho = 1e-7, # prob of recombination
#                              k = 1, # number of generations
#                              p1 = 0.25, # AF
#                              p2 = 0.9,
#                              noise = 0.1,
#                              seqerr = 0.125){
#
#   prob = function(d=d,t=t,k=k){1 - exp(-rho * (d) * k)}
#
#   # init
#   ret <- matrix(NA, n, 2)
#   t <- 0
#   pos <- sort(sample(seq(1:CHROMlen), size=n, replace = F))
#
#   # first state
#   ret[1, 1] <- sample( c("M","F"), size = 1, prob = c(0.5, 0.5) )
#
#   # draw subsequent states of recombination block
#   for (i in 2:n) {
#     d <- pos[i] - pos[i-1]
#
#     if(ret[i-1,1] == "M"){
#       ret[i,1] <- sample( c("M","F"), size = 1, prob = c(1-prob(d=d,t=t,k=k), prob(d=d,t=t,k=k)) )
#     } else {
#       ret[i,1] <- sample(c("F","M"), size = 1, prob = c(1-prob(d=d,t=t,k=k), prob(d=d,t=t,k=k)) )
#     }
#
#     if(ret[i-1,1] != ret[i,1]){
#       t <- d
#     }
#
#
#   }
#
#   ret[,2] <- sapply(ret[,1],
#                     function(x){switch(x,
#                                        "M"={
#                                          if(rbinom(n=1,size=1, p = seqerr) == 1){
#                                            runif(n=1, min=0, max=1)
#                                          } else{
#                                            p1 + runif(n=1, min= -noise, max=noise) # some noise
#                                          }
#
#                                        },
#
#                                        "F"={
#                                          if(rbinom(n=1,size=1, p = seqerr) == 1){
#                                            runif(n=1, min=0, max=1)
#                                          } else{
#                                            p2 + runif(n=1, min= -noise, max=noise) # some noise
#                                          }
#
#                                        }
#                     )
#                     }, simplify = T
#   )
#   ret <- tibble::tibble(POS = pos, parent = ret[,1], AF = as.numeric(ret[,2]))
#
#
#   # end of overall function
#   return(ret)
# }
#
#
# ##---------------------------------------------------
# ##           Run Code
# ##---------------------------------------------------
#
# p1 <- 0.5
# p2 <- 0.9
#
# ret <- recombo_blocker(n = 5e3, # number of loci
#                        rho = 1e-7, # prob of recombination
#                        CHROMlen = 1e7, # length of chromosome (directly interacts with rho...obviously)
#                        p1 = p1, # AF for M
#                        p2 = p2, # AF for F
#                        seqerr = 0.1, # frequency that we hit a sequencing error
#                        noise = 0.1 # error that can occur around true AF for M and F
# )
#
# # find recombo breakpoints
# brkpts <- ret %>%
#   dplyr::mutate(leadpar = dplyr::lead(parent)) %>%
#   dplyr::filter(parent != leadpar) %>%
#   dplyr::select(POS) %>%
#   cbind(., fct = letters[1:nrow(.)]) %>%
#   dplyr::mutate(fct = factor(fct))
#
# blockdf <- ret %>%
#   dplyr::left_join(x=., y=brkpts, by=c("POS")) %>%
#   tidyr::fill(fct, .direction = c("up")) %>%
#   dplyr::mutate(fct = ifelse(is.na(fct), "end", fct)) %>%
#   dplyr::group_by(fct) %>%
#   summarise(parent = unique(parent), minpos = min(POS), maxpos=max(POS))
#
#
#
#
# PlotObj <- ggplot() +
#   geom_point(data=ret, aes(x=POS, y=AF, color=parent)) +
#   geom_rect(data=blockdf, aes(xmin = minpos, xmax=maxpos, ymin=1.01, ymax=1.1, fill=parent)) +
#   geom_hline(yintercept = p1, color="#de2d26", linetype="dashed") +
#   geom_hline(yintercept = p2, color="#de2d26", linetype="dashed") +
#   scale_color_manual("Parental Strain", values = brewer.pal(n = length(levels(factor(ret$parent))), name = "Set2")) +
#   scale_fill_manual("Parental Strain", values = brewer.pal(n = length(levels(factor(ret$parent))), name = "Set2")) +
#   ggtitle("Simulation of Recombination with \n Error in AF like Pv Cross") +
#   ylab("Allele Frequency") + xlab("POS") +
#   theme(plot.title = element_text(hjust = 0.5, family="Arial", size=15, face = "bold"),
#         axis.text.x = element_blank(),
#         axis.title.x = element_text(family="Arial", size=13, face = "bold"),
#         axis.text.y = element_text(family="Arial", size=11, face = "bold", hjust=0.5),
#         axis.title.y = element_text(family="Arial", size=13, face = "bold"),
#         legend.title = element_text(family="Arial", size=10, face = "bold"),
#         legend.text = element_text(family="Arial", size=8, face = "bold"),
#         legend.title.align = 0.5,
#         panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.line = element_line(colour = "black", size=2),
#         axis.ticks = element_blank()
#   )
#
# plot(PlotObj)
#
#
#
#
# ##---------------------------------------------------
# ##           Testing RunMed
# ##---------------------------------------------------
#
#
# ret <- ret %>%
#   dplyr::mutate(AF_smoothed = stats::runmed(x=AF, k=41, endrule=c("median"), algorithm=c("Turlach")))
#
#
#
#
# SmoothPlotObj <- ggplot() +
#   geom_point(data=ret, aes(x=POS, y=AF_smoothed, color=parent)) +
#   geom_rect(data=blockdf, aes(xmin = minpos, xmax=maxpos, ymin=1.01, ymax=1.1, fill=parent)) +
#   geom_hline(yintercept = p1, color="#de2d26", linetype="dashed") +
#   geom_hline(yintercept = p2, color="#de2d26", linetype="dashed") +
#   scale_color_manual("Parental Strain", values = brewer.pal(n = length(levels(factor(ret$parent))), name = "Set2")) +
#   scale_fill_manual("Parental Strain", values = brewer.pal(n = length(levels(factor(ret$parent))), name = "Set2")) +
#   ggtitle("Simulation of Recombination with \n Error in AF like Pv Cross") +
#   ylab("Allele Frequency") + xlab("POS") +
#   theme(plot.title = element_text(hjust = 0.5, family="Arial", size=15, face = "bold"),
#         axis.text.x = element_blank(),
#         axis.title.x = element_text(family="Arial", size=13, face = "bold"),
#         axis.text.y = element_text(family="Arial", size=11, face = "bold", hjust=0.5),
#         axis.title.y = element_text(family="Arial", size=13, face = "bold"),
#         legend.title = element_text(family="Arial", size=10, face = "bold"),
#         legend.text = element_text(family="Arial", size=8, face = "bold"),
#         legend.title.align = 0.5,
#         panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.line = element_line(colour = "black", size=2),
#         axis.ticks = element_blank()
#   )
#
# plot(SmoothPlotObj)
#
#
#
# # That's very nice...
#
#
