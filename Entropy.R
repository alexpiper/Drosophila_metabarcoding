

#Script defines function to find relative entropy of "data.set"
#load seqinr package
library(seqinr)

seqs <- seqinr::read.fasta("reference/fwh_insecta_aligned_curated.fasta")

data.set <- seqs[1:100]



#Supply fastas to be analyzed in form of list or vector (read.fasta() from seqinr).
# Define desired window size (default = 500). If quick plot is desired, set plot = T. 

#data.set: A list or vector of a gene sequence as read in as a fasta.
#window: set size of window over which to calculate relative entropy.
#slide: define how to slide the window. If completely distinct windows are desired set: slide = window size.
#plot: setting plot = True generates a simple plot in addition to vector output
EntropySlide <- function(data.set, window = 500, slide = 100, plot = F){
  #change list to a vector if necessary
  if(!(class(data.set) %in% c("list","vector"))){
    stop("Function requires a list or vector")
  }
  if(class(data.set) == "list") data.set <- unlist(data.set)
  #don't want slide bigger than window
  if(slide > window) stop("Slide must be less than size of window!")
  #set non-nucleotide reads to NA and drop
  data.set[!(data.set %in% c("a", "g", "c", "t"))] <- NA
  data.set <- na.omit(data.set)
  #calculate proportions of nucleotides for entire genome
  prop.a <- length(which(data.set == "a")) / length(data.set)
  prop.g <- length(which(data.set == "g")) / length(data.set)
  prop.c <- length(which(data.set == "c")) / length(data.set)
  prop.t <- length(which(data.set == "t")) / length(data.set)
  list.props <- list(a = prop.a, g = prop.g, c = prop.c, t = prop.t)
  #calculate number of element needed in each vector (total number of windows)
  size.vectors <- ceiling((length(data.set) - window) / slide)
  window.margins <- rep(0, size.vectors)
  window.margins <- cbind(window.margins, window.margins)
  window.margins[,1] <- seq(1, by = slide, length.out = size.vectors)
  window.margins[,2] <- seq(window + 1, by = slide, length.out = size.vectors)
  #empty list to hold entropy data by base pair, section
  window.list <- list(a = rep(0, size.vectors), g = rep(0, size.vectors), c = rep(0, size.vectors), t = rep(0, size.vectors))
  pb <- txtProgressBar(min = 0, max = length(window.margins[,1]), initial = 0, char = "=", width = NA, title, label, style = 1)
  for(i in 1:length(window.margins[,1])){
    window.list$a[i] <- length(which(data.set[window.margins[i,1]:window.margins[i,2]] == "a")) / window
    window.list$g[i] <- length(which(data.set[window.margins[i,1]:window.margins[i,2]] == "g")) / window
    window.list$c[i] <- length(which(data.set[window.margins[i,1]:window.margins[i,2]] == "c")) / window
    window.list$t[i] <- length(which(data.set[window.margins[i,1]:window.margins[i,2]] == "t")) / window
    setTxtProgressBar(pb, i)
  }
  close(pb)
  #set new list to hold entropy values
  list.entropy <- list(a = rep(0,size.vectors), g = rep(0,size.vectors), c = rep(0,size.vectors), t = rep(0,size.vectors))
  #vectorized entropy function:
  list.entropy$a <- window.list$a * log(window.list$a / prop.a)
  list.entropy$g <- window.list$g * log(window.list$g / prop.g)
  list.entropy$c <- window.list$c * log(window.list$c / prop.c)
  list.entropy$t <- window.list$t * log(window.list$t / prop.t)
  #entropy for each base pair must be summed to give total window entropy
  entropy <- list.entropy$a + list.entropy$g + list.entropy$c + list.entropy$t
  if(plot == T) {
    tick.positions <- c(1, ceiling(length(window.margins[,1]) * (1 / 4)), ceiling(length(window.margins[,1]) * (1/2)), ceiling(length(window.margins[,1]) * (3 / 4)), length(window.margins[,1]))
    tick.names <- as.character(window.margins[tick.positions])
    plot(entropy, xlab = "Position", ylab="Relative Entropy", xaxt = "n", type='l')
    axis(side = 1, at = tick.positions, labels = tick.names)
  }
  return(entropy)
}