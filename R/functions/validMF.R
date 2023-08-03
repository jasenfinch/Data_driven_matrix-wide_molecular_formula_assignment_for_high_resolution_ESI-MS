
validMF <- function(MF){
  grepl("^(([CHONPS]*)(([0-9])*))+$",MF) 
}
