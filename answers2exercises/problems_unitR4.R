#Answers to R4 Problems in https://r4biochemists.github.io/
#Victor Mateo, 2024

#load genetic code
genetic_code <- read.table("genetic_code.csv", #ALTER: read.csv2
                           sep=";", 
                           header=TRUE,
                           row.names=1)

#fun:
translate_dna <- function(file_name) {
    
    #Read info
    nuc_seq <- readLines(file_name)[-1] #[-1] to remove header (what if actual fa?)
    nuc_n <- nchar(nuc_seq)
    
    #Warning if len is not divisible by 3 
    if (nuc_n%%3 != 0) {
      warning(paste0("Sequence has ", nuc_n," nucleotides and ", nuc_n%/%3, " codons." 
                     ,"Last ", nuc_n%%3," nucleotides will not be used."))
    }
    
    #Translation alg.
    aa_seq <- ""
    for (i in seq(from=1, to=nuc_n ,by=3)) { 
      tmp_codon <- substr(nuc_seq,i,i+2)
      tmp_AA <- genetic_code[tmp_codon,] #vector, colnames lost
      aa_seq <- paste0(aa_seq,tmp_AA)
    }
    return(aa_seq)
    
}
lacz_aa <- translate_dna("lacz.fa")
print(lacz_aa)

pipolb_partial_aa <- translate_dna("piPolB.fasta") #Does not work

#Fix function for actual fasta files
translate_dna_fasta <- function(file_name) {
  
  #Read info
  nuc_seq <- paste(readLines(file_name)[-1], collapse="") #add paste + params and fixed
  nuc_n <- nchar(nuc_seq)

  #Warning if len is not divisible by 3 
  if (nuc_n%%3 != 0) {
    warning(paste("Sequence has ", nuc_n," nucleotides and ", nuc_n%/%3, " codons." 
                  ,"Last ", nuc_n%%3," nucleotides will not be used."),
                  sep = "")
  }
  
  #Translation alg.
  aa_seq <- ""
  for (i in seq(from=1, to=nuc_n ,by=3)) { #nuc_n
    tmp_codon <- substr(nuc_seq,i,i+2)
    tmp_AA <- genetic_code[tmp_codon,] #vector, colnames lost
    aa_seq <- paste0(aa_seq,tmp_AA)
  }
  
  return(list("lenght"=nuc_n, 
              "codons"=nuc_n%/%3, 
              "remainder"=nuc_n%%3,
              "amino_acids"=aa_seq))
  
}

pipolb_partial_aa <- translate_dna_fasta("piPolB_partial.fasta") #Does not work
pipolb_partial_aa

### price_calculator (One argument version) ###

price_calculator <- function(input) {
  
  #Get samples and category from input object/variable/whatever
  #(It must contain "samples" and "category" inside)
  samples <- input$samples
  category <- input$category

  categories <- c(1, 1.15, 2)
  names(categories) = c("normal", "priority", "urgent")
  if (samples < 10) {
    price <- 19 * samples * categories[which(names(categories) ==
                                               category)]
  } else if (samples < 50) {
    price <- 14 * samples * categories[which(names(categories) ==
                                               category)]
  } else if (samples >= 50) {
    price <- 10 * samples * categories[which(names(categories) ==
                                               category)]
  }
  paste("El precio es de", price, "euros.")
}
test_list <- list(samples = 10, category = "normal")
test_df <- data.frame(samples = 10, category = "normal")
price_calculator(test_list)
