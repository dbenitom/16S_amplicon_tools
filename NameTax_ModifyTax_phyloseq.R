# Modifications for taxonomic tables in phyloseq objects
require(phyloseq)

# Taxonomy table of the phyloseq object ps
tax.tab <- data.frame(tax_table(ps))

# NameTax function. It adds a prefix in front of each taxonomical level. Modified from https://github.com/joey711/phyloseq/issues/773
NameTax <- function(x, ind){
  if(is.na(x[ind])){
    x[ind] <- x[ind]
  } else {
    if(ind==1){x[ind] <- paste("d", x[ind], sep="_")} else{             # Domain
      if(ind==2){x[ind] <- paste("p", x[ind], sep="_")} else{           # Phylum
        if(ind==3){x[ind] <- paste("c", x[ind], sep="_")} else{         # Class
          if(ind==4){x[ind] <- paste("o", x[ind], sep="_")} else{       # Order
            if(ind==5){x[ind] <- paste("f", x[ind], sep="_")} else{     # Family
              if(ind==6){x[ind] <- paste("g", x[ind], sep="_")} else{   # Genus
                if(ind==7){x[ind] <- paste("s", x[ind], sep="_")}       # Species. Still needs to be polished to merge genus and species names in the single binomial nomenclature. 
              }
            }
          }
        }
      }
    }
  }
}

# ModifyTax function. It replaces NA levels with the lowest taxonomic rank that could be classified.
ModifyTax <- function(x,ind){
  #   xth row in the dataframe
  #   ind taxonomy level to change
  if(is.na(x[ind])){
    nonNa <- which(!is.na(x[-ind])) # which taxa are not NA excepting the one we're interested in.
    maxNonNa <- max(nonNa)
    x[ind] <- x[maxNonNa]
  }else{x[ind] <- x[ind]}
}

# Apply the functions. The order is important, first the NameTax function and second the ModifyTax function.
for (i in 1:7) {
  tax_table(ps)[,i] <- apply(tax.tab,1,NameTax,ind=i)
}

for (i in 1:7) {
  tax_table(dada2_phyloseq_trial)[,i] <- apply(tax.tab,1,ModifyTax,ind=i)
}
