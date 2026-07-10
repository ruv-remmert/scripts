#!/usr/bin/env Rscript

pathSTD <- getwd()

mito     <- data.frame(chr="Mito", start=0,      end=85779)
rDNA     <- data.frame(chr="XII",  start=451000,  end=469500)
irregular  <- data.frame(chr="X",   start=472500,  end=484000)
irregular2 <- data.frame(chr="VIII",start=85200,   end=95160)
irregular3 <- data.frame(chr="VIII",start=212300,  end=216200)
irregular4 <- data.frame(chr="VII", start=0,       end=8900)
irregular5 <- data.frame(chr="VIII",start=531400,  end=553500)
irregular6 <- data.frame(chr="II",  start=29700,   end=36200)

blacklist <- rbind(mito, rDNA, irregular,
                   irregular2, irregular3, irregular4,
                   irregular5, irregular6)

write.table(blacklist, row.names = F, col.names = F,
          file="data/blacklist.bed", quote = F, sep="\t")
