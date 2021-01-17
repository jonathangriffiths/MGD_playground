#!/bin/R
pal = EmbryoCelltypeColours
df = data.frame(name = names(pal), col = pal)
write.table(df, file="../data/palette.tsv", sep="\t",
    row.names=FALSE, col.names=FALSE, quote=FALSE)
