library(CoDaSeq)
library(ALDEx2)


f = read.delim ("otutab.txt", header=T, row.names= 1, sep="\t")

f <- codaSeq.filter(otu, min.reads=10000, min.prop=0.01, min.occurrence = 0.1, samples.by.row=FALSE)

conds <- read.delim ("condition.txt", header=T, row.names= , sep="\t")

f.x <- aldex.clr(f, as.vector(conds$Region), mc.samples=128, verbose=FALSE, useMC = TRUE)



f.e <- aldex.effect(f.x, as.vector(conds$Region),
                    include.sample.summary=TRUE,
                    verbose=FALSE,
                    useMC = TRUE)

write.csv(f.e,"7f-e-diff.csv")

E.E.clr <- t(f.e[,grep("rab.sample", colnames(f.e))])
rownames(E.E.clr) <- gsub("rab.sample.", "", rownames(E.E.clr))
exp <- apply(E.E.clr, 1, function(x) 2^x)
E.clr <- t(apply(exp, 2, function(x) log2(x) - mean(log2(x))))
write.csv(E.clr,"8-E-clr.csv")


#################
f.t <- aldex.ttest(f.x, as.vector(conds$Region))
low.p <- which(f.t$we.eBH < 0.05)
write.csv(f.t,"10-f-t.csv")


