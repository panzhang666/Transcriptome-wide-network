library(glasso)
library(preprocessCore)
library("ppcor")
library(corpcor)
library("MBESS")
library(network)
library(QUIC)
library("ROCR")
library("hash")

setwd("/Users/panzhang/Desktop/AS/TWN")
tee=read.csv("filtered_TG_geneTPM.csv", row.names = 1)
tee_sd <- apply(tee, 1, sd)
tee <- cbind(tee,tee_sd)
#Sortbysd Select top 1000 genes
tee <- tee[order(tee[,"tee_sd"],decreasing=TRUE),]
tee <- tee[1:2000,-21]
tee <- t(tee)
scale_te=apply(tee,2,scale)

#scale_te <- t(scale_te)

#new_te=normalize.quantiles(scale_te)
new_te = scale_te
colnames(new_te) <- colnames(tee)
rownames(new_te) <- rownames(tee)
colnames(new_te) <- toupper(colnames(new_te))
#write.table(new_te,file="NA_TE_100_twn.txt",sep = "\t",quote = FALSE)



irr=read.csv("filtered_TG_isoRatio.csv", row.names = 1)
#geneName = as.data.frame(table(tee$new_id))
irr_sd <- apply(irr, 1, sd)
irr <- cbind(irr,irr_sd)
#Sortbysd Select top 2000 isoforms
irr <- irr[order(irr[,"irr_sd"],decreasing=TRUE),]
irr <- irr[1:3000,-21]
irr = t(irr)
scale_ir=apply(irr,2,scale)

new_ir = scale_ir
#new_ir=normalize.quantiles(scale_ir)
colnames(new_ir) <- colnames(irr)
rownames(new_ir) <- rownames(irr)

#write.table(new_ir,file="NA_IR_200_twn.txt",sep = "\t",quote = FALSE)


geneAnnot <- read.csv("geneAnnot.csv")
transcriptAnnot <- read.csv("transcriptAnnot.csv")

penalty_parameter <- function(te, ir, transcriptAnnot, p.te, p.ir, p.te.ir, p.precision=0.05){
  n.te <- dim(te)[2]
  n.ir <- dim(ir)[2]
  feature <- append(colnames(te), colnames(ir))
  n.total = n.te + n.ir
  r <- matrix(0, nrow = n.total , ncol = n.total)
  #store transcript annotation information
  h <- hash(transcriptAnnot$transcript_id, transcriptAnnot$gene_id)
  for (i in seq(2, n.total)){
    for (j in seq(i-1)){
      if ((i <= n.te) & (j <= n.te)){
        r[i,j] = p.te
      }
      else if ((i > n.te) & (j > n.te)){
        if (h[[feature[i]]] == h[[feature[j]]]){
          r[i,j] = p.precision
        }
        else{
          r[i,j] = p.ir
        }
      }
      else{
        if(h[[feature[i]]] == feature[j]){
          r[i,j] = p.precision
        }
        else{
          r[i,j] = p.te.ir
        }
      }
    }
  }
  colnames(r) <- feature
  rownames(r) <- feature
  return(r)
}

systematric <- function(x){
  ncx <- ncol(x)
  for (i in seq(ncx-1)){
    for (j in seq(i+1,ncx)){
      x[i,j] <- x[j,i]
    }
  }
  return(x)
}

te<- new_te
ir <- new_ir
te.ir <- cbind(te,ir)
S <- cov(te.ir)
min(S)
max(S)


rho.QUIC <- penalty_parameter(new_te, new_ir, transcriptAnnot, 0.8, 0.65, 0.6)
rho<- systematric(rho.QUIC)
isSymmetric(rho)
#change rho
#rho[which(rho == 0.8)] <- 0.85
#rho[which(rho == 0.75)] <- 0.8
#rho[which(rho == 0.6)] <- 0.65
#table(rho)


#colnames(rho) <- append(colnames(new_te), colnames(new_ir))
#rownames(rho) <- append(colnames(new_te), colnames(new_ir))



rho_twn_run <- function(te, ir,transcriptAnnot, rho){
  te.ir <- cbind(te,ir)
  S <- cov(te.ir)
  a <- glasso(S,rho)
  P <- a$wi
  A <- ifelse(P!=0 & row(P)!=col(P),1,0)
  n.te <- ncol(te)
  n.ir <- ncol(ir)
  n.total = n.te + n.ir
  genes <- colnames(te)
  isoforms <- colnames(ir)
  genes.isoforms <- colnames(te.ir)
  h <- hash(transcriptAnnot$transcript_id, transcriptAnnot$gene_id)
  edge <- c("node1","node2","type")
  for (i in seq(2,n.total)) {
    for (j in seq(i-1)) {
      if (A[i,j] == 1){
        if ((i <= n.te) & (j <= n.te)){
          edge <- rbind(edge, c(genes.isoforms[i],genes.isoforms[j],1))
        }
        else if ((i > n.te) & (j > n.te)) {
          if (h[[genes.isoforms[i]]] != h[[genes.isoforms[j]]] ){
            edge <- rbind(edge, c(genes.isoforms[i],genes.isoforms[j],3))
          }
        }
        else{
          if (h[[genes.isoforms[i]]] != genes.isoforms[j]){
            edge <- rbind(edge, c(genes.isoforms[i],genes.isoforms[j],2))
          }
        }
      }
    }
  }
  out <- as.data.frame(edge)
  colnames(P) <- genes.isoforms
  rownames(P) <- genes.isoforms
  return(list(out,P))  
}

results <- rho_twn_run(new_te, new_ir, transcriptAnnot, rho)
p <- as.data.frame(results[2])
twn.edge <- as.data.frame(results[1])
colnames(twn.edge) <-  c("node1","node2","type")
twn.edge <- twn.edge[-1,]
table(twn.edge[,"type"])
scale_free <- function(twn.edge){
  te_nodes<-append(as.character(twn.edge$node1),as.character(twn.edge$node2))
  node_degree <- as.data.frame(table(te_nodes))
  degree_freq <- as.data.frame(table(node_degree$Freq))
  k <- as.numeric(as.character(degree_freq[,1]))
  pk <- as.numeric(as.character(degree_freq[,2]))
  pk <- pk/sum(pk)
  cor.v <- cor(log(k), log(pk))
  return(cor.v )
}
scale_free(twn.edge)  
te.te <-  twn.edge[(twn.edge$type == 1),]
te.ir <-  twn.edge[(twn.edge$type == 2),]
ir.ir <-  twn.edge[(twn.edge$type == 3),]
scale_free(te.te)
scale_free(te.ir)
scale_free(ir.ir)
write.table(twn.edge,file="20sample_TG_TPM_2000_3000_0.8_0.65_0.65.txt",sep = "\t",quote = FALSE)