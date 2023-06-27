############################
### LIKELIHOOD FUNCTIONS ###
############################
require(mvtnorm)
require(graphics)
require(MASS)
require(mixAK)
require(MCMCpack)
require(plotrix)
require(fields)
require(matrixcalc)
require(nlme)
require(gdata)
require(grid)

# FUNCOES
distancia <- function(x1,y1,x2,y2) {
  dist = sqrt((x1-x2)^2+(y1-y2)^2)
  return(dist)
}

matriz.cheia <- function(sigma,rho,H,phi){
    p <- dim(sigma)[1]
    n <- dim(H)[1]
    MAT.COV <- matrix(NA, ncol = n*p, nrow=n*p)
    sapply(1:p, 
         function(i){
           sapply(i:p, 
                  function(j){
                    # Basta gerar sigma pelo produto de sigma_i vezes sigma_j
                    # Constroi os 4 blocos
                    MAT.COV[(((i-1)*n)+1):(i*n),(((j-1)*n)+1):(j*n)] <<- sigma[i,j]*rho[i,j]*((1+(H/phi[i,j]))^(-1))
                    # Repete p banda inferior
                    MAT.COV[lower.tri(MAT.COV)] <<- t(MAT.COV)[lower.tri(MAT.COV)]
                  })
         })
    return(MAT.COV)}

sub.matriz <- function(M,r,c){
  rg <- (row(M)-1)%/%r+1
  cg <- (col(M)-1)%/%c+1
  rci <- (rg-1)*max(cg) + cg
  N <- prod(dim(M))/r/c
  cv <- unlist(lapply(1:N,function(x) M[rci==x]))
  dim(cv) <- c(r,c,N)
  cv
} 

dmvnorm.chol <- function(x,mean=rep(0,(n*p)),A=diag(p),R=diag(n),log=FALSE){

    if (is.vector(x)) {
    x <- matrix(x, ncol = length(x))
  }
  dec.1 <- tryCatch(chol(A), error = function(e) e)
  dec.2 <- tryCatch(chol(R), error = function(e) e)
  
  tmp <- backsolve(kronecker(dec.1,dec.2), t(x) - mean, transpose = TRUE)
  rss <- colSums(tmp^2)
  logretval <- -n*sum(log(diag(dec.1)))-p*sum(log(diag(dec.2)))-0.5*n*p*log(2*pi)-0.5*rss
  
  names(logretval) <- rownames(x)
  if (log) 
    logretval
  else exp(logretval)
}

vero.aprox= function(mat_cheia,n,y_vet,mu_vetor){
  VERO.APROX.CH.DMV <- NULL
  v <- sub.matriz(mat_cheia,n,n)
  # MATRIZ PERMUTADA
  R_SIGMA <- matrix(c(v),ncol=(p^2))
  # DECOMPOSICAO EM VALORES SINGULARES
  mat.sing <- svd(R_SIGMA)
  vec.R <- sqrt(mat.sing$d[1])*mat.sing$u[,1]
  vec.A <- sqrt(mat.sing$d[1])*mat.sing$v[,1]
  # DA DECOMPOSICAO EXTRAIMOS A E R
  A <- matrix(vec.A,ncol=p)
  R <- matrix(vec.R,ncol=n)
  VERO.APROX.CH.DMV <- dmvnorm.chol(y_vet,mu_vetor,-A,-R,log=TRUE)
  return(VERO.APROX.CH.DMV)
}

ler_arquivo_csv = function(nome, cenario_var, nseq_usado = 250){
  #if (cenario_var == 1){
    nome_arquivo = paste0(nome,"_20220910_cenario_",cenario_var,"nseq_",nseq_usado,".csv")
    return(as.matrix(read.csv(nome_arquivo, sep = ";", dec = ",")[,-1]))
  #} else{
 #   print(NULL)
 # }
}

escolha_cenario = function(x){
  cenario <<- x
  if (x == 1) {
    phi.11 <<- 0.1
    phi.22 <<- 0.1
  } else if (x == 2) {
    phi.11 <<- 0.1
    phi.22 <<- 0.2
  }else if (x == 3) {
    phi.11 <<- 0.1
    phi.22 <<- 0.5
  }else if (x == 4) {
    phi.11 <<- 0.1
    phi.22 <<- 0.7
  } else{
    print('Insira valor de 1 a 4')
  }
  
}

cenario = 0

# Numero de variaveis
p <- 2
# Numero de localizacoes
n <- 100

escolha_cenario(1)
#### LAT, LONG, H ####

lat_long = read.csv2(paste0("Localizacoes_contorno_20220910_cenario_",cenario,"n_",n,".csv"))

lat = lat_long[,2]
long = lat_long[,3]

plot(long,lat,pch=16)


## Criando a Matriz de Distâncias H
# ij que varia de acordo com a componente
# h varia: matriz das distancias do espaço
H<- matrix(0,nrow=n,ncol=n)

sapply(1:(n-1), 
       function(k){
         sapply((k+1):n, 
                function(l){
                  H[k,l] <<- distancia(long[k],lat[k],long[l],lat[l])
                  H[l,k] <<- H[k,l]
                })
       })
dim(H)

# Media usada foi 0: latitude e longitude nao estao entrando como covariaveis
MU.vetor <- rep(0,(n*p))
### PHI, RHO, SIGMA ####
#Cenario
phi.12 <- (phi.11+phi.22)/2
# Na funcao do artigo, ha a restricao do c12: sera adaptado

# delta <- matrix(rep(0,(p^2)),ncol=p)
phi.vet <- c(phi.11,phi.12,phi.22)

# tentar fazer sem FOR
phi <- matrix(rep(0,(p^2)),ncol=p)
phi[lower.tri(phi,diag=TRUE)] <- phi.vet
# Tem que usar o transposto pq havia motivo para usa-lo
phi[upper.tri(phi)] <- t(phi)[upper.tri(phi)]
# phi[upper.tri(phi)] <- phi[lower.tri(phi)]

# b_ij no artigo = phi
# phi <- 0.3
# restricao correlacao: rho.12 = phi.11*phi.22/(phi.12^2)
rho.12 <- 0.4 # a mesma corr para todos os cenários
# rho.13 <- 0.6 
# rho.23 <- 0.1

rho.vet <- c(rho.12) #,rho.13, rho.23) 

# tentar fazer sem FOR
rho <- matrix(rep(1,(p^2)),ncol=p)
rho[lower.tri(rho)] <- rho.vet
rho[upper.tri(rho)] <- t(rho)[upper.tri(rho)]

# SIGMA
sigma.1 <- 1
sigma.2 <- 1.5

# Refazer esse produto cruzado de forma generalizada
sigma.vet <- c(sigma.1^2,sigma.1*sigma.2,sigma.2^2)

# tentar fazer sem FOR
sigma <- matrix(rep(NA,(p^2)),ncol=p)
sigma[lower.tri(sigma, diag=TRUE)] <- sigma.vet
sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(sigma)]


#### Análise da Verossimilhança ####
L <- 1 #(replicas)

# Fixar todos os params e variar 1 de cada vez (inicialmente)
n.seq <- 250# fazer 500 para os demais, 100 é para contorno

#### SEQUENCIAS - cenario 1 ####

seq.rho.12 <- seq(0.05,.5,l=n.seq)

seq.phi.11 <- seq(0.015,0.4,l=n.seq)
seq.phi.22 <- seq(0.015,0.4,l=n.seq)

seq.s.1 <- seq(.05,1.8,l=n.seq)
seq.s.2 <- seq(.05,4.5,l=n.seq) #sigmas

#### SEQUENCIAS - cenario 2 ####

seq.rho.12 <- seq(0.05,.5,l=n.seq)

seq.phi.11 <- seq(0.015,0.4,l=n.seq)
seq.phi.22 <- seq(0.015,0.4,l=n.seq)

seq.s.1 <- seq(.05,1.8,l=n.seq)
seq.s.2 <- seq(.05,4.5,l=n.seq) #sigmas

#### SEQUENCIAS - cenario 3 ####

seq.rho.12 <- seq(0.05,.5,l=n.seq)

seq.phi.11 <- seq(0.015,0.4,l=n.seq)
seq.phi.22 <- seq(0.015,0.4,l=n.seq)

seq.s.1 <- seq(.05,1.8,l=n.seq)
seq.s.2 <- seq(.05,4.5,l=n.seq) #sigmas

#### SEQUENCIAS - cenario 4 ####

seq.rho.12 <- seq(0.05,.5,l=n.seq)

seq.phi.11 <- seq(0.015,0.4,l=n.seq)
seq.phi.22 <- seq(0.015,.9,l=n.seq)

seq.s.1 <- seq(.05,1.8,l=n.seq)
seq.s.2 <- seq(.05,5,l=n.seq) #sigmas


#### AGRUPANDO DADOS ####
setwd("C:/Users/Caroline/Dropbox/Projeto_AproxSep/Scripts/Caroline/IC_TCC/S4_contorno/S4_contorno_dados_gerados")

vero.mat.phi.11.sigma.1 <- ler_arquivo_csv("vero.mat.phi.11.sigma.1", cenario)
vero.mat.phi.11.sigma.2 <- ler_arquivo_csv("vero.mat.phi.11.sigma.2", cenario)
vero.mat.phi.11.rho.12 <-  ler_arquivo_csv('vero.mat.phi.11.rho.12', cenario)
vero.mat.phi.22.sigma.1<-  ler_arquivo_csv('vero.mat.phi.22.sigma.1', cenario)
vero.mat.phi.22.sigma.2 <-  ler_arquivo_csv('vero.mat.phi.22.sigma.2', cenario)
vero.mat.phi.22.rho.12 <-  ler_arquivo_csv('vero.mat.phi.22.rho.12', cenario)
vero.mat.sigma.1.sigma.2 <-  ler_arquivo_csv('vero.mat.sigma.1.sigma.2', cenario)
vero.mat.sigma.1.rho.12<- ler_arquivo_csv('vero.mat.sigma.1.rho.12', cenario)
vero.mat.sigma.2.rho.12<-  ler_arquivo_csv('vero.mat.sigma.2.rho.12', cenario)

aprox.mat.phi.11.sigma.1 <- ler_arquivo_csv('aprox.mat.phi.11.sigma.1', cenario)
aprox.mat.phi.11.sigma.2 <- ler_arquivo_csv('aprox.mat.phi.11.sigma.2', cenario)
aprox.mat.phi.11.rho.12 <- ler_arquivo_csv('aprox.mat.phi.11.rho.12', cenario)
aprox.mat.phi.22.sigma.1<- ler_arquivo_csv('aprox.mat.phi.22.sigma.1', cenario)
aprox.mat.phi.22.sigma.2 <- ler_arquivo_csv('aprox.mat.phi.22.sigma.2', cenario)
aprox.mat.phi.22.rho.12 <- ler_arquivo_csv('aprox.mat.phi.22.rho.12', cenario)
aprox.mat.sigma.1.sigma.2 <- ler_arquivo_csv('aprox.mat.sigma.1.sigma.2', cenario)
aprox.mat.sigma.1.rho.12<- ler_arquivo_csv('aprox.mat.sigma.1.rho.12', cenario)
aprox.mat.sigma.2.rho.12<- ler_arquivo_csv('aprox.mat.sigma.2.rho.12', cenario)



m.phi.11.sigma.1 <- max(vero.mat.phi.11.sigma.1)
aprox.mat.m.phi.11.sigma.1 <- max(aprox.mat.phi.11.sigma.1)

m.phi.11.sigma.2 <- max(vero.mat.phi.11.sigma.2) 
aprox.mat.m.phi.11.sigma.2 <- max(aprox.mat.phi.11.sigma.2)

m.phi.11.rho.12 <- max(vero.mat.phi.11.rho.12) 
aprox.mat.m.phi.11.rho.12 <- max(aprox.mat.phi.11.rho.12)

m.phi.22.sigma.1 <- max(vero.mat.phi.22.sigma.1) 
aprox.mat.m.phi.22.sigma.1 <- max(aprox.mat.phi.22.sigma.1)

m.phi.22.sigma.2 <- max(vero.mat.phi.22.sigma.2) 
aprox.mat.m.phi.22.sigma.2 <- max(aprox.mat.phi.22.sigma.2)

m.phi.22.rho.12 <- max(vero.mat.phi.22.rho.12) 
aprox.mat.m.phi.22.rho.12 <- max(aprox.mat.phi.22.rho.12)

m.sigma.1.sigma.2 <- max(vero.mat.sigma.1.sigma.2) 
aprox.mat.m.sigma.1.sigma.2 <- max(aprox.mat.sigma.1.sigma.2)

m.sigma.1.sigma.2 <- max(vero.mat.sigma.1.sigma.2) 
aprox.mat.m.sigma.1.sigma.2 <- max(aprox.mat.sigma.1.sigma.2)

m.sigma.1.rho.12 <- max(vero.mat.sigma.1.rho.12) 
aprox.mat.m.sigma.1.rho.12 <- max(aprox.mat.sigma.1.rho.12)

m.sigma.2.rho.12 <- max(vero.mat.sigma.2.rho.12) 
aprox.mat.m.sigma.2.rho.12 <- max(aprox.mat.sigma.2.rho.12)

#### ARQUIVOS ####

setwd("C:/Users/Caroline/Dropbox/Projeto_AproxSep/Scripts/Caroline/IC_TCC/S4_contorno/Figuras")

pdf(paste0("curvas_contorno_agrupado_20220910_v3_cenario_",cenario,"_n_",n,".pdf",sep=""),
    width=18,
    height=15)

jpeg(file = paste0("curvas_contorno_agrupado_20220910_v3_cenario_",cenario,"_n_",n,".jpeg",sep=""), 
     quality=100, 
     width=550, 
     height=500)

#### GRAFICOS AGRUPADOS - cenario 1 ####

par(mfrow=c(3,3),mar=c(5,5,0.8,0.1),cex=1)# fazer 1.2 para jpeg e 2.8 pdf
seq.phi.11 <- seq(0.02,0.5,l=dim(exp(vero.mat.phi.11.sigma.1-m.phi.11.sigma.1))[1])
seq.s.1 <- seq(.55,2,l=dim(exp(vero.mat.phi.11.sigma.1-m.phi.11.sigma.1))[1])
contour(seq.phi.11, 
        seq.s.1, 
        exp(vero.mat.phi.11.sigma.1-m.phi.11.sigma.1),
        xlim=c(0.01,0.28),
        ylim=c(0.58,1.3),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(sigma[1]),
        xlab = expression(phi[11]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.phi.11,
        seq.s.1,
        exp(aprox.mat.phi.11.sigma.1-aprox.mat.m.phi.11.sigma.1),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=sigma.vet[1],lty = 2, lwd = 2, col="gray40")
abline(v=phi.11,lty = 2, lwd = 2, col="gray40")


# seq.phi.11 <- seq(0.01,2,l=dim(exp(vero.mat.phi.11.sigma.2-m.phi.11.sigma.2))[1])
# seq.s.2 <- seq(1.4,2.1,l=dim(exp(vero.mat.phi.11.sigma.1-m.phi.11.sigma.1))[1])

contour(seq.phi.11, 
        seq.s.2, 
        exp(vero.mat.phi.11.sigma.2-m.phi.11.sigma.2),
        xlim=c(0.05,.4),
        ylim=c(.2,4.5),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(sigma[2]),
        xlab = expression(phi[11]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.phi.11,
        seq.s.2,
        exp(aprox.mat.phi.11.sigma.2-aprox.mat.m.phi.11.sigma.2),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=sigma.2,lty = 2, lwd = 3, col="gray40")
abline(v=phi.11,lty = 2, lwd = 3, col="gray40")

# seq.phi.11 <- seq(0.02, 0.6, l = dim(exp(vero.mat.phi.11.rho.12-m.phi.11.rho.12))[1])
# seq.rho.12 <- seq(0.19,1,l= dim(exp(vero.mat.phi.11.rho.12-m.phi.11.rho.12))[1])

contour(seq.phi.11, 
        seq.rho.12, 
        exp(vero.mat.phi.11.rho.12-m.phi.11.rho.12),
        xlim=c(0.01,0.37),
        # ylim=c(0.25,0.65),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(rho[12]),
        xlab = expression(phi[11]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.phi.11,
        seq.rho.12,
        exp(aprox.mat.phi.11.rho.12-aprox.mat.m.phi.11.rho.12),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=rho.12,lty = 2, lwd = 3, col="gray40")
abline(v=phi.11,lty = 2, lwd = 3, col="gray40")

# 
# seq.phi.22 <- seq(0.01,.03,l= dim(exp(vero.mat.phi.22.sigma.1-m.phi.22.sigma.1))[1])
# seq.s.1 <- seq(.07,1.95,l=dim(exp(vero.mat.phi.22.sigma.1-m.phi.22.sigma.1))[1])
contour(seq.phi.22, 
        seq.s.1, 
        exp(vero.mat.phi.22.sigma.1-m.phi.22.sigma.1),
        xlim=c(0.05,.35),
        ylim=c(.7,2),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(sigma[1]),
        xlab = expression(phi[22]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.phi.22,
        seq.s.1,
        exp(aprox.mat.phi.22.sigma.1-aprox.mat.m.phi.22.sigma.1),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=sigma.1,lty = 2, lwd = 3, col="gray40")
abline(v=phi.22,lty = 2, lwd = 3, col="gray40")


# seq.phi.22 <- seq(0.02, 0.55, l = dim(exp(vero.mat.phi.22.sigma.2-m.phi.22.sigma.2))[1])
# seq.s.2 <- seq(0.3,2.5,l= dim(exp(vero.mat.phi.22.sigma.2-m.phi.22.sigma.2))[1])

contour(seq.phi.22, 
        seq.s.2, 
        exp(vero.mat.phi.22.sigma.2-m.phi.22.sigma.2),
        xlim=c(0.05,.4),
        ylim=c(.7,4.6),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(sigma[2]),
        xlab = expression(phi[22]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.phi.22,
        seq.s.2,
        exp(aprox.mat.phi.22.sigma.2-aprox.mat.m.phi.22.sigma.2),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=sigma.2,lty = 2, lwd = 3, col="gray40")
abline(v=phi.22,lty = 2, lwd = 3, col="gray40")


contour(seq.phi.22, 
        seq.rho.12, 
        exp(vero.mat.phi.22.rho.12-m.phi.22.rho.12),
        xlim=c(0.01,.2),
        ylim=c(0.05,.45),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(rho[12]),
        xlab = expression(phi[22]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.phi.22,
        seq.rho.12,
        exp(aprox.mat.phi.22.rho.12-aprox.mat.m.phi.22.rho.12),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=rho.12,lty = 2, lwd = 3, col="gray40")
abline(v=phi.22,lty = 2, lwd = 3, col="gray40")



contour(seq.s.1, 
        seq.s.2, 
        exp(vero.mat.sigma.1.sigma.2-m.sigma.1.sigma.2),
        xlim=c(0.7,1.7),
        ylim=c(.3,4),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(sigma[2]),
        xlab = expression(sigma[1]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.s.1,
        seq.s.2,
        exp(aprox.mat.sigma.1.sigma.2-aprox.mat.m.sigma.1.sigma.2),
        col='gray55',
        lwd=2,
        
        drawlabels = FALSE,
        add = TRUE)
abline(h=sigma.2,lty = 2, lwd = 3, col="gray40")
abline(v=sigma.1,lty = 2, lwd = 3, col="gray40")


contour(seq.s.1, 
        seq.rho.12, 
        exp(vero.mat.sigma.1.rho.12-m.sigma.1.rho.12),
        xlim=c(0.8,1.9),
        ylim=c(0.2,0.5),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(rho[12]),
        xlab = expression(sigma[1]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.s.1,
        seq.rho.12,
        exp(aprox.mat.sigma.1.rho.12-aprox.mat.m.sigma.1.rho.12),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=rho.12,lty = 2, lwd = 3, col="gray40")
abline(v=sigma.1,lty = 2, lwd = 3, col="gray40")



contour(seq.s.2, 
        seq.rho.12, 
        exp(vero.mat.sigma.2.rho.12-m.sigma.2.rho.12),
        xlim=c(1.4,4.2),
        ylim=c(0.2,.5),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(rho[12]),
        xlab = expression(sigma[2]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.s.2,
        seq.rho.12,
        exp(aprox.mat.sigma.2.rho.12-aprox.mat.m.sigma.2.rho.12),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=rho.12,lty = 2, lwd = 3, col="gray40")
abline(v=sigma.2,lty = 2, lwd = 3, col="gray40")

######
dev.off()
#### GRAFICOS AGRUPADOS - cenario 2 ####

par(mfrow=c(3,3),mar=c(5,5,0.15,0.15),cex=1)# fazer 1.2 para jpeg e 2.5 pdf
seq.phi.11 <- seq(0.02,0.5,l=dim(exp(vero.mat.phi.11.sigma.1-m.phi.11.sigma.1))[1])
seq.s.1 <- seq(.55,2,l=dim(exp(vero.mat.phi.11.sigma.1-m.phi.11.sigma.1))[1])
contour(seq.phi.11, 
        seq.s.1, 
        exp(vero.mat.phi.11.sigma.1-m.phi.11.sigma.1),
        xlim=c(0.03,0.4),
        ylim=c(0.8,1.75),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(sigma[1]),
        xlab = expression(phi[11]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.phi.11,
        seq.s.1,
        exp(aprox.mat.phi.11.sigma.1-aprox.mat.m.phi.11.sigma.1),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=sigma.vet[1],lty = 2, lwd = 2, col="gray40")
abline(v=phi.11,lty = 2, lwd = 2, col="gray40")


# seq.phi.11 <- seq(0.01,2,l=dim(exp(vero.mat.phi.11.sigma.2-m.phi.11.sigma.2))[1])
# seq.s.2 <- seq(1.4,2.1,l=dim(exp(vero.mat.phi.11.sigma.1-m.phi.11.sigma.1))[1])

contour(seq.phi.11, 
        seq.s.2, 
        exp(vero.mat.phi.11.sigma.2-m.phi.11.sigma.2),
        xlim=c(0.02,.48),
        ylim=c(.4,4),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(sigma[2]),
        xlab = expression(phi[11]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.phi.11,
        seq.s.2,
        exp(aprox.mat.phi.11.sigma.2-aprox.mat.m.phi.11.sigma.2),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=sigma.2,lty = 2, lwd = 3, col="gray40")
abline(v=phi.11,lty = 2, lwd = 3, col="gray40")

# seq.phi.11 <- seq(0.02, 0.6, l = dim(exp(vero.mat.phi.11.rho.12-m.phi.11.rho.12))[1])
# seq.rho.12 <- seq(0.19,1,l= dim(exp(vero.mat.phi.11.rho.12-m.phi.11.rho.12))[1])

contour(seq.phi.11, 
        seq.rho.12, 
        exp(vero.mat.phi.11.rho.12-m.phi.11.rho.12),
        xlim=c(0.05,0.47),
        # ylim=c(0.25,0.65),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(rho[12]),
        xlab = expression(phi[11]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.phi.11,
        seq.rho.12,
        exp(aprox.mat.phi.11.rho.12-aprox.mat.m.phi.11.rho.12),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=rho.12,lty = 2, lwd = 3, col="gray40")
abline(v=phi.11,lty = 2, lwd = 3, col="gray40")

# 
# seq.phi.22 <- seq(0.01,.03,l= dim(exp(vero.mat.phi.22.sigma.1-m.phi.22.sigma.1))[1])
# seq.s.1 <- seq(.07,1.95,l=dim(exp(vero.mat.phi.22.sigma.1-m.phi.22.sigma.1))[1])
contour(seq.phi.22, 
        seq.s.1, 
        exp(vero.mat.phi.22.sigma.1-m.phi.22.sigma.1),
        xlim=c(0.01,.4),
        ylim=c(.5,2),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(sigma[1]),
        xlab = expression(phi[22]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.phi.22,
        seq.s.1,
        exp(aprox.mat.phi.22.sigma.1-aprox.mat.m.phi.22.sigma.1),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=sigma.1,lty = 2, lwd = 3, col="gray40")
abline(v=phi.22,lty = 2, lwd = 3, col="gray40")


# seq.phi.22 <- seq(0.02, 0.55, l = dim(exp(vero.mat.phi.22.sigma.2-m.phi.22.sigma.2))[1])
# seq.s.2 <- seq(0.3,2.5,l= dim(exp(vero.mat.phi.22.sigma.2-m.phi.22.sigma.2))[1])

contour(seq.phi.22, 
        seq.s.2, 
        exp(vero.mat.phi.22.sigma.2-m.phi.22.sigma.2),
        xlim=c(0.01,.37),
        ylim=c(.2,4.5),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(sigma[2]),
        xlab = expression(phi[22]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.phi.22,
        seq.s.2,
        exp(aprox.mat.phi.22.sigma.2-aprox.mat.m.phi.22.sigma.2),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=sigma.2,lty = 2, lwd = 3, col="gray40")
abline(v=phi.22,lty = 2, lwd = 3, col="gray40")


contour(seq.phi.22, 
        seq.rho.12, 
        exp(vero.mat.phi.22.rho.12-m.phi.22.rho.12),
        xlim=c(0.03,.3),
        ylim=c(0.1,.55),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(rho[12]),
        xlab = expression(phi[22]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.phi.22,
        seq.rho.12,
        exp(aprox.mat.phi.22.rho.12-aprox.mat.m.phi.22.rho.12),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=rho.12,lty = 2, lwd = 3, col="gray40")
abline(v=phi.22,lty = 2, lwd = 3, col="gray40")



contour(seq.s.1, 
        seq.s.2, 
        exp(vero.mat.sigma.1.sigma.2-m.sigma.1.sigma.2),
        xlim=c(0.7,1.7),
        ylim=c(.95,4.5),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(sigma[2]),
        xlab = expression(sigma[1]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.s.1,
        seq.s.2,
        exp(aprox.mat.sigma.1.sigma.2-aprox.mat.m.sigma.1.sigma.2),
        col='gray55',
        lwd=2,
        
        drawlabels = FALSE,
        add = TRUE)
abline(h=sigma.2,lty = 2, lwd = 3, col="gray40")
abline(v=sigma.1,lty = 2, lwd = 3, col="gray40")


contour(seq.s.1, 
        seq.rho.12, 
        exp(vero.mat.sigma.1.rho.12-m.sigma.1.rho.12),
        xlim=c(0.55,1.85),
        ylim=c(0.09,0.45),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(rho[12]),
        xlab = expression(sigma[1]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.s.1,
        seq.rho.12,
        exp(aprox.mat.sigma.1.rho.12-aprox.mat.m.sigma.1.rho.12),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=rho.12,lty = 2, lwd = 3, col="gray40")
abline(v=sigma.1,lty = 2, lwd = 3, col="gray40")



contour(seq.s.2, 
        seq.rho.12, 
        exp(vero.mat.sigma.2.rho.12-m.sigma.2.rho.12),
        xlim=c(.55,3.8),
        ylim=c(0.08,.45),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(rho[12]),
        xlab = expression(sigma[2]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.s.2,
        seq.rho.12,
        exp(aprox.mat.sigma.2.rho.12-aprox.mat.m.sigma.2.rho.12),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=rho.12,lty = 2, lwd = 3, col="gray40")
abline(v=sigma.2,lty = 2, lwd = 3, col="gray40")

######
dev.off()
#### GRAFICOS AGRUPADOS - cenario 3 ####

par(mfrow=c(3,3),mar=c(5,5,0.15,0.15),cex=2.5)# fazer 1 para jpeg e 2.5 pdf
seq.phi.11 <- seq(0.02,0.5,l=dim(exp(vero.mat.phi.11.sigma.1-m.phi.11.sigma.1))[1])
seq.s.1 <- seq(.55,2,l=dim(exp(vero.mat.phi.11.sigma.1-m.phi.11.sigma.1))[1])
contour(seq.phi.11, 
        seq.s.1, 
        exp(vero.mat.phi.11.sigma.1-m.phi.11.sigma.1),
        xlim=c(0.03,0.5),
        ylim=c(0.6,2),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(sigma[1]),
        xlab = expression(phi[11]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.phi.11,
        seq.s.1,
        exp(aprox.mat.phi.11.sigma.1-aprox.mat.m.phi.11.sigma.1),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=sigma.vet[1],lty = 2, lwd = 2, col="gray40")
abline(v=phi.11,lty = 2, lwd = 2, col="gray40")


# seq.phi.11 <- seq(0.01,2,l=dim(exp(vero.mat.phi.11.sigma.2-m.phi.11.sigma.2))[1])
# seq.s.2 <- seq(1.4,2.1,l=dim(exp(vero.mat.phi.11.sigma.1-m.phi.11.sigma.1))[1])

contour(seq.phi.11, 
        seq.s.2, 
        exp(vero.mat.phi.11.sigma.2-m.phi.11.sigma.2),
        xlim=c(0.02,.47),
        ylim=c(.3,4.5),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(sigma[2]),
        xlab = expression(phi[11]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.phi.11,
        seq.s.2,
        exp(aprox.mat.phi.11.sigma.2-aprox.mat.m.phi.11.sigma.2),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=sigma.2,lty = 2, lwd = 3, col="gray40")
abline(v=phi.11,lty = 2, lwd = 3, col="gray40")

# seq.phi.11 <- seq(0.02, 0.6, l = dim(exp(vero.mat.phi.11.rho.12-m.phi.11.rho.12))[1])
# seq.rho.12 <- seq(0.19,1,l= dim(exp(vero.mat.phi.11.rho.12-m.phi.11.rho.12))[1])

contour(seq.phi.11, 
        seq.rho.12, 
        exp(vero.mat.phi.11.rho.12-m.phi.11.rho.12),
        xlim=c(0.01,0.4),
        # ylim=c(0.25,0.65),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(rho[12]),
        xlab = expression(phi[11]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.phi.11,
        seq.rho.12,
        exp(aprox.mat.phi.11.rho.12-aprox.mat.m.phi.11.rho.12),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=rho.12,lty = 2, lwd = 3, col="gray40")
abline(v=phi.11,lty = 2, lwd = 3, col="gray40")

# 
# seq.phi.22 <- seq(0.01,.03,l= dim(exp(vero.mat.phi.22.sigma.1-m.phi.22.sigma.1))[1])
# seq.s.1 <- seq(.07,1.95,l=dim(exp(vero.mat.phi.22.sigma.1-m.phi.22.sigma.1))[1])
contour(seq.phi.22, 
        seq.s.1, 
        exp(vero.mat.phi.22.sigma.1-m.phi.22.sigma.1),
        xlim=c(0.01,.55),
        ylim=c(.52,1.95),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(sigma[1]),
        xlab = expression(phi[22]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.phi.22,
        seq.s.1,
        exp(aprox.mat.phi.22.sigma.1-aprox.mat.m.phi.22.sigma.1),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=sigma.1,lty = 2, lwd = 3, col="gray40")
abline(v=phi.22,lty = 2, lwd = 3, col="gray40")


# seq.phi.22 <- seq(0.02, 0.55, l = dim(exp(vero.mat.phi.22.sigma.2-m.phi.22.sigma.2))[1])
# seq.s.2 <- seq(0.3,2.5,l= dim(exp(vero.mat.phi.22.sigma.2-m.phi.22.sigma.2))[1])

contour(seq.phi.22, 
        seq.s.2, 
        exp(vero.mat.phi.22.sigma.2-m.phi.22.sigma.2),
        xlim=c(0.01,.53),
        ylim=c(.2,4),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(sigma[2]),
        xlab = expression(phi[22]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.phi.22,
        seq.s.2,
        exp(aprox.mat.phi.22.sigma.2-aprox.mat.m.phi.22.sigma.2),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=sigma.2,lty = 2, lwd = 3, col="gray40")
abline(v=phi.22,lty = 2, lwd = 3, col="gray40")


contour(seq.phi.22, 
        seq.rho.12, 
        exp(vero.mat.phi.22.rho.12-m.phi.22.rho.12),
        xlim=c(0.01,.53),
        ylim=c(0.05,.52),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(rho[12]),
        xlab = expression(phi[22]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.phi.22,
        seq.rho.12,
        exp(aprox.mat.phi.22.rho.12-aprox.mat.m.phi.22.rho.12),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=rho.12,lty = 2, lwd = 3, col="gray40")
abline(v=phi.22,lty = 2, lwd = 3, col="gray40")



contour(seq.s.1, 
        seq.s.2, 
        exp(vero.mat.sigma.1.sigma.2-m.sigma.1.sigma.2),
        xlim=c(0.5,2.1),
        ylim=c(.05,4.5),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(sigma[2]),
        xlab = expression(sigma[1]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.s.1,
        seq.s.2,
        exp(aprox.mat.sigma.1.sigma.2-aprox.mat.m.sigma.1.sigma.2),
        col='gray55',
        lwd=2,
        
        drawlabels = FALSE,
        add = TRUE)
abline(h=sigma.2,lty = 2, lwd = 3, col="gray40")
abline(v=sigma.1,lty = 2, lwd = 3, col="gray40")


contour(seq.s.1, 
        seq.rho.12, 
        exp(vero.mat.sigma.1.rho.12-m.sigma.1.rho.12),
        xlim=c(0.55,2),
        ylim=c(0.05,0.53),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(rho[12]),
        xlab = expression(sigma[1]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.s.1,
        seq.rho.12,
        exp(aprox.mat.sigma.1.rho.12-aprox.mat.m.sigma.1.rho.12),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=rho.12,lty = 2, lwd = 3, col="gray40")
abline(v=sigma.1,lty = 2, lwd = 3, col="gray40")



contour(seq.s.2, 
        seq.rho.12, 
        exp(vero.mat.sigma.2.rho.12-m.sigma.2.rho.12),
        xlim=c(.05,4.5),
        ylim=c(0.05,.55),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(rho[12]),
        xlab = expression(sigma[2]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.s.2,
        seq.rho.12,
        exp(aprox.mat.sigma.2.rho.12-aprox.mat.m.sigma.2.rho.12),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=rho.12,lty = 2, lwd = 3, col="gray40")
abline(v=sigma.2,lty = 2, lwd = 3, col="gray40")

######
dev.off()
#### GRAFICOS AGRUPADOS - cenario 4 ####

par(mfrow=c(3,3),mar=c(5,5,0.15,0.15),cex=2.5)# fazer 1 para jpeg e 2.5 pdf
seq.phi.11 <- seq(0,0.5,l=dim(exp(vero.mat.phi.11.sigma.1-m.phi.11.sigma.1))[1])
seq.s.1 <- seq(.55,2,l=dim(exp(vero.mat.phi.11.sigma.1-m.phi.11.sigma.1))[1])
contour(seq.phi.11, 
        seq.s.1, 
        exp(vero.mat.phi.11.sigma.1-m.phi.11.sigma.1),
        xlim=c(0.02,0.48),
        ylim=c(0.5,2),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(sigma[1]),
        xlab = expression(phi[11]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.phi.11,
        seq.s.1,
        exp(aprox.mat.phi.11.sigma.1-aprox.mat.m.phi.11.sigma.1),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=sigma.vet[1],lty = 2, lwd = 2, col="gray40")
abline(v=phi.11,lty = 2, lwd = 2, col="gray40")


# seq.phi.11 <- seq(0.01,2,l=dim(exp(vero.mat.phi.11.sigma.2-m.phi.11.sigma.2))[1])
# seq.s.2 <- seq(1.4,2.1,l=dim(exp(vero.mat.phi.11.sigma.1-m.phi.11.sigma.1))[1])

contour(seq.phi.11, 
        seq.s.2, 
        exp(vero.mat.phi.11.sigma.2-m.phi.11.sigma.2),
        xlim=c(0.03,.45),
        ylim=c(.3,4.8),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(sigma[2]),
        xlab = expression(phi[11]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.phi.11,
        seq.s.2,
        exp(aprox.mat.phi.11.sigma.2-aprox.mat.m.phi.11.sigma.2),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=sigma.2,lty = 2, lwd = 3, col="gray40")
abline(v=phi.11,lty = 2, lwd = 3, col="gray40")

# seq.phi.11 <- seq(0.02, 0.6, l = dim(exp(vero.mat.phi.11.rho.12-m.phi.11.rho.12))[1])
# seq.rho.12 <- seq(0.19,1,l= dim(exp(vero.mat.phi.11.rho.12-m.phi.11.rho.12))[1])

contour(seq.phi.11, 
        seq.rho.12, 
        exp(vero.mat.phi.11.rho.12-m.phi.11.rho.12),
        xlim=c(0,.48),
        ylim=c(0.09,0.45),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(rho[12]),
        xlab = expression(phi[11]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.phi.11,
        seq.rho.12,
        exp(aprox.mat.phi.11.rho.12-aprox.mat.m.phi.11.rho.12),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=rho.12,lty = 2, lwd = 3, col="gray40")
abline(v=phi.11,lty = 2, lwd = 3, col="gray40")

# 
# seq.phi.22 <- seq(0.01,.03,l= dim(exp(vero.mat.phi.22.sigma.1-m.phi.22.sigma.1))[1])
# seq.s.1 <- seq(.07,1.95,l=dim(exp(vero.mat.phi.22.sigma.1-m.phi.22.sigma.1))[1])
contour(seq.phi.22, 
        seq.s.1, 
        exp(vero.mat.phi.22.sigma.1-m.phi.22.sigma.1),
        xlim=c(0.01,.73),
        ylim=c(.52,2),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(sigma[1]),
        xlab = expression(phi[22]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.phi.22,
        seq.s.1,
        exp(aprox.mat.phi.22.sigma.1-aprox.mat.m.phi.22.sigma.1),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=sigma.1,lty = 2, lwd = 3, col="gray40")
abline(v=phi.22,lty = 2, lwd = 3, col="gray40")


# seq.phi.22 <- seq(0.02, 0.55, l = dim(exp(vero.mat.phi.22.sigma.2-m.phi.22.sigma.2))[1])
# seq.s.2 <- seq(0.3,2.5,l= dim(exp(vero.mat.phi.22.sigma.2-m.phi.22.sigma.2))[1])

contour(seq.phi.22, 
        seq.s.2, 
        exp(vero.mat.phi.22.sigma.2-m.phi.22.sigma.2),
        xlim=c(0.01,.73),
        ylim=c(.09,4.8),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(sigma[2]),
        xlab = expression(phi[22]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.phi.22,
        seq.s.2,
        exp(aprox.mat.phi.22.sigma.2-aprox.mat.m.phi.22.sigma.2),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=sigma.2,lty = 2, lwd = 3, col="gray40")
abline(v=phi.22,lty = 2, lwd = 3, col="gray40")


contour(seq.phi.22, 
        seq.rho.12, 
        exp(vero.mat.phi.22.rho.12-m.phi.22.rho.12),
        xlim=c(0.01,.73),
        ylim=c(0.05,.52),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(rho[12]),
        xlab = expression(phi[22]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.phi.22,
        seq.rho.12,
        exp(aprox.mat.phi.22.rho.12-aprox.mat.m.phi.22.rho.12),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=rho.12,lty = 2, lwd = 3, col="gray40")
abline(v=phi.22,lty = 2, lwd = 3, col="gray40")



contour(seq.s.1, 
        seq.s.2, 
        exp(vero.mat.sigma.1.sigma.2-m.sigma.1.sigma.2),
        xlim=c(0.57,1.85),
        ylim=c(.05,5),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(sigma[2]),
        xlab = expression(sigma[1]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.s.1,
        seq.s.2,
        exp(aprox.mat.sigma.1.sigma.2-aprox.mat.m.sigma.1.sigma.2),
        col='gray55',
        lwd=2,
        
        drawlabels = FALSE,
        add = TRUE)
abline(h=sigma.2,lty = 2, lwd = 3, col="gray40")
abline(v=sigma.1,lty = 2, lwd = 3, col="gray40")


contour(seq.s.1, 
        seq.rho.12, 
        exp(vero.mat.sigma.1.rho.12-m.sigma.1.rho.12),
        xlim=c(0.65,1.8),
        ylim=c(0.05,0.45),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(rho[12]),
        xlab = expression(sigma[1]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.s.1,
        seq.rho.12,
        exp(aprox.mat.sigma.1.rho.12-aprox.mat.m.sigma.1.rho.12),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=rho.12,lty = 2, lwd = 3, col="gray40")
abline(v=sigma.1,lty = 2, lwd = 3, col="gray40")



contour(seq.s.2, 
        seq.rho.12, 
        exp(vero.mat.sigma.2.rho.12-m.sigma.2.rho.12),
        xlim=c(.05,4.6),
        ylim=c(0.05,.5),
        lwd=2,
        drawlabels = FALSE, 
        ylab = expression(rho[12]),
        xlab = expression(sigma[2]),
        cex.lab =1.5,
        cex.axis = 1.5)
contour(seq.s.2,
        seq.rho.12,
        exp(aprox.mat.sigma.2.rho.12-aprox.mat.m.sigma.2.rho.12),
        col='gray55',
        lwd=2,
        drawlabels = FALSE,
        add = TRUE)
abline(h=rho.12,lty = 2, lwd = 3, col="gray40")
abline(v=sigma.2,lty = 2, lwd = 3, col="gray40")

######
dev.off()



