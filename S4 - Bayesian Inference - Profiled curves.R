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


matriz.cheia <- function(sigma,rho,H_perfilada,phi){
    p <- dim(sigma)[1]
    n <- dim(H_perfilada)[1]
    MAT.COV_perfilada <- matrix(NA, ncol = n*p, nrow=n*p)
    sapply(1:p, 
         function(i){
           sapply(i:p, 
                  function(j){
                    # Basta gerar sigma pelo produto de sigma_i vezes sigma_j
                    # Constroi os 4 blocos
                    MAT.COV_perfilada[(((i-1)*n)+1):(i*n),(((j-1)*n)+1):(j*n)] <<- sigma[i,j]*rho[i,j]*((1+(H_perfilada/phi[i,j]))^(-1))
                    # Repete p banda inferior
                    MAT.COV_perfilada[lower.tri(MAT.COV_perfilada)] <<- t(MAT.COV_perfilada)[lower.tri(MAT.COV_perfilada)]
                  })
         })
    return(MAT.COV_perfilada)}

sub.matriz <- function(M,r,c){
  rg <- (row(M)-1)%/%r+1
  cg<- (col(M)-1)%/%c+1
  rci <- (rg-1)*max(cg) + cg
  N <- prod(dim(M))/r/c
  cv <- unlist(lapply(1:N,function(x) M[rci==x]))
  dim(cv) <- c(r,c,N)
  cv
} 

#help(backsolve)
#aprox.phi.11
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
   L <-1
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
  sapply(1:L, 
         function(dd){
           VERO.APROX.CH.DMV[dd] <<- dmvnorm.chol(y_vet[,dd],mu_vetor,-A,-R,log=TRUE)
         })  
  return(VERO.APROX.CH.DMV)
}

ler_arquivo_csv = function(nome, cenario_var, nseq_usado = 500){
  if (cenario_var == 1){
    nome_arquivo = paste0(nome,"_cenario_1nseq_",nseq_usado,".csv")
    return(as.matrix(read.csv(nome_arquivo, sep = ";", dec = ",")[,-1]))
  } else{
    print(NULL)
  }
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

#### Geração do Dado ####
# Quero comparar o tempo de todas as execucoes da verossimilhanca cheia e aproximada
# Para todas as configuracoes (combinacoes) de n e p
# qq = tamanho do vetor p
# qq=1
# ww = tamanho do vetor n
# ww=1

cenario = 0
# Cenario 
setwd(paste0("C:/Users/Caroline/Dropbox/Projeto_AproxSep/Scripts/Caroline/IC_TCC/S4_perfiladas/Dados_gerados"))
escolha_cenario(1)
# Numero de variaveis
p <- 2
# Numero de localizacoes
n <- 100
# Numero de amostras que quero --replicas
# L <- 1

# spatial range - alcance espacial

# lat <- runif(n,0,1)
# long <- runif(n,0,1)	
#write.csv2(cbind(lat,long),file =  paste0("Localizacoes_contorno_20220910_cenario_",cenario,"n_",n,".csv"))

# ou usar os salvos

lat_long = read.csv2(paste0("Localizacoes_perfilada_cenario_2n_100.csv"))

lat_perfilada = lat_long[,2]
long_perfilada = lat_long[,3]

plot(long_perfilada,lat_perfilada,pch=16)
### MATRIZ H, PHI, RHO, SIGMA ####

## Criando a Matriz de Distâncias H
# ij que varia de acordo com a componente
# h varia: matriz das distancias do espaço
H_perfilada<- matrix(0,nrow=n,ncol=n)

sapply(1:(n-1), 
       function(k){
         sapply((k+1):n, 
                function(l){
                  H_perfilada[k,l] <<- distancia(long_perfilada[k],lat_perfilada[k],long_perfilada[l],lat_perfilada[l])
                  H_perfilada[l,k] <<- H_perfilada[k,l]
                })
       })
dim(H_perfilada)

# Media usada foi 0: latitude e longitude nao estao entrando como covariaveis
#MU.vetor <- rep(0,(n*p)) # anterior

# Media usada foi 0: latitude e longitude nao estao entrando como covariaveis
MU.vetor <- rep(0,(n*p))

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

## MATRIZ DE COVARIÂNCIA

#MAT.COV_perfilada <- matrix(0,nrow=(n*p),ncol=(n*p)) 

MAT.COV_perfilada <- matrix(0,nrow=(n*p),ncol=(n*p)) 

# Como é simétrica, construo só a banda superior e a repito na banda inferior 
# sigma_i*sigma_j*rho_ij*exp(-h/phi_ij) -- funcao de correlacao espacial = exp
# para uma variavel c ela msm: vai dar 1 tds os elementos do produto
# ww = tamanho do vetor n
#ww=1

sapply(1:p, 
       function(i){
         sapply(i:p, 
                function(j){
                  # Basta gerar sigma pelo produto de sigma_i vezes sigma_j
                  # Constroi os 4 blocos
                  MAT.COV_perfilada[(((i-1)*n)+1):(i*n),(((j-1)*n)+1):(j*n)] <<- sigma[i,j]*rho[i,j]*((1+(H_perfilada/phi[i,j]))^(-1))
                  # Repete p banda inferior
                  MAT.COV_perfilada[lower.tri(MAT.COV_perfilada)] <<- t(MAT.COV_perfilada)[lower.tri(MAT.COV_perfilada)]
                })
       })

# image(MAT.COV_perfilada)
is.positive.definite(MAT.COV_perfilada)
#isSymmetric(MAT.COV_perfilada)


# São negativas definidas pois todos os autovalores sao negativos
# eigen(A)
# Error in eigen(R) : matriz não quadrada em 'eigen'
#eigen(R)$values 

# Geramos o Y verdadeiro utilizando matriz cheia, nao aproximada
# erro <- rmvnorm(1,rep(0,n*p),MAT.COV_perfilada)
# Y.vetor_perfilada <- MU.vetor + t(erro)

# sempre salvar o dado
#write.csv2(Y.vetor_perfilada,file = paste0("Y_contorno_20220910_cenario_",cenario,"n_",n,".csv"))

# ou usar o salvo

Y.vetor_perfilada <- read.csv2(file = 'Y_perfilada_cenario_v22n_100.csv')[,2]#paste("Y_cenario_",cenario,"n_",n,".csv")
hist(Y.vetor_perfilada[1:100])
hist(Y.vetor_perfilada[101:200])

#### Análise da Verossimilhança ####
L <- 1 #(replicas)

# Fixar todos os params e variar 1 de cada vez (inicialmente)
n.seq <-  10000# fazer 250 para os demais, 100 é para contorno

#### ANOTACOES ####
# Vimos que reduz o tempo conforme aumenta a dimensão (proporcionalmente)
# Mas será a qualidade boa? É o que avaliaremos.
# Desenho da vero: curva de contorno. De 2 em 2 parametros

# Gerar lat e long
# Definir valores verdadeiros dos params
# Calcular matrix de cov
# Gerar o y
# Calcular a vero
# Faremos sem replicas

# rho12, phi11, phi12, phi22, sigma1, sigma2
# fazer para a cheia e a aproximada
# dentro do msm for consigo fazer para cada item, calcular a matriz cheia e a aproximada
# e a vero cheia e a vero aproximada (dmvnorm, dmvnorm.chol)

# Avaliamos nos intervalos estrategicos: ja temos os valores verdadeiros
# Queremos ver que de fato estao centradas no vdd valor
# Queremos que a aproximada e a cheia sejam centradas no real

# Serao 6 graficos, colocar em uma figura só -- individuais

# Na hora dos parametros 2 a 2 será um for duplo, apenas. O resto permanece igual

# Primeiro: fazer o do rho
# mudando de 0.5 para 0.6 (cauchy)
# restricao funcao de correlacao cauchy
# rho.12 <= phi.11*phi.22/(phi.12^2) = 0.75 no cenario 3

#### SEQUENCIAS ####

#-0.5
seq.rho.12 <- seq(0.1,(phi.11*phi.22/(phi.12^2)),l=n.seq)

seq.phi.11 <- seq(0.02,0.4,l=n.seq)
seq.phi.22 <- seq(0.02,0.4,l=n.seq)
#seq.phi.12 <- seq(0.05,0.5,l=n.seq)

seq.s.1 <- seq(.05,1.8,l=n.seq)
seq.s.2 <- seq(.05,1.8,l=n.seq) #sigmas

vero.phi.11 <- vero.phi.22 <- vero.phi.12 <- rep(0,n.seq); 
vero.rho.12 <- rep(0,n.seq)
vero.s.1 <- vero.s.2 <- rep(0,n.seq)

#aproximacoes
aprox.phi.11 <- aprox.phi.22 <- aprox.phi.12 <- rep(0,n.seq); 
aprox.s.1 <- aprox.s.2 <- rep(0,n.seq)
aprox.rho.12 <- rep(0,n.seq);

#### ANALISE PERFILADA ####

#### PHI.11 ####
i = 50
seq.phi.11 <- seq(0.02,.22,l=n.seq)
for(i in 1:n.seq){
  phi.p <- matrix(rep(0,(p^2)),ncol=p) # phi do intervalo
  # ver restricoes de phi.12, rho.12
  phi.p[lower.tri(phi.p,diag=TRUE)] <- c(seq.phi.11[i],(seq.phi.11[i] + phi.22)/2,phi.22)
  phi.p[upper.tri(phi.p)] <- t(phi.p)[upper.tri(phi.p)]
  cov <- matriz.cheia(sigma,rho,H_perfilada,phi.p) #funcao que preenche MAT.COV_perfilada
  vero.phi.11[i] <- dmvnorm(as.vector(Y.vetor_perfilada),MU.vetor,cov,log = TRUE)#;print(i)
  aprox.phi.11[i] <- vero.aprox(cov,n,as.matrix(Y.vetor_perfilada),MU.vetor);print(i)
  } 
# fazer log-verossimilhanca
(m.phi.11 <- max(vero.phi.11))# transformacao para mudar a escala da densidade de 0 a 1
# corrige tbm o problema da funcao explodir, que daria infinito no R

# APROX
(aprox.m.phi.11 <- max(aprox.phi.11))

par(mar=c(4,4,0.5,0.5),cex=1)
# padronizacao, densidade no max = 1, sem explodir e ter problemas de escala 
plot(seq.phi.11,
     exp(vero.phi.11-m.phi.11),
     type='l',
     xlab=expression(phi[11]),
     ylab='density',
     #xlim=c(0.00001,.65),
     ylim=c(0,1),
     cex.lab=2, cex.axis=1.6)
lines(seq.phi.11,exp(aprox.phi.11-aprox.m.phi.11),type='l', col='gray55',
      lwd=2)
abline(v=phi.11,lty=2,lwd=2)
abline(h=0)
dev.off()

phi_11_perfilada = cbind(vero.phi.11, aprox.phi.11)
write.csv2(phi_11_perfilada,file =paste0("phi.11_perfiladas_cenario_",cenario,"_nseq",n.seq,".csv"))
getwd()
#### PHI.22 ####
seq.phi.22 <- seq(0.05,1.05,l=n.seq)
for(i in 1:n.seq){
  phi.p <- matrix(rep(0,(p^2)),ncol=p) # phi do intervalo
  # ver restricoes de phi.12, rho.12
  phi.p[lower.tri(phi.p,diag=TRUE)] <- c(phi.11,(phi.11 + seq.phi.22[i])/2,seq.phi.22[i])
  phi.p[upper.tri(phi.p)] <- t(phi.p)[upper.tri(phi.p)]
  cov <- matriz.cheia(sigma,rho,H_perfilada,phi.p) #funcao que preenche MAT.COV_perfilada
  vero.phi.22[i] <- dmvnorm(as.vector(Y.vetor_perfilada),MU.vetor,cov,log = TRUE)
  aprox.phi.22[i] <- vero.aprox(cov,n,Y.vetor_perfilada,MU.vetor);print(i)
} 
m.phi.22 <- max(vero.phi.22)
aprox.m.phi.22 <- max(aprox.phi.22)
par(mar=c(4,4,0.5,0.5),cex=1)
plot(seq.phi.22,
     exp(vero.phi.22-m.phi.22),
     type='l',
     #xlim=c(0.1,0.8),
     xlab=expression(phi[22]),
     ylab='density',
     cex.lab=1.8, 
     cex.axis=1.8)
lines(seq.phi.22,exp(aprox.phi.22-aprox.m.phi.22),type='l', col='gray55',
      lwd=2)
abline(v=phi.22,lty=2,lwd=2)
abline(h=0)


#### SIGMA 1 ####
seq.s.1 <- seq(.75,2,l=n.seq)
i = 1
for(i in 1:n.seq){
  s.vet <- c(seq.s.1[i]^2,seq.s.1[i]*sigma.2,sigma.2^2)
  sigma.p <- matrix(rep(NA,(p^2)),ncol=p)
  sigma.p[lower.tri(sigma.p, diag=TRUE)] <- s.vet
  sigma.p[upper.tri(sigma.p)] <- t(sigma.p)[upper.tri(sigma.p)]
  cov <- matriz.cheia(sigma.p,rho,H_perfilada,phi)
  vero.s.1[i] <- dmvnorm(as.vector(Y.vetor_perfilada),MU.vetor,cov,log = TRUE)
  aprox.s.1[i] <- vero.aprox(cov,n,Y.vetor_perfilada, MU.vetor);print(i)
  
}
(m.vero.s.1 <- max(vero.s.1))
(aprox.m.s.1 <- max(aprox.s.1))


par(mar=c(4,4,0.5,0.5),cex=1)
plot(seq.s.1,
     exp(vero.s.1-m.vero.s.1),
     type='l',
     xlab=expression(sigma[1]),
     ylab='density',
     cex.lab=1.8, cex.axis=1.8)
lines(seq.s.1,exp(aprox.s.1-aprox.m.s.1),type='l', col='gray55',
      lwd=2)
abline(v=sigma.vet[1],lty=2,lwd=2)
abline(h=0)


#### SIGMA 2 ####
seq.s.2 <- seq(0.45,1.9,l=n.seq)
for(i in 1:n.seq){
  s.vet <- c(sigma.1^2,sigma.1*seq.s.2[i],seq.s.2[i]^2)
  sigma.p <- matrix(rep(NA,(p^2)),ncol=p)
  sigma.p[lower.tri(sigma.p, diag=TRUE)] <- s.vet
  sigma.p[upper.tri(sigma.p)] <- t(sigma.p)[upper.tri(sigma.p)]
  cov <- matriz.cheia(sigma.p,rho,H_perfilada,phi)
  vero.s.2[i] <- dmvnorm(as.vector(Y.vetor_perfilada),MU.vetor,cov,log = TRUE)
  
  aprox.s.2[i] <- vero.aprox(cov,n,Y.vetor_perfilada, MU.vetor);print(i)
  
}
(m.vero.s.2 <- max(vero.s.2))
(aprox.m.s.2 <- max(aprox.s.2))
par(mar=c(4,4,0.5,0.5),cex=1)
plot(seq.s.2,
     exp(vero.s.2-m.vero.s.2),
     type='l',
     xlab=expression(sigma[2]),
     ylab='density',
     cex.lab=1.8,
     cex.axis=1.8)
lines(seq.s.2,exp(aprox.s.2-aprox.m.s.2),type='l', col='gray55',
      lwd=2)
abline(v=sigma.vet[2],lty=2,lwd=2)
abline(h=0)

#### RHO 12 ####
seq.rho.12 <- seq(-0.01,(phi.11*phi.22/(phi.12^2)+.2),l=n.seq)
for(i in 1:n.seq){
  rho.vet <- seq.rho.12[i]
  rho.p <- matrix(rep(1,(p^2)),ncol=p)
  rho.p[lower.tri(rho.p)] <- rho.vet
  rho.p[upper.tri(rho.p)] <- t(rho.p)[upper.tri(rho.p)]
  cov <- matriz.cheia(sigma,rho.p,H_perfilada,phi)
  vero.rho.12[i] <- dmvnorm(as.vector(Y.vetor_perfilada),MU.vetor,cov,log = TRUE)
  aprox.rho.12[i] <- vero.aprox(cov,n,Y.vetor_perfilada, MU.vetor);print(i)
}
m.vero.rho.12 <- max(vero.rho.12)
aprox.m.rho.12 <- max(aprox.rho.12)
par(mar=c(4,4,0.5,0.5),cex=1)
plot(seq.rho.12,
     exp(vero.rho.12-m.vero.rho.12),
     type='l',
     xlab=expression(rho[12]),
     #xlim=c(-0.1,(phi.11*phi.22/(phi.12^2)-.3)),
     ylab='density',
     cex.lab=1.8, 
     cex.axis=1.8)
lines(seq.rho.12,exp(aprox.rho.12-aprox.m.rho.12),type='l', col='gray55',
      lwd=2)
abline(v=rho[1,2],lty=2,lwd=2)
abline(h=0)

dev.off()
####### GRAFICOS PERFILADA ##########

salvar_grafico = function(idioma, formato){
  setwd(paste0("C:/Users/Caroline/Dropbox/Projeto_AproxSep/Scripts/Caroline/IC_TCC/S4_perfiladas/Figuras/", idioma))
  
  if(idioma == 'pt'){
    densidade = 'densidade'
  }else{
    densidade = 'density'
  }
  #densidade'
  
  if(formato == 'pdf'){
    pdf(paste0("curvas_perfiladas_cenario_",cenario,"_n_",n,"_" ,idioma, ".pdf",sep=""),
        width=30,
        height=15)
    tamanho = 2.8
  }else{
    jpeg(file = paste0("curvas_perfiladas_cenario", cenario,
                       "_n_",n,"_" ,idioma,".jpeg",sep=""), 
                quality=100, width=1024, height=500) 
    tamanho = 1.4
  }

  
  #### GRAFICOS AGRUPADOS ####
  
  par(mfrow=c(2,3),mar=c(5,5.5,0.5,0.6),cex=tamanho)# fazer 1.4 para jpeg 2.8 pdf
  
  plot(seq.phi.11,
       exp(vero.phi.11-m.phi.11),
       cex.axis=1.7,
       cex.lab=1.7,
       type='l',
       lwd=6,
       xlab=expression(phi[11]),
       ylab=densidade)
  
  lines(seq.phi.11,
        exp(aprox.phi.11-aprox.m.phi.11),
        type='l', 
        col='gray55',
        lwd=2)
  abline(v=phi.11,lty=2,lwd=2)
  abline(h=0)
  
  plot(seq.phi.22,
       exp(vero.phi.22-m.phi.22),
       cex.axis=1.7,
       cex.lab=1.7,
       type='l',
       lwd=6,
       xlab=expression(phi[22]),
       ylab=densidade)
  
  lines(seq.phi.22,
        exp(aprox.phi.22-aprox.m.phi.22),
        type='l',
        col='gray55',
        lwd=2)
  abline(v=phi.22,lty=2,lwd=2)
  abline(h=0)
  
  
  plot(seq.s.1,
       exp(vero.s.1-m.vero.s.1),
       cex.axis=1.7,
       cex.lab=1.7,
       type='l',
       lwd=6,
       xim = c(0.7, 1.2),
       xlab=expression(sigma[11]),
       ylab=densidade)
  lines(seq.s.1,
        exp(aprox.s.1-aprox.m.s.1),
        type='l',
        col='gray55',
        lwd=2)
  abline(v=sigma.vet[1],lty=2,lwd=2)
  abline(h=0)
  
  plot(seq.s.2,
       exp(vero.s.2-m.vero.s.2),
       cex.axis=1.7,
       cex.lab=1.7,
       type='l',
       xim = c(1.3, 2.8),
       lwd=6,
       xlab=expression(sigma[22]),
       ylab=densidade)
  
  lines(seq.s.2,
        exp(aprox.s.2-aprox.m.s.2),
        type='l',
        col='gray55',
        lwd=2)
  abline(v=sigma.vet[2],lty=2,lwd=2)
  abline(h=0)
  
  plot(seq.rho.12,
       exp(vero.rho.12-m.vero.rho.12),
       cex.axis=1.7,
       cex.lab=1.7,
       type='l',
       #xlim=c(-0.1,(phi.11*phi.22/(phi.12^2)-.4)),
       ylim=c(0,1),
       lwd=6,
       xlab=expression(rho[12]),
       ylab=densidade)
  
  lines(seq.rho.12,
        exp(aprox.rho.12-aprox.m.rho.12),
        type='l',
        col='gray55',
        lwd=2)
  abline(v=rho[1,2],lty=2,lwd=2)
  abline(h=0)
  
  #######  ##########
  dev.off()
}


salvar_grafico('pt', 'jpeg')
# devido as restricoes da funcao de correlacao cauchy, 
# nao faremos as curvas perfiladas para phi.12
