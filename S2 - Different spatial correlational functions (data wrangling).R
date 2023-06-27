##-------------------------------------- ##
## ----- FUNÇÃO DMVNORM MODIFICADA ----- ##
##-------------------------------------- ##

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
require(ggplot2)#function keep

#############################################################


rm(list=ls())

# Dividindo a Matriz Cheia em submatrizes
## MATRIZ DE COVARIÂNCIA PERMUTADA
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

distancia <- function(x1,y1,x2,y2) {
  dist = sqrt((x1-x2)^2+(y1-y2)^2)
  return(dist)
}

save(sub.matriz, dmvnorm.chol, distancia, 
     file="VeroEnvFunctions_teste.RData")


ww = 1
ij = 1
#incremento = -0.01
#standardize.error.index = 0

####################################################
# REGERANDO USANDO A FUNÇÃO DE CORRELCAO EXPONENCIAL
####################################################

setwd("C:/Users/Caroline/Dropbox/Projeto_AproxSep/Scripts/Caroline/IC_TCC/S2")

for(ww in 1:1){
  for(ij in 1:4){
  
  # Cleans environment
  rm(list=ls()[!(ls() %in% c("ww","ij" 
                            # ,"incremento", "standardize.error.index", "H"
                             ))])
    
    
    corr.11 = list()
    corr.12 = list()
    corr.22 = list()
  
  # Loads necessary functions
  load("RData/VeroEnvFunctions_teste.RData")
  
  # incremento para variar cenarios
  # phi
  if(ij == 1){
    # cenario 1 - modelo separavel
    cenario <- 1
    phi.11 <- 0.1
    phi.22 <- 0.1
  }else if(ij == 2){
    # cenario 2 - modelo pouco nao separavel
    cenario <- 2
    phi.11 <- 0.1
    phi.22 <- 0.2
  }else if(ij == 3){
    # cenario 3 - modelo muito nao separavel
    cenario <- 3
    phi.11 <- 0.1
    phi.22 <- 0.5
  } else if(ij == 4){
    # cenario 3 - modelo muito nao separavel
    cenario <- 4
    phi.11 <- 0.1
    phi.22 <- 0.7
  }  # Para todas as configuracoes (combinacoes) de n e p
  # qq = tamanho do vetor p
  qq=1
  # ww = tamanho do vetor n
  ww=1
  #set.seed(117054003)
  # Numero de variaveis
  p <- 2
  # Numero de localizacoes
  n <- c(100)
  # Numero de amostras que quero --replicas
  L <- 1000 # repetir 1000 vezes
  
  
  lat <- runif(n[ww],0,1)
  long <- runif(n[ww],0,1)		
  
  plot(long,lat,pch=16)
  
  ## Criando a Matriz de Distâncias H
  # ij que varia de acordo com a componente
  # h varia: matriz das distancias do espaço
  H <- matrix(0,nrow=n[ww],ncol=n[ww])
  
  sapply(1:(n[ww]-1), 
         function(k){
           sapply((k+1):n[ww], 
                  function(l){
                    H[k,l] <<- distancia(long[k],lat[k],long[l],lat[l])
                    H[l,k] <<- H[k,l]
                  })
         })
  dim(H)
  
  # spatial range - alcance espacial
  
  # Media usada foi 0: latitude e longitude nao estao entrando como covariaveis
  MU.vetor <- rep(0,(n[ww]*p))
  
  phi.12 <- (phi.11+phi.22)/2
  # phi.12 <- sqrt((phi.11^2+phi.22^2)/2)
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
  
  
  # VETOR DE DISTANCIAS
  
  d = seq(0,2, l = L)
  
  
  # CORRELACAO
  
  # DUVIDAS
  
  # A corr deve ser uma >lista< com a funcao de correlacao
  # para cada cenario?
  
  # por que o vetor de distancias teria um indice? o objetivo
  # eh fazer todas as distancias de 0 a 2, certo?
  
  corr.11[[ij]]= exp(-d/phi.11)
  
  corr.12[[ij]] = exp(-d/phi.12)
  
  corr.22[[ij]] = exp(-d/phi.22)
  

  rho.12 <- 0.4 # fixar este valor, baseado na restricao do rho na funcao Cauchy
 
  #rho.12 <- sqrt((phi.11)*(phi.22))/(phi.12) # calculo feito da restricao do rho na funcao Cauchy
  
  # sqrt(0.5*((phi.11^2)+(phi.22^2)))
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
  
  MAT.COV <- matrix(0,nrow=(n[ww]*p),ncol=(n[ww]*p)) 
  
  # ---------- APOS AQUI ELA SERA POSITIVA DEFINIDA ---------------#
  # Como é simétrica, construo só a banda superior e a repito na banda inferior 
  # sigma_i*sigma_j*rho_ij*exp(-h/phi_ij) -- funcao de correlacao espacial = exp
  # para uma variavel c ela msm: vai dar 1 tds os elementos do produto
  sapply(1:p, 
         function(i){
           sapply(i:p, 
                  function(j){
                    # Basta gerar sigma pelo produto de sigma_i vezes sigma_j
                    # Constroi os 4 blocos
                    MAT.COV[(((i-1)*n[ww])+1):(i*n[ww]),(((j-1)*n[ww])+1):(j*n[ww])] <<- sigma[i,j]*rho[i,j]*exp(-H/phi[i,j])
                    # Repete p banda inferior
                    MAT.COV[lower.tri(MAT.COV)] <<- t(MAT.COV)[lower.tri(MAT.COV)]
                  })
         })
  #image(MAT.COV)
  is.positive.definite(MAT.COV)
  isSymmetric(MAT.COV)
  
  # Ate aqui é padrão: modelo não separavel
  # Agora vamos aproximar a matriz cheia: matriz permutada
  
  ## MATRIZ DE COVARIÂNCIA PERMUTADA
  
  v <- sub.matriz(MAT.COV,n[ww],n[ww])
  
  # MATRIZ PERMUTADA
  R_SIGMA <- matrix(c(v),ncol=(p^2))
  # n^2 x p^2
  
  # DECOMPOSICAO EM VALORES SINGULARES
  mat.sing <- svd(R_SIGMA)
  vec.R <- sqrt(mat.sing$d[1])*mat.sing$u[,1]
  vec.A <- sqrt(mat.sing$d[1])*mat.sing$v[,1]
  
  # DA DECOMPOSICAO EXTRAIMOS A E R
  A <- matrix(vec.A,ncol=p)
  R <- matrix(vec.R,ncol=n[ww])
  

  
  # Data frame com tempos de execução
  df_exp <- data.frame(#TempoCheia = as.numeric(tempo.cheia), 
    #TempoAprox = as.numeric(tempo.aprox.CH.DMV),
    #ErroPadronizado = as.numeric(standardize.error.index),
    phis = rbind(as.matrix(rep("phi.11", L)),as.matrix(rep("phi.22", L))),

    
    vetor_distancias = rep(d, 2),
    #corrphi11 = unlist(corr.12[ij]),
    #corrphi22 = unlist(corr.22[ij]),
    
    funcao_cor_exp = rbind(as.matrix(unlist(corr.11[ij])), as.matrix(unlist(corr.22[ij]))),
    
    phi22_valor = rep(phi.22, 2)
    
    #incre = incremento,
   # Positiva_definida = is.positive.definite(MAT.COV)
   )
  
  #Salva data frame para cada caso
  
  save(df_exp, 
      file=paste0(#"diff_n", n[ww],
                  "RData/vetor_dist_20221220_v2_cenario_", cenario, "_exponencial.RData"))
  
  paste(print(is.positive.definite(MAT.COV)))
  
 }
}


##### GRÁFICO

#### VETOR DE DISTANCIAS ####

rm(list=ls())


files <- list.files(pattern="vetor_dist_20221220cenario_")
for (i in 1:length(files)) {
  
  # Carregando dataframes
  load(files[i])
  
  # Linha dos graficos
  lcols = c('black', 'black')
  linetype = c('solid', 'dashed')
  
  ########################################################################## 
  # VERSAO PORTUGUES
  ########################################################################## ]
  
  cenario = paste0("Cenário ", i)
  
  titulo_fill = "Alcance espacial"
  legenda_fill = c('Umidade', 'Temperatura')
  
  grafico_pt =  ggplot(data = df, aes(x = vetor_distancias, y = funcao_cor, group = phis))+
   
                  geom_line(aes(linetype = phis, col = phis))  +
                  
                  labs( x = '\nVetor de distâncias',
                        y = "Função de correlação\n",
                        title = 'Análise da correlação das variáveis no espaço',
                        subtitle =  cenario)+
                  
                  
                  scale_color_manual(name = titulo_fill, values = lcols,
                                     labels = legenda_fill) +
                  
                  
                  scale_linetype_manual(name = titulo_fill, values = linetype,
                                        labels = legenda_fill) +     
                  
                  
                  theme(axis.text.x =element_text(size=12),
                        axis.title=element_text(size=14,face="bold"),
                        axis.title.x = element_text(size=30),
                        axis.title.y = element_text(size=18),
                        legend.text.align = 0.5,
                        legend.text = element_text(size = 20),
                        plot.title = element_text(hjust = 0.5))+
                
                  theme_minimal() 
  
  # salva os graficos gerados
  ggsave(
    filename = paste0(cenario, '_v2.jpeg'),
    plot = last_plot(),
    device = 'jpeg',
    #path = NULL,
    scale = 1,
    width = 15,
   height = 9,
    units = c("cm"#,"in",  "mm", "px"
              ),
    dpi = 300,
    limitsize = FALSE,
    bg = NULL)
  
  
  ########################################################################## 
  # VERSAO INGLES
  ########################################################################## 
  
  scenario = paste0("Scenario ", i)
  
  titulo_fill = "Spacial range"
  legenda_fill = c('Humidity', 'Temperature')
  
  grafico_eng = ggplot(data = df, aes(x = vetor_distancias, y = funcao_cor, group = phis))+
    
                  geom_line(aes(linetype=phis, col = phis))  +
                  
                  labs( x = '\nDistance vector',
                        y = "Correlation function\n",
                        title = 'Analysis of the correlation of variables in space',
                        subtitle =  scenario)+
                  
                  
                  scale_color_manual(name = titulo_fill, values = lcols,
                                     labels = legenda_fill) +
                  
                  
                  scale_linetype_manual(name = titulo_fill, values = linetype,
                                        labels = legenda_fill) +     
                  
                  
                  theme(axis.text.x =element_text(size=12),
                        axis.title=element_text(size=14,face="bold"),
                        axis.title.x = element_text(size=30),
                        axis.title.y = element_text(size=18),
                        #legend.text.align = 0.5,
                        legend.text = element_text(size = 20),
                        plot.title = element_text(hjust = 0.5))+
                  
                  theme_minimal() 
  
  
  # salva os graficos gerados
  ggsave(
    filename = paste0(scenario, '_eng.jpeg'),
    plot = last_plot(),
    device = 'jpeg',
    #path = NULL,
    scale = 1,
    width = 15,
    height = 9,
    units = c("cm"#,"in",  "mm", "px"
    ),
    dpi = 300,
    limitsize = TRUE,
    bg = NULL)
  
}

dev.off()
