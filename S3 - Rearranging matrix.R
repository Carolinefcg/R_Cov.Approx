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
require(gdata) #function keep

#############################################################

setwd("C:/Users/Caroline/Dropbox/Projeto_AproxSep/Scripts/Caroline/IC_TCC/S3/Rdata")

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


# Retirar o FOR de fora -- OK
# Tirar o p que tava com colchete [qq] -- OK
# Substituir FORs por APPLYs -- OK
# Para a comparacao ser justa, o FOR final nao deveria ser um apos o outro
# Fazer um FOR só pra aproximada e outro para a cheia
# Pois o R deveria estar igual em termos de memoria
# Para ele calcular de fato só o tempo da verossimilhanca

# Ler a media, o R, o A, a MAT.COV
# Calcular pra cheia
# Limpar TUDO
# Ler de novo
# Calcular da aproximada

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

## Funçao que define a distância entre dois pontos no espaco

distancia <- function(x1,y1,x2,y2) {
  dist = sqrt((x1-x2)^2+(y1-y2)^2)
  return(dist)
}

save(sub.matriz, dmvnorm.chol, distancia, 
     file="VeroEnvFunctions_teste.RData")
ww = 1
ij = 1
incremento = -0.01
standardize.error.index = 0

# Numero de localizacoes
n <- c(100)

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

while( standardize.error.index >= 0){
  
  #print(standardize.error.index)
  
  # Cleans environment
  rm(list=ls()[!(ls() %in% c("ww",#"ij", 
                             "incremento", "standardize.error.index", "H"))])
  
  # Loads necessary functions
  load("VeroEnvFunctions_teste.RData")
  
  # incremento para variar cenarios
  
  incremento = incremento + 0.01
  #incremento = 0.4
  
  print(incremento)
  
  phi.11 <- 0.1
  phi.22 <- phi.11 + incremento
  cenario <- phi.11 + incremento
  # phi
  if(phi.22  > 0.71){
    break
     }
  # Quero comparar o tempo de todas as execucoes da verossimilhanca cheia e aproximada
  # Para todas as configuracoes (combinacoes) de n e p
  # qq = tamanho do vetor p
   qq=1
  # ww = tamanho do vetor n
   ww=1
  #set.seed(117054003)
  # Numero de variaveis
  p <- 2
  # Numero de localizacoes
  n <- c(200)
  # Numero de amostras que quero --replicas
  L <- 1
  
  # spatial range - alcance espacial
  
  # Media usada foi 0: latitude e longitude nao estao entrando como covariaveis
  MU.vetor <- rep(0,(n[ww]*p))
  
  phi.12 <- (phi.11+phi.22)/2
  # phi.12 <- sqrt((phi.11^2+phi.22^2)//2)
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

  
  rho.12 <- 0.4 # fixar este valor, baseado na restricao do rho na funcao Cauchy
  #rho.12 <- sqrt((phi.11)*(phi.22))//(phi.12)
    
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
  # sigma_i*sigma_j*rho_ij*exp(-h//phi_ij) -- funcao de correlacao espacial = exp
  # para uma variavel c ela msm: vai dar 1 tds os elementos do produto
  sapply(1:p, 
         function(i){
           sapply(i:p, 
                  function(j){
                    # Basta gerar sigma pelo produto de sigma_i vezes sigma_j
                    # Constroi os 4 blocos
                    MAT.COV[(((i-1)*n[ww])+1):(i*n[ww]),(((j-1)*n[ww])+1):(j*n[ww])] <<- sigma[i,j]*rho[i,j]*((1+(H/phi[i,j]))^(-1))
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
  
  
  
  standardize.error.index <- (norm(MAT.COV - kronecker(A,R), type= "F")/norm(MAT.COV, type= "F"))/sqrt(1 - (1/(p^2)))
  
  
  
  
  # São negativas definidas pois todos os autovalores sao negativos
  # eigen(A)
  # Error in eigen(R) : matriz não quadrada em 'eigen'
  #eigen(R)$values 
  
  # Geramos o Y verdadeiro utilizando matriz cheia, nao aproximada
  is.positive.definite(MAT.COV)
  
  
  
 # erro <- rmvnorm(L,rep(0,n[ww]*p),MAT.COV)
  #Y.vetor <- MU.vetor + t(erro)
  # Y é n*p x 200 
  # No nosso caso 100 x 200 pois temos 50 localizacoes
  # n primeiras linhas é a primera var 
  # e as n ultimas sao as demais
  
  # 2n x 2n
  # Vetor Repetido L veze
  
  
  # Data frame com tempos de execução
  df <- data.frame(#TempoCheia = as.numeric(tempo.cheia), 
                       #TempoAprox = as.numeric(tempo.aprox.CH.DMV),
                       ErroPadronizado = as.numeric(standardize.error.index),
                       phi22 = phi.22,
                       incre = incremento,
                       Positiva_definida = is.positive.definite(MAT.COV))
  
  # Salva data frame para cada caso
  save(df, 
       file=paste0("diff_n", n[ww], "_cenario_phi", cenario, "_cauchy.RData"))
 
   paste(print(incremento),print(is.positive.definite(MAT.COV)))
 
}
#}


##### GRÁFICO

#### ERRO PADRONIZADO ####

rm(list=ls())


files <- list.files(pattern="diff_n100_cenario")

cenarios_std_error <- c()
#cenarios_phi22 <- c()
cenarios_incre <- c()


sapply(1:length(files), function(x){
  load(files[x])
  cenarios_std_error<<- c(cenarios_std_error, df$ErroPadronizado)
  #cenarios_phi22 <<- c(cenarios_phi22, tempos$phi22)
  cenarios_incre<<- c(cenarios_incre, df$incre)
})

(df2 = data.frame(cenarios_std_error, cenarios_incre
                 #cenarios_phi22
                 )) 

require(ggplot2)
require(hrbrthemes)


theme_set(theme_bw())

grafico_pt = ggplot(data = df2, aes(x = cenarios_incre, y = cenarios_std_error))+
  geom_line(color="black", size=2, alpha=0.9) +
  labs( x = '\nIncremento no alcance espacial'#(expression(phi[22]))
        , y = "Erro padronizado\n")+
  theme(axis.text.x =element_text(size=20),
        axis.text.y =element_text(size=20),
        #axis.title=element_text(size=14,face="bold"),
          axis.title.x = element_text(size=23),
         axis.title.y = element_text(size=23)
        #legend.text.align = 0.5,
        )
#+ggtitle("Erro padrozinado com incremento em phi|Distribuição Cauchy")
#standard error

grafico_pt

grafico_eng = ggplot(data = df2, aes(x = cenarios_incre, y = cenarios_std_error))+
  geom_line(color="black", size=2, alpha=0.9) +
  labs( x = 'Increment'#(expression(phi[22]))
        , y = "Standard error")+
  theme(#axis.text.x =element_text(size=12),
    #axis.title=element_text(size=14,face="bold"),
    axis.title.x = element_text(size=17.5),
    axis.title.y = element_text(size=17.5)
    #legend.text.align = 0.5,
  )

grafico_eng

salvar_grafico = function(grafico, idioma, formato){
  ggsave(
    filename = (paste0('erro_incremento_', idioma,'.', formato)),
    plot = grafico,
    device = formato,
    path = paste0("C:/Users/Caroline/Dropbox/Projeto_AproxSep/Scripts/Caroline/IC_TCC/S3/Figuras/", idioma),
    scale = 0.65,
    width = 25,
    height = 20,
    units = c("cm"#,"in",  "mm", "px"
    ),
    dpi = 600,
    limitsize = TRUE,
    bg = NULL)
  #dev.off()
}

salvar_grafico(grafico_eng, 'eng', 'jpeg')

salvar_grafico_v2 = function(grafico, idioma, formato){
  ggsave(
    filename = (paste0('erro_incremento_', idioma,'_v2.', formato)),
    plot = grafico,
    device = formato,
    path = paste0("C:/Users/Caroline/Dropbox/Projeto_AproxSep/Scripts/Caroline/IC_TCC/S3/Figuras/", idioma),
    scale = 0.65,# 55 para jpeg
    width = 45,
    height = 22,
    units = c("cm"#,"in",  "mm", "px"
    ),
    dpi = 600,
    limitsize = TRUE,
    bg = NULL)
  #dev.off()
}

salvar_grafico_v2(grafico_pt, 'pt', 'pdf')
