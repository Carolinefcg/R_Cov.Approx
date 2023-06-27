##### GRÁFICO #####

# v <- read.csv2(file="C:/Users/rafae_000/Dropbox/UFRJ/PROJETO/Reuniões/2016/5- Maio/24-05-2016/variacoes.csv",header=TRUE)
rm(list=ls())

#rodar fora do Rproject
setwd("C:/Users/Caroline/Dropbox/Projeto_AproxSep/Scripts/Caroline/IC_TCC/S1/Rdata")
files <- list.files(pattern="tempos_diff_20220908_n")

cenario = c()
n_localizacoes = c()
tempo_cheia = c()
tempo_aprox = c()
reducao = c()

nome_arquivo = c()

sapply(1:length(files), function(x){
  load(files[x])
  
  inicio_cenario = unlist(gregexpr('cenario', files[x]))
  fim_cenario = unlist(gregexpr('cenario', files[x]))+7
  
  inicio_n = unlist(gregexpr('n', files[x]))[1]+1
  fim_n = unlist(gregexpr('0', files[x]))[length(unlist(gregexpr('0', files[x])))]
  
  cenario<<- c(cenario, substr((files[x]), start = inicio_cenario, stop = fim_cenario))
  n_localizacoes<<- c(n_localizacoes, as.numeric( substr((files[x]), start = inicio_n, stop = fim_n)))
  tempo_cheia <<- c(tempo_cheia, tempos$TempoCheia)
  tempo_aprox <<- c(tempo_aprox, tempos$TempoAprox)
  reducao <<- c(reducao, as.numeric((1-(tempos$TempoAprox/tempos$TempoCheia))*100))
  nome_arquivo <<- c(nome_arquivo, files[x])
  
})

v = data.frame(reducao, nome_arquivo, cenario, n_localizacoes,tempo_cheia, tempo_aprox)


v = v %>%
  mutate(reducao_cenario1 = ifelse(cenario == 'cenario1', reducao, NA),
         reducao_cenario2 = ifelse(cenario == 'cenario2', reducao, NA),
         reducao_cenario3 = ifelse(cenario == 'cenario3', reducao, NA),
         reducao_cenario4 = ifelse(cenario == 'cenario4', reducao, NA),
         n_c1 = ifelse(cenario == 'cenario1', n_localizacoes, NA),
         n_c2 = ifelse(cenario == 'cenario2', n_localizacoes, NA),
         n_c3 = ifelse(cenario == 'cenario3', n_localizacoes, NA),
         n_c4 = ifelse(cenario == 'cenario4', n_localizacoes, NA)) %>%
  arrange(n_localizacoes)


plot(v$n_localizacoes, v$reducao)


dev.off()

library(dplyr)

c1 = v %>%
  dplyr::select(tempo_cheia, tempo_aprox,reducao_cenario1, n_c1) %>%
  filter(!is.na(n_c1))

c2 = v %>%
  dplyr::select(tempo_cheia, tempo_aprox,reducao_cenario2, n_c2) %>%
  filter(!is.na(n_c2))

c3 = v %>%
  dplyr::select(tempo_cheia, tempo_aprox,reducao_cenario3, n_c3) %>%
  filter(!is.na(n_c3))

c4 = v %>%
  dplyr::select(tempo_cheia, tempo_aprox,reducao_cenario4, n_c4) %>%
  filter(!is.na(n_c4))

#### GRAFICOS #####
setwd("C:/Users/Caroline/Dropbox/Projeto_AproxSep/Scripts/Caroline/IC_TCC/S1/Figuras")


pdf(paste0("tempo_cenario_pt_v2.pdf",sep=""),
    width=35, # era 25
    height=23)


jpeg(file =  paste0("tempo_cenario_pt_v2.jpeg",sep=""),
     quality=100, width=1100, height=710)

par(mar=c(4.5,4.5,0.5,0.5),cex=2.8) # jpeg 2.8 pdf 6
#
#### GRAFICO PT ####
plot(c1$reducao_cenario1~c1$n_c1, type="b", bty="l",
     ylim=c(0,100),xlab="Localizações",ylab="Redução do tempo (em %)",
     lwd=2 , pch=20, col = '#000000', xaxt='n', axes = FALSE)
axis(1, # para x
     seq(0,1000,by=100),
     labels = seq(0,1000,by=100),
     cex.axis=0.85)
axis(2, # para y
     seq(0,100,by=5),
     labels = seq(0,100,by=5),
     cex.axis=0.85)
lines(c2$reducao_cenario2~c2$n_c2,type="b",pch=8,col = 'gray20',lwd=2)
lines(c3$reducao_cenario3~c3$n_c3,type="b",pch=17,col = 'gray60',lwd=2)
lines(c4$reducao_cenario4~c4$n_c4,type="b",pch=15,col = 'gray40',lwd=2)
legend("bottomright",legend=expression('Cenário 1','Cenário 2',
                                       'Cenário 3', 'Cenário 4'),
       pch=c(20,8,17,15), 
       col = c('#000000','gray20','gray60', 'gray40'))

dev.off()


#
#### GRAFICO ENG ####

plot(c1$reducao_cenario1~c1$n_c1, type="b", bty="l",
     ylim=c(1,100),xlab="Spatial locations",ylab="Time reduction (%)",
     lwd=2 , pch=20, col = '#000000', xaxt='n', axes = FALSE)
axis(1, # para x
     seq(0,1000,by=100),
     labels = seq(0,1000,by=100),
     cex.axis=0.85)
axis(2, # para y
     seq(0,100,by=5),
     labels = seq(0,100,by=5),
     cex.axis=0.85)
lines(c2$reducao_cenario2~c2$n_c2,type="b",pch=8,col = 'gray20',lwd=2)
lines(c3$reducao_cenario3~c3$n_c3,type="b",pch=17,col = 'gray60',lwd=2)
lines(c4$reducao_cenario4~c4$n_c4,type="b",pch=15,col = 'gray40',lwd=2)
legend("bottomright",legend=expression('Scenario 1','Scenario 2', 'Scenario 3', 'Scenario 4'),
       pch=c(20,8,17,15), 
       col = c('#000000','gray20','gray60', 'gray40'))

dev.off()