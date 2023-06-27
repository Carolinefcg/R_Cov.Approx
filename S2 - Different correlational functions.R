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
require(ggplot2)
library(gridExtra)
library(grid)
library(ggpubr)
library(dplyr)

setwd("C:/Users/Caroline/Dropbox/Projeto_AproxSep/Scripts/Caroline/IC_TCC/S2/Rdata")

# cauchy
files <- list.files(pattern="vetor_dist_20220909_cenario")

# exponencial
files_exp <- list.files(pattern="vetor_dist_20221220_v2_cenario_")

cenario = c(rep(1, 2000), rep(2, 2000), rep(3, 2000), rep(4, 2000))


df_cauchy = tibble(NULL)
for (i in 1:4) {
  df_cauchy = rbind(df_cauchy, get(load(files[i])))
}

df_cauchy = cbind(df_cauchy, cenario)


df_exp = tibble(NULL)
for (i in 1:4) {
  df_exp = rbind(df_exp, get(load(files_exp[i])))
}

df_exp = cbind(df_exp, cenario)

names(df_cauchy)
names(df_exp)


dados_exp_cauchy = inner_join(df_cauchy, df_exp, 
                              by = c('phis', 'vetor_distancias', 
                                     'phi22_valor', 'cenario'))

# 1 - 1 2000
# 2 - 2001 4000
# 3- 4001 6000
# 4 -6001 8000

for(i in 1:4){
  print(paste0('inicio ', ((i*2000)-1999), ' fim', i*2000))
}


########################################
# GERAR GRAFICO 
########################################

gerar_grafico = function(idioma){
  theme_set(theme_bw())
  legenda_fill = c('Variável 1', 'Variável 2')
  if (idioma == 'pt'){
    cenario = 'Cenário'
    xlabel ='\nDistância'
    ylabel = "Função de correlação espacial\n"
   
  }else {
      cenario = 'Scenario'
      xlabel ='\nDistance'
      ylabel = 'Spatial correlation function'    }
  for (i in 1:4) {
    # Carregando dataframes
   # load(files[i])
    df = dados_exp_cauchy[((i*2000)-1999):(i*2000),]
    
    # Linha dos graficos
    lcols = c('black', 'black')
    linetype = c('solid', 'dashed')
    cenario_title = paste(cenario, i)
    
    titulo_fill = ""
    
    grafico_pt =  ggplot(data = df, aes(x = vetor_distancias, y = funcao_cor, group = phis))+
      
      geom_line(aes(linetype = phis, col = phis)#, show.legend = FALSE
      )   +
      
      labs( x = xlabel,
            y = ifelse(i==1, ylabel, ""),
            fill = NULL,
            group = NULL,
            title = cenario_title
      )+
      
      scale_color_manual(name = titulo_fill, 
                         values = c('black', 'gray40'),
                         breaks =  c('phi.11', 'phi.22'),
                         labels =legenda_fill
      ) +
      scale_linetype_manual(name = titulo_fill, 
                            values = linetype,
                            breaks =  c('phi.11', 'phi.22'),
                            labels =legenda_fill
      ) +   
      
      theme(axis.text.x =element_text(size=13),
            axis.title=element_text(size=13#,face="bold"
                                    ),
            axis.title.x = element_text(size=13),
            axis.title.y = element_text(size=13),
            # legend.text.align = 0.5,
            legend.text = element_text(size = 14),
            #plot.title = element_text(hjust = 0.5)
            #,legend.position = 'none'
            )
    
    # Seprando para o grid
    if (i == 1) {
      g1 = grafico_pt
      
    }else if(i == 2) {
      g2 = grafico_pt
      
    }else if(i == 3){
      g3 = grafico_pt
      
    }else{
      g4 = grafico_pt
    }
  }  
  return(ggarrange(g1, g2, g3, g4, common.legend = TRUE, legend="bottom"))
}

(fig = gerar_grafico('pt'))

########################################
# GERAR GRAFICO 1X4
########################################

gerar_grafico_colunas = function(idioma){
  theme_set(theme_bw())
  legenda_fill = c('Variável 1', 'Variável 2')
  if (idioma == 'pt'){
    cenario = 'Cenário'
    xlabel =''#'\nDistância'
    ylabel =''#"Função de correlação espacial\n"
    
  }else {
    cenario = 'Scenario'
    xlabel ='\nDistance'
    ylabel = 'Spatial correlation function'    }
  for (i in 1:length(files)) {
    # Carregando dataframes
    load(files[i])
    
    # Linha dos graficos
    lcols = c('black', 'black')
    linetype = c('solid', 'dashed')
    cenario_title = paste(cenario, i)
    
    titulo_fill = ""
    
    grafico_pt =  ggplot(data = df, aes(x = vetor_distancias,   y = funcao_cor, group = phis))+
      
      geom_line(aes(linetype = phis
                    , col = phis)#, show.legend = FALSE
      ,size = 1.5)   +
      
      labs( x = xlabel,
            y =ifelse(i ==1,  ylabel, ''),
            fill = NULL,
            group = NULL,
            title = cenario_title
      )+
      
      scale_color_manual(name = titulo_fill, 
                         values = c('black', 'gray40'),
                         breaks =  c('phi.11', 'phi.22'),
                         labels =legenda_fill
      ) +
      scale_linetype_manual(name = titulo_fill, 
                            values = linetype,
                            breaks =  c('phi.11', 'phi.22'),
                            labels =legenda_fill
      ) +   
      
      theme(axis.text.x =element_text(size=15),
            axis.text.y =element_text(size=15),
            axis.title.y = element_text(size=18),
            axis.title=element_text(size=18),#,face="bold",
            plot.title = element_text(size=18),#,face="bold"
           
            # legend.text.align = 0.5,
            legend.text = element_text(size = 18)
            #plot.title = element_text(hjust = 0.5)
            #,legend.position = 'none'
      )
    
    # Seprando para o grid
    if (i == 1) {
      g1 = grafico_pt
      
    }else if(i == 2) {
      g2 = grafico_pt
      
    }else if(i == 3){
      g3 = grafico_pt
      
    }else{
      g4 = grafico_pt
    }
  }  
  return(annotate_figure(ggarrange(g1, g2, g3, g4, 
                          common.legend = TRUE,
                          legend="bottom", ncol=4),
                          left = text_grob("Função de correlação espacial", 
                          rot = 90,# vjust = 2,
                          size =20,
                          hjust = .4
                          ), #gp = gpar(cex = 1.6)
                          bottom = textGrob("Distância", 
                          gp = gpar(cex = 1.6), vjust = -3)))
}
#                                  top = textGrob("Análise da função de correlação em cada cenário",gp=gpar(fontsize=15,font=1))))
(fig = gerar_grafico_colunas('pt'))
(fig = gerar_grafico_colunas('pt'))

########################################
# GERAR GRAFICO SEPARADO
########################################

gerar_grafico_separado = function(idioma, cenario_val){
  theme_set(theme_bw())
  if (idioma == 'pt'){
    cenario = 'Cenário'
    xlabel ='\nDistância'
    ylabel = ifelse(cenario_val ==1, "Função de correlação espacial\n", "")
    legenda_fill = c('Variável 1', 'Variável 2')
    
  }else {
    cenario = 'Scenario'
    xlabel ='\nDistance'
    ylabel = 'Spatial correlation function' 
    legenda_fill = c('Variable 1', 'Variable 2')}
  
  k= cenario_val 
  for (i in k:length(files)) {
    # Carregando dataframes
    #df = get(load(files_exp[i]))
    df = dados_exp_cauchy[((i*2000)-1999):(i*2000),]
    
    # Linha dos graficos
    lcols = c('black', 'black')
    linetype = c('solid', 'dashed')
    cenario_title = paste(cenario, i)
    
    titulo_fill = ""
    
    grafico_pt =  ggplot(data = df, aes(x = vetor_distancias, 
                                       # y = funcao_cor_exp,
                                        group = phis))+
      
      geom_line(aes(y = funcao_cor, linetype = phis#, col = phis
                    ,color='black'
                    )#, show.legend = FALSE
      )   +

      
      geom_line(aes(y = funcao_cor_exp,  linetype = phis,
                    #col =  phis,

                    color='gray30'))   +
      
      
      labs( x = xlabel,
            y = ylabel,
            fill = NULL,
            group = NULL,
            title = cenario_title
      )+
      
      scale_color_manual(name = titulo_fill,
                        values = c('black', 'gray30'),
                         #breaks =  c('phi.11', 'phi.22', 'phi.11', 'phi.22')
                         labels = c('Exp','Cauchy')
      ) +
      scale_linetype_manual(name = titulo_fill,
                            values = linetype,
                            breaks =  c('phi.11', 'phi.22'),
                            labels = legenda_fill
                            #,colours = c('a', 'b')
      ) +
      
      theme(axis.text.x =element_text(size=13),
            axis.title=element_text(size=13#,face="bold"
            ),
            axis.title.x = element_text(size=13),
            axis.title.y = element_text(size=13),
            # legend.text.align = 0.5,
            legend.text = element_text(size = 14),
            #plot.title = element_text(hjust = 0.5)
            legend.position = 'bottom'
      )
  return(grafico_pt)
  }
}
#dev.off()

(fig = gerar_grafico_separado('pt', 2))

########################################
# GERAR GRAFICO 2X2
########################################


gerar_grafico_2_2 = function(idioma){
  theme_set(theme_bw())
  legenda_fill = c('Variável 1', 'Variável 2')
  if (idioma == 'pt'){
    cenario = 'Cenário'
    xlabel =''#'\nDistância'
    ylabel =''#"Função de correlação espacial\n"
    
  }else {
    cenario = 'Scenario'
    xlabel ='\nDistance'
    ylabel = 'Spatial correlation function'    }
  for (i in 1:length(files)) {
    # Carregando dataframes
    load(files[i])
    
    # Linha dos graficos
    lcols = c('black', 'black')
    linetype = c('solid', 'dashed')
    cenario_title = paste(cenario, i)
    
    titulo_fill = ""
    
    grafico_pt =  ggplot(data = df, aes(x = vetor_distancias,   y = funcao_cor, group = phis))+
      
      geom_line(aes(linetype = phis
                    , col = phis)#, show.legend = FALSE
                ,size = 1.5)   +
      
      labs( x = xlabel,
            y =ifelse(i ==1,  ylabel, ''),
            fill = NULL,
            group = NULL,
            title = cenario_title
      )+
      
      scale_color_manual(name = titulo_fill, 
                         values = c('black', 'gray40'),
                         breaks =  c('phi.11', 'phi.22'),
                         labels =legenda_fill
      ) +
      scale_linetype_manual(name = titulo_fill, 
                            values = linetype,
                            breaks =  c('phi.11', 'phi.22'),
                            labels =legenda_fill
      ) +   
      
      theme(axis.text.x =element_text(size=15),
            axis.text.y =element_text(size=15),
            axis.title.y = element_text(size=18),
            axis.title=element_text(size=18),#,face="bold",
            plot.title = element_text(size=18),#,face="bold"
            
            # legend.text.align = 0.5,
            legend.text = element_text(size = 18)
            #plot.title = element_text(hjust = 0.5)
            #,legend.position = 'none'
      )
    
    # Seprando para o grid
    if (i == 1) {
      g1 = grafico_pt
      
    }else if(i == 2) {
      g2 = grafico_pt
      
    }else if(i == 3){
      g3 = grafico_pt
      
    }else{
      g4 = grafico_pt
    }
  }  
  return(annotate_figure(ggarrange(g1, g2, g3, g4, 
                                   common.legend = TRUE,
                                   legend="bottom", ncol=2, nrow =2),
                         left = text_grob("Função de correlação espacial", 
                                          rot = 90,# vjust = 2,
                                          size =20,
                                          hjust = .4
                         ), #gp = gpar(cex = 1.6)
                         bottom = textGrob("Distância", 
                                           gp = gpar(cex = 1.6), vjust = -3)))
}
(fig = gerar_grafico_2_2('pt'))
(fig = gerar_grafico_colunas('pt'))

########################################
# SALVAR GRAFICO COLUNAS
########################################

salvar_grafico_colunas = function(idioma, formato){
  ggsave(
    filename = (paste0('vetor_correlacao_', idioma,'_v4.', formato)),
    plot = gerar_grafico_colunas(idioma),
    device = formato,
    path = paste0("C:/Users/Caroline/Dropbox/Projeto_AproxSep/Scripts/Caroline/IC_TCC/S2/Figuras/", idioma),
    scale = 0.85,
    width = 30,
    height = 16,
    units = c("cm"#,"in",  "mm", "px"
    ),
    dpi = 600,
    limitsize = TRUE,
    bg = NULL)
  dev.off()
}
salvar_grafico_colunas('pt', 'jpeg')

salvar_grafico_2_2 = function(idioma, formato){
  ggsave(
    filename = (paste0('vetor_correlacao_', idioma,'_v5.', formato)),
    plot = gerar_grafico_2_2(idioma),
    device = formato,
    path = paste0("C:/Users/Caroline/Dropbox/Projeto_AproxSep/Scripts/Caroline/IC_TCC/S2/Figuras/", idioma),
    scale = 0.85,
    width = 20,
    height = 20,
    units = c("cm"#,"in",  "mm", "px"
    ),
    dpi = 600,
    limitsize = TRUE,
    bg = NULL)
  dev.off()
}
salvar_grafico_2_2('pt', 'pdf')

########################################
# SALVAR GRAFICO 
########################################

salvar_grafico('pt', 'jpeg')

salvar_grafico = function(idioma, formato){
  ggsave(
    filename = (paste0('vetor_correlacao_', idioma,'_v3.', formato)),
    plot = gerar_grafico(idioma),
    device = formato,
    path = paste0("C:/Users/Caroline/Dropbox/Projeto_AproxSep/Scripts/Caroline/IC_TCC/S2/Figuras/", idioma),
    scale = 0.85,
    width = 18,
    height = 18.5,
    units = c("cm"#,"in",  "mm", "px"
    ),
    dpi = 600,
    limitsize = TRUE,
    bg = NULL)
  dev.off()
}

salvar_grafico('pt', 'pdf')


########################################
# SALVAR GRAFICO SEPARADO
########################################

salvar_grafico_separado = function(idioma, formato, cenario){
  ggsave(
    filename = (paste0('vetor_correlacao_cauchy_exp_cenario', cenario,'_', idioma,'.', formato)),
    plot = gerar_grafico_separado(idioma, cenario),
    device = formato,
    path = paste0("C:/Users/Caroline/Dropbox/Projeto_AproxSep/Scripts/Caroline/IC_TCC/S2/Figuras/", idioma),
    scale = .6,
    width = 23,
    height = 17,
    units = c("cm"#,"in",  "mm", "px"
    ),
    dpi = 600,
    limitsize = TRUE,
    bg = NULL)
  dev.off()
  gerar_grafico_separado('pt', cenario)
}

salvar_grafico_separado('pt', 'pdf',1)
salvar_grafico_separado('pt', 'pdf',2)
salvar_grafico_separado('pt', 'pdf',3)
salvar_grafico_separado('pt', 'pdf',4)

