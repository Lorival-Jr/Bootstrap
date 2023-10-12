# 2° Trabalho - Parte 1 ---------------------------------------------------

# Método Bootstrap -------------
# Ex.1 Para esta atividade será considerado a base de dados "Abalone" 

# Determine a média e mediana das variáveis "Whole weight", "Diameter" e "Lenght". 
# Determine o estimativa e o erro-padrão da média e mediana via método Bootstrap considerando
# - beta = 1000000 de replicações 
# - Tamanhos das amostras bootstrap, 20,40,...,200
# - Grafico histograma

# Determinar o vício e o eqm das estimativas bootstrap

# Determinar os intervalos de 90%, 95%, 99% via bootstrap usando os tipos:
# Normal, Básico (percentil), t de Student, BCA

# Determinar o tempo computacional dos processo
# Escrever relatório em Rmarkdown

# Obs: não deve ser usado pacotes

# Colunas da base de dados: Importa o banco e ve

#sex, length, height, whole weight, shucked weight, viscera weight, shell weight



abalone <- read.csv('Abalone.txt')
str(abalone)

names(abalone) <- c('sex', 'length', 'diameter', 'height', 'whole_weight', 'shucked_weight', 'viscera_weight', 'shell_weight', 'rings')

head(abalone)



a <- c(1,2,10,50,1,1,90,100,1)
b <- quantile(a, c(0.5,0.95))



# Função pika -------------------------------------------------------------

# calcular nessa função também os ics e retornar uma lista
Bootstrap <- function(N,enes,var,banco= abalone, conf.level=0.95, seed=NULL)
{
  est_bootstrap <- data.frame('media' = rep(NA, 10), 'mediana' = NA, 'erro_media' = NA, 'erro_mediana' = NA, 'vicio_media' = NA, 'vicio_mediana' = NA, 'EQM_media' = NA, 'EQM_mediana' = NA,  'n' = NA)
  
  ICs_media   <-  data.frame('Tipo_IC' = rep(NA,40*length(conf.level)), 'LI'= NA, 'LS' = NA, n = NA)
  ICs_mediana <-  data.frame('Tipo_IC' = rep(NA,40*length(conf.level)), 'LI'= NA, 'LS' = NA, n = NA)
    
    
  media_global   <- mean(banco[[var]])
  mediana_global <- median(banco[[var]])
  
  if(!is.null(seed)){set.seed(seed)}
  
  
  amostras <- lapply(enes, banco = banco, var = var, N = N,
                     FUN = function(x, banco, var,N){
                       index <- round(runif(x*N, 0.5, 0.4 + length(banco[[var]])))
                       amostra <- banco[[var]][index]
                       return(amostra)
                     }
  )
  
  for(i in enes)
  {
    
    estatistiscas <- lapply(split(amostras[[i/enes[1]]], ceiling(seq_along(amostras[[i/enes[1]]])/i)), 
                            FUN = function(x)
                            {return(c(mean(x), median(x)))})

    estatistiscas <- do.call(rbind, estatistiscas)

    valores <-  c(mean(estatistiscas[,1]), median(estatistiscas[,2]), # media e mediana
                  sd(estatistiscas[,1]), # Erro media
                  sd(estatistiscas[,2]), # Erro mediana
                  (mean(estatistiscas[,1]) - media_global),           # Vicio media
                  (median(estatistiscas[,2]) - mediana_global),       # Vicio mediana
                  (mean(estatistiscas[,1]) - media_global)^2,         # Eqm media
                 (median(estatistiscas[,2]) - mediana_global)^2,       # Eqm mediana
                  i
                  
    )

    est_bootstrap[i/enes[1], ] <- valores
    
    # Intervalos de confiança ---
    for(conf_index in 1:length(conf.level)){
      
      conf <- conf.level[conf_index]
        # Normal
        IC_norm_mean   <- c(valores[1] - qnorm(conf)*valores[3], valores[1] + qnorm(conf)*valores[3])
        IC_norm_median <- c(valores[2] - qnorm(conf)*valores[4], valores[2] + qnorm(conf)*valores[4])
        
        # Percentil Reverso (básico)
        Q1 <- quantile(estatistiscas[,1], c(1-conf, conf))
        Q2 <- quantile(estatistiscas[,2], c(1-conf, conf))
        
        IC_perc_mean <-   c(2*valores[1] - Q1[[2]], 2*valores[1] - Q1[[1]])
        IC_perc_median <- c(2*valores[2] - Q2[[2]], 2*valores[2] - Q2[[1]])
        
        # t-Student
        
        IC_t_mean   <- c(valores[1] - qt(conf, i -1 )*valores[3],
                         valores[1] + qt(conf, i -1)*valores[3])
        IC_t_median <- c(valores[2] - qt(conf, i -1)*valores[4],
                         valores[2] + qt(conf, i -1)*valores[4])   
        
        # BCA
        IC_bca_mean   <- c(NA, NA)
        IC_bca_median <- c(NA, NA)
        
        # jutando todos
        
        norm_name <- paste0('Normal_', conf)
        perc_name <- paste0('Percentil_', conf)
        t_name <- paste0('t_', conf)
        bca_name <- paste0('BCA_', conf)
        
        ICs_media[i/enes[1] + (conf_index -1)*40, ]       <- c(norm_name,IC_norm_mean,i)
        ICs_media[i/enes[1] +10 + (conf_index -1)*40, ]   <- c(perc_name,IC_perc_mean,i)   
        ICs_media[i/enes[1] +20 + (conf_index -1)*40, ]   <- c(t_name,IC_t_mean,i)
        ICs_media[i/enes[1] +30 + (conf_index -1)*40, ]   <- c(bca_name,IC_bca_mean,i)
        
        ICs_mediana[i/enes[1] + (conf_index -1)*40, ]       <- c(norm_name,IC_norm_median,i)
        ICs_mediana[i/enes[1] +10 + (conf_index -1)*40, ]   <- c(perc_name,IC_perc_median,i)   
        ICs_mediana[i/enes[1] +20 + (conf_index -1)*40, ]   <- c(t_name,IC_t_median,i)
        ICs_mediana[i/enes[1] +30 + (conf_index -1)*40, ]   <- c(bca_name,IC_bca_median,i)
    }
    
    jpeg(filename = paste0('img/tend_central_', i, '.jpg' ), width = 1080, height = 720, quality = 100)
    par(mfrow = c(1,2), mar = c(5, 5, 4, 6))
    hist(estatistiscas[,1], xlab = paste('Valores da média com n=', i), 
         ylab = 'Frequência', main = '',  prob = TRUE)
    abline(v = mean(estatistiscas[,1]), col = 'red', lty = 2)
    hist(estatistiscas[,2], xlab = paste('Valores da mediana com n=', i), 
         ylab = 'Frequência', main = '', prob = TRUE)
    abline(v = median(estatistiscas[,2]), col = 'blue', lty = 2)
    
    legend(x = 0, y= 0,
           legend = c('Estimativa da Média (Bootstrap)', 'Estimativa da Mediana (Bootstrap)'),
           col = c('red','blue'),
           lty = 2,
           bty = 'n',
           xpd = T)
    
    par(mfrow = c(1,1))
    
    dev.off()
    
    cat('Repetição:',i/enes[1], '\n')
    
    
  }
  
  output <- list('estimativas' = est_bootstrap, 'IC_media' = ICs_media, 'IC_mediana' = ICs_mediana)
  
  return(output)
  
}


# diametrer ----------------------------------------------------------------

# Não foi possível usar N= 1.000.000 por falta de RAM, foi utilizado 700.000 repetições
diameter <- Bootstrap(10000, seq(20,200,20), 'diameter', conf.level =  c(0.9, 0.95, 0.99))


as.numeric(diameter$IC_media$LI[1:40]) < mean(abalone$diameter) & mean(abalone$diameter) < as.numeric(diameter$IC_media$LS[1:40])
# OU seja com 0.9 de confiança todos intervalos contém o parametro verdadeiro da media

as.numeric(diameter$IC_media$LI[41:80]) < mean(abalone$diameter) & mean(abalone$diameter) < as.numeric(diameter$IC_media$LS[41:80])
# OU seja com 0.95 de confiança todos intervalos contém o parametro verdadeiro da media

as.numeric(diameter$IC_media$LI[81:120]) < mean(abalone$diameter) & mean(abalone$diameter) < as.numeric(diameter$IC_media$LS[81:120])
# OU seja com 0.99 de confiança todos intervalos contém o parametro verdadeiro da media

as.numeric(diameter$IC_mediana$LI[1:40]) < median(abalone$diameter) &
  mean(abalone$diameter) < as.numeric(diameter$IC_mediana$LS[1:40])
# OU seja com 0.9 de confiança todos intervalos contém o parametro verdadeiro da mediana

as.numeric(diameter$IC_mediana$LI[41:80]) < median(abalone$diameter) &
  mean(abalone$diameter) < as.numeric(diameter$IC_mediana$LS[41:80])
# OU seja com 0.95 de confiança todos intervalos contém o parametro verdadeiro da mediana

as.numeric(diameter$IC_mediana$LI[81:120]) < median(abalone$diameter) &
  mean(abalone$diameter) < as.numeric(diameter$IC_mediana$LS[81:120])
# OU seja com 0.99 de confiança todos intervalos contém o parametro verdadeiro da mediana


# Whole.weight ------------------------------------------------------------



whole.weight <- Bootstrap(10000, seq(20,200,20), 'whole_weight', conf.level =  c(0.9, 0.95, 0.99))


as.numeric(whole.weight$IC_media$LI[1:40]) < mean(abalone$whole_weight) & mean(abalone$whole_weight) < as.numeric(whole.weight$IC_media$LS[1:40])
# OU seja com 0.9 de confiança todos intervalos contém o parametro verdadeiro da media

as.numeric(whole.weight$IC_media$LI[41:80]) < mean(abalone$whole_weight) & mean(abalone$whole_weight) < as.numeric(whole.weight$IC_media$LS[41:80])
# OU seja com 0.95 de confiança todos intervalos contém o parametro verdadeiro da media

as.numeric(whole.weight$IC_media$LI[81:120]) < mean(abalone$whole_weight) & mean(abalone$whole_weight) < as.numeric(whole.weight$IC_media$LS[81:120])
# OU seja com 0.99 de confiança todos intervalos contém o parametro verdadeiro da media



as.numeric(whole.weight$IC_mediana$LI[1:40]) < median(abalone$whole_weight) &
  mean(abalone$whole_weight) < as.numeric(whole.weight$IC_mediana$LS[1:40])
# OU seja com 0.9 de confiança todos intervalos contém o parametro verdadeiro da mediana

as.numeric(whole.weight$IC_mediana$LI[41:80]) < median(abalone$whole_weight) &
  mean(abalone$whole_weight) < as.numeric(whole.weight$IC_mediana$LS[41:80])
# OU seja com 0.95 de confiança todos intervalos contém o parametro verdadeiro da mediana

as.numeric(whole.weight$IC_mediana$LI[81:120]) < median(abalone$whole_weight) &
  mean(abalone$whole_weight) < as.numeric(whole.weight$IC_mediana$LS[81:120])
# OU seja com 0.99 de confiança todos intervalos contém o parametro verdadeiro da mediana


# lenght ------------------------------------------------------------------

length <- Bootstrap(10000, seq(20,200,20), 'length', conf.level =  c(0.9, 0.95, 0.99))


as.numeric(length$IC_media$LI[1:40]) < mean(abalone$length) &
  mean(abalone$length) < as.numeric(length$IC_media$LS[1:40])
# OU seja com 0.9 de confiança todos intervalos contém o parametro verdadeiro da media

as.numeric(length$IC_media$LI[41:80]) < mean(abalone$length) &
  mean(abalone$length) < as.numeric(length$IC_media$LS[41:80])
# OU seja com 0.95 de confiança todos intervalos contém o parametro verdadeiro da media

as.numeric(length$IC_media$LI[81:120]) < mean(abalone$length) &
  mean(abalone$length) < as.numeric(length$IC_media$LS[81:120])
# OU seja com 0.99 de confiança todos intervalos contém o parametro verdadeiro da media

as.numeric(length$IC_mediana$LI[1:40]) < median(abalone$length) &
  mean(abalone$length) < as.numeric(length$IC_mediana$LS[1:40])
# OU seja com 0.9 de confiança todos intervalos contém o parametro verdadeiro da mediana

as.numeric(length$IC_mediana$LI[41:80]) < median(abalone$length) &
  mean(abalone$length) < as.numeric(length$IC_mediana$LS[41:80])
# OU seja com 0.95 de confiança todos intervalos contém o parametro verdadeiro da mediana

as.numeric(length$IC_mediana$LI[81:120]) < median(abalone$length) &
  mean(abalone$length) < as.numeric(length$IC_mediana$LS[81:120])
# OU seja com 0.99 de confiança todos intervalos contém o parametro verdadeiro da mediana


# Análise das medidas Dianóstico ------------------------------------------

diameter$estimativas
whole.weight$estimativas
length$estimativas

grafico_diag <- function(var, ylab, ylim=NULL, linha_0 = TRUE, arquivo = 'lixo')
{
  
  medida <- sub('.*_', '',var)
  
  
  jpeg(paste0('img/', arquivo, '.jpg'), 1080, 720, quality = 100,)
  
  if(!is.null(ylim)){
    plot(x= whole.weight$estimativas$n, whole.weight$estimativas[[var]],
         type = 'b', col = '#004586', lty = 2, lwd = 2, ylim = ylim,
         main = paste(ylab, 'Bootstrap para',medida, 'por tamanho de amostra'),
         xlab = 'Tamanho de amostra',
         ylab = ylab)}
  else{
    plot(x= whole.weight$estimativas$n, whole.weight$estimativas[[var]],
         type = 'b', col = '#004586', lty = 2, lwd = 2,
         main = paste(ylab, 'Bootstrap para',medida, 'por tamanho de amostra'),
         xlab = 'Tamanho de amostra',
         ylab = ylab)}
  
    lines(x= diameter$estimativas$n,y=diameter$estimativas[[var]], type = 'b', col = '#FFD320', lty = 2, lwd = 2)
    lines(x= length$estimativas$n,y=length$estimativas[[var]],     type = 'b', col = '#579D1C', lty = 2, lwd = 2)
    if(linha_0 == TRUE) abline(h = 0, lty = 2, lwd = 2, col = 'red')
    legend(x = "topright",          # Position
           legend = c('Whole.Weight',
                      'Length',
                      'Diameter'),  
           lty = 2,           # Line types
           col = c('#004586', '#579D1C', '#FFD320'),           # Line colors
           lwd = 2)  
  

    dev.off()
    if(!is.null(ylim)){
      plot(x= whole.weight$estimativas$n, whole.weight$estimativas[[var]],
           type = 'b', col = '#004586', lty = 2, lwd = 2, ylim = ylim,
           main = paste(ylab, 'Bootstrap para',medida, 'por tamanho de amostra'),
           xlab = 'Tamanho de amostra',
           ylab = ylab)}
    else{
      plot(x= whole.weight$estimativas$n, whole.weight$estimativas[[var]],
           type = 'b', col = '#004586', lty = 2, lwd = 2,
           main = paste(ylab, 'Bootstrap para',medida, 'por tamanho de amostra'),
           xlab = 'Tamanho de amostra',
           ylab = ylab)}
    
    lines(x= diameter$estimativas$n,y=diameter$estimativas[[var]], type = 'b', col = '#FFD320', lty = 2, lwd = 2)
    lines(x= length$estimativas$n,y=length$estimativas[[var]],     type = 'b', col = '#579D1C', lty = 2, lwd = 2)
    if(linha_0 == TRUE) abline(h = 0, lty = 2, lwd = 2, col = 'red')
    legend(x = "topright",          # Position
           legend = c('Whole.Weight',
                      'Length',
                      'Diameter'),  
           lty = 2,           # Line types
           col = c('#004586', '#579D1C', '#FFD320'),           # Line colors
           lwd = 2)  
}

whole.weight$estimativas
grafico_diag(var = 'vicio_media',ylab = 'Vício', arquivo = 'vicio_media')
# O vício para media é baixo, e conforme o n aumenta ele se aproxima de zero
grafico_diag(var = 'vicio_mediana',ylab = 'Vício', ylim = c(-0.0025, 0.0008), arquivo = 'vicio_mediana')
# O vício para mediana é baixo, e conforme o n aumenta ele se aproxima de zero

grafico_diag(var = 'erro_media',ylab = 'Erro', arquivo = 'erro_media', ylim = c(0,0.11))
# Erro para media é baixo, e conforme o n aumenta ele diminui
grafico_diag(var = 'erro_mediana',ylab = 'Erro', ylim = c(0,0.15), arquivo = 'erro_mediana')
# Erro para mediana é baixo, e conforme o n aumenta ele diminui
grafico_diag(var = 'EQM_media',ylab = 'EQM', arquivo = 'EQM_media')
# O EQM da media é minusculo, aproximadamente 0.0000006 para whole.weight em n = 20, conforme n aumenta ele estabiliza em 0
grafico_diag(var = 'EQM_mediana',ylab = 'EQM', ylim = c(0, 8e-6), arquivo = 'EQM_mediana')
# O EQM da mediana é minusculo, aproximadamente 0.000006, mesmo para n =20, conforme n aumenta ele estabiliza em 0

