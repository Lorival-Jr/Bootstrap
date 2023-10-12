library(microbenchmark)


microbenchmark('laço' = for(i in 1:length(n_amostra))
{
  n <- n_amostra[i]
  vetor_medias   <- c(rep(NA,N))
  vetor_medianas <- c(rep(NA,N))
  
  for(j in 1:N){
    
    media   <- mean(sample(abalone$diameter, n, T))
    mediana <- median(sample(abalone$diameter, n, T))
    vetor_medias[j]   <- media
    vetor_medianas[j] <- mediana
    
    print(j)
  }
  
  est_bootstrap[i,] <- c(mean(vetor_medias), median(vetor_medianas), sd(vetor_medias), sd(vetor_medianas))
  
  jpeg(filename = paste0('img/tend_central_', n, '.jpg' ), width = 1080, height = 720, quality = 100)
  par(mfrow = c(1,2))
  hist(vetor_medias, xlab = paste('Valores da média com n=', n), 
       ylab = 'Frequência', main = '')
  abline(v = mean(media), col = 'red', lty = 2)
  hist(vetor_medianas, xlab = paste('Valores da mediana com n=', n), 
       ylab = 'Frequência', main = '')
  abline(v = median(mediana), col = 'blue', lty = 2)
  
  legend('bottom',
         legend = c('Estimativa da Média (Bootstrap)', 'Estimativa da Mediana (Bootstrap)'),
         col = c('red','blue'),
         lty = 2,
         bty = 'n')
  
  par(mfrow = c(1,1))
  
  dev.off()
},
'apply' = Bootstrap(10000, seq(20,200,20), 'diameter'), times = 10)


install.packages('tictoc')
library(tictoc)




Bootstrap <- function(N,enes,var,banco= abalone, conf.level=0.95, seed=NULL)
{
  est_bootstrap <- data.frame('media' = rep(NA, 10), 'mediana' = NA, 'erro_media' = NA, 'erro_mediana' = NA, 'vicio_media' = NA, 'vicio_mediana' = NA, 'EQM_media' = NA, 'EQM_mediana' = NA,  'n' = NA)
  
  ICs_media   <-  data.frame('Tipo_IC' = rep(NA,40*length(conf.level)), 'LI'= NA, 'LS' = NA, n = NA)
  ICs_mediana <-  data.frame('Tipo_IC' = rep(NA,40*length(conf.level)), 'LI'= NA, 'LS' = NA, n = NA)
  
  tempo_IC_norm       <- list('tic'=0, 'toc'=0)
  tempo_IC_perc       <- list('tic'=0, 'toc'=0)
  tempo_IC_t       <- list('tic'=0, 'toc'=0)
  tempo_IC_bca       <- list('tic'=0, 'toc'=0)
  tempo_graficos  <- list('tic'=0, 'toc'=0)
  tempo_bootstrap <- list('tic'=0, 'toc'=0)
  
  
  media_global   <- mean(banco[[var]])
  mediana_global <- median(banco[[var]])
  
  if(!is.null(seed)){set.seed(seed)}
  
  tic()
  
  amostras <- lapply(enes, banco = banco, var = var, N = N,
                     FUN = function(x, banco, var,N){
                       index <- round(runif(x*N, 0.5, 0.4 + length(banco[[var]])))
                       amostra <- banco[[var]][index]
                       return(amostra)
                     }
  )
  T_bootstrap <- toc() 
  tempo_bootstrap$tic <- append(T_bootstrap$tic, tempo_bootstrap$tic)
  tempo_bootstrap$toc <- append(T_bootstrap$toc, tempo_bootstrap$toc)
  for(i in enes)
  {
    tic()
    estatistiscas <- lapply(split(amostras[[i/enes[1]]], ceiling(seq_along(amostras[[i/enes[1]]])/i)), 
                            FUN = function(x)
                            {return(c(mean(x), median(x)))})
    
    estatistiscas <- do.call(rbind, estatistiscas)
    
    T_bootstrap <- toc() 
    tempo_bootstrap$tic <- append(T_bootstrap$tic, tempo_bootstrap$tic)
    tempo_bootstrap$toc <- append(T_bootstrap$toc, tempo_bootstrap$toc)
    
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
      
      tic()
      IC_norm_mean   <- c(valores[1] - qnorm(conf)*valores[3], valores[1] + qnorm(conf)*valores[3])
      IC_norm_median <- c(valores[2] - qnorm(conf)*valores[4], valores[2] + qnorm(conf)*valores[4])
      T_norm <- toc() 
      
      tempo_IC_norm$tic <- append(T_norm$tic, tempo_IC_norm$tic)
      tempo_IC_norm$toc <- append(T_norm$toc, tempo_IC_norm$toc)
      
      # Percentil Reverso (básico)
      tic()
      Q1 <- quantile(estatistiscas[,1], c(1-conf, conf))
      Q2 <- quantile(estatistiscas[,2], c(1-conf, conf))
      
      IC_perc_mean <-   c(2*valores[1] - Q1[[2]], 2*valores[1] - Q1[[1]])
      IC_perc_median <- c(2*valores[2] - Q2[[2]], 2*valores[2] - Q2[[1]])
      
      T_perc <- toc() 
      
      tempo_IC_perc$tic <- append(T_norm$tic, tempo_IC_perc$tic)
      tempo_IC_perc$toc <- append(T_norm$toc, tempo_IC_perc$toc)
      
      # t-Student
      
      tic()
      IC_t_mean   <- c(valores[1] - qt(conf, i -1 )*valores[3],
                       valores[1] + qt(conf, i -1)*valores[3])
      IC_t_median <- c(valores[2] - qt(conf, i -1)*valores[4],
                       valores[2] + qt(conf, i -1)*valores[4])   
      T_t <- toc() 
      
      tempo_IC_t$tic <- append(T_norm$tic, tempo_IC_t$tic)
      tempo_IC_t$toc <- append(T_norm$toc, tempo_IC_t$toc)
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
    
    tic()
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
    
    T_graf <- toc() 
    
    tempo_graficos$tic <- append(as.numeric(T_graf$tic), tempo_graficos$tic)
    tempo_graficos$toc <- append(as.numeric(T_graf$toc), tempo_graficos$toc)
    
    cat('Repetição:',i/enes[1], '\n')
    
    
  }
  
  output <- list('estimativas' = est_bootstrap, 'IC_media' = ICs_media, 'IC_mediana' = ICs_mediana)
  output <- list('bootstrap'= tempo_bootstrap, 'graficos' = tempo_graficos, 'IC_norm' = tempo_IC_norm,
                 'IC_t' = tempo_IC_t, 'IC_perc' = tempo_IC_perc )
  
  
  return(output)
  
}
tic()
tempos <- Bootstrap(700000, seq(20,200,20), 'diameter', conf.level =  c(0.9, 0.95, 0.99))
tempo_tudo <- toc()

tempo_tudo <- as.numeric(tempo_tudo$toc - tempo_tudo$tic)
tempo_graficos <- sum(tempos$graficos$toc - tempos$graficos$tic)
tempo_bootstrap <- sum(tempos$bootstrap$toc - tempos$bootstrap$tic)
tempo_IC_norm <- sum(tempos$IC_norm$toc - tempos$IC_norm$tic)
tempo_IC_t <- sum(tempos$IC_t$toc - tempos$IC_t$tic)
tempo_IC_perc <- sum(tempos$IC_perc$toc - tempos$IC_perc$tic)

cat(' Tempo total do processo:',  tempo_tudo,'\n',
    'Tempo do método Bootstrap:', tempo_bootstrap, '\n',
    'Tempo do IC Normal:',        tempo_IC_norm, '\n',
    'Tempo do IC Percentil:',     tempo_IC_perc, '\n',
    'Tempo do IC T Student:',     tempo_IC_t, '\n',
    'Tempo do IC BCA: NA \n',
    'Tempo dos gráficos:', tempo_graficos, '\n')
