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

b[[2]

# Função pika -------------------------------------------------------------

# calcular nessa função também os ics e retornar uma lista
Bootstrap <- function(N,enes,var,banco= abalone, conf.level=0.95)
{
  est_bootstrap <- data.frame('media' = rep(NA, 10), 'mediana' = NA, 'erro_media' = NA, 'erro_mediana' = NA, 'vicio_media' = NA, 'vicio_mediana' = NA, 'EQM_media' = NA, 'EQM_mediana' = NA,  'n' = NA)
  
  media_global   <- mean(banco[[var]])
  mediana_global <- median(banco[[var]])
  
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
    
    # Intervalos de confiança
    quantils <- quantile(estatistiscas, c(1-conf.level, conf.level))
    
    
    
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
  
  return(est_bootstrap)
  
}

dados <- Bootstrap(100, seq(20,200,20), 'diameter')
dados


# Intervalo de confiança --------------------------------------------------

# Determinar os intervalos de 90%, 95% e 99% via bootstrap usando os tipos "normal", "percentil", "t de Student" e "BCA". Todos devem ser implementados na mão, sem uso de pacotes.
help("t.test")
IC <- function(dados, var, IC = c('normal', 'percentil', 't', 'BCA'), conf.level= 0.95)
{
  
  if(var == 'media')
  {
    variavel <- dados['media']
    erro     <- dados['erro_media']
    vicio    <- dados['vicio_media']
    eqm      <- dados['EQM_media']
  }
  else if(var == 'mediana')
  {
    variavel <- dados['mediana']
    erro     <- dados['erro_mediana']
    vicio    <- dados['vicio_mediana']
    eqm      <- dados['EQM_mediana']
  }
  
  if(IC == 'normal')
  {
   IC <- c(variavel - qnorm(conf.level)*erro, variavel + qnorm(conf.level)*erro) 
  }  
  
  
}
