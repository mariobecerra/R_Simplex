library(dplyr)

# Mario Becerra Contreras
# 2016-08-31

simplex <- function(c, A, b, verbose = 0) {
  # Salida: Lista con las siguientes entradas 
  #   sol: valor  optimo del problema, 
  #   fx: valor de la funcion objetivo en x, 
  #   flag: indica los siguientes casos: 
  #     flag == 0 funcion objetivo no acotada superiormente, 
  #     flag == 1 se encontró solucion óptima,
  #     flag == 2 conjunto factible vacío.
  #   verbose: El nivel de verbosidad de impresión (0 o 1)
  
  # Para probar
  # simplex(c(3,2), matrix(c(-1, 1, 2, 3, 1, -1), ncol = 2), c(12, 8, 10))
  # simplex(c(1,1), matrix(c(4, -1, 2, 1, -5, 2), ncol = 2, byrow = T), c(8, 10, 2))
  
  m <- dim(A)[1]
  n <- dim(A)[2]
  
  
  A <- cbind(A, matrix(rep(0, m*m), ncol = m))
  A <- t(cbind(matrix(rep(0, n*(n + m)), ncol = n), t(A)))
  
  B=(n+1):(n+m); #Básicas
  N=1:n; #No básicas
  
  if(sum( b < 0 ) > 0) {
    flag <- -1
    cat('Conjunto factible vacío')
  } else  {
    flag <- 2
  }
  
  
  fx <- 0
  b <- c(rep(0, n), b)
  c <- c(c, rep(0, m))
  iter <- 0
  r <- expand.grid(1:(m+n), 1:(m+n)) %>% arrange(Var1) %>% as.matrix()
  while(flag == 2 & iter < 5000) {
    iter = iter + 1
    
    if(verbose) cat("\n\n\n", iter, "\n\n")
    if(verbose) cat("A\n")
    if(verbose) print(A[B, N])
    if(verbose) cat("Ax: ", (A %*% b), "\n")
    if(verbose) cat("b:", b[B], "\n")
    if(verbose) cat("c: ", c[N], "\n")
    if(verbose) cat("B: ", B, "\n")
    if(verbose) cat("N: ", N, "\n")
    if(verbose) cat("fx: ", fx, "\n")
    
    if (sum(c > 0) > 0) {
      k <- which(c > 0)[1] # columna pivote
      temp <- rep(Inf, m+n) 
      temp[A[,k] > 0] <- 0
      temp <- temp + b/A[,k]
      l <- which(temp == min(temp, na.rm = T))[1]
      if(verbose) cat("l: ", l, "\n")
      if(verbose) cat("k: ", k, "\n")
      if(verbose) cat("A[l,k]: ", A[l,k], "\n")
      
      if(sum(l == 0, na.rm = T) | sum(temp, na.rm = T)==0 | sum(temp == Inf, na.rm = T) == m ) {
        flag = 0;
        print('El problema no es acotado');
      } else {
        fx <- fx + c[k]*b[l]/A[l,k] #Actualizar func obj
        temp1 <- N[N != k] # Conjunto de las no básicas menos k
        temp2 <- B[B != l] # Conjunto de las básicas menos l
        
        c[temp1] <- c[temp1] - c[k]*t(A[l,temp1])/A[l,k] #Act cj para todo j No básica menos k
        c[l] <- -c[k]/A[l,k] #Actualizar ck
        
        # Para las básicas menos l
        b[temp2] <- b[temp2] - A[temp2, k]*b[l]/A[l,k] #Actualizar bi pa todo i Básica menos l
        A[temp2, temp1] <- A[temp2, temp1] - (1/A[l,k]) * A[temp2,k] %*% t(A[l,temp1]) #Actualizar aij pa todo j No básico menos k y todo i Básica menos l
        A[temp2, l] = - A[temp2, k] / A[l,k] #Actualizar ail para todo i Básica menos l
        
        b[k] = b[l]/A[l,k] #Actualizar bk
        
        # Para las no básicas menos k
        A[k,temp1] = A[l,temp1] / A[l,k] #Actualizar akj
        
        A[k,l] = 1/(A[l,k]) #Actualizar akl
        
        # Actualizar variables básicas
        B <- sort(c(temp2, k)) # Nuevo conjunto de variables básicas 
        N <- sort(c(temp1, l)) # Nuevo conjunto de variables no básicas
        c[B] <- 0
        b[N] <- 0
        temp3 <- expand.grid(N, B) %>% arrange(Var1) %>% as.matrix()
        
        # Estúpidamente ineficiente, pero no hubo tiempo de optimizar esto
        idx <- rep(F, nrow(r))
        for(i in seq_along(idx)) {
          sumas <- rep(0, nrow(temp3))
          for(j in seq_along(sumas)){
            aaa <- r[i,] == temp3[j,]
            sumas[j] <- sum(aaa)
          }
          
          idx_loop <- which(sumas == 2)
          if(length(idx_loop) > 0) idx[i] <- T
        }
        
        q <- r[!idx,]
        
        for(i in 1:nrow(q)) {
          A[q[i,2], q[i,1]] <- 0;
        } # end for
        
      } # end if
      
    } else {
      flag <- 1
      print('Solución óptima encontrada')
    }
    
  } # end while
  x <- b[1:n];
  return(list(sol = x, fx = fx, flag = flag))
}