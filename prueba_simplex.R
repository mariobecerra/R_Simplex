source("simplex.R")

# Prueba Simplex

# Problema 1
# fx = 260.7

A1 = matrix(c(2, 3, 4, 1, 2, 9), ncol = 2, byrow = T)
c1 = c(21, 31)
b1 = c(25, 32, 54)

simplex(c1, A1, b1)



# Problema 2
# fx = 290.18

A2 = matrix(c(4, 5, 3, 7, 2, 9, 9, 11, 7), ncol = 3, byrow = T)
b2 = c(19, 25, 36)
c2 = c(51, 62, 76)

simplex(c2, A2, b2)



# Problema 3
# fx = 8

A2 = matrix(c(4, -1, 2, 1, -5, 2), ncol = 2, byrow = T)
b2 = c(8, 10, 2)
c2 = c(1, 1)

simplex(c2, A2, b2)



# Problema 4
# fx = 22

A4 = matrix(c(-1, 3, 1, 1, 2, -1), ncol = 2, byrow = T)
b4 = c(12, 8, 10)
c4 = c(3, 2)

simplex(c4, A4, b4)

