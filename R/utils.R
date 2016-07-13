belong.to.disk  <- function(x, center, halfDiameter.power2) t(dif <- x - center) %*% dif  < halfDiameter.power2
