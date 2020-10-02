# Here we are going to use the package coupled_gauss_logistic.R
source("coupled_gauss_logistic.R")

library(ggplot2)
library(plotly)

# Figuras do Anderson 
out <- parameterSpace(delta=1.8, 
          gama=seq(0,7,length=1000), 
          beta=seq(0,1,length=1000),
          x0=0.2, transient=300, maxiter=300)
plt <- plotParSpace(out, zmin=-7, zmax=+7)
plt

a <- highlightExtreme(delta=1.8, gama=seq(0,7,length=100), beta=seq(0,1,length=100),
                 from="xm",to="xr",maxiter=1)
a
ggplotly(a)

# Extreme curves
e1mr <- extreme(delta=c(1.8,1.8),gama=c(0,7),beta=c(0,1),
                p1=c(1.8,7,0.25), p2=c(1.8,7,0.75),direction="left",
                from="xm", to="xr", maxiter=1, step=0.001)

e2rm <- extreme(delta=c(1.8,1.8),gama=c(0,7),beta=c(0,1),
          p1=c(1.8,7,0.01), p2=c(1.8,7,0.25),direction="left",
          from="xr", to="xm", maxiter=2, step=0.001)

e2mr_1 <- extreme(delta=c(1.8,1.8),gama=c(0,7),beta=c(0,1),
                p1=c(1.8,7,0.25), p2=c(1.8,7,0.5),direction="left",
                from="xm", to="xr", maxiter=2, step=0.001)
e2mr_2 <- extreme(delta=c(1.8,1.8),gama=c(0,7),beta=c(0,1),
                  p1=c(1.8,7,0.5), p2=c(1.8,7,0.75),direction="left",
                  from="xm", to="xr", maxiter=2, step=0.001)

plt + geom_path(data=e1mr, aes(gama,beta), size=0.001, colour="cyan") +
      geom_path(data=e2rm, aes(gama,beta), size=0.001, colour="black") +
      geom_path(data=e2mr_1, aes(gama,beta), size=0.001, colour="blue") +
      geom_path(data=e2mr_2, aes(gama,beta), size=0.001, colour="blue")

save.image(file="testes_anderson.RData")
