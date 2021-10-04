################################################################
#   Lab 1 Algunos Métodos de Evaluación con Datos Limitados
#   Método: LBSPR
#   Prof: Luis A. Cubillos
#   Curso: V SIBECORP
#   Fecha: 4 Octubre 2021
#==============================================================


# Carga paquetes necesarios -----------------------------------------------
#install.packages("LBSPR")
library(LBSPR)

# Datos de frecuencia de tallas -------------------------------------------

# Datos de la pesquería artesanal de reineta Brama australis en Chile

lfd_reineta = read.csv("Datos/LFD_Brama_artesanal_total.csv")
head(lfd_reineta)
summary(lfd_reineta[-1])


# Preparación y estimación del estatus ------------------------------------

MyPars <- new("LB_pars")
MyPars@Species <- "Reineta"
#NOTA:La maxima longitud no puede ser menor que Linf
range(lfd_reineta[,1]) #rango de tallas
MyPars@Linf <- 56.9 # parametros Moyano
# Madurez
MyPars@L50 <- 37.7 #LEAL et al.2017
MyPars@L95 <- 43
#calcular M/K Langostino
MyPars@MK <- 1.94 #IFOP M=0.35, K=0.18
MyPars@L_units <- "cm"
Brama_LFD1 <- new("LB_lengths",
                  LB_pars=MyPars,
                  file="Datos/LFD_Brama_artesanal_total.csv",
                  dataType="freq",header=TRUE)
Brama_M1 <- LBSPRfit(MyPars, Brama_LFD1)

#Ajuste a los datos
plotSize(Brama_M1)
#Selectividad comparada con madurez
plotMat(Brama_M1)
#Diagnóstico
plotEsts(Brama_M1)

Brama_M1@Ests
Brama_M1@YPR
Brama_M1@Vars

MyPars@SPR <- 0.4
yr <- 14
MyPars@SL50 <- Brama_M1@SL50[yr]
MyPars@SL95 <- Brama_M1@SL95[yr]
plotTarg(MyPars, Brama_LFD1, yr=yr)


