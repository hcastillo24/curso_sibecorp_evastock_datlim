################################################################
#   Lab 1 Algunos Métodos de Evaluación con Datos Limitados
#   Método: CMSY2, OCOM, zBRT
#   Prof: Luis A. Cubillos
#   Curso: V SIBECORP
#   Fecha: 4 Octubre 2021
#==============================================================


# Carga paquetes necesarios -----------------------------------------------
#devtools::install_github("cfree14/datalimited2")
library(datalimited2)

# Datos: serie de tiempo de capturas --------------------------------------
reineta <- structure(list(
  yr = 1994:2020,
  ct = c(1186, 3930, 5585, 5998, 6332, 6828, 8159, 15156, 4429, 2645, 3764, 12707, 2517, 3743, 6160, 15199, 16977, 28814, 23079, 11955, 35975, 34218, 27586, 25267, 28175, 44288, 38109)),
  .Names = c("yr", "ct"),
  class = "data.frame",
  row.names = 1:27)

# Prepara gráfico clásico -------------------------------------------------

yt <- reineta$ct
range(reineta$yr)
captura_reineta <- ts(yt,start = 1994,frequency = 1)
plot(captura_reineta,ty="b",pch=19,ylim=c(0,max(reineta$ct)*1.2),las=1,ylab="Captura (toneledas)",xlab="Años")


# CMSY2 - Froese et al.  2017 ---------------------------------------------

catch = reineta$ct
yr = reineta$yr
OCm1 = cmsy2(year=yr, catch=catch, resilience="High")
str(OCm1)
old.op <- par()
plot_dlm(OCm1)

par(old.op)
m1.Bt = OCm1$ref_ts$b
m1.Bt_l = OCm1$ref_ts$b_lo   
m1.Bt_h = OCm1$ref_ts$b_hi
plot(yr, m1.Bt, type="b", bty="n", pch=19,xlab="",ylab="Biomasa",ylim=c(0,max(m1.Bt_h)))
lines(yr, m1.Bt_l, lty=2)
lines(yr, m1.Bt_h, lty=2)

CMSY1_ref_pts = OCm1[["ref_pts"]]
CMSY1_ref_pts


# OCOM Zhou et al. (2017a) ------------------------------------------------

OCm2 = ocom(year = yr, catch = catch, m = 0.35)
plot_dlm(OCm2)

par(old.op)
m2.Bt = OCm2$ref_ts$b
m2.Bt_l = OCm2$ref_ts$b_lo   
m2.Bt_h = OCm2$ref_ts$b_hi
plot(yr, m2.Bt, type="b", bty="n", pch=19,xlab="",ylab="Biomasa",ylim=c(0,max(m2.Bt_h)))
lines(yr, m2.Bt_l, lty=2)
lines(yr, m2.Bt_h, lty=2)

OCOM_ref_pts = OCm2[["ref_pts"]]
OCOM_ref_pts


# zBRT de Zhou et al. (2017b) ---------------------------------------------

OCm3 <- zbrt(year=yr, catch=catch)
plot_dlm(OCm3)





