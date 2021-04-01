
# Estimación de la función S

# Datos de hepatitis

# libreria necesaria

require(survival)

# Grupo Control

tpo1 <- c(1,2,3,3,3,5,5,16,16,16,16,16,16,16,16) # Datos hepatitis pag 65
cens1 <- c(0,0,1,1,0,0,0,0,0,0,0,0,0,0,0) # Censurados 0, Observado 1

# Grupo Esteroide

tpo2 <- c(1,1,1,1,4,5,7,8,10,10,12,16,16,16)
cens2 <- c(1,1,1,0,0,1,1,1,1,0,0,0,0,0)

# Datos Total

tpo <- c(1,2,3,3,3,5,5,16,16,16,16,16,16,16,16,1,1,1,1,4,5,7,8,10,10,12,16,16,16)
cens <- c(0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,1,1,1,1,0,0,0,0,0)
grupo <- c(rep(1,15),rep(2,14)) # Crea un vector para diferenciar los grupos
                                # Se es de grupo control o de grupo esteroides

# Estimador K-M Grupo Control

ekm1 <- survfit(Surv(tpo1,cens1)~1)
summary(ekm1)
plot(ekm1,xlab="Tiempo (semanas)", ylab= "S(t) estimada")

# Estimador K-M Grupo Esteroide

ekm2 <- survfit(Surv(tpo2,cens2)~1)
summary(ekm2)
plot(ekm2,xlab="Tiempo (semanas)", ylab= "S(t) estimada")


# Estimador K-M Grupos

ekm <- survfit(Surv(tpo,cens)~grupo)
summary(ekm)
plot(ekm,lty=c(2,1),xlab="Tiempo (semanas)", ylab= "S(t) estimada")
legend(1,0.3,lty=c(2,1), c("Control", "Esteroide"), lwd=1, bty="n")

# Estimador NA Grupo de Control, Fleming and Harrington = fh


ena1 <- survfit(Surv(tpo1,cens1)~1, type = "fh")
summary(ena1)
plot(ena1,xlab="Tiempo (semanas)", ylab= "S(t) estimada")



# Estimador NA Grupo de Eseroide, Fleming and Harrington = fh


ena2 <- survfit(Surv(tpo2,cens2)~1, type = "fh")
summary(ena2)
plot(ena2,xlab="Tiempo (semanas)", ylab= "S(t) estimada")


 # Estimador NA Grupos, Fleming and Harrington = fh


ena <- survfit(Surv(tpo,cens)~grupo, type = "fh")
summary(ena)
plot(ena,xlab="Tiempo (semanas)", ylab= "S(t) estimada")




























