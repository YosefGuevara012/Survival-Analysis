
# Lectura de los datos de cancer de ovarios

can_ovario <- read.csv("stress.csv", header=T, sep=";")
attach(can_ovario)

# Cargar los paquetes survival y muhaz

require(survival)
require(muhaz)

# Estimación de la funcion de riesgo


riesgo.can <- muhaz(Tiempo, Censura)
plot(riesgo.can)
