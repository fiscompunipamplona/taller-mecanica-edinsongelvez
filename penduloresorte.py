##Péndulo resorte

import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
import matplotlib.animation as anim

# Parámetros del sistema
# Longitud del resorte en reposo (m)
l = 1.0
# Constante de resorte (N / m)
k = 15.0
# La masa unida al resorte (kg)
m = 1.0
# # Aceleración gravitacional
g = 9.81

# Valores iniciales
t0 = 0.0
theta0 = np.pi / 3.0
omega0 = 0.0
x0 = 0.25
v0 = 0.0

# ODE ajustes
dt = 1.0 / 30.0
t_max = 60.0

# Los límites de la cuadrícula
plt.xlim(-2.0 * l, 2.0 * l)
plt.ylim(-3.0 * l, 0.0)
plt.gca().set_aspect('equal')
plt.gca().grid()

#origin = plt.plot([ 0.0 ], [ 0.0 ], 'ro')
spring, = plt.plot([], [], 'b-', linewidth = 2)
mass, = plt.plot([], [], 'ro', markersize = 10)

frames = []

# Esta función dibuja el sistema dados los valores de la generalizada.
# coordenadas.
def draw(n):
  spr_x = (l + frames[n][2]) * np.sin(frames[n][0])
  spr_y = -1.0 * (l + frames[n][2]) * np.cos(frames[n][0])
  spring.set_data([ 0.0, spr_x ], [ 0.0, spr_y ])
  mass.set_data([ spr_x ], [ spr_y ])
  return spring, mass

def deriv(t, y):
  dth = y[1]
  dx = y[3]

  #  Aceleración angular
  d2th = -1.0 * (g * np.sin(y[0]) + 2.0 * dx * dth) / (l + y[2])
  # Aceleración de la masa respecto al muelle.
  d2x = (l + y[2]) * dth ** 2 + g * np.cos(y[0]) - k / m * y[2]
  
  return [ dth, d2th, dx, d2x ]

solver = ode(deriv).set_integrator('vode', method = 'bdf')
solver.set_initial_value([ theta0, omega0, x0, v0 ], t0)

while solver.successful() and solver.t < t_max:
  solver.integrate(solver.t + dt)
  frames.append(solver.y)

ani = anim.FuncAnimation(plt.gcf(), draw, np.arange(0, len(frames)),
  interval = 1e3 * dt)

#ani.save('spring_pendulum.mp4', fps = 30)
plt.show()

# Se observa el movimiento de la partícula, y su trayectoria, se sitúa en la posición del extremo libre del muelle
# sin deformar, L0=1 m por debajo del punto de sujeción del muelle, mientras que en los fundamentos físicos, se ha
# situado en la posición de equilibrio de la partícula, ye=mg/k por debajo del origen, o Le=L0+mg/k por debajo del
# punto de sujeción.

# En este caso se resuelve el sistema de dos ecuaciones diferencial de segundo orden
#por el método de Runge-Kutta. Este procedimiento numérico produce, en general, buenos resultados,
#que podemos confirmar observando que la energía se mantiene constante.

