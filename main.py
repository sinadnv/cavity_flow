## Init
import numpy as np
import matplotlib.pyplot as plt

Lx = 2
Ly = 2
nx = 41
ny = 41
nt = 500
nit = 50
c = 1
dx = Lx / (nx - 1)
dy = Ly / (ny - 1)
x = np.linspace(0, Lx, nx)
y = np.linspace(0, Ly, ny)
X, Y = np.meshgrid(x, y)

rho = 1
nu = .1
dt = .001

u = np.zeros((ny, nx))
v = np.zeros((ny, nx))
p = np.zeros((ny, nx))
b = np.zeros((ny, nx))


def pCalc(u,v,p, dx, dy):
    error = 1
    counter = 0
    while error > 1e-6:
        error = 0
        pn = p.copy()
        for i in range(1,ny-1):
            for j in range(1,nx-1):
                p[i, j] = ((pn[i + 1, j] + pn[i - 1, j]) * (dy * dy) + (pn[i, j + 1] + p[i, j - 1] * (dx * dx))) /\
                          (2 * (dx * dx + dy * dy)) \
                          + (rho * dx * dx * dy * dy) / (2 * (dx * dx + dy * dy)) * (
                                  ((u[i + 1, j] - u[i - 1, j]) / (2 * dx)) ** 2
                                  + 2 * ((u[i, j + 1] - u[i, j - 1]) / (2 * dy)) * ((v[i + 1, j] - v[i - 1, j]) / (2 * dx))
                                  + ((v[i, j + 1] - v[i, j - 1]) / (2 * dy)) ** 2)
                error = error + np.abs(pn[i,j]-p[i,j])
        p[:, -1] = p[:, -2]
        p[0, :] = p[1, :]
        p[:, 0] = p[:, 1]
        p[-1, :] = 0
        counter = counter + 1
        print(error)
    return p


##
def uvCalc(u,v,p,dx,dy,dt):
    un = u.copy()
    vn = v.copy()
    for i in range(1,ny-1):
        for j in range(1,nx-1):
            u[i, j] = un[i, j] \
                      - un[i, j] * (dt / dx) * (un[i, j] - un[i - 1, j]) \
                      - vn[i, j] * (dt / dy) * (un[i, j] - un[i, j - 1]) \
                      - (dt / (rho * 2 * dx)) * (p[i + 1, j] - p[i - 1, j]) \
                      + (nu * dt / (dx * dx)) * (un[i + 1, j] - 2 * un[i, j] + un[i - 1, j]) \
                      + (nu * dt / (dy * dy)) * (un[i, j + 1] - 2 * un[i, j] + un[i, j - 1])

            v[i, j] = vn[i, j] \
                      - un[i, j] * (dt/dx) * (vn[i, j] - vn[i - 1, j]) \
                      - vn[i, j] * (dt / dy) * (un[i, j] - un[i, j - 1]) \
                      - (dt / (rho * 2 * dy)) * (p[i, j + 1] - p[i, j - 1]) \
                      + (nu * dt / (dx * dx)) * (vn[i + 1, j] - 2 * vn[i, j] + vn[i - 1, j]) \
                      + (nu * dt / (dy * dy)) * (vn[i, j + 1] - 2 * vn[i, j] + vn[i, j - 1])

    u[0, :]  = 0
    u[:, 0]  = 0
    u[:, -1] = 0
    u[-1, :] = 1
    v[0, :]  = 0
    v[-1, :] = 0
    v[:, 0]  = 0
    v[:, -1] = 0


    return u, v


for t in range(2):
    print('Printing results at t = ', t)
    p = pCalc(u,v,p, dx, dy)
    u,v = uvCalc(u,v,p,dx,dy,dt)


plt.contourf(X,Y, p)
plt.show()