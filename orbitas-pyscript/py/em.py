import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from scipy.integrate import quad
import math
#from PIL import Image
import matplotlib.animation as animation
from matplotlib.patches import Circle

def gerarOrbita1():
    
    from js import x, v

    x0 = float(x)
    v0 = float(v)

    def v(u, l):
        v = -u + (l ** 2) * (u ** 2) / 2 - (l ** 2) * (u ** 3)
        return v

    rs_sun = 3  # Hipótese sobre o raio de Schwarzschild do corpo central (em km). Ligeiramente maior que o do Sol.

    rst = 2.0 * (x0 / rs_sun)
    ust = 1 / rst
    l = rst * v0
    E = v(ust, l) + 0.000000001  # Hipótese: dr/dt = 0 inicialmente.
    norbit = 10

    npoints = 500

    if l > math.sqrt(12):
        coef = [- 3 * (l ** 2), l ** 2, -1]
        a = np.roots(coef)
        umax = sorted(a)[1].real
        umin = sorted(a)[0].real
        vmin = v(umin, l)
        vmax = v(umax, l)

    coef = [- (l ** 2), (l ** 2) / 2, -1, -E]
    roots = np.roots(coef)
    tp1 = roots[2]
    tp2 = roots[1]
    tp3 = roots[0]

    eps = 0.00000001
    if l > math.sqrt(12):
        if E < 0 and ust < tp2.real:
            u1 = tp1.real * (1 + eps)
            u2 = tp2.real * (1 - eps)
        elif 0 < E < vmax and ust < tp2.real:
            u1 = ust / 20
            u2 = ust
            norbit = 0.5
        elif E < vmax and ust > tp3.real:
            u1 = 0.5
            u2 = tp3.real * (1 + eps)
            norbit = 0.5
        elif E > vmax:
            u1 = ust
            u2 = 0.5
            norbit = 0.5
    else:
        if E >= 0:
            u1 = ust
            u2 = 0.5
            norbit = 0.5
        else:
            u1 = tp1.real * (1 + eps)  # COM TP3 NÃO GERA ORBITAS COM X0 PEQUENOS
            u2 = 0.5
            norbit = 0.5

    w = sp.Symbol('w')


    def tau_integrand(w):
        tau_integrand = w ** (-2) * (2.0 * (E - v(w, l))) ** (-1 / 2)
        return tau_integrand


    Ttotal, erroT = quad(tau_integrand, u1, u2)  # Computes total time to go from u1 to u2.
    dt = Ttotal / npoints  # Sets the time step as 1/100 of the total time.
    ud = [u2]
    for i in range(npoints + 20):
        ud.append(ud[i] - dt * ud[i] ** 2 * (2.0 * (E - v(ud[i], l))) ** (1 / 2.0))
        if ud[-1].imag != 0 or math.isnan(ud[-1]):
            ud = ud[:-2]
            break
    uc = ud[::-1]
    n = len(uc)


    def theta(w):
        theta = l * (2.0 * (E - v(w, l))) ** (-1 / 2)
        return theta


    delphi, erro = quad(theta, u1, u2)

    if abs(u1 - ust) < abs(u2 - ust):
        phi1 = []
        for i in range(len(uc)):
            a = quad(theta, u1, uc[i])
            phi1.append(abs(a[0]))

        phi2 = []
        for j in range(len(ud)):
            b = quad(theta, u2, ud[j])
            phi2.append(abs(b[0]))

        if norbit == 0.5:
            utotal = uc
        else:
            utotal = np.concatenate([uc, ud] * (norbit))
    else:
        phi2 = []
        for i in range(len(uc)):
            a = quad(theta, u1, uc[i])
            phi2.append(abs(a[0]))

        phi1 = []
        for j in range(len(ud)):
            b = quad(theta, u2, ud[j])
            phi1.append(abs(b[0]))

        if norbit == 0.5:
            utotal = ud
        else:
            utotal = np.concatenate([ud, uc] * (norbit))

    accphi = [0] * (len(utotal))

    if norbit == 0.5:
        accphi = phi1
        x = [0] * (len(utotal))
        y = [0] * (len(utotal))
        for i in range(len(utotal)):
            x[i] = (math.cos(accphi[i])) / utotal[i] * (rs_sun / 2.0)
            y[i] = (math.sin(accphi[i])) / utotal[i] * (rs_sun / 2.0)
    else:
        for i in range(norbit):
            for j in range(n):
                accphi[j + (2 * i * n)] = 2 * i * delphi + phi1[j]
                accphi[j + ((2 * i + 1) * n)] = ((2 * i) + 1) * delphi + phi2[j]
        x = [0] * (2 * norbit * n)
        y = [0] * (2 * norbit * n)
        for i in range(2 * norbit * n):
            x[i] = (math.cos(accphi[i])) / utotal[i] * (rs_sun / 2.0)
            y[i] = (math.sin(accphi[i])) / utotal[i] * (rs_sun / 2.0)

    fig = plt.figure()

    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.gca().set_aspect('equal')
    ax = plt.gca()
    ax.spines['bottom'].set_color('white')
    ax.tick_params(axis='x', colors='white')
    ax.tick_params(axis='y', colors='white')
    ax.spines['top'].set_color('white')
    ax.spines['right'].set_color('white')
    ax.spines['left'].set_color('white')
    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')
    fig.patch.set_facecolor('#0E1117')
    ax.set_facecolor("black")
    circle = Circle((0, 0), rs_sun, color='dimgrey', linewidth=0)
    plt.gca().add_patch(circle)
    if x0 <= 3.8:
        plt.axis([x0 - 2, x0 + 2, -2, 2])
    elif x0 >= 15:
        plt.axis([x0 - 100, x0 + 100, -100, 100])
    else:
        plt.axis([- (rs_sun / 2.0) * x0, (rs_sun / 2.0) * x0, - (rs_sun / 2.0) * x0, (rs_sun / 2.0) * x0])

    # Montagem do gif

    graph, = plt.plot([], [], color="gold", markersize=3, label='Tempo: 0 s')
    L = plt.legend(loc=1)
    # o_plot = st.pyplot(plt)

    plt.close()  # Não mostra a imagem de fundo


    def animate(i):
        lab = 'Tempo: ' + str(round(dt * i * (rs_sun / 2.0) * 3e-5,-int(math.floor(math.log10(abs(dt * (rs_sun / 2.0) * 3e-5)))))) + ' s'
        graph.set_data(x[:i], y[:i])
        L.get_texts()[0].set_text(lab)  # Atualiza a legenda a cada frame
        return graph,


    skipframes = int(len(x) / 200)
    if skipframes == 0:
        skipframes = 1

    # for i in range(len[x]):
    # animate(i)
    # time.sleep(0.01)

    ani1 = animation.FuncAnimation(fig, animate, frames=range(0, len(x), skipframes), interval=30, blit=True,repeat=False)
    display(ani1, target="graph", append=True)

    # HTML(ani1.to_jshtml())
    # components.html(ani1.to_jshtml(),height=800)

    # st.pyplot(fig)
    # width = st.sidebar.slider("plot width", 1, 25, 3)