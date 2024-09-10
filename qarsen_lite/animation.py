import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt


def animate_1D_dynamics(x_grid: np.ndarray,
                        dyn: list,
                        t_grid: list,
                        V_array: list,
                        magnitude=True,
                        real=False,
                        imaginary=False,
                        V_ylim=(0,1),
                        psi_ylim=(-1,1),
                        frame_delay=200,
                        **options):
    """Generate an animation of the wavefunction evolving with time overlaped with the potential.
    
    Mandatory arguments
    x: 1D array of equally spaced points in space.
    dyn: 2D array with wavefunction values at several instants in time.
    dt: time step in atomic units.
    V_array: 1D array with the values of the potential a the points x.
    
    Optional keyword arguments
    magnitude = bool: show magnitude of the wavefunction (default:True).
    real = bool: show real part of the wavefunction (default=True).
    imaginary = bool: show imaginary part of the wavefunction (default=False).
    V_ylim = (num,num): set y axis limits for the potential plot.
    psi_ylim = (num,num): set y axis limits for the wavefunction plot.
    frame_delay = num: set the delay between frames in milliseconds (default:200).
    
    Any other keywords will be passed as plot options (you can use for example
    xlim=(num,num)).
    """
    
    fig = plt.figure()
    #create axes for potential plot
    potplot = plt.gca()
    potplot.plot(x_grid, V_array, color="gray")
    potplot.set(xlabel='Real space $/a_0$',
                #xlabel='Computational vector of length $2^{}={}$ $(|n\\rangle)$'.format(round(np.log(len(x))/np.log(2)), len(x)),
                ylabel="$V /E_h$",
                ylim=V_ylim,
                **options)
    
    # create another axes set for the wavefunction sharing the x axis
    psiplot = potplot.twinx()
    
    def dyn_plot(step, sim_array):
        psiplot.cla() # clear current axes
        for particle in sim_array:
            if magnitude:
                psiplot.plot(x_grid, np.abs(particle[step]), color="black", label="$|\\Psi|$", linestyle=":", marker=".")
            if real:
                psiplot.plot(x_grid, np.real(particle[step]), color="blue", label="$Re(\\Psi)$", linestyle=":", marker=".")
            if imaginary:
                psiplot.plot(x_grid, np.imag(particle[step]), color="red", label="$Im(\\Psi)$", linestyle=":", marker=".")
                
        if type(t_grid[1])==float:
            psiplot.text(0.1, 0.9, "t="+str(np.round(t_grid[step]))+"au", transform=psiplot.transAxes)
        else:
            psiplot.text(0.1, 0.9, "t="+str(np.round(t_grid[step]))+"au", transform=psiplot.transAxes)

        psiplot.set(ylabel="$\\psi$", ylim=psi_ylim, **options)
        psiplot.legend(loc="upper right")
        
        return psiplot

    return animation.FuncAnimation(fig, dyn_plot, frames=len(dyn[0]), interval=frame_delay, fargs=(dyn,))