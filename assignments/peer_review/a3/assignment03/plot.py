from numpy import loadtxt
from matplotlib.pyplot import subplots

# load the raw data
x_data, y_data = loadtxt("DATA_FILES/DATA.dat", unpack=True)

# load the fit data
x_fit, y_fit = loadtxt("q5_output.dat", unpack=True)

# begin defining figure
fig, ax = subplots(1, 1, figsize=(10,10))

# plot the raw data from the periodic process
ax.scatter(x_data, y_data, label="Raw Data")

# plot the fit data from the chi squared minimization
ax.plot(x_fit, y_fit, color="black", lw=3.0, label="Fit")

ax.legend(loc=0, prop={"size": 18})
ax.set_ylabel(r"$y$", fontsize=20)
ax.set_xlabel(r"$x$", fontsize=20)

fig.tight_layout()
fig.savefig("q5_plot.pdf")
