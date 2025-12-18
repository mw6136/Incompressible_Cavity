import numpy as np
import matplotlib.pyplot as plt
import glob
from tqdm import tqdm
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argparse

field_key = {
    "Density": "Density",
    "VelX1": "X-Velocity",
    "VelX2": "Y-Velocity",
    "Press": "Pressure",
    "Temp": "Temperature",
    "VelMag": "Velocity Magnitude",
    "dpdx": "dpdx",
    "dpdy": "dpdy",
}

parser = argparse.ArgumentParser()

parser.add_argument("-d", "--directory", help="Directory containing plotting files")

parser.add_argument("-f", "--field", help="Field to plot")

parser.add_argument("-n", "--ngrid", help="Number of grid cells of simulation")

parser.add_argument(
    "--dt", help="Number of seconds (in omega * t) between data output saves"
)

parser.add_argument(
    "-g", "--ng", help="Include ghost zones in plot", action="store_true"
)

parser.add_argument(
    "-m", "--multi", help="Generate multipanel plot of Mach number", action="store_true"
)
args = parser.parse_args()

N = int(args.ngrid)
# print(N)

gmax, gmin = 0.0, 1e20

for f in tqdm(
    sorted(glob.glob(f"{args.directory}/*.h5")), desc="Getting max and min vals :"
):

    strnum = f.split(".")[-2]
    # print(float(strnum))
    ds = h5py.File(f, "r")

    field = args.field

    if args.field == "VelMag":
        u = np.rot90(ds["VelX1"][:].reshape(N + 2, N + 2)[1:-1, 1:-1], 3)
        v = np.rot90(ds["VelX2"][:].reshape(N + 2, N + 2)[1:-1, 1:-1], 3)

        plotarray = np.sqrt(u**2 + v**2)
    else:
        plotarray = np.rot90(ds[field][:].reshape(N + 2, N + 2)[1:-1, 1:-1], 3)

    if args.field == "Press":
        print(
            (plotarray[1, :] - plotarray[0, :]).min(),
            (plotarray[:, 1] - plotarray[:, 0]).min(),
        )
    pmax = np.amax(plotarray)
    pmin = np.amin(plotarray)
    if pmax > gmax:
        gmax = pmax
    else:
        gmax = gmax

    if pmin < gmin:
        gmin = pmin
    else:
        gmin = gmin


for f in tqdm(
    sorted(glob.glob(f"{args.directory}/*.h5")),
    desc=f"Plotting {field_key[args.field]} :",
):

    strnum = f.split(".")[-2]
    ds = h5py.File(f, "r")

    field = args.field

    if args.field == "VelMag":
        if args.ng:
            u = np.rot90(ds["VelX1"][:].reshape(N + 2, N + 2)[:, :], 3)
            v = np.rot90(ds["VelX2"][:].reshape(N + 2, N + 2)[:, :], 3)
        else:
            u = np.rot90(ds["VelX1"][:].reshape(N + 2, N + 2)[1:-1, 1:-1], 3)
            v = np.rot90(ds["VelX2"][:].reshape(N + 2, N + 2)[1:-1, 1:-1], 3)

        plotarray = np.sqrt(u**2 + v**2)
    else:
        if args.ng:
            plotarray = np.rot90(ds[field][:].reshape(N + 2, N + 2)[:, :], 3)
        else:
            plotarray = np.rot90(ds[field][:].reshape(N + 2, N + 2)[1:-1, 1:-1], 3)

    fig, ax = plt.subplots(figsize=(8, 8))

    x, y = np.linspace(0, 1, len(plotarray)), np.linspace(0, 1, len(plotarray))
    X, Y = np.meshgrid(x, y)

    plt.title(
        rf"{field_key[field]}, $\omega t$ = {round(float(strnum) * float(args.dt),7)}",
        fontsize=16,
    )

    im = ax.pcolormesh(X, Y, plotarray, cmap="jet", vmax=gmax, vmin=gmin)

    # if args.field == "VelMag":
    #     ax.quiver(X,Y,u,v)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.02)
    fig.add_axes(cax)
    fig.colorbar(im, cax=cax)

    plt.savefig(f"{field}_{strnum}_lid.png", bbox_inches="tight", dpi=1200)
    plt.show()
    plt.close()

if args.multi:

    strnum = "0100"
    strnum_08 = "0020"

    f_08 = f"C:\Users\micha\OneDrive\Documents\CTRFL\MAE557\Project2\Ma_08\test.output.{strnum_08}.h5"
    f_04 = f"C:\Users\micha\OneDrive\Documents\CTRFL\MAE557\Project2\Ma_04\test.output.{strnum}.h5"
    f_02 = f"C:\Users\micha\OneDrive\Documents\CTRFL\MAE557\Project2\Ma_02\test.output.{strnum}.h5"
    f_01 = f"C:\Users\micha\OneDrive\Documents\CTRFL\MAE557\Project2\Ma_01\test.output.{strnum}.h5"

    ds_08 = h5py.File(f_08, "r")
    ds_04 = h5py.File(f_04, "r")
    ds_02 = h5py.File(f_02, "r")
    ds_01 = h5py.File(f_01, "r")

    field = args.field

    if args.field == "VelMag":
        u_08 = np.rot90(ds_08["VelX1"][:].reshape(N + 2, N + 2)[1:-1, 1:-1], 3)
        v_08 = np.rot90(ds_08["VelX2"][:].reshape(N + 2, N + 2)[1:-1, 1:-1], 3)

        plotarray_08 = np.sqrt(u_08**2 + v_08**2)

        u_04 = np.rot90(ds_04["VelX1"][:].reshape(N + 2, N + 2)[1:-1, 1:-1], 3)
        v_04 = np.rot90(ds_04["VelX2"][:].reshape(N + 2, N + 2)[1:-1, 1:-1], 3)

        plotarray_04 = np.sqrt(u_04**2 + v_04**2)

        u_02 = np.rot90(ds_02["VelX1"][:].reshape(N + 2, N + 2)[1:-1, 1:-1], 3)
        v_02 = np.rot90(ds_02["VelX2"][:].reshape(N + 2, N + 2)[1:-1, 1:-1], 3)

        plotarray_02 = np.sqrt(u_02**2 + v_02**2)

        u_01 = np.rot90(ds_01["VelX1"][:].reshape(N + 2, N + 2)[1:-1, 1:-1], 3)
        v_01 = np.rot90(ds_01["VelX2"][:].reshape(N + 2, N + 2)[1:-1, 1:-1], 3)

        plotarray_01 = np.sqrt(u_01**2 + v_01**2)
    else:
        plotarray = np.rot90(ds[field][:].reshape(N + 2, N + 2)[1:-1, 1:-1], 3)

    fig, ax = plt.subplots(2, 2, figsize=(8, 8), sharex=True, sharey=True)

    im_08 = ax[0, 0].pcolormesh(
        X,
        Y,
        plotarray_08,
        cmap="jet",
        vmax=np.amax(plotarray_08),
        vmin=np.amin(plotarray_08),
    )
    ax[0, 0].set_title(r"$\mathrm{Ma} = 0.8$")

    im_04 = ax[0, 1].pcolormesh(
        X,
        Y,
        plotarray_04,
        cmap="jet",
        vmax=np.amax(plotarray_04),
        vmin=np.amin(plotarray_04),
    )
    ax[0, 1].set_title(r"$\mathrm{Ma} = 0.4$")

    im_02 = ax[1, 0].pcolormesh(
        X,
        Y,
        plotarray_02,
        cmap="jet",
        vmax=np.amax(plotarray_02),
        vmin=np.amin(plotarray_02),
    )
    ax[1, 0].set_title(r"$\mathrm{Ma} = 0.2$")

    im_01 = ax[1, 1].pcolormesh(
        X,
        Y,
        plotarray_01,
        cmap="jet",
        vmax=np.amax(plotarray_01),
        vmin=np.amin(plotarray_01),
    )
    ax[1, 1].set_title(r"$\mathrm{Ma} = 0.1$")

    divider = make_axes_locatable(ax[0, 1])
    cax = divider.append_axes("right", size="5%", pad=0.02)
    fig.add_axes(cax)
    fig.colorbar(im_04, cax=cax)

    divider = make_axes_locatable(ax[1, 1])
    cax = divider.append_axes("right", size="5%", pad=0.02)
    fig.add_axes(cax)
    fig.colorbar(im_01, cax=cax)

    divider = make_axes_locatable(ax[0, 0])
    cax = divider.append_axes("right", size="5%", pad=0.02)
    fig.add_axes(cax)
    fig.colorbar(im_08, cax=cax)

    divider = make_axes_locatable(ax[1, 0])
    cax = divider.append_axes("right", size="5%", pad=0.02)
    fig.add_axes(cax)
    fig.colorbar(im_02, cax=cax)

    plt.savefig("multiplot.png", bbox_inches="tight", dpi=1200)
    plt.show()
    plt.close()
