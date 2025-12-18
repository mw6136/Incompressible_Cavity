import numpy as np
import matplotlib.pyplot as plt
import glob
from tqdm import tqdm

ylo, yhi = 0, 0

for f in tqdm(sorted(glob.glob("*.csv")), "Generating plots:"):

    savenum = f.split(".")[-2]
    spat = f.split(".")[1]
    time = f.split(".")[0]

    ds = np.genfromtxt(f, delimiter=",")

    if savenum == "0000":
        ylo, yhi = ds[:-1].min() - 0.1 * abs(ds[:-1].min()), ds[:-1].max() + 0.1 * abs(
            ds[:-1].max()
        )

    fig = plt.figure(figsize=(10, 6))

    plt.plot(np.linspace(0, 2 * np.pi, ds[:-1].size), ds[:-1])
    plt.ylim(ylo, yhi)

    plt.xlabel("x", fontsize=16)
    plt.ylabel("u", fontsize=16)
    plt.title(f"Temporal : {time}, Spatial : {spat}\nt = {ds[-1]}")

    plt.savefig(f"plot.{savenum}.png", bbox_inches="tight", dpi=1200)
    plt.close()
	