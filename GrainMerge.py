import matplotlib.pyplot as plt
import numpy as np
from mtex import *

def select_grains(grains, matName):
    global stage, selected_grains, selected_parents
    stage = 0
    ax = plt.gca()
    ax.set_title("Select mis-indexed grains (RED) - ENTER to continue")
    fig = plt.gcf()
    selected_grains = []
    selected_parents = []
    rounds = []

    def onclick(event):
        global stage, selected_grains, selected_parents
        if event.inaxes is ax:
            x, y = event.xdata, event.ydata
            g = grains(x, y)
            if len(g) == 0:
                return

            gid = int(np.asarray(g.id).reshape(-1)[0])   # grain ID
            ind = int(np.asarray(grains.id2ind(gid)).reshape(-1)[0])  # Position for displaying grain

            if stage == 0:
                if gid not in selected_grains:
                    selected_grains.append(gid)
                    hold(True)
                    plot(grains[ind].boundary, lineWidth=4, lineColor='red')
                    hold(False)
                    fig.canvas.draw_idle()

            elif stage == 1:
                if gid not in selected_parents:
                    selected_parents.append(gid)
                    hold(True)
                    plot(grains[ind].boundary, lineWidth=4, lineColor='green')
                    hold(False)
                    fig.canvas.draw_idle()   

    def onkey(event):
        global stage, selected_grains, selected_parents
        if stage == 0:
            if event.key in ('enter', 'return'):
                stage = 1
                ax.set_title("Select parent grains (GREEN) - ENTER to finish")
                print(f"\nChildren: {selected_grains}\n")
                fig.canvas.draw_idle()

        elif stage == 1:
            if event.key in ('enter', 'return'):
                rounds.append((selected_grains[:], selected_parents[:]))
                print(f"Parents: {selected_parents}\n")
                stage = 2
                fig.canvas.draw_idle()

        elif stage == 2:
            fig.canvas.mpl_disconnect(cid)
            plt.close(fig)

    fig.canvas.mpl_connect('button_press_event', onclick)
    plt.gcf().canvas.manager.set_window_title(f'Select Grains - {matName}')
    cid = fig.canvas.mpl_connect('key_press_event', onkey)
    plt.show()

    return selected_grains, selected_parents

def find_boundary_index(boundaries, grain_a, grain_b):
    """Return the row index where [grain_a, grain_b] (or [grain_b, grain_a]) appears, else -1."""
    mask = ((boundaries[:, 0] == grain_a) & (boundaries[:, 1] == grain_b)) | \
           ((boundaries[:, 0] == grain_b) & (boundaries[:, 1] == grain_a))
    idx = np.where(mask)[0]
    return idx[0] if len(idx) else -1

def grainmerge(ebsd, grains, matName, ipfkey):
    plot(ebsd[matName], ebsd[matName].orientations, micronbar='on')
    hold(True)
    largeGrains = grains[grains.numPixel > 20]
    text(largeGrains, largeGrains.id)
    children, parents = select_grains(grains, matName) # Pull children and parent grains out to combine them.
    print(children, parents) # We will just call the function again if user wants to make more selections.
    pairs = grains.neighbors() # This is important for validating user selections.

    idx = find_boundary_index(pairs, children, parents) # Finds specific boundary between selected grains.
    print(f"\nFound at index {idx}\n" if idx != -1 else "\n Not found\n")

    # Find the largest selected parent grain and grab its mean orientation
    L_parent = 0
    L_parent_meanO = None
    for i in range(len(parents)):
        parent_area = float(np.asarray(grains['id', parents[i]].area).item())
        if parent_area > L_parent:
            L_parent = parent_area
            L_parent_meanO = grains['id', parents[i]].meanOrientation

    ind = np.isin(ebsd.grainId, children) # Index child grains
    ebsd.rotations.data[ind] = L_parent_meanO.data # Set children orientations to the parent's mean orientation.
    grains = calcGrains(ebsd, angle=5 * degree, minPixel=10) # Recalculate grains

    plot(ebsd[matName], ipfkey.orientation2color(ebsd[matName].orientations), micronbar='on')
    hold(True)
    plot(grains.boundary, lineWidth=1)
    hold(False)
    plt.gcf().canvas.manager.set_window_title(f'Merging Grains - {matName}')
    plt.show()

    return ebsd