import numpy as np
import cv2 
import matplotlib.pyplot as plt

def rotation_angle(R):
    return np.rad2deg(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1)))

def angle(v1,v2):
    return np.rad2deg(np.arccos(np.clip(v1.dot(v2) / np.linalg.norm(v1) / np.linalg.norm(v2), -1, 1)))

def skew(v):
    return np.array([[0,-v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])

def read_grayscale_image(path):
    image = cv2.imread(path)
    image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    return image


def savefig(filepath, fig=None, dpi=None, remove_axis=True):
    # TomNorway - https://stackoverflow.com/a/53516034
    if not fig:
        fig = plt.gcf()

    plt.subplots_adjust(0, 0, 1, 1, 0, 0)
    if remove_axis:
        for ax in fig.axes:
            ax.axis('off')
            ax.margins(0, 0)
            ax.xaxis.set_major_locator(plt.NullLocator())
            ax.yaxis.set_major_locator(plt.NullLocator())

    fig.savefig(filepath, pad_inches=0, bbox_inches='tight', dpi=dpi)