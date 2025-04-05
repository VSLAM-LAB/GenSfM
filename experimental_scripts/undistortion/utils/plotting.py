import matplotlib.pyplot as plt


def plot_1d_radial_reprojs(R,t,pp,points3D,width,height,**kwargs):
    
    for X in points3D:
        z = (R @ X + t)[0:2]

        # pp + alpha * z = d
        # lambda = (d-pp)/z
        alphas = [-pp[0]/z[0], -pp[1]/z[1], (width-pp[0])/z[0], (height-pp[1])/z[1]]
        alpha = min([alpha for alpha in alphas if alpha > 0])

        plt.plot([pp[0], pp[0]+alpha*z[0]], [pp[1], pp[1]+alpha*z[1]],**kwargs)
        