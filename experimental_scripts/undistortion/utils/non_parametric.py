import numpy as np
import pyimplicitdist
import poselib
import random
import torch

def project_non_parametric_intrinsic(R,t,X,calib):
    Z = (R @ X.T).T + t
    x = pyimplicitdist.distort(Z, calib)
    return x


def ransac_pose(p2d,p3d,calib,thr=12.0,max_iters=100):
    N = len(p2d)
    best_score = np.inf
    best_model = None
    best_inl = []
    thr_sq = thr * thr

    p2d_undist = np.array(pyimplicitdist.undistort(p2d, calib))

    for iter in range(max_iters):
        sample = random.sample(range(N), 3)

        models = poselib.p3p(p2d_undist[sample], p3d[sample])

        for p in models:            
            proj_p2d = project_non_parametric_intrinsic(p.R,p.t,p3d,calib)
            err = np.sum((proj_p2d-p2d)**2,axis=1)
            score = np.sum(np.clip(err,0,thr_sq))

            if score < best_score:
                best_score = score
                best_inl = err < thr_sq
                best_model = p
    return best_model, best_inl
        

def undistort_image(img, calib, width, height, hfov=70, rotation=np.eye(3)):
    # compute focal length
    focal = width / (2.0 * np.tan(np.deg2rad(hfov)/2))
    print(f'focal = {focal}')

    if len(img.shape) == 2:
        img = img[None]
    else:
        # We have three channels - permute
        img = np.transpose(img,axes=[2,0,1])

    if img.dtype == np.uint8:
        img = img.astype(np.float64)
        img /= 255.0
    
    # Coordinates in the original image
    ind_x = range(width)
    ind_y = range(height)
    xx, yy = np.meshgrid(ind_x, ind_y, indexing='xy')

    xx = (xx - width/2) / focal
    yy = (yy - height/2) / focal
    

    X = xx.flatten()
    Y = yy.flatten()
    Z = np.ones(X.shape[0])

    points3D = np.r_[X[None],Y[None],Z[None]].T

    points3D = (rotation @ points3D.T).T

    points2D = np.array(pyimplicitdist.distort(points3D, calib))
   
    x = np.reshape(points2D[:,0],(height,width))
    y = np.reshape(points2D[:,1],(height,width))

    #import ipdb
    #ipdb.set_trace()

    c, h, w = img.shape
    x = (x - w/2) / (w/2)
    y = (y - h/2) / (h/2)

    # Note that we flip x,y here since we index by [y,x] in the image
    points2D = np.stack((x,y),axis=-1)

    # Add batch axis
    img = img[None]
    points2D = points2D[None]

    # and interpolate
    im_undist = torch.nn.functional.grid_sample(torch.from_numpy(img), torch.from_numpy(points2D))
    im_undist = im_undist.cpu().detach().numpy()[0]
    im_undist = np.transpose(im_undist, axes=[1,2,0])
    return im_undist
    
    
