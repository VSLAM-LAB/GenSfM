import numpy as np

def colmap_camera_to_dict(camera):
    camera_dict = {}
    camera_dict['model'] = camera.model
    camera_dict['width'] = camera.width
    camera_dict['height'] = camera.height
    camera_dict['params'] = camera.params
    return camera_dict

def project_kb(R,t,X,K,params):    
    Xcam = [R @ U + t for U in X]
    return apply_kb(Xcam,K,params)

def apply_kb(Xcam,K,params):
    x = []
    for Z in Xcam:        
        r = np.linalg.norm(Z[0:2])
        z = Z[2]
                
        if(r > 1e-8):
            theta = np.arctan2(r,z)
            theta2 = theta * theta
            theta4 = theta2 * theta2
            theta6 = theta4 * theta2
            theta8 = theta6 * theta2        
            theta_d = theta * (1.0 + params[0] * theta2 + params[1] * theta4 + params[2] * theta6 + params[3] * theta8)
            x_proj = np.r_[(theta_d/r) * Z[0:2], 1.0]
        else:
            x_proj = np.r_[Z[0:2], 1.0]

        x_proj = K @ x_proj
        x.append(x_proj[0:2] / x_proj[2])        
    return np.array(x)


def rms_kb(x,R,t,X,K,params):
    x_proj = project_kb(R,t,X,K,params)
    return np.sqrt(np.mean([np.linalg.norm(z1-z2)**2 for z1,z2 in zip(x,x_proj)]))

def rms_radial(x,R,t,X,pp):
    res2 = []

    for (z,Z) in zip(x,X):
        Z = R @ Z + t
        zz = Z[0:2] / np.linalg.norm(Z[0:2])

        res2.append(np.linalg.norm(zz.dot(z-pp) * zz - (z-pp))**2)
    
    return np.sqrt(np.mean(res2))