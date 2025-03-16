import glob
import scipy.io
# import pytsamp
import numpy as np
from itertools import product

def unproject(pts, camera_dict):
    """
    Convert 2D image points to 3D bearing vectors based on camera calibration parameters.

    Args:
        pts (np.array): An Nx2 array of image points in pixel coordinates.
        camera_dict (dict): Dictionary containing camera intrinsic parameters and model:
                            - 'model': The camera model ('OPENCV_FISHEYE' or 'SIMPLE_RADIAL').
                            - 'params': List of model-specific parameters.

    Returns:
        np.array: An Nx3 array of 3D unit-bearing vectors corresponding to the 2D points.
    """
    pts = np.array(pts)
    model = camera_dict['model']
    params = camera_dict['params']

    if model == 'OPENCV_FISHEYE':
        # 'kb' model (Kannala-Brandt)
        fx, fy, cx, cy, k1, k2, k3, k4 = params
        # Normalize image coordinates to the camera frame
        x = (pts[:, 0] - cx) / fx
        y = (pts[:, 1] - cy) / fy
        r = np.sqrt(x**2 + y**2)

        # Compute theta (the undistorted angle)
        theta_d = r
        theta = theta_d / (1 + k1 * theta_d**2 + k2 * theta_d**4 + k3 * theta_d**6 + k4 * theta_d**8)
        scale = np.tan(theta) / (r + 1e-9)  # Add small epsilon to avoid division by zero
        # scale = np.sin(theta) / (r + 1e-9)
        # Convert to bearing vectors
        bearing_vecs = np.vstack((x * scale, y * scale, np.ones_like(x))).T
        bearing_vecs /= np.linalg.norm(bearing_vecs, axis=1, keepdims=True)

    elif model == 'SIMPLE_RADIAL':
        # 'bc' model (Brown-Conrady)
        f, cx, cy, k1 = params
        # Normalize image coordinates to the camera frame
        x = (pts[:, 0] - cx) / f
        y = (pts[:, 1] - cy) / f
        r = np.sqrt(x**2 + y**2)

        # Apply inverse radial distortion
        radial_distortion = 1 + k1 * r**2
        x /= radial_distortion
        y /= radial_distortion

        # Convert to bearing vectors
        bearing_vecs = np.vstack((x, y, np.ones_like(x))).T
        bearing_vecs /= np.linalg.norm(bearing_vecs, axis=1, keepdims=True)

    else:
        raise ValueError(f"Unsupported camera model: {model}")

    return bearing_vecs
class CameraPose:
    def __init__(self):
        self.R = np.eye(3)
        self.t = np.zeros(3)
                

def calib_to_dict(calib):
    if calib['proj_model'] == 'kb':
        return {
            'model': 'OPENCV_FISHEYE',
            'width': -1, 
            'height': -1,
            'params': [calib['K'][0,0], calib['K'][1,1], calib['K'][0,2], calib['K'][1,2], calib['proj_params'][0], calib['proj_params'][1], calib['proj_params'][2], calib['proj_params'][3]]
        }
    elif calib['proj_model'] == 'bc':
        return {
            'model': 'SIMPLE_RADIAL',
            'width': -1, 
            'height': -1,
            'params': [0.5*(calib['K'][0,0]+calib['K'][1,1]),calib['K'][0,2], calib['K'][1,2], calib['proj_params']]
        }
    else:
        print('Unsupported camera model!')
        return {}
    

def project(points_3d, camera_dict):
    """
    Projects 3D points in the camera coordinate system onto the 2D image plane
    using different projection models.

    Args:
        points_3d (np.array): Nx3 array of 3D points in the camera coordinate system.
        camera_dict (dict): Dictionary containing camera model and intrinsic parameters:
            - 'model': The camera model ('OPENCV_FISHEYE' or 'SIMPLE_RADIAL').
            - 'params': List of model-specific parameters.

    Returns:
        np.array: Nx2 array of 2D points in pixel coordinates.
    """
    model = camera_dict['model']
    params = camera_dict['params']

    # Normalize by Z to get normalized image coordinates
    points_2d = points_3d[:, :2] / points_3d[:, 2][:, np.newaxis]  # (N, 2)
    r = np.linalg.norm(points_2d, axis=1)  # Radius in normalized coordinates

    if model == 'OPENCV_FISHEYE':
        # 'kb' model (Kannala-Brandt)
        fx, fy, cx, cy, k1, k2, k3, k4 = params
        theta = np.arctan(r)
        theta_d = theta * (1 + k1 * theta**2 + k2 * theta**4 + k3 * theta**6 + k4 * theta**8)
        scaling = theta_d / (r + 1e-9)  # Add small epsilon to avoid division by zero
        points_2d *= scaling[:, np.newaxis]
        points_2d[:, 0] = fx * points_2d[:, 0] + cx
        points_2d[:, 1] = fy * points_2d[:, 1] + cy

    elif model == 'SIMPLE_RADIAL':
        # 'bc' model (Brown-Conrady)
        f, cx, cy, k1 = params
        radial_distortion = 1 + k1 * r**2
        points_2d *= radial_distortion[:, np.newaxis]
        points_2d[:, 0] = f * points_2d[:, 0] + cx
        points_2d[:, 1] = f * points_2d[:, 1] + cy

    else:
        raise ValueError(f"Unsupported camera model: {model}")

    return points_2d

# Note, this does not work for all cameras as their distortion model does not fit the image edges
def normalize_calib_data(data, gt_poses):
    # Find diagonal of the largest board
    boards = data['boards']
    if type(boards) != list:
        boards = [boards]

    diags = []
    for b in boards:
        dists = [np.linalg.norm(a-b) for a,b in product(b['X'].T,b['X'].T)]        
        diags.append(np.max(dists))            
    scale = np.max(diags)
    
    # normalize 3D points w.r.t. board size
    for b in boards:
        b['X'] = b['X'].astype(np.float64)     
        b['X'] /= scale
        b['Rt'] = b['Rt'].astype(np.float64)        
        b['Rt'][:,3] /= scale
    
    # normalize camera poses
    for k in range(len(gt_poses)):
        gt_poses[k].t /= scale

def compute_fov(calib, imgsize):
    camera_dict = calib_to_dict(calib)
    pp = [calib['K'][0,2], calib['K'][1,2]]
    pts = np.array([[pp[0], pp[1]], [imgsize[1], pp[1]], [0.0, pp[1]]])
    bearing_vecs = unproject(pts, camera_dict)
    fov1 = np.arccos(np.dot(bearing_vecs[0], bearing_vecs[1]))
    fov2 = np.arccos(np.dot(bearing_vecs[0], bearing_vecs[2]))

    return (fov1+fov2) * 180 / np.pi    

# def compute_effective_fov(calib, pts):
#     camera_dict = calib_to_dict(calib)
#     pp = [calib['K'][0,2], calib['K'][1,2]]
#     all_pts = np.concatenate([x[0] for x in pts])

#     center_vec = np.array(unproject([pp], camera_dict))
#     bearing_vecs = np.array(unproject(all_pts, camera_dict))

#     cos_theta = bearing_vecs @ center_vec.T
#     theta = np.arccos(cos_theta)
#     fov = 2 * np.max(theta) * 180.0 / np.pi

#     return fov   

def compute_effective_fov(calib, pts):
    camera_dict = calib_to_dict(calib)
    
    # Ensure principal point is formatted correctly
    pp = np.array([[calib['K'][0,2], calib['K'][1,2]]])  

    # Ensure points are properly stacked
    all_pts = np.vstack([x[0] for x in pts])

    # Compute bearing vectors
    center_vec = np.array(unproject(pp, camera_dict)).squeeze()
    bearing_vecs = np.array(unproject(all_pts, camera_dict))

    # Compute angles
    cos_theta = np.dot(bearing_vecs, center_vec)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)  # Ensure stability in arccos

    # theta = np.arccos(cos_theta)
    theta = np.arctan2(np.linalg.norm(np.cross(bearing_vecs, center_vec), axis=1),
                   np.dot(bearing_vecs, center_vec))
    fov = 2 * np.max(theta) * 180.0 / np.pi

    return fov

def load_data(query_file, max_error=3.0):
    import glob
    import scipy.io
    import numpy as np

    dataset = scipy.io.loadmat(query_file, simplify_cells=True)
    # split the dataset into train_set and test_set based on dataset['img'] name arrays, select the components of the dataset where 'img' have 'train' in their name for train_set and 'test' for test_set
    dataset_train = {'img': [], 'orpc': [], 'corners': [], 'boards': dataset['boards'], 'imgsize': dataset['imgsize']}
    dataset_test = {'img': [], 'orpc': [], 'corners': [], 'boards': dataset['boards'], 'imgsize': dataset['imgsize']}
    # breakpoint()
    



    dataset_name = query_file[query_file.rfind('/') + 1:-4]
    query_path = query_file[:query_file.rfind('/') + 1:]

    # Check if each calibration file exists
    bc_path = query_path + 'calib/' + dataset_name + '_calib_bc.mat'
    kb_path = query_path + 'calib/' + dataset_name + '_calib_kb.mat'
    bc_exists = glob.glob(bc_path)
    kb_exists = glob.glob(kb_path)

    if not (bc_exists or kb_exists):
        print('No calibration files found!')
        return None

    # Load available calibration files
    if bc_exists and kb_exists:
        calibs = {
            'bc': scipy.io.loadmat(bc_path, simplify_cells=True),
            'kb': scipy.io.loadmat(kb_path, simplify_cells=True)
        }
        # Choose the simpler model if it performs similarly
        if calibs['bc']['calib_ir'] >= calibs['kb']['calib_ir'] * 0.99 and calibs['bc']['calib_rms'] <= calibs['kb']['calib_rms'] * 1.01:
            calib = calibs['bc']
        else:
            calib = calibs['kb']
    elif bc_exists:
        calib = scipy.io.loadmat(bc_path, simplify_cells=True)
    elif kb_exists:
        calib = scipy.io.loadmat(kb_path, simplify_cells=True)
    print(calib['Rt'].shape)
    
    calib_train = {'proj_params': calib['proj_params'], 'K': calib['K'], 'Rt': [], 'calib_rms': calib['calib_rms'], 'calib_ir': calib['calib_ir'], 'proj_model': calib['proj_model']}
    calib_test = {'proj_params': calib['proj_params'], 'K': calib['K'], 'Rt': [], 'calib_rms': calib['calib_rms'], 'calib_ir': calib['calib_ir'], 'proj_model': calib['proj_model']}
    # Split the calibration data into train_set and test_set based on the dataset['img'] name arrays
    # Collect correspondences and pose data
    corrs = collect_correspondences(dataset)
    # poses = []
    for i, img_path in enumerate(dataset['img']):
        if 'train' in img_path:
            # Populate dataset_train
            dataset_train['img'].append(img_path)
            dataset_train['orpc'].append(dataset['orpc'][i])
            dataset_train['corners'].append(dataset['corners'][i])
            
            # Populate calib_train with corresponding calibration data
            # calib_train['proj_params'].append(calib['proj_params'])
            # calib_train['K'].append(calib['K'])
            calib_train['Rt'].append(calib['Rt'][:, :, i])
            # calib_train['calib_rms'].append(calib['calib_rms'])
            # calib_train['calib_ir'].append(calib['calib_ir'])
        else:
            # Populate dataset_test
            dataset_test['img'].append(img_path)
            dataset_test['orpc'].append(dataset['orpc'][i])
            dataset_test['corners'].append(dataset['corners'][i])
            
            # Populate calib_test with corresponding calibration data
            # calib_test['proj_params'].append(calib['proj_params'])
            # calib_test['K'].append(calib['K'])
            calib_test['Rt'].append(calib['Rt'][:, :, i])
            # calib_test['calib_rms'].append(calib['calib_rms'])
            # calib_test['calib_ir'].append(calib['calib_ir'])

    # Convert lists to numpy arrays
    calib_train['Rt'] = np.stack(calib_train['Rt'], axis=2)
    calib_test['Rt'] = np.stack(calib_test['Rt'], axis=2)
    dataset_train['img'] = np.array(dataset_train['img'], dtype=object)
    dataset_train['orpc'] = np.array(dataset_train['orpc'], dtype=object)
    # dataset_train['corners'] = np.array(dataset_train['corners'], dtype=object)
    
    dataset_test['img'] = np.array(dataset_test['img'], dtype=object)
    dataset_test['orpc'] = np.array(dataset_test['orpc'], dtype=object)
    # dataset_test['corners'] = np.array(dataset_test['corners'], dtype=object)
    print(calib_train['Rt'].shape)
    
    # Convert calibration data lists to numpy arrays for train and test sets
    for key in ['proj_params', 'K', 'Rt', 'calib_rms', 'calib_ir']:
        calib_train[key] = np.array(calib_train[key])
        calib_test[key] = np.array(calib_test[key])

    corrs_train = collect_correspondences(dataset_train)
    corrs_test = collect_correspondences(dataset_test)
    poses_train = []
    poses_test = []
    for k in range(calib_train['Rt'].shape[2]):
        p = CameraPose()
        p.R = calib_train['Rt'][:, 0:3, k]
        p.t = calib_train['Rt'][:, 3, k]
        poses_train.append(p)
    for k in range(calib_test['Rt'].shape[2]):
        p = CameraPose()
        p.R = calib_test['Rt'][:, 0:3, k]
        p.t = calib_test['Rt'][:, 3, k]
        poses_test.append(p)
    poses = []
    for k in range(calib['Rt'].shape[2]):
        pose = CameraPose()
        pose.R = calib['Rt'][:, 0:3, k]
        pose.t = calib['Rt'][:, 3, k]
        poses.append(pose)
    # Filter correspondences
    corrs_train = filter_corrs(corrs_train, poses_train, calib_train, max_error)
    corrs_test = filter_corrs(corrs_test, poses_test, calib_test, max_error)


    # Compute field of view
    # fov = compute_fov(calib, dataset['imgsize'])
    fov = compute_effective_fov(calib, corrs)

    # return dataset_train, dataset_test, calib_train, calib_test, dataset, calib, poses_train, poses_test, poses
    return dataset_train, dataset_test, corrs_train, corrs_test, poses_train, poses_test, calib_train, calib_test, fov

# Grabs all 2D-3D correspondences from the data
def collect_correspondences(data):    
    boards = data['boards']
    if type(boards) != list:
        boards = [boards]

    num_pts = [b['X'].shape[1] for b in boards]
    offset = np.cumsum([0] + num_pts)
    corrs = []

    for im in data['corners']:
        p2d = im['x'].transpose()
        idx = im['cspond'].transpose()    
        idx = idx - 1
        p3d = []
        full_idx = idx[:,0] + offset[idx[:,1]]
        for (pt_idx, b_idx) in idx:
            b = boards[b_idx]                        
            
            X = b['Rt'] @ np.r_[b['X'][:,pt_idx], 0.0, 1.0]
            
            p3d.append(X)
            
        corrs.append((np.array(p2d), np.array(p3d), full_idx, idx[:,1]))
    return corrs

def filter_corrs(corrs, poses, calib, max_error):
    camera_dict = calib_to_dict(calib)

    for k in range(len(corrs)):
        pose = poses[k]
        (p2d, p3d, full_idx, idx) = corrs[k]
        
        X_trans = p3d @ pose.R.T + pose.t.T
        p2d_proj = np.array(project(X_trans, camera_dict))
        err = np.sqrt(np.sum((p2d_proj - p2d)**2,axis=1))

        inliers = np.array(err < max_error)

        corrs[k] = (p2d[inliers], p3d[inliers], full_idx[inliers], idx[inliers])

        #print(f"orig = {np.mean(err)}")
        # debug
        #min_err = 1000
        #for p in poses:
        #    X_trans = p3d @ p.R.T + p.t.T
        #    p2d_proj = np.array(pytsamp.project(X_trans, camera_dict))
        #    err = np.sqrt(np.sum((p2d_proj - p2d)**2,axis=1))
        #    err = np.mean(err)
        #    print(err)
        #    if err < min_err:
        #        min_err = err

    return corrs

