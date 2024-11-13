# %%
import sys
# sys.path.append('/home/yihan/cvg/implicit_radial_sfm')
sys.path.append('../')
from calib_reader import *
from python.database import *
# from ``scripts/python/database.py`` import everything
# from scripts.python.database import *
import os
import pyimplicitdist
import matplotlib.pyplot as plt




# %%
base_path = './data/'
output_db_base = './data/databases'
output_calib_base = './calibrations'

# dataset_train, dataset_test, calib_train, calib_test, dataset, calib, poses_train, poses_test, poses=load_data('/home/yihan/cvg/implicit_radial_sfm/experimental_scripts/babelcalib/data/OCamCalib_Fisheye190deg.mat')
# load data and split to train and test
# dataset_train, dataset_test, corrs_train, corrs_test, poses_train, poses_test, calib_train, calib_test, fov = load_data('/home/yihan/cvg/implicit_radial_sfm/experimental_scripts/babelcalib/data/UZH_DAVIS_indoor_45.mat')

dataset_train, dataset_test, corrs_train, corrs_test, poses_train, poses_test, calib_train, calib_test, fov = load_data('./data/UZH_DAVIS_indoor_45.mat')


# %%
def rotation_to_quaternion(R):
    q = np.zeros(4)
    q[0] = np.sqrt(1 + R[0, 0] + R[1, 1] + R[2, 2]) / 2
    q[1] = (R[2, 1] - R[1, 2]) / (4 * q[0])
    q[2] = (R[0, 2] - R[2, 0]) / (4 * q[0])
    q[3] = (R[1, 0] - R[0, 1]) / (4 * q[0])
    return q

# %%
# Get intrinsic calibration from implicit distortion model
# get 2d points from corrs_train
def get_intrinsics(dataset_train,corrs_train, poses_train, output_path):
   
    points2d = []
    points3d = []
    for corr in corrs_train:
        points2d.append(corr[0])
        points3d.append(corr[1])
    # points2d = np.array(points2d)
    # points3d = np.array(points3d)
    cm_opt = pyimplicitdist.CostMatrixOptions()
    image_size = dataset_train['imgsize']
    pp = np.array([image_size[1]/2, image_size[0]/2])

    
    
    poses = []
    for pose in poses_train:
        # transform pose.R to quaternion
        q  = rotation_to_quaternion(pose.R)
        p = pyimplicitdist.CameraPose(q, pose.t)
        poses.append(p)
    # get intrinsic calibration
    opt_1dradial = pyimplicitdist.PoseRefinement1DRadialOptions()
    opt_1dradial.verbose = True
    opt_1dradial.optimize_pose = False
    print("original pp", pp)
    out = pyimplicitdist.joint_pose_refinement_1D_radial(points2d, points3d, poses, pp, opt_1dradial)
    pp = out['pp']
    print("updated pp", pp)
    cost_matrix =  pyimplicitdist.build_cost_matrix_multi(points2d, cm_opt, pp)
    
    pyimplicitdist.estimate
    print('Cost matrix built')
    # if not os.path.exists(output_path):
    
    total_num_points = np.sum([x.shape[0] for x in points2d])
    print("Start calibrating", total_num_points, "points")
    if not os.path.exists(output_path):
        intrinsics = pyimplicitdist.calibrate_multi(points2d, points3d, cost_matrix,  pp, poses)
        np.savetxt(output_path, intrinsics.theta_r)
        print('Intrinsic calibration saved to {}'.format(output_path))
    
    return out


# %%


# %%
# # plot the intrisics.theta_r
# import matplotlib.pyplot as plt
# write intrinsics.theta_r to a file
def spline_fit(theta_r_path):
    theta = []
    r = []
    # read theta_r from file
    theta_r = np.loadtxt(theta_r_path)
    for i in range(0, len(theta_r)):
        theta.append(theta_r[i][0])
        r.append(theta_r[i][1]) 
    # fit a cubic spline to theta_r
    import scipy.interpolate
    thetas = np.array(theta)
    rs = np.array(r)
    theta_r_spline = scipy.interpolate.CubicSpline(thetas, rs, bc_type='natural')
    return thetas, rs, theta_r_spline


# %%
# calculate the mean reprojection error for test set
def calculate_mean_reproj_error(corrs_test, poses_test,theta_r_spline, thetas, out):
    mean_reprojection_error = 0
    point_count = 0
    total_points = 0
    
    thetas_full = []
    rs_full = []
    errors_full = []
    for i in range(0, len(corrs_test)):
        corr = corrs_test[i]
        points2d = corr[0]
        points3d = corr[1]
        # get the camera pose
        pose = poses_test[i]
        # for each point, calculate the reprojection error
        for j in range(0, len(points2d)):
            point2d = points2d[j]
            point3d = points3d[j]
            # project 3d point into camera coordinate using pose.R and pose.t
            point3d_cam = pose.R @ point3d + pose.t
            # get the viewing angle theta from the point3d_cam
            rho = np.sqrt(point3d_cam[0]**2 + point3d_cam[1]**2)
            theta = np.arctan2(rho, point3d_cam[2])
            if theta >= thetas[0] and theta <= thetas[-1]:
            # get r from theta_r_spline
                r = theta_r_spline(theta)
                # f = r / (np.tan(theta) + 1e-8)
                # get the projection of 3d point
                r_ori = np.sqrt(point3d_cam[0]**2 + point3d_cam[1]**2)
                u_proj = r * point3d_cam[0] / r_ori + out['pp'][0]
                v_proj = r * point3d_cam[1] / r_ori + out['pp'][1]
                # u_proj = f * point3d_cam[0] / point3d_cam[2] + out['pp'][0]
                # v_proj = f * point3d_cam[1] / point3d_cam[2] + pp[1]
                # calculate the reprojection error
                error = (u_proj - point2d[0])**2 + (v_proj - point2d[1])**2
                # print(error)
                
                if error > 10:
                    continue
                mean_reprojection_error += error
                point_count += 1
                
                thetas_full.append(theta)
                rs_full.append(theta_r_spline(theta))
                errors_full.append(error)
            total_points += 1
    
    if point_count == 0:
        breakpoint()
    
    errors_full = np.array(errors_full)
    inliers = errors_full < 10
    plt.hist(errors_full[inliers], bins=100)
    plt.show()
    
    breakpoint()
    mean_reprojection_error = np.sqrt(mean_reprojection_error / point_count)
    calibrated_points = point_count / total_points
    return mean_reprojection_error, calibrated_points

         

    

# %%
results = {'OV_corner': [], 'OV_cube': [], 'OV_single_plane': [], 'Kalibr': [], 'OCamCalib': [], 'UZH_DAVIS': [], 'UZH_Snapdragon': []}
base_path = './data/'
output_base = './eval/'
spline_base = './spline/'
query_lists = sorted(os.listdir(base_path))[:1]
out_full = {}
# query_lists = [query_file for query_file in query_lists if query_file.find("OV_corner") >= 0][:1]
for query_file in query_lists:
    if query_file.endswith('.mat'):
        # print(query_file)
        dataset_train, dataset_test, corrs_train, corrs_test, poses_train, poses_test, calib_train, calib_test, fov = load_data(base_path + query_file)
        output_path = output_base + query_file.split('.')[0] + '_theta_r.txt'
    
        out = get_intrinsics(dataset_train, corrs_train, poses_train, output_path)
        out_full[query_file] = out

cmd = "../../build/src/get_intrinsics"
os.system(cmd)
for query_file in query_lists:
    if query_file.endswith('.mat'):
        dataset_train, dataset_test, corrs_train, corrs_test, poses_train, poses_test, calib_train, calib_test, fov = load_data(base_path + query_file)
        output_path = output_base + query_file.split('.')[0] + '_theta_r.txt'
    
        # out = get_intrinsics(dataset_train, corrs_train, poses_train, output_path)
        out = out_full[query_file]
        spline_path = spline_base + query_file.split('.')[0] + '_spline.txt'
        if  os.path.exists(spline_path):
            # thetas, theta_r_spline = spline_fit(spline_path)
            thetas, rs, theta_r_spline = spline_fit(output_path)
            mean_reprojection_error, calibrated_points = calculate_mean_reproj_error(corrs_test, poses_test,theta_r_spline, thetas,out)
            breakpoint()
            if query_file.startswith('OV_corner'):
                results['OV_corner'].append([mean_reprojection_error, calibrated_points])
            elif query_file.startswith('OV_cube'):
                results['OV_cube'].append([mean_reprojection_error, calibrated_points])
            elif query_file.startswith('OV_single_plane'):
                results['OV_single_plane'].append([mean_reprojection_error, calibrated_points])
            elif query_file.startswith('Kalibr'):
                results['Kalibr'].append([mean_reprojection_error, calibrated_points])
            elif query_file.startswith('OCamCalib'):
                results['OCamCalib'].append([mean_reprojection_error, calibrated_points])
            elif query_file.startswith('UZH_DAVIS'):
                results['UZH_DAVIS'].append([mean_reprojection_error, calibrated_points])
            elif query_file.startswith('UZH_Snapdragon'):
                results['UZH_Snapdragon'].append([mean_reprojection_error, calibrated_points])

breakpoint()
# calculate the mean reprojection error and calibrated points for each dataset
for key in results:
    mean_reprojection_error = 0
    calibrated_points = 0
    for result in results[key]:
        mean_reprojection_error += result[0]
        calibrated_points += result[1]
    mean_reprojection_error /= len(results[key])
    calibrated_points /= len(results[key])
    # keep 4 decimal places
    mean_reprojection_error = round(mean_reprojection_error, 3)
    calibrated_points = round(calibrated_points * 100, 3)
    print(key, mean_reprojection_error, calibrated_points)


# %%
results

# %%
def add_correspondences_to_database(database_path, corrs, poses, calib):
    db = COLMAPDatabase.connect(database_path)
    db.create_tables()

    # Define camera parameters and add camera to database
    camera_model = 12 # OPENCV_FISHEYE or appropriate model ID
    width, height = calib.get('width', 1280), calib.get('height', 800)
    # params = np.array(calib['params'])
    params = np.array([width, height, 350.000000, 700.000000, 1050.000000, 1400.000000, 1750.000000, 2100.000000, 2450.000000, 2800.000000, 3150.000000, 3500.000000, 3400.000000, 3380.000000, 3360.000000, 3340.000000, 3330.000000, 3310.000000, 3290.000000, 3280.000000, 3280.000000, 3300.000000])
    camera_id = db.add_camera(camera_model, width, height, params)

    # Add images and keypoints to database
    image_ids = []
    for i, (p2d, _, _, _) in enumerate(corrs):
        image_name = f"image_{i}.png"  # Replace with actual image names if available
        image_id = db.add_image(image_name, camera_id)
        image_ids.append(image_id)
        db.add_keypoints(image_id, np.array(p2d, dtype=np.float32))

    # Find matches and add two-view geometries with the correct pair_id
    for i in range(len(corrs)):
        for j in range(i + 1, len(corrs)):
            matches = []

            # Identify common 3D points between images i and j
            _, p3d_i, full_idx_i, _ = corrs[i]
            _, p3d_j, full_idx_j, _ = corrs[j]
            common_indices = {idx: (np.where(full_idx_i == idx)[0][0],
                                    np.where(full_idx_j == idx)[0][0])
                              for idx in set(full_idx_i) & set(full_idx_j)}

            for idx, (pt_i, pt_j) in common_indices.items():
                matches.append([pt_i, pt_j])

            if matches:
                # Convert matches list to numpy array
                matches_array = np.array(matches, dtype=np.uint32)

                # Calculate pair_id for the image pair
                if image_ids[i] < image_ids[j]:
                    # pair_id = image_ids_to_pair_id(image_ids[i], image_ids[j])
                # pair_id = image_ids_to_pair_id(image_ids[i], image_ids[j])

                # Add matches to the database using the computed pair_id
                    db.add_matches(image_ids[i], image_ids[j], matches_array)
                    db.add_two_view_geometry(image_ids[i],image_ids[j], matches_array)
                    print(f"Added {len(matches)} matches between image {i} and image {j}")
        

    db.commit()
    db.close()

# %%
# list all the files in the directory ending with .mat
for query_file in os.listdir(base_path):
    if query_file.endswith('.mat'): 

        print('Processing', query_file)
        corrs, poses, calib, fov = load_data(base_path + query_file)
        # check the name of the file
        calib_dict = calib_to_dict(calib)
        if 'OV_' in query_file:
            calib_dict['width'] = 1280
            calib_dict['height'] = 800
        elif 'Kalibr_BF2M2020S23' in query_file:
            calib_dict['width'] = 1280
            calib_dict['height'] = 1024
        elif 'Kalibr_BM2820' in query_file:
            calib_dict['width'] = 1280
            calib_dict['height'] = 1024
        elif 'Kalibr_BF5M13720' in query_file:
            calib_dict['width'] = 1280
            calib_dict['height'] = 1024
        elif 'Kalibr_BM4018S118' in query_file: 
            calib_dict['width'] = 1280
            calib_dict['height'] = 1024
        elif 'Kalibr_ENTANIYA' in query_file:
            calib_dict['width'] = 1680
            calib_dict['height'] = 1680
        elif 'Kalibr_EUROC' in query_file:
            calib_dict['width'] = 752
            calib_dict['height'] = 480
        elif 'Kalibr_GOPRO' in query_file:
            calib_dict['width'] = 1680
            calib_dict['height'] = 1680
        elif 'Kalibr_TUMVI' in query_file:
            calib_dict['width'] = 512
            calib_dict['height'] = 512
        elif 'UZH_DAVIS' in query_file:
            calib_dict['width'] = 346
            calib_dict['height'] = 260
        elif 'OCamCalib_Fisheye1' in query_file:
            calib_dict['width'] = 1032
            calib_dict['height'] = 778
        elif 'OCamCalib_Fisheye2' in query_file:
            calib_dict['width'] = 748
            calib_dict['height'] = 480
        elif 'OCamCalib_Fisheye190deg' in query_file:
            calib_dict['width'] = 752
            calib_dict['height'] = 480
        elif 'OCamCalib_GOPR' in query_file:
            calib_dict['width'] = 3840
            calib_dict['height'] = 2880
        elif 'OCamCalib_KaidanOmni' in query_file:
            calib_dict['width'] = 640
            calib_dict['height'] = 480
        elif 'OCamCalib_Ladybug' in query_file:
            calib_dict['width'] = 1024
            calib_dict['height'] = 768
        elif 'OCamCalib_MiniOmni' in query_file:
            calib_dict['width'] = 752
            calib_dict['height'] = 480
        elif 'OCamCalib_Omni' in query_file:
            calib_dict['width'] = 1280
            calib_dict['height'] = 960
        elif 'OCamCalib_VMRImage' in query_file:
            calib_dict['width'] = 1024
            calib_dict['height'] = 768
        else:
            print('Unknown image size for', query_file) 
            continue

        # create the database
        output_db_dir = os.path.join(output_db_base, query_file.split('.')[0]) + '.db'
        print(output_db_dir)
        # add_correspondences_to_database(output_db_dir, corrs, poses, calib_dict)
        # save the calibration to the output directory
        output_calib_dir = os.path.join(output_calib_base, query_file.split('.')[0]) + '.txt'
        # write calib_dict['params'] to the output directory
        with open(output_calib_dir, 'w') as f:
            for item in calib_dict['params']:
                f.write("%s\n" % item)
       

# %%
calib_dict

# %%
calib_dict['width'] = 1280
calib_dict['height'] = 800

# %%
# OV: 1280, 800
# Kalibr_BF2M2020S23: 1280, 1024
# Kalibr_BM2820: 1280, 1024
# Kalibr_BF5M13720: 1280, 1024
# Kalibr_BM4018S118: 1280, 1024
# Kalibr_ENTANIYA: 1680, 1680
# Kalibr_EUROC: 752, 480
# Kalibr_GoPro: 1680, 1680
# Kalibr_TUMVI: 512, 512
# UZH_DAVIS: 346, 260
# UZH_Snapdragon:  640, 480
# OCamCalib_Fisheye_1: 1032, 778
# OCamCalib_Fisheye_2: 748, 480
# OCamCalib_Fisheye190deg: 752, 480
# OCamCalib_GOPR: 3840, 2880
# OcamCaliB_KaidanOmni: 640, 480
# OCamCalib_Ladybug: 1024, 768
# OCamCalib_MiniOmni: 752, 480
# OcamCalib_Omni: 1280, 960
# OCamCalib_VMRImage: 1024, 768







# %%
# read a pgm image and determine the image size
import cv2
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
image_path = '/home/yihan/Downloads/OCamCalib/OCamCalib/VMRImage/train/VMRImage0.jpg'
image = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
plt.imshow(image, cmap='gray')
# turn off axis
plt.axis('off')
# size of the image
image.shape


# %%
len(corrs)

# %%


# %%
add_correspondences_to_database('database_cube.db', corrs, poses, calib_dict)


