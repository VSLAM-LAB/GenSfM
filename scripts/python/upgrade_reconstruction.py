import pyimplicitdist
import pycolmap
import numpy as np

# datasets = [
#     "courtyard",
#     "delivery_area",
#     "electro",
#     "facade",
#     "kicker",
#     "meadow",
#     "office",
#     "pipes",
#     "playground",
#     "relief",
#     "relief_2",
#     "terrace",
#     "terrains"
# ]

# base = "dslr"
# mode = "training"

# # TODO: Change the directory here to mactch the data
# input_dir = "sparse/0"
# output_dir = "test"

# for dataset in datasets:
# TODO: Also change the directory here
# dir_base = "../../../datasets/ETH3D/mvs/multi_view_" + mode + "_" + base + "_jpg/" + dataset + "/"

# First, read in the reconstruction.
# Here, we assume that the reconstruction has been converted into standard COLMAP format.
reconstruction = pycolmap.Reconstruction()
reconstruction.read('/home/yihan/cvg/datasets/cata_mainbuilding/1120_radial/txt_2')
output_dir = '/home/yihan/cvg/datasets/cata_mainbuilding/1120_radial/upgraded/' 

cm_opt = pyimplicitdist.CostMatrixOptions()
refinement_opt = pyimplicitdist.PoseRefinementOptions()
refinement_opt.verbose = True

ba_opt = pyimplicitdist.BundleAdjustmentOptions()
ba_opt.upgrade_result = True
ba_opt.min_curr_num = 3
ba_opt.filter_thres = 50


cam_id_to_image_ids = {}
for image_id in reconstruction.reg_image_ids():
    image = reconstruction.images[image_id]
    cam_id = image.camera_id
    if cam_id not in cam_id_to_image_ids:
        cam_id_to_image_ids[cam_id] = []
    cam_id_to_image_ids[cam_id].append(image_id)
    

# Then, prepare the structure for the camera upgrade
for camera in reconstruction.cameras.values():
    p2d_multi = []
    p3d_multi = []
    
    poses_multi = []
    pts_index_multi = []
    p3d_all = []
    
    p3d_id_to_index = {}
    p3d_index_to_id = {}
    for i, image_id in enumerate(cam_id_to_image_ids[camera.camera_id]):
        
        image = reconstruction.images[image_id]
        p2d = []
        pts_index = [] 
        points2D = image.points2D
        for point2D in points2D:
            if point2D.has_point3D():
                point3D = reconstruction.points3D[point2D.point3D_id]
                p2d.append(point2D.xy)
                
                if point2D.point3D_id not in p3d_id_to_index:
                    p3d_id_to_index[point2D.point3D_id] = len(p3d_all)
                    p3d_index_to_id[len(p3d_all)] = point2D.point3D_id
                    p3d_all.append(point3D.xyz)
                    
                pts_index.append(p3d_id_to_index[point2D.point3D_id])
        

        p2d_multi.append(p2d)
        pts_index_multi.append(pts_index)
        
        initial_pose = pyimplicitdist.CameraPose()
        initial_pose.q_vec = image.cam_from_world.rotation.quat[[3, 0, 1, 2]]
        initial_pose.t = image.cam_from_world.translation
        
        poses_multi.append(initial_pose)
    
    
    # Perform pose upgrade
    pp = np.array([camera.principal_point_x, camera.principal_point_y])
    cost_matrix = pyimplicitdist.build_cost_matrix_multi(p2d_multi, cm_opt, pp)
    
    print("pose refinement done")
    
    
    # Perform implicit bundle adjustment
    # Curently, camera pose are first upgraded, then the points3D are upgraded
    out = pyimplicitdist.bundle_adjustment(p2d_multi, p3d_all, pts_index_multi, cost_matrix, cm_opt, poses_multi, pp, ba_opt)
    
    poses_new = out["poses"]
    point3D_new = out["points3D"]
    
    # Write the new camera poses, points3d to the reconstruction
    for i, image_id in enumerate(cam_id_to_image_ids[camera.camera_id]):
        image = reconstruction.images[image_id]
        image.cam_from_world.rotation.quat = poses_new[i].q_vec[[1, 2, 3, 0]]
        image.cam_from_world.translation = poses_new[i].t
        
    for j in range(len(p3d_all)):
        reconstruction.points3D[p3d_index_to_id[j]].xyz = point3D_new[j]
    
    
    reconstruction.write(output_dir)
    
    






