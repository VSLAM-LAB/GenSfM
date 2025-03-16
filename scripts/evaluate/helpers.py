import numpy as np
import math
import poselib
from scipy.spatial.transform import Rotation

M_PI = 3.14159265359


def calc_angle(pose1, pose2):
    cos_r = ((pose1.R.T @ pose2.R).trace() - 1) / 2
    cos_r = min(max(cos_r, -1.), 1.)

    return math.acos(cos_r) * 180 / M_PI

def calc_trans(pose1, pose2):
    return np.linalg.norm((pose1.R.T @ pose1.t - pose2.R.T @ pose2.t))

def calc_trans_angle(pose1, pose2):
    # If the translation is too small, we say the error is zero
    if np.linalg.norm(pose1.t) < 1e-10 or np.linalg.norm(pose2.t) < 1e-10:
        return 0
    cos_r = (pose1.t.dot(pose2.t)) / (np.linalg.norm(pose1.t) * np.linalg.norm(pose2.t))
    cos_r = min(max(cos_r, -1.), 1.)
    return math.acos(cos_r) * 180 / M_PI

def calc_relative_pose(pose1, pose2):
    R =  pose2.R @ pose1.R.T
    quat = Rotation.from_matrix(R).as_quat()[[3, 0, 1, 2]]
    t = -R @ pose1.t + pose2.t
    pose = poselib.CameraPose()
    pose.q = quat
    pose.t = t

    return pose

# Return errors_r, error_t
# R = R_j * R_i.T
# t = R_j * (c_i - c_j)
def calc_pairwise_relative_pose_batch(poses):
    R_batch = np.zeros([len(poses) * 3, 3])
    c_batch = np.zeros([len(poses), 3])

    for i in range(len(poses)):
        if poses[i] is None:
            R_batch[i * 3:i * 3 + 3, :] = np.eye(3)
            c_batch[i] = np.zeros([3, ])
        else:
            R_batch[i * 3:i * 3 + 3, :] = poses[i].R
            c_batch[i] = poses[i].center()
    
    # Pre allocate the memory
    R_rel = np.zeros([3, len(poses) * len(poses) * 3])
    t_rel = np.zeros([3, len(poses) * len(poses)])
    for j in range(len(poses)):
        R_rel[:, len(poses) * 3 * j : len(poses) * 3 * (j + 1)] = R_batch[j * 3:j * 3 + 3, :] @ R_batch.T
        t_rel[:, len(poses) * j : len(poses) * (j + 1)] = R_batch[j * 3:j * 3 + 3, :] @ (c_batch - c_batch[j]).T
    
    # Perform column-wise normalization
    t_rel[:, np.linalg.norm(t_rel, axis=0) < 1e-10] = 0
    t_rel = t_rel / (1e-10 + np.linalg.norm(t_rel, axis=0))

    return R_rel, t_rel

def calc_pairwise_relative_error_batch(poses_calc, poses_gt):
    R_rel_calc, t_rel_calc = calc_pairwise_relative_pose_batch(poses_calc)
    R_rel_gt, t_rel_gt = calc_pairwise_relative_pose_batch(poses_gt)
    R_diff = R_rel_calc - R_rel_gt
    t_diff = t_rel_calc - t_rel_gt


    R_error = np.zeros([R_diff.shape[1] // 3, ])
    for i in range(3):
        for j in range(3):
            R_error += R_diff[i, j::3] * R_diff[i, j::3]
    R_error = 1 - R_error / 4
    R_error = np.arccos(np.clip(R_error, -1, 1)) * 180 / M_PI

    t_error = np.sum(t_diff * t_diff, axis=0)
    t_error = 1 - t_error / 2
    t_error = np.arccos(np.clip(t_error, -1, 1)) * 180 / M_PI

    for i in range(len(poses_gt)):
        if poses_calc[i] is None:
            t_error[i * len(poses_calc) : (i + 1) * len(poses_calc)] = 180
            R_error[i * len(poses_calc) : (i + 1) * len(poses_calc)] = 180
            t_error[i::len(poses_calc)] = 180
            R_error[i::len(poses_calc)] = 180
        if poses_gt[i] is None:
            t_error[i * len(poses_gt) : (i + 1) * len(poses_gt)] = 180
            R_error[i * len(poses_gt) : (i + 1) * len(poses_gt)] = 180
            t_error[i::len(poses_gt)] = 180
            R_error[i::len(poses_gt)] = 180

    return R_error, t_error


def establish_accum(data, times, bins=20):
    count, bins_count = np.histogram(data, bins=bins, range=[0, times])
    
    # finding the PDF of the histogram using count values
    pdf = count / data.shape[0]
    
    # using numpy np.cumsum to calculate the CDF
    # We can also find using the PDF values by looping and adding
    cdf = np.cumsum(pdf)
        
    return cdf

def compute_recall(errors):
    num_elements = len(errors)
    sort_idx = np.argsort(errors)
    errors = np.array(errors.copy())[sort_idx]
    recall = (np.arange(num_elements) + 1) / num_elements
    return errors, recall


def compute_auc(errors, thresholds):
    errors, recall = compute_recall(errors)

    recall = np.r_[0, recall]
    errors = np.r_[0, errors]

    aucs = []
    for t in thresholds:
        last_index = np.searchsorted(errors, t, side="right")
        r = np.r_[recall[:last_index], recall[last_index-1]]
        e = np.r_[errors[:last_index], t]
        auc = np.trapz(r, x=e)/t
        aucs.append(auc*100)
    return aucs

def parse_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    image_data = []
    for line in lines:
        if line.startswith("#") or line.strip() == "":
            # Skip comments and empty lines
            continue

        # Split the line into components
        parts = line.strip().split(" ")
        if len(parts) != 10:
            continue

        # Extract the values
        image_id = int(parts[0])
        qw, qx, qy, qz = map(float, parts[1:5])
        tx, ty, tz = map(float, parts[5:8])
        camera_id = int(parts[8])
        image_name = parts[9]
        # only keep the image_name after the last '/'
        image_name = image_name.split('/')[-1]

        # Calculate the camera center using quaternion and translation
        # Camera center C = -R^T * T, where R is the rotation matrix from quaternion
        q = np.array([qw, qx, qy, qz])
        t = np.array([tx, ty, tz])

        # Compute the rotation matrix from the quaternion
        R = quaternion_to_rotation_matrix(q)
        camera_center = -np.dot(R.T, t)

        image_data.append((image_name, camera_center))

    return image_data

def quaternion_to_rotation_matrix(q):
    # Normalize the quaternion
    q = q / np.linalg.norm(q)
    qw, qx, qy, qz = q

    # Construct the rotation matrix
    R = np.array([
        [1 - 2*qy**2 - 2*qz**2, 2*qx*qy - 2*qz*qw, 2*qx*qz + 2*qy*qw],
        [2*qx*qy + 2*qz*qw, 1 - 2*qx**2 - 2*qz**2, 2*qy*qz - 2*qx*qw],
        [2*qx*qz - 2*qy*qw, 2*qy*qz + 2*qx*qw, 1 - 2*qx**2 - 2*qy**2]
    ])
    return R

def parse_file_pose(file_path):
     with open(file_path, 'r') as file:
        lines = file.readlines()

        image_data = []
        for line in lines:
            if line.startswith("#") or line.strip() == "":
                # Skip comments and empty lines
                continue

            # Split the line into components
            parts = line.strip().split(" ")
            if len(parts) != 10:
                continue

            # Extract the values
            image_id = int(parts[0])
            qw, qx, qy, qz = map(float, parts[1:5])
            tx, ty, tz = map(float, parts[5:8])
            camera_id = int(parts[8])
            image_name = parts[9]
            # only keep the image_name after the last '/'
            image_name = image_name.split('/')[-1]

            image_data.append((image_name, (qw, qx, qy, qz, tx, ty, tz)))

        return image_data