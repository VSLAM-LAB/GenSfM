import numpy as np
import sqlite3
import os
import sys

# CAMERA_MODEL = 1
CAMERA_MODEL = 12


def filter_duplicates(all_keyp, all_desc):
    for im_k in range(0,8):
        if(im_k == 7):
            print("Filtering pair ", 7, 0)
            bad_ind = filter_pair(np.asarray(all_keyp[7]), np.asarray(all_keyp[0]))
        else:
            print("Filtering pair ", im_k, im_k+1)
            bad_ind = filter_pair(np.asarray(all_keyp[im_k]), np.asarray(all_keyp[im_k+1]))

        print(all_keyp[im_k].shape)
        mask = np.ones(all_keyp[im_k].shape[0], dtype=bool)
        mask[bad_ind] = False
        all_keyp[im_k] = np.asarray(all_keyp[im_k])[mask]
        all_desc[im_k] = np.asarray(all_desc[im_k])[mask]
        print(all_keyp[im_k].shape)

    return all_keyp, all_desc
    
def filter_pair(keyp1,keyp2):
    
    # Calculate pairwise distances between all points in keyp1 and keyp2
    diff = keyp1[:, None, 0:2] - keyp2[None, :, 0:2]
    dist = np.linalg.norm(diff, axis=2)

    # Check for any distances below the threshold
    min_dist = np.min(dist, axis=1)

    # Find indices where the distance is less than 2.0
    bad_ind = np.where(min_dist < 2.0)[0]

    return bad_ind.tolist()
    



if len(sys.argv) < 3:
    print('Needs three arguments!')
    exit()

db_path = sys.argv[1]
if db_path[-3:] == ".db":
    filename_split_db = db_path
else:
    filename_split_db = db_path + '/database_split.db'

print("Opening database: " + filename_split_db)
if not os.path.exists(filename_split_db):
    print('Error db does not exist!')
    exit()

connection_split = sqlite3.connect(filename_split_db)
cursor_split = connection_split.cursor()

db_path = sys.argv[2]
if db_path[-3:] == ".db":
    filename_db = db_path
else:
    filename_db = db_path + '/database.db'

print("Opening database: " + filename_db)
if not os.path.exists(filename_db):
    print('Error db does not exist!')
    exit()

connection = sqlite3.connect(filename_db)
cursor = connection.cursor()

params = sys.argv[3].split(',')
if len(params) != 7:
    print('Third argument must be list of 7 parameters (width,height,pp(1),pp(2),r_inner,r_outer,k)')
    exit()

width, height, ppx, ppy, r_inner, r_outer, k_dist = params
r_inner = float(r_inner)
r_outer = float(r_outer)
k_dist = float(k_dist)
ppx = float(ppx)
ppy = float(ppy)

print("width=",width)
print("height=",height)
print("pp=",ppx,ppy)
print("radius=",r_inner,r_outer)
print("k_dist=",k_dist)
pp = np.asarray([ppx,ppy], np.float64)

# Make sure new database is clean
cursor.execute("DELETE FROM cameras;")
cursor.execute("DELETE FROM images;")
cursor.execute("DELETE FROM keypoints;")
cursor.execute("DELETE FROM descriptors;")


# Find original images names (before split)
im_names = []
print("Listing images in database...\n")
cursor_split.execute('SELECT image_id, name FROM images WHERE name LIKE "%_0.png";')
for row in cursor_split:
    image_idx, name = row        
    name = name[0:-6]
    print("Image: " + name + ", id = " + str(image_idx))
    im_names.append(name)
    
# Insert camera

print("Creating camera...")
if CAMERA_MODEL == 1:
    pp = np.asarray([1,1,ppx,ppy], np.float64)
else:
    pp = np.asarray([ppx,ppy], np.float64)
print(pp)
cursor.execute('INSERT INTO cameras VALUES (?,?,?,?,?,?)', (1, CAMERA_MODEL, width, height, pp.tobytes(), 0))

image_id = 1

for im_name in im_names:
    print("Processing: ",im_name)
    all_keyp = np.zeros((0,6),dtype=np.float32)
    all_desc = np.zeros((0,128),dtype=np.uint8)

    print("Adding image to new database.")
    cursor.execute("INSERT INTO images VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (image_id, im_name, 1, None, None, None, None, None, None, None))
    all_keyp = []
    all_desc = []
    for im_k in range(0,8):
        # Fetch keypoints
        cursor_split.execute("SELECT rows, cols, data FROM keypoints LEFT JOIN images ON keypoints.image_id = images.image_id WHERE images.name = '" + im_name + "_" + str(im_k)  + ".png';")
        n_rows, n_columns, raw_data = cursor_split.fetchone()
        keyp = np.asarray(np.frombuffer(raw_data, dtype=np.float32).reshape(n_rows, n_columns)).copy()
        print("Read ",n_rows, n_columns)

        # Transform keypoint back into original image
        r = keyp[:,1] + r_inner
        phi = keyp[:,0]/k_dist + np.pi/4 * im_k

        keyp[:,0] = ppx + r * np.cos(phi)
        keyp[:,1] = ppy + r * np.sin(phi)
        
        # Fetch descriptors
        cursor_split.execute("SELECT rows, cols, data FROM descriptors LEFT JOIN images ON descriptors.image_id = images.image_id WHERE images.name = '" + im_name + "_" + str(im_k)  + ".png';")
        n_rows, n_columns, raw_data = cursor_split.fetchone()
        desc = np.asarray(np.frombuffer(raw_data, dtype=np.uint8).reshape(n_rows, n_columns)).copy()

        all_keyp.append(keyp)
        all_desc.append(desc)

    all_keyp, all_desc = filter_duplicates(all_keyp, all_desc)

    all_keyp = np.concatenate(all_keyp, axis=0)
    all_desc = np.concatenate(all_desc, axis=0)
    
    print("In total:", all_keyp.shape, all_desc.shape)
    # all_keyp = np.asnumpy(all_keyp)
    all_keyp = np.asarray(all_keyp)
    #
    # all_desc = np.asnumpy(all_desc)
    all_desc = np.asarray(all_desc)
    cursor.execute("INSERT INTO keypoints VALUES (?, ?, ?, ?)",
            (image_id,) + all_keyp.shape + (all_keyp.tostring(),))

    cursor.execute("INSERT INTO descriptors VALUES (?, ?, ?, ?)",
            (image_id,) + all_desc.shape + (all_desc.tostring(),))

    
    cursor.execute("SELECT image_id,rows,cols FROM descriptors")
    for row in cursor:
        print(row)

    image_id = image_id + 1

    
    
connection.commit()
cursor.close()
connection.close()


cursor_split.close()
connection_split.close()
