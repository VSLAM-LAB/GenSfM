import sqlite3
import numpy as np


MAX_IMAGE_ID = 2**31 - 1

def image_ids_to_pair_id(image_id1, image_id2):
    if image_id1 > image_id2:
        image_id1, image_id2 = image_id2, image_id1
    return image_id1 * MAX_IMAGE_ID + image_id2


def pair_id_to_image_ids(pair_id):
    image_id2 = int(pair_id % MAX_IMAGE_ID)
    image_id1 = int((pair_id - image_id2) / MAX_IMAGE_ID)
    return image_id1, image_id2


def array_to_blob(array):
    return array.tostring()


def blob_to_array(blob, dtype, shape=(-1,)):
    return np.fromstring(blob, dtype=dtype).reshape(*shape)


def read_matches(matches_path):
    matches = {}
    points2D = np.zeros((0,2))
    points3D = np.zeros((0,3))

    with open(matches_path) as f:
        lines = f.read().splitlines()

        for line in lines:
            row = line.split(' ')
            points2D = np.r_[ points2D, np.expand_dims(np.array(row[0:2],dtype=float),axis=0)]
            points3D = np.r_[ points3D, np.expand_dims(np.array(row[2:5],dtype=float),axis=0)]
            
    return points2D, points3D
        

def read_matches_list(dataset_path, max_queries=None, image_list='image_list.txt'):
    queries = []
    with open(dataset_path + '/' + image_list) as f:
        lines = f.read().splitlines()
        for line in lines:
            row = line.split(' ')      
            query = {}  
            query['name'] = row[0]
            query['camera'] = {}
            query['camera']['model'] = row[1]
            query['camera']['width'] = int(row[2])
            query['camera']['height'] = int(row[3])
            n_params = len(row) - 11
            query['camera']['params'] = np.array(row[4:4+n_params],dtype=float) 

            query['qvec'] = np.array(row[4+n_params:8+n_params],dtype=float)
            query['tvec'] = np.array(row[8+n_params:],dtype=float)

            p2d, p3d = read_matches(dataset_path + '/' + query['name'] + '.matches.txt')
            query['matches'] = (p2d,p3d)

            queries.append(query)

            if max_queries is not None:
                if len(queries) >= max_queries:
                    break
    return queries
