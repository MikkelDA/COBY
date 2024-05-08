import numpy as np

class rotation_matrix_from_vectors:
    def rotation_matrix_from_vectors(self, vector1, vector2):
        """
        Math explananation:
        https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/476311#476311
        Code explanation:
        https://stackoverflow.com/questions/45142959/calculate-rotation-matrix-to-align-two-vectors-in-3d-space
        Find the rotation matrix that aligns vector1 to vector2
        :param vector1: A 3d "source" vector
        :param vector2: A 3d "destination" vector
        :return rotation_matrix: A transform matrix (3x3) which when applied to vector1, aligns it with vector2.
        """
        a = (vector1 / np.linalg.norm(vector1)).reshape(3)
        b = (vector2 / np.linalg.norm(vector2)).reshape(3)
        v = np.cross(a, b) # cross product
        if any(v): # If not all zeros then calculate rotation matrix
            ### rotation_matrix = I + [v]_x + [v]_x^2 * ((1-c)/s^2)
            c = np.dot(a, b) # cosine of angle
            ### skew-symmetric cross-product (SSCP) matrix of v
            SSCP = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
            ### (1-c) / s^2 = (1 - c) / (1 - c^2) = 1 / (1 + c)
            last_bit = 1 / (1 + c)
            rotation_matrix = np.eye(3) + SSCP + SSCP.dot(SSCP) * last_bit
            return rotation_matrix
        else: # Rotation matrix is simply the identity matrix
            return numpy.eye(3)
