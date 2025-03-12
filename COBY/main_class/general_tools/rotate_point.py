import math

class rotate_point:
    def rotate_point(self, posx, posy, posz, x_ang_deg, y_ang_deg, z_ang_deg):
        '''
        Rotates a point in 3D space using its original position and angles in degrees.
        '''
        ### Convert angles to radians
        x_ang = math.radians(x_ang_deg)
        y_ang = math.radians(y_ang_deg)
        z_ang = math.radians(z_ang_deg)

        ### Calculate rotation matrices
        rx = [[1               , 0               , 0               ],
              [0               , math.cos(x_ang) , -math.sin(x_ang)],
              [0               , math.sin(x_ang) , math.cos(x_ang) ]]
        
        ry = [[math.cos(y_ang) , 0               , math.sin(y_ang) ],
              [0               , 1               , 0               ],
              [-math.sin(y_ang), 0               , math.cos(y_ang) ]]
        
        rz = [[math.cos(z_ang) , -math.sin(z_ang), 0               ],
              [math.sin(z_ang) , math.cos(z_ang) , 0               ],
              [0               , 0               , 1               ]]
        
#         rx = [[1, 0, 0]                             , [0, math.cos(x_ang) , -math.sin(x_ang)], [0, math.sin(x_ang), math.cos(x_ang)]]
#         ry = [[math.cos(y_ang), 0, math.sin(y_ang)] , [0, 1, 0]                              , [-math.sin(y_ang)  , 0, math.cos(y_ang)]]
#         rz = [[math.cos(z_ang), -math.sin(z_ang), 0], [math.sin(z_ang) , math.cos(z_ang), 0] , [0, 0, 1]]

        ### Combine rotation matrices
        rxy  = [[sum(a * b for a, b in zip(row, col)) for col in zip(*ry)] for row in rx]
        rxyz = [[sum(a * b for a, b in zip(row, col)) for col in zip(*rz)] for row in rxy]

        ### Calculate new position of point
        new_posx, new_posy, new_posz = [rxyz[ax][0] * posx + rxyz[ax][1] * posy + rxyz[ax][2] * posz for ax in [0, 1, 2]]

        return new_posx, new_posy, new_posz

