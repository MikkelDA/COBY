def center_coords(self, coords):
    '''
    Calculate lipid-spicific x/y-center
    '''
    coord_diff = (max(coords) + min(coords)) / 2
    centered_coords = [coord - coord_diff for coord in coords]
    return centered_coords

