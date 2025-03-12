class get_geoms_list:
    def get_geoms_list(self, geoms):
        '''
        Takes a geometry object from shapely and returns a list of all geometries/subgeometries
        '''
        geoms_list = []
        if geoms.geom_type in ["Point", "LineString", "LinearRing", "Polygon"]:
            self.print_term("(get_geoms_list)", "Geometry", geoms.geom_type, debug = True)
            geoms_list.extend([geoms])
        elif geoms.geom_type in ["MultiPoint", "MultiLineString", "MultiPolygon"]:
            self.print_term("(get_geoms_list)", "Geometry", geoms.geom_type, debug = True)
            geoms_list.extend(list(geoms.geoms))
        elif geoms.geom_type == "GeometryCollection":
            ### Recursively check GeometryCollections as they can contain "Multi"-geometries
            self.print_term("(get_geoms_list)", "Geometry", geoms.geom_type, debug = True)
            for sub_geoms in list(geoms.geoms):
                geoms_list.extend(list(geoms.geoms))
        elif type(geoms) in [list, tuple]:
            ### Recursively check nested geometries inside lists and tuples
            self.print_term("(get_geoms_list)", "List/tuple of geometries", debug = True)
            for sub_geoms in geoms:
                geoms_list.extend(self.get_geoms_list(sub_geoms))
        elif type(geoms) == dict:
            ### Recursively check nested geometries inside dictionaries
            self.print_term("(get_geoms_list)", "Dict of geometries", debug = True)
            for sub_geoms in geoms.values():
                geoms_list.extend(self.get_geoms_list(sub_geoms))
        geoms_list = [geom for geom in geoms_list if not geom.is_empty]
        return geoms_list
    