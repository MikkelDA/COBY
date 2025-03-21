class get_number_from_string:
    def get_number_from_string(self, s):
        isnumber, isint = self.is_number(s)
        if isnumber and isint:
            return int(s)
        elif isnumber:
            return float(s)
        else:
            return False
