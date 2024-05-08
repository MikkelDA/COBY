class is_number:
    def is_number(self, s):
        
        if type(s) == str:
            if s.startswith("-"):
                s = s[1:]
        else:
            s = str(s)
        
        try:
            test = float(s) # Errors if not floatable and returns exception
            number = True
            if s.isdigit():
                integer = True
            else:
                integer = False
            return number, integer
        
        except ValueError:
            return False, False
    
