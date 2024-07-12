### Recursive prints dictionary keys and values
### Appropriate number of whitespaces added to front by recursive calls to the function
### Only used for debugging during development
def dictprinter(d, r = 0):
    if type(d) == dict:
        for key, vals in d.items():
            print(r*"    " + str(key))
            dictprinter(vals, r+1)
    elif type(d) in [str, int, float, list, tuple]:
        print(r*"    " + str(d))
