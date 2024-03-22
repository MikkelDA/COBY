### General tools
def flatten(matrix):
    ### https://realpython.com/python-flatten-list/
    ### Fastest one shown
    flat_list = []
    for row in matrix:
        flat_list += row
    return flat_list

