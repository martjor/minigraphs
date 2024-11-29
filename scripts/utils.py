from scipy.sparse import load_npz, save_npz
from networkx import from_scipy_sparse_array, to_scipy_sparse_array
from yaml import dump, safe_load

def load_graph(file):
    return from_scipy_sparse_array(load_npz(file))

def save_graph(file,graph):
    save_npz(file,to_scipy_sparse_array(graph))
    
def load_dict(file):
    with open(file,'r') as f: 
        return safe_load(f)

def save_dict(file,dict):
    with open(file,'w') as f:
        dump(dict,f)
        

    
