# Functions to extract protein features for learning-based HP-PPI predictions

from scipy import sparse

def domain_features(protein, domain_dict, domain_set):
    '''Generate a feature vector of protein domain occurences in a protein.
    
    Parameters
    ----------
    protein : string
        Protein identifier (must be contained in the keys of `domain_dict`)
    
    domain_dict : dict
        Mapping of protein identifiers to a list of its domains
    
    domain_set : set, list, or list-like objects
        Set of unique domains contained in datasets
    
    Return
    ------
    features : 1xlen(domain_set) sparse matrix
        Feature vector of domain occurences in `protein` (1 for presence and 0 for absence) represented as a sparse matrix
    '''
    
    domain_profile = []
    for domain in domain_set:
        if domain in domain_dict[protein]:
            domain_profile.append(1)
        else:
            domain_profile.append(0)
    
    features = sparse.csr_matrix(domain_profile)
    return features