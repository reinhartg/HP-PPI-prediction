# Functions to extract protein features for learning-based HP-PPI predictions

import numpy as np

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
    features : array, shape=[len(domain_set)]
        Feature vector of domain occurences in `protein` (1 for presence and 0 for absence)
    '''
    
    domain_counts = []
    for domain in domain_set:
        c = domain_dict[protein].count(domain)
        domain_counts.append(c)
    
    features = np.array(domain_counts)
    return features