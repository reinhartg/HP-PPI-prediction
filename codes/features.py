# Functions to extract protein features for learning-based HP-PPI predictions
from scipy.sparse import csr_matrix

def domain_features(protein_pair, domain_dict, domain_set):
    '''Generate a feature vector of protein domain occurences in a protein.
    
    Parameters
    ----------
    protein_pair : iterable
        A pair of protein identifiers (must be contained in the keys of `domain_dict`)
    
    domain_dict : dict
        Mapping of protein identifiers to a list of its domains
    
    domain_set : set, list, or list-like objects
        Set of unique domains contained in datasets
    
    Return
    ------
    features : list, shape=[len(pfam_set)]
        Feature vector of domain occurences in `protein` (1 for presence and 0 for absence)
    '''
    
    features = []
    for domain in domain_set:
        
        # Extract domains from pair
        domains = domain_dict[protein_pair[0]] + domain_dict[protein_pair[1]]
        
        # Check existence of a domain in the pair
        if domain in domains:
            features.append(1)
        else:
            features.append(0)
    
    return csr_matrix(features)