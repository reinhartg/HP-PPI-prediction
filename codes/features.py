# Functions to extract protein features for learning-based HP-PPI predictions
from scipy.sparse import csr_matrix

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
    features : sparse matrix, shape=[len(pfam_set)]
        Feature vector of domain occurences in `protein` (1 for presence and 0 for absence)
    '''
    
    domain_profile = []
    domain_counts = domain_dict[protein] # domain counts of the input protein
    
    for domain in domain_set:
        
        # Check existence of each domain in the protein
        if domain in domain_counts:
            domain_profile.append(domain_counts[domain])
        else:
            domain_profile.append(0)
    
    # Get features as a sparse matrix
    features = csr_matrix(domain_profile)
    return features