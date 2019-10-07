"""
uses scikit-bio and pandas


"""


import skbio
assert skbio.__version__>= '0.5.1' , "Need scikit-bio >= 0.5.1"

import pandas as pd
from skbio.stats import composition

def min_filer(D,min_count=1,min_N=2):

    """

        Filter matrix or dataset (samples x features) by filter criteria

        if min_N is a fraction, the fraction of the nuber of samples is taken

    """


    if (0< min_N< 1):
        min_N= np.round(D.shape[0]*min_N)

    subset= D.loc[:,(D>=min_count).sum()>=min_N]

    print(f"Remove {D.shape[1]-subset.shape[1]} ({subset.shape[1]/ D.shape[1]})"
          f" features because below the filter criteria (> {min_count} in at least {min_N} samples) ")

    return subset


def normalize_clr(data):

    "replace zeros and apply clr"

    assert data.shape[0]< data.shape[1], "samples should be indexes, I don't think you have"

    normalized=composition.clr(composition.multiplicative_replacement(data))
    normalized= pd.DataFrame(normalized,
                             index= data.index,columns= data.columns)

    return normalized

if __name__ == '__main__':

    Nmetadata_columns=2
    excel_file="/Volumes/m-phyme/GTrajkovs/Database/E3&E5_template.xlsx"
    sheet_name="16S"
    new_sheet_name="16S_normalized"
    D= pd.read_excel(excel_file,sheet_name=sheet_name,index_col=tuple(range(Nmetadata_columns)))

    #transpose to get samples x features
    D= D.T
    normalized = normalize_clr(D)
    #transpose back
    normalized= normalized.T

    normalized.to_excel(excel_file,sheet_name=new_sheet_name)
