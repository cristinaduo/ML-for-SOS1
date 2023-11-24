# Further Optimisation: Multiprocessing might be suitable for accelerating the execution of the 2 fuctions defined in this module.

import re

import numpy   as np
import pandas  as pd
import utilies as U

from rdkit      import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import DataStructs

def SMILES_DataSet(dataset_df = None, smiles_list = [], smiles_loc = None, fp_radius = 3, fp_bits = 512):

    '''
    Use this function to generate the dataframe of fingerprint
    dataset_df: the input dataset should be a dataframe
    smiles_list: Or directly input in the form of series without inputing dataframe and column name.
    smiles_loc: the column name that consists of SMILEs strings
    fp_radius = the radius of Morgan fingerprint
    fp_bits = the number of fingerprint bits of Morgan fingerprint
    '''
    
    if len(smiles_list) == 0:

        smiles = dataset_df[smiles_loc]
        smiles_list = np.array(smiles).tolist()

    else:
        
        smiles_list = smiles_list
    
    mols = []
    Error_01 = []
    
    counter = 0
    BackTracker = []
    
    for counter in range (0, len(smiles_list)):

        try:
            mols.append(Chem.MolFromSmiles(smiles_list[counter]))
            BackTracker.append(counter)
        except:
            Error_01.append([smiles_list[counter], counter])
    
# Note that adding Hs would change the fingerprint albeit the same molecule

    Mols = []
    Error_02 = []
    counter = 0
    
    for counter in range (0, len(mols)):
        try:
            Mols.append(AllChem.AddHs(mols[counter]))
        except:
            Error_02.append([smiles_list[BackTracker[counter]], BackTracker[counter]])

    morgans = [AllChem.GetMorganFingerprintAsBitVect(Mol, radius = fp_radius,
                nBits= fp_bits, useChirality = True) for Mol in Mols]
    morgan_bits =  [morgan.ToBitString() for morgan in morgans] 

# find every single digit
    
    pattern = re.compile('.{1}') 
    morgan_bits = [','.join(pattern.findall(morgan)) for morgan in morgan_bits]

    fp_list = []
    for bit in morgan_bits:
        single_fp = bit.split(',')   # split the string by commas
        single_fp = [float(fp) for fp in single_fp] # transfer string to float32
        fp_list.append(single_fp)

    fp_df = pd.DataFrame(np.array(fp_list))
    fp_df.columns = fp_df.columns.astype(str)

# rename the columns

    for i in range(fp_df.columns.shape[0]):
        fp_df.rename(columns = {fp_df.columns[i]: fp_df.columns[i] + "pChEMBL"}, inplace = True)

    return fp_df, Error_01, Error_02

def SMILES_DataSet_Lightway(use_mols = False, mols_list = [], dataset_df = None, smiles_list = [], smiles_loc = None, division = 1, file = None, fp_radius = 0, fp_bits = 0):

    '''
    Use this function to generate the dataframe of fingerprint
    dataset_df: the input dataset should be a dataframe
    smiles_list: Or directly input in the form of series without inputing dataframe and column name.
    smiles_loc: the column name that consists of SMILEs strings
    division: separation rate of the whole computation to save RAM 
    fp_radius: the radius of Morgan fingerprint
    fp_bits: the number of fingerprint bits of Morgan fingerprint
    '''
    
    if use_mols == False:
        if len(smiles_list) == 0:
            smiles = dataset_df[smiles_loc]
            smiles_list = np.array(smiles).tolist()
        else:
            smiles_list = smiles_list
        unit = int(len(smiles_list) / division)
    else:
        unit = int(len(mols_list) / division)
    
    D = 0
    After = 0
    Before = 0
    
    Error_01 = []
    Error_02 = []
    BackTracker = []
    
    for D in range (0, (division - 1)):
    
        mols = []
    
        for counter_01 in range ((0 + unit * D), (unit + unit * D)):

            try:
                if use_mols == False:
                    mols.append(Chem.MolFromSmiles(smiles_list[counter_01]))
                else:
                    mols.append(mols_list[counter_01])
                BackTracker.append(counter_01)
            except:
                Error_01.append([mols_list[counter_01], counter_01])
    
# Note that adding Hs would change the fingerprint albeit the same molecule.

        Mols = []
        Before = After
    
        for counter_02 in range (0, len(mols)):
            try:
                Mols.append(AllChem.AddHs(mols[counter_02]))
            except:
                Error_02.append([mols_list[BackTracker[Before + counter_02]], BackTracker[Before + counter_02]])
            After = Before + counter_02

        morgans = [AllChem.GetMorganFingerprintAsBitVect(Mol, radius = fp_radius, nBits= fp_bits, useChirality = True) for Mol in Mols]
        morgan_bits =  [morgan.ToBitString() for morgan in morgans] 

# find every single digit
    
        pattern = re.compile('.{1}') 
        morgan_bits = [','.join(pattern.findall(morgan)) for morgan in morgan_bits]

        fp_list = []
        for bit in morgan_bits:
            single_fp = bit.split(',')   
            single_fp = [float(fp) for fp in single_fp]
            fp_list.append(single_fp)

        fp_df = pd.DataFrame(np.array(fp_list))
        fp_df.columns = fp_df.columns.astype(str)

# rename the columns

        for i in range(fp_df.columns.shape[0]):
            fp_df.rename(columns = {fp_df.columns[i]:fp_df.columns[i] + "pChEMBL"}, inplace = True)
            
        U.save_dataset(fp_df, path = ('./datasets/' + file + '/'), file_name = ('FP_3_512_D' + str(D)), idx = False)
        
        print('\nDivision ', D, ' Completed.')
        
# last iteration

    for D in range (1):
    
        mols = []
    
        for counter_01 in range ((unit * (division - 1)), len(smiles_list)):

            try:
                mols.append(Chem.MolFromSmiles(smiles_list[counter_01]))
                BackTracker.append(counter_01)
            except:
                Error_01.append([smiles_list[counter_01], counter_01])

        Mols = []
        
        Before = After
    
        for counter_02 in range (0, len(mols)):
            try:
                Mols.append(AllChem.AddHs(mols[counter_02]))
            except:
                Error_02.append([smiles_list[BackTracker[Before + counter_02]], BackTracker[Before + counter_02]])
            After = Before + counter_02

        morgans = [AllChem.GetMorganFingerprintAsBitVect(Mol, radius = fp_radius, nBits= fp_bits, useChirality = True) for Mol in Mols]
        morgan_bits =  [morgan.ToBitString() for morgan in morgans] 
    
        pattern = re.compile('.{1}') 
        morgan_bits = [','.join(pattern.findall(morgan)) for morgan in morgan_bits]

        fp_list = []
        for bit in morgan_bits:
            single_fp = bit.split(',')   # split the string by commas
            single_fp = [float(fp) for fp in single_fp] # transfer string to float32
            fp_list.append(single_fp)

        fp_df = pd.DataFrame(np.array(fp_list))
        fp_df.columns = fp_df.columns.astype(str)

# rename the columns

        for i in range(fp_df.columns.shape[0]):
            fp_df.rename(columns = {fp_df.columns[i]:fp_df.columns[i] + "pChEMBL"}, inplace = True)
            
        U.save_dataset(fp_df, path = ('./datasets/' + file + '/'), file_name = ('FP_' + str(fp_radius) + '_' + str(fp_bits) + '_D' + str(division - 1)), idx = False)
        
        print('\nCongratulation and say thanks to your computer! All Completed.')
        
    return Error_01, Error_02
