import pandas as pd
import numpy as np
import re

def convert_epicv2_to_v1(input_file, output_file):
    print("Loading data...")
    betas = pd.read_csv(input_file, index_col=0)
    
    print("Extracting base names...")
    base_names = betas.index.str.replace(r'_[BT][CO]\d+$', '', regex=True)
    
    print("Converting to EPIC v1 format...")
    betas_v1 = betas.groupby(base_names).mean()
    
    print("Saving results...")
    betas_v1.to_csv(output_file)
    
    print(f"Conversion complete!")
    print(f"Original probes: {len(betas)}")
    print(f"EPIC v1 probes: {len(betas_v1)}")
    print(f"Output saved to: {output_file}")
    
    return None
