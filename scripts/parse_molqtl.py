#!/usr/bin/env python3

# general libraries
import os
import importlib
import subprocess
import logging
import sys
import json

# time and date
import time
from datetime import datetime, timedelta

# argument parsing
import argparse

# data manipulation
import pandas as pd
import polars as pl

# Change log:
# * v1.2.4 2024-12-04: Fixed dynamic creation of the log-file name to include the analysis type, and ensure it is saved in the correct output directory.
# * v1.2.3 2024-10-03: Fixed dynamic annotation of CpGs for gene-to-variant lookup in mQTLs.
# * v1.2.2 2024-10-03: Updated the dynamic Ensembl ID annotation function to work with `mygene` instead of `pybiomart`.
# * v1.2.1 2024-10-03: Added the `--type` to annotate the target list with Ensembl IDs if the type is 'gene'; changed the lookup accordingly.
# * v1.1.1 2024-10-03: Added the `--analysis` argument to the parser. Added logger.info statements to print the analysis type.
# * v1.1.0 2024-05-31: Initial version. 
# Version and license information 
VERSION_NAME = 'Parse molQTL results'
VERSION = '1.2.4'
VERSION_DATE = '2024-12-04'
COPYRIGHT = 'Copyright 1979-2024. Sander W. van der Laan | s.w.vanderlaan [at] gmail [dot] com | https://vanderlaanand.science.'
COPYRIGHT_TEXT = '''
The MIT License (MIT).

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and 
associated documentation files (the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial 
portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT 
LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Reference: http://opensource.org.
'''

# Define some global directories
# molQTL data directories 
# eQTL data directories
NOM_CIS_EQTL_loc = "/Users/slaan3/git/CirculatoryHealth/molqtl/results/version1_aernas1_firstrun/nom_cis_eqtl"
NOM_CIS_EQTL_SEX_loc = "/Users/slaan3/git/CirculatoryHealth/molqtl/results/version1_aernas1_firstrun/eqtl_gender"
NOM_CIS_EQTL_SEXINT_loc = "/Users/slaan3/git/CirculatoryHealth/molqtl/results/version1_aernas1_firstrun/sex_int_annot"
NOM_CIS_EQTL_SMOKING_loc = "/Users/slaan3/git/CirculatoryHealth/molqtl/results/version1_aernas1_firstrun/eqtl_smoking"
NOM_CIS_EQTL_SMOKINGINT_loc = "/Users/slaan3/git/CirculatoryHealth/molqtl/results/version1_aernas1_firstrun/smoking_int_annot"
PERM_TRANS_EQTL_loc = "/Users/slaan3/git/CirculatoryHealth/molqtl/results/version1_aernas1_firstrun/perm_trans_eqtl"

# mQTL data directories 
PERM_CIS_MQTL_loc = "/Users/slaan3/git/CirculatoryHealth/molqtl/results/perm_cis_mqtl"
PERM_TRANS_MQTL_loc = "/Users/slaan3/git/CirculatoryHealth/molqtl/results/perm_trans_mqtl"

# Reference and GWAS data directories 
REF_loc = "/Users/slaan3/PLINK/references"
GWAS_loc = "/Users/slaan3/PLINK/_GWAS_Datasets"

# Define the UtrechtSciencePark color scheme as a dictionary
# Website to convert HEX to RGB: http://hex.colorrrs.com.
# For some functions you should divide these numbers by 255.
UTRECHT_COLOR_SCHEME = {
    1: "#FBB820",   # yellow
    2: "#F59D10",   # gold
    3: "#E55738",   # salmon
    4: "#DB003F",   # darkpink
    5: "#E35493",   # lightpink
    6: "#D5267B",   # pink
    7: "#CC0071",   # hardpink
    8: "#A8448A",   # lightpurple
    9: "#9A3480",   # purple
    10: "#8D5B9A",  # lavendel
    11: "#705296",  # bluepurple
    12: "#686AA9",  # purpleblue
    13: "#6173AD",  # lightpurpleblue
    14: "#4C81BF",  # seablue
    15: "#2F8BC9",  # skyblue
    16: "#1290D9",  # azurblue
    17: "#1396D8",  # lightazurblue
    18: "#15A6C1",  # greenblue
    19: "#5EB17F",  # seaweedgreen
    20: "#86B833",  # yellowgreen
    21: "#C5D220",  # lightmossgreen
    22: "#9FC228",  # mossgreen
    23: "#78B113",  # lightgreen (for X chromosome)
    24: "#49A01D",  # green (for Y chromosome)
    25: "#595A5C",  # grey (for XY)
    26: "#A2A3A4",  # lightgrey (for MT)
    # Additional colors if needed
    27: "#D7D8D7",  # midgrey
    28: "#ECECEC",  # very lightgrey
    29: "#FFFFFF",  # white
    30: "#000000"   # black
}

# Set up logging
def setup_logger(script_name, verbose, analysis_type=None, output_dir=None):
    """Setup the logger to log to a file with the date and script name, and also log to the console."""
    # Get current date in the format YYYYMMDD
    date_str = datetime.now().strftime('%Y%m%d')

    # Construct log file name
    # log_file = f"{date_str}.{script_name}.log"
    log_file = f"{date_str}.{script_name}"
    if analysis_type:
        log_file += f".{analysis_type}"
    log_file += ".log"

    # Ensure the log file is in the specified output directory
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)  # Ensure the directory exists
        log_file = os.path.join(output_dir, log_file)

    # Initialize the logger
    logger = logging.getLogger(script_name)
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)

    # Create file handler to log to a file
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG)

    # Create console handler to print to console
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG if verbose else logging.INFO)

    # Formatter for the log messages
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    # Add both handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)

    return logger

# Function to check if a package is installed and install it if it is not
def check_install_package(package_name):
    """Check if a package is installed and install it if it is not."""
    try:
        importlib.import_module(package_name)
    except ImportError:
        logger.info(f'{package_name} is not installed. Installing it now...')
        subprocess.check_call(['pip', 'install', package_name])

# Function to map Ensembl IDs to gene information using MyGeneInfo
def map_gene_symbols_ids(gene_symbols, verbose=False, debug=False, logger=None):
    """
    Map GeneSymbols to Ensembl IDs and expanded gene information using MyGeneInfo.
    
    Parameters:
    - gene_symbols: List of GeneSymbols to map.
    
    Returns:
    - A DataFrame with expanded gene information.
    """
    
    # Import MyGeneInfo package
    check_install_package('mygene')
    logger.debug(f'  > Importing the MyGeneInfo package.')
    import mygene
    
    mg = mygene.MyGeneInfo()
    
    # Query MyGeneInfo with the list of GeneSymbols
    logger.debug(f'  > Querying MyGeneInfo with the list of GeneSymbols.')
    gene_info = mg.querymany(gene_symbols, 
                             scopes='symbol',  # Scope set to 'symbol' for GeneSymbol query
                             fields='ensembl.gene,symbol,name,alias,type_of_gene,map_location,genomic_pos,entrezgene,HGNC,summary', # Fields to retrieve
                             species='human')

    # Convert the query response into a DataFrame
    logger.debug(f'  > Converting the MyGeneInfo response to a DataFrame.')
    gene_info_df = pd.DataFrame(gene_info)

    # Check if the 'notfound' column exists
    if 'notfound' in gene_info_df.columns:
        missing_symbols = gene_info_df[gene_info_df['notfound'] == True]['query'].to_list()
        if missing_symbols:
            logger.warning(f"Some GeneSymbols were not found in MyGeneInfo: {missing_symbols}")
    else:
        logger.info("No missing GeneSymbols in the MyGeneInfo response.")
    
    # Filter out unnecessary columns and keep relevant information
    logger.debug(f'  > Filtering out unnecessary columns and keeping relevant information.')
    gene_info_df = gene_info_df[['query', 'ensembl', 'symbol', 'name', 'alias', 'type_of_gene', 
                                 'map_location', 'genomic_pos', 'entrezgene', 'HGNC', 'summary']]

    # Ensure no missing data breaks downstream code
    gene_info_df.fillna('NA', inplace=True)

    # Flatten the Ensembl IDs, handling cases where it might be a list, dict, or None
    logger.debug(f'  > Flattening the Ensembl IDs.')
    def extract_ensembl_id(ensembl):
        if isinstance(ensembl, dict):
            return ensembl.get('gene', None)  # Extract 'gene' field if present
        elif isinstance(ensembl, list):
            return ensembl[0]['gene'] if ensembl and isinstance(ensembl[0], dict) else None  # Take the first gene if it's a list of dicts
        return None
    logger.debug(f'  > Flattened.')
    gene_info_df['EnsemblID'] = gene_info_df['ensembl'].apply(extract_ensembl_id)

    # Flatten alias field (convert list to a single string)
    if 'alias' in gene_info_df.columns:
        logger.debug(f'  > Flattening the alias field.')
        gene_info_df['alias'] = gene_info_df['alias'].apply(lambda x: ', '.join(x) if isinstance(x, list) else x)

    # Split genomic_pos into multiple columns
    if 'genomic_pos' in gene_info_df.columns:
        gene_info_df['genomic_pos_chr'] = gene_info_df['genomic_pos'].apply(lambda x: x.get('chr', 'NA') if isinstance(x, dict) else 'NA')
        gene_info_df['genomic_pos_start'] = gene_info_df['genomic_pos'].apply(lambda x: x.get('start', 'NA') if isinstance(x, dict) else 'NA')
        gene_info_df['genomic_pos_end'] = gene_info_df['genomic_pos'].apply(lambda x: x.get('end', 'NA') if isinstance(x, dict) else 'NA')
        gene_info_df['genomic_pos_strand'] = gene_info_df['genomic_pos'].apply(lambda x: x.get('strand', 'NA') if isinstance(x, dict) else 'NA')
        
        # Drop the original genomic_pos column
        logger.debug(f'  > Dropping the original genomic_pos column.')
        gene_info_df.drop(columns=['genomic_pos'], inplace=True)

    # Drop the original 'ensembl' column
    logger.debug(f'  > Dropping the original Ensembl column.')
    gene_info_df.drop(columns=['ensembl'], inplace=True)

    # Show the first 5 rows of the DataFrame plus the shape for debugging after flattening
    logger.debug(f'  > Shape of the DataFrame: {gene_info_df.shape[0]} rows, {gene_info_df.shape[1]} columns.')
    logger.debug(f'  > First 5 rows of the DataFrame:\n{gene_info_df.head()}\n')
    
    logger.info(f"GeneSymbol annotations with Ensembl IDs and gene info extracted from MyGeneInfo.")
    return gene_info_df

# Function to map Gene Symbol to CpG IDs
def map_genes_to_cpg(gene_symbols):
    """
    Map gene symbols to CpG IDs using the IlluminaHumanMethylation450kanno.ilmn12.hg19 R package.
    
    Parameters:
    - gene_symbols: List of gene symbols to map.
    
    Returns:
    - A dictionary mapping gene symbols to CpG IDs.
    """
    
    # Import rpy2 and activate pandas-R conversion only when the function is called
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    
    # Activate pandas-R data conversion
    pandas2ri.activate()
    
    # Import R's base package
    base = ro.r['base']
    
    # Load the IlluminaHumanMethylation450kanno.ilmn12.hg19 package in R
    ro.r('library(IlluminaHumanMethylation450kanno.ilmn12.hg19)')
    
    # Access the manifest and extract CpG mappings
    ro.r('data(IlluminaHumanMethylation450kanno.ilmn12.hg19)')
    
    # The following command extracts the CpG to gene mapping from the annotation package
    ro.r('annotation_data <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)')
    
    # Convert the annotation data to a pandas DataFrame
    cpg_data_df = ro.r('annotation_data')
    cpg_data_df = pandas2ri.ri2py(cpg_data_df)
    
    # Select relevant columns for mapping: "UCSC_RefGene_Name" for genes and "Name" for CpG IDs
    cpg_data_df = cpg_data_df[['UCSC_RefGene_Name', 'Name']]
    
    # Create a dictionary where GeneSymbol is the key and CpG_ID is the value
    gene_to_cpg = {}
    for gene_symbol in gene_symbols:
        cpg_ids = cpg_data_df[cpg_data_df['UCSC_RefGene_Name'].str.contains(gene_symbol)]['Name'].tolist()
        if cpg_ids:
            gene_to_cpg[gene_symbol] = ','.join(cpg_ids)
        else:
            gene_to_cpg[gene_symbol] = 'NA'  # If no mapping found
    
    return gene_to_cpg

# Function to update target list with expanded gene information
def update_target_list(target_data, gene_info_df, output_file, verbose=False, debug=False, logger=None):
    """
    Update the target gene list with expanded gene information (e.g., symbol, name, alias, etc.)
    and save it as an Excel file.
    
    Parameters:
    - target_data: Polars DataFrame containing the original target gene data.
    - gene_info_df: DataFrame containing the expanded gene information from MyGeneInfo.
    - output_file: Path to save the updated file.
    
    Returns:
    - Updated target DataFrame with expanded gene information.
    """

    # Convert the Polars DataFrame to a pandas DataFrame
    logger.debug(f'  > Converting the target DataFrame to a pandas DataFrame.')
    target_data_pd = target_data.to_pandas()

    # Merge the target list with the expanded gene information on EnsemblID
    logger.debug(f'  > Merging the target list with the gene information.')
    updated_data = pd.merge(target_data_pd, gene_info_df, how='left', left_on='EnsemblID', right_on='query')
    
    # Convert the merged DataFrame back to Polars DataFrame
    logger.debug(f'  > Converting the updated DataFrame back to Polars.')
    updated_data_pl = pl.DataFrame(updated_data)
    
    # Write the updated data to an Excel file
    logger.info(f"Saving the updated target list to '{output_file}'.")
    with fastexcel.Workbook(output_file) as workbook:
        workbook.write_sheet("Genes", updated_data)
    
    logger.info(f"GeneSymbol annotations with Ensembl IDs and gene info saved to '{output_file}'.")
    return updated_data_pl  # Return the updated target DataFrame for downstream processing

# Function to merge and export molQTL results
def molqtl_merge_and_export(target_variants, sumstats, left_col, right_col, sort_column, output_csv, verbose=False, debug=False, logger=None):
    """
    Merges target variants with molQTL summary statistics and exports the result to a CSV file.

    Parameters:
    - target_variants: Polars DataFrame containing target variants.
    - sumstats: Polars DataFrame containing molQTL summary statistics.
    - left_col: Column name in target_variants to merge on.
    - right_col: Column name in sumstats to merge on.
    - sort_column: Column name to sort the merged DataFrame.
    - output_csv: Output CSV file path.
    - verbose: Boolean flag to control verbosity.
    - debug: Boolean flag to control debug logging.

    Returns:
    - None
    """
    # Debug: log shape and first 5 rows of the input QTL datasets
    if debug:
        logger.debug(f'  > Target variants DataFrame shape: {target_variants.shape[0]} rows, {target_variants.shape[1]} columns.')
        logger.debug(f'  > First 5 rows of target variants:\n{target_variants.head()}\n')
        logger.debug(f'  > Summary statistics DataFrame shape: {sumstats.shape[0]} rows, {sumstats.shape[1]} columns.')
        logger.debug(f'  > First 5 rows of summary statistics:\n{sumstats.head()}\n')

    logger.info(f'  > Starting merge of target variants or genes with molQTL summary statistics.')
    temp = target_variants.join(sumstats, left_on=left_col, right_on=right_col, how="inner")

    logger.info(f'  > Sorting the DataFrame by column "{sort_column}" in descending order.')
    result = temp.sort(sort_column)

    # Filter to keep only the top variant for each EnsemblID
    if left_col == "EnsemblID":  # Ensure this is for gene-type lookups
        logger.info(f'  > Dropping duplicates to keep only the top variant for each EnsemblID.')
        result = result.unique(subset="EnsemblID", keep="first")  # Keep the first occurrence of each EnsemblID

    # Debug: log shape and first 5 rows after the merge
    if debug:
        logger.debug(f'  > Merged DataFrame shape: {result.shape[0]} rows, {result.shape[1]} columns.')
        logger.debug(f'  > First 5 rows of merged DataFrame:\n{result.head()}\n')

    if debug:
        logger.debug(f'  > Removing temporary DataFrame from memory.')
        del temp

    logger.info(f'  > Exporting the Polars DataFrame to a CSV file.')
    result.write_csv(output_csv)

# Main function
def main():
    """Main function to parse molQTL results."""

# Parse command-line arguments
    parser = argparse.ArgumentParser(description=f'''
+ {VERSION_NAME} v{VERSION} +

This script extracts molQTL results from the Athero-Express Biobank Study. 
The phenotype of interest is given by `--trait` and the project name is given by `--project_name`. 
As an example and by default, the script uses `PCSK9` as `--trait` and `--project_name`.

The script requires the genome build of the results (`--build`), which can be either b37 or b38. Currently, 
the script only supports b37. 
The script also requires the type of lookup to perform (`--analysis`), which can be either:
- `cise` (cis-eQTL, nominal), which includes cis-eQTL results from the nominal analysis,
- `cism` (cis-mQTL, permuted), which includes cis-mQTL results from the permuted analysis (as the nominal 
   analysis yields too many results and is computationally expensive),
- `transe` (trans-eQTL, permuted), which includes the trans-eQTL results from the permuted analysis,
- `transm` (trans-mQTL, permuted), which includes the trans-mQTL results from the permuted analysis, or
- `all` (get all results), this is the default. 
Currently, the script does not support the analysis types `cisesmoking`, `cisesex`, `cisesmokingint`, and 
`cisesexint`.
The script also requires the type of lookup to perform (`--type`), which can be either:
- `variant` for variant-to-gene(s) lookup, or
- `gene` for gene-to-variant(s) lookup. 
In the `--type variant` case for cis- and trans-mQTL results it will automatically map 
EnsemblIDs to CpGs and save these IDs in a file with the CpG IDs called 'targets/targets.annot.xlsx'.
In the `--type gene` case it will automatically map EnsemblIDs to gene names 
and save these IDs in a file with the gene names called 'targets/targets.annot.xlsx'.

The script can print extra information using the `--verbose` argument. 

The `--debug` argument prints debug information. Note: this creates a lot of output - think carefully.

The `--version` argument prints the version and exits the script.

Example usage:
    python parse_molqtl.py --trait PCSK9 --project_name PCSK9_project --build b37 --analysis cise --verbose
        ''',
        epilog=f'''
+ {VERSION_NAME} v{VERSION}. {COPYRIGHT} \n{COPYRIGHT_TEXT}+''', 
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--trait', '-t', type=str, default='PCSK9', help='Phenotype of interest (default: PCSK9). Required.')
    parser.add_argument('--project_name', '-p', type=str, default='PCSK9', help='Name of the project (default: PCSK9). Required.')
    parser.add_argument('--analysis', '-a', type=str, default='all',
    choices=['cise', 'transe', 'cism', 'transm', 'all', 'cisesmoking', 'cisesex', 'cisesmokingint', 'cisesexint'], 
    help='Type of lookup to perform. (default: cise). Choices:\n'
        '  cise (cis-eQTL, nominal),\n'
        '  transe (trans-eQTL, permuted),\n'
        '  cism (cis-mQTL, permuted),\n'
        '  transm (trans-mQTL, permuted),\n'
        '  cisesmoking (cis-eQTL, smoking) -- not available yet,\n'
        '  cisesex (cis-eQTL, sex) -- not available yet,\n'
        '  cisesmokingint (cis-eQTL, smoking interaction) -- not available yet,\n'
        '  cisesexint (cis-eQTL, sex interaction) -- not available yet,\n'
        '  all (get all results), this the default.\n'
        '  Required.')
    parser.add_argument('--build', '-b', type=str, default='b37', choices=['b37', 'b38'], help='Genome build of the results. (default: b37). Choices: b37, b38. Required.')
    parser.add_argument('--type', '-y', type=str, default='variant', choices=['variant', 'gene'], 
        help="Type of lookup to perform. Choices: 'variant' for variant-to-gene(s), 'gene' for gene-to-variant(s) lookup. Required.")
    parser.add_argument('--verbose', '-v', action='store_true', help='Print extra information. Optional.')
    parser.add_argument('--debug', '-d', action='store_true', help='Print debug information. Note: this creates a lot of output - think carefully. Optional.')
    parser.add_argument('--version', '-V', action='version', version=f'%(prog)s {VERSION} ({VERSION_DATE}).')
    args = parser.parse_args()

    # Start the timer
    start_time = time.time()
    
    # Set some general defaults
    TRAIT_OF_INTEREST = args.trait
    PROJECTNAME = args.project_name
    BUILD = args.build

    # Output directory for the molQTL results
    molQTL_loc = "molQTL_results"

    # Get today's date
    FORMATTED_TODAY = today_date.strftime("%Y%m%d")

    # Setup logger after parsing arguments
    global logger
    # logger = setup_logger(parse_molqtl_results, args.verbose)
    logger = setup_logger('parse_molqtl_results', args.verbose, analysis_type=args.analysis, output_dir=molQTL_loc)

    # Start the script
    logger.info(f"+ {VERSION_NAME} v{VERSION} ({VERSION_DATE}) +")
    logger.info(f"\nStarting extraction job {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.")

    # Check if required arguments are provided
    required_args = ['trait', 'project_name', 'analysis', 'build']
    missing_args = [arg for arg in required_args if not getattr(args, arg)]
    if missing_args:
        logger.error(f"Error. The following required arguments are missing: {', '.join(missing_args)}.\n")
        parser.print_help()
        sys.exit(1)
    if args.build == 'b38':
        logger.error(f"Error. The script only supports b37 at the moment. Please provide the correct genome build.\n")
        parser.print_help()
        sys.exit(1)
    unsupported_analyses = ['cisesmoking', 'cisesex', 'cisesmokingint', 'cisesexint']
    if args.analysis in unsupported_analyses:
        logger.error(
            f"The analysis type '{args.analysis}' is not supported at the moment. "
            f"Supported types are: cise, transe, cism, transm, all."
        )
        parser.print_help()
        sys.exit(1)

    # Print extra information 
    logger.info(f"Trait................... {args.trait}")
    logger.info(f"Project name............ {args.project_name}")
    logger.info(f"Analysis type........... {args.analysis}")
    logger.info(f"Genome build............ {args.build}")
    logger.info(f"Verbose mode............ {'On' if args.verbose else 'Off (default)'}")
    logger.info(f"Debug mode.............. {'On' if args.debug else 'Off (default)'}")
    logger.info(f"Running version......... {VERSION} ({VERSION_DATE})")
    logger.info(f"Running on.............. {today_date}\n\n")

    # Create directories if they don't exist
    for directory in [molQTL_loc]:
        logger.info(f'Creating the directory (if it did not exist): {directory}')
        if not os.path.exists(directory):
            os.makedirs(directory)

    # Check contents of reference and GWAS directories
    if args.debug:
        logger.debug(f'Checking contents of the reference directory:')
        logger.debug(os.listdir(REF_loc))
        logger.debug(f'Checking contents of the GWAS directory:')
        logger.debug(os.listdir(GWAS_loc))

    # Check contents of molQTL data directories
    if args.debug:
        logger.debug(f'Checking contents of the molQTL data directories.')
        for directory in [NOM_CIS_EQTL_loc, PERM_CIS_MQTL_loc, PERM_TRANS_EQTL_loc, PERM_TRANS_MQTL_loc]:
            logger.debug(f'Checking contents of the directory: {directory}')
            logger.debug(os.listdir(directory))

    # Load target variants or genes based on --type argument
    logger.info(f'> Loading the target list from the Excel-file \'targets.xlsx\' based on type [{args.type}].')
    if args.type == 'variant':
        target_sheet = 'Variants'
    else:
        target_sheet = 'Genes'

    try:
        target_variants = pl.read_excel(
            source=os.path.join("targets/targets.xlsx"),
            sheet_name=target_sheet
        )
    except FileNotFoundError:
        logger.error(f"The target variants file was not found in sheet {target_sheet}.")
        sys.exit(1)
    except Exception as e:
        logger.error(f"An error occurred while reading the target variants file: {e}")
        sys.exit(1)

    # Annotate GeneSymbols in the Genes sheet with Ensembl IDs and gene info
    if args.type == 'gene' and args.analysis in ['cise', 'transe', 'all']:
        logger.info(f"Annotating GeneSymbol in the 'Genes' sheet with Ensembl IDs and gene information.")
        
        # Extract the list of GeneSymbols from the Genes sheet
        logger.debug(f'Extracting GeneSymbols from the target list.')
        gene_symbols = target_variants['GeneSymbol'].unique().to_list()
        
        # Use MyGeneInfo to map GeneSymbols to expanded gene information
        logger.debug(f'Mapping GeneSymbols to Ensembl IDs and gene information.')
        gene_info_df = map_gene_symbols_ids(gene_symbols, verbose=args.verbose, debug=args.debug, logger=logger)
        
        # Merge the annotations into the target list without overwriting existing data
        logger.debug(f'Merging gene information with the target list.')
        target_variants_pd = target_variants.to_pandas()

        # Only merge if the target DataFrame doesn't already contain EnsemblID info
        if 'EnsemblID' not in target_variants_pd.columns:
            logger.debug(f'Merging gene information with the target list.')
            target_variants_pd = pd.merge(target_variants_pd, gene_info_df, how='left', left_on='GeneSymbol', right_on='symbol')

        # Convert the updated DataFrame back to a Polars DataFrame
        logger.debug(f'Converting the updated DataFrame back to Polars.')
        # Replace 'NA' with None before converting back to Polars
        target_variants_pd = target_variants_pd.replace('NA', None)

        # Ensure columns with numeric data are of numeric type (e.g., converting strings to ints or floats where necessary)
        logger.info(f"Converting columns with numeric data to the appropriate type. If it throws a warning, the column will be left as is. This is expected.")
        for col in target_variants_pd.columns:
            if target_variants_pd[col].dtype == 'object':  # Check if the column is of object type (strings)
                try:
                    target_variants_pd[col] = pd.to_numeric(target_variants_pd[col])  # Attempt conversion to numeric
                except ValueError:
                    logger.warning(f"Could not convert column '{col}' to numeric. Retaining original values.")

        updated_target_variants = pl.DataFrame(target_variants_pd)

        # Save to the output Excel file (append mode if the file already exists)
        logger.debug(f'Saving the updated target list to an Excel file.')
        output_file = os.path.join(molQTL_loc, "targets.annot.xlsx")
        # Write the DataFrame to an Excel file
        target_variants_pd.to_excel(output_file, sheet_name='Genes', index=False)

        logger.info(f"GeneSymbol annotations with Ensembl IDs and gene info saved to '{output_file}'.")

    # Annotate CpG symbols in mQTLs for 'cism' and 'transm' analyses
    if args.type == 'variant' and args.analysis in ['cism', 'transm', 'all']:
        logger.info(f"Annotating CpG symbols in the mQTL summary statistics.")
        
        # Load the mQTL data (assuming CpG symbols are in a column called 'CpG')
        if args.analysis in ['cism', 'all']:
            sumstats_perm_cis_mqtl = pl.read_csv(os.path.join(PERM_CIS_MQTL_loc, "tensormqtl.perm_cis_mqtl.txt"), has_header=True, separator="\t")
        if args.analysis in ['transm', 'all']:
            sumstats_perm_trans_mqtl = pl.read_parquet(os.path.join(PERM_TRANS_MQTL_loc, "tensormqtl_perm_trans_qtl_pairs.annot.parquet"))
        
        # Extract CpG symbols from the mQTL data
        cpg_symbols = sumstats_perm_cis_mqtl['CpG'].unique().to_list() if 'CpG' in sumstats_perm_cis_mqtl.columns else []
        
        # Map CpG symbols to gene information
        cpg_mapping = map_genes_to_cpg(cpg_symbols)
        
        # Convert cpg_mapping to a DataFrame and merge with the mQTL data
        cpg_df = pd.DataFrame(list(cpg_mapping.items()), columns=['GeneSymbol', 'CpG_ID'])
        sumstats_perm_cis_mqtl = sumstats_perm_cis_mqtl.to_pandas()

        # Ensure CpG_ID is not already in the DataFrame to prevent overwriting
        if 'CpG_ID' not in sumstats_perm_cis_mqtl.columns:
            sumstats_perm_cis_mqtl = pd.merge(sumstats_perm_cis_mqtl, cpg_df, how='left', left_on='CpG', right_on='CpG_ID')
        
        # Save the annotated mQTL data to a file (append mode if needed)
        annotated_mqtl_file = os.path.join(molQTL_loc, "annotated_mqtl_results.csv")
        sumstats_perm_cis_mqtl.to_csv(annotated_mqtl_file, index=False)

        logger.info(f"Annotated mQTL data saved to '{annotated_mqtl_file}'.")

    # Load nominal cis-eQTL data
    if args.analysis == 'cise' or args.analysis == 'all':
        logger.info(f'> Loading nominal cis-eQTL data.')
        try:
            sumstats_nom_cis_eqtl = pl.read_parquet(source=os.path.join(NOM_CIS_EQTL_loc, "tensorqtl_nominal_cis_qtl_pairs.annot.parquet"))
        except FileNotFoundError:
            logger.error("The nominal cis-eQTL data file was not found.")
            sys.exit(1)
        except Exception as e:
            logger.error(f"An error occurred while reading the nominal cis-eQTL data: {e}")
            sys.exit(1)

        # Merge and export nominal cis-eQTL data
        # ┌─────────────────┬───────────┬──────────────┬──────────┬───┬────────┬───────────┬───────────────────────┬──────────┐
        # │ EnsemblID       ┆ VariantID ┆ tss_distance ┆ CAF_eQTL ┆ … ┆ TotalN ┆ MAF       ┆ MissingDataProportion ┆ HWE_P    │
        # │ ---             ┆ ---       ┆ ---          ┆ ---      ┆   ┆ ---    ┆ ---       ┆ ---                   ┆ ---      │
        # │ str             ┆ str       ┆ i32          ┆ f32      ┆   ┆ i64    ┆ f64       ┆ f64                   ┆ f64      │
        # ╞═════════════════╪═══════════╪══════════════╪══════════╪═══╪════════╪═══════════╪═══════════════════════╪══════════╡
        # │ ENSG00000187634 ┆ 1:693731  ┆ -240519      ┆ 0.134185 ┆ … ┆ 2124   ┆ 0.122957  ┆ 0.000075              ┆ 0.544841 │
        # │ ENSG00000187634 ┆ 1:714596  ┆ -219654      ┆ 0.038339 ┆ … ┆ 2124   ┆ 0.0322165 ┆ 0.000024              ┆ 0.272504 │
        # │ ENSG00000187634 ┆ 1:715367  ┆ -218883      ┆ 0.038339 ┆ … ┆ 2124   ┆ 0.0347528 ┆ 0.000019              ┆ 0.176551 │
        # │ ENSG00000187634 ┆ 1:717485  ┆ -216765      ┆ 0.038339 ┆ … ┆ 2124   ┆ 0.0343883 ┆ 0.000022              ┆ 0.174817 │
        # │ ENSG00000187634 ┆ 1:720381  ┆ -213869      ┆ 0.038339 ┆ … ┆ 2124   ┆ 0.0355205 ┆ 0.00002               ┆ 0.110452 │
        # └─────────────────┴───────────┴──────────────┴──────────┴───┴────────┴───────────┴───────────────────────┴──────────┘
        logger.info(f'> Merging and exporting nominal cis-eQTL data.')
        nom_cis_eqtl_out_file=os.path.join(molQTL_loc, FORMATTED_TODAY + "." + PROJECTNAME + "." + BUILD + "." + TRAIT_OF_INTEREST + ".nom_cis_eqtl.csv")
        if args.type == 'gene':
            molqtl_merge_and_export(updated_target_variants, sumstats_nom_cis_eqtl, "EnsemblID", "EnsemblID", "pval_nominal", nom_cis_eqtl_out_file, verbose=args.verbose, debug=args.debug, logger=logger)
        else:
            molqtl_merge_and_export(target_variants, sumstats_nom_cis_eqtl, "VariantID", "VariantID", "pval_nominal", nom_cis_eqtl_out_file, verbose=args.verbose, debug=args.debug, logger=logger)

        # Clean up
        if args.debug:
            logger.debug(f'Removing the nominal cis-eQTL data from memory.')
        del sumstats_nom_cis_eqtl

    # Load nominal cis-mQTL data
    if args.analysis == 'cism' or args.analysis == 'all':
        logger.info(f'> Loading nominal cis-mQTL data.')
        file_path_perm_cis_mqtl = os.path.join(PERM_CIS_MQTL_loc, "tensormqtl.perm_cis_mqtl.txt")
        try:
            sumstats_perm_cis_mqtl = pl.read_csv(file_path_perm_cis_mqtl, has_header=True, separator="\t", ignore_errors=True)
        except FileNotFoundError:
            logger.error("The nominal cis-mQTL data file was not found.")
            sys.exit(1)
        except Exception as e:
            logger.error(f"An error occurred while reading the nominal cis-mQTL data: {e}")
            sys.exit(1)

        # Merge and export nominal cis-mQTL data
        # ┌──────────────┬─────────┬─────────────┬─────────────┬───┬───────────┬───────────┬──────────┬────────────────────────┐
        # │ phenotype_id ┆ num_var ┆ beta_shape1 ┆ beta_shape2 ┆ … ┆ pval_perm ┆ pval_beta ┆ qval     ┆ pval_nominal_threshold │
        # │ ---          ┆ ---     ┆ ---         ┆ ---         ┆   ┆ ---       ┆ ---       ┆ ---      ┆ ---                    │
        # │ str          ┆ i64     ┆ f64         ┆ f64         ┆   ┆ f64       ┆ f64       ┆ f64      ┆ f64                    │
        # ╞══════════════╪═════════╪═════════════╪═════════════╪═══╪═══════════╪═══════════╪══════════╪════════════════════════╡
        # │ cg21870274   ┆ 309     ┆ 1.01037     ┆ 15.0882     ┆ … ┆ 0.531247  ┆ 0.530703  ┆ 0.468602 ┆ 0.001232               │
        # │ cg08258224   ┆ 731     ┆ 1.03658     ┆ 40.9219     ┆ … ┆ 0.984102  ┆ 0.986523  ┆ 0.605861 ┆ 0.000509               │
        # │ cg18147296   ┆ 750     ┆ 1.02594     ┆ 43.8044     ┆ … ┆ 0.0017    ┆ 0.001693  ┆ 0.006203 ┆ 0.000454               │
        # │ cg13938959   ┆ 787     ┆ 1.00503     ┆ 41.991      ┆ … ┆ 0.49525   ┆ 0.48971   ┆ 0.45302  ┆ 0.000433               │
        # │ cg12445832   ┆ 787     ┆ 1.01013     ┆ 42.2714     ┆ … ┆ 0.118388  ┆ 0.118428  ┆ 0.211842 ┆ 0.000439               │
        # └──────────────┴─────────┴─────────────┴─────────────┴───┴───────────┴───────────┴──────────┴────────────────────────┘
        logger.info(f'> Merging and exporting nominal cis-mQTL data.')
        perm_cis_mqtl_out_file=os.path.join(molQTL_loc, FORMATTED_TODAY + "." + PROJECTNAME + "." + BUILD + "." + TRAIT_OF_INTEREST + ".perm_cis_mqtl.csv")
        molqtl_merge_and_export(target_variants, sumstats_perm_cis_mqtl, "VariantID", "variant_id", "pval_nominal", perm_cis_mqtl_out_file, verbose=args.verbose, debug=args.debug, logger=logger)

        # Clean up
        if args.debug:
            logger.debug(f'Removing the nominal cis-mQTL data from memory.')
        del sumstats_perm_cis_mqtl

    # Load nominal trans-eQTL data
    if args.analysis == 'transe' or args.analysis == 'all':
        logger.info(f'> Loading the permuted trans-eQTL data.')
        file_path_perm_trans_eqtl = os.path.join(PERM_TRANS_EQTL_loc, "tensorqtl_trans_full.trans_qtl_pairs.parquet")
        try:
            sumstats_perm_trans_eqtl = pl.read_parquet(file_path_perm_trans_eqtl)
        except FileNotFoundError:
            logger.error("The permuted trans-eQTL data file was not found.")
            sys.exit(1)
        except Exception as e:
            logger.error(f"An error occurred while reading the permuted trans-eQTL data: {e}")
            sys.exit(1)

        # Merge and export nominal trans-eQTL data
        # ┌────────────┬─────────────────┬──────────┬───────────┬──────────┬──────────┬───────────────────┐
        # │ variant_id ┆ phenotype_id    ┆ pval     ┆ b         ┆ b_se     ┆ af       ┆ __index_level_0__ │
        # │ ---        ┆ ---             ┆ ---      ┆ ---       ┆ ---      ┆ ---      ┆ ---               │
        # │ str        ┆ str             ┆ f64      ┆ f32       ┆ f32      ┆ f32      ┆ i64               │
        # ╞════════════╪═════════════════╪══════════╪═══════════╪══════════╪══════════╪═══════════════════╡
        # │ 1:730087   ┆ ENSG00000131408 ┆ 0.000009 ┆ -0.419756 ┆ 0.093878 ┆ 0.010383 ┆ 0                 │
        # │ 1:752307   ┆ ENSG00000108039 ┆ 0.000006 ┆ -0.308768 ┆ 0.067392 ┆ 0.019169 ┆ 1                 │
        # │ 1:752593   ┆ ENSG00000108039 ┆ 0.000006 ┆ -0.308768 ┆ 0.067392 ┆ 0.019169 ┆ 2                 │
        # │ 1:752617   ┆ ENSG00000108039 ┆ 0.000006 ┆ -0.308768 ┆ 0.067392 ┆ 0.019169 ┆ 3                 │
        # │ 1:754121   ┆ ENSG00000108039 ┆ 0.000006 ┆ -0.308768 ┆ 0.067392 ┆ 0.019169 ┆ 4                 │
        # └────────────┴─────────────────┴──────────┴───────────┴──────────┴──────────┴───────────────────┘

        logger.info(f'> Merging and exporting the permuted trans-eQTL data.')
        perm_trans_eqtl_out_file=os.path.join(molQTL_loc, FORMATTED_TODAY + "." + PROJECTNAME + "." + BUILD + "." + TRAIT_OF_INTEREST + ".perm_trans_eqtl.csv")
        if args.type == 'gene':
            molqtl_merge_and_export(updated_target_variants, sumstats_perm_trans_eqtl, "EnsemblID", "phenotype_id", "pval", perm_trans_eqtl_out_file, verbose=args.verbose, debug=args.debug, logger=logger)
        else:
            molqtl_merge_and_export(target_variants, sumstats_perm_trans_eqtl, "VariantID", "variant_id", "pval", perm_trans_eqtl_out_file, verbose=args.verbose, debug=args.debug, logger=logger)

        # Clean up
        if args.debug:
            logger.debug(f'Removing the permuted trans-eQTL data from memory.')
        del sumstats_perm_trans_eqtl

    # Load nominal trans-mQTL data
    if args.analysis == 'transm' or args.analysis == 'all':
        logger.info(f'> Loading the permuted trans-mQTL data.')
        file_path_perm_trans_mqtl = os.path.join(PERM_TRANS_MQTL_loc, "tensormqtl_perm_trans_qtl_pairs.annot.parquet")
        try: 
            sumstats_perm_trans_mqtl = pl.read_parquet(file_path_perm_trans_mqtl)
        except FileNotFoundError:
            logger.error("The permuted trans-mQTL data file was not found.")
            sys.exit(1)
        except Exception as e:
            logger.error(f"An error occurred while reading the permuted trans-mQTL data: {e}")
            sys.exit(1)

        # Merge and export nominal trans-mQTL data
        # ┌───────────┬────────────┬───────────┬───────────┬───┬────────┬───────────┬───────────────────────┬──────────┐
        # │ VariantID ┆ CpG        ┆ pval_perm ┆ Beta      ┆ … ┆ TotalN ┆ MAF       ┆ MissingDataProportion ┆ HWE_P    │
        # │ ---       ┆ ---        ┆ ---       ┆ ---       ┆   ┆ ---    ┆ ---       ┆ ---                   ┆ ---      │
        # │ str       ┆ str        ┆ f64       ┆ f32       ┆   ┆ i64    ┆ f64       ┆ f64                   ┆ f64      │
        # ╞═══════════╪════════════╪═══════════╪═══════════╪═══╪════════╪═══════════╪═══════════════════════╪══════════╡
        # │ 1:693731  ┆ cg23922040 ┆ 0.000007  ┆ -0.356294 ┆ … ┆ 2124   ┆ 0.122957  ┆ 0.000075              ┆ 0.544841 │
        # │ 1:714596  ┆ cg26079250 ┆ 4.7736e-7 ┆ 0.601047  ┆ … ┆ 2124   ┆ 0.0322165 ┆ 0.000024              ┆ 0.272504 │
        # │ 1:715367  ┆ cg26079250 ┆ 4.7736e-7 ┆ 0.601047  ┆ … ┆ 2124   ┆ 0.0347528 ┆ 0.000019              ┆ 0.176551 │
        # │ 1:717485  ┆ cg26079250 ┆ 4.7736e-7 ┆ 0.601047  ┆ … ┆ 2124   ┆ 0.0343883 ┆ 0.000022              ┆ 0.174817 │
        # │ 1:720381  ┆ cg26079250 ┆ 4.7736e-7 ┆ 0.601047  ┆ … ┆ 2124   ┆ 0.0355205 ┆ 0.00002               ┆ 0.110452 │
        # └───────────┴────────────┴───────────┴───────────┴───┴────────┴───────────┴───────────────────────┴──────────┘
        logger.info(f'> Merging and exporting the permuted trans-mQTL data.')
        perm_trans_mqtl_out_file=os.path.join(molQTL_loc, FORMATTED_TODAY + "." + PROJECTNAME + "." + BUILD + "." + TRAIT_OF_INTEREST + ".perm_trans_mqtl.csv")
        molqtl_merge_and_export(target_variants, sumstats_perm_trans_mqtl, "VariantID", "VariantID", "pval_perm", perm_trans_mqtl_out_file, verbose=args.verbose, debug=args.debug, logger=logger)

        # Clean up
        if args.debug:
            logger.debug(f'Removing the permuted trans-mQTL data from memory.')
        del sumstats_perm_trans_mqtl

    # Calculate and print execution time
    elapsed_time = time.time() - start_time
    time_delta = timedelta(seconds=elapsed_time)
    formatted_time = str(time_delta).split('.')[0]
    logger.info(f"Script executed on {today_date.strftime('%Y-%m-%d')}. Total execution time: {formatted_time}.")

    # Print the version and license information
    logger.info(f"{VERSION_NAME} v{VERSION} ({VERSION_DATE}). {COPYRIGHT}")
    logger.info(COPYRIGHT_TEXT.strip())

if __name__ == "__main__":
    main()
# End of file