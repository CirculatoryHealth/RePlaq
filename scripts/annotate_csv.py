#!/usr/bin/env python3

# general libraries
import argparse
import logging
import sys
import os
import gzip

# data manipulation
import polars as pl

# time and date
import time
from datetime import datetime, timedelta

# Change log:
# * v1.1.0 2024-10-18: Initial version. 
# Version and license information 
VERSION_NAME = 'AnnotateCSV'
VERSION = "1.0.0"
VERSION_DATE = '2024-10-18'
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

# Set up logging
def setup_logger(output_base, verbose):
    """Setup the logger to log to a file in the same directory as the output file and also log to the console."""
    # Get current date in the format YYYYMMDD
    date_str = datetime.now().strftime('%Y%m%d')
    
    # Create log file path based on the output base with date prefix
    log_file = f"{os.path.dirname(output_base)}/{date_str}.{os.path.basename(output_base)}.log"

    # Ensure the directory for the log file exists
    log_dir = os.path.dirname(log_file)
    if log_dir and not os.path.exists(log_dir):
        os.makedirs(log_dir)

    # Initialize the logger
    logger = logging.getLogger(output_base)
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


# Function to annotate genes
def annotate_genes(dataframe, scope_type):
    # gene annotation
    import mygene
    # Initialize MyGene client
    mg = mygene.MyGeneInfo()

    '''Annotate a list of genes with biotype, description, and summary.'''
    # Extract gene symbols or Ensembl IDs from the relevant column
    genes = dataframe["gene"].to_list()
    
    logger.info(f"> Start annotating {len(genes)} genes using {scope_type}...\n")

    # Check for empty or missing gene entries
    if any(not gene or gene == "NA" for gene in genes):
        logger.warning(f"Some entries in the gene column are empty or missing. They will be annotated with 'NA'.")
    
    # Set the scope based on whether it's a gene symbol or Ensembl ID
    query_scope = "symbol" if scope_type == "symbol" else "ensembl.gene"
    
    # Query MyGene.info for annotations
    annotations = mg.querymany([gene for gene in genes if gene and gene != "NA"], scopes=query_scope, fields='symbol,name,type_of_gene,map_location,entrezgene,HGNC,summary,disease', species='human')
    
    # Possible queries mygene for additional information
    # annotations = mg.querymany([gene for gene in genes if gene and gene != "NA"], scopes=query_scope, fields='symbol,name,alias,type_of_gene,map_location,genomic_pos,entrezgene,HGNC,summary,pathway,go', species='human')
    # query	symbol	name	alias	type_of_gene	map_location	genomic_pos	entrezgene	HGNC	summary
    # ENSG00000187634	SAMD11	Sterile Alpha Motif...	None	protein-coding	1q25.1	{"start": 92312, ...}	148398	37102	Involved in signaling...
    # ENSG00000079974	CENPM	Centromere Protein M	None	protein-coding	22q12.1	{"start": 44922049, ...}	55880	17675	Part of centromere complex...

    # Process results into a dictionary
    annotations_dict = {ann['query']: ann for ann in annotations if 'symbol' in ann}
    
    # Add annotation columns to the dataframe
    descriptions = []
    biotypes = []
    map_location = []
    entrezgene = []
    hgnc = []
    summaries = []
    diseases = []
    
    logger.debug("> Adding annotations to the dataframe for the following genes:")
    for gene in genes:
        if not gene or gene == "NA": # Append 'NA' if missing or empty
            descriptions.append('NA')
            biotypes.append('NA')
            map_location.append('NA')
            entrezgene.append('NA')
            hgnc.append('NA')
            summaries.append('NA')
            diseases.append('NA')
        else:
            logger.debug(f"  > [{gene}]")
            ann = annotations_dict.get(gene, {})
            descriptions.append(ann.get('name', 'NA'))
            biotypes.append(ann.get('type_of_gene', 'NA'))
            map_location.append(ann.get('map_location', 'NA'))
            entrezgene.append(ann.get('entrezgene', 'NA'))
            hgnc.append(ann.get('HGNC', 'NA'))
            summaries.append(ann.get('summary', 'NA'))
            diseases.append(ann.get('disease', 'NA'))
    
    dataframe = dataframe.with_columns([
        pl.Series('gene_name', descriptions),
        pl.Series('biotype', biotypes),
        pl.Series('map_location', map_location),
        pl.Series('entrezgene', entrezgene),
        pl.Series('HGNC', hgnc),
        pl.Series('summary', summaries),
        pl.Series('disease', diseases)
    ])
    
    return dataframe

def main():

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Annotate a list of genes with biotype, description, and summary.")
    """Main function to parse molQTL results."""

# Parse command-line arguments
    parser = argparse.ArgumentParser(description=f'''
+ {VERSION_NAME} v{VERSION} +

Annotate a list of genes with biotype, description, and summary.

The scripts will annotate the given CSV-file (`--input`) which at the minimum contains the following column:
- gene: GeneSymbol or Ensembl ID. It can either use Ensembl gene IDs or gene symbols and use this information 
to query MyGene.info return a CSV file with the annotations. By default, the script will try to detect the gene
column name automatically. If the column name is different, it can be specified using the `--gene` argument.

Annotations include:
- biotype: The biotype of the gene.
- description: The description of the gene.
- summary: A summary of the gene.
- entrezgene: The Entrez Gene ID.
- HGNC: The HGNC gene symbol.
- map_location: The chromosomal location of the gene.

The script will output a CSV file with the annotated data. By default, the output file will be saved in the same
directory as the input file with the suffix '_annotated'. The output file can be compressed with gzip using the
`--compression` argument.

The `--verbose` argument enables verbose output. The `--version` argument prints the version and exits the script.

Example usage:
    python annotate.csv.py --input genes.csv --gene gene --output annotated_genes.csv --compression gzip --verbose
        ''',
        epilog=f'''
+ {VERSION_NAME} v{VERSION}. {COPYRIGHT} \n{COPYRIGHT_TEXT}+''', 
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--input", required=True, help="Input CSV file containing the gene data. For example: input.csv Required.")
    parser.add_argument("-g", "--gene", required=False, help="Name of the gene column, e.g. gene, symbol, ensemblid, etc. By default it will be determined on the spot. Optional.")
    parser.add_argument("-o", "--output", required=False, help="Output CSV file to save the annotated data. Default: input_annotated.csv. Optional.")
    parser.add_argument("-c", "--compression", choices=["gzip", "regular"], default="regular", help="Specify output compression: gzip or regular. Default: regular. Optional.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose output.")
    parser.add_argument("-V", "--version", action="version", version=f"%(prog)s {VERSION}", help="Show program's version number and exit.")
    
    args = parser.parse_args()
    
    # Get the directory where the input file is located
    input_dir = os.path.dirname(args.input)
    
    # Get current date in the format YYYYMMDD
    date_str = datetime.now().strftime('%Y%m%d')

    # Set output base name: if --output is not given, take the basename of the --input and append "_annotated"
    if args.output:
        output_base = os.path.join(input_dir, f"{date_str}.{os.path.basename(args.output)}")
    else:
        input_basename = os.path.splitext(os.path.basename(args.input))[0]  # Get input file basename without extension
        output_base = os.path.join(input_dir, f"{date_str}.{input_basename}_annotated")  # Use the input directory for the output

    # Setup logging in the same directory as input
    global logger
    logger = setup_logger(output_base, args.verbose)
    
    # Start the timer
    start_time = time.time()

    # Log script start
    today_date = datetime.now()
    logger.info(f"+ {VERSION_NAME} v{VERSION} ({VERSION_DATE}) +")
    logger.info(f"Starting job at {today_date.strftime('%Y-%m-%d %H:%M:%S')}.\n")

    # Print extra information 
    logger.info(f"Input................... {args.input}")
    logger.info(f"Gene.................... {args.gene if args.gene else 'Auto-detect'}")
    logger.info(f"Output base name........ {output_base}")
    logger.info(f"Compression............. {args.compression}")
    logger.info(f"Verbose mode............ {'On' if args.verbose else 'Off'}\n")

    # Read the input CSV file
    logger.debug(f"> Reading input CSV file: {args.input}")
    try:
        df = pl.read_csv(args.input)
    except Exception as e:
        logger.error(f"Error reading CSV file: {e}")
        sys.exit(1)

    # Check for gene symbol or EnsemblID column names
    logger.debug("> Checking for gene or EnsemblID column.")
    symbol_columns = ["gene", "Gene", "GENE", "GeneSymbol", "Symbol", "symbol"]
    ensembl_columns = ["EnsemblID", "ENSlID", "ensembleid"]

    # If the user provided the column name via -g/--gene, use that column.
    if args.gene:
        if args.gene in df.columns:
            gene_column = args.gene
        else:
            logger.error(f"Specified column '{args.gene}' not found in the input file.")
            sys.exit(1)
    else:
        # Automatically find the first matching column from symbol or ensembl lists
        gene_column = next((col for col in df.columns if col in symbol_columns + ensembl_columns), None)

    if not gene_column:
        logger.error(f"Missing required column. Expected one of {symbol_columns + ensembl_columns}.")
        sys.exit(1)

    # If a column named 'gene' already exists, rename it to avoid conflicts
    if 'gene' in df.columns and gene_column != 'gene':
        logger.warning(f"'gene' column already exists. Renaming existing 'gene' column to avoid conflict.")
        df = df.rename({'gene': 'original_gene'})

    # Rename the identified gene column to 'gene'
    df = df.rename({gene_column: "gene"})

    # Determine the scope type (symbol or ensembl)
    scope_type = "symbol" if gene_column in symbol_columns else "ensembl"

    # Annotate the genes
    logger.debug(f"> Annotating genes using {scope_type}.\n")
    df = annotate_genes(df, scope_type)

    # Save the annotated data in the same directory as the input
    print(f"")
    logger.debug(f"> Saving annotated data to {output_base}.csv")
    output_file = f"{output_base}.csv"
    if args.compression == "gzip":
        logger.debug(f"> Compressing output file with gzip.")
        output_file += ".gz"
        with gzip.open(output_file, 'wt') as f:
            df.write_csv(f)
    else:
        df.write_csv(output_file)

    print(f"")
    logger.info(f"All done. Annotated data saved to [{output_file}].")

    # Calculate and print execution time
    elapsed_time = time.time() - start_time
    time_delta = timedelta(seconds=elapsed_time)
    formatted_time = str(time_delta).split('.')[0]
    logger.info(f"Total execution time: {formatted_time}.\n")

    # Print version and license information
    logger.info(f"{VERSION_NAME} v{VERSION} ({VERSION_DATE}). {COPYRIGHT}")
    logger.info(COPYRIGHT_TEXT.strip())

if __name__ == "__main__":
    main()
