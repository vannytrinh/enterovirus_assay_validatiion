import pandas as pd
import numpy as np
from ete3 import NCBITaxa

ncbi = NCBITaxa()

# --- FUNCTIONS TO PARSE DATA ---

# date: date value of datasets file 
def datasets_get_year(date):
    '''
    Retrieve year from datasets file date elements 
    ex. '2000-01-09T00:00:00Z' --> '2000'
    '''
    if (date is not np.NaN):
        return int(str(date).split('-')[0])
    else:
        return np.NaN

# filepath: filepath to simulate_PCR file 
def get_accessions(filepath): 
    '''
    Returns list of hit accessions in simulate_PCR file 
    '''
    # read sim_pcr data
    data = pd.read_csv(filepath, sep='\t')

    # get list of pcr hits
    acc = pd.DataFrame()
    acc['Accession'] = data['Full_Hit_ID'].apply(lambda r: str(r).split(' ')[0])
    acc = acc.drop_duplicates(subset=["Accession"], keep='first')

    return acc["Accession"].tolist()

def filter_data(file):
    # read in datasets file
    data = pd.read_csv(file, sep='\t')

    # filter for only human hosts
    data = data[data['Host Taxonomic ID'] == 9606]
    # filter out missing collection dates
    data = data[~data['Isolate Collection date'].isna()]
    # filter for complete sequences
    data = data[data['Completeness'] == 'COMPLETE']
    
    return data.copy()

# return df with accession number, taxid, and if the strain was hit by assay for all potential assay targets
# df: datatsets dataframe
# pcr_acc: list of hit accessions
def assess_data(data, pcr_acc):
    '''
    Parse dataset file and indicate if accession was hit by assay 
    '''

    # make new copy of df 
    df = data.copy()

    # assign boolean if accession was hit
    df['Hit'] = df['Accession'].apply(lambda x: x in pcr_acc)

    # retreive collection year
    df['Collection year'] = df['Isolate Collection date'].apply(datasets_get_year)
    df['Release year'] = df['Release date'].apply(datasets_get_year)

    # return df w/ relevant columns
    return df[['Accession', 'Virus Taxonomic ID', 'Collection year', 'Release year', 'Hit']].copy()

# --- FUNCTIONS TO ANALYZE DATA ---

# return df summarizing the hits and totals for each taxid from input df
# df: df from load_data, with rows of accession, taxid, and hit boolean
def count_data(df):
    '''
    Summarize total and hit data by taxid 
    '''
    counts = df.groupby(['Virus Taxonomic ID'])['Hit'] \
                  .agg([lambda x: (x == True).sum(), 'count']).reset_index()
    counts.columns = ['Virus Taxonomic ID','Hits', 'Total']
    return counts

def collapse_data(counts, taxonomy_groups):
    df = counts.merge(taxonomy_groups, on='Virus Taxonomic ID', how='left')
    df = df.drop(columns=['Virus Taxonomic ID'])
    df = df.groupby(['Collapse Name', 'Collapse TaxId']).sum().reset_index()
    df = df.rename(columns={'Collapse Name': 'Virus Group', 'Collapse TaxId': 'TaxID'})
    return df.copy()

def summarize_assay(assay_data, collapse_info):
    counts = count_data(assay_data)
    counts_collapsed = collapse_data(counts, collapse_info)
    return counts_collapsed

def summarize_assays(assays, assay_names, collapse_info):
    summaries = []

    for assay, assay_name in zip(assays, assay_names):
        summary = summarize_assay(assay, collapse_info)
        summary = summary.rename(columns={'Hits':f'{assay_name} Hits'})
        summaries.append(summary)

    combined_summary = summaries[0]
    for summary in summaries[1:]:
        combined_summary = pd.merge(combined_summary, summary, on=['Virus Group', 'TaxID', 'Total'], how='outer')

    column_order = ['Virus Group', 'TaxID'] + [f'{assay_name} Hits' for assay_name in assay_names] + ['Total']

    return combined_summary[column_order].copy()

# --- FUNCTIONS FOR HEATMAP --- 

def aggregate_counts(data, years, target_taxid):
    '''
    Expands given data so that it contains data from each given year 
    Data for a given year includes data collected up to that year 
    '''
    # collect data for each time point
    all_data = []

    # iter through each input year 
    for yr in years:
        # count data for each taxid collected before or in given year
        counts = count_data(data[data["Collection year"] <= yr])
        # collapse data by species level
        counts = collapse_data(counts, target_taxid)
        counts.rename(columns={'Hits': 'Hits '+str(yr), 'Total': 'Total '+str(yr)}, inplace=True)
    
        all_data.append(counts)

    # combine data for each timepoint into a combined dataframe
    combined_df = all_data[0]
    for df in all_data[1:]:
        combined_df = pd.merge(combined_df, df, on=['Virus Group', 'TaxID'], how='outer')
    combined_df = combined_df.fillna(0)
    
    combined_df.set_index('Virus Group', inplace=True)

    return combined_df

def calc_totals_ratios(aggregated_counts):
    hits = aggregated_counts.filter(regex="^Hits")
    totals = aggregated_counts.filter(regex="^Total")
    
    hits.columns = hits.columns.str.replace('Hits ', '')
    totals.columns = totals.columns.str.replace('Total ', '')
    
    ratios = hits / totals
    
    return totals, ratios

def separate_year_types(combined_df):
    collection = combined_df.filter(like="Collection").copy()
    release = combined_df.filter(like="Release").copy()
    
    collection.columns = collection.columns.str.replace(' Collection', '')
    release.columns = release.columns.str.replace(' Release', '')
    
    collection.index.name = 'Collection Year'
    release.index.name = 'Release Year'
    
    return collection, release

# --- FUNCTIONS FOR TIMEPLOT ---

def count_years(df):
    # collect all types of counts here
    # Total Collection, Total Release, Hit Collection, Hit Release
    all_df = []

    # df of only hit accessions 
    hit_df = df[df['Hit']]

    all_df.append(df[['Collection year']].value_counts().sort_index().to_frame()\
                  .rename(columns={'count': 'Total Collection'}))
    all_df.append(df[['Release year']].value_counts().sort_index().to_frame()\
                  .rename(columns={'count':'Total Release'}))
    all_df.append(hit_df[['Collection year']].value_counts().sort_index().to_frame()\
                  .rename(columns={'count':'Hit Collection'}))
    all_df.append(hit_df[['Release year']].value_counts().sort_index().to_frame()\
                  .rename(columns={'count':'Hit Release'}))

    counts = pd.concat(all_df, axis=1)

    counts = counts.reset_index()
    counts.set_index('level_0', inplace=True)
    counts.index.name = 'Year'

    counts.index = counts.index.astype(int)
    min_year = counts.index[0]
    max_year = counts.index[-1]
    new_index = range(min_year, max_year+1)
    counts = counts.reindex(new_index)

    return counts

def make_cumulative(counts):
    cumulative = counts.fillna(0)
    cumulative = cumulative.cumsum()
    cumulative = cumulative.replace(0, np.nan)
    return cumulative