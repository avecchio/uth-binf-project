import json
import xmltodict
import ensembl_rest
import requests
import shutil
import gzip
import psycopg2
import datetime
import matplotlib.pyplot as plt
import liftover
import shutil
import urllib.request as reques        
import os.path
import wget
import numpy as np
from scipy.stats import chisquare

from liftover import get_lifter
from Bio import SeqIO
from contextlib import closing
from os import path


import secrets

from scipy.stats import chi2_contingency

from collections import defaultdict, OrderedDict
import warnings
import gffutils
import pybedtools
import pandas as pd
import copy
import os
import re
from gffutils.pybedtools_integration import tsses
from copy import deepcopy
from collections import OrderedDict, Callable
import errno



def make_directory(directory):
    try:
        os.makedirs(directory)
    except OSError as e:
        pass

def unzip_file(local_filepath):
    unzipped_filename = local_filepath.replace(".gz","")
    with gzip.open(local_filepath, 'rb') as f_in:
        with open(unzipped_filename, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return unzipped_filename


def sync_databases(working_directory, file_name, url, unzip):
    if path.exists(f'./{working_directory}/{file_name}') == False:
        local_filename = file_name
        local_filepath = f'./{working_directory}/{local_filename}'
        if unzip:
            local_filepath = local_filepath + ".gz"

        wget.download(url, local_filepath)

        filename = local_filepath
        if unzip:
            filename = unzip_file(local_filepath)
        return filename

def myconverter(o):
    if isinstance(o, datetime.datetime):
        return o.__str__()

def db_cache(file_name, callback, callback_params):
    if path.exists(file_name):
        with open(file_name) as json_file:
            data = json.load(json_file)
            return data
    else:
        data = callback(callback_params)
        with open(file_name, 'w') as outfile:
            json_str = json.dumps(data, default = myconverter)
            outfile.write(json_str)
        return data

def plot_bar_chart(frequencies_dict, xlabel, ylabel, title, filename):
    keys = list(frequencies_dict.keys())
    values = list(frequencies_dict.values())

    plt.figure(figsize=(15, 15))
    plt.bar(keys, values, color='green') #, yerr=variance)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.xticks(rotation = 45)
    
    plt.savefig(filename, dpi=100)

def convert_coordinate(chromosome, coordinate):
    converter = get_lifter('hg19', 'hg38')
    return converter[chromosome][coordinate][0][1]

def extract_content(response, restype):
    if response.status_code != 200:
        print('error')
        return None
    else:
        if restype == 'json':
            return json.loads(response.text)
        else:
            return xmltodict.parse(response.text)

def is_gwas_snp_associated(snpid, condition):
    study_url = f'https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{snpid}/studies'
    response = requests.get(study_url)
    content = extract_content(response, 'json')
    studies = content['_embedded']['studies']
    for study in studies:
        trait = study['diseaseTrait']['trait']
        if condition in trait:
            return True
    return False

def get_gwas_alleles(snpid):
    server = "https://rest.ensembl.org"
    ext = f"/variation/human/{snpid}?"
    
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
        
    decoded = r.json()
    return repr(decoded)

def query_gwas(params):
    gene_name, condition = params
    associated_snps = []
    query_url = f'https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/search/findByGene?geneName={gene_name}'
    response = requests.get(query_url)
    content = extract_content(response, 'json')
    snps = content['_embedded']['singleNucleotidePolymorphisms']
    for snp in snps:
        snpid = snp['rsId']
        alleles = get_gwas_alleles(snpid)
        print('===========================')
        is_associated = is_gwas_snp_associated(snpid, condition)
        
        if is_associated:
            for allele in alleles['mappings']:
                associated_snps.append({
                    'reference_allele': allele['ancestral_allele'],
                    'alternate_allele': allele['allele_string'],
                    'start': int(allele['start']),
                    'stop': int(allele['end']),
                })
    return associated_snps

def query_ncbi(query_url, restype):
    response = requests.get(query_url)
    content = extract_content(response, restype)
    ids = content['esearchresult']['idlist']
    return ids

def query_clinvar(gene_name, is_test):
    retmax = 20 if is_test == True else 100000 
    query_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={gene_name}&retmax={retmax}&retmode=json'
    return query_ncbi(query_url, 'json')

def get_clinvar_entry(id, condition):
    query_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id={id}&from_esearch=true&retmode=json'
    res = requests.get(query_url)
    content = extract_content(res, 'xml')
    variants = []

    record = content['ClinVarResult-Set']['VariationArchive']['InterpretedRecord']
    try:        
        if 'TraitMappingList' in record:
            is_condition = False
            traits = record['TraitMappingList']['TraitMapping']
            if type(traits).__name__ == 'list':
                for trait in traits:
                    trait_name = trait['MedGen']['@Name']
                    is_condition = is_condition or condition in trait_name
            else:
                trait_name = traits['MedGen']['@Name']
                is_condition = is_condition or condition in trait_name
            locations = record['SimpleAllele']['Location']['SequenceLocation']

            if is_condition:
                for location in locations:                
                        if location['@Assembly'] == 'GRCh38':
                            variants.append({
                                'start': int(location['@start']),
                                'stop': int(location['@stop']),
                                'reference_allele': location['@referenceAlleleVCF'],
                                'alternate_allele': location['@alternateAlleleVCF']
                            })
    except Exception as e:
        print(e)
        print('-----------------')
    return variants

def pg_query(conn_string, sql_query):
    try:
        print('connecting')
        conn = psycopg2.connect(conn_string)
        cur = conn.cursor()
        print('executing query')
        cur.execute(sql_query)
        row = cur.fetchone()
        print('fetching results')

        results = []
        while row is not None:
            results.append(list(row))
            row = cur.fetchone()
        cur.close()
        return results

    except (Exception, psycopg2.DatabaseError) as error:
        print(error)
        return []
    finally:
        if conn is not None:
            conn.close()


def query_rnacentral(params):
    gene_name = params[0]
    """ query data from the vendors table """
    conn = None
    
    columns = '''
    acc.id AS acc_id,
    acc.accession AS acc_accession,
    acc.anticodon AS acc_anticodon,
    acc.common_name AS acc_common_name,
    acc.experiment AS acc_experiment,
    acc.gene AS acc_gene,
    acc.species AS acc_species,
    acc.standard_name AS acc_standard_name,
    r.has_coordinates AS r_has_coordinates,
    r.is_active AS r_is_active,
    sr.assembly_id AS sr_assembly,
    sr.chromosome AS sr_chromosome,
    sr.exon_count AS sr_exon_count,
    sr.region_name AS sr_region_name,
    sr.region_start AS sr_region_start,
    sr.region_stop AS sr_region_stop,
    sr.strand AS sr_strand,
    r.description as description,
    r.rna_type
    '''

    sql_query = f'''    
    SELECT
    {columns}
    from rnacen.rnc_accessions acc
    LEFT JOIN rnacen.xref x on (acc.accession = x.ac)
    left join rnacen.rnc_rna_precomputed r on (x.upi = r.upi)
    left join rnacen.rnc_sequence_regions sr on (r.id = sr.urs_taxid)
    '''

    conn_string = "host='hh-pgsql-public.ebi.ac.uk' dbname='pfmegrnargs' user='reader' password='NWDMCE5xdipIjRrp'"

    results = []
    where_clauses = [f"where gene like '%{gene_name}%'", f"where r.description like '%{gene_name}%'"]

    for clause in where_clauses:
        results = results + pg_query(conn_string, sql_query + clause)

    return [ 
        {"acc_id": x[0],
        "acc_accession": x[1],
        "acc_anticodon": x[2],
        "acc_common_name": x[3],
        "acc_experiment": x[4],
        "acc_gene": x[5],
        "acc_species": x[6],
        "acc_standard_name": x[7],
        "r_has_coordinates": x[8],
        "r_is_active": x[9],
        "sr_assembly": x[10],
        "sr_chromosome": x[11],
        "sr_exon_count": x[12],
        "sr_region_name": x[13],
        "sr_region_start": x[14],
        "sr_region_stop": x[15],
        "sr_strand": x[16],
        "description": x[17],
        "rna_type": x[18]}
    for x in results]

def get_ensembl_data(params):
    gene_name = params[0]
    results = ensembl_rest.symbol_lookup(
        species='homo sapiens',
        symbol=gene_name,
        params={'expand': True}
    )
    return results

def sync_gene_enhancers(working_directory):
    enhancer_paths = []
    cell_type_enhancers = [
        'A375','A549','AML_blast','Astrocyte','BJ','Bronchia_epithelial',
        'Caco-2','Calu-3','CD14+','CD19+','CD20+','CD34+',
        'CD36+','CD4+','CD8+','Cerebellum','CUTLL1','DOHH2',
        'ECC-1','ESC_neuron','Esophagus','Fetal_heart','Fetal_kidney','Fetal_muscle_leg',
        'Fetal_placenta','Fetal_small_intestine','Fetal_spinal_cord','Fetal_stomach','Fetal_thymus','FT246',
        'FT33','GM10847','GM12878','GM12891','GM12892','GM18505',
        'GM18526','GM18951','GM19099','GM19193','GM19238','GM19239',
        'GM19240','H1','H9','HCC1954','HCT116','HEK293T',
        'HEK293','Hela-S3','Hela','HepG2','HFF','HL-60',
        'hMADS-3','HMEC','hNCC','HSMM','HT1080','HT29',
        'HUVEC','IMR90','Jurkat','K562','Kasumi-1','KB',
        'Keratinocyte','Left_ventricle','LHCN-M2','Liver','LNCaP-abl','LNCaP',
        'Lung','MCF-7','MCF10A','ME-1','Melanocyte','melanoma',
        'Mesendoderm','MS1','Myotube','Namalwa','NB4','NHDF',
        'NHEK','NHLF','NKC','OCI-Ly7','Osteoblast','Ovary',
        'PANC-1','Pancreas','Pancreatic_islet','PBMC','PC3','PrEC',
        'SGBS_adipocyte','SK-N-SH','SK-N-SH_RA','Skeletal_muscle','Small_intestine','Sperm',
        'Spleen','T47D','T98G','th1','Thymus','U2OS',
        'VCaP','ZR75-30'
    ]

    for cell_type in cell_type_enhancers:
        url = f'http://www.enhanceratlas.org/data/AllEPs/hs/{cell_type}_EP.txt'
        try:
            cell_type_enhancer_name = f'{cell_type}.hg19.txt'
            sync_databases(working_directory, cell_type_enhancer_name, url, False)
            enhancer_paths.append(cell_type_enhancer_name)
        except:
            print('unable to download: ' + cell_type)

    return enhancer_paths

def extract_circular_rnas(file_path, gene_name, chromosome):
    circular_rnas = []
    with open(file_path) as fp:
        lines = fp.readlines()
        counter = 0
        for line in lines:
            entry = line.replace("\n", "").split("\t")
            if (len(entry) == 13):
                start = entry[1]
                end = entry[2]
                strand = entry[3]
                identifier = entry[4]
                if entry[11] == gene_name:
                    circular_rnas.append({
                        'identifier': identifier,
                        'start': convert_coordinate(chromosome, int(start)),
                        'end': convert_coordinate(chromosome, int(end)),
                        'type': 'circular_rna',
                        'strand': strand
                    })
            counter = counter + 1
    return circular_rnas

def get_coordinates(coord_str):
    chr, coordinates = coord_str.split(":")
    start, end = coordinates.split("-")
    return chr, start, end

def extract_enhancers(working_directory, file_path, gene_identifier, chromosome):
    enhancers = {}
    counter = 0
    with open(f'./{working_directory}/{file_path}') as fp:
        lines = fp.readlines()
        for line in lines:
            entry = line.replace("\n","").split("\t")
            if gene_identifier in line:
                location_str = entry[0].split("$")[0]
                coord_str, gene_id = location_str.split("_")
                if (gene_id == gene_identifier):
                    chr, start, end = get_coordinates(coord_str)
                    coordinate_id = f'{start}:{end}'
                    if coordinate_id not in enhancers:
                        enhancers[coordinate_id] = {
                            'identifier': f'Enhancer{counter}',
                            'start': convert_coordinate(chromosome, int(start)),
                            'end': convert_coordinate(chromosome, int(end)),
                            'type': 'enhancer',
                            'strand': '+'
                        }
            counter = counter + 1
    return list(enhancers.values())

def extract_promoters(file_path, gene_name):
    promoters = []
    counter = 0
    with open(file_path) as fp:
        lines = fp.readlines()
        for line in lines:
            entry = line.split(" ")
            if len(entry) == 8:
                gene = entry[3]
                start = entry[1]
                end = entry[2]
                strand = entry[5]
                if gene_name in gene:
                    promoters.append({
                        'identifier': f'promoter{counter}',
                        'start': int(start),
                        'end': int(end),
                        'type': 'promoter',
                        'strand': strand
                    })
            else:
                print('not a proper entry for file: ' + file_path)
            counter = counter + 1
    return promoters

def extract_insulators(file_path, gene_name, chromosome):
    known_insulators = []
    insulators = []
    counter = 0
    with open(file_path) as fp:
        lines = fp.readlines()
        for line in lines:
            if counter > 0 and len(line.split("\t")) > 4:
                entry = line.split("\t")
                identifier = entry[0]
                species = entry[1]
                coord_str = entry[2]

                five_prime_gene = ''
                three_prime_gene = ''

                if (len(coord_str.split(":")) == 2):
                    five_prime_gene = entry[3]
                    three_prime_gene = entry[4]
                else:
                    coord_str = entry[3]
                    five_prime_gene = entry[4]
                    three_prime_gene = entry[5]

                chr, start, end = get_coordinates(coord_str)

                flanking = (gene_name in five_prime_gene or gene_name in three_prime_gene)
                if flanking and (species == 'Human'):
                    insulators.append({
                        'identifier': identifier,
                        'start': convert_coordinate(chromosome, int(start)),
                        'end': convert_coordinate(chromosome, int(end)),
                        'type': 'insulator',
                        'strand': '+'
                    })
            counter = counter + 1
    return insulators

def extract_non_coding_rnas(non_coding_rnas):
    rna_regions = []
    for non_coding_rna in non_coding_rnas:
        identifier = non_coding_rna['ID']
        start = non_coding_rna['start']
        end = non_coding_rna['end']
        biotype = non_coding_rna['biotype']
        rna_regions.append()
    return rna_regions

def extract_sno_rnas(file_path, chromosome, gene_name):
    rnas = []
    with open(file_path) as fp:
        lines = fp.readlines()
        for line in lines:
            metadata = line.split("\t")
            information = metadata[8].strip()[1:-2].split(",")
            start = metadata[3]
            end = metadata[4]
            host_gene = information[4].split(":")[1]
            if gene_name in host_gene:
                rnas.append({
                    'identifier': information[2],
                    'start': convert_coordinate(chromosome, int(start)),
                    'end': convert_coordinate(chromosome, int(end)),
                    'type': 'snornas',
                    'strand': '?'
                })
    return rnas

def extract_genecode_features(file_path, ensembl_gene_id):
    features_dict = {}
    features = []
    with open(file_path) as fp:
        lines = fp.readlines()
        for line in lines:
            metadata = line.split("\t")
            if len(metadata) > 5:
                chr = metadata[0]
                start = metadata[3]
                end = metadata[4]
                biotype = metadata[2]
                information = metadata[8].split(";")
                identifier = information[0]
                gene_id = information[2]
                
                accepted_biotypes = [
                    'exon', 'three_prime_UTR', 'five_prime_UTR',
                ]

                is_type = (biotype in accepted_biotypes)
                if is_type and (ensembl_gene_id in gene_id):        
                    if biotype not in features_dict:
                        features_dict[biotype] = 0
                    features_dict[biotype] += 1            
                    features.append({
                        'identifier': identifier,
                        'start': int(start),
                        'end': int(end),
                        'type': biotype,
                        'strand': '+'
                    })

    print(features_dict)
    return features


def get_ncbi_clinical_variants(params):
    gene_name, condition = params
    variants = []
    clinical_variants = query_clinvar(gene_name, False)
    for clinical_variant in clinical_variants:
        variant_coordinates = get_clinvar_entry(clinical_variant, condition)
        variants = variants + variant_coordinates
    return variants

def dedup_regions(regions):
    unique_regions = []
    unique_region_coordinates = {}
    for region in regions:
        start = region['start']
        end = region['end']
        coordinate = f'{start}-{end}'
        if coordinate not in unique_region_coordinates:
            unique_region_coordinates[coordinate] = 0
            unique_regions.append(region)
    return unique_regions

def filter_nc_rnas(chromosome, nc_rnas):
    filtered_nc_rnas = {}
    rna_ids = []
    for rna in nc_rnas:
        start = rna['sr_region_start']
        end = rna['sr_region_stop']
        if start is not None and end is not None and rna['acc_species'] == 'Homo sapiens':
            coord_key = f'{start}-{end}'
            if coord_key not in filtered_nc_rnas:
                filtered_nc_rnas[coord_key] = rna

    nc_rnas_list = []
    for rna in filtered_nc_rnas:
        start = filtered_nc_rnas[rna]['sr_region_start']
        end = filtered_nc_rnas[rna]['sr_region_stop']
        biotype = filtered_nc_rnas[rna]['rna_type']
        identifier = filtered_nc_rnas[rna]['sr_region_name']
        assembly = filtered_nc_rnas[rna]['sr_assembly']
        if assembly == 'GRCh19':
            nc_rnas_list.append({
                'identifier': identifier,
                'start': convert_coordinate(chromosome, int(start)),
                'end': convert_coordinate(chromosome, int(end)),
                'type': biotype,
                'strand': strand
            })
        elif assembly == 'GRCh38':
            nc_rnas_list.append({
                'identifier': identifier,
                'start': int(start),
                'end': int(end),
                'type': biotype,
                'strand': strand
            })
        else:
            print(f'Error: Unknown assembly type [{assembly}] for {identifier}')

    return nc_rnas_list

def write_regions_to_gff3(gene_name, chromosome, regions):
    entries = []
    for region in regions:
        region_type = region['type']
        region_start = region['start']
        region_end = region['end']
        region_id = region['identifier']
        entry = f'{gene_name}_regions . {region_type} {region_start} {region_end} . + . ID={region_id};'
        entries.append(entry)

        gff3_contents = '\n'.join(entries)

        file1 = open(f"{gene_name}_regions.gff3","w")
        file1.write(gff3_contents) 
        file1.close()

def get_intron_counts(regions):
    intron_count = regions['transcript']
    for region_type in ['exon', 'three_prime_UTR', 'five_prime_UTR']:
        intron_count = intron_count - regions[region_type]
    del regions['transcript']
    return regions

def get_sequence_from_ensembl(chromosome, start, end):
    server = "https://rest.ensembl.org"
    ext = f"/sequence/region/human/{chromosome}:{start}..{end}:1?coord_system_version=GRCh38"
    
    r = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"})
    return r.text

def dna_to_rna(dna):
    return dna.replace("T", "U")

def generate_rna_structure(identifier, directory, rna_type, rna):

    make_directory(directory)
    tmp_path = secrets.token_hex(nbytes=16)
    structure_path = f'{directory}/{tmp_path}'
    make_directory(structure_path)

    fasta_name = f'{structure_path}/{identifier}.fasta'
    with open(fasta_name, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")
    is_circular = " -c " if rna_type == "circularRNA" else ""

    dbn_file = f'{structure_path}/{identifier}_{structure_path}.dbn'
    st_file = dbn_file.replace(".dbn", ".st")
    os.system(f'RNAfold {is_circular} {fasta_name} > {dbn_file}')
    os.system(f'perl bpRNA.pl {dbn_file} ')
    #os.system(f'mv {st_file} {structure_path}')

def calculate_random_chance_statistics(regions):
    total_number_of_variants = 0
    regional_length = 0
    affected_regional_length = 0

    for region in regoins:
        rna_start = rna_region['region']['start']
        rna_end = rna_region['region']['end']
        rna_length = end - start + 1
        regional_length += rna_length

        variants = rna_region['variants']
        number_of_variants = len(variants)
        if number_of_variants > 0:
            total_number_of_variants += number_of_variants
            affected_regional_length += rna_length

    expected_variants_in_region = total_number_of_variants / regional_length
    expected_variants_in_affected_regions = total_number_of_variants / affected_regional_length

    return expected_variants_in_region, expected_variants_in_affected_regions


# these are zero based. Should be -1 to all coordinates
# will need to restructure these and to account for all variants
def delete_sequence(dna, start, end):
    before = dna[0:start]
    after = dna[end+1:]
    return before + after

def insert_sequence(dna, start, end, sequence):
    before = dna[0:start+1]
    after = dna[end-1:]
    return before + sequence + after

def modify_sequence(dna, start, end, sequence):
    end = start + len(sequence)
    before = dna[0:start]
    after = dna[end:]
    return before + sequence + after

def mutate_dna(start, end, mutation, dna):
    return dna[0:start-1] + mutation + dna[end+1:]

def generate_structural_variants(regions):
    rna_regions = regions
    for rna_region in rna_regions:
        rna_identifier = rna_region['region']['identifier']
        rna_start = rna_region['region']['start']
        rna_end = rna_region['region']['end']
        rna_type = rna_region['region']['type']

        dna = get_sequence_from_ensembl(chromosome, start, end)
        rna = dna_to_rna(dna)
        #generate_rna_structure(identifier, rna_type, rna)

        for variant in rna_region['variants']:
            mutated_dna = mutate_dna(variant['start'], variant['stop'], variant['alternate_allele'], dna)
            mutated_rna = dna_to_rna(mutated_dna)
            #generate_rna_structure(identifier, rna_type, rna)

def count_overlapping_variant_frequencies(regions, variants):
    print(f'regions: ' + str(len(regions)))
    print(f'variants: ' + str(len(variants)))
    unique_variant_regions = {}
    for variant in variants:
        counter = 0
        for region in regions:
            after_start = region['start'] <= variant['start'] and region['end'] >= variant['start']
            before_end = region['start'] <= variant['stop'] and region['end'] >= variant['stop']
            if after_start or before_end:
                counter += 1
        if str(counter) not in unique_variant_regions:
            unique_variant_regions[str(counter)] = 0
        unique_variant_regions[str(counter)] += 1
    
    return unique_variant_regions

def count_emperical_regional_variant_frequencies(regions, variants):
    print(f'regions: ' + str(len(regions)))
    print(f'variants: ' + str(len(variants)))
    regional_frequencies = {}            
    region_type_lengths = {}

    for region in regions:
        if region['type'] not in regional_frequencies:
            regional_frequencies[region['type']] = {}
        if region['identifier'] not in regional_frequencies[region['type']]:
            regional_frequencies[region['type']][region['identifier']] = {
                'region': region,
                'variants': []
            }
        if region['type'] not in region_type_lengths:
            region_type_lengths[region['type']] = 0
        region_type_lengths[region['type']] += abs(region['end'] - region['start']) + 1

        for variant in variants:
            after_start = region['start'] <= variant['start'] and region['end'] >= variant['start']
            before_end = region['start'] <= variant['stop'] and region['end'] >= variant['stop']
            if after_start or before_end:
                regional_frequencies[region['type']][region['identifier']]['variants'].append(variant)

    return regional_frequencies, region_type_lengths

def count_statistical_regional_variant_frequencies(regions, variants):
    variant_regions = {}
    for variant in variants:
        regions_impacted = []
        counter = 0
        for region in regions:
            if region['type'] not in variant_regions:
                variant_regions[region['type']] = 0
            after_start = region['start'] <= variant['start'] and region['end'] >= variant['start']
            before_end = region['start'] <= variant['stop'] and region['end'] >= variant['stop']
            if after_start or before_end:
                region_type = region['type']
                if region_type not in regions_impacted:
                    regions_impacted.append(region_type)  
        num_regions_impacted = len(regions_impacted)
        if num_regions_impacted > 0:
            frequency_score = 1 / num_regions_impacted
            for region in regions_impacted:
                variant_regions[region] += frequency_score
    
    return variant_regions

def double_bar_chart(regional_frequencies, expected_regional_frequencies):

    regional_frequency_values = list(regional_frequencies.values())
    N = len(regional_frequency_values)
    ind = np.arange(N)
    width = 0.27

    fig = plt.figure()
    ax = fig.add_subplot(111)

    yvals =  regional_frequency_values
    rects1 = ax.bar(ind, yvals, width, color='b')
    zvals = list(expected_regional_frequencies.values())
    rects2 = ax.bar(ind+width, zvals, width, color='g')

    ax.set_ylabel('Scores')
    ax.set_xticks(ind+width)
    ax.set_xticklabels( list(regional_frequencies.keys()) )
    ax.legend( (rects1[0], rects2[0]), ('Actual distribution', 'Expected distribution') )

    def autolabel(rects):
        for rect in rects:
            h = rect.get_height()
            ax.text(rect.get_x()+rect.get_width()/2., 1.05*h, '%d'%int(h),
                    ha='center', va='bottom')

    autolabel(rects1)
    autolabel(rects2)

    plt.xticks(rotation = 45)
    
    plt.savefig('expected_actual_distribution.png', dpi=100)


def read_genecode_bed_file(biotype, gene, file_path):
    genecode_regions = []
    counter = 0
    with open(file_path) as fp:
        lines = fp.readlines()
        for line in lines:
            metadata = line.strip().split("\t")
            start = int(metadata[1])
            end = int(metadata[2])
            gene_name = metadata[3]
            strand = metadata[5]

            if gene_name in gene:
                genecode_regions.append({
                    'identifier': f'{biotype}-{str(counter)}',
                    'start': start,
                    'end': end,
                    'type': biotype,
                    'strand': strand
                })
        counter += 1
    return genecode_regions

def create_gene_dict(db):
    '''
    Store each feature line db.all_features() as a dict of dicts
    '''
    gene_dict = DefaultOrderedDict(lambda: DefaultOrderedDict(lambda: DefaultOrderedDict(list)))
    for line_no, feature in enumerate(db.all_features()):
        gene_ids = feature.attributes['gene_id']
        feature_type = feature.featuretype
        if feature_type == 'gene':
            if len(gene_ids)!=1:
                logging.warning('Found multiple gene_ids on line {} in gtf'.format(line_no))
                break
            else:
                gene_id = gene_ids[0]
                gene_dict[gene_id]['gene'] = feature
        else:
            transcript_ids = feature.attributes['transcript_id']

            for gene_id in gene_ids:
                for transcript_id in transcript_ids:
                    gene_dict[gene_id][transcript_id][feature_type].append(feature)
    return gene_dict

class DefaultOrderedDict(OrderedDict):
    # Source: http://stackoverflow.com/a/6190500/562769
    def __init__(self, default_factory=None, *a, **kw):
        if (default_factory is not None and
           not isinstance(default_factory, Callable)):
            raise TypeError('first argument must be callable')
        OrderedDict.__init__(self, *a, **kw)
        self.default_factory = default_factory

    def __getitem__(self, key):
        try:
            return OrderedDict.__getitem__(self, key)
        except KeyError:
            return self.__missing__(key)

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value

    def __reduce__(self):
        if self.default_factory is None:
            args = tuple()
        else:
            args = self.default_factory,
        return type(self), args, None, None, self.items()

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        return type(self)(self.default_factory, self)

    def __deepcopy__(self, memo):
        import copy
        return type(self)(self.default_factory,
                          copy.deepcopy(self.items()))

    def __repr__(self):
        return 'OrderedDefaultDict(%s, %s)' % (self.default_factory,
                                               OrderedDict.__repr__(self))

def get_gene_list(gene_dict):
    return list(set(gene_dict.keys()))

def get_UTR_regions(gene_dict, gene_id, transcript, cds):
    if len(cds)==0:
        return [], []
    utr5_regions = []
    utr3_regions = []
    utrs = gene_dict[gene_id][transcript]['UTR']
    first_cds = cds[0]
    last_cds = cds[-1]
    for utr in utrs:
        ## Push all cds at once
        ## Sort later to remove duplicates
        strand = utr.strand
        if strand == '+':
            if utr.stop < first_cds.start:
                utr.feature_type = 'five_prime_UTR'
                utr5_regions.append(utr)
            elif utr.start > last_cds.stop:
                utr.feature_type = 'three_prime_UTR'
                utr3_regions.append(utr)
            else:
                raise RuntimeError('Error with cds')
        elif strand == '-':
            if utr.stop < first_cds.start:
                utr.feature_type = 'three_prime_UTR'
                utr3_regions.append(utr)
            elif utr.start > last_cds.stop:
                utr.feature_type = 'five_prime_UTR'
                utr5_regions.append(utr)                
            else:
                raise RuntimeError('Error with cds')    
    return utr5_regions, utr3_regions

def create_bed(regions, bedtype='0'):
    '''Create bed from list of regions
    bedtype: 0 or 1
        0-Based or 1-based coordinate of the BED
    '''
    bedstr = ''
    for region in regions:
        if len(region.attributes['gene_id']) == 1:
            ## GTF start is 1-based, so shift by one while writing 
            ## to 0-based BED format
            if bedtype == '0':
                start = region.start - 1
            else:
                start = region.start
            bedstr += '{}\t{}\t{}\t{}\t{}\t{}\n'.format(region.chrom,
                                                start,
                                                region.stop,
                                                re.sub('\.\d+', '', region.attributes['gene_id'][0]),
                                                '.',
                                                region.strand)
    return bedstr

def rename_regions(regions, gene_id):
    #print(len(regions))
    regions = list(regions)
    if len(regions) == 0:
        return []
    for region in regions:
        region.attributes['gene_id'] = gene_id
    return regions

def merge_regions(db, regions):
    if len(regions) == 0:
        return []
    merged = db.merge(sorted(list(regions), key=lambda x: x.start))
    return merged

def merge_regions_nostrand(db, regions):
    if len(regions) == 0:
        return []
    merged = db.merge(sorted(list(regions), key=lambda x: x.start), ignore_strand=True)
    return merged

from bioinfokit.analys import gff

def create_genecode_db(gff_file):
    gtf_file = gff_file.replace(".gff3", ".gtf").split("/")[2]
    if path.exists(gtf_file) == False:
        gff.gff_to_gtf(file=gff_file)

    gtf_db_file = gff_file.replace(".gff3", ".gtf.db")
    print(gtf_db_file)
    print(path.exists(gtf_db_file))
    if path.exists(gtf_db_file) == False:
        print('uh...')
        db = gffutils.create_db(gtf_file, dbfn=gtf_db_file, merge_strategy='merge', force=True, disable_infer_transcripts=True, disable_infer_genes=True)

    db = gffutils.FeatureDB(gtf_db_file, keep_order=True)
    gene_dict = create_gene_dict(db)
    return gene_dict, db

def extract_features(working_directory, gff_file):
    data_files = [
        f'./{working_directory}/utr5.bed.gz',
        f'./{working_directory}/utr3.bed.gz',
        f'./{working_directory}/exon.bed.gz',
        f'./{working_directory}/intron.bed.gz'
    ]

    existing_data = True
    for data_file in data_files:
        existing_data = existing_data and path.exists(data_file)


    if existing_data == False:
        gene_dict, db = create_genecode_db(gff_file)
        utr5_bed = ''
        utr3_bed = ''
        gene_bed = ''
        exon_bed = ''
        intron_bed = ''
        start_codon_bed = ''
        stop_codon_bed = ''
        cds_bed = ''

        gene_list = []

        for gene_id in get_gene_list(gene_dict):
            gene_list.append(gene_dict[gene_id]['gene'])
            
            utr5_regions, utr3_regions = [], []
            exon_regions, intron_regions = [], []
            star_codon_regions, stop_codon_regions = [], []
            cds_regions = []
            
            for feature in gene_dict[gene_id].keys():
                if feature == 'gene':
                    continue
                cds = list(gene_dict[gene_id][feature]['CDS'])
                exons = list(gene_dict[gene_id][feature]['exon'])
                merged_exons = merge_regions(db, exons)
                introns = db.interfeatures(merged_exons)
                #utr5_region, utr3_region = get_UTR_regions(gene_dict, gene_id, feature, cds)
                utr5_region = list(gene_dict[gene_id][feature]['five_prime_utr'])
                utr3_region = list(gene_dict[gene_id][feature]['three_prime_utr'])
                utr5_regions += utr5_region
                utr3_regions += utr3_region
                exon_regions += exons
                intron_regions += introns
                cds_regions += cds
                
            merged_utr5 = merge_regions(db, utr5_regions)
            renamed_utr5 = rename_regions(merged_utr5, gene_id)
            
            #print(merged_utr5)

            merged_utr3 = merge_regions(db, utr3_regions)
            renamed_utr3 = rename_regions(merged_utr3, gene_id)
            
            #print(merged_utr3)

            #merged_exons = merge_regions(db, exon_regions)
            renamed_exons = rename_regions(exon_regions, gene_id)
            
            #print(merged_exons)

            merged_introns = merge_regions(db, intron_regions)
            renamed_introns = rename_regions(merged_introns, gene_id)
            
            #print(merged_introns)

            merged_cds = merge_regions(db, cds_regions)
            #print(merged_cds)
            renamed_cds = rename_regions(merged_cds, gene_id)
            
            utr3_bed += create_bed(renamed_utr3)
            utr5_bed += create_bed(renamed_utr5)
            exon_bed += create_bed(renamed_exons)
            intron_bed += create_bed(renamed_introns)
            cds_bed += create_bed(renamed_cds)

        gene_bed = create_bed(gene_list)
        gene_bedtool = pybedtools.BedTool(gene_bed, from_string=True)
        utr5_bedtool = pybedtools.BedTool(utr5_bed, from_string=True)
        utr3_bedtool = pybedtools.BedTool(utr3_bed, from_string=True)
        exon_bedtool = pybedtools.BedTool(exon_bed, from_string=True)
        intron_bedtool = pybedtools.BedTool(intron_bed, from_string=True)
        cds_bedtool = pybedtools.BedTool(cds_bed, from_string=True)

        utr5_cds_subtracted = utr5_bedtool.subtract(cds_bedtool)
        utr3_cds_subtracted = utr3_bedtool.subtract(cds_bedtool)
        utr5_cds_subtracted.remove_invalid().sort().saveas(os.path.join(working_directory, 'utr5.bed.gz'))
        utr3_cds_subtracted.remove_invalid().sort().saveas(os.path.join(working_directory, 'utr3.bed.gz'))
        gene_bedtool.remove_invalid().sort().saveas(os.path.join(working_directory, 'gene.bed.gz'))
        exon_bedtool.remove_invalid().sort().saveas(os.path.join(working_directory, 'exon.bed.gz'))
        intron_bedtool.remove_invalid().sort().saveas(os.path.join(working_directory, 'intron.bed.gz'))
        cds_bedtool.remove_invalid().sort().saveas(os.path.join(working_directory, 'cds.bed.gz'))

        for data_file in data_files:
            unzip_file(data_file)

    return data_files

def main():
    working_directory = 'data'
    gene_name = 'FTO'

    make_directory(working_directory)

    condition = 'Growth retardation'

    condition_filename = '_'.join(condition.split(" "))
    ensemble_cds_metadata = db_cache(f'./{working_directory}/{gene_name}_{condition_filename}_ensembl.json', get_ensembl_data, [gene_name])

    gene_id = ensemble_cds_metadata['id']
    region_name = ensemble_cds_metadata['seq_region_name']
    chromosome = f'chr{region_name}'
    #print(gene_id)
    #associated_enhancer_paths = sync_gene_enhancers(working_directory)
    #sync_databases(working_directory, 'Hs_EPDnew.bed', 'ftp://ccg.epfl.ch/epdnew/H_sapiens/current/Hs_EPDnew.bed', False)    
   
    #sync_databases(working_directory, 'gencode.v37.chr_patch_hapl_scaff.annotation.gff3', 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.chr_patch_hapl_scaff.annotation.gff3.gz', True)
    #gff3_file = f'./{working_directory}/gencode.v37.chr_patch_hapl_scaff.annotation.gff3'
    #print('extracting features')
    #feature_files = extract_features(working_directory, gff3_file)
    #sync_databases(working_directory, 'Supplementary_Dataset_S5.gff', 'http://snoatlas.bioinf.uni-leipzig.de/Supplementary_Dataset_S5.gff', False)
    #sync_databases(working_directory, 'human-circdb.hg19.txt', 'http://www.circbase.org/download/hsa_hg19_circRNA.txt', False)
    #sync_databases(working_directory, 'insulators-experimental.hg19.txt', 'https://insulatordb.uthsc.edu/download/CTCFBSDB1.0/allexp.txt.gz', True)
    #sync_databases(working_directory, 'insulators-computational.hg19.txt', 'https://insulatordb.uthsc.edu/download/allcomp.txt.gz', True)

    variants = []
    gwas_snps = db_cache(f'./{working_directory}/{gene_name}_{condition_filename}_gwas_variants.json', query_gwas, [gene_name, condition])
    variants = variants + gwas_snps

    ncbi_clinical_variants = db_cache(f'./{working_directory}/{gene_name}_{condition_filename}_clinical_variants.json', get_ncbi_clinical_variants, [gene_name, condition])
    variants = variants + ncbi_clinical_variants

    print(len(variants))
    for variant in variants:
        print(variant)
    #regions = []

    #generate_structural_variants([])

    #introns = read_genecode_bed_file('intron', gene_id, f'./data/intron.bed')
    #regions = regions + introns
    #print("introns: " + str(len(introns)))
    #exons = read_genecode_bed_file('exon', gene_id, f'./data/exon.bed')
    #print("exons: " + str(len(exons)))
    #regions = regions + exons
    #utr_three_primes = read_genecode_bed_file('UTR3', gene_id, f'./data/utr3.bed')
    #print("three prime: " + str(utr_three_primes))
    #regions = regions + utr_three_primes
    #utr_five_primes = read_genecode_bed_file('UTR5', gene_id, f'./data/utr5.bed')
    #print("five prime: " + str(utr_five_primes))
    #regions = regions + utr_five_primes


    #nc_rnas = db_cache(f'./{working_directory}/{gene_name}_{condition_filename}_rna_central.json', query_rnacentral, [gene_name])
    #filtered_nc_rnas = filter_nc_rnas(chromosome, nc_rnas)
    #regions = regions + filtered_nc_rnas

    #sno_rnas = extract_sno_rnas(f'./{working_directory}/Supplementary_Dataset_S5.gff', chromosome, gene_name)
    #regions = regions + sno_rnas

    features = extract_genecode_features(f'./gencode.v37.chr_patch_hapl_scaff.annotation.gff3', gene_id)
    #print(features)
    #regions = regions + features

    #promoters = extract_promoters(f'./{working_directory}/Hs_EPDnew.bed', gene_name)
    #regions = regions + promoters

    #circular_rnas = extract_circular_rnas(f'./{working_directory}/human-circdb.hg19.txt', gene_name, chromosome)
    #regions = regions + dedup_regions(circular_rnas)

    #computational_insulators = extract_insulators(f'./{working_directory}/insulators-computational.hg19.txt', gene_name, chromosome)
    #regions = regions + dedup_regions(computational_insulators)

    #experimental_insulators = extract_insulators(f'./{working_directory}/insulators-experimental.hg19.txt', gene_name, chromosome)
    #regions = regions + dedup_regions(experimental_insulators)

    #enhancers = []
    #for enhancer_path in associated_enhancer_paths:
    #    enhancers = enhancers + extract_enhancers(working_directory, enhancer_path, gene_id, chromosome)
    #regions = regions + dedup_regions(enhancers)
    

    #unique_variant_regions = count_overlapping_variant_frequencies(regions, variants)
    #regional_frequencies, region_type_lengths = count_emperical_regional_variant_frequencies(regions, variants)
    #variant_regions = count_statistical_regional_variant_frequencies(regions, variants)

    #write_regions_to_gff3(gene_name, chromosome, regions)

    #variant_region_expected_frequencies = {}
    #total_lengths = sum(list(region_type_lengths.values()))
    #for region in region_type_lengths:
    #    length_proportion = region_type_lengths[region] / total_lengths
    #    variant_region_expected_frequencies[region] = len(variants) * length_proportion

    #regional_frequency_counts = {}
    #for region_type in regional_frequencies:
    #    if region_type not in regional_frequency_counts:
    #        regional_frequency_counts[region_type] = 0
    #    for key in list(regional_frequencies[region_type].keys()):
    #        counts = len(regional_frequencies[region_type][key]['variants'])
    #        regional_frequency_counts[region_type] += counts

    #plot_bar_chart(regional_frequency_counts, 'Regions', 'Region Frequency', 'Frequency of Variants per Genomic Region', 'variant_frequencies.png')
    #plot_bar_chart(unique_variant_regions, 'Variants', 'Overlapping Region Frequency', 'Variants in overlapping regions', 'unique_variants.png')
    #double_bar_chart(variant_regions, variant_region_expected_frequencies)

    #print('chisquare')

    #chi_square_data = [
    #    list(variant_regions.values()),
    #    list(variant_region_expected_frequencies.values())
    #]

    #stat, p, dof, expected = chi2_contingency(chi_square_data)
    
    # interpret p-value
    #alpha = 0.05
    #print('Stats: ' + stat)
    #print("Dof" + dof)
    #print("p value is " + str(p))
    #if p <= alpha:
    #    print('Dependent (reject H0)')
    #else:
    #    print('Independent (H0 holds true)')
main()
