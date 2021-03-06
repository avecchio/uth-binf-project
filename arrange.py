import requests
import json
import xmltodict
#import ssl
#import urllib.request
import ensembl_rest
import requests
import shutil
import gzip

import os.path
from os import path

import wget

def make_working_directory():
    try:
        os.makedirs('work')
    except OSError as e:
        pass

def sync_databases(file_name, url, unzip):
    if path.exists(f'./work/{file_name}') == False:
        local_filename = file_name # if body != None else url.split('/')[-1]
        local_filepath = f'./work/{local_filename}'

        wget.download(url, local_filepath)

        if unzip:
            unzipped_filename = local_filepath[:-3]
            with gzip.open(local_filepath, 'rb') as f_in:
                with open(unzipped_filename, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            return unzipped_filename
        return local_filepath

def db_cache(file_name, callback, callback_params):
    if path.exists(file_name):
        with open(file_name) as json_file:
            data = json.load(json_file)
            return data
    else:
        data = callback(callback_params)
        with open(file_name, 'w') as outfile:
            json.dump(data, outfile)
        return data

def locus(item):
    return item['locus']

def arrange(variants, regions):
    frequencies = {}
    frequencies['unknown'] = 0

    items = []
    for variant in variants:
        items.append(variant)
        if variant['region'] not in frequencies:
            frequencies[variant['region']] = 0
    for region in regions:
        items.append(region)
    items.sort(key=locus)

    region = 'unknown'
    
    for item in items:
        if item['type'] == 'variant':
            frequencies[region] = frequencies[region] + 1
        elif item['type'] == 'region':
            if item['position'] == 'start':
                region = item['name']
            else:
                region = 'unknown'
    return frequencies

def construct_rna_variants(start, rna, variants):
    rnas = [rna]
    for variant in variants:
        s = rna[variant['coordinate'] - start] = variant['mutation']
        rnas.append(s)
    return rnas

def extract_id_results(content):
    return content['esearchresult']['idlist']

def extract_content(response, restype):
    if response.status_code != 200:
        print('error')
        return None
    else:
        #print(response.text)
        if restype == 'json':
            return json.loads(response.text)
        else:
            return xmltodict.parse(response.text)

def generate_rna_structure(rna):
    os.system(f'')

def extract_rna_features(structural_path):
    os.system(f'')

def query_gwas(gene_name):
    query_url = f'https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/search/findByGene?geneName={gene_name}'
    response = requests.get(query_url)

def query_ncbi(query_url, restype):
    response = requests.get(query_url)
    content = extract_content(response, restype)
    ids = extract_id_results(content)
    return ids

def query_dbsnp(gene_name, is_test):
    retmax = 20 if is_test == True else 100000 
    query_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=snp&term={gene_name}&retmax={retmax}&retmode=json'
    return query_ncbi(query_url, 'json')

def query_clinvar(gene_name, condition, is_test):
    retmax = 20 if is_test == True else 100000 
    query_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={gene_name}&retmax={retmax}&retmode=json'
    return query_ncbi(query_url, 'json')

def query_dbvar(gene_name, condition):
    retmax = 20 if is_test == True else 100000 
    query_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=dbvar&term=%28{gene_name}%5BGene%20Name%5D%29%&retmax={retmax}&retmode=json'
    return query_ncbi(query_url)

def get_clinvar_entry(id):
    query_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id={id}&from_esearch=true&retmode=json'
    print(query_url)
    res = requests.get(query_url)
    content = extract_content(res, 'xml')
    print(content)
    record = content['ClinVarResult-Set']['VariationArchive']['InterpretedRecord']
    print(record['ClinicalAssertionList'])

#    name = record['TraitMappingList']['TraitMapping']['MedGen']['@Name']

def get_dbsnp_coords(id):
    query_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&id={id}&rettype=json&retmode=text'
    res = requests.get(query_url)
    content = extract_content(res)
    snapshot_data = content['primary_snapshot_data']
    print(snapshot_data)

def query_gwas(gene_name, condition):
    pass

def query_rnacentral(params):
    chromosome, gene_start, gene_end = params    
    query_url = f'https://rnacentral.org/api/v1/overlap/region/homo_sapiens/{chromosome}:{gene_start}-{gene_end}'
    res = requests.get(query_url, verify=False)
    content = extract_content(res, 'json')
    return content


def extract_sequence(dna, start, end):
    pass

def mutate_sequence(sequence, pos, mutation):
    pass

def generate_structure(sequence):
    pass

def get_ensembl_data(params):
    gene_name = params[0]
    results = ensembl_rest.symbol_lookup(
        species='homo sapiens',
        symbol=gene_name,
        params={'expand': True}
    )
    return results

def sync_enhancers():
    cell_type_enhancers = [
    '786-O','A375','A549','Adipocyte','AML_blast',
    'Astrocyte','BE2C','BJ','Bronchia_epithelial','B_cell_blood',
    'C4-2','Caco-2','Calu-3','CCRF-CEM','CD133+',
    'CD14+','CD14+_monocyte','CD19+','CD20+','CD34+',
    'CD36+','CD4+','CD8+','Cerebellum','CMK',
    'CNCC','Colo320','Colo829','CUTLL1','CyT49',
    'Denditric_cell','DLD1','DOHH2','ECC-1','Endometrial_stromal_cell',
    'endometrioid_adenocarcinoma','ESC','ESC_neuron','ESC_NPC','Esophagus',
    'EWS502','Fetal_brain','Fetal_heart','Fetal_kidney','Fetal_lung',
    'Fetal_muscle_leg','Fetal_placenta','Fetal_small_intestine','Fetal_spinal_cord','Fetal_stomach',
    'Fetal_thymus','Fibroblast_foreskin','Foreskin_keratinocyte','FT246','FT33',
    'GC_B_cell','Gliobla','GM10847','GM12878','GM12891',
    'GM12892','GM18486','GM18505','GM18507','GM18508',
    'GM18516','GM18522','GM18526','GM18951','GM19099',
    'GM19141','GM19193','GM19238','GM19239','GM19240',
    'H1','H128','H2171','H54','H9',
    'HACAT','HCASMC','HCC1954','HCT116','Heart',
    'HEK293','HEK293T','HeLa-S3','Hela','Hepatocyte',
    'HepG2','HFF','HL-60','hMADS-3','HMEC',
    'hNCC','HSC','HSMM','HSMMtube','HT1080',
    'HT29','HuCCT1','HUES64','HUVEC','IMR90',
    'iPSC','Jurkat','K562','Kasumi-1','KATO3',
    'KB','KELLY','Keratinocyte','Kidney','Kidney_cortex',
    'Left_ventricle','LHCN-M2','Liver','LNCaP-1F5','LNCaP-abl',
    'LNCaP','LoVo','LP-1','LS174T','Lung',
    'LY1','Macrophage','MCF-7','MCF10A','MDA-MB-231',
    'ME-1','Melanocyte','melanoma','Mesendoderm','MM1S',
    'Monocyte','MS1','MSC_BM','Myotube','Namalwa',
    'NB4','NCCIT','NGP','NH-A','NHBE',
    'NHDF','NHEK','NHLF','NKC','NT2-D1',
    'OCI-LY1','OCI-Ly7','Osteobl','Osteoblast','Ovary',
    'P493-6','PANC-1','Pancreas','Pancreatic_islet','PBMC',
    'PC3','Plasma_cell_myeloma','PrEC','Raji','Ramos',
    'REH','Retina','RPE','RPTEC','RS4-11',
    'SEM','SGBS_adipocyte','SH-SY5Y','SK-MEL-5','SK-N-MC',
    'SK-N-SH','SK-N-SH_RA','Skeletal_muscle','SkMC','Small_intestine',
    'Sperm','Spleen','T47D-MTVL','T47D','T98G',
    'TC-797','Th1','Th2','Thymus','Treg_cell',
    'Trophoblast','U2OS','U87','Urothelial_cell','VCaP',
    'ZR75-1','ZR75-30'
    ]

    for cell_type in cell_type_enhancers:
        url = f'http://www.enhanceratlas.org/data/download/enhancer/hs/{cell_type}.bed'
        try:
            sync_databases(f'{cell_type}.bed', url, False)
        except:
            print('unable to download: ' + cell_type)

def extract_circular_rnas(file_path, identifiers, start, end):
    with open(file_path) as fp:
        lines = fp.readlines()
        for line in lines:
            if len(line.split("\t")) == 4:
                identifier, start, end, tf, = line.split("\t")
                if identifier.split(".")[0] in identifiers:
                    print(line)

def extract_genecode_features(file_path, chromosome, start, end):
    with open(file_path) as fp:
        lines = fp.readlines()
        for line in lines:
            pass

def extract_enhancers(path_to_beds, chromosome, start, end):
    pass

def extract_insulators(file_path, chromosome, start, end):
    with open(file_path) as fp:
        lines = fp.readlines()
        for line in lines:
            pass

def main():
    gene_name = 'FTO'

    make_working_directory()
    sync_enhancers()
    sync_databases('gencode.v37.chr_patch_hapl_scaff.annotation.gff3', 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.chr_patch_hapl_scaff.annotation.gff3.gz', True)    
    sync_databases('human-circdb.txt', 'http://www.circbase.org/download/hscore/hscores_human_gencode27.tar.gz', True)
    #sync_databases('insulators-experimental.txt', 'https://insulatordb.uthsc.edu/download/allexp.txt.gz', True)
    sync_databases('insulators-computational.txt', 'https://insulatordb.uthsc.edu/download/allcomp.txt.gz', True)

    ensemble_cds_metadata = db_cache('ensembl.json', get_ensembl_data, ('FTO'))

    id = ensemble_cds_metadata['id']
    chromosome = ensemble_cds_metadata['seq_region_name']
    gene_start = ensemble_cds_metadata['start']
    gene_end = ensemble_cds_metadata['end']

    #print(id)

    parent_identifiers = []
    for entry in ensemble_cds_metadata['Transcript']:
        if 'Translation' in entry:
            parent_identifiers.append(entry['Translation']['Parent'])

    extract_circular_rnas(f'./work/human-circdb', parent_identifiers, gene_start, gene_end)

    #rna_central_rnas = db_cache('rna_central.json', query_rnacentral, (chromosome, gene_start, gene_end))
    #print(len(rna_central_rnas))

#query_clinvar(gene_name, '', True)
#ncbi_snps = query_dbsnp(gene_name, True)
#print(ncbi_snps)
#for ncbi_snp in ncbi_snps:
#    pass

#get_clinvar_entry('1599393422')

#get_dbsnp_coords(268)
#query_ncbi('FTO')

main()





