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

def sync_gene_enhancers():
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
            cell_type_enhancer_name = f'{cell_type}.txt'
            sync_databases(cell_type_enhancer_name, url, False)
            enhancer_paths.append(cell_type_enhancer_name)
        except:
            print('unable to download: ' + cell_type)

    return enhancer_paths


def sync_enhancers():
    enhancer_paths = []
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
            cell_type_enhancer_name = f'{cell_type}.bed'
            sync_databases(cell_type_enhancer_name, url, False)
            enhancer_paths.append(cell_type_enhancer_name)
        except:
            print('unable to download: ' + cell_type)

    return enhancer_paths


def extract_circular_rnas(file_path, gene_name, chromosome, start, end, is_associated):
    circular_rnas = []
    with open(file_path) as fp:
        lines = fp.readlines()
        counter = 0
        for line in lines:
            entry = line.replace("\n", "").split("\t")
            if (len(entry) == 13):
                start = entry[1]
                end = entry[2]
                is_locale = (chr == f'chr{chromosome}' and int(gene_start)<int(start) and int(end)<int(gene_end))
                if is_associated and entry[11] == gene_name:
                    circular_rnas.append({
                        'identifier': f'CircularRna{counter}',
                        'start': start,
                        'end': end,
                        'type': 'circular_rna'
                    })
                elif is_locale and (is_associated == False):
                    circular_rnas.append({
                        'identifier': f'CircularRna{counter}',
                        'start': start,
                        'end': end,
                        'type': 'circular_rna'
                    })
            counter = counter + 1
    return circular_rnas

def extract_enhancers(file_path, chromosome, gene_start, gene_end, is_associated):
    enhancers = []
    counter = 0
    with open(f'./work/{file_path}') as fp:
        lines = fp.readlines()
        for line in lines:
            entry = line.replace("\n","").split("\t")
            if (len(entry) == 4):
                chr, start, end, tf = entry
                is_locale = (chr == f'chr{chromosome}' and int(gene_start)<int(start) and int(end)<int(gene_end))
                if is_locale:
                    enhancers.append({
                        'identifier': f'Enhancer{counter}',
                        'start': start,
                        'end': end,
                        'type': 'enhancer'
                    })
                    print(entry)
                    counter = counter + 1
    return enhancers

def extract_insulators(file_path, gene_name, chromosome, gene_start, gene_end, is_associated):
    insulators = []
    counter = 0
    with open(file_path) as fp:
        lines = fp.readlines()
        for line in lines:
            if counter > 0 and len(line.split("\t")) > 4:
                entry = line.split("\t")
                identifier = entry[0]
                species = entry[1]
                locus = entry[2].split(":")

                five_prime_gene = ''
                three_prime_gene = ''

                if (len(locus) == 2):
                    five_prime_gene = entry[3]
                    three_prime_gene = entry[4]
                else:
                    locus = entry[3].split(":")
                    five_prime_gene = entry[4]
                    three_prime_gene = entry[5]

                chr = locus[0]
                coordinates = locus[1].split("-")
                start = coordinates[0]
                end = coordinates[1]

                is_locale = (chr == f'chr{chromosome}' and int(gene_start)<int(start) and int(end)<int(gene_end))
                if (gene_name == five_prime_gene or gene_name == three_prime_gene) and (is_associated):
                    print(five_prime_gene, three_prime_gene)
                    insulators.append({
                        'identifier': identifier,
                        'start': start,
                        'end': end,
                        'type': 'insulator'
                    })
                elif is_locale and (is_associated == False):
                    insulators.append({
                        'identifier': identifier,
                        'start': start,
                        'end': end,
                        'type': 'insulator'
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

def extract_genecode_features(file_path, chromosome, gene_start, gene_end):
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
                identifier = metadata[8].split(";")[0][3:]
                
                is_locale = (chr == f'chr{chromosome}' and int(gene_start)<int(start) and int(end)<int(gene_end))
                is_type = (biotype in ['exon', 'CDS', 'three_prime_UTR', 'five_prime_UTR', 'transcript'])
                if is_locale and is_type:
                    features.append({
                        'identifier': identifier,
                        'start': start,
                        'end': end,
                        'type': biotype
                    })
    return features

def does_overlap_exist(regions):
    stack = []

    coordinates = []
    for region in regions:
        coordinates.append({
            'identifier': region['identifier'],
            'coordinate': int(region['start'])
        })
        coordinates.append({
            'identifier': region['identifier'],
            'coordinate': int(region['end'])
        })
    coordinates.sort()
    for coordinate in coordinates:
        if coordinate['identifier'] == stack[-1]['identifier']:
            pass
            # remove
        else:
            pass
    print(stack)
    return len(stack) == 0


def main():
    gene_name = 'FTO'

    make_working_directory()
    non_associated_enhancer_paths = sync_enhancers()
    associated_enhancer_paths = sync_gene_enhancers()
    #sync_databases('gencode.v37.chr_patch_hapl_scaff.annotation.gff3', 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.chr_patch_hapl_scaff.annotation.gff3.gz', True)    
    sync_databases('human-circdb.txt', 'http://www.circbase.org/download/hsa_hg19_circRNA.txt', False)
    sync_databases('insulators-experimental.txt', 'https://insulatordb.uthsc.edu/download/CTCFBSDB1.0/allexp.txt.gz', True)
    sync_databases('insulators-computational.txt', 'https://insulatordb.uthsc.edu/download/allcomp.txt.gz', True)

    ensemble_cds_metadata = db_cache('ensembl.json', get_ensembl_data, ('FTO'))

    id = ensemble_cds_metadata['id']
    chromosome = ensemble_cds_metadata['seq_region_name']
    gene_start = ensemble_cds_metadata['start']
    gene_end = ensemble_cds_metadata['end']

    is_associated = True
    regions = []

    #for enhancer_path in enhancer_paths:
    #    print(enhancer_path)
    #    extract_enhancers(enhancer_path, chromosome, gene_start, gene_end)

    #features = extract_genecode_features(f'./work/gencode.v37.chr_patch_hapl_scaff.annotation.g', chromosome, gene_start, gene_end)
    #regions = regions + features

    #circular_rnas = extract_circular_rnas(f'./work/human-circdb.txt', gene_name, chromosome, gene_start, gene_end, is_associated)
    #regions = regions + circular_rnas

    #computational_insulators = extract_insulators(f'./work/insulators-computational', gene_name, chromosome, gene_start, gene_end, is_associated)
    #regions = regions + computational_insulators

    #experimental_insulators = extract_insulators(f'./work/insulators-experimental', gene_name, chromosome, gene_start, gene_end, is_associated)
    #regions = regions + experimental_insulators


    print(regions)
    print(len(regions))

#    non_coding_rnas = db_cache('rna_central.json', query_rnacentral, (chromosome, gene_start, gene_end))
#    stats = {
#        'lncRNA': 0
#    }
#    print(non_coding_rnas[0])
#    for rna in non_coding_rnas:
#        print(rna['start'], rna['end'])
#        if 'Parent' in rna:
#            print(rna['Parent'])
#        if rna['biotype'] == 'lncRNA':
#            stats['lncRNA'] = stats['lncRNA'] + 1
#        print(rna['Parent'])
#    print(stats)
#    print(len(non_coding_rnas))

#    print(non_coding_rnas)
    #regions = regions + extract_non_coding_rnas(non_coding_rnas)


#    print(regions)

#    parent_identifiers = []
#    for entry in ensemble_cds_metadata['Transcript']:
#        if 'Translation' in entry:
#            parent_identifiers.append(entry['Translation']['Parent'])
#
#    

main()





