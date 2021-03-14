import requests
import json
import xmltodict
#import ssl
#import urllib.request
import ensembl_rest
import requests
import shutil
import gzip
import psycopg2
import datetime

import shutil
import urllib.request as request
from contextlib import closing
        
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
        local_filename = file_name 
        local_filepath = f'./work/{local_filename}'

        wget.download(url, local_filepath)

        if unzip:
            unzipped_filename = local_filepath[:-3]
            with gzip.open(local_filepath, 'rb') as f_in:
                with open(unzipped_filename, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            return unzipped_filename
        return local_filepath

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

def arrange(variants, regions):
    #frequencies = {}
    #frequencies['unknown'] = 0

    items = []
    for variant in variants:
        items.append({
            'name': 'variant',
            'position': variant['start'],
        })
        items.append({
            'name': 'variant',
            'position': variant['stop'],
        })
    for region in regions:
        items.append({
            'name': region['identifier'] + '_start',
            'position': region['start'],
            'type': region['type']
        })
        items.append({
            'name': region['identifier'] + '_end',
            'position': region['end'],
            'type': region['type']
        })
    
    items = sorted(items, key=lambda item: item['position'])

    #    items.append(region)
    #items.sort(key=locus)

    #region = 'unknown'
    
    #for item in items:
    #    if item['type'] == 'variant':
    #        frequencies[region] = frequencies[region] + 1
    #    elif item['type'] == 'region':
    #        if item['position'] == 'start':
    #            region = item['name']
    #        else:
    #            region = 'unknown'
    #return frequencies
    for item in items:
        print(item)

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
        if restype == 'json':
            return json.loads(response.text)
        else:
            return xmltodict.parse(response.text)

def generate_rna_structure(rna):
    os.system(f'')

def extract_rna_features(structural_path):
    os.system(f'')

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

def query_gwas(params):
    gene_name, condition = params
    associated_snps = []
    query_url = f'https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/search/findByGene?geneName={gene_name}'
    response = requests.get(query_url)
    content = extract_content(response, 'json')
    snps = content['_embedded']['singleNucleotidePolymorphisms']
    for snp in snps:
        snpid = snp['rsId']
        is_associated = is_gwas_snp_associated(snpid, condition)
        if is_associated:
            associated_snps.append({
                'identifier': f'gwas_var{counter}',
                'coordinate': int(start),
            })
    return associated_snps

def query_ncbi(query_url, restype):
    response = requests.get(query_url)
    content = extract_content(response, restype)
    ids = extract_id_results(content)
    return ids

def query_dbsnp(gene_name, is_test):
    retmax = 20 if is_test == True else 100000 
    query_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=snp&term={gene_name}&retmax={retmax}&retmode=json'
    return query_ncbi(query_url, 'json')

def query_clinvar(gene_name, is_test):
    retmax = 20 if is_test == True else 100000 
    query_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={gene_name}&retmax={retmax}&retmode=json'
    return query_ncbi(query_url, 'json')

def query_dbvar(gene_name, condition):
    retmax = 20 if is_test == True else 100000 
    query_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=dbvar&term=%28{gene_name}%5BGene%20Name%5D%29%&retmax={retmax}&retmode=json'
    return query_ncbi(query_url)

def get_clinvar_entry(id, condition):
    query_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id={id}&from_esearch=true&retmode=json'
    res = requests.get(query_url)
    content = extract_content(res, 'xml')
    variants = []

    try:
        record = content['ClinVarResult-Set']['VariationArchive']['InterpretedRecord']
        trait_name = record['TraitMappingList']['TraitMapping']['MedGen']['@Name']
        locations = record['SimpleAllele']['Location']['SequenceLocation']

        for location in locations:
            if condition in trait_name:
                variants.append({
                    'start': int(location['@start']),
                    'stop': int(location['@stop']),
                    'reference_allele': location['@referenceAlleleVCF'],
                    'alternate_allele': location['@alternateAlleleVCF']
                })
            else:
                print(condition)
                print(trait_name)
                print('--------------')
    except:
        print('error')
    return variants


def get_dbsnp_coords(id):
    query_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&id={id}&rettype=json&retmode=text'
    res = requests.get(query_url)
    content = extract_content(res)
    snapshot_data = content['primary_snapshot_data']
    print(snapshot_data)

def query_rnacentral(params):
    gene_name = params[0]
    """ query data from the vendors table """
    conn = None

    #--accession, feature_start, feature_end, feature_name, species,
    #--external_id, anticodon, sr.chromosome, gene, gene_synonym, locus_tag,
    #--product, r.rna_type, ac, r.taxid, acc.description, r.description,
    #--region_name, strand, region_start, region_stop, exon_count
    
    sql_query = f'''
    SELECT * from rnacen.rnc_accessions acc
    where gene like '%{gene_name}%'
    '''
    t_host = "hh-pgsql-public.ebi.ac.uk"
    t_port = "5432"
    t_dbname = "pfmegrnargs"
    t_user = "reader"
    t_pw = "NWDMCE5xdipIjRrp"
    try:
        conn = psycopg2.connect(host=t_host, port=t_port, dbname=t_dbname, user=t_user, password=t_pw)
        cur = conn.cursor()
        cur.execute(sql_query)
        row = cur.fetchone()

        results = []
        while row is not None:
            results.append(row)
            row = cur.fetchone()
        cur.close()
        return results
    except (Exception, psycopg2.DatabaseError) as error:
        print(error)
        return []
    finally:
        if conn is not None:
            conn.close()

def is_locale(gene_start, gene_end, start, end, strand):
    if strand == "+":
        return int(gene_start)<int(start) and int(end)<int(gene_end)
    else:
        return int(gene_start)<int(start) and int(end)<int(gene_end)

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

def extract_circular_rnas(file_path, gene_name):
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
                        'start': int(start),
                        'end': int(end),
                        'type': 'circular_rna',
                        'strand': strand
                    })
            counter = counter + 1
    return circular_rnas

def get_coordinates(coord_str):
    chr, coordinates = coord_str.split(":")
    start, end = coordinates.split("-")
    return chr, start, end

def extract_enhancers(file_path, gene_identifier):
    enhancers = {}
    counter = 0
    with open(f'./work/{file_path}') as fp:
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
                            'start': int(start),
                            'end': int(end),
                            'type': 'enhancer'
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
                        'type': 'promoter'
                        'strand': strand
                    })
            else:
                print('not a proper entry for file: ' + file_path)
            counter = counter + 1
    return promoters

def extract_insulators(file_path, gene_name):
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
                        'start': int(start),
                        'end': int(end),
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

def extract_genecode_features(file_path, ensembl_gene_id):
    features = []
    exon_counter = 0
    five_utr = 0
    three_utr = 0
    transcript = 0
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
                
                is_type = (biotype in ['exon', 'CDS', 'three_prime_UTR', 'five_prime_UTR', 'transcript'])
                if is_type and (ensembl_gene_id in gene_id):
                    if (biotype == 'exon'):
                        exon_counter = exon_counter + 1
                    elif biotype == 'three_prime_UTR':
                        three_utr = three_utr + 1
                    elif biotype == 'five_prime_UTR':
                        five_utr = five_utr + 1
                    elif biotype == 'transcript':
                        transcript = transcript + 1
                    features.append({
                        'identifier': identifier,
                        'start': int(start),
                        'end': int(end),
                        'type': biotype
                    })

    return features

def extract_sequence(dna, start, end):
    pass

def mutate_sequence(sequence, pos, mutation):
    pass

def generate_structure(sequence):
    pass

def get_ncbi_clinical_variants(params):
    gene_name, condition = params
    variants = []
    clinical_variants = query_clinvar(gene_name, False)
    for clinical_variant in clinical_variants:
        variant_coordinates = get_clinvar_entry(clinical_variant, condition)
        variants = variants + variant_coordinates
    return variants

def gff_writer(regions):
    for region in regions:
        pass

def dedup_regions(regions):
    print(len(regions))
    unique_regions = []
    unique_region_coordinates = {}
    for region in regions:
        start = region['start']
        end = region['end']
        coordinate = f'{start}-{end}'
        if coordinate not in unique_region_coordinates:
            unique_region_coordinates[coordinate] = 0
            unique_regions.append(region)
    print(unique_region_coordinates)
    print(len(unique_regions))
    return unique_regions

def main():
    gene_name = 'FTO'
    #rnas = db_cache('rna_central.json', query_rnacentral, [gene_name])

    condition = 'Growth retardation'
    make_working_directory()
    associated_enhancer_paths = sync_gene_enhancers()
    sync_databases('Hs_EPDnew.bed', 'ftp://ccg.epfl.ch/epdnew/H_sapiens/current/Hs_EPDnew.bed', False)    
   
    sync_databases('gencode.v37.chr_patch_hapl_scaff.annotation.gff3', 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.chr_patch_hapl_scaff.annotation.gff3.gz', True)    
    sync_databases('human-circdb.txt', 'http://www.circbase.org/download/hsa_hg19_circRNA.txt', False)
    sync_databases('insulators-experimental.txt', 'https://insulatordb.uthsc.edu/download/CTCFBSDB1.0/allexp.txt.gz', True)
    sync_databases('insulators-computational.txt', 'https://insulatordb.uthsc.edu/download/allcomp.txt.gz', True)

    ensemble_cds_metadata = db_cache('ensembl.json', get_ensembl_data, [gene_name])

    gene_id = ensemble_cds_metadata['id']

    variants = []
    gwas_snps = db_cache('gwas_variants.json', query_gwas, [gene_name, condition])
    #variants = variants + gwas_snps

    ncbi_clinical_variants = db_cache('clinical_variants.json', get_ncbi_clinical_variants, [gene_name, condition])
    #variants = variants + ncbi_clinical_variants

    regions = []

    features = extract_genecode_features(f'./work/gencode.v37.chr_patch_hapl_scaff.annotation.g', gene_id)
    #regions = regions + features

    promoters = extract_promoters('./work/Hs_EPDnew.bed', gene_name)
    regions = regions + promoters

    circular_rnas = extract_circular_rnas(f'./work/human-circdb.txt', gene_name)
    #regions = regions + dedup_regions(circular_rnas)

    computational_insulators = extract_insulators(f'./work/insulators-computational', gene_name)
    #regions = regions + dedup_regions(computational_insulators)

    experimental_insulators = extract_insulators(f'./work/insulators-experimental', gene_name)
    #regions = regions + dedup_regions(experimental_insulators)

    enhancers = []
    for enhancer_path in associated_enhancer_paths:
        enhancers = enhancers + extract_enhancers(enhancer_path, gene_id)
    #regions = regions + dedup_regions(enhancers)

    arrange(variants, regions)

    #non_coding_rnas = db_cache('rna_central.json', query_rnacentral, (chromosome, gene_start, gene_end))

#    stats = {
#        'lncRNA': 0
#    }
    #features = {}
    #biotypes = {}
    #print(non_coding_rnas[0])
    #for rna in non_coding_rnas:
    #    print(rna['start'], rna['end'])
    #    if rna['biotype'] not in biotypes:
    #        biotypes[rna['biotype']] = 0
    #    biotypes[rna['biotype']] = biotypes[rna['biotype']] + 1
    #
    #    if rna['feature_type'] not in features:
    #        features[rna['feature_type']] = 0
    #    features[rna['feature_type']] = features[rna['feature_type']] + 1

main()





