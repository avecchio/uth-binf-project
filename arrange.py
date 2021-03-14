import requests
import json
import xmltodict
import ensembl_rest
import requests
import shutil
import gzip
import psycopg2
import datetime
#import matplotlib.pyplot as plt

import shutil
import urllib.request as reques        
import os.path
import wget

from contextlib import closing
from os import path

def make_working_directory(directory):
    try:
        os.makedirs(directory)
    except OSError as e:
        pass

def sync_databases(working_directory, file_name, url, unzip):
    if path.exists(f'./{working_directory}/{file_name}') == False:
        local_filename = file_name 
        local_filepath = f'./{working_directory}/{local_filename}'

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

def plot_variant_frequencies():
    pass
    #fig = plt.figure()
    #ax = fig.add_axes([0,0,1,1])
    #langs = ['C', 'C++', 'Java', 'Python', 'PHP']
    #students = [23,17,35,29,12]
    #ax.bar(langs,students)
    #plt.show()

def plot_overlapping_variants():
    pass

    # An "interface" to matplotlib.axes.Axes.hist() method
    #n, bins, patches = plt.hist(x=d, bins='auto', color='#0504aa',
    #                            alpha=0.7, rwidth=0.85)
    #plt.grid(axis='y', alpha=0.75)
    #plt.xlabel('Value')
    #plt.ylabel('Frequency')
    #plt.title('My Very Own Histogram')
    #plt.text(23, 45, r'$\mu=15, b=3$')
    #maxfreq = n.max()
    # Set a clean upper y-axis limit.
    #plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)

def arrange(variants, regions):
    #frequencies = {}
    #frequencies['unknown'] = 0

    items = []
    region_counters = {}
    variant_counters = {}

    for region in regions:
        region_counters[region['identifier']] = {
            'region': region,
            'variants': []
        }

    for variant in variants:
        unique_regions = []
        for region in regions:
            in_region = False
            if in_region:
                region['type']


#    for variant in variants:
#        items.append({
#            'name': 'variant',
#            'position': variant['start'],
#        })
#        items.append({
#            'name': 'variant',
#            'position': variant['stop'],
#        })
#    for region in regions:
#        items.append({
#            'name': region['identifier'] + '_start',
#            'position': region['start'],
#            'type': region['type']
#        })
#        items.append({
#            'name': region['identifier'] + '_end',
#            'position': region['end'],
#            'type': region['type']
#        })
    
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
    
    columns = '''
    acc.id AS acc_id,
    acc.accession AS acc_accession,
    acc.anticodon AS acc_anticodon,
    acc.classification AS acc_classification,
    acc.common_name AS acc_common_name,
    acc.db_xref AS acc_acc_db_xref,
    acc.description AS acc_description,
    acc.division AS acc_division,
    acc.experiment AS acc_experiment,
    acc.external_id AS acc_external_id,
    acc.feature_end AS acc_feature_end,
    acc.feature_name AS acc_feature_name,
    acc.feature_start AS acc_feature_start,
    acc.function AS acc_function,
    acc.gene AS acc_gene,
    acc.gene_synonym AS acc_gene_synonym,
    acc.inference AS acc_inference,
    acc.is_composite AS acc_is_composite,
    acc.keywords AS acc_keywords,
    acc.locus_tag AS acc_locus_tag,
    acc.mol_type AS acc_mol_type,
    acc.ncrna_class AS acc_ncrna_class,
    acc.non_coding_id AS acc_non_coding_id,
    acc.note AS acc_note,
    acc.old_locus_tag AS acc_old_locus_tag,
    acc.optional_id AS acc_acc_optional_id,
    acc.ordinal AS acc_ordinal,
    acc.organelle AS acc_organelle,
    acc.parent_ac AS acc_parent_ac,
    acc.product AS acc_product,
    acc.project AS acc_project,
    acc.seq_version AS acc_seq_version,
    acc.species AS acc_species,
    acc.standard_name AS acc_standard_name,
    x.id AS xref_id,
    x.created AS xref_created,
    x.last AS xref_last,
    x.upi AS xref_upi,
    x.deleted AS xref_deleted,
    x.taxid AS xref_taxid,
    x.timestamp AS xref_timestamp,
    x.userstamp AS xref_userstamp,
    x.version AS xref_version,
    x.version_i AS xref_version_i,
    r.id AS r_id,
    r.upi AS r_upi,
    r.databases AS r_databases,
    r.description AS r_description,
    r.has_coordinates AS r_has_coordinates,
    r.is_active AS r_is_active,
    r.last_release AS r_last_release,
    r.rfam_problems AS r_rfam_problems,
    r.rna_type AS r_rna_type,
    r.short_description AS r_short_description,
    r.taxid AS r_taxid,
    r.update_date AS r_update_date,
    sr.id AS sr_id,
    sr.assembly_id AS sr_assembly,
    sr.urs_taxid AS sr_urs_taxid,
    sr.chromosome AS sr_chromosome,
    sr.exon_count AS sr_exon_count,
    sr.identity AS sr_identity,
    sr.providing_databases AS sr_providing_databases,
    sr.region_name AS sr_region_name,
    sr.region_start AS sr_region_start,
    sr.region_stop AS sr_region_stop,
    sr.strand AS sr_strand,
    sr.was_mapped AS sr_was_mapped
    '''
    sql_query = f'''
    
    SELECT
    {columns}
    from rnacen.rnc_accessions acc
    LEFT JOIN rnacen.xref x on (acc.accession = x.ac)
    left join rnacen.rnc_rna_precomputed r on (x.upi = r.upi)
    left join rnacen.rnc_sequence_regions sr on (r.id = sr.urs_taxid)    

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
            results.append(list(row))
            row = cur.fetchone()
        cur.close()
        result_dictionaries = []
        return [ 
            {"acc_id,": x[0],
            "acc_accession,": x[1],
            "acc_anticodon,": x[2],
            "acc_classification,": x[3],
            "acc_common_name,": x[4],
            "acc_acc_db_xref,": x[5],
            "acc_description,": x[6],
            "acc_division,": x[7],
            "acc_experiment,": x[8],
            "acc_external_id,": x[9],
            "acc_feature_end,": x[10],
            "acc_feature_name,": x[11],
            "acc_feature_start,": x[12],
            "acc_function,": x[13],
            "acc_gene,": x[14],
            "acc_gene_synonym,": x[15],
            "acc_inference,": x[16],
            "acc_is_composite,": x[17],
            "acc_keywords,": x[18],
            "acc_locus_tag,": x[19],
            "acc_mol_type,": x[20],
            "acc_ncrna_class,": x[21],
            "acc_non_coding_id,": x[22],
            "acc_note,": x[23],
            "acc_old_locus_tag,": x[24],
            "acc_acc_optional_id,": x[25],
            "acc_ordinal,": x[26],
            "acc_organelle,": x[27],
            "acc_parent_ac,": x[28],
            "acc_product,": x[29],
            "acc_project,": x[30],
            "acc_seq_version,": x[31],
            "acc_species,": x[32],
            "acc_standard_name,": x[33],
            "xref_id,": x[34],
            "xref_created,": x[35],
            "xref_last,": x[36],
            "xref_upi,": x[37],
            "xref_deleted,": x[38],
            "xref_taxid,": x[39],
            "xref_timestamp,": x[40],
            "xref_userstamp,": x[41],
            "xref_version,": x[42],
            "xref_version_i,": x[43],
            "r_id,": x[44],
            "r_upi,": x[45],
            "r_databases,": x[46],
            "r_description,": x[47],
            "r_has_coordinates,": x[48],
            "r_is_active,": x[49],
            "r_last_release,": x[50],
            "r_rfam_problems,": x[51],
            "r_rna_type,": x[52],
            "r_short_description,": x[53],
            "r_taxid,": x[54],
            "r_update_date,": x[55],
            "sr_id,": x[56],
            "sr_assembly,": x[57],
            "sr_urs_taxid,": x[58],
            "sr_chromosome,": x[59],
            "sr_exon_count,": x[60],
            "sr_identity,": x[61],
            "sr_providing_databases,": x[62],
            "sr_region_name,": x[63],
            "sr_region_start,": x[64],
            "sr_region_stop,": x[65],
            "sr_strand,": x[66],
            "sr_was_mapped": x[67]}
        for x in results]
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
            cell_type_enhancer_name = f'{cell_type}.txt'
            sync_databases(working_directory, cell_type_enhancer_name, url, False)
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

def extract_enhancers(working_directory, file_path, gene_identifier):
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
                        'type': 'promoter',
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

def convert_coordinates(file_path):
    os.system('CrossMap.py bed hg18ToHg19.over.chain.gz test.hg18.bed')

def main():
    working_directory = 'work'
    gene_name = 'FTO'

    make_working_directory(working_directory)
    rnas = db_cache(f'./{working_directory}/rna_central.json', query_rnacentral, [gene_name])

    condition = 'Growth retardation'
    associated_enhancer_paths = sync_gene_enhancers(working_directory)
    sync_databases(working_directory, 'Hs_EPDnew.bed', 'ftp://ccg.epfl.ch/epdnew/H_sapiens/current/Hs_EPDnew.bed', False)    
   
    sync_databases(working_directory, 'gencode.v37.chr_patch_hapl_scaff.annotation.gff3', 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.chr_patch_hapl_scaff.annotation.gff3.gz', True)    
    sync_databases(working_directory, 'human-circdb.txt', 'http://www.circbase.org/download/hsa_hg19_circRNA.txt', False)
    sync_databases(working_directory, 'insulators-experimental.txt', 'https://insulatordb.uthsc.edu/download/CTCFBSDB1.0/allexp.txt.gz', True)
    sync_databases(working_directory, 'insulators-computational.txt', 'https://insulatordb.uthsc.edu/download/allcomp.txt.gz', True)

    ensemble_cds_metadata = db_cache(f'./{working_directory}/ensembl.json', get_ensembl_data, [gene_name])

    gene_id = ensemble_cds_metadata['id']

    variants = []
    gwas_snps = db_cache(f'./{working_directory}/gwas_variants.json', query_gwas, [gene_name, condition])
    #variants = variants + gwas_snps

    ncbi_clinical_variants = db_cache(f'./{working_directory}/clinical_variants.json', get_ncbi_clinical_variants, [gene_name, condition])
    #variants = variants + ncbi_clinical_variants

    #regions = []

    features = extract_genecode_features(f'./{working_directory}/gencode.v37.chr_patch_hapl_scaff.annotation.g', gene_id)
    #regions = regions + features

    promoters = extract_promoters(f'./{working_directory}/Hs_EPDnew.bed', gene_name)
    #regions = regions + promoters

    circular_rnas = extract_circular_rnas(f'./{working_directory}/human-circdb.txt', gene_name)
    #regions = regions + dedup_regions(circular_rnas)

    computational_insulators = extract_insulators(f'./{working_directory}/insulators-computational', gene_name)
    #regions = regions + dedup_regions(computational_insulators)

    experimental_insulators = extract_insulators(f'./{working_directory}/insulators-experimental', gene_name)
    #regions = regions + dedup_regions(experimental_insulators)

    enhancers = []
    for enhancer_path in associated_enhancer_paths:
        enhancers = enhancers + extract_enhancers(enhancer_path, gene_id)
    #regions = regions + dedup_regions(enhancers)

    #arrange(variants, regions)

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





