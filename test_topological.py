
variants = [
    {
        'identifier': 'test1',
        'coordinate': 2
    },
    {
        'identifier': 'test2',
        'coordinate': 3
    },
    {
        'identifier': 'test3',
        'coordinate': 5
    }
]

regions = [
    {
        'start': 2,
        'end': 6,
        'type': 'circrna',
        'identifier': 'circrna1'
    },
    {
        'start': 5,
        'end': 10,
        'type': 'lincrna',
        'identifier': 'lincrna1'
    },
    {
        'start': 14,
        'end': 18,
        'type': 'exon',
        'identifier': 'exon1'
    },
]


def count_regional_variant_frequencies(regions, variants):
    regional_frequencies = {}
    unique_variant_regions = {}
    for variant in variants:
        counter = []
        for region in regions:
            if region['type'] not in regional_frequencies:
                regional_frequencies[region['type']] = 0
            after_start = region['start'] <= variant['coordinate']
            before_end = region['end'] >= variant['coordinate']
            print(after_start, before_end)
            if after_start and before_end:
                regional_frequencies[region['type']] += 1
                if region['type'] not in counter:
                    counter.append(region['type'])
        unique_variant_region = len(counter)
        if str(unique_variant_region) not in unique_variant_regions:
            unique_variant_regions[str(unique_variant_region)] = 0
        unique_variant_regions[str(unique_variant_region)] += 1
    return regional_frequencies, unique_variant_regions

def test_can_count():
    regional_frequency_results = {
        'exon': 0,
        'circrna': 3,
        'lincrna': 1
    }

    unique_variant_region_results = {
        '1': 2,
        '2': 1
    }

    regional_frequencies, unique_variant_regions = count_regional_variant_frequencies(regions, variants)
    assert regional_frequency_results == regional_frequencies
    assert unique_variant_region_results == unique_variant_regions



def mutate_dna(start, end, mutation, dna):
    if (end == start):
        return dna[0:start-1] + mutation + dna[end:]
    else:
        return dna[0:start-1] + mutation + dna[end+1:]

def test_deletion_mutation():
    sequence = "AAACTTGCAAAA"
    start = 5
    end = 8
    new_sequence = "T"
    m_seq = mutate_dna(start, end, new_sequence, sequence)
    #assert m_seq == "AAACTGC"
    #print(m_seq)

def test_snp_mutation():
    sequence = "GGG"
    start = 2
    end = 2
    new_sequence = "T"
    m_seq = mutate_dna(start, end, new_sequence, sequence)
    print(m_seq)

def test_insertion_mutation():
    sequence = "GGG"
    start = 2
    end = 2
    new_sequence = "TT"
    m_seq = mutate_dna(start, end, new_sequence, sequence)
    print(m_seq)

def dedup_regions(region_edges):
    sorted_region_edges = sorted(region_edges, key=lambda edge: edge['coordinate'])  
    edges = []
    stack = 0
    start = None

    for sorted_region_edge in sorted_region_edges:
        if sorted_region_edge['edge'] == 'start' and stack == 0:
            start = sorted_region_edge['coordinate']
            stack += 1
        elif sorted_region_edge['edge'] == 'start' and stack > 0:
            stack += 1
        elif sorted_region_edge['edge'] == 'end':
            stack -= 1
            if stack == 0:
                edges.append({
                    'start': start,
                    'end': sorted_region_edge['coordinate']
                })
                start = None
    return edges


def test_merge_regions_case_alpha():
    regions = [
        {'edge': 'start', 'coordinate': 3},
        {'edge': 'end', 'coordinate': 8},

        {'edge': 'start', 'coordinate': 10},
        {'edge': 'end', 'coordinate': 18},
        {'edge': 'start', 'coordinate': 12},
        {'edge': 'end', 'coordinate': 16},

        {'edge': 'start', 'coordinate': 20},
        {'edge': 'end', 'coordinate': 27},
    ]
    print(dedup_regions(regions))

def test_merge_regions_case_beta():
    regions = [
        {'edge': 'start', 'coordinate': 5},
        {'edge': 'end', 'coordinate': 9},
        {'edge': 'start', 'coordinate': 7},
        {'edge': 'end', 'coordinate': 10},

        {'edge': 'start', 'coordinate': 6},
        {'edge': 'end', 'coordinate': 12},
    ]
    print(dedup_regions(regions))


def LCSSubStr(X: str, Y: str,
                   m: int, n: int):
    # Create a table to store lengths of
    # longest common suffixes of substrings.
    # Note that LCSuff[i][j] contains length
    # of longest common suffix of X[0..i-1] and
    # Y[0..j-1]. The first row and first
    # column entries have no logical meaning,
    # they are used only for simplicity of program
    LCSuff = [[0 for i in range(n + 1)]
                 for j in range(m + 1)]
 
    # To store length of the
    # longest common substring
    length = 0
 
    # To store the index of the cell
    # which contains the maximum value.
    # This cell's index helps in building
    # up the longest common substring
    # from right to left.
    row, col = 0, 0
 
    # Following steps build LCSuff[m+1][n+1]
    # in bottom up fashion.
    for i in range(m + 1):
        for j in range(n + 1):
            if i == 0 or j == 0:
                LCSuff[i][j] = 0
            elif X[i - 1] == Y[j - 1]:
                LCSuff[i][j] = LCSuff[i - 1][j - 1] + 1
                if length < LCSuff[i][j]:
                    length = LCSuff[i][j]
                    row = i
                    col = j
            else:
                LCSuff[i][j] = 0
 
    # if true, then no common substring exists
    if length == 0:
        print("No Common Substring")
        return
 
    # allocate space for the longest
    # common substring
    resultStr = ['0'] * length
 
    # traverse up diagonally form the
    # (row, col) cell until LCSuff[row][col] != 0
    while LCSuff[row][col] != 0:
        length -= 1
        resultStr[length] = X[row - 1] # or Y[col-1]
 
        # move diagonally up to previous cell
        row -= 1
        col -= 1
 
    # required longest common substring
    return ''.join(resultStr)
 
# This code is contributed by
# sanjeev2552
import requests

#def calculate_positions(dna, dna_subset):
#    subset_length = len(dna_subset)
#    start = dna.index(dna_subset) + 1
#    end = start + subset_length - 1
#    return start, end

#dna = 'ATACGGTGAGGTAGCCGATATAGATATAGCGCTAAGGAGATAGTGCTAGACG'
#subset = 'CGATATAGATATAGCGCTAAGGAG'

#start, end = calculate_positions(dna, subset)
#print(start, end)


def order_search(coordinates, index):
    for i in range(len(coordinates)):
        print('index')
        print(i)
        print(coordinates[i]['order'])
        print(coordinates[i]['order'] == index)
        if int(coordinates[i]['order']) == int(index):
            return coordinates[i]
    return None

def remap(real_rna, coordinates):
    remapped_coordinates = []
    orders = real_rna[1:-1].split("||")

    print(orders)
    for counter in range(len(coordinates)):
        print(counter)
        order = orders[counter]
        print(order)
        coordinate = order_search(coordinates, order)
        print(coordinate)
        print(counter)
        coordinate['order'] = counter
        remapped_coordinates.append(coordinate)

    return remapped_coordinates

coordinates = [
    {'start': 9000, 'end': 9090, 'order': 1},
    {'start': 0, 'end': 30, 'order': 2},
    {'start': 3000, 'end': 3030, 'order': 3},
    {'start': 2000, 'end': 2020, 'order': 4},
]

print(remap('|4||3||1||2|', coordinates))

def get_sequence_from_ensembl(chromosome, start, end):
    server = "https://rest.ensembl.org"
    ext = f"/sequence/region/human/{chromosome}:{start}..{end}:1?coord_system_version=GRCh37"
    
    r = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"})
    return r.text

print(get_sequence_from_ensembl('chr1', 5, 8))

#location = 'chr1:-'
#data = 2	5 ,6 (7)	135, 428	0, 5936

#def coordinate_map(chromosome, genomic_start, genomic_end, actual_rna):
#    fasta_dna = get_sequence_from_ensembl(chromosome, genomic_start, genomic_end)
#    dna_string = ''.join(fasta_dna.split("\n")[1:])
#    print(dna_string)
#    print(actual_rna)
#    #processing = True
#    sub_string = LCSSubStr(actual_rna, dna_string, len(actual_rna), len(dna_string))
#    print(sub_string)
#while processing:
#    
#    print(len(sub_string))
#    if sub_string == None:
#        processing = False
#    else:
#        actual_rna = actual_rna.replace(sub_string, "-")

# Driver Code
#if __name__ == "__main__":
#    raw_sequence = 'GTCCCACCCG AAAGATGCCC CCCAGCGCCA GTGCCGTGGA CTTCTTCCAG CTCTTTGTCC CAGACAACGT CCTCAAGAAC ATGGTGGTGC AGACAAACAT GTATGCCAAG AAGTTCCAGG AGCGGTTTGG GAGCGACGGA GCCTGGGTGG AGGTGACGCT GACGGAGATG AAGGCGTTCC TGGGCTACAT GATCTCCACC AGCATCTCCC ACTGCGAGTC CGTCCTCAGC ATCTGGAGCG GAGGCTTCTA CAGCAACCGC AGCCTCGCCC TCGTCATGAG CCAGGCCCGC TTCGAGAAGA TCCTCAAGTA CTTCCACGTC GTGGCCTTCC GCTCCAGCCA GACCACGCAC GGGCTCTACA AGGTCCAGCC CTTCCTCGAC TCCCTGCAGA ACAGCTTCGA CTCTGCCTTC AGGCCTTCCC AAACCCAGGT GCTACATGAA CCCCTGATCG ATGAGGATCC TGTATTCATT GCCACGTGCA CAGAGCGGGA GCTGCGAAAG AGGAAAAAGC GGAAATTCAG CCTCTGGGTC AGACAATGTT CTTCCACTGG CTTCATCATC CAG'
#    seq = raw_sequence.replace(" ", "")
#    coordinate_map('chr1', 230486703, 230493067, seq)


#test_merge_regions_case_alpha()
#test_merge_regions_case_beta()
#test_deletion_mutation()
#test_snp_mutation()
#test_insertion_mutation()



def calculate_positions(dna, dna_subset):
    subset_length = len(dna_subset)
    start = dna.index(dna_subset) + 1
    end = start + subset_length - 1
    return start, end

print(calculate_positions('TGGGGCTCCCACCCGCCGTCCTGTTGGGGAACTGCGGAGATTCACCCCAGC', 'CGCCGTCCTGTTGGGGAACTGCG'))