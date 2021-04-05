
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

test_merge_regions_case_alpha()
test_merge_regions_case_beta()
#test_deletion_mutation()
#test_snp_mutation()
#test_insertion_mutation()