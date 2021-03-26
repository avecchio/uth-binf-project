
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
