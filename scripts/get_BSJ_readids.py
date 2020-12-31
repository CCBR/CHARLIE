import argparse

def star_parse(fusion, out):
    '''
    Parse fusion junctions from STAR aligner
    '''
    print('Start parsing fusion junctions from STAR...')
    junc = defaultdict(int)
    with open(fusion, 'r') as junc_f:
        for line in junc_f:
            flag = int(line.split()[6])
            if flag < 0:
                continue
            chr1, site1, strand1, chr2, site2, strand2 = line.split()[:6]
            if chr1 != chr2 or strand1 != strand2:
                continue
            if strand1 == '+':
                start = int(site2)
                end = int(site1) - 1
            else:
                start = int(site1)
                end = int(site2) - 1
            if start > end:
                continue
            junc_id = '%s\t%d\t%d' % (chr1, start, end)
            junc[junc_id] += 1
    total = 0
    with open(out, 'w') as outf:
        for i, j in enumerate(junc):
            outf.write('%s\tFUSIONJUNC_%d/%d\t0\t+\n' % (j, i, junc[j]))
            total += junc[j]
    print('Converted %d fusion reads!' % total)
