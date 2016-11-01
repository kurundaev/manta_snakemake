import collections


def parseVcfInfoField(infofield):
    parts = [x.split('=') for x in infofield.split(';')]
    for i in range(len(parts)):
        if len(parts[i]) == 1:
            parts[i] = [parts[i][0], None]
    return dict(parts)


def parseVcfPRSR(format, normal, tumour):
    # Tuple order: ref, alt
    data = { 'N': { 'PR': (0, 0), 'SR': (0, 0) }, 'T': {'PR': (0, 0), 'SR': (0, 0) } }
    formatparts = format.split(':')
    normalparts = [[int(y) for y in x.split(',')] for x in normal.split(':')]
    tumourparts = [[int(y) for y in x.split(',')] for x in tumour.split(':')]
    for formati, normali, tumouri in zip(formatparts, normalparts, tumourparts):
        data['N'][formati] = normali
        data['T'][formati] = tumouri
    return data


def filterVcf(infile, outfile, args):
    """
    Filter the VCF accessed by file-like object (with a line-by-line
    iterator) infile according to the settings in args, writing the
    passed variants to the file-like object outfile.
    """
    normal_idx = None
    tumour_idx = None
    bnd_cache = {}

    tally = collections.Counter()

    for line in infile:
        # Emit header lines without modification
        if line.startswith('##'):
            outfile.write(line)
            continue

        lineparts = line.rstrip().split('\t')

        if line.startswith('#'):
            # Parse the header to identify which VCF columns
            # contain the normal and tumour genotypes.
            try:
                normal_idx = 9 + lineparts[9:].index(args.normal)
            except ValueError:
                raise ValueError('Normal sample ID {} not found in VCF'.format(args.normal))

            try:
                tumour_idx = 9 + lineparts[9:].index(args.tumour)
            except ValueError:
                raise ValueError('Tumour sample ID {} not found in VCF'.format(args.tumour))

            outfile.write(line)
            continue

        if normal_idx == None or tumour_idx == None:
            raise ValueError('Sample header line not found -- input is not in VCF format.')

        lineinfo = parseVcfInfoField(lineparts[7])

        # If requested, skip non-passing lines
        if not args.allownonpass and lineparts[6] != 'PASS':
            tally['FailNotPass'] += 1
            continue

        # If requested, skip imprecise variants
        if not args.allowimprecise and 'IMPRECISE' in lineinfo:
            tally['FailImprecise'] += 1
            continue

        # Parse the FORMAT and genotype lines to get PR and SR counts
        allele_counts = parseVcfPRSR(lineparts[8], lineparts[normal_idx], lineparts[tumour_idx])

        # Calculate depths
        depth_normal = sum([sum(allele_counts['N'][t]) for t in ['PR', 'SR']])
        depth_tumour = sum([sum(allele_counts['T'][t]) for t in ['PR', 'SR']])

        # Check depths
        if depth_normal < args.minnormaldepth or depth_tumour < args.mintumourdepth:
            tally['FailDepth'] += 1
            continue

        # Calculate VAFs
        vaf_normal = (allele_counts['N']['PR'][1] + allele_counts['N']['SR'][1]) / float(depth_normal)
        vaf_tumour = (allele_counts['T']['PR'][1] + allele_counts['T']['SR'][1]) / float(depth_tumour)

        # Check VAFs
        if vaf_normal > args.maxnormalvaf or vaf_tumour < args.mintumourvaf:
            tally['FailVAF'] += 1
            continue

        # If this line is not a BND, emit it directly
        if lineinfo['SVTYPE'] != 'BND':
            tally['Pass{}'.format(lineinfo['SVTYPE'])] += 1
            outfile.write(line)
        else:
            # This line is a BND.  It needs to be considered together with its mate,
            # as they must pass or fail together.  If we've already seen the line's
            # mate (which means that it passed), emit.  If we haven't seen this 
            # line's mate, add it to the cache for later.
            if lineinfo['MATEID'] in bnd_cache:
                outfile.write(bnd_cache[lineinfo['MATEID']])
                del bnd_cache[lineinfo['MATEID']]
                outfile.write(line)
                tally['PassBND'] += 2
            else:
                bnd_cache[lineparts[2]] = line

    tally['FailPartner'] += len(bnd_cache)

    sys.stdout.write('Filter summary\n' + '\n'.join(['{:14} {}'.format(k, v) for k, v in tally.items()]) + '\n')


if __name__ == '__main__':
    import sys
    import argparse

    def parseCommandArgs():
        parser = argparse.ArgumentParser(description = 'Filter variants in a Manta VCF to those suitable for breakend primer design.', epilog = 'Mark Pinese')

        parser_io = parser.add_argument_group('input/output arguments')
        parser_filters = parser.add_argument_group('filter arguments')

        parser_io.add_argument('--input', '-i', dest = 'input', type = str, required = False, default = '-', help = 'input Manta VCF (if -, use stdin)', metavar = 'p')
        parser_io.add_argument('--output', '-o', dest = 'output', type = str, required = False, default = '-', help = 'output filtered Manta breakend VCF (if -, use stdout)', metavar = 'p')
        parser_io.add_argument('--normal' '-n', dest = 'normal', type = str, required = False, default = 'NORMAL', help = 'Normal (non-tumour) DNA VCF sample ID', metavar = 's')
        parser_io.add_argument('--tumour' '-t', dest = 'tumour', type = str, required = False, default = 'TUMOR', help = 'Tumour DNA VCF sample ID', metavar = 's')

        parser_filters.add_argument('--transonly', action = 'store_true', dest = 'transonly', default = False, help = 'consider trans rearrangements only')
        parser_filters.add_argument('--allowimprecise', action = 'store_true', dest = 'allowimprecise', default = False, help = 'allow variants with the IMPRECISE flag')
        parser_filters.add_argument('--allownonpass', action = 'store_true', dest = 'allownonpass', default = False, help = 'allow variants without the PASS flag')
        parser_filters.add_argument('--mintumourvaf', dest = 'mintumourvaf', type = float, default = 0.3, help = 'minimum variant frequency in the tumour sample', metavar = 'f')
        parser_filters.add_argument('--maxnormalvaf', dest = 'maxnormalvaf', type = float, default = 0.0, help = 'maximum variant frequency in the normal sample', metavar = 'f')
        parser_filters.add_argument('--mintumourdepth', dest = 'mintumourdepth', type = int, default = 40, help = 'minimum total depth in the tumour sample', metavar = 'i')
        parser_filters.add_argument('--minnormaldepth', dest = 'minnormaldepth', type = int, default = 40, help = 'minimum total depth in the normal sample', metavar = 'i')

        return parser.parse_args()


    args = parseCommandArgs()

    if args.input == '-':
        infile = sys.stdin
    else:
        infile = open(args.input, 'rt')

    if args.output == '-':
        outfile = sys.stdout
    else:
        outfile = open(args.output, 'wt')


    filterVcf(infile, outfile, args)
