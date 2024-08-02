import cyvcf2
import argparse
import csv

def parse_arguments():
    parser = argparse.ArgumentParser(description='McDonald-Kreitman test on VCF data')
    parser.add_argument('-c', '--chromosome', required=True, help='Chromosome')
    parser.add_argument('-s', '--start', type=int, required=True, help='Start position')
    parser.add_argument('-e', '--end', type=int, required=True, help='End position')
    parser.add_argument('-n', '--name', required=True, help='Gene name')
    parser.add_argument('-v', '--vcf', required=True, help='VCF file')
    parser.add_argument('-o', '--output', help='Output file', default=None)
    return parser.parse_args()

def is_fixed(genotypes):
    """Check if all genotypes are the same (fixed) and return the fixed genotype, excluding missing genotypes.
    """
    # Exclude missing genotypes, correctly represented by 3 in cyvcf2
    non_missing_genotypes = [gt for gt in genotypes if gt != 3]
    
    if not non_missing_genotypes:
        return False, None  # All genotypes are missing
    
    # Check if all remaining genotypes are the same
    first_genotype = non_missing_genotypes[0]
    if all(gt == first_genotype for gt in non_missing_genotypes):
        return True, first_genotype  # All non-missing genotypes are the same
    else:
        return False, None  # There are different genotypes present


def main():
    args = parse_arguments()
    
    vcf_reader = cyvcf2.VCF(args.vcf)

    Pn, Ps, Dn, Ds = 0, 0, 0, 0

    for record in vcf_reader(f"{args.chromosome}:{args.start}-{args.end}"):
        genotypes = record.gt_types 

        # The last 16 samples are A. suecica and the rest are the parents accessions
        arenosa_genotypes = genotypes[:-16]
        asuecica_genotypes = genotypes[-16:]
        #print(arenosa_genotypes)
        #print(asuecica_genotypes)

        fixed_in_asuecica, asuecica_allele = is_fixed(asuecica_genotypes)
        fixed_in_arenosa, arenosa_allele = is_fixed(arenosa_genotypes)

        if fixed_in_asuecica and fixed_in_arenosa and asuecica_allele != arenosa_allele:
            if 'ANN' in record.INFO:
                for ann in record.INFO['ANN']:
                    ann_fields = ann.split('|')
                    impact = ann_fields[2]
                    if impact == 'LOW':
                        Ds += 1
                    elif impact in ['MODERATE', 'HIGH']:
                        Dn += 1

        is_polymorphic_in_suecica = any(gt > 0 for gt in asuecica_genotypes)
        if 'ANN' in record.INFO:
            for ann in record.INFO['ANN']:
                ann_fields = ann.split('|')
                impact = ann_fields[2]
                if impact == 'LOW' and is_polymorphic_in_suecica:
                    Ps += 1
                elif impact in ['MODERATE', 'HIGH'] and is_polymorphic_in_suecica:
                    Pn += 1

    pn_ps_ratio = Pn / Ps if Ps != 0 else 'undefined'
    dn_ds_ratio = Dn / Ds if Ds != 0 else 'undefined'

    output_file = args.output if args.output else f"{args.name}_MKresult.csv"
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Gene', 'Pn', 'Ps', 'Dn', 'Ds', 'Pn/Ps', 'Dn/Ds'])
        writer.writerow([args.name, Pn, Ps, Dn, Ds, pn_ps_ratio, dn_ds_ratio])

    ##print(f"MK results saved to {output_file}")

if __name__ == '__main__':
    main()
