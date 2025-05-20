# Code construit Par Felix-Antoine

import pysam
import os
import numpy as np
import pickle
import pyfaidx
from multiprocessing import Pool
from skbio import TabularMSA, DNA
from skbio.alignment import local_pairwise_align_ssw, AlignmentStructure
from math import log
import pandas as pd
import argparse

def load_exon_data(parquet_path):
    df = pd.read_parquet(parquet_path)

    # exon_dict: transcript_id -> list of (chrom, start, end)
    exon_dict = {}
    for tid, grp in df.groupby("transcript_id"):
        exon_dict[tid] = [(row.chromosome, row.start, row.end) for row in grp.itertuples()]

    # transcripts_dict: gene_name -> list of transcript_ids
    transcripts_dict = df.groupby("gene_name")["transcript_id"].unique().apply(list).to_dict()

    # gene_exons_dict: gene_name -> {chrom, strand, exons: [{start,end}, ...]}
    gene_exons_dict = {}
    for gene, grp in df.groupby("gene_name"):
        gene_exons_dict[gene] = {
            "chrom": grp.chromosome.iloc[0],
            "strand": grp.strand.iloc[0],
            "exons": [{"start": row.start, "end": row.end} for row in grp.itertuples()]
        }
    return gene_exons_dict, exon_dict, transcripts_dict

# Get gene coords for remapping
def get_gene_coords_dict(gene_exons_dict):
    gene_coords = dict()
    for gene_id, coords in gene_exons_dict.items():
        gene_coords[gene_id] = {
            'chrom':coords['chrom'],
            'start':min(x['start'] for x in coords['exons']),
            'end':max(x['end'] for x in coords['exons']),
            'strand':coords['strand']
            }
    return gene_coords

# Get reads_stats
# 'bamfile'= contain the aligned reads and 'bamfile_disc'= discordant reads
def get_reads_stats(bamfile, bamfile_disc):
    # Remove all the reads with a length <1
    insert_lengths = [r.tlen for r in bamfile.fetch() if r.tlen>1]
    # Sort them by ascending order
    insert_lengths = sorted(insert_lengths)
    # calculate the index of the 99.5th percentile of the insert lengths 
    idx = int(np.ceil(len(insert_lengths)*0.995))
    # max insert size for disc reads
    insert995 = insert_lengths[idx]

    # Dict of differents stats for reads and discordant reads
    return {
        'mean insert len': np.mean(insert_lengths),
        'median insert len': np.median(insert_lengths),
        'insert995': insert995,
        'disc_read_lengths': [r.tlen for r in bamfile_disc.fetch() if r.tlen>1],
        'good_reads':bamfile.count(),
        'disc_reads':bamfile_disc.count()
    }

def create_reads_stats(reads_stats_path, bamfile, bamfile_disc):
    if os.path.exists(reads_stats_path):
        reads_stats = pickle.load(open(reads_stats_path, 'rb'))
        if 'insert995' not in reads_stats:
            reads_stats = get_reads_stats(bamfile, bamfile_disc)
            pickle.dump(reads_stats, open(reads_stats_path, 'wb'))
        else:
            reads_stats = pickle.load(open(reads_stats_path, 'rb'))
        
    else:
        reads_stats = get_reads_stats(bamfile, bamfile_disc)
        pickle.dump(reads_stats, open(reads_stats_path, 'wb'))
    return reads_stats

def sequence_to_repvec(sequence, N):
    """
    Computes the repetition vector (as seen in Wooton, 1993) from a
    given sequence of a biopolymer with `N` possible residues.

    :param sequence: the nucleotide or protein sequence to generate a repetition vector for.
    :param N: the total number of possible residues in the biopolymer `sequence` belongs to.
    """
    encountered_residues = set()
    repvec = []
    for residue in sequence:
        if residue not in encountered_residues:
            residue_count = sequence.count(residue)
            repvec.append(residue_count)
            encountered_residues.add(residue)

        if len(encountered_residues) == N:
            break

    while len(repvec) < N:
        repvec.append(0)

    return sorted(repvec, reverse=True)

def sequence_entropy(sequence, N=4):
    repvec = sequence_to_repvec(sequence, N)
    L = len(sequence)
    entropy = sum([-1*(n/L)*log((n/L), N) for n in repvec if n != 0])  # formule de l'entropie de Shannon
    return entropy

# Filter DRPs for remapping
def filtering_DRPs(samfile_disc, insert995, entropy_min=0.7):
    done_reads = []
    good_pairs = []
    
    gene_good_read = [read for read in samfile_disc.fetch() 
                       if np.abs(read.tlen) > insert995 and read.is_paired and not read.mate_is_unmapped and not read.is_duplicate and read.is_read1]
    for ra in gene_good_read:
        rb = samfile_disc.mate(ra)
        if ra not in done_reads and rb not in done_reads:
            if sequence_entropy(ra.seq)>=entropy_min and sequence_entropy(rb.seq)>=entropy_min:
                good_pairs.append((ra,rb))
        done_reads.extend((ra, rb))

    return good_pairs

def get_reverse_gene_coords(gene_coords):
    reverse_gene_coords = {}
    for gene in gene_coords:
        if gene_coords[gene]['chrom'] in reverse_gene_coords:
            reverse_gene_coords[gene_coords[gene]['chrom']].update({(gene_coords[gene]['start'], gene_coords[gene]['end'], gene_coords[gene]['strand']):gene})
        else:
            reverse_gene_coords.update({gene_coords[gene]['chrom'] : {(gene_coords[gene]['start'], gene_coords[gene]['end'], gene_coords[gene]['strand']):gene}})
    return reverse_gene_coords

def get_possible_genes(reverse_gene_coords, chr, start, end):
    possible_gene = []
    for gene_coords in reverse_gene_coords[chr]:
        if start in range(gene_coords[0], gene_coords[1]) or end in range(gene_coords[0], gene_coords[1]):
                possible_gene.append(reverse_gene_coords[chr][gene_coords])
    return possible_gene

def find_overlap(possible_gene, gene_coords_dict, read):
    genes_strand = [gene_coords_dict[gene]['strand'] for gene in possible_gene]
    if read.is_forward == True:
        good_strand = '+'
    elif read.is_reverse == True:
        good_strand = '-'
    good_indexs = [i for i, x in enumerate(genes_strand) if x == good_strand]
    if len(good_indexs) == 1:
        return possible_gene[good_indexs[0]]
    else:
        return possible_gene

def get_gene(gene_coords_dict, reverse_gene_coords, drp, samfile_disc):
    r1, r2 = drp
    chr1 = samfile_disc.get_reference_name(r1.reference_id)
    start1 = r1.reference_start + 1
    end1 = r1.reference_end
    chr2 = samfile_disc.get_reference_name(r2.reference_id)
    start2 = r2.reference_start + 1
    end2 = r2.reference_end
    chroms = set([str(x) for x in range(1, 23)] + ['X', 'Y', 'MT'])
    
    if chr1 in chroms:
        possible_gene_1 = get_possible_genes(reverse_gene_coords, chr1, start1, end1)
    elif chr1 not in chroms:
        possible_gene_1 = []

    if chr2 in chroms:
        possible_gene_2 = get_possible_genes(reverse_gene_coords, chr2, start2, end2)
    elif chr2 not in chroms:
        possible_gene_2 = []


    if len(possible_gene_1) == 1:
        gene1 = possible_gene_1[0]     
    elif len(possible_gene_1) == 0:
        gene1 = (chr1, start1, end1)       
    elif len(possible_gene_1) > 1:
        gene1 = find_overlap(possible_gene_1, gene_coords_dict, drp[0])

    
    if len(possible_gene_2) == 1:
        gene2 = possible_gene_2[0]   
    elif len(possible_gene_2) == 0:
        gene2 = (chr2, start2, end2)
    elif len(possible_gene_2) > 1:
        gene2 = find_overlap(possible_gene_2 ,gene_coords_dict, drp[1])

    return gene1, gene2

def get_reads_gene_dict(top_notch_DRP, gene_coords_dict, samfile_disc):
    reads_gene_dict = {}
    reverse_gene_coords = get_reverse_gene_coords(gene_coords_dict)
    for drp in top_notch_DRP:
        gene1, gene2 = get_gene(gene_coords_dict, reverse_gene_coords, drp, samfile_disc)
        if drp[0] in reads_gene_dict: 
            print(f"{drp[0]} duplicate", flush=True)
        if drp[1] in reads_gene_dict:
            print(f"{drp[1]} duplicate", flush=True)
        reads_gene_dict.update({drp[0]:gene1, drp[1]:gene2})

    return reads_gene_dict

def separation(top_notch_DRP, reads_gene_dict, output_path):
    same_gene_DRP = []
    diff_genes_DRP = []
    both_intergenic_DRP = []
    one_intergenic_DRP = []
    too_much_genes = []
    
    for (r1, r2) in top_notch_DRP:
        if isinstance(reads_gene_dict[r1], list) or isinstance(reads_gene_dict[r2], list):
            too_much_genes.append((r1.to_dict(), r2.to_dict()))
        elif isinstance(reads_gene_dict[r1], tuple) and isinstance(reads_gene_dict[r2], tuple):
            both_intergenic_DRP.append((r1.to_dict(), r2.to_dict()))
        elif isinstance(reads_gene_dict[r1], tuple) or isinstance(reads_gene_dict[r2], tuple):
            one_intergenic_DRP.append((r1.to_dict(), r2.to_dict()))
        elif reads_gene_dict[r1] == reads_gene_dict[r2]:
            same_gene_DRP.append((r1, r2))
        elif reads_gene_dict[r1] != reads_gene_dict[r2]:
            diff_genes_DRP.append((r1, r2))

    pickle.dump(both_intergenic_DRP, open(f'{output_path}/both_intergenic_DRP.pkl', 'wb'))
    pickle.dump(one_intergenic_DRP, open(f'{output_path}/one_intergenic_DRP.pkl', 'wb'))
    pickle.dump(too_much_genes, open(f'{output_path}/too_much_genes.pkl', 'wb'))
    
    return same_gene_DRP, diff_genes_DRP

def get_insert_size(start_end_positions_1, start_end_positions_2): # Insert size include reads size (+300 nt)
    (start_1, end_1) = start_end_positions_1
    (start_2, end_2) = start_end_positions_2

    if start_2 > start_1:
        insert_size = int(np.abs(end_2-start_1)) # insert size after re-mapping rb in same gene as ra

    elif start_2 < start_1:
        insert_size = int(np.abs(end_1-start_2)) # insert size after re-mapping rb in same gene as ra

    elif start_2 == start_1: # and end_1 == end_2: # voir ce qu'on fait avec ceux là      ÉGALE AU START OU ÉGALE AU END
        #insert_size = int(np.abs(end_1-start_2))   IDÉE À VOIR SI OK, ON POURRAIT ALLER VOIR LA OVERLAPPING REGION
        insert_size = 'EQUAL'

    return insert_size

# Remapping with SSW
def exon_remap(transcripts_dict, gene, exon_dict, hg38_genome, rb_seq):
    
    exons = [] # une liste qui contiendra tous les exons concaténés (donc les 10 transcripts alternatifs)
    for transcript in transcripts_dict[gene]:
        exon_seq = ""
        for exon in exon_dict[transcript]:
            exon_seq += hg38_genome[exon[0]][exon[1]:exon[2]] # concaténation
        exons.append(exon_seq) # ajout de l'exon concaténé (transcript) à la liste
        
    best_score = 0
    for i in range(len(exons)):  # on trouve le transcripts (exon concat) avec le meilleurs score pour la séquence selon ssw
        new_score = local_pairwise_align_ssw(DNA(rb_seq),DNA(exons[i]))[1] #local_pairwise_align_ssw fonction de skbio DNA fait partie de la fonction
        if new_score > best_score:
            best_score = new_score
            index_exons = i
    
    pos_exons = local_pairwise_align_ssw(DNA(rb_seq),DNA(exons[index_exons]))[2][1]  # on récupère la position de l'alignement dans l'exon concaténé
    
    exon_interet = exon_dict[transcripts_dict[gene][index_exons]]  # on recupère la liste des exons qui composent l'exon concaténé d'interet

    exon_cc = [0]  # cette liste contiendra les position des exons composant l'exon concaténé d'interet , on saura quels exons exactement s'aligne avec ssw
    for i in range(len(exon_interet)):
        exon_cc.append(exon_interet[i][2]-exon_interet[i][1]+exon_cc[i])

    for i in range(len(exon_cc)):  # on récupère le numéro des exons qui s'aligne avec le read selon ssw
        if (exon_cc[i] > pos_exons[0]):
            for j in range(i,len(exon_cc)):
                if (exon_cc[j] > pos_exons[1]):
                    i = i - 1
                    j = j - 1 # la liste commencait avec un 0 donc on ajuste avec -1
                    break
            break
            
    if i == j:
        start = exon_interet[i][2] - (exon_cc[i+1]-pos_exons[0]) + 1 # on récupère la position génomique de la fin d'exon - le petit bout calculé entre l'alignement et la coordonne
        end = exon_interet[i][2] - (exon_cc[i+1]-pos_exons[1]) + 1
        start_end = (start,end)
    else:
        start = exon_interet[i][2] - (exon_cc[i+1]-pos_exons[0]) + 1 # on récupère la position génomique de la fin d'exon - le petit bout calculé entre l'alignement et la coordonne
        end = exon_interet[j][1] + (pos_exons[1] - exon_cc[j]) + 1
        start_end = (start,end)
        
    return best_score,start_end

def remap_same_gene(same_gene_DRP, reads_gene_dict, hg38_genome, transcripts_dict, exon_dict, gene_coords_dict):
    same_gene_remap = []
    for r1, r2 in same_gene_DRP:
        gene_id = reads_gene_dict[r1]
        coords = gene_coords_dict[gene_id]
        
        align_1, score_orig_1, start_end_positions_orig_1 = local_pairwise_align_ssw(DNA(hg38_genome[coords['chrom']][coords['start']:coords['end']]), DNA(r1.query_sequence))
        align_2, score_orig_2, start_end_positions_orig_2 = local_pairwise_align_ssw(DNA(hg38_genome[coords['chrom']][coords['start']:coords['end']]), DNA(r2.query_sequence))

        orig_insert = get_insert_size(start_end_positions_orig_1[0], start_end_positions_orig_2[0]) # start_end_positions_orig_1 is mapping coords on gene [0] (start, end) and on read [1] (start, end) 
        if orig_insert == 'EQUAL' : continue
        
        score_exon_1, start_end_positions_exon_1 = exon_remap(transcripts_dict, gene_id, exon_dict, hg38_genome, r1.query_sequence)
        score_exon_2, start_end_positions_exon_2 = exon_remap(transcripts_dict, gene_id, exon_dict, hg38_genome, r2.query_sequence)

        remap_2_exon_insert = get_insert_size((r1.reference_start, r1.reference_end), start_end_positions_exon_2) # insert size if read 2 remapped on exon concat
        if remap_2_exon_insert == 'EQUAL' : continue # voir ce qu'on fait avec ceux là
        
        remap_1_exon_insert = get_insert_size((r2.reference_start, r2.reference_end), start_end_positions_exon_1) # insert size if read 2 remapped on exon concat
        if remap_1_exon_insert == 'EQUAL' : continue # voir ce qu'on fait avec ceux là

        remap_both_exon_concat_insert = get_insert_size(start_end_positions_exon_2, start_end_positions_exon_1) # insert size between remap of both genes on exon concat
        if remap_both_exon_concat_insert == 'EQUAL' : continue # voir ce qu'on fait avec ceux là
        
        # read to dict
        r1 = r1.to_dict()
        r1.update({
            'score_orig' : score_orig_1,
            'start_end_positions_orig' : start_end_positions_orig_1,
            'orig_insert_size' : orig_insert + 1,  # Insert size between the two original map (same in r2 dict)
            'score_exon' : score_exon_1,
            'start_end_positions_exon' : start_end_positions_exon_1,
            'remap_exon_insert' : remap_1_exon_insert +1,  # Insert size between r1 remap on exon gene 1 and original map of r2
            'remap_both_exon_insert' : remap_both_exon_concat_insert +1, # Insert size between r1 and r2 remap on exon concat gene
            'gene':gene_id
            })

        r2 = r2.to_dict()
        r2.update({
            'score_orig' : score_orig_2,
            'start_end_positions_orig' : start_end_positions_orig_2,
            'orig_insert_size' : orig_insert + 1,  # Insert size between the to original map (same in r1 dict)
            'score_exon' : score_exon_2,
            'start_end_positions_exon' : start_end_positions_exon_2,
            'remap_exon_insert' : remap_2_exon_insert +1,  # Insert size between r2 remap on exon gene 2 and original map of r1
            'remap_both_exon_insert' : remap_both_exon_concat_insert +1, # Insert size between r1 and r2 remap on exon concat gene
            'gene':gene_id
            })
        
        same_gene_remap.append((r1, r2))
    return same_gene_remap

def remap_diff_gene(diff_genes_DRP, reads_gene_dict, hg38_genome, transcripts_dict, exon_dict, gene_coords_dict):
    diff_gene_remap = []
    for r1, r2 in diff_genes_DRP:
        gene_1 = reads_gene_dict[r1]
        gene_2 = reads_gene_dict[r2]
        coords_1 = gene_coords_dict[gene_1]
        coords_2 = gene_coords_dict[gene_2]

        
# SSW MAP ON ORIGINAL GENE
        align_orig_1, score_orig_1, start_end_positions_orig_1 = local_pairwise_align_ssw(DNA(hg38_genome[coords_1['chrom']][coords_1['start']:coords_1['end']]), DNA(r1.query_sequence))
        align_orig_2, score_orig_2, start_end_positions_orig_2 = local_pairwise_align_ssw(DNA(hg38_genome[coords_2['chrom']][coords_2['start']:coords_2['end']]), DNA(r2.query_sequence))

        orig_insert = get_insert_size(start_end_positions_orig_1[0], start_end_positions_orig_2[0]) # start_end_positions_orig is mapping coords on gene [0] (start, end) and on read [1] (start, end) 
        if orig_insert == 'EQUAL' : continue

        
# SSW MAP ON MATE GENE
        align_mate_1, score_mate_1, start_end_positions_mate_1 = local_pairwise_align_ssw(DNA(hg38_genome[coords_2['chrom']][coords_2['start']:coords_2['end']]), DNA(r1.query_sequence))
        align_mate_2, score_mate_2, start_end_positions_mate_2 = local_pairwise_align_ssw(DNA(hg38_genome[coords_1['chrom']][coords_1['start']:coords_1['end']]), DNA(r2.query_sequence))
        
        mate_insert_1 = get_insert_size((r2.reference_start, r2.reference_end), start_end_positions_mate_1[0]) # start_end_positions_mate is mapping coords on gene [0] (start, end) and on read [1] (start, end) 
        if mate_insert_1 == 'EQUAL' : continue

        mate_insert_2 = get_insert_size((r1.reference_start, r1.reference_end), start_end_positions_mate_2[0])
        if mate_insert_2 == 'EQUAL' : continue

        
# SSW MAP ON ORIGINAL GENE EXON CONCAT 
        score_exon_1, start_end_positions_exon_1 = exon_remap(transcripts_dict, gene_1, exon_dict, hg38_genome, r1.query_sequence)
        score_exon_2, start_end_positions_exon_2 = exon_remap(transcripts_dict, gene_2, exon_dict, hg38_genome, r2.query_sequence)

        remap_2_exon_insert = get_insert_size((r1.reference_start, r1.reference_end), start_end_positions_exon_2) # Insert size between r2 remap on exon gene 2 and original map of r1
        if remap_2_exon_insert == 'EQUAL' : continue # voir ce qu'on fait avec ceux là
        
        remap_1_exon_insert = get_insert_size((r2.reference_start, r2.reference_end), start_end_positions_exon_1) # Insert size between r1 remap on exon gene 1 and original map of r2
        if remap_1_exon_insert == 'EQUAL' : continue # voir ce qu'on fait avec ceux là

        
# SSW MAP ON MATE GENE EXON CONCAT 
        score_exon_mate_1, start_end_positions_exon_mate_1 = exon_remap(transcripts_dict, gene_2, exon_dict, hg38_genome, r1.query_sequence)
        score_exon_mate_2, start_end_positions_exon_mate_2 = exon_remap(transcripts_dict, gene_1, exon_dict, hg38_genome, r2.query_sequence)

        remap_2_exon_mate_insert = get_insert_size((r1.reference_start, r1.reference_end), start_end_positions_exon_mate_2) # Insert size between r2 remap on exon gene 1 and original map of r1
        if remap_2_exon_mate_insert == 'EQUAL' : continue # voir ce qu'on fait avec ceux là
        
        remap_1_exon_mate_insert = get_insert_size((r2.reference_start, r2.reference_end), start_end_positions_exon_mate_1) # Insert size between r1 remap on exon gene 2 and original map of r2
        if remap_1_exon_mate_insert == 'EQUAL' : continue # voir ce qu'on fait avec ceux là

        
# read to dict
        r1 = r1.to_dict()
        r1.update({
            'score_orig' : score_orig_1,
            'start_end_positions_orig' : start_end_positions_orig_1,
            'orig_insert_size' : orig_insert + 1,  # Insert size between the to original map (same in r2 dict)
            'score_mate' : score_mate_1,
            'start_end_positions_mate' : start_end_positions_mate_1,
            'mate_insert_size' : mate_insert_1 +1,  # Insert size between r1 remap on gene 2 and original map of r2
            'score_exon' : score_exon_1,
            'start_end_positions_exon' : start_end_positions_exon_1,
            'exon_insert_size' : remap_1_exon_insert +1,  # Insert size between r1 remap on exon gene 1 and original map of r2
            'score_exon_mate' : score_exon_mate_1,
            'start_end_positions_exon_mate' : start_end_positions_exon_mate_1,
            'exon_mate_insert_size' : remap_1_exon_mate_insert +1,  # Insert size between r1 remap on exon gene 2 and original map of r2
            'gene_orig' : gene_1,
            'gene_mate' : gene_2
            })

        r2 = r2.to_dict()
        r2.update({
            'score_orig' : score_orig_2,
            'start_end_positions_orig' : start_end_positions_orig_2,
            'orig_insert_size' : orig_insert + 1,  # Insert size between the to original map (same in r1 dict)
            'score_mate' : score_mate_2,
            'start_end_positions_mate' : start_end_positions_mate_2,
            'mate_insert_size' : mate_insert_2 + 1,  # Insert size between r2 remap on gene 1 and original map of r1gene_coords
            'score_exon' : score_exon_2,
            'start_end_positions_exon' : start_end_positions_exon_2,
            'exon_insert_size' : remap_2_exon_insert + 1,  # Insert size between r2 remap on exon gene 2 and original map of r1
            'score_exon_mate' : score_exon_mate_2,
            'start_end_positions_exon_mate' : start_end_positions_exon_mate_2,
            'exon_mate_insert_size' : remap_2_exon_mate_insert +1,  # Insert size between r2 remap on exon gene 1 and original map of r1
            'gene_orig' : gene_2,
            'gene_mate' : gene_1
            })
        
        diff_gene_remap.append((r1, r2))
    return diff_gene_remap

def running_remapping(bamfile_path, bamfile_disc_path, genome_fasta_path, gene_exons_dict, exon_dict, transcripts_dict):
    # Set variables
    exp_path = '/'.join(bamfile_path.split('/')[:-1])
    exp_name = bamfile_path.split('/')[-1]
    output_path = f'{exp_path}/pipeline_output'

    print(f'starting {exp_name}...', flush=True)

    bamfile = pysam.AlignmentFile(bamfile_path, "rb")
    bamfile_disc = pysam.AlignmentFile(bamfile_disc_path, "rb")
    genome_fasta_path = pyfaidx.Fasta(genome_fasta_path, as_raw=True)

    # Create outpath
    if not os.path.isdir(output_path):
        os.makedirs (output_path)
        print(f'{output_path} dir created', flush=True)

    # Get gene_coords
    gene_coords_dict = get_gene_coords_dict(gene_exons_dict)
    print(f'{exp_name} gene_coords_dict done', flush=True)
    
    # Get reads_stats
    reads_stats_path = f'{output_path}/reads_stats.pkl'
    reads_stats = create_reads_stats(reads_stats_path, bamfile, bamfile_disc)
    print(f'{exp_name} reads_stats done', flush=True)

    # Top Notch list
    top_notch_DRP = filtering_DRPs(bamfile_disc, reads_stats['insert995'], entropy_min=0.7)
    # pickle.dump(reads_gene_dict, open(f'{output_path}{ex_id}/reads_gene.pkl', 'wb')) # Cant pickle the parsing object 'read'
    print(f'{exp_name} Filtration done', flush=True)

    # Get reads_gene_dict
    reads_gene_dict = get_reads_gene_dict(top_notch_DRP, gene_coords_dict, bamfile_disc)
    print(f'{exp_name} reads_gene_dict done', flush=True)

    # Separation
    same_gene_DRP, diff_gene_DRP = separation(top_notch_DRP, reads_gene_dict, output_path)
    print(f'{exp_name} Separation done', flush=True)

    # Remapping same gene
    remapped_same_gene = remap_same_gene(same_gene_DRP, reads_gene_dict, genome_fasta_path, transcripts_dict, exon_dict, gene_coords_dict)
    pickle.dump(remapped_same_gene, open(f'{output_path}/remapped_same_gene.pkl', 'wb'))
    print(f'{exp_name} remapped_same_gene done', flush=True)

    # Remapping diff gene
    remapped_diff_gene = remap_diff_gene(diff_gene_DRP, reads_gene_dict, genome_fasta_path, transcripts_dict, exon_dict, gene_coords_dict)
    pickle.dump(remapped_diff_gene, open(f'{output_path}/remapped_diff_gene.pkl', 'wb'))
    print(f'{exp_name} remapped_diff_gene done', flush=True)

    print(f'{exp_name} done!', flush=True)

if __name__ == '__main__':
    # Snakemake met cet objet à dispo
    bamfile           = snakemake.input.bam
    bamfile_disc      = snakemake.input.disc_bam
    exon_parquet      = snakemake.input.genome_dict
    genome_fasta_path = snakemake.input.genome

    out_same = snakemake.output.same_gene_remap
    out_diff = snakemake.output.diff_gene_remap

    # 1) charge vos dictionnaires à partir du parquet
    gene_exons_dict, exon_dict, transcripts_dict = load_exon_data(exon_parquet)

    # 2) ouvre les BAM
    bam   = pysam.AlignmentFile(bamfile,      "rb")
    disc  = pysam.AlignmentFile(bamfile_disc, "rb")
    hg38  = pyfaidx.Fasta(genome_fasta_path, as_raw=True)

    # 3) calcule vos stats et DRPs
    stats       = create_reads_stats(f"{os.path.dirname(bamfile)}/reads_stats.pkl", bam, disc)
    top_notch   = filtering_DRPs(disc, stats['insert995'])
    rgd         = get_reads_gene_dict(top_notch, get_gene_coords_dict(gene_exons_dict), disc)
    same_DRP, diff_DRP = separation(top_notch, rgd, os.path.dirname(bamfile))

    # 4) remappage
    remapped_same = remap_same_gene(same_DRP, rgd, hg38, transcripts_dict, exon_dict, get_gene_coords_dict(gene_exons_dict))
    remapped_diff = remap_diff_gene(diff_DRP, rgd, hg38, transcripts_dict, exon_dict, get_gene_coords_dict(gene_exons_dict))

    # 5) écriture des deux résultats finaux
    pickle.dump(remapped_same, open(out_same, 'wb'))
    pickle.dump(remapped_diff, open(out_diff, 'wb'))
