#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 10:15:59 2018

@author: bogao
"""

import gzip
import os
from pymongo import MongoClient

biosamples_dict = {}
individuals_dict = {}
variants_dict = {}
callsets_dict = {}


def parse_cnvgz(filepath):
    global biosamples_dict, individuals_dict, variants_dict, callsets_dict
    with gzip.open(filepath, 'rt') as fin:
        for i in range(6):
            next(fin)
        for line in fin:
            line = line.strip().split('\t')
            chro = line[4]
            start = line[5]
            end = line[6]
            strand = line[7]
#            variant_classification = line[8]
            variant_type = line[9]
            ref_allele = line[10]
            tumor_allele_1 = line[11]
            tumor_allele_2 = line[12]
            tumor_sample_barcode = line[15]
            normal_sample_barcode = line[16]
            
            biosample_id_tumor = tumor_sample_barcode
            biosample_id_normal = normal_sample_barcode
            individual_id = tumor_sample_barcode[0:12]
            variant_id = '{}_{}_{}_{}'.format(chro,start,end,strand)
            callset_id = '{}_{}'.format(tumor_sample_barcode, normal_sample_barcode)
            
            if biosample_id_tumor not in biosamples_dict:
                biosamples_dict[biosample_id_tumor] = {
                        'id': biosample_id_tumor,
                        'name': biosample_id_tumor,
                        'description': 'Tumor sample',
                        'individual_id': individual_id}
            if biosample_id_normal not in biosamples_dict:
                biosamples_dict[biosample_id_normal] = {
                        'id': biosample_id_normal,
                        'name': biosample_id_normal,
                        'description': 'Normal sample',
                        'individual_id': individual_id}
                
            if individual_id not in individuals_dict:
                individuals_dict[individual_id] = {
                        'id': individual_id
                        }
            
            callsets_dict[callset_id] = {
                    'id': callset_id,
                    'biosample_id': biosample_id_tumor,
                    'variantset_id': 'GRCh38'}


            if variant_id not in variants_dict:
                
#                all_bases = {ref_allele:0}
#                if tumor_allele_1 not in all_bases:
#                    all_bases[tumor_allele_1] = len(all_bases)
#                if tumor_allele_2 not in all_bases:
#                    all_bases[tumor_allele_2] = len(all_bases)
#                
#                alt_bases = []
#                for i in range(len(all_bases)):
#                    alt_bases.append(all_bases.
#                
#                alt_bases = [tumor_allele_1,tumor_allele_2]
#                if ref_allele in alt_bases:
#                    alt_bases.remove(ref_allele)
#                all_bases = alt_bases.copy()
#                all_bases.insert(0,ref_allele)
                
                all_bases = []
                all_bases.append(ref_allele)
                if tumor_allele_1 not in all_bases:
                    all_bases.append(tumor_allele_1)
                if tumor_allele_2 not in all_bases:
                    all_bases.append(tumor_allele_2)
                
                variants_dict[variant_id] = {
                        'id': variant_id,
                        'reference_name': chro,
                        'start': start,
                        'end': end,
                        'variant_set': 'GRCh38',
                        'reference_bases': ref_allele,
                        'alternate_bases': all_bases[1:],
                        'variant_type': variant_type,
                        'calls': [{'callset_id': callset_id,
                                   'genotype': [all_bases.index(tumor_allele_1), all_bases.index(tumor_allele_2)],
                                   'info': {'biosample_id': biosample_id_tumor}
                                   }]
                    }
            else:
                call_ids = []
                for call in variants_dict[variant_id]['calls']:
                    call_ids.append(call['callset_id'])
                if callset_id not in call_ids:
                    all_bases = variants_dict[variant_id]['alternate_bases']
                    all_bases.insert(0, ref_allele)
                    if tumor_allele_1 not in all_bases:
                        all_bases.append(tumor_allele_1)
                    if tumor_allele_2 not in all_bases:
                        all_bases.append(tumor_allele_2)
                    variants_dict[variant_id]['alternate_bases'] = all_bases[1:]
                    variants_dict[variant_id]['calls'].append({
                            'callset_id': callset_id,
                            'genotype': [all_bases.index(tumor_allele_1), all_bases.index(tumor_allele_2)],
                            'info': {'biosample_id': biosample_id_tumor} 
                                })
                    
                
def parse_count(filepath):
    with gzip.open(filepath, 'rt') as fin:
        for i in range(6):
            next(fin)
        for line in fin:
            line = line.strip().split('\t')
#            if len(line[10]) > 1:
#                print(filepath)
#                print(line)
            if line[10] != line[11]:
                print(filepath)
                print(line)
         
def write2db():
    db_name = 'tcga_openCNV_ga4gh'
    
    db = MongoClient()[db_name]
    db['biosamples'].drop()
    db['biosamples'].insert_many(biosamples_dict.values())
    db['individuals'].drop()
    db['individuals'].insert_many(individuals_dict.values())
    db['callsets'].drop()
    db['callsets'].insert_many(callsets_dict.values())
    db['variants'].drop()
    db['variants'].insert_many(variants_dict.values())      
       


fpath = '/Volumes/originalData/TCGA/2018-01-08_TCGA_OpenMaskedSNV'
for root,subdirs,files in os.walk(fpath):
    for f in files:
        if f.lower().endswith('.gz'):
            parse_cnvgz(os.path.join(root,f))
#            parse_count(os.path.join(root,f))
    




#fpath = '/Volumes/originalData/TCGA/2018-01-08_TCGA_OpenMaskedSNV/03652df4-6090-4f5a-a2ff-ee28a37f9301/TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.gz'
#parse_cnvgz(fpath)

write2db()