from pymongo import MongoClient
import json
import os
import datetime
import sys


#######################################################################
## This script convert TCGA masked copy number variation files into ###
## GA4GH data schema, and store it into a database                  ###
#######################################################################


#######################################################################
# Requirements:
# 1. a directiory which has the TCGA data in the orginal folder & file  
#   structure.
# 2. a JSON file downloaded through the TCGA API, which contains the 
#   meta-info of each file
#######################################################################


# Check arguments
if len(sys.argv) != 3:
    print("Usage: >python tcga_cnv_2_ga4gh.py [path of data directiory] [path of meta-file]")
    sys.exit()
else:
    datapath = sys.argv[1]
    metafile = sys.argv[2]



# Init variables

# datapath = '/Volumes/arraymapIncoming/TCGA/MaskedCNV_2017_10_03'
# datapath = '/Volumes/originalData/TCGA/MaskedCNV_2017_10_03'
# metafile = '/Users/bogao/Desktop/projects/TCGA/tcga_masked_cnv_meta_171205.json'
db_name = 'tcga_cnv_masked_ga4gh'
with open(metafile, 'r') as fi:
    data = json.load(fi)

# hits = data['data']['hits'][0:10]
hits = data['data']['hits']

biosamples = {}
individuals = {}
callsets = {}
variants = {}


# walk through each file in the meta-file
for file in hits:
    file_path = '{}/{}/{}'.format(datapath,file['id'], file['file_name'])
    
    # init ids
    biosample_id = file['cases'][0]['samples'][0]['sample_id']
    individual_id = file['cases'][0]['case_id']
    project_id = file['cases'][0]['project']['project_id']
    callset_id = file['id']


    # biosample
    biosamples[biosample_id] = {
        'id': biosample_id,
        'name': project_id,
        'description': file['cases'][0]['samples'][0]['sample_type'],
        'bio_characteristics': [
            {
                'description': file['cases'][0]['samples'][0]['sample_type'],
                'ontology_terms': [
                    {
                        'term_id': 'idcom:_{}'.format(file['cases'][0]['diagnoses'][0]['morphology']),
                        'term_label': None
                    },
                    {
                        'term_id': 'icdot:_{}'.format(file['cases'][0]['diagnoses'][0]['tissue_or_organ_of_origin']),
                        'term_label': None
                    }
                ],
                'negated_ontology_terms': []
            }
        ],
        'created': datetime.datetime.utcnow(),
        'individual_id': individual_id,
        'individual_age_at_collection': None,
        'external_identifiers': [
            {
                'database': 'TCGA',
                'identifier': biosample_id
            }
        ],
        'location': None,
        'attributes':{
            'tumor_stage': file['cases'][0]['diagnoses'][0]['tumor_stage']
        }
    }
    
    
    # individual
    if individual_id not in individuals:
        sex = {'term_id': 'PATO:0020000', 'term_label': 'genotypic sex' }
        if file['cases'][0]['demographic']['gender'] == 'male':
            sex = {'term_id': 'PATO:0020001', 'term_label': 'male genotypic sex' }
        elif file['cases'][0]['demographic']['gender'] == 'female':
            sex = {'term_id': 'PATO:0020002', 'term_label': 'female genotypic sex' }
        
        individuals[individual_id] = {
            'id': individual_id,
            'species': {'term_id': 'NCBITaxon:9606', 'term_label': 'Homo sapiens' },
            'sex': sex,
            'external_identifiers': [
                {
                    'database': 'TCGA',
                    'identifier': individual_id
                }
            ],
            'created': datetime.datetime.utcnow(),
            'attributes':{
                'death': file['cases'][0]['diagnoses'][0]['vital_status'],
                'year_of_birth':file['cases'][0]['demographic']['year_of_birth'],
                'age_at_diagnosis':file['cases'][0]['diagnoses'][0]['age_at_diagnosis'],
                'days_to_death':file['cases'][0]['diagnoses'][0]['days_to_death'],
                'race':file['cases'][0]['demographic']['race'],
                'ethnicity':file['cases'][0]['demographic']['ethnicity']
            }
        }
    

    # callset
    if callset_id not in callsets:
        callsets[callset_id] = {
            'id': callset_id,
            'biosample_id': biosample_id,
            'variantset_id': 'GRCH38',
            'created': datetime.datetime.utcnow(),
            'info': {
                'data_type':file['data_type'],
                'file_name':file['file_name']
            }
        }
    

    # variants
    if os.path.isfile(file_path):
        with open(file_path, 'r') as fin:
            next(fin)
            variants_cnv = []
            for line in fin:
                line = line.split()
                chro = str(line[1])
                start = int(float(line[2]))
                end = int(float(line[3]))
                probes = int(line[4])
                value = float(line[5])
                
                if value >0:
                    variant_type = 'DUP'
                    alternate_bases = '<DUP>'
                    cipos = [-1000,1000]
                    ciend = [-1000,1000]
                else:
                    variant_type = 'DEL'
                    alternate_bases = '<DEL>'
                    cipos = [-1000,1000]
                    ciend = [-1000,1000]
                
                tag = '{}_{}_{}_{}'.format(chro, start, end, variant_type)
                
                call = {
                    'call_set_id': callset_id, 
                    'genotype': ['.', '.'], 
                    'info': {
                        'segvalue': value, 
                        'biosample_id': biosample_id,
                        'probes':probes
                    }
                }
                
                if tag in variants:
                    variants[tag]['updated'] = datetime.datetime.utcnow()
                    variants[tag]['calls'].append(call)
                else:
                    variants[tag] = {
                        'id': tag,
                        'reference_name': chro,
                        'start':start,
                        'end':end,
                        'variant_set': 'GRCH38',
                        'updated': datetime.datetime.utcnow(),
                        'created': datetime.datetime.utcnow(),
                        'reference_bases': '.',
                        'alternate_bases': alternate_bases,
                        'variant_type': variant_type,
                        'cipos': cipos,
                        'ciend': ciend,
                        'calls': [call]
                    }
                
    else:
        print(file['file_name'])

# write to database
db = MongoClient()[db_name]
db['biosamples'].drop()
db['biosamples'].insert_many(biosamples.values())
db['individuals'].drop()
db['individuals'].insert_many(individuals.values())
db['callsets'].drop()
db['callsets'].insert_many(callsets.values())
db['variants'].drop()
db['variants'].insert_many(variants.values())