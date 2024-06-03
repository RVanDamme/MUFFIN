# db_utils.py

import sqlite3
import csv
import glob
import os
from sql_commands import *


def connect_db(path):
    return sqlite3.connect(path)


def create_tables(cur):
    for table, sql in CREATE_TABLES.items():
        cur.execute(sql)
    # conn.commit()


def drop_all_tables(cur, tables):
    for table in tables:
        cur.execute(DROP_TABLE.format(table_name=table))
    # conn.commit()


def drop_table(conn, table_name):
    cur = conn.cursor()
    cur.execute(DROP_TABLE.format(table_name=table_name))
    conn.commit()


def insert_taxonomy(cur, tax_data):
    # Vérifier si l'entrée existe déjà
    cur.execute(SELECT_TAXONOMY_ID, tax_data[2:9])
    result = cur.fetchone()
    if result:
        return result[0]
    else:
        # Insérer la nouvelle entrée taxonomique
        cur.execute(INSERT_TAXONOMY, tax_data[2:9])
        return cur.lastrowid


def insert_bin(cur, bin_name, taxonomic_id):
    cur.execute(INSERT_BIN, (bin_name, taxonomic_id))


def load_pathways(pathway_file):
    """Charge les données des pathways depuis un fichier et les stocke dans un dictionnaire."""
    pathways_dict = {}
    with open(pathway_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            map_id = row[0].strip()
            pathway_description = row[1].strip() if len(row) > 1 else None
            pathways_dict[map_id] = pathway_description
    return pathways_dict


def load_ko_descriptions(ko_file):
    """Charge les descriptions des KO depuis un fichier et les stocke dans un dictionnaire."""
    ko_dict = {}
    with open(ko_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            ko_id = row[0].strip()
            description = row[1].split(";") if len(row) > 1 else ["Unknown"]
            ko_dict[ko_id] = description
    return ko_dict


def process_pathways_and_kos(cur, pathway, maps_dict, pathways_dict, kos, bin_map_list):
    if pathway not in bin_map_list:
        bin_map_list.append(pathway)

    if pathway not in maps_dict:
        # Tente d'insérer le pathway dans la base de données si ce n'est pas déjà fait
        description = pathways_dict.get(pathway, "Description not available")  # Utilise une description par défaut si non trouvée
        cur.execute(INSERT_MAP, (description, pathway, pathway))
        maps_dict[pathway] = kos
    else:
        map_ko_list = maps_dict.get(pathway)
        original_length = len(map_ko_list)
        map_ko_list.extend([ko for ko in kos if ko not in map_ko_list])
        if len(map_ko_list) > original_length:
            maps_dict[pathway] = map_ko_list


def link_bins_to_pathways(cur, bin_name, bin_map_list):
    for map_number in bin_map_list:
        # Récupérer l'ID du bin
        cur.execute(SELECT_BIN_ID, (bin_name,))
        bin_id = cur.fetchone()
        if bin_id:
            bin_id = bin_id[0]
            # Récupérer l'ID de la map
            cur.execute(SELECT_MAP_ID, (map_number,))
            map_id = cur.fetchone()
            if map_id:
                map_id = map_id[0]
                # Insérer la jointure dans bin_map si elle n'existe pas déjà
                cur.execute(INSERT_BIN_MAP, (bin_id, map_id, bin_id, map_id))


def link_maps_to_kos(cur, maps_dict):
    for map_name in maps_dict:
        cur.execute(SELECT_MAP_ID, (map_name,))
        map_id = cur.fetchone()
        if map_id:
            map_id = map_id[0]
            for ko in maps_dict.get(map_name):
                if not ko == "":
                    ko_entry = ko.split(":")[1]
                    cur.execute(SELECT_KEGG_ID, (ko_entry,))
                    ko_id = cur.fetchone()
                    if ko_id:
                        ko_id = ko_id[0]
                        # Insérer la jointure dans bin_map si elle n'existe pas déjà
                        cur.execute(INSERT_MAP_KEGG, (map_id, ko_id, map_id, ko_id))


def add_bin_quality_values(cur, quality_report_path):
    with open(quality_report_path, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            # Get the data for each bin
            bin_name = row['Name'] + ".fa"
            completeness = float(row['Completeness'])
            contamination = float(row['Contamination'])
            # Update the Bin table with the new values
            cur.execute(UPDATE_BIN_QUALITY, (completeness, contamination, bin_name))


def set_ko(cur, kegg_translate_table, maps_dict):
    unique_kos = set()
    for kos in maps_dict.values():
        unique_kos.update(kos)
    unique_ko_list = list(unique_kos)
    for ko in unique_ko_list:
        if not ko == "":
            ko_id = ko.split(":")[1]
            kegg_info = kegg_translate_table.get(ko_id)
            if kegg_info:
                kegg_name = kegg_info[0]
                kegg_full_name = kegg_info[1]
                cur.execute(INSERT_KEGG, (ko_id, kegg_name, kegg_full_name, ko_id))

            else:
                cur.execute(INSERT_KEGG, (ko_id, "null", "null", ko_id))


def process_taxonomy_files(cur, file_path):
    files = glob.glob(file_path)
    for file_path in files:
        with open(file_path, 'r') as file:
            next(file)  # Sauter la première ligne (entête)
            line = next(file).strip()  # Lire la deuxième ligne
            data = line.split(',')
            assigned_satus = data[1]
            bin_name = data[0]
            if not assigned_satus == "nomatch":
                tax_id = insert_taxonomy(cur, data)
                insert_bin(cur, bin_name, tax_id)
            else:
                insert_bin(cur, bin_name, "null")


def process_annotation_file(cur, annotation_files_path, pathways_dict):
    maps_dict = {}
    for filename in glob.glob(annotation_files_path):
        with open(filename, 'r', newline='') as file:
            base_name = os.path.basename(filename)
            name_part = base_name.split('.')[0] + '.' + base_name.split('.')[1]
            print("Fichier traité : ", name_part)

            reader = csv.DictReader(file, delimiter='\t')
            bin_map_list = []
            for row in reader:
                # Extraire la colonne "KEGG_Pathway"
                kegg_pathways = row['KEGG_Pathway']
                kegg_kos = row['KEGG_ko']
                kos = kegg_kos.split(',')
                if kegg_pathways:
                    for pathway in kegg_pathways.split(','):
                        pathway = pathway.strip()
                        if pathway.startswith('map'):
                            process_pathways_and_kos(cur, pathway, maps_dict, pathways_dict, kos, bin_map_list)

                else:
                    pathway = name_part + "_not_mapped"
                    process_pathways_and_kos(cur, pathway, maps_dict, pathways_dict, kos, bin_map_list)

            bin_name = name_part + ".fa"
            link_bins_to_pathways(cur, bin_name, bin_map_list)

    return maps_dict


def link_full_line_with_kos(cur, bin_id, kegg_gos, kegg_kos, kegg_free_desc):
    cur.execute(INSERT_BIN_EXTRA, (bin_id, kegg_gos, kegg_kos, kegg_free_desc))
    kegg_extra_line_id = cur.lastrowid
    if kegg_extra_line_id:
        kos = kegg_kos.split(',')
        for ko in kos:
            ko_entry = ko.split(":")[1]
            cur.execute(SELECT_KEGG_ID, (ko_entry,))
            ko_id = cur.fetchone()
            if ko_id:
                ko_id = ko_id[0]
                cur.execute(INSERT_BIN_EXTRA_KEGG, (kegg_extra_line_id, ko_id, kegg_extra_line_id, ko_id))


def set_full_annot_table(cur, annotation_files_path):
    for filename in glob.glob(annotation_files_path):
        with open(filename, 'r', newline='') as file:
            base = os.path.basename(filename)
            name_part = base.split('.')[0] + '.' + base.split('.')[1]
            bin_name = name_part + ".fa"
            cur.execute(SELECT_BIN_ID, (bin_name,))
            bin_id = cur.fetchone()
            if bin_id:
                bin_id = bin_id[0]
                reader = csv.DictReader(file, delimiter='\t')
                for row in reader:
                    # print(row)
                    kegg_kos = row['KEGG_ko']
                    kegg_gos = row['GOs']
                    kegg_free_desc = row['eggNOG free text desc.']

                    if kegg_kos:
                        link_full_line_with_kos(cur, bin_id, kegg_gos, kegg_kos, kegg_free_desc)
