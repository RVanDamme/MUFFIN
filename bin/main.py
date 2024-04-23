# main.py
from db_utils import *

# Chemin de base vers la base de données
DATABASE_PATH = "./pankegg2.db"

# Chemins des fichiers de données
PATHWAY_FILE = '../pathway.txt'
KO_FILE = '../ko.txt'
ANNOTATION_PATH = '../annotate/bin_annotation/*.annotations.tsv'
SOURMASH_TXT_DIR = '../classify/sourmash/*.txt'
CHECKM2_DIR = '../classify/checkm2/checkm2_dir/quality_report.tsv'

TABLES_TO_DROP = [
    'taxonomy', 'bin', 'map', 'kegg',
    'bin_map', 'map_kegg', 'bin_extra', 'bin_extra_kegg'
]


def main():
    conn = connect_db(DATABASE_PATH)
    cur = conn.cursor()

    drop_all_tables(cur, TABLES_TO_DROP)
    create_tables(cur)
    process_taxonomy_files(cur, SOURMASH_TXT_DIR)
    add_bin_quality_values(cur, CHECKM2_DIR)
    ko_translate_table = load_ko_descriptions(KO_FILE)
    map_translate_table = load_pathways(PATHWAY_FILE)
    map_ko_dict = process_annotation_file(cur, ANNOTATION_PATH, map_translate_table)
    set_ko(cur, ko_translate_table, map_ko_dict)
    link_maps_to_kos(cur, map_ko_dict)
    set_full_annot_table(cur, ANNOTATION_PATH)

    conn.commit()
    cur.close()
    conn.close()


if __name__ == "__main__":
    main()
