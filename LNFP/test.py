import os
import cobra
import pandas as pd
def add_sink(complementary_sink, folder_to_save):
    model = cobra.io.read_sbml_model(complementary_sink)
    metabolites_data = []
    # 遍历模型中的所有代谢物
    df = pd.DataFrame(columns=['name', 'inchi'])
    for metabolite in model.metabolites:
        # 获取代谢物名称
        name = metabolite.name

        # 获取代谢物的InChIKey
        annotation_keys = metabolite.annotation.keys()
        if 'inchi_key' in annotation_keys:
            inchikey = metabolite.annotation['inchi_key'] if 'inchi_key' in metabolite.annotation else ''
            print(inchikey)
            # inchi = inchikey2inchi(inchikey)
            # df.loc[len(df.index)] = [name, inchi]

    # sink_file_path = os.path.join(folder_to_save, 'add_sink.csv')
    # df.to_csv(sink_file_path, index=False)
    # organism = get_organism_from_dataframe(df, add_Hs=True)

    # return organism
add_sink('lnfp.xml', '.')
