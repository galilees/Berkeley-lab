
import numpy as np
import pandas as pd
import json
import os

def open_file():
    path_to_json = '/net/cci-filer3/home/galilees/pdb_survey_gol/'
    data_list = []
    for file_name in [file for file in os.listdir(path_to_json) if file.endswith('.json')]:
        with open(path_to_json + file_name) as json_file:
            data = json.load(json_file)
            data_list.append(data)
    return data_list
   
def data_fail_count(data_list):
    count = 0 
    for data in data_list:
        success = data['success']
        if success == False:
            count += 1
    return count 
    

def data_check(data_list):
    checked_data = []
    for data in data_list:
    
        success = data['success']
        if success == True:
            checked_data.append(data)

    return checked_data

def get_data(checked_data, resname):

     for nested_dict in data_list:
       data = nested_dict[resname]
     #print(nested_dict.items())
     return data 

def get_nearby_res(data):
    nearby_list =[]
    for k,v  in data.items():
        if k == "near_by_res":
            nearby_list.append(v)
    print(nearby_list)
    return nearby_list

def get_res_count(nearby_list):
    res_count_list = []
    for i in nearby_list:
        for j in i:
            for dict in j:
              print("inner structure", dict,type(dict))
            #   for k,v in dict:
            #       res_count_list.append(v)
    return res_count_list

def data_frame(res_count_list):
    df = pd.DataFrame(res_count_list)
    return df

data_list = open_file()
#print("all data from json files: ", data_list)
print("Total number of files read:", len(data_list))

data_fail = data_fail_count(data_list)
print("The number of files that failed check: ", data_fail)

list_checked_data = data_check(data_list)
print("number of files that passed the fail check: ",  len(list_checked_data))
#print("Data from files that that passed the fail check: ", list_checked_data )

gol_data = get_data(list_checked_data,'GOL')
#print("GOL data from files" , gol_data )
nearby_res_gol= get_nearby_res(gol_data)
#print("nearby residues:", nearby_res_gol)

res_count_gol = get_res_count(nearby_res_gol)
# # print("count of residues near gol:", res_count_gol)

# hoh_data = get_data(data_list, 'HOH')
# # print("GOL data from files" , gol_data )
# nearby_res_hoh= get_nearby_res(hoh_data)
# # print("nearby residues:", nearby_res_gol)
# res_count_hoh = get_res_count(nearby_res_hoh)
# # print("count of residues near gol:", res_count_gol)

# df_gol = data_frame(res_count_gol)
# print("original dataframe", df_gol)
# df_hoh = data_frame(res_count_hoh)
# print(df_hoh)


# df_gol_0 = df_gol.fillna(0)
# df_hoh_0 = df_hoh.fillna(0)
# print(df_gol_0)
# print(df_hoh_0)
# print("oiginal dataframe columns", df_gol_0.columns)
# print("oiginal dataframe columns", df_hoh_0.columns)

# # final_table_columns =['HOH','ARG','HIS','LYS','ASP','GLU','SER','THR','ASN','GLN','CYS','GLY','PRO','ALA','VAL','ILE','LEU','MET','PHE','TYR','TRP']
# # df_gol_sample_aa = df_gol_0[df_gol_0.columns.intersection(final_table_columns)]

# # df_hoh_sample_aa = df_hoh_0[df_hoh_0.columns.intersection(final_table_columns)]
# # print(df_gol_sample_aa)
# # print(df_gol_sample_aa.columns)


