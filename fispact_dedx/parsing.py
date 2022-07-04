import os 
import pandas as pd 
import numpy as np

def get_flux(root): 
    fluxes = {}
    #for item in os.listdir(root):
    os.chdir(root)
    directories=[d for d in os.listdir(root) if os.path.isdir(d)]
    #print(directories)
    for item in directories:
        #print(item)
        for file in os.listdir(f"{root}/{item}"):
            if file == 'fluxes':
                try:
                    fluxes[item] = np.genfromtxt(f"{root}/{item}/fluxes", skip_footer=1).flatten()
                except: 
                    fluxes[item] = np.append(np.genfromtxt(f"{root}/{item}/fluxes", skip_footer=2).flatten(),0)
    return fluxes
                
def extract_nuclide(filepath):
    data = {}
    irrad_list = []
    cool_list = []
    with open(filepath, "r") as file: 
        counter = 0
        time_interval = 0
        for line in file:
            if "TIME INTERVAL" in line: 
                counter = 1
                time_interval = line[23:25].strip()
                flux_val = float(line.split("FLUX AMP IS")[1].split("/cm^2/s")[0].strip())
                time_val = line.split("ELAPSED TIME IS")[1].split('*')[0]
                if flux_val == 0: 
                    cool_list.append(time_val)
                else: 
                    irrad_list.append(time_val)
                val = []
            elif counter >= 1 and line[0] == '0': 
                #print(val)
                try: 
                    df = pd.DataFrame(val, columns = ["ELEMENT", "MASS", "STABILITY", "ATOMS","GRAMS","Bq","b-Energy","a-Energy","g-Energy","DOSE RATE"])
                    #df["ATOMS","GRAMS","Bq","b-Energy","a-Energy","g-Energy","DOSE RATE"] = pd.to_numeric(df["ATOMS","GRAMS","Bq","b-Energy","a-Energy","g-Energy","DOSE RATE"])
                except: 
                    df = pd.DataFrame(val, columns = ["ELEMENT", "MASS", "STABILITY", "ATOMS","GRAMS","Bq","b-Energy","a-Energy","g-Energy","DOSE RATE", "HALF-LIFE"])
                    #df["ATOMS","GRAMS","Bq","b-Energy","a-Energy","g-Energy","DOSE RATE", "HALF-LIFE"] = pd.to_numeric(df["ATOMS","GRAMS","Bq","b-Energy","a-Energy","g-Energy","DOSE RATE", "HALF-LIFE"])
                data[time_interval] = df
                counter = 0
            elif counter >= 1:
                counter += 1
            if counter > 4 and line != '': 
                line = [line[2:4]] + [line[4:7]] + [line[11:13]] + [float(x) if x != 'Stable' else x for x in line[13:].split()]
                val.append(line)
                continue
    return data, irrad_list, cool_list


def convert_time(irrad_list, cool_list): 
    def inner_convert(time_list):
        for idx in range(len(time_list)): 
            item = time_list[idx].strip()
            if 'ps' in item: 
                num = item.strip('ps')
                mul = 1e9
            elif 's' in item: 
                num = item.strip('s')
                mul = 1/(3600*24) 
            elif 'm' in item: 
                num = item.strip('m')
                mul = 1/(60*24)
            elif 'h' in item: 
                num = item.strip('h')
                mul = 1/24
            elif 'd' in item: 
                num = item.strip('d')
                mul = 1 
            elif 'y' in item: 
                num = item.strip('y')
                mul = 365.25
            time_list[idx] = float(num) * mul
        return time_list
    irrad_list = inner_convert(irrad_list)
    cool_list = inner_convert(cool_list)
    end_irrad = irrad_list[-1]
    irrad_list = [x-end_irrad for x in irrad_list]
    full_list = irrad_list + cool_list
    return full_list, irrad_list[-1]

def extract_data(root, data_type, top_energy): 
    full_data = {}
    for item in os.listdir(root):
        if 'DOSE' not in item and 'Bq' not in item and 'ATOMS' not in item:
            val = float(item)
            if top_energy >= val:
                for file in os.listdir(f"{root}/{item}"):
                    if file == 'inventory.out':
                        if data_type == "ATOMS" or data_type == "Bq" or data_type == "DOSE RATE":
                            full_data[item], irrad_list, cool_list = extract_nuclide(f"{root}/{item}/inventory.out")
                            time_list, irrad_finish = convert_time(irrad_list, cool_list)
                        else: 
                            raise ValueError(f"Data type must be 'ATOMS', 'Bq' or 'DOSE' - not {data_type}")
    return full_data, time_list

def add_total_files(dataframe):
    element_list = []
    for product in dataframe['PRODUCT']:
        if 'prod' not in product:
            element = product.split('-')[0]
            if element not in element_list:
                element_list.append(element)
    dataframe = dataframe.set_index('PRODUCT')
    for element in element_list:
        dataframe.loc[f"{element}-total", :] = dataframe[
            (dataframe.index.str.contains(element)) & (
                    dataframe.index.str.contains('tot') == False)].sum(axis=0,
                                                                                             skipna=True)
    dataframe.loc['total', :] = dataframe[
        (dataframe.index.str.contains("tot") == False) ].sum(axis=0,
                                                                                         skipna=True)
    dataframe = dataframe.reset_index()
    #dataframe.loc['tot' in dataframe['PRODUCT'], 'STABILITY'] = 10
    for i, row in dataframe.iterrows(): 
        if 'tot' in row['PRODUCT']: 
            dataframe.loc[i, 'STABILITY'] = ''
    return dataframe

def concat_extracted_data(full_data, time_list, data_type): 
    concat_dic = {}
    time_int_dic = full_data[list(full_data.keys())[0]].keys()
    for time in time_int_dic:
        new_data = []
        for run in full_data:
            new_data.append(full_data[run][time])
        try: 
            df = pd.concat(new_data).groupby(by=["ELEMENT", "MASS", "STABILITY", "HALF-LIFE"]).sum().reset_index()
        except: 
            df = pd.concat(new_data).groupby(by=["ELEMENT", "MASS", "STABILITY"]).sum().reset_index()
        concat_dic[time] = df

    cols = ["ELEMENT", "MASS", "STABILITY", f"{data_type}"]
    df_list = []
    for time_step in concat_dic.keys(): 
        concat_dic[time_step] = concat_dic[time_step][cols]
        concat_dic[time_step].set_index(['ELEMENT', 'MASS', 'STABILITY'], inplace=True)
        concat_dic[time_step].columns = [time_step]
        df_list.append(concat_dic[time_step])
    df3 = pd.concat(df_list)#, join="outer", sort=True, axis=0, ignore_index=True)
    df3 = df3.reset_index() 
    df3 = df3.groupby(by=["ELEMENT", "MASS", "STABILITY"]).sum().reset_index()
    df3.insert(0, 'PRODUCT', df3['ELEMENT'].str.strip()+'-'+df3['MASS'].str.strip())
    del df3['ELEMENT'] 
    del df3['MASS']
    #df3.set_index(['PRODUCT'], inplace=True)
    df3 = add_total_files(df3)
    df3.columns = ["PRODUCT", "STABILITY"] + time_list
    #del df3['STABILITY']
    return df3