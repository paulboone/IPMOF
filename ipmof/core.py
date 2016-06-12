import os

import xlrd


def mof_properties(core_file_path):
    """
    Reads properties of MOFs in the CoRE database from 'CoRE.xlsx' file in 'doc' folder.
    Properties are exported in dictonary format with following inputs:
    - mof_names - metal_type - density - pld - lcd - vsa - gsa
    - void_fraction - pore_volume - unmofied_mofs
    """
    # Read Excel File
    core_data = xlrd.open_workbook(excel_file_path)
    # Read columns to acquire force field parameters
    mof_names = core_data.sheets()[0].col_values(1)[9:]
    metal_type = core_data.sheets()[0].col_values(2)[9:]
    density = core_data.sheets()[0].col_values(3)[9:]
    pld = core_data.sheets()[0].col_values(4)[9:]
    lcd = core_data.sheets()[0].col_values(5)[9:]
    vsa = core_data.sheets()[0].col_values(6)[9:]
    gsa = core_data.sheets()[0].col_values(7)[9:]
    void_fraction = core_data.sheets()[0].col_values(8)[9:]
    pore_volume = core_data.sheets()[0].col_values(9)[9:]
    unmofied_mofs = core_data.sheets()[9].col_values(1)[3:]

    mof_properties = {'mof_names': mof_names,
                      'metal_type': metal_type,
                      'density': density,
                      'pld': pld,
                      'lcd': lcd,
                      'vsa': vsa,
                      'gsa': gsa,
                      'void_fraction': void_fraction,
                      'pore_volume': pore_volume,
                      'unmofied_mofs': unmofied_mofs
                      }

    return mof_properties


def high_void_fraction(mof_properties, vf_limit):
    """
    Selects MOFs with void fraction higher than given void fraction limit.
    """
    void_fraction = mof_properties['void_fraction']
    high_vf_mofs = {'vf': [], 'name': []}

    # vf_limit = 0.9
    for vf_index, vf in enumerate(void_fraction):
        if vf > vf_limit:
            high_vf_mofs['vf'].append(vf)
            high_vf_mofs['name'].append(mof_names[vf_index])

    print('Gathered a total of', len(high_vf_mofs['name']), 'MOFs')
    print('With void fraction >', vf_limit)
