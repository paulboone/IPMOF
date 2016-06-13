import os

import xlrd


def core_mof_properties(core_file_path):
    """
    Reads properties of MOFs in the CoRE database from 'CoRE.xlsx' file in 'doc' folder.
    Properties are exported in dictonary format with following inputs:
    - mof_names - metal_type - unmofied_mofs
    - density (g/cm3) - pore_volume (cm3/g) - void_fraction
    - pld (pore limiting diamater)
    - lcd (largest cavity diamater)
    - vsa (volumetric surface area m2/cm3)
    - gsa (gravimetric surface area m2/g)
    CoRE Database: 10.1021/cm502594j
    """
    # Read Excel File
    core_data = xlrd.open_workbook(core_file_path)
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

    mof_properties = {'mof_name': mof_names,
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


def core_mof_sort(mof_properties, sort='void_fraction', limit='0.9'):
    """
    Selects MOFs with given property higher than given limit.
    Default property is void fraction with a limit of 0.9.
    Available properties are:
        - density (g/cm3) - pore_volume (cm3/g) - void_fraction
        - pld (pore limiting diamater) - lcd (largest cavity diamater)
        - vsa (volumetric surface area m2/cm3)
        - gsa (gravimetric surface area m2/g)
    """
    sorted_mofs = {sort: [], 'name': []}
    sort_prop = mof_properties[sort]

    for prop_index, prop in enumerate(sort_prop):
        if prop > limit:
            sorted_mofs[sort].append(prop)
            sorted_mofs['name'].append(mof_properties['mof_name'][prop_index])

    print('Gathered a total of', len(sorted_mofs['name']), 'MOFs')
    print('With', sort, '>', limit)

    return sorted_mofs


def core_mof_dir(sorted_mofs, mol2_dir):
    """
    Collects directories for the given MOF list (sorted_mofs).
    The directories are taken from the database directory (mol2_dir) provided.
    Collected directories are returned in list format.
    """
    all_mofs = os.listdir(mol2_dir)
    mof_dirs = []
    missing_mofs = []

    for mof in sorted_mofs['name']:
        mof_found = False
        for mof_file_name in all_mofs:
            mof_name = mof_file_name.split('.')[0].split('_')[0]
            if mof == mof_name:
                mof_dir = os.path.join(mol2_dir, mof_file_name)
                mof_dirs.append(mof_dir)
                mof_found = True
        if not mof_found:
            missing_mofs.append(mof)

    print(len(missing_mofs), 'mofs are missing')
    print(len(mof_dirs), 'total mofs found')

    return mof_dirs
