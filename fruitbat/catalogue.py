import os
import pandas as pd

from e13tools import docstring_substitute

from fruitbat import Frb, methods, cosmologies

__all__ = ["create_analysis_catalogue", "read_frb_row"]


@docstring_substitute(meth=methods.avaliable_methods(),
                      cosmo=cosmologies.avaliable_cosmologies())
def create_analysis_catalogue(filename="fruitbat_analysis_catalogue",
                              dataset='default', method='Inoue2004',
                              cosmology='Planck18'):
    """
    Analyses an FRB dataset and produces a CSV file containing the
    estimated redshift, fluence, energy and luminosity for each FRB in
    additional to its measured quantities.

    Parameters
    ----------
    filename: str, optional
        The output file name. Default: 'fruitbat_analysis_catalogue'

    dataset: str, optional
        The path to the FRBCAT dataset. The dataset is required to have
        the following columns: 'frb_name', 'utc', 'telescope',
        'rop_raj', 'rop_decj', 'rop_gl', 'rop_gb', 'rop_bandwidth',
        'rop_centre_frequency', 'rmp_dm', 'rmp_width', 'rmp_snr',
        'rmp_flux'. If ``dataset='default'`` then the builtin dataset
        will be used. The builtin dataset was last updated: 2019-04-08.
        Default: 'default'

    method: str, optional
        The dispersion  measure - redshift relation to use when
        calculating the redshift. Avaliable methods: %(meth)s.
        Default: 'Inoue2004'

    cosmology: str, optional
        The cosmology to assume when calculating redshift.
        Avaliable cosmologies: %(cosmo)s. Default: 'Planck18'

    Generates
    ---------
    A csv file with the output of the of the analysis.
    """
    if dataset == 'default':
        dataset = os.path.join('data', 'frbcat_database_20190408.csv')

    df = pd.read_csv(dataset)

    columns = ["Name", "Telescope", "RAJ", "DECJ", "gl", "gb",
               "DM (pc/cm3)", "Width (ms)", "Bandwidth (MHz)",
               "Centre_Frequency (MHz)", "Flux (Jy)", "SNR",
               "DM_galaxy (pc/cm3)", "z",  "Fluence (Jy ms)", "Energy (erg)",
               "Luminosity (erg/s)", "Method", "Cosmology"]

    df_out = pd.DataFrame(index=range(0, len(df)), columns=columns)
    for item, row in df.iterrows():
        data = read_frb_row(row)

        frb = Frb(data["dm"], raj=data["raj"], decj=data["decj"],
                  gl=data["gl"], gb=data["gb"], width=data["width"],
                  peak_flux=data["flux"], obs_bandwidth=data["bw"],
                  obs_freq_central=data["centre_freq"], name=data["name"])

        # Calculate FRB properties
        frb.calc_dm_galaxy()
        frb.calc_redshift()
        energy = frb.calc_energy()
        luminosity = frb.calc_luminosity()

        df_out.iloc[item] = [frb.name, data["telescope"], frb.raj,
                             frb.decj, frb.gl, frb.gb, frb.dm.value,
                             frb.width.value, frb.obs_bandwidth.value,
                             frb.obs_freq_central.value, frb.peak_flux.value,
                             frb.snr, frb.dm_galaxy.value, frb.z.value,
                             frb.fluence.value, energy.value,
                             luminosity.value, method, cosmology]

        output_name = "".join([filename, ".csv"])
        df_out.to_csv(output_name)


def read_frb_row(row):
    """
    Reads the row of a :obj:`~pandas.DataFrame` and retrieves the 
    data in the correct format.

    Parameters
    ----------
    row: :obj:`~pandas.core.series.Series`
        The series containing FRB data.

    Returns
    -------
    A dictionary containing the FRB paramemters
    """
    drow = {"name": row['frb_name'],
            "utc": row['utc'],
            "telescope": row['telescope'],
            "dm": float(row['rmp_dm'].split('&plusmn')[0]),
            "gl": row['rop_gl'],
            "gb": row['rop_gb'],
            "raj": row['rop_raj'],
            "decj": row['rop_decj'],
            "bw": float(row['rop_bandwidth']),
            "width": float(row['rmp_width']),
            "snr": float(row['rmp_snr']),
            "flux": float(row['rmp_flux']),
            "centre_freq": float(row['rop_centre_frequency'])
            }
    return drow
