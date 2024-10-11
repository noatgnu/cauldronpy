import os
import numpy as np
import pandas as pd
import click


def create_contrast_matrix(annotation_file: pd.DataFrame) -> pd.DataFrame:
    # Load the annotation file


    # Get unique conditions
    unique_conditions = annotation_file['Condition'].unique()

    # Create an empty DataFrame for the contrast matrix
    contrast_matrix = pd.DataFrame(0, index=[], columns=unique_conditions)

    # Create comparisons
    for i, condition1 in enumerate(unique_conditions):
        for condition2 in unique_conditions[i + 1:]:
            comparison_name = f"{condition1}_vs_{condition2}"
            contrast_matrix.loc[comparison_name] = 0
            contrast_matrix.loc[comparison_name, condition1] = 1
            contrast_matrix.loc[comparison_name, condition2] = -1

    return contrast_matrix

class MSstats:
    def __init__(self, r_home: str):
        self.setup_r_home(r_home)
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.packages import importr
        self.msstats = importr("MSstats")
        pandas2ri.activate()
        self.pandas2ri = pandas2ri
        self.ro = ro


    def setup_r_home(self, r_home: str):
        if not os.path.exists(r_home):
            raise FileNotFoundError(f"R home directory {r_home} does not exist")
        os.environ["R_HOME"] = r_home

    def load_data(self, file_path: str) -> pd.DataFrame:
        if file_path.endswith(".csv"):
            return pd.read_csv(file_path)
        elif file_path.endswith(".tsv") or file_path.endswith(".txt"):
            return pd.read_csv(file_path, sep="\t")
        else:
            raise ValueError("Unsupported file format")

    def run_msstats(self, annotation_file: str, report_file: str, matrix_file: str = None, input_type: str = "DIA-NN"):
        report_file = pd.read_csv(report_file, sep="\t")
        annotation_file = pd.read_csv(annotation_file, sep="\t")

        # Create contrast matrix if not provided
        if matrix_file is None:
            contrast_matrix = create_contrast_matrix(annotation_file)
        else:
            contrast_matrix = pd.read_csv(matrix_file, sep="\t")

        if input_type == "DIA-NN":
            msstats_data = self.msstats.DIANNtoMSstatsFormat(
                report_file,
                annotation=annotation_file,
            )
        elif input_type == "Skyline":
            msstats_data = self.msstats.SkylinetoMSstatsFormat(
                report_file,
                annotation=annotation_file,
            )
        elif input_type == "Spectronaut":
            msstats_data = self.msstats.SpectronautoMSstatsFormat(
                report_file,
                annotation=annotation_file,
            )
        else:
            raise ValueError("Unsupported input type")

        with (self.ro.default_converter + self.pandas2ri.converter).context():
            msstats_data_df = self.ro.conversion.rpy2py(msstats_data)

        print("msstats_data dimensions:", msstats_data_df.shape)
        print("msstats_data sample:", msstats_data_df.head())

        if msstats_data_df.shape[0] == 0 or msstats_data_df.shape[1] == 0:
            raise ValueError("msstats_data is empty or has invalid dimensions")

        processed_data = self.msstats.dataProcess(msstats_data_df, use_log_file=True)
        print("Processed data dimensions:", processed_data.shape)

        if processed_data.shape[0] == 0 or processed_data.shape[1] == 0:
            raise ValueError("Processed data is empty or has invalid dimensions")
        result, *_ = self.msstats.groupComparison(constrast_matrix=contrast_matrix, data=processed_data, use_log_file=True)
        return result

@click.command()
@click.option("--annotation", help="Annotation file path", required=True)
@click.option("--report", help="Report file path", required=True)
@click.option("--matrix", help="Matrix file path", required=True)
@click.option("--input_type", help="Input type", required=True)
@click.option("--r_home", help="R home directory", required=True)
def main(annotation, report, matrix, input_type, r_home):
    msstats = MSstats(r_home)
    result = msstats.run_msstats(annotation, report, matrix, input_type)
    print(result)



