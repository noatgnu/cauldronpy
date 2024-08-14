import os

import pandas as pd
import click

class MSstats:
    def __init__(self, r_home: str):
        self.setup_r_home(r_home)
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.packages import importr
        self.msstats = importr("MSstats")
        pandas2ri.activate()


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

    def run_msstats(self, annotation_file: str, report_file: str, matrix_file: str, input_type: str):
        report_file = pd.read_csv(report_file, sep="\t")
        annotation_file = pd.read_csv(annotation_file, sep="\t")
        matrix_file = pd.read_csv(matrix_file, sep="\t")
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

        processed_data = self.msstats.dataProcess(msstats_data, use_log_file=True)
        result, *_ = self.msstats.groupComparison(constrast_matrix=matrix_file, data=processed_data, use_log_file=True)
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



