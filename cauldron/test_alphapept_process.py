import traceback
from unittest import TestCase
from cauldron.alphapept_process import preprocess_data, load_file, main, differential_analysis, \
    differential_analysis_limma
import pandas as pd
from click.testing import CliRunner

from cauldron.differential_analysis import set_up_R_HOME


class Test(TestCase):
    def test_alphapept_process(self):
        annotation = pd.read_csv("test/annotation.txt", sep="\t")
        data = load_file(file_path="test/imputed.data.txt", metadata_path="test/annotation.txt", engine="generic", index="Precursor.Id", intensity_column=list(annotation["Sample"]))
        data = preprocess_data(data, impute="knn", normalize="zscore")
        print(data.mat.columns)
        print(data.mat.index)
        result = differential_analysis(data, "BCA_LT-IP", "BCA_LT-MockIP")

        print(result)
        differential_analysis(data, "BCA_LT-IP", "BCA_LT-MockIP", method="sam")

class TestLimma(TestCase):
    def test_alphapept_process(self):
        annotation = pd.read_csv("test/annotation.txt", sep="\t")
        set_up_R_HOME(r'D:\WebstormProjects\cauldron\app\bin\win\R-Portable\App\R-Portable')
        data = load_file(file_path="test/imputed.data.txt", metadata_path="test/annotation.txt", engine="generic", index="Precursor.Id", intensity_column=list(annotation["Sample"]))
        data = preprocess_data(data, impute="knn", normalize="quantile")
        comparison = "D:\\cauldrontest\\d2f93641-a596-4771-8854-3ef821069456\\comparison.txt"
        comparison_matrix = pd.read_csv(comparison, sep="\t")
        results = differential_analysis_limma(data, comparison_matrix, r_home=r'D:\WebstormProjects\cauldron\app\bin\win\R-Portable\App\R-Portable', log2=True)

class TestMainFunction(TestCase):
    def setUp(self):
        annotation = pd.read_csv("test/annotation.txt", sep="\t")
        self.runner = CliRunner()
        self.file_path = 'test/imputed.data.txt'
        self.metadata_path = 'test/annotation.txt'
        self.engine = 'generic'
        self.index = 'Protein.Group'
        self.intensity_column = ",".join(list(annotation["Sample"]))
        self.impute = 'knn'
        self.normalize = 'zscore'
        self.data_completeness = 0.3
        self.group_A = 'BCA_LT-IP'
        self.group_B = 'BCA_LT-MockIP'
        self.method = 'ttest'
        self.p_value = 0.05
        self.permutation = 10
        self.output_folder = 'test'
        self.merge_column = 'Genes'

    def test_main(self):
        result = self.runner.invoke(main, [
            '--file_path', self.file_path,
            '--metadata_path', self.metadata_path,
            '--engine', self.engine,
            '--index', self.index,
            '--intensity_column', self.intensity_column,
            '--impute', self.impute,
            '--normalize', self.normalize,
            '--data_completeness', self.data_completeness,
            '--group_a', self.group_A,
            '--group_b', self.group_B,
            '--method', self.method,
            '--p_value', self.p_value,
            '--permutation', self.permutation,
            '--output_folder', self.output_folder,
            '--merge_columns', self.merge_column,
            '--log2'
        ])
        print(result.output)
        print(result.exception)
        print(result.exit_code)

    def test_main2(self):
        result = self.runner.invoke(main, [
        "--file_path",
        "D:\\WebstormProjects\\angular-electron-main\\app\\examples\\diann\\imputed.data.txt",
        "--metadata_path",
        "D:\\WebstormProjects\\angular-electron-main\\app\\examples\\differential_analysis\\annotation.txt",
        "--index",
        "Precursor.Id",
        "--comparison_matrix",
        "D:\\cauldrontest\\d2f93641-a596-4771-8854-3ef821069456\\comparison.txt",
        "--output_folder",
        "D:\\cauldrontest\\d2f93641-a596-4771-8854-3ef821069456",
        "--data_completeness",
        "0.3",
        "--impute",
        "knn",
        "--normalize",
        "quantile",
        "--method",
        "welch-ttest",
        "--engine",
        "generic",
        "--merge_columns",
        "Protein.Ids,Genes",
        "--log2"])
        if result.exception:
            traceback.print_exception(type(result.exception), result.exception, result.exception.__traceback__)

    def test_main_limma(self):
        result = self.runner.invoke(main, [
        "--file_path",
        "D:\\WebstormProjects\\angular-electron-main\\app\\examples\\diann\\imputed.data.txt",
        "--metadata_path",
        "D:\\WebstormProjects\\angular-electron-main\\app\\examples\\differential_analysis\\annotation.txt",
        "--index",
        "Precursor.Id",
        "--comparison_matrix",
        "D:\\cauldrontest\\d2f93641-a596-4771-8854-3ef821069456\\comparison.txt",
        "--output_folder",
        "D:\\cauldrontest\\d2f93641-a596-4771-8854-3ef821069456",
        "--data_completeness",
        "0.3",
        "--impute",
        "knn",
        "--normalize",
        "quantile",
        "--method",
        "limma",
        "--engine",
        "generic",
        "--merge_columns",
        "Protein.Ids,Genes",
        "--log2", "--r_home", r'D:\WebstormProjects\cauldron\app\bin\win\R-Portable\App\R-Portable'])
        if result.exception:
            traceback.print_exception(type(result.exception), result.exception, result.exception.__traceback__)