"""
Module with functionalities to implement 2D VCM
"""
from PyVCM2D.ProblemSetup import ReadYaml

class PyVCM(ReadYaml):

    def __init__(self, case=None):

        print('\n### PyVCM initialized ###')

        super().__init__(case)
