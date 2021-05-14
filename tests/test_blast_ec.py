# LICENSE
# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import re
import sys
import asyncio
import os, csv
import unittest

from gws.file import File
from gws.model import Protocol, Study, Experiment, Resource
from gws.settings import Settings
settings = Settings.retrieve()

from omix.blast import BlastEC
from gws.unittest import GTest

#from gena.network import Network
#from gena.context import Context
#from gena.biomodel import Biomodel
#from gena.fba import FluxAnalyzer


class TestBlastEC(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        Experiment.drop_table()
        Resource.drop_table()
        File.drop_table()        
        Protocol.drop_table()
        Study.drop_table()
        BlastEC.drop_table()
        GTest.init()
        pass

    @classmethod
    def tearDownClass(cls):
        Experiment.drop_table()
        Protocol.drop_table()
        Resource.drop_table()
        File.drop_table()
        Study.drop_table()
        BlastEC.drop_table()
        pass

    def test_BlastEC(self):
        
        data_dir = settings.get_dir("omix:testdata_dir")
       
        file_path = os.path.join(data_dir, "fungi.uniprotKB.faa") 
        with open(file_path, newline='') as f:
            fastadb = csv.reader(f, delimiter='\t')
            
        file_path = os.path.join(data_dir, "fungi.uniprotKB.tab") 
        with open(file_path, newline='') as f:
            tab = csv.reader(f, delimiter='\t')

        file_path = os.path.join(data_dir, "fasta_subject.faa") 
        with open(file_path, newline='') as f:
            fastasubject = csv.reader(f, delimiter='\t')

        file = File()
        file.path = file_path
        file.save()
        blt = BlastEC()
        blt.input["fasta_file"] = file
#        blt.input["fasta_file"] = file_path       

        def _on_end(*args, **kwargs):
            f = blt.output["filtered_blast_file"]
            
            result_content = f.read()
            
            
            print("----")
            print(result_content)
            
            print("----")
            
            file_path = os.path.join(data_dir, "blast.output.pipe_sep.txt")
            with open(file_path) as fp:
                expected_result_content = fp.read()
                print(expected_result_content)
                self.assertEqual( result_content, expected_result_content  )
            
        
        e = blt.create_experiment(study=GTest.study, user=GTest.user)
        e.on_end(_on_end)
        asyncio.run( e.run() )