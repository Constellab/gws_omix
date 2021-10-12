# LICENSE
# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
import json
from gws_core import File, Settings, BaseTestCase, TaskTester, GTest
from gws_omix import BlastEC

class TestBlastEC(BaseTestCase):
    
    async def test_BlastEC(self):
        settings = Settings.retrieve()
        data_dir = settings.get_variable("gws_omix:testdata_dir")
        file = File()
        file.path = os.path.join(data_dir, "fasta_subject.faa")

        settings = Settings.retrieve()
        testdata_dir = settings.get_variable("gws_omix:testdata_dir")

        # Running BlastEC.py
        tester = TaskTester(
            params = {'taxonomy': 'fungi', "uniprot_db_dir": os.path.join(testdata_dir, "uniprot_kb")},
            inputs = {'fasta_file': file},
            task_type = BlastEC
        )
        outputs = await tester.run()
        blast_ec_file = outputs['filtered_blast_ec_file']
        result_content = blast_ec_file.read()


        # Get the expected file output       
        expected_file_path = os.path.join(data_dir, "results.test.blast.csv")

        with open(expected_file_path, 'r') as fp:
            expected_result_content = fp.read()
            print("----")
            print(result_content)
            print("----")
            print(expected_result_content)
            print("----")
            # Comparing results
            self.assertEqual( result_content, expected_result_content  )


        text_view = blast_ec_file.view_head_as_raw_text()
        text_view_dict = text_view.to_dict(page=1, page_size=500)
        print(json.dumps(text_view_dict, indent=2))
        self.assertEqual( text_view_dict["total_number_of_pages"], 5 )


### Oldies ###


#     def test_BlastEC(self):
        
#         data_dir = settings.get_dir("omix:testdata_dir")
       
#         file_path = os.path.join(data_dir, "fungi.uniprotKB.faa") 
#         with open(file_path, newline='') as f:
#             fastadb = csv.reader(f, delimiter='\t')
            
#         file_path = os.path.join(data_dir, "fungi.uniprotKB.tab") 
#         with open(file_path, newline='') as f:
#             tab = csv.reader(f, delimiter='\t')

#         file_path = os.path.join(data_dir, "fasta_subject.faa") 
#         with open(file_path, newline='') as f:
#             fastasubject = csv.reader(f, delimiter='\t')

#         file = File()
#         file.path = file_path
#         file.save()
#         blt = BlastEC()
#         blt.input["fasta_file"] = file
# #        blt.input["fasta_file"] = file_path       

#         def _on_end(*args, **kwargs):
#             f = blt.output["filtered_blast_file"]
            
#             result_content = f.read()
            
            
#             print("----")
#             print(result_content)
            
#             print("----")
            
#             file_path = os.path.join(data_dir, "blast.output.pipe_sep.txt")
#             with open(file_path) as fp:
#                 expected_result_content = fp.read()
#                 print(expected_result_content)
#                 self.assertEqual( result_content, expected_result_content  )
            
        
#         e = blt.create_experiment(study=GTest.study, user=GTest.user)
#         e.on_end(_on_end)
#         asyncio.run( e.run() )