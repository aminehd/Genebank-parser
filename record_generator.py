from os import truncate
from sys import path
from typing import Iterator
from pathlib import Path
import seq_parser
import re
import pandas as pd

import os
# import records either sequence type or feature type
# this class can be used to spit records like a cds record and save the records of interest in a data frame
# can pass column names that you want in data frame
# example extra_cols = ['coded_by'] ( db_xref is a defualt col)
#      CDS             1..406
#                      /gene="PSMC5"
#                      /coded_by="NM_001293084.1:1..1221"
#                      /db_xref="GeneID:101445299"
# returns datafreame of minimum of batch size and the remained items in the file
# iterate over all data frame by generator style iterating on Recorde_generator.next_batch
import gzip
class Recorde_generator:
    def __init__(self, files_name, batch_size_features, batch_size_origins, extra_cols = []):
        # self.files_dir = files_dir
        self.batch_size_features = batch_size_features
        self.batch_size_origins = batch_size_origins
        self.origins = []
        self.features = []
        self.current_accession = None
        self.file_name = files_name
        self.extra_cols = extra_cols
        self.num_of_features = 0
        self.num_of_records = 0

    def get_lines_from_files(self):
        pathlist = Path(self.files_dir).rglob("*.seq")


        for file in pathlist:
            f = open(file, "r")
            lines = f.readlines()
            f.close()
            yield  (lines, file)

    def next_batch(self):
        if ( self.file_name.endswith("gz")):
            openner = gzip.open
            stringer = lambda x : str(x, 'utf-8')
        else:
            openner = open

            stringer = lambda x : str(x, 'latin1')

        with openner(self.file_name, 'rb') as f:
            lines = f.readlines()
            file_name = self.file_name

            parser = seq_parser.Parser_HEAD(self.extra_cols)
            for line in lines:
                line = stringer(line)
                self.sanity_check(line)
                for el in parser.parse_line(line):
                    tag, record = el
                    if tag == "not compelete" and record and "ACCESSION" in record.keys():
                        self.current_accession = record["ACCESSION"]
                    if tag == "feature":
                        record["ACCESSION"] = self.current_accession
                        record["file_name"] = re.sub(r'.*[^a-z0-9]([a-z0-9]+)\.seq.*', r'\1',  self.file_name)
                        self.features.append(record)
                        if(len(self.features) == self.batch_size_features):
                            yield  ("features", self.convert_to_df(self.features))
                            self.features = []
                    if tag == "origin":
                        self.origins.append(record)
                        record["file_name"] = re.sub(r'.*[^a-z0-9]([a-z0-9]+)\.seq.*', r'\1',  self.file_name)
                        if(len(self.origins) == self.batch_size_origins):
                            yield  ("origins", self.origins)
                            self.origins = []


        if(len(self.features) > 0):
            yield  ("features", self.convert_to_df(self.features))
        if(len(self.origins) > 0):
            yield  ("origins", self.origins)

        yield ("sanity", (self.num_of_features, self.num_of_records, re.sub(r'.*[^a-z0-9]([a-z0-9]+)\.seq\.gz', r'\1',  self.file_name) ))
    def sanity_check(self, line):
        pattern = r"     (\S+)(\s+)(\S+)(.*)"
        x = False
        m = re.match(pattern, line)
        if m:
            if(len(m.group(1) + m.group(2)) == 16 ):
                x = True
                self.num_of_features += 1

        if(re.match("ACCESSION(\s*)", line)):
            self.num_of_records += 1

        return x
    def convert_to_df(self, batch):
        pd_batch = pd.DataFrame(batch)

        try:
            # reorder columns

            #process the columns
            pd_batch["translation"] = pd_batch["translation"].apply(self.combined_translations)
            pd_batch["db_xref"] = pd_batch["db_xref"].apply(lambda x: "" if len(x) == 0 else x[0])
            pd_batch["ACCESSION"] = pd_batch["ACCESSION"].apply(lambda x: x.strip())
            pd_batch["interval"] = pd_batch["interval"].apply(lambda x: x.strip())

            if('note' in pd_batch):
                print( "note is there")

                pd_batch["note"] = pd_batch["note"].apply(self.combined_notes)
            return pd_batch
        except:
            print("failed")
    def combined_translations(self, translations):
        if (len(translations) > 0):
            pattern = r"([^A-Z]*)([A-Z]+)([^A-Z]*)"
            captal_getter = lambda x: "" if re.match( pattern, x) == None else re.match(pattern, x).group(2)
            tc = list(map(captal_getter, translations))
            a = "".join(tc)
            return a
        return  ""

    def combined_notes(self, notes):
        if (len(notes) > 0):
            a = "".join(notes)
            return a
        return  ""

