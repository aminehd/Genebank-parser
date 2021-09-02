import re
import pprint
pp = pprint.PrettyPrinter(indent=4)
import finders
import sys
# TO DO
# pep formating
# test strategy
# error handling (logs line number, current line)
# stat collection while parsing -- putting checks and balances while testing

# reads line by line, yields a record if it's the end of a record
# record can be a feature with accession id or a sequence with accession id 

class Parser_HEAD:
    def __init__(self, extra_cols = []):
        self.extra_cols = extra_cols
        self.num_line = 0
        self.command = ""
        self.records = []
        self.sequences = []
        self.monitor = {"sequence_starts_NAN": [], "featuresempty": [], "originEmpty": [], "line_number": 0, "visited_area": [], "unexpected_area_end": {"ids_end": [], "features_end": [], "sequence_end": []}}
        dae_instance = finders.Date_Accession_Version_Extractor(self.monitor)
        dae_instance.start_over("dummy")
        ff_instance = finders.FeatureFinder(self.monitor, self.extra_cols)
        seq_instance = finders.Find_Seuence(self.monitor)
        self.region_data_extractor = {"get_id": dae_instance , "get_features":  ff_instance, "sequence": seq_instance}
        self.command = "get_id"
        self.data_extractor = None
        self.currentRecord = None
    def parse_line(self, line):

        self.monitor["line_number"] += 1

        # we have three extractors- id, features, and origin
        self.data_extractor = self.get_extractor()

        # use same extractor
        if self.data_extractor.is_next_same(line):
            for record in self.data_extractor.process_next_line(line):
                yield record
        else:

            for record in self.wrap_up():
                yield record

            # find the instance of extractor class based on the command
            self.data_extractor = self.get_extractor()
            self.data_extractor.start_over(line)
            pass

    def get_extractor(self):
        return self.region_data_extractor[self.command]
    def wrap_up(self):
        item = None
        for x in self.data_extractor.wrap_up():
            item = x


        if self.command == "get_id":
            self.command = "get_features"
            x, self.currentRecord = item
        elif self.command == "get_features":
            if not item or len(item) <= 1 or item[0] != "feature":
                self.monitor["featuresempty"].append(self.monitor["line_number"])

            self.currentSeq = []
            self.command = "sequence"

        elif self.command == "sequence":
            if not item or len(item) <= 1 or item[0] != "origin":
                self.monitor["originEmpty"].append(self.monitor["line_number"])
            else:
                self.currentRecord["ORIGIN"] = item[1]
                item = (item[0], self.currentRecord)
            self.command = "get_id"

        yield item

    def store_dates_as_pickle(region, result):
        pass
