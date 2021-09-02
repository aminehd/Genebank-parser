import re

class Region_Data_Extractor:
    def __init__(self) -> None:
        pass

    def start_over(self, line):
        pass

    def is_next_same(self, line):
        indent_size, dummy1, dummy2= self.get_indent_and_parts(line)
        return indent_size != 0

    def wrap_up(self):
        pass
    def process_next_line(self, line):
        pass
    def get_indent_and_parts(self, line):
        indenting_group = re.match( r"(\s*)(\S*)(.*)", line)
        return  (len(indenting_group.group(1)), indenting_group.group(2), indenting_group.group(3))

class Date_Accession_Version_Extractor(Region_Data_Extractor):
    def __init__(self, monitor) -> None:
        self.monitor = monitor
        pass
    def start_over(self, line):
        self.monitor["visited_area"].append(1)
        self.record_info = {"ACCESSION": "", "LOCUS": "", "VERSION": ""}
        self.record = {}
    def is_next_same(self, line):
        b = super(Date_Accession_Version_Extractor, self).is_next_same(line)
        m = re.search(r"\A(\s*)\bFEATURES\b(\s)*Location", line)
        if m and b :
            self.monitor["unexpected_area_end"]["ids_end"].append(self.monitor["line_number"])
        return not m
    def process_next_line(self, line):
        indent, part1, part2  = super(Date_Accession_Version_Extractor, self).get_indent_and_parts(line)
        if( indent == 0 ):
            if part1 in self.record_info.keys():
                self.record_info[part1] = part2
        return
        yield

    def wrap_up(self):
        yield ("not compelete", self.record_info)

class Find_Seuence(Region_Data_Extractor):
    def __init__(self, monitor):
        self.monitor = monitor
    def start_over(self, line):
        self.monitor["visited_area"].append(3)

        self.sequences = []
        if( not re.match(r"(\s*)\bORIGIN\b(\s*)", line)):
            self.sequences.append(line)

    def is_next_same(self, line):
        b = super(Find_Seuence, self).is_next_same(line)

        return b
    def process_next_line(self, line):
        m = re.search(r"(\d+).*", line)
        if not m:
            self.monitor["sequence_starts_NAN"].appned(self.monitor["line_number"])
        self.sequences.append(line)
        return
        yield
    def wrap_up(self):
        yield ("origin", self.sequences)

class FeatureFinder(Region_Data_Extractor):
    def __init__(self, monitor, extra_cols = []) -> None:
        self.features = []
        self.monitor = monitor
        self.extra_cols = extra_cols

        pass

    def start_over(self, line):
        self.features = []
        self._next_is_feature = True
        self.indent = 10000
        self.last_feature = "non_existing"
        self.monitor["visited_area"].append(2)
        self.currfeautre = None
        self.just_started = True

    def is_next_same(self, line):
        b = super(FeatureFinder, self).is_next_same(line)
        if not b:
            m = re.match(r"(\s*)\bORIGIN\b(\s*)", line)
            if m == None:
                self.monitor["unexpected_area_end"]["features_end"].append(self.monitor["line_number"])
        return  b

    def wrap_up(self):
        if self.currfeautre:
            yield ("feature", self.currfeautre)
        # return self.features

    def process_next_line(self, line):
        if(self.just_started):
            self.just_started = False
            self.level_one_indent, xx, xxx = self.get_indent_and_parts(line)
        indent_size, part1, part2 = super(FeatureFinder, self).get_indent_and_parts(line)
        m2 = None
        if part2:


            m2 = re.search(r"(\d+)", part2)
        if indent_size == self.level_one_indent and m2:
            if self.currfeautre:
                yield ("feature", self.currfeautre)

            set_dict = {
                "type": part1,
                "interval":part2,
                "protein_id": "",
                "translation": [],
                "gene": "",
                "db_xref": [],
                "product": "",
                "locus_tag": ""
            }
            input_dict = {k[1:]: [] if k[0] == 'M' else "" for k in self.extra_cols }
            self.currfeautre = {**set_dict, **input_dict}
            self.indent = indent_size
            self.features.append(self.currfeautre)
            self.linesoftranslation = None
            return

        self.indent = indent_size

        feature_line = re.match(r"(\s*)(\/)([^=]*)(\=)(.*)", line)
        if(feature_line):
            feature_name, feature_value = feature_line.group(3), feature_line.group(5)
            self.last_feature = feature_name
            if feature_line.group(3) in self.currfeautre.keys():
                if hasattr(self.currfeautre[feature_name] , 'append'):
                    self.currfeautre[feature_name].append(feature_value)
                else:
                    self.currfeautre[feature_name] = feature_value
        elif self.last_feature in self.currfeautre.keys():
            if hasattr(self.currfeautre[self.last_feature],"append"):
                self.currfeautre[self.last_feature].append(line)
            else:
                self.currfeautre[self.last_feature] = line

        self.features[-1] = self.currfeautre
        ##search line to a tree, each leaf has an action for current and next command

