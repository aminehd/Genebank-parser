# use record generator to read a 100 features and call it
import record_generator
import pprint
pp = pprint.PrettyPrinter( indent=4)
import pandas as pd

if __name__ == "__main__":
    pd.set_option('display.max_rows', 500)
    pd.set_option('display.max_columns', 4)
    pd.set_option('display.width', 1000)
    record_gen = record_generator.Recorde_generator("../data/gbbct1.seq", 100, 100, [])
    for tag, batch in record_gen.next_batch():
        if (tag == "features"):
            with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
                print(batch)
            break