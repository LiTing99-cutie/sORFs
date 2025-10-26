import sys
import ast

if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} <offdict_file>")
    sys.exit(1)

input_file = sys.argv[1]

with open(input_file, 'r') as f:
    for line in f:
        if 'offdict' in line:
            dict_str = line.split('=', 1)[1].strip()
            offdict = ast.literal_eval(dict_str)
            for k, v in offdict.items():
                if isinstance(k, int):
                    print(f"{k}\t{v}") 