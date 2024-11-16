import sys

input_file = sys.argv[1]
output_file= sys.argv[2]

with open(input_file, "r") as input_f, open(output_file, "w") as output_f:
    for line in input_f:
        if line.startswith(">"):
            header_components = line.strip().split(".")
            header_components = header_components[:1]
            modified_header =" ".join(header_components)
            output_f.write(modified_header + "\n")
        else:
            output_f.write(line)
