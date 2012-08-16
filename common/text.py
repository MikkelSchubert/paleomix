import string

def padded_table(table, min_padding = 4):
    sizes = [0] * len(table[0])
    for row in table:
        for (ii, field) in enumerate(row):
            sizes[ii] = max(sizes[ii], len(str(field)))
    sizes = [(size + min_padding) for size in sizes]

    for row in table:
        padded = []
        for (ii, field) in enumerate(row):
            padded.append(string.ljust(str(field), sizes[ii]))
        yield "".join(padded).rstrip()

